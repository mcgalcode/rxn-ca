"""Jobflow jobs for rxn-ca simulations.

These jobs can be used standalone or composed into flows for running
rxn-ca simulations in workflow frameworks.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional
import json

from jobflow import job
from monty.json import MontyEncoder, MontyDecoder

from .schemas import ReactionLibraryData, SimulationOutput


def _build_reaction_library(
    chemical_system: str,
    temperatures: List[float],
    ensure_phases: List[str] = None,
    metastability_cutoff: float = 0.1,
    exclude_theoretical: bool = True,
) -> tuple:
    """Build phase set and reaction library for a chemical system.

    Args:
        chemical_system: Chemical system string (e.g., "Ba-Ti-O")
        temperatures: List of temperatures (K) to score reactions at
        ensure_phases: List of phase formulas that MUST be included
        metastability_cutoff: Energy above hull cutoff for phases
        exclude_theoretical: Whether to exclude theoretical phases

    Returns:
        Tuple of (phase_set, reaction_lib)
    """
    from rxn_ca.utilities.get_entries import get_entries
    from rxn_ca.phases import SolidPhaseSet
    from rxn_ca.utilities.get_scored_rxns import get_scored_rxns
    from rxn_network.enumerators.basic import BasicEnumerator
    from rxn_network.enumerators.minimize import MinimizeGibbsEnumerator
    from rxn_network.enumerators.utils import run_enumerators

    # Get entries from Materials Project
    entries = get_entries(
        chem_sys=chemical_system,
        metastability_cutoff=metastability_cutoff,
        ensure_phases=ensure_phases or [],
        exclude_theoretical_phases=exclude_theoretical,
    )
    print(f"Got {len(entries)} entries for {chemical_system}")
    if ensure_phases:
        print(f"  (ensured inclusion of: {ensure_phases})")

    # Create phase set
    phase_set = SolidPhaseSet.from_entry_set(entries)

    # Enumerate reactions
    enumerators = [MinimizeGibbsEnumerator(), BasicEnumerator()]
    rxn_set = run_enumerators(enumerators, entries)
    print(f"Enumerated {len(rxn_set)} reactions")

    # Compute reactions at all temperatures
    temp_rxn_mapping = rxn_set.compute_at_temperatures(temperatures)

    # Score reactions
    reaction_lib = get_scored_rxns(
        rxn_set,
        temps=temperatures,
        phase_set=phase_set,
        rxns_at_temps=temp_rxn_mapping,
        parallel=True,
    )

    return phase_set, reaction_lib


@job
def setup_reaction_library(
    chemical_system: str,
    temperatures: List[float],
    ensure_phases: List[str] = None,
    metastability_cutoff: float = 0.1,
    exclude_theoretical: bool = True,
    save_to_file: bool = True,
) -> ReactionLibraryData:
    """Set up phase set and reaction library for a chemical system.

    This job fetches thermodynamic data from Materials Project, enumerates
    possible reactions, and scores them at each temperature. This is typically
    the most expensive step and can be shared across multiple simulations
    in the same chemical system.

    Args:
        chemical_system: Chemical system string (e.g., "Ba-Ti-O")
        temperatures: List of temperatures (K) to score reactions at
        ensure_phases: List of phase formulas that MUST be included
            even if they would otherwise be filtered out (e.g., phases
            known to exist from experimental observations)
        metastability_cutoff: Energy above hull cutoff for phases
        exclude_theoretical: Whether to exclude theoretical phases
        save_to_file: If True, save reaction library to a JSON file

    Returns:
        ReactionLibraryData with phase set, reaction library, and metadata
    """
    phase_set, reaction_lib = _build_reaction_library(
        chemical_system,
        temperatures,
        ensure_phases,
        metastability_cutoff,
        exclude_theoretical,
    )

    # Optionally save reaction library to file
    reaction_library_path = None
    if save_to_file:
        filename = f"reaction_library_{chemical_system.replace('-', '_')}.json"
        reaction_library_path = str(Path.cwd() / filename)
        with open(reaction_library_path, "w") as f:
            json.dump(reaction_lib.as_dict(), f, cls=MontyEncoder)
        print(f"Saved reaction library to {reaction_library_path}")

    return ReactionLibraryData(
        phase_set_dict=phase_set.as_dict(),
        reaction_library_dict=reaction_lib.as_dict(),
        chemical_system=chemical_system,
        temperatures=temperatures,
        phases_available=list(phase_set.phases),
        reaction_library_path=reaction_library_path,
    )


@job
def run_simulation(
    recipe: "ReactionRecipe",
    reaction_library_data: ReactionLibraryData = None,
    chemical_system: str = None,
    ensure_phases: List[str] = None,
    metastability_cutoff: float = 0.1,
    save_to_file: bool = True,
    metadata: Dict[str, Any] = None,
    compress_freq: int = 100,
    num_frames: int = None,
    analysis_only: bool = False,
) -> SimulationOutput:
    """Run an rxn-ca simulation.

    Args:
        recipe: ReactionRecipe specifying reactants, heating schedule, etc.
        reaction_library_data: Pre-computed reaction library from setup_reaction_library.
            If not provided, will build one from scratch (requires chemical_system).
        chemical_system: Chemical system string, required if reaction_library_data
            not provided
        ensure_phases: Phases to ensure are included, used if building library
            from scratch
        metastability_cutoff: Energy above hull cutoff, used if building library
            from scratch
        save_to_file: If True, save the full result doc to a JSON file
        metadata: Optional user-provided metadata for tagging/provenance
        compress_freq: If set, store a full state snapshot every compress_freq
            steps instead of per-step diffs. Cannot be combined with num_frames.
            In analysis_only mode this sets the observation cadence instead.
        num_frames: If set, store roughly this many full state snapshots over
            the run. Cannot be combined with compress_freq.
        analysis_only: If True, retain only the reduced per-phase-volume view of
            each frame (via a PhaseVolumeObserver). The saved result doc is much
            smaller, but only volume-derived analyses are available afterward.

    Returns:
        SimulationOutput with analyzed results and file references
    """
    from rxn_ca.phases import SolidPhaseSet
    from rxn_ca.core.recipe import ReactionRecipe
    from rxn_ca.utilities.parallel_sim import run_sim_parallel
    from rxn_ca.utilities.single_sim import run_single_sim
    from rxn_ca.analysis import BulkReactionAnalyzer
    from rxn_ca.analysis.reaction_step_analyzer import AnalysisMode, AnalysisQuantity
    from rxn_ca.analysis.visualization.phase_trace_calculator import (
        PhaseTraceCalculator,
        PhaseTraceConfig,
    )
    from rxn_ca.reactions import ReactionLibrary

    recipe_dict = recipe.as_dict()

    # Set up phase set and reaction library
    if reaction_library_data is not None:
        if reaction_library_data.reaction_library_path:
            # The library was saved to a file (e.g. grid runs ship an empty-dict
            # stub to keep the FireWorks spec small). Load the full library from
            # disk; it carries the phase set as `reaction_lib.phases`.
            reaction_lib = ReactionLibrary.from_file(reaction_library_data.reaction_library_path)
            phase_set = reaction_lib.phases
        elif reaction_library_data.reaction_library_dict:
            phase_set = SolidPhaseSet.from_dict(reaction_library_data.phase_set_dict)
            reaction_lib = ReactionLibrary.from_dict(reaction_library_data.reaction_library_dict)
        else:
            raise ValueError(
                "reaction_library_data must include reaction_library_path or reaction_library_dict"
            )
        chem_sys = reaction_library_data.chemical_system
        reaction_library_path = reaction_library_data.reaction_library_path
    else:
        if chemical_system is None:
            raise ValueError("chemical_system required when reaction_library_data not provided")
        all_temps = recipe.heating_schedule.all_temps
        phase_set, reaction_lib = _build_reaction_library(
            chemical_system, all_temps, ensure_phases, metastability_cutoff
        )
        chem_sys = chemical_system
        reaction_library_path = None

    # Reference the library by path (when one exists) so the saved result doc
    # does not embed the full library.
    if num_frames is not None:
        compress_freq = None

    if recipe.num_realizations > 1:
        result_doc = run_sim_parallel(
            recipe=recipe,
            reaction_lib=reaction_lib,
            phase_set=phase_set,
            compress_freq=compress_freq,
            num_frames=num_frames,
            analysis_only=analysis_only,
            reaction_lib_path=reaction_library_path,
        )
    else:
        result_doc = run_single_sim(
            recipe=recipe,
            reaction_lib=reaction_lib,
            phase_set=phase_set,
            compress_freq=compress_freq,
            num_frames=num_frames,
            analysis_only=analysis_only,
            reaction_lib_path=reaction_library_path,
        )

    # Optionally save result doc to file
    result_doc_path = None
    if save_to_file:
        filename = f"result_doc_{chem_sys.replace('-', '_')}.json"
        result_doc_path = str(Path.cwd() / filename)
        with open(result_doc_path, "w") as f:
            json.dump(result_doc.as_dict(), f, cls=MontyEncoder)
        print(f"Saved result doc to {result_doc_path}")

    # Analyze results
    analyzer = BulkReactionAnalyzer.from_result_doc(result_doc)

    # Get final molar amounts
    final_molar_amounts = analyzer.get_all_absolute_molar_amounts(
        analyzer.last_loaded_step_idx
    )

    # Build trajectory: molar amounts at each step.
    #
    # get_all_absolute_molar_amounts(step) returns ONLY the phases present at that
    # step (absent phases are omitted). Appending per-step like that gives each
    # phase a list whose length is "number of steps it happened to exist for",
    # with NO record of WHICH steps -- so the per-phase trajectories lose their
    # common step alignment and a consumed precursor becomes indistinguishable
    # from a late-forming product. That produced nonsense mass fractions
    # downstream (e.g. a single phase reading 1.0 at step 0).
    #
    # Use the same aligned builder the plotter uses (PhaseTraceCalculator): it
    # walks every step and 0-fills phases absent at that step, keeping all phases
    # on the common step axis. matter_phases=None keeps gas phases (e.g. CO2),
    # matching the previous get_all_absolute_molar_amounts behaviour;
    # minimum_required_prevalence=0.0 keeps every phase that ever appears.
    step_indices: List[int] = list(analyzer.loaded_step_idxs)
    temp_trajectory: List[float] = [
        recipe.heating_schedule.temp_at(step_idx) for step_idx in step_indices
    ]

    trace_calc = PhaseTraceCalculator(
        analyzer.loaded_step_groups, analyzer.step_analyzer
    )
    molar_traces = trace_calc.get_general_traces(
        PhaseTraceConfig(minimum_required_prevalence=0.0),
        AnalysisQuantity.MOLES,
        AnalysisMode.ABSOLUTE,
        matter_phases=None,
    )
    molar_trajectory: Dict[str, List[float]] = {
        t.name: [float(y) for y in t.ys] for t in molar_traces
    }

    return SimulationOutput(
        final_molar_amounts=final_molar_amounts,
        molar_amounts_trajectory=molar_trajectory,
        temperature_trajectory=temp_trajectory,
        step_indices=step_indices,
        phase_set_dict=phase_set.as_dict(),
        reaction_library_path=reaction_library_path,
        result_doc_path=result_doc_path,
        recipe_dict=recipe_dict,
        chemical_system=chem_sys,
        num_realizations=recipe.num_realizations,
        metadata=metadata,
    )
