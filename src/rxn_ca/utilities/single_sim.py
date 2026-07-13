from .heating_schedule_runner import MeltAndRegrindMultiRunner, HeatingScheduleRunner
from ..core.recipe import ReactionRecipe
from ..reactions import ReactionLibrary

from ..computing.schemas.ca_result_schema import RxnCAResultDoc

from pylattica.core import Simulation
from rxn_network.reactions.reaction_set import ReactionSet

from ..phases import SolidPhaseSet
from ..core.reaction_controller import ReactionController
from ..core.liquid_swap_controller import LiquidSwapController
from ..core.reaction_calculator import ReactionCalculator
from ..core.pairwise_reaction_calculator import PairwiseReactionCalculator

_UPDATE_SCHEMES = {
    "independent": ReactionCalculator,
    "pairwise": PairwiseReactionCalculator,
}

from .get_scored_rxns import get_scored_rxns
from .setup_reaction import setup_reaction, setup_noise_reaction


def run_single_sim(recipe: ReactionRecipe,
                   base_reactions: ReactionSet = None,
                   reaction_lib: ReactionLibrary = None,
                   initial_simulation: Simulation = None,
                   phase_set: SolidPhaseSet = None,
                   existing_lib: ReactionLibrary = None,
                   compress_freq: int = None,
                   num_frames: int = None,
                   analysis_only: bool = False,
                   reaction_lib_path: str = None) -> RxnCAResultDoc:
    """Run a single simulation.

    Args:
        recipe: The reaction recipe to simulate
        base_reactions: Base reactions for scoring (required if reaction_lib not provided)
        reaction_lib: Pre-computed complete reaction library (skips all scoring)
        initial_simulation: Initial simulation state (optional)
        phase_set: Phase set for the system
        existing_lib: Existing library with some temps already scored.
            New temps will be added to this library incrementally.
        compress_freq: If set, store a full state snapshot every compress_freq
            steps instead of per-step diffs. Avoids slow reconstruction during
            analysis. Cannot be combined with num_frames.
        num_frames: If set, store roughly this many full state snapshots over
            the run. Cannot be combined with compress_freq.
        analysis_only: If True, do not retain full simulation states. Instead
            record only the reduced per-phase-volume view (via a
            PhaseVolumeObserver) at the compress_freq / num_frames cadence, and
            store it under RxnCAResultDoc.observed_results. This is dramatically
            lighter on memory and disk, but only volume-derived analyses are
            available afterward (no per-site or reaction-choice analyses).
        reaction_lib_path: If set, the result doc references the library at this
            path instead of embedding it (the library is loaded lazily on
            access). When reaction_lib is not provided, the library is loaded
            from this path to run the simulation.

    Returns:
        RxnCAResultDoc with simulation results
    """
    if reaction_lib is None and reaction_lib_path is not None:
        reaction_lib = ReactionLibrary.from_file(reaction_lib_path)

    if base_reactions is None and reaction_lib is None:
        raise ValueError("Must provide either base_reactions or reaction_lib")

    if reaction_lib is None:
        # Get required temperatures
        required_temps = recipe.heating_schedule.all_temps
        cached_temps = existing_lib.temps if existing_lib else []
        new_temps = [t for t in required_temps if t not in cached_temps]

        if new_temps:
            print("================= RETRIEVING AND SCORING REACTIONS =================")
            if existing_lib:
                print(f"    (Reusing {len(cached_temps)} cached temps, scoring {len(new_temps)} new temps)")

        reaction_lib: ReactionLibrary = get_scored_rxns(
            base_reactions,
            heating_sched=recipe.heating_schedule,
            scorer_class=recipe.get_score_class(),
            phase_set=phase_set,
            existing_lib=existing_lib,
        )

    print()
    print()
    print()

    if len(recipe.exclude_phases) > 0:
        reaction_lib = reaction_lib.exclude_phases(recipe.exclude_phases)

    if recipe.exact_phase_set is not None:
        reaction_lib = reaction_lib.limit_phase_set(recipe.exact_phase_set)

    if initial_simulation is None:

        print("================= SETTING UP SIMULATION =================")

        initial_simulation = setup_noise_reaction(
            reaction_lib.phases,
            precursor_mole_ratios = recipe.reactant_amounts,
            size = recipe.simulation_size,
            packing_fraction = recipe.packing_fraction
        )

    print(f'================= RUNNING SIMULATION =================')

    calculator_class = _UPDATE_SCHEMES.get(recipe.update_scheme)
    if calculator_class is None:
        raise ValueError(
            f"Unknown update_scheme {recipe.update_scheme!r}; "
            f"expected one of {sorted(_UPDATE_SCHEMES)}"
        )

    rxn_calculator = calculator_class(
        LiquidSwapController.get_neighborhood_from_structure(initial_simulation.structure),
        atmospheric_species=recipe.atmospheric_phases
    )

    controller = LiquidSwapController(
        initial_simulation.structure,
        rxn_calculator=rxn_calculator,
    )

    runner = HeatingScheduleRunner()

    result = runner.run_multi(
        initial_simulation,
        reaction_lib,
        recipe.heating_schedule,
        controller=controller,
        compress_freq=compress_freq,
        num_frames=num_frames,
        analysis_only=analysis_only,
    )

    # In analysis-only mode the heavy full result is discarded; only the
    # lightweight per-phase-volume observations are persisted. The final state
    # is still available from the retained live state.
    final_simulation = Simulation(result.last_step, initial_simulation.structure)
    if analysis_only:
        result_doc = RxnCAResultDoc(
            recipe=recipe,
            observed_results=result.observed_results,
            reaction_library=reaction_lib,
            reaction_library_path=reaction_lib_path,
            phases=reaction_lib.phases,
            final_simulation=final_simulation,
        )
    else:
        result_doc = RxnCAResultDoc(
            recipe=recipe,
            results=[result],
            reaction_library=reaction_lib,
            reaction_library_path=reaction_lib_path,
            phases=reaction_lib.phases,
            final_simulation=final_simulation,
        )

    return result_doc

