from rxn_network.entries.entry_set import GibbsEntrySet
from rxn_network.entries.utils import process_entries
from mp_api.client import MPRester
from pymatgen.core import Composition
from rxn_ca.phases import SolidPhaseSet
from rxn_ca.phases.utils import remove_phases_from_entry_set, remove_theoretical_phases


from typing import List

def get_entries(chem_sys: str,
                metastability_cutoff: float = 0.1,
                ensure_phases: List[str] = [],
                custom_entries: List = [],
                thermo_types: List[str] = ["GGA_GGA+U"],
                filter_at_temperature: int = 300,
                exclude_theoretical_phases: bool = False,
                exclude_phases: List[str] = [],
                **kwargs) -> GibbsEntrySet:
    # Note: custom entries should be retrieved from MP using the
    # additional_criteria={"thermo_types": ["GGA_GGA+U"]} option
    # First we enumerate entries
    with MPRester() as mpr:
        mp_computed_struct_entries = mpr.get_entries_in_chemsys(
            elements=chem_sys,
            additional_criteria={"thermo_types": thermo_types},
        )

    # exclude_phases are dropped (every entry of the formula) before the hull is
    # built, as if Materials Project had none: a phase that DFT places too low
    # no longer pushes other phases past the metastability cutoff. An excluded
    # phase wins over ensure_phases.
    drop = {Composition(f).reduced_formula for f in exclude_phases}
    ensure_phases = [f for f in ensure_phases if Composition(f).reduced_formula not in drop]
    all_entries = [e for e in [*custom_entries, *mp_computed_struct_entries]
                   if e.composition.reduced_formula not in drop]
    entry_set = process_entries(
        all_entries,
        300, 
        e_above_hull=1.0,
        formulas_to_include=ensure_phases,
        calculate_e_above_hulls=True,
        filter_at_temperature=filter_at_temperature,
        **kwargs
    )

    stable_entries = entry_set.filter_by_stability(metastability_cutoff)

    if exclude_theoretical_phases:
        solid_phase_set = SolidPhaseSet.from_entry_set(stable_entries)
        filtered = remove_theoretical_phases(stable_entries, solid_phase_set)

        # Re-add any ensure_phases that were removed by the theoretical filter
        # (they may be marked theoretical in MP but we know they're real)
        if ensure_phases:
            # Normalize to reduced formulas for consistent matching (e.g., NH4Cl -> H4NCl)
            ensure_set = {Composition(f).reduced_formula for f in ensure_phases}
            filtered_formulas = {e.composition.reduced_formula for e in filtered}
            for e in stable_entries:
                formula = e.composition.reduced_formula
                if formula in ensure_set and formula not in filtered_formulas:
                    filtered.add(e)

        return filtered
    else:
        return stable_entries


