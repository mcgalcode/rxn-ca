from rxn_network.entries.entry_set import GibbsEntrySet
from rxn_network.entries.utils import process_entries
from mp_api.client import MPRester
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
                **kwargs) -> GibbsEntrySet:
    # Note: custom entries should be retrieved from MP using the
    # additional_criteria={"thermo_types": ["GGA_GGA+U"]} option
    # First we enumerate entries
    with MPRester() as mpr:
        mp_computed_struct_entries = mpr.get_entries_in_chemsys(
            elements=chem_sys,
            additional_criteria={"thermo_types": thermo_types},
        )
        
    all_entries = [*custom_entries, *mp_computed_struct_entries]
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
        return remove_theoretical_phases(stable_entries, solid_phase_set)
    else:
        return stable_entries


