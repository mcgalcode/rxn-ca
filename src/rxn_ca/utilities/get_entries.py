from rxn_network.entries.entry_set import GibbsEntrySet
from rxn_network.entries.utils import process_entries
from mp_api.client import MPRester


from typing import List

def get_entries(chem_sys: str,
                stability_cutoff: float = 0.1,
                ensure_phases: List[str] = [],
                custom_entries: List = [],
                **kwargs) -> GibbsEntrySet:
    # Note: custom entries should be retrieved from MP using the
    # additional_criteria={"thermo_types": ["GGA_GGA+U"]} option
    # First we enumerate entries
    with MPRester() as mpr:
        mp_computed_struct_entries = mpr.get_entries_in_chemsys(
            elements=chem_sys,
            additional_criteria={"thermo_types": ["GGA_GGA+U"]},
        )
        
    all_entries = [*custom_entries, *mp_computed_struct_entries]
    entry_set = process_entries(
        all_entries,
        300, 
        stability_cutoff,
        formulas_to_include=ensure_phases,
        calculate_e_above_hulls=True,
        **kwargs
    )
    return entry_set