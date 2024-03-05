from rxn_network.entries.entry_set import GibbsEntrySet
from .solid_phase_set import SolidPhaseSet

from typing import List

def remove_phases_from_entry_set(eset: GibbsEntrySet, phases: List[str]):
    eset_copy = eset.copy()
    for p in phases:
        entry = eset_copy.get_min_entry_by_formula(p)
        if entry is not None:
            eset_copy.discard(entry)
    
    return eset_copy

def remove_theoretical_phases(eset: GibbsEntrySet, solid_phase_set: SolidPhaseSet) -> GibbsEntrySet:
    filtered_entry_set = eset.copy()
    for p in solid_phase_set.phases:
        if p != 'Free Space':
            e = filtered_entry_set.get_min_entry_by_formula(p)
            if solid_phase_set.is_theoretical(p):
                filtered_entry_set.discard(e)
    
    return filtered_entry_set