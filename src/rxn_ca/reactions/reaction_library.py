from __future__ import annotations

from .scored_reaction_set import ScoredReactionSet
from .scored_reaction import ScoredReaction
from ..phases.solid_phase_set import SolidPhaseSet

from monty.json import MSONable

import json
from typing import List, Dict


class ReactionLibrary(MSONable):
    """Contains a mapping of temperatures to ScoredReactionSet objects
    which contain reactions scored at the given temperature. Used in multi-stage
    reaction simulations where different stages in the reaction are run at
    different temperatures.
    """

    @classmethod
    def from_dict(cls, d):
        library = cls(
            phases = SolidPhaseSet.from_dict(d['phases'])
        )

        for t, scored_rxns in d.get('lib').items():
            rxns = [ScoredReaction.from_dict(r) for r in scored_rxns["reactions"]]
            library.add_rxns_at_temp(
                ScoredReactionSet(rxns, phase_set=library.phases),
                t
            )

        return library
    
    @classmethod
    def from_file(cls, fpath):

        with open(fpath, 'r+') as f:
            d = json.load(f)
            return cls.from_dict(d)


    def __init__(self, phases: SolidPhaseSet):
        self.lib: Dict[int, ScoredReactionSet] = {}
        self.phases = phases
        self.metadata = {}

    def add_rxns_at_temp(self, rxns: ScoredReactionSet, temp: int) -> int:
        self.lib[int(temp)] = rxns
        return temp
    
    def exclude_phases(self, phases) -> ReactionLibrary:
        lib = ReactionLibrary(self.phases)
        for t, rxns in self.lib.items():
            lib.add_rxns_at_temp(rxns.exclude_phases(phases), t)
        
        return lib
    
    def get_rxns_at_temp(self, temp: int) -> ScoredReactionSet:
        return self.lib[temp]
    
    def get_lib_from_ids(self, rxn_ids: List[int]) -> ReactionLibrary:
        deduped_ids = list(set(rxn_ids))
        lib = ReactionLibrary(self.phases)
        for t, rxns in self.lib.items():
            pruned_rxn_set = ScoredReactionSet([], lib.phases)
            for rxn_id in deduped_ids:
                rxn = rxns.get_rxn_by_id(rxn_id)
                pruned_rxn_set.add_rxn(rxn, rxn_id)
                
            lib.add_rxns_at_temp(pruned_rxn_set, t)
        return lib
    
    def add_metadata(self, rxn_id, metadata):
        if rxn_id not in self.metadata:
            self.metadata[rxn_id] = metadata
        else:
            self.metadata[rxn_id] = {**self.metadata[rxn_id], **metadata}
    
    def limit_phase_set(self, phases) -> ReactionLibrary:
        lib = ReactionLibrary(self.phases)
        for t, rxns in self.lib.items():
            lib.add_rxns_at_temp(rxns.limit_phases(phases), t)
        return lib
    
    @property
    def temps(self):
        return list(self.lib.keys())
    

    def as_dict(self):
        sup = {"@module": self.__class__.__module__, "@class": self.__class__.__name__}

        lib = {
            temp: rset.as_dict()
            for temp, rset in self.lib.items()
        }

        return {
            **sup,
            "phases": self.phases.as_dict(),
            "lib": lib,
        }

    def to_file(self, fpath):
        with open(fpath, 'w+') as f:
            json.dump(self.as_dict(), f)