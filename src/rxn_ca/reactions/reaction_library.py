from .scored_reaction_set import ScoredReactionSet
from ..phases.solid_phase_set import SolidPhaseSet
from rxn_network.reactions.reaction_set import ReactionSet

from monty.json import MSONable
from dataclasses import dataclass

from rxn_network.reactions.reaction_set import ReactionSet

import json


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
            library.add_rxns_at_temp(
                ScoredReactionSet.from_dict(scored_rxns),
                t
            )

        return library
    
    @classmethod
    def from_file(cls, fpath):

        with open(fpath, 'r+') as f:
            d = json.load(f)
            return cls.from_dict(d)



    def __init__(self, rxn_set: ReactionSet = None, phases: SolidPhaseSet = None):
        self.lib = {}
        if phases is None:
            self.phases = SolidPhaseSet.from_rxn_set(rxn_set)
        else:
            self.phases = phases

    def add_rxns_at_temp(self, rxns: ScoredReactionSet, temp: int) -> int:
        self.lib[int(temp)] = rxns
        return temp
    
    def get_rxns_at_temp(self, temp: int) -> ScoredReactionSet:
        return self.lib[temp]
    
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