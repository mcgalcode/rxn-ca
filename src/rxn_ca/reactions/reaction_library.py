from .scored_reaction_set import ScoredReactionSet
from ..core.solid_phase_set import SolidPhaseSet
from rxn_network.reactions.reaction_set import ReactionSet


class ReactionLibrary():
    """Contains a mapping of temperatures to ScoredReactionSet objects
    which contain reactions scored at the given temperature. Used in multi-stage
    reaction simulations where different stages in the reaction are run at
    different temperatures.
    """

    def __init__(self, rxn_set: ReactionSet, phases: SolidPhaseSet):
        self._lib = {}
        self.original_rxns = rxn_set
        self.phases = phases

    def add_rxns_at_temp(self, rxns: ScoredReactionSet, temp: int) -> int:
        self._lib[temp] = rxns
        return temp
    
    def get_rxns_at_temp(self, temp: int) -> ScoredReactionSet:
        return self._lib[temp]