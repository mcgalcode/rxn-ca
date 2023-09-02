from .scored_reaction_set import ScoredReactionSet
from ..core.solid_phase_set import SolidPhaseSet
from rxn_network.reactions.reaction_set import ReactionSet


from ..core.heating import HeatingSchedule
from ..reactions.scored_reaction import ScoredReaction
from ..reactions.scorers import score_rxns, TammanHuttigScore

from typing import List
from tqdm import tqdm


class ReactionLibrary():
    """Contains a mapping of temperatures to ScoredReactionSet objects
    which contain reactions scored at the given temperature. Used in multi-stage
    reaction simulations where different stages in the reaction are run at
    different temperatures.
    """

    def __init__(self, rxn_set: ReactionSet):
        self._lib = {}
        self.original_rxns = rxn_set
        self.phases = SolidPhaseSet.from_rxn_set(rxn_set)

    def add_rxns_at_temp(self, rxns: ScoredReactionSet, temp: int) -> int:
        self._lib[temp] = rxns
        return temp
    
    def get_rxns_at_temp(self, temp: int) -> ScoredReactionSet:
        return self._lib[temp]
    
    @property
    def temps(self):
        return list(self._lib.keys())
    
    def score_rxns_for_heating_schedule(self, heating_sched: HeatingSchedule, score_class = TammanHuttigScore):
        """Returns a ReactionLibrary with reactions calculated at every
        temperature in this heating schedule.

        Args:
            rxns (ReactionSet): The ReactionSet containing the reactions to 
            score at each temperature.

        Returns:
            ReactionLibrary:
        """
        sched_temps = heating_sched.all_temps

        new_temps = list(set(sched_temps) - set(self.temps))
        if len(new_temps) == 0:
            return
        
        scorers = [score_class(temp=t, phase_set=self.phases) for t in new_temps]

        rsets: List[ReactionSet] = []
        for t in tqdm(new_temps, desc="Calculating reaction energies at temperatures..."):
            rsets.append(self.original_rxns.set_new_temperature(t))

        for t, scorer, rset in zip(new_temps, scorers, rsets):
            scored_rxns: List[ScoredReaction] = score_rxns(rset, scorer, phase_set=self.phases)
            rxn_set = ScoredReactionSet(scored_rxns, self.phases)
            self.add_rxns_at_temp(rxn_set, t)