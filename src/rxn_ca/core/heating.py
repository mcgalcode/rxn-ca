from typing import List

from rxn_network.reactions.reaction_set import ReactionSet

from ..reactions.reaction_library import ReactionLibrary
from ..reactions.scorers import score_rxns_many_temps

class HeatingStep():
    """Captures the information about a single step during a heating schedule.
    """

    def __init__(self, duration, temp):
        self.duration = duration
        self.temp = temp

class HeatingSchedule():
    """Captures the information of a heating schedule, e.g. ramping up
    to a particular temperature, holding, and then cooling back down
    """

    @classmethod
    def from_const_temp(cls, temp, duration):
        return cls([(duration, temp)])

    def __init__(self, schedule):
        # schedules at first can just be a series of steps, e.g.:
        # 2000 steps, 300k
        # 5000 steps, 900k
        # 2000 steps 500k
        self.steps: List[HeatingStep] = []

        for step in schedule:
            heating_step = HeatingStep(*step)
            self.steps.append(heating_step)

    def calculate_rxns(self, rxns: ReactionSet) -> ReactionLibrary:
        """Returns a ReactionLibrary with reactions calculated at every
        temperature in this heating schedule.

        Args:
            rxns (ReactionSet): The ReactionSet containing the reactions to 
            score at each temperature.

        Returns:
            ReactionLibrary:
        """
        temps = [s.temp for s in self.steps]
        return score_rxns_many_temps(rxns, temps)