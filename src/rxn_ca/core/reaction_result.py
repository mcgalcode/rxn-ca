from pylattica.core import SimulationState, SimulationResult
from ..reactions import ScoredReactionSet
from .heating import HeatingSchedule

class ReactionResult(SimulationResult):
    """A class that stores the result of running a simulation. Keeps track of all
    the steps that the simulation proceeded through, and the set of reactions that
    was used in the simulation.
    """

    @classmethod
    def from_dict(cls, res_dict):
        diffs = res_dict["diffs"]
        res = ReactionResult(
            SimulationState.from_dict(res_dict["initial_state"]),
            ScoredReactionSet.from_dict(res_dict["rxn_set"])
        )        
        for diff in diffs:
            res.add_step(diff)
        return res

    def __init__(self,
                 starting_state: SimulationState,
                 rxn_set: ScoredReactionSet,
                 heating_schedule: HeatingSchedule = None):
        """Initializes a ReactionResult with the reaction set used in the simulation

        Args:
            rxn_set (ScoredReactionSet):
        """
        super().__init__(starting_state)
        self.rxn_set: ScoredReactionSet = rxn_set
        self.heating_schedule = heating_schedule
    
    def as_dict(self):
        return {
            **super().as_dict(),
            "rxn_set": self.rxn_set.as_dict(),
        }