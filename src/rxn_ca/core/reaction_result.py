from pylattica.core import SimulationState, SimulationResult
from pylattica.core.constants import SITES, GENERAL

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
        )
        for diff in diffs:
            if SITES in diff:
                diff[SITES] = { int(k): v for k, v in diff[SITES].items() }
            if GENERAL not in diff and SITES not in diff:
                diff = { int(k): v for k, v in diff.items() }
            res.add_step(diff)
        return res

    def __init__(self,
                 starting_state: SimulationState):
        """Initializes a ReactionResult with the reaction set used in the simulation

        Args:
            rxn_set (ScoredReactionSet):
        """
        super().__init__(starting_state)
    