import matplotlib.pyplot as plt
import plotly.graph_objects as go
import numpy as np

from pylattica.core import SimulationState, SimulationResult
from .solid_phase_set import SolidPhaseSet
from ..reactions import ScoredReactionSet



class ReactionResult(SimulationResult):
    """A class that stores the result of running a simulation. Keeps track of all
    the steps that the simulation proceeded through, and the set of reactions that
    was used in the simulation.
    """

    @classmethod
    def from_dict(cls, res_dict):
        res = ReactionResult(
            ScoredReactionSet.from_dict(res_dict["rxn_set"]),
            SolidPhaseSet.from_dict(res_dict["phase_set"])
        )
        for step_dict in res_dict["steps"]:
            res.add_step(SimulationState.from_dict(step_dict))

        return res

    def __init__(self, starting_state: SimulationState, rxn_set: ScoredReactionSet, phase_set: SolidPhaseSet):
        """Initializes a ReactionResult with the reaction set used in the simulation

        Args:
            rxn_set (ScoredReactionSet):
        """
        super().__init__(starting_state)
        self.phase_set = phase_set
        self.rxn_set: ScoredReactionSet = rxn_set

    def as_dict(self):
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "steps": [s.as_dict() for s in self.steps],
            "rxn_set": self.rxn_set.as_dict(),
            "phase_set": self.phase_set.as_dict()
        }