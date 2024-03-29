from ..reactions import ScoredReactionSet
from .heating import HeatingSchedule
from pylattica.core import Simulation

class ReactionSimulation():


    def __init__(self,
                 reactions: ScoredReactionSet,
                 simulation: Simulation,
                 heating_schedule: HeatingSchedule = None):
        self.state = simulation.state
        self.structure = simulation.structure
