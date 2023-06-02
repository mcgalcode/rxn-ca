from ..reactions import ScoredReactionSet
from pylattica.core import Simulation

class ReactionSimulation():


    def __init__(self, reactions: ScoredReactionSet, simulation: Simulation):
        self.reactions = reactions
        self.state = simulation.state
        self.structure = simulation.structure
