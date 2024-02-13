from pylattica.core.periodic_structure import PeriodicStructure
from pylattica.core.simulation_state import SimulationState
from pylattica.structures.square_grid.neighborhoods import VonNeumannNbHood2DBuilder, VonNeumannNbHood3DBuilder
from pylattica.core.basic_controller import BasicController

from .reaction_result import ReactionResult
from ..reactions import ScoredReactionSet
from .reaction_calculator import ReactionCalculator

NB_HOOD_RADIUS = 5

class ReactionController(BasicController):

    @classmethod
    def get_neighborhood_builder(cls, nb_builder = VonNeumannNbHood2DBuilder):
        neighborhood_radius = NB_HOOD_RADIUS
        print(f'Using neighborhood of size {neighborhood_radius}')
        return nb_builder(neighborhood_radius)

    @classmethod
    def get_neighborhood_builder_from_structure(cls, structure: PeriodicStructure, nb_builder = None):
        if nb_builder is None:
            if structure.dim == 3:
                nb_builder = VonNeumannNbHood3DBuilder
            else:
                nb_builder = VonNeumannNbHood2DBuilder
        return cls.get_neighborhood_builder(nb_builder=nb_builder)

    @classmethod
    def get_neighborhood_from_structure(cls, structure: PeriodicStructure):
        return cls.get_neighborhood_builder_from_structure(structure).get(structure)

    def __init__(self,
        structure: PeriodicStructure,
        rxn_calculator: ReactionCalculator,
    ) -> None:
        self.reaction_calculator = rxn_calculator
        self.structure = structure

    def set_rxn_set(self, rxn_set: ScoredReactionSet):
        self.reaction_calculator.set_rxn_set(rxn_set)

    def instantiate_result(self, starting_state: SimulationState):
        return ReactionResult(starting_state)

    def get_state_update(self, site_id: int, prev_state: SimulationState):
        return self.reaction_calculator.get_state_update(site_id, prev_state)