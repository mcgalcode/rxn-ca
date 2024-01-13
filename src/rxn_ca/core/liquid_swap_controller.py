import math
import numpy as np
import random
from typing import Dict, List, Tuple
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY
from pylattica.core.periodic_structure import PeriodicStructure
from pylattica.core.constants import GENERAL, SITE_ID, SITES
from pylattica.core.simulation_state import SimulationState
from pylattica.structures.square_grid.neighborhoods import VonNeumannNbHood2DBuilder, VonNeumannNbHood3DBuilder
from pylattica.core.basic_controller import BasicController

from ..phases.solid_phase_set import SolidPhaseSet
from .reaction_result import ReactionResult
from .constants import VOLUME, GASES_EVOLVED, GASES_CONSUMED
from ..reactions import ScoredReactionSet
from .reaction_calculator import ReactionCalculator

class LiquidSwapController(BasicController):

    @classmethod
    def get_neighborhood(cls, nb_builder = VonNeumannNbHood2DBuilder):
        neighborhood_radius = 1
        print(f'Using neighborhood of size {neighborhood_radius}')
        return nb_builder(neighborhood_radius)

    @classmethod
    def get_neighborhood_from_structure(cls, structure: PeriodicStructure, nb_builder = None):
        if nb_builder is None:
            if structure.dim == 3:
                nb_builder = VonNeumannNbHood3DBuilder
            else:
                nb_builder = VonNeumannNbHood2DBuilder
        return cls.get_neighborhood(nb_builder=nb_builder)

    def __init__(self,
        structure: PeriodicStructure,
        scored_rxns: ScoredReactionSet = None,
        inertia = 1,
        open_species = {}
    ) -> None:
        self.rxn_set = scored_rxns
        self.structure = structure
        nb_hood_builder = LiquidSwapController.get_neighborhood_from_structure(structure)
        self.nb_graph = nb_hood_builder.get(structure)
        self.reaction_calculator = ReactionCalculator(self.nb_graph, scored_rxns, inertia, open_species)
        # Defines the atmosphere
    def set_rxn_set(self, rxn_set: ScoredReactionSet):
        self.rxn_set = rxn_set
    
    def set_temperature(self, temp: int):
        self.temperature = temp

    def instantiate_result(self, starting_state: SimulationState):
        return ReactionResult(starting_state)

    def get_state_update(self, site_id: int, prev_state: SimulationState):
        np.random.seed(None)

        site_state = prev_state.get_site_state(site_id)
        species = site_state[DISCRETE_OCCUPANCY]
        updates = {}

        if species == SolidPhaseSet.FREE_SPACE:
            return updates
        
        if self.rxn_set.phases.is_melted(species, self.temperature):
            pass
            # melt and swap
        else:
            updates = self.reaction_calculator.get_state_update(site_id, prev_state)

        return updates
