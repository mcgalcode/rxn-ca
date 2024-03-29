import numpy as np
import random
import math
from typing import Dict, List, Tuple
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY
from pylattica.core.periodic_structure import PeriodicStructure
from pylattica.core.constants import GENERAL, SITE_ID, SITES
from pylattica.core.simulation_state import SimulationState
from pylattica.core.neighborhoods import Neighborhood
from pylattica.structures.square_grid.neighborhoods import VonNeumannNbHood3DBuilder
from pylattica.core.basic_controller import BasicController

from ..phases.solid_phase_set import SolidPhaseSet
from .reaction_result import ReactionResult
from .constants import VOLUME, REACTION_CHOSEN
from ..reactions import ScoredReactionSet
from .reaction_calculator import ReactionCalculator

def swap_chance(tm_frac):
    num = (20*tm_frac - 18.5)
    den = 2 * math.sqrt((20*tm_frac - 18.5) ** 2 + 1)
    return num / den + 1 / 2

class LiquidSwapController(BasicController):

    @classmethod
    def get_neighborhood_from_structure(cls, structure: PeriodicStructure):
        return VonNeumannNbHood3DBuilder(1).get(structure)
    
    def __init__(self,
        structure: PeriodicStructure,
        rxn_calculator: ReactionCalculator,
    ) -> None:
        self.structure = structure
        self.reaction_calculator = rxn_calculator
        self.temperature = None

    def set_rxn_set(self, rxn_set: ScoredReactionSet):
        self.reaction_calculator.set_rxn_set(rxn_set)
    
    def set_temperature(self, temp: int):
        self.temperature = temp

    def instantiate_result(self, starting_state: SimulationState):
        return ReactionResult(starting_state)

    def get_state_update(self, site_id: int, prev_state: SimulationState):
        np.random.seed(None)

        site_state = prev_state.get_site_state(site_id)
        species = site_state[DISCRETE_OCCUPANCY]
        updates = {}
        updates[GENERAL] = {
            REACTION_CHOSEN: None
        }

        if species == SolidPhaseSet.FREE_SPACE:
            return updates
        
        diff = self.temperature / self.reaction_calculator.rxn_set.phases.get_melting_point(species)

        if species == SolidPhaseSet.FREE_SPACE or random.random() < swap_chance(diff):
            nb_ids = self.reaction_calculator.neighborhood_graph.neighbors_of(site_id)

            other_id = random.choice(nb_ids)
            other_state = prev_state.get_site_state(other_id)

            updates[SITES] = {
                site_id: {
                    DISCRETE_OCCUPANCY: other_state[DISCRETE_OCCUPANCY],
                    VOLUME: other_state[VOLUME]
                },
                other_id: {
                    DISCRETE_OCCUPANCY: site_state[DISCRETE_OCCUPANCY],
                    VOLUME: site_state[VOLUME]                    
                }
            }
        else:
            updates = self.reaction_calculator.get_state_update(site_id, prev_state)

        return updates
