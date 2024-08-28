from typing import Dict
import random

from ..analysis.reaction_step_analyzer import ReactionStepAnalyzer
from pylattica.core import SimulationState
from pylattica.core import BasicController
from pylattica.core.neighborhood_builders import NeighborhoodBuilder
from pylattica.core.periodic_structure import PeriodicStructure
from pylattica.core.simulation_state import SimulationState
from pylattica.discrete import PhaseSet
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY, VACANT
from pylattica.structures.square_grid.neighborhoods import MooreNbHoodBuilder

class PhaseGrowthController(BasicController):

    def __init__(
        self,
        phase_set: PhaseSet,
        periodic_struct: PeriodicStructure,
        desired_phase_vols: Dict,
        background_phase: str = VACANT,
        nb_builder: NeighborhoodBuilder = None,
    ) -> None:
        
        self.background_phase = background_phase
        self.phase_set: PhaseSet = phase_set
        self.analyzer = ReactionStepAnalyzer(self.phase_set)
        self.desired_phase_vols_abs = desired_phase_vols

        self.known_empty_ids = periodic_struct.site_ids.copy()

        if nb_builder is None:
            self.nb_builder = MooreNbHoodBuilder(1, dim=periodic_struct.dim)
        else:
            self.nb_builder = nb_builder

        self.nb_graph = self.nb_builder.get(periodic_struct)
    
    def get_random_site(self, _):
        if len(self.known_empty_ids) > 0:
            return random.choice(self.known_empty_ids)
        else:
            return 1

    def get_state_update(self, site_id: int, prev_state: SimulationState):
        
        curr_state = prev_state.get_site_state(site_id)
        if curr_state[DISCRETE_OCCUPANCY] == self.background_phase:
            nb_phases = set()
            for nb_id in self.nb_graph.neighbors_of(site_id):
                nb_phase = prev_state.get_site_state(nb_id)[DISCRETE_OCCUPANCY]
                if nb_phase != self.background_phase:
                    nb_phases.add(nb_phase)

            if len(nb_phases) > 0:
                chosen_phase = None
                volume_breakdown = self.analyzer.set_step_group(prev_state).get_all_absolute_phase_volumes()

                deficient_candidates = []

                for phase in nb_phases:
                    diff = volume_breakdown[phase] - self.desired_phase_vols_abs[phase]
                    deficient_candidates.append((phase, diff))
                
                most_deficient = min(deficient_candidates, key = lambda c: c[1])[0]

                if len(deficient_candidates) > 0:
                    chosen_phase = most_deficient
                    self.known_empty_ids.remove(site_id)
                    return {DISCRETE_OCCUPANCY: chosen_phase}
                else:
                    return {}
            else:
                return {}
        else:
            if site_id in self.known_empty_ids:
                self.known_empty_ids.remove(site_id)
            return {}