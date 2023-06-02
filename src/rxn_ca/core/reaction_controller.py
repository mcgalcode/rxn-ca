import math
import numpy as np
import random
from typing import Dict, List, Tuple
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY
from pylattica.core.periodic_structure import PeriodicStructure
from pylattica.core.simulation_state import GENERAL, SITES, SimulationState
from pylattica.square_grid.neighborhoods import MooreNbHoodBuilder, VonNeumannNbHood2DBuilder, VonNeumannNbHood3DBuilder
from pylattica.core.basic_controller import BasicController

from .normalizers import normalize
from .solid_phase_set import SolidPhaseSet
from .reaction_result import ReactionResult
from .reaction_setup import VOLUME
from ..reactions import ScoredReactionSet, ScoredReaction

DEFAULT_GASES = [
    "H2",
    "O2",
    "H2O",
    "CO",
    "CO2",
    "N2"
]

class ReactionController(BasicController):

    @classmethod
    def get_neighborhood_from_size(cls, size, nb_builder = VonNeumannNbHood2DBuilder):
        neighborhood_radius = cls.nb_radius_from_size(size)
        print(f'Using neighborhood of size {neighborhood_radius}')
        return nb_builder(neighborhood_radius)

    @classmethod
    def get_neighborhood_from_structure(cls, structure: PeriodicStructure, nb_builder = None):
        if nb_builder is None:
            if structure.dim == 3:
                nb_builder = VonNeumannNbHood3DBuilder
            else:
                nb_builder = VonNeumannNbHood2DBuilder
        return cls.get_neighborhood_from_size(structure.bounds[0], nb_builder=nb_builder)

    @classmethod
    def nb_radius_from_size(cls, size: int) -> int:
        """Provided the side length of the simulation stage, generate the
        appropriate filter size for the simulation. This routine chooses the smaller of:
            1. The filter size that will cover the entire reaction stage
            2. 25, which is a performance compromise
        If the filter size is 25, then we are looking at reactant pairs up to 12 squares
        away from eachother. Given that probability of proceeding scales as 1/d^2, this is equivalent
        to cutting off possibility after the scaling drops below 1/144.

        Args:
            size (int): The side length of the simulation stage.

        Returns:
            int: The side length of the filter used in the convolution step.
        """
        return math.floor(min((size - 1) * 2 + 1, 21) / 2)

    def __init__(self,
        structure: PeriodicStructure,
        scored_rxns: ScoredReactionSet = None,
        inertia = 1,
        open_species = {},
        free_species = None,
        temperature = None,
    ) -> None:
        
        if free_species is None:
            self.free_species = DEFAULT_GASES
        else:
            self.free_species = free_species

        print(f'Interpreting {self.free_species} as gaseous in this reaction')

        if scored_rxns is not None:
            self.rxn_set = scored_rxns
        
        self.structure = structure
        self.temperature = temperature
        self.phase_set: SolidPhaseSet = scored_rxns.phases

        nb_hood_builder = ReactionController.get_neighborhood_from_structure(structure)

        self.nb_graph = nb_hood_builder.get(structure)
        self.nucleation_nb_graph = MooreNbHoodBuilder(dim = structure.dim).get(structure)
        self.inertia = inertia

        # proxy for partial pressures
        self.effective_open_distances = {}
        for specie, strength in open_species.items():
            self.effective_open_distances[specie] = strength
    
    def get_random_site(self):
        return random.randint(0,len(self.structure.site_ids) - 1)

    def instantiate_result(self, starting_state: SimulationState):
        return ReactionResult(starting_state, self.rxn_set)

    def get_state_update(self, site_id: int, prev_state: SimulationState):
        np.random.seed(None)

        if random.random() < self.inertia:
            return {}

        center_site_state = prev_state.get_site_state(site_id)

        center_species = center_site_state[DISCRETE_OCCUPANCY]

        if center_species == self.phase_set.FREE_SPACE:
            return {}
        else:
            possible_reactions = self.rxns_at_site(site_id, prev_state)
            chosen_rxn = self.choose_reaction(possible_reactions)
            rxn: ScoredReaction = chosen_rxn["reaction"]
            if rxn.is_identity:
                return {}

            site_updates = {
                site_id: self.get_cell_updates(center_site_state, rxn),
            }
            
            if not chosen_rxn['open_el']:
                other_state = chosen_rxn['other_site_state']
                other_site_id = other_state["_site_id"]
                site_updates[other_site_id] = self.get_cell_updates(other_state, rxn)

            return site_updates


    def immediate_neighbors(self, site_id: int, state: SimulationState):
        neighbor_phases = []
        for nb_id in self.nucleation_nb_graph.neighbors_of(site_id):
            nb_state = state.get_site_state(nb_id)
            neighbor_phases.append(nb_state[DISCRETE_OCCUPANCY])

        return neighbor_phases

    def get_rxns_from_step(self, simulation_state, coords):
        site_id = self.structure.site_at(coords)['_site_id']
        return self.rxns_at_site(site_id, simulation_state)

    def rxns_at_site(self, site_id: int, state: SimulationState):
        curr_state = state.get_site_state(site_id)
        this_phase = curr_state[DISCRETE_OCCUPANCY]

        neighbor_phases = self.immediate_neighbors(site_id, state)

        # Look through neighborhood, enumerate possible reactions
        rxns = []

        for nb_id, distance in self.nb_graph.neighbors_of(site_id, include_weights=True):

            site_state = state.get_site_state(nb_id)
            neighbor_phase = site_state[DISCRETE_OCCUPANCY]
            rxn, score = self.get_rxn_and_score([neighbor_phase, this_phase], distance, neighbor_phases, this_phase)
            rxns.append({
                'other_site_state': site_state,
                'reaction': rxn,
                'score': score,
                'other_phase': neighbor_phase,
                'open_el': False
            })

        # Readd open species reactions later
        possible_reactions = [*rxns, *self.rxns_with_open_species(this_phase, neighbor_phases)]

        return possible_reactions

    def rxns_with_open_species(self, this_phase, neighbor_phases):
        rxns = []
        for specie, dist in self.effective_open_distances.items():
            rxn, score = self.get_rxn_and_score([this_phase, specie], dist, neighbor_phases, this_phase)
            rxns.append({
                'reaction': rxn,
                'score': score,
                'other_site_state': None,
                'other_phase': None,
                'open_el': True
            })

        return rxns

    def choose_reaction(self, rxns_and_scores: List[Dict]) -> Tuple[ScoredReaction, int, str]:
        scores: list[float] = [
            rxn['score'] for rxn in rxns_and_scores
        ]

        scores: np.array = np.array(scores)
        normalized: np.array = normalize(scores)
        choices: np.array = np.array(range(0,len(rxns_and_scores)))

        chosen_idx = np.random.choice(choices, p=normalized)

        return rxns_and_scores[chosen_idx]

    def get_cell_updates(self, cell_state, reaction: ScoredReaction):
        updates = {}
        
        current_spec = cell_state[DISCRETE_OCCUPANCY]
        current_vol = cell_state[VOLUME]

        # The current site should be replaced if a randomly chosen reactant is the current species
        current_phase_replacement = self.get_phase_replacement_from_reaction(reaction, current_spec, current_vol)
        if current_phase_replacement is not None:
            volume_ratio = reaction.product_reactant_stoich_ratio
            updates[DISCRETE_OCCUPANCY] = current_phase_replacement
            updates[VOLUME] = volume_ratio * cell_state[VOLUME]
        
        return updates

    
    def get_phase_replacement_from_reaction(self, rxn: ScoredReaction, reactant_phase: str, reactant_vol: float) -> Dict:
        stoich_fraction = rxn.reactant_stoich_fraction(reactant_phase)
        # IMPORTANT: This division ensures that the likelihood of consuming a particular cell decreases with
        # the size of that cell - it should take twice as many "tries" to consume twice as much
        # volume
        adjusted = stoich_fraction / reactant_vol
        if random.random() < adjusted:
            prod_sf = np.array([rxn.product_stoich(p) for p in rxn.products])
            likelihoods: np.array = prod_sf / prod_sf.sum()
            new_phase_name = str(np.random.choice(list(rxn.products), p=likelihoods))

            if new_phase_name in self.free_species:
                new_phase_name = self.phase_set.FREE_SPACE
            
            return new_phase_name
        else:   
            return None

    def get_rxn_and_score(self, reactants, distance, neighbor_phases, replaced_phase):

        possible_reaction = self.rxn_set.get_reaction(reactants)
        if possible_reaction is not None:
            score = self.get_score_contribution(possible_reaction.competitiveness, distance)
            if len(possible_reaction.reactants) > 1:
                score = self.adjust_score_for_nucleation(score, neighbor_phases, possible_reaction.products, possible_reaction.reactants, replaced_phase)
            return (possible_reaction, score)
        else:
            # Utilize the self reaction
            rxn = self.rxn_set.get_reaction([replaced_phase])
            score = self.get_score_contribution(self.inertia, distance)
            return (rxn, score)


    def adjust_score_for_nucleation(self, score, neighbors, products, reactants, current_phase):
        return score
        # new_score = score * 1
        # for neighbor in neighbors:
        #     if neighbor in products:
        #         new_score = new_score * 3
        
        # return new_score

    # def adjust_score_for_nucleation(self, score, neighbors, products, reactants):
    #     # penalize a new phase not neighboring any products
    #     if not any([n in products for n in neighbors]):
    #         return score * 0.15
    #     else:
    #         return score

    def get_score_contribution(self, weight, distance):
        return weight * 1 / distance ** 3
