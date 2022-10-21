import math
import numpy as np
import random
from typing import Dict, List, Tuple
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY
from pylattica.core.periodic_structure import PeriodicStructure
from pylattica.core.simulation_state import GENERAL, SITES, SimulationState
from pylattica.square_grid.neighborhoods import MooreNbHoodBuilder, VonNeumannNbHood2DBuilder, VonNeumannNbHood3DBuilder
from pylattica.core.basic_controller import BasicController

from .solid_phase_set import SolidPhaseSet
from .reaction_result import ReactionResult
from ..reactions import ScoredReactionSet, ScoredReaction

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
        phase_set: SolidPhaseSet,
        structure: PeriodicStructure,
        scored_rxns: ScoredReactionSet = None,
        inertia = 1,
        open_species = {},
        free_species = [],
        temperature = None,
    ) -> None:

        if scored_rxns is not None:
            self.rxn_set = scored_rxns
        
        self.structure = structure
        self.temperature = temperature
        self.phase_set: SolidPhaseSet = phase_set

        nb_hood_builder = ReactionController.get_neighborhood_from_structure(structure)

        self.nb_graph = nb_hood_builder.get(structure)
        self.nucleation_nb_graph = MooreNbHoodBuilder(1, dim=structure.dim).get(structure)
        self.inertia = inertia
        self.free_species = free_species
        # proxy for partial pressures
        self.effective_open_distances = {}
        for specie, strength in open_species.items():
            self.effective_open_distances[specie] = strength
    
    def get_random_site(self):
        return random.randint(0,len(self.structure.site_ids) - 1)

    def instantiate_result(self, starting_state: SimulationState):
        return ReactionResult(starting_state, self.rxn_set, self.phase_set)

    def get_state_update(self, site_id: int, prev_state: SimulationState):
        np.random.seed(None)

        curr_state = prev_state.get_site_state(site_id)

        curr_species = curr_state[DISCRETE_OCCUPANCY]
        if curr_species == self.phase_set.FREE_SPACE:
            return {}
        else:
            possible_reactions = self.rxns_at_site(site_id, prev_state)
            chosen_rxn, other_site_id, other_phase = self.choose_reaction(possible_reactions)

            site_updates = self.get_updates_from_reaction(
                chosen_rxn,
                site_id,
                other_site_id,
                curr_species,
                other_phase
            )

            return {
                GENERAL: {
                    'rxn': str(chosen_rxn)
                },
                SITES: site_updates
            }


    def immediate_neighbors(self, site_id: int, state: SimulationState):
        neighbor_phases = []
        for nb_id in self.nucleation_nb_graph.neighbors_of(site_id):
            nb_state = state.get_site_state(nb_id)
            neighbor_phases.append(nb_state[DISCRETE_OCCUPANCY])

        return neighbor_phases

    def get_rxns_from_step(self, simulation_state, coords):
        site_id = self.structure.site_at(coords)['id']
        return self.rxns_at_site(site_id, simulation_state)

    def rxns_at_site(self, site_id: int, state: SimulationState):
        curr_state = state.get_site_state(site_id)
        this_phase = curr_state[DISCRETE_OCCUPANCY]

        neighbor_phases = self.immediate_neighbors(site_id, state)

        # Look through neighborhood, enumerate possible reactions
        rxns = []

        for nb_id, distance in self.nb_graph.neighbors_of(site_id, include_weights=True):

            neighbor_phase = state.get_site_state(nb_id)[DISCRETE_OCCUPANCY]
            rxn, score = self.get_rxn_and_score([neighbor_phase, this_phase], distance, neighbor_phases, this_phase)
            rxns.append({
                'site_id': nb_id,
                'other_phase': neighbor_phase,
                'reaction': rxn,
                'score': score
            })

        # Readd open species reactions later
        # possible_reactions = {**rxns, **self.rxns_with_open_species(this_phase, neighbor_phases)}

        return rxns

    def rxns_with_open_species(self, this_phase, neighbor_phases):
        rxns = {}
        for specie, dist in self.effective_open_distances.items():
            rxn, score = self.get_rxn_and_score([this_phase, specie], dist, neighbor_phases, this_phase)
            if str(rxn) in rxns:
                rxns[str(rxn)]['score'] += score
            else:
                rxns[str(rxn)] = {
                    'reaction': rxn,
                    'score': score
                }

        return rxns

    def choose_reaction(self, rxns_and_scores: List[Dict]) -> Tuple[ScoredReaction, int, str]:
        rxns: list[ScoredReaction] = [
            rxn['reaction'] for rxn in rxns_and_scores
        ]
        scores: list[float] = [
            rxn['score'] for rxn in rxns_and_scores
        ]
        site_ids: list[float] = [
            rxn['site_id'] for rxn in rxns_and_scores
        ]
        phases: list[float] = [
            rxn['other_phase'] for rxn in rxns_and_scores
        ]

        scores: np.array = np.array(scores)
        normalized: np.array = normalize(scores)
        choices: np.array = np.array(range(0,len(rxns)))

        chosen_idx = np.random.choice(choices, p=normalized)

        # chosen_rxn: ScoredReaction = rxns[np.argmax(scores)]

        return rxns[chosen_idx], site_ids[chosen_idx], phases[chosen_idx]

    def get_updates_from_reaction(self,
                                  rxn: ScoredReaction,
                                  current_site_id,
                                  other_site_id,
                                  current_spec,
                                  other_spec) -> Tuple[int, ScoredReaction]:
        # Assume at this point that the reaction is proceeding. It has
        # X chance of consuming A and Y chance of consuming B

        updates = {
            current_site_id: {},
            other_site_id: {}
        }

        # The current site should be replaced if a randomly chosen reactant is the current species
        current_phase_replacement = self.get_phase_replacement_from_reaction(rxn, current_spec)
        if current_phase_replacement is not None:
            updates[current_site_id][DISCRETE_OCCUPANCY] = current_phase_replacement
        
        other_phase_replacement = self.get_phase_replacement_from_reaction(rxn, other_spec)
        if other_phase_replacement is not None:
            updates[other_site_id][DISCRETE_OCCUPANCY] = other_phase_replacement

        return updates
    
    def get_phase_replacement_from_reaction(self, rxn: ScoredReaction, phase: str) -> Dict:
        chosen_reactant: str = self.choose_reactant_by_volume(rxn)        
        if chosen_reactant == phase:
            new_phase_name: str = self.choose_product_by_volume(rxn)

            if new_phase_name in self.free_species:
                new_phase_name = self.phase_set.FREE_SPACE
            
            return new_phase_name
        else:
            return None

    def choose_product_by_volume(self, rxn):
        product_stoichs = [rxn.product_stoich(p) for p in rxn.products]
        return self.choose_phase_from_volume_weighted_stoich(rxn.products, product_stoichs)

    def choose_reactant_by_volume(self, rxn):
        reactant_stoichs = [rxn.reactant_stoich(p) for p in rxn.reactants]
        return self.choose_phase_from_volume_weighted_stoich(rxn.reactants, reactant_stoichs)

    def choose_phase_from_volume_weighted_stoich(self, phases, stoichs):
        volumes = self.rxn_set.volumes
        likelihoods: np.array = np.array([stoichs[idx] * volumes.get(phase, 1) for idx, phase in enumerate(phases)])
        likelihoods: np.array = likelihoods / likelihoods.sum()
        return str(np.random.choice(phases, p=np.array(likelihoods)))


    def get_rxn_and_score(self, reactants, distance, neighbor_phases, replaced_phase):

        possible_reaction = self.rxn_set.get_reaction(reactants)

        if possible_reaction is not None:
            score = self.get_score_contribution(possible_reaction.competitiveness, distance)
            score = self.adjust_score_for_nucleation(score, neighbor_phases, possible_reaction.products, possible_reaction.reactants)
            return (possible_reaction, score)
        else:
            # Utilize the self reaction
            rxn = self.rxn_set.get_reaction([replaced_phase])
            score = self.get_score_contribution(self.inertia, distance)
            return (rxn, score)


    # def adjust_score_for_nucleation(self, score, neighbors, products, reactants):
    #     total_friendly_neighbors = len(neighbors)

    #     for neighbor in neighbors:
    #         if neighbor not in products:
    #             score = score * 0.8
    #         if neighbor not in reactants and neighbor not in products:
    #             total_friendly_neighbors -= 1

    #     if total_friendly_neighbors == 0:
    #         return 0
    #     else:
    #         return score

    def adjust_score_for_nucleation(self, score, neighbors, products, reactants):
        # penalize a new phase not neighboring any products
        if not any([n in products for n in neighbors]):
            return score * 0.15
        else:
            return score

    def get_score_contribution(self, weight, distance):
        return weight * 5 / distance
