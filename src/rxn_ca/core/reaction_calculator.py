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

from .normalizers import normalize
from ..phases.solid_phase_set import SolidPhaseSet
from .reaction_result import ReactionResult
from .constants import VOLUME, GASES_EVOLVED, GASES_CONSUMED
from ..reactions import ScoredReactionSet, ScoredReaction

from dataclasses import dataclass, field
from copy import copy

def choose_from_list(choices, scores):
    scores: np.array = np.array(scores)
    normalized: np.array = normalize(scores)
    idxs: np.array = np.array(range(0,len(choices)))

    chosen_idx = np.random.choice(idxs, p=normalized)

    return choices[chosen_idx]    


def scale_score_by_distance(score, distance):
    return score * 1 / distance ** 2

NB_HOOD_RADIUS = 5

@dataclass
class SiteInteraction:

    score: float
    site_states: List[Dict] = field(default_factory=list)
    reactions: List[ScoredReaction] = field(default_factory=list)
    atmosphere_reactant: str = None
    is_no_op: bool = False

class ReactionCalculator():

    def __init__(self,
        neighborhood_graph,
        scored_rxns: ScoredReactionSet = None,
        inertia = 1,
        open_species = {},
    ) -> None:
        self.rxn_set = scored_rxns
        self.inertia = inertia
        self.neighborhood_graph = neighborhood_graph
        self.effective_open_distances = copy(open_species)

    def set_rxn_set(self, rxn_set: ScoredReactionSet):
        self.rxn_set = rxn_set

    def get_state_update(self, site_id: int, prev_state: SimulationState):
        updates = {}

        # Get the set of possible interactions - cell-cell reactions,cell-gas reactions and no-ops
        possible_interactions = self.possible_interactions_at_site(site_id, prev_state)
        selected_interaction = self.choose_interaction(possible_interactions)

        if selected_interaction.is_no_op:
            return updates
        
        updates[GENERAL] = {}
        updates[SITES] = {}

        # Select a reaction - recall the convex reaction hull: there are often
        # many possible reactions between two precursors
        rxns: List[ScoredReaction] = selected_interaction.reactions
        selected_reaction: ScoredReaction = choose_from_list(rxns, [rxn.competitiveness for rxn in rxns])

        # Proceed this reaction at all relevant site states
        for site_state in selected_interaction.site_states:
            site_species = site_state[DISCRETE_OCCUPANCY]
            site_vol     = site_state[VOLUME]
            site_id      = site_state[SITE_ID]

            if not self.should_reaction_proceed(selected_reaction, site_species, site_vol):
                continue

            product_phase  = self.get_product_from_reaction(selected_reaction)
            product_volume = selected_reaction.convert_reactant_amt_to_product_amt(site_species, site_vol, product_phase)

            # If it's a gaseous product, do some accounting to maintain mass balance
            # then replace the phase with empty space
            if self.rxn_set.phases.is_gas(product_phase):
                gas_amts = prev_state.get_general_state().get(GASES_EVOLVED)
                if product_phase in gas_amts:
                    gas_amts[product_phase] = gas_amts[product_phase] + product_volume
                else:
                    gas_amts[product_phase] = product_volume

                updates[GENERAL][GASES_EVOLVED] = gas_amts
                # updates[SITES][site_id] = {
                #     DISCRETE_OCCUPANCY: SolidPhaseSet.FREE_SPACE,
                # }
            
            # Otherwise, if there is a selected product it must be solid, so just replace
            # the old phase with the new one and the updated volume
            elif product_phase is not None:
                updates[SITES][site_id] = {
                    DISCRETE_OCCUPANCY: product_phase,
                    VOLUME: product_volume
                }
            
        return updates

    def possible_interactions_at_site(self, site_one_id: int, state: SimulationState):
        site_one_state = state.get_site_state(site_one_id)
        site_one_phase = site_one_state[DISCRETE_OCCUPANCY]

        # Look through neighborhood, enumerate possible reactions
        possible_interactions = []

        for nb_id, distance in self.neighborhood_graph.neighbors_of(site_one_id, include_weights=True):
            site_two_state = state.get_site_state(nb_id)
            site_two_phase = site_two_state[DISCRETE_OCCUPANCY]
            possible_reactions = self.rxn_set.get_reactions([site_two_phase, site_one_phase])

            if len(possible_reactions) > 0 and not possible_reactions[0].is_identity:
                interaction_score = self.adjust_score_for_distance(possible_reactions[0].competitiveness, distance)
                interaction = SiteInteraction(
                    site_states=[site_one_state, site_two_state],
                    reactions=possible_reactions,
                    atmosphere_reactant=None,
                    score=interaction_score
                )
            else:
                interaction = SiteInteraction(
                    is_no_op=True,
                    score=self.inertia
                )
            
            possible_interactions.append(interaction)

        # Add interactions with atmosphere
        possible_interactions = [*possible_interactions, *self.interactions_with_open_species(site_one_state)]

        return possible_interactions

    def interactions_with_open_species(self, site_state: Dict):
        site_phase = site_state[DISCRETE_OCCUPANCY]
        interactions = []

        for specie, dist in self.effective_open_distances.items():
            rxns = self.rxn_set.get_reactions([site_phase, specie])
            if len(rxns) > 0:
                interaction_score = self.adjust_score_for_distance(rxns[0].competitiveness, dist)
                interactions.append(SiteInteraction(
                    site_states=[site_state],
                    reactions=rxns,
                    atmosphere_reactant=specie,
                    score=interaction_score
                ))

        return interactions

    def choose_interaction(self, interactions: List[SiteInteraction]) -> SiteInteraction:
        scores: list[float] = [
            interaction.score for interaction in interactions
        ]
        return choose_from_list(interactions, scores)
    
    def should_reaction_proceed(self, rxn: ScoredReaction, reactant_phase: str, reactant_vol: float) -> Dict:
        stoich_fraction = rxn.solid_reactant_stoich_fraction(reactant_phase)
        # IMPORTANT: This division ensures that the likelihood of consuming a particular cell decreases with
        # the size of that cell - it should take twice as many "tries" to consume twice as much
        # volume
        adjusted = stoich_fraction / reactant_vol
        return random.random() < adjusted
    
    def get_product_from_reaction(self, rxn: ScoredReaction) -> Dict:
        products = list(rxn.products)
        product_stoich_coeffs = np.array([rxn.product_stoich(p) for p in products])
        likelihoods: np.array = product_stoich_coeffs / product_stoich_coeffs.sum()
        new_phase_name = str(np.random.choice(products, p=likelihoods))
        return new_phase_name

    def adjust_score_for_distance(self, score, distance):
        return score * 1 / distance ** 3
