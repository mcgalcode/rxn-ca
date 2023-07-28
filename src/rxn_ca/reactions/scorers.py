import math
from tqdm import tqdm
from abc import ABC, abstractmethod
from .scored_reaction import ScoredReaction
from .scored_reaction_set import ScoredReactionSet
from .reaction_library import ReactionLibrary

from ..core.solid_phase_set import SolidPhaseSet
from .utils import get_phase_vols
from typing import List

from rxn_network.reactions.reaction_set import ReactionSet
from rxn_network.reactions.computed import ComputedReaction

def arrhenius_score(energy, temp):
    return math.exp(-energy / (8.6e-5 * temp * 5))


class BasicScore(ABC):

    @abstractmethod
    def score(self, rxn: ComputedReaction):
        pass

class ArrheniusScore(BasicScore):

    def __init__(self, temp):
        self.temp = temp

    def score(self, rxn):
        return math.exp(-rxn.energy_per_atom / (8.6e-5 * self.temp * 8))

class ConstantScore(BasicScore):

    def __init__(self, score):
        self.score = score

    def score(self, _):
        return self.score


def score_rxns(reactions: list[ComputedReaction], scorer: BasicScore):
    rxn_list = list(reactions.get_rxns())
    all_phases = list(set([ c.reduced_formula for r in rxn_list for c in r.compositions]))
    vols = get_phase_vols(all_phases)

    phase_set = SolidPhaseSet(all_phases, volumes=vols)

    scored_reactions = []

    for rxn in tqdm(reactions, desc="Scoring reactions..."):
        scored_rxn = ScoredReaction.from_rxn_network(scorer.score(rxn), rxn, vols)
        scored_reactions.append(scored_rxn)

    return scored_reactions, phase_set

def score_rxns_many_temps(reactions: ReactionSet, temps: List[int]):
    scorers = [ArrheniusScore(t) for t in temps]
    all_phases = list(set([e.composition.reduced_formula for e in reactions.entries]))
    vols = get_phase_vols(all_phases)
    phase_set = SolidPhaseSet(all_phases, volumes=vols)

    rsets: List[ReactionSet] = []
    for t in tqdm(temps, desc="Calculating reaction energies at temperatures..."):
        rsets.append(reactions.set_new_temperature(t))

    rxn_library = ReactionLibrary(reactions, phase_set)
    for t, scorer, rset in zip(temps, scorers, rsets):
        scored_rxns: List[ScoredReaction] = []
        for rxn in tqdm(rset.get_rxns(), desc=f'Calculating reaction scores at {t}'):
            scored_rxn = ScoredReaction.from_rxn_network(scorer.score(rxn), rxn, vols)
            scored_rxns.append(scored_rxn)

        rxn_set = ScoredReactionSet(scored_rxns, phase_set)
        rxn_library.add_rxns_at_temp(rxn_set, t)
    
    return rxn_library


def score_rxn_network_rxn_set(reactions: ReactionSet, scorer):
    scores = [scorer.score(rxn) for rxn in reactions]
    rxn_scores = zip(reactions, scores)
    scored_reactions = []
    for rxn_score in rxn_scores:
        rxn = rxn_score[0]
        score = rxn_score[1]
        scored_rxn = ScoredReaction.from_rxn_network(score, rxn)
        scored_reactions = scored_reactions + [scored_rxn]

    return scored_reactions