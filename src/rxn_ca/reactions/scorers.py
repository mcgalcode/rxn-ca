import math
from .scored_reaction import ScoredReaction

from rxn_network.reactions.reaction_set import ReactionSet

def arrhenius_score(energy, temp):
    return math.exp(-energy / (8.6e-5 * temp * 5))

class ArrheniusScore():

    def __init__(self, temp):
        self.temp = temp

    def score(self, rxn):
        return math.exp(-rxn.energy_per_atom / (8.6e-5 * self.temp * 8))

class ConstantScore():

    def __init__(self, score):
        self.score = score

    def score(self, _):
        return self.score


def score_rxns(reactions: list[ScoredReaction], scorer):
    scores = [scorer.score(rxn) for rxn in reactions]
    rxn_scores = zip(reactions, scores)
    scored_reactions = []
    for rxn_score in rxn_scores:
        rxn = rxn_score[0]
        score = rxn_score[1]
        scored_rxn = ScoredReaction.from_rxn_network(score, rxn)
        scored_reactions = scored_reactions + [scored_rxn]

    return scored_reactions

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