import math
from tqdm import tqdm
from abc import ABC, abstractmethod
from .scored_reaction import ScoredReaction

from ..phases.solid_phase_set import SolidPhaseSet
from typing import List

from rxn_network.reactions.reaction_set import ReactionSet
from rxn_network.reactions.computed import ComputedReaction

from ..phases.gasses import DEFAULT_GASES

def softplus(x):
    return 1/3 * math.log(1 + math.exp(3*x))

def tamman_score_exp(t_tm_ratio):
    return math.exp(4.82*(t_tm_ratio) - 3.21)

def tamman_score_softplus(t_tm_ratio):
    return math.log(1 + math.exp(14 * (t_tm_ratio - 0.63)))

def huttig_score_exp(t_tm_ratio):
    return math.exp(2.41*(t_tm_ratio) - 0.8)

def huttig_score_softplus(t_tm_ratio):
    return math.log(1 + math.exp(8 * (t_tm_ratio - 0.33)))

def erf(x):
    return 0.5 * (1 + math.erf(-35 * (x + 0.03)))

def tamman_erf_score(tm_ratio, delta_g):
    return tamman_score_softplus(tm_ratio) * erf(delta_g)

def huttig_erf_score(tm_ratio, delta_g):
    return huttig_score_softplus(tm_ratio) * erf(delta_g)
         
class BasicScore(ABC):

    def __init__(self, phase_set: SolidPhaseSet, temp: int):
        self.phases = phase_set
        self.temp = temp

    @abstractmethod
    def score(self, rxn: ComputedReaction):
        pass

class TammanHuttigScoreExponential(BasicScore):
    # https://en.wikipedia.org/wiki/Tammann_and_H%C3%BCttig_temperatures


    def score(self, rxn: ComputedReaction):
        phases = [c.reduced_formula for c in rxn.reactants]
        non_gasses = [p for p in phases if p not in DEFAULT_GASES]
        mps = [self.phases.get_melting_point(p) for p in non_gasses]
        min_mp = min(mps)

        # Softplus adjustment
        # delta_g_adjustment = softplus(-rxn.energy_per_atom)
        delta_g_adjustment = softplus(-(2*rxn.energy_per_atom + 0.8))
        

        if len(non_gasses) < len(phases):
            # Huttig
            return huttig_score_exp(self.temp / min_mp) * delta_g_adjustment
        else:
            # Tamman
            return tamman_score_exp(self.temp / min_mp) * delta_g_adjustment

class TammanHuttigScoreSoftplus(BasicScore):
    # https://en.wikipedia.org/wiki/Tammann_and_H%C3%BCttig_temperatures


    def score(self, rxn: ComputedReaction):
        phases = [c.reduced_formula for c in rxn.reactants]
        non_gasses = [p for p in phases if p not in DEFAULT_GASES]
        mps = [self.phases.get_melting_point(p) for p in non_gasses]
        min_mp = min(mps)

        # Softplus adjustment
        delta_g_adjustment = softplus(-rxn.energy_per_atom)
        delta_g_adjustment = softplus(-(2*rxn.energy_per_atom + 0.8))

        if len(non_gasses) < len(phases):
            # Huttig
            return huttig_score_softplus(self.temp / min_mp) * delta_g_adjustment
        else:
            # Tamman
            return tamman_score_softplus(self.temp / min_mp) * delta_g_adjustment


class TammanHuttigScoreErf(BasicScore):
    # https://en.wikipedia.org/wiki/Tammann_and_H%C3%BCttig_temperatures


    def score(self, rxn: ComputedReaction):
        phases = [c.reduced_formula for c in rxn.reactants]
        non_gasses = [p for p in phases if p not in self.phases.gas_phases]
        mps = [self.phases.get_melting_point(p) for p in non_gasses]
        min_mp = min(mps)

        # Softplus adjustment
        delta_g_adjustment = erf(rxn.energy_per_atom)

        if len(non_gasses) < len(phases):
            # Huttig
            return huttig_score_softplus(self.temp / min_mp) * delta_g_adjustment
        else:
            # Tamman
            return tamman_score_softplus(self.temp / min_mp) * delta_g_adjustment
    


def score_rxns(reactions: ReactionSet, scorer: BasicScore, phase_set: SolidPhaseSet = None):
    scored_reactions = []

    for rxn in tqdm(reactions.get_rxns(), desc=f"Scoring reactions... at temp {scorer.temp}"):
        reactants = [r.reduced_formula for r in rxn.reactants]
        non_gases = [r for r in reactants if r not in phase_set.gas_phases]
        if len(non_gases) > 0:
            scored_rxn = ScoredReaction.from_rxn_network(scorer.score(rxn), rxn, phase_set.volumes)
            scored_reactions.append(scored_rxn)

    return scored_reactions