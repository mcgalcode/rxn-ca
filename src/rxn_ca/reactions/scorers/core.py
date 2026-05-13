from __future__ import annotations

import math
from abc import ABC, abstractmethod
from typing import Optional

from rxn_network.reactions.computed import ComputedReaction

from ...phases.gasses import DEFAULT_GASES
from ...phases.solid_phase_set import SolidPhaseSet

KB = 8.6173303e-5  # Boltzmann constant in eV/K
Na = 6.02214076e23  # Avogadro's number


def softplus(x):
    return 1 / 3 * math.log(1 + math.exp(3 * x))


def tamman_score_exp(t_tm_ratio):
    return math.exp(4.82 * (t_tm_ratio) - 3.21)


def tamman_score_softplus(t_tm_ratio):
    return math.log(1 + math.exp(14 * (t_tm_ratio - 0.8)))


def huttig_score_exp(t_tm_ratio):
    return math.exp(2.41 * (t_tm_ratio) - 0.8)


def huttig_score_softplus(t_tm_ratio):
    return 0.25 * math.log(1 + math.exp(30 * (t_tm_ratio - 0.33)))


def erf(x):
    return 0.5 * (1 + math.erf(-35 * (x + 0.03)))


class BasicScore(ABC):
    def __init__(self, phase_set: SolidPhaseSet, temp: Optional[int] = None):
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

        delta_g_adjustment = softplus(-(2 * rxn.energy_per_atom + 0.8))

        if len(non_gasses) < len(phases):
            return huttig_score_exp(self.temp / min_mp) * delta_g_adjustment
        return tamman_score_exp(self.temp / min_mp) * delta_g_adjustment


class TammanHuttigScoreSoftplus(BasicScore):
    # https://en.wikipedia.org/wiki/Tammann_and_H%C3%BCttig_temperatures

    def score(self, rxn: ComputedReaction):
        phases = [c.reduced_formula for c in rxn.reactants]
        non_gasses = [p for p in phases if p not in DEFAULT_GASES]
        mps = [self.phases.get_melting_point(p) for p in non_gasses]
        min_mp = min(mps)

        delta_g_adjustment = softplus(-(2 * rxn.energy_per_atom + 0.8))

        if len(non_gasses) < len(phases):
            return huttig_score_softplus(self.temp / min_mp) * delta_g_adjustment
        return tamman_score_softplus(self.temp / min_mp) * delta_g_adjustment


class TammanHuttigScoreErf(BasicScore):
    # https://en.wikipedia.org/wiki/Tammann_and_H%C3%BCttig_temperatures

    def score(self, rxn: ComputedReaction):
        phases = [c.reduced_formula for c in rxn.reactants]
        non_gasses = [p for p in phases if p not in self.phases.gas_phases]
        mps = [self.phases.get_melting_point(p) for p in non_gasses]
        min_mp = min(mps)

        delta_g_adjustment = erf(rxn.energy_per_atom)

        if len(non_gasses) == 1:
            return huttig_score_softplus(self.temp / min_mp) * delta_g_adjustment
        return tamman_score_softplus(self.temp / min_mp) * delta_g_adjustment


class GibbsErfScore(BasicScore):
    def score(self, rxn: ComputedReaction):
        return erf(rxn.energy_per_atom)


class TammanScore(BasicScore):
    # https://en.wikipedia.org/wiki/Tammann_and_H%C3%BCttig_temperatures

    def score(self, rxn: ComputedReaction):
        phases = [c.reduced_formula for c in rxn.reactants]
        non_gasses = [p for p in phases if p not in self.phases.gas_phases]
        mps = [self.phases.get_melting_point(p) for p in non_gasses]
        min_mp = min(mps)

        delta_g_adjustment = erf(rxn.energy_per_atom)
        return tamman_score_softplus(self.temp / min_mp) * delta_g_adjustment


class ConstantScore(BasicScore):
    def score(self, _):
        return 1.0


class TammanTightLinear(BasicScore):
    def score(self, rxn: ComputedReaction):
        phases = [c.reduced_formula for c in rxn.reactants]
        non_gasses = [p for p in phases if p not in self.phases.gas_phases]
        mps = [self.phases.get_melting_point(p) for p in non_gasses]
        min_mp = min(mps)

        def _score(x):
            return 1 / 2 * (1 + math.erf(20 * (x - 0.6))) * (1 / 0.6 * x)

        delta_g_adjustment = erf(rxn.energy_per_atom)
        return _score(self.temp / min_mp) * delta_g_adjustment
