from __future__ import annotations

import math
from abc import ABC, abstractmethod
from typing import Optional

import numpy as np
from pymatgen.core import Composition
from tqdm import tqdm

from rxn_network.entries.entry_set import GibbsEntrySet
from rxn_network.reactions.computed import ComputedReaction
from rxn_network.reactions.reaction_set import ReactionSet
from rxn_network.thermo.chempot_diagram import ChemicalPotentialDiagram

from ..phases.gasses import DEFAULT_GASES
from ..phases.solid_phase_set import SolidPhaseSet
from ..utilities.transport import TransportDatabase
from .scored_reaction import ScoredReaction

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


class DiffusionScorer(BasicScore):
    """
    A scorer that uses the diffusion coefficients of the species in the reaction to score the reaction.

    Args:
        phase_set: Phase set for volumes and melting points.
        entry_set: Gibbs entries.
        temp: Temperature (K).
        precursor_size: Precursor size in meters (default 1e-8).
        scale_factor: Scales flux contribution (default 5e12).
        chemsys: Chemical system string for transport data.
        self_diffusion: If True, use self Onsager tensor for fluxes.
        tamman_factor: Multiplier on T/Tm for onset factor.
        mode: How to pick chemical-potential distance between reactant domains ('min'|'max'|'mean').
    """

    # https://arxiv.org/abs/2501.08560 (Section 4.5)

    def __init__(
        self,
        phase_set: SolidPhaseSet,
        entry_set: GibbsEntrySet,
        temp: int,
        precursor_size: float = 1e-8,
        scale_factor: float = 5e12,
        chemsys: Optional[list[str] | str] = None,
        self_diffusion: bool = False,
        tamman_factor: float = 1.0,
        mode: str = "min",
    ):
        self.phase_set = phase_set
        self.temp = temp
        self.chem_pot_diagram = ChemicalPotentialDiagram(
            entry_set.get_entries_with_new_temperature(self.temp)
        )
        self.max_flux = 1e20 * KB * self.temp
        self.min_flux = 1e8 * KB * self.temp
        self.chemsys = chemsys if chemsys else self.chem_pot_diagram.chemical_system
        self.precursor_size = precursor_size
        self.scale_factor = scale_factor
        self.mode = mode
        self.tamman_factor = tamman_factor
        self.self_diffusion = self_diffusion
        self.transport_db: TransportDatabase = TransportDatabase.from_chemsys(self.chemsys)

    def score(self, reaction: ComputedReaction):
        dG = reaction.energy_per_atom
        
        phases = [c.reduced_formula for c in reaction.reactants]
        non_gasses = [p for p in phases if p not in DEFAULT_GASES]
        mps = [self.phase_set.get_melting_point(p) for p in non_gasses]
        min_mp = np.min(mps)

        system_factor = (
            1 / self.precursor_size**2 / self.scale_factor / Na / KB / self.temp
        )

        onset_factor = tamman_score_softplus(self.temp / min_mp * self.tamman_factor)
        
        if dG > 0.01:
            return erf(dG)*onset_factor

        if len(reaction.reactants) > 1:
            reactants_dict = {c.reduced_formula: v for c, v in reaction.reactant_coeffs.items()}
            products_dict = {c.reduced_formula: v for c, v in reaction.product_coeffs.items()}
        else:
            reactants_dict = {c.reduced_formula: -v for c, v in reaction.product_coeffs.items()}
            products_dict = {c.reduced_formula: -v for c, v in reaction.reactant_coeffs.items()}

        mu_dict = mu_distance(list(reactants_dict.keys()), self.chem_pot_diagram, self.mode)

        all_fluxes = []
        mol_normalizer = sum(products_dict.values())

        for p, coeff in products_dict.items():
            el_fraction = {el: 0.0 for el in self.chem_pot_diagram.chemical_system.split("-")}
            for sp in Composition(p).elements:
                el_fraction[sp.symbol] = Composition(p).get_atomic_fraction(sp)

            vol = self.phase_set.get_vol(p)

            if p not in DEFAULT_GASES:
                melt_pt = self.phase_set.get_melting_point(p)
                if self.temp > melt_pt:
                    K_d = sum(
                        abs(self.max_flux / el_fraction[sp] * vol)
                        for sp in el_fraction
                        if el_fraction[sp] > 0
                    )
                    all_fluxes.append(K_d * coeff / mol_normalizer)
                    continue

            fluxes = self.transport_db.compute_fluxes(
                formula=Composition(p).reduced_formula,
                temperature=float(self.temp),
                mu=mu_dict,
                self_transport=self.self_diffusion,
            )
            K_d = sum(
                abs(fluxes[sp] / el_fraction[sp] / vol)
                for sp in fluxes
                if el_fraction.get(sp, 0.0) > 0
            )
            all_fluxes.append(K_d * coeff / mol_normalizer)

        selected_flux = np.max(all_fluxes)

        selectivity_factor = selected_flux * system_factor * erf(dG)
        selectivity_factor = (
            softplus(selectivity_factor)
            if selectivity_factor < 2
            else selectivity_factor
        )
        score = onset_factor * selectivity_factor
        return score if score < 20 else 20


def mu_distance(
    phases: list, chempot: ChemicalPotentialDiagram, mode: str = "min"
) -> dict[str, float]:
    """
    Compute the chemical potential distance vector between the stability domains
    of two phases and return it as a mapping of element symbol -> mu distance (eV).
    """
    if len(phases) > 2:
        phases = [p for p in phases if p not in DEFAULT_GASES]
        if len(phases) > 2:
            raise ValueError(
                "Only works for reactions between two solids or a solid and a gas"
            )

    domains = []
    for phase in phases:
        if phase in list(chempot.domains.keys()):
            domains.append(chempot.domains[phase])
        else:
            domains.append(chempot.metastable_domains[phase])

    if mode == "max":
        max_distance = 0
        ind_1, ind_2 = 0, 0
        for node1 in range(len(domains[0])):
            for node2 in range(len(domains[1])):
                distance = np.linalg.norm(domains[0][node1] - domains[1][node2])
                if distance > max_distance:
                    max_distance = distance
                    ind_1, ind_2 = node1, node2
        mu_vec = domains[0][ind_1] - domains[1][ind_2]

    elif mode == "min":
        min_distance = 1000
        ind_1, ind_2 = 0, 0
        for node1 in range(len(domains[0])):
            for node2 in range(len(domains[1])):
                distance = np.linalg.norm(domains[0][node1] - domains[1][node2])
                if distance < min_distance:
                    min_distance = distance
                    ind_1, ind_2 = node1, node2
        mu_vec = domains[0][ind_1] - domains[1][ind_2]

    elif mode == "mean":
        mu_vec = np.mean(domains[0], axis=0) - np.mean(domains[1], axis=0)

    else:
        raise ValueError(f"Unknown mode={mode!r}")

    return {str(el): float(val) for el, val in zip(chempot.elements, mu_vec)}


def score_rxns(
    reactions: ReactionSet, scorer: BasicScore, phase_set: Optional[SolidPhaseSet] = None
):
    scored_reactions = []

    for rxn in tqdm(
        reactions.get_rxns(), desc=f"Scoring reactions... at temp {scorer.temp}"
    ):
        reactants = [r.reduced_formula for r in rxn.reactants]
        non_gases = [r for r in reactants if r not in phase_set.gas_phases]
        if len(non_gases) > 0:
            scored_rxn = ScoredReaction.from_rxn_network(
                scorer.score(rxn), rxn, phase_set.volumes
            )
            scored_reactions.append(scored_rxn)

    return scored_reactions
