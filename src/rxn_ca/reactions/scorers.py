import math
from tqdm import tqdm
from abc import ABC, abstractmethod
from .scored_reaction import ScoredReaction

import numpy as np
import pandas as pd
import warnings
import pkg_resources

from ..phases.solid_phase_set import SolidPhaseSet
from pymatgen.core import Composition
from typing import List

from rxn_network.reactions.reaction_set import ReactionSet
from rxn_network.reactions.computed import ComputedReaction
from rxn_network.thermo.chempot_diagram import ChemicalPotentialDiagram
from rxn_network.entries.entry_set import GibbsEntrySet

from ..phases.gasses import DEFAULT_GASES

KB = 8.6173303e-5 # Boltzmann constant in eV/K
Na = 6.02214076e23 # Avogadro's number

def softplus(x):
    return 1/3 * math.log(1 + math.exp(3*x))

def tamman_score_exp(t_tm_ratio):
    return math.exp(4.82*(t_tm_ratio) - 3.21)

def tamman_score_softplus(t_tm_ratio):
    return math.log(1 + math.exp(14 * (t_tm_ratio - 0.8)))

def huttig_score_exp(t_tm_ratio):
    return math.exp(2.41*(t_tm_ratio) - 0.8)

def huttig_score_softplus(t_tm_ratio):
    return 0.25 * math.log(1 + math.exp(30 * (t_tm_ratio - 0.33)))

def erf(x):
    return 0.5 * (1 + math.erf(-35 * (x + 0.03)))

def tamman_erf_score(tm_ratio, delta_g):
    return tamman_score_softplus(tm_ratio) * erf(delta_g)

def huttig_erf_score(tm_ratio, delta_g):
    return huttig_score_softplus(tm_ratio) * erf(delta_g)
         
class BasicScore(ABC):

    def __init__(self, phase_set: SolidPhaseSet, temp: int = None):
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

        if len(non_gasses) == 1:
            # Huttig
            return huttig_score_softplus(self.temp / min_mp) * delta_g_adjustment
        else:
            # Tamman
            return tamman_score_softplus(self.temp / min_mp) * delta_g_adjustment

class GibbsErfScore(BasicScore):
    # https://en.wikipedia.org/wiki/Tammann_and_H%C3%BCttig_temperatures

    def score(self, rxn: ComputedReaction):
        return erf(rxn.energy_per_atom)

class TammanScore(BasicScore):
    # https://en.wikipedia.org/wiki/Tammann_and_H%C3%BCttig_temperatures

    def score(self, rxn: ComputedReaction):
        phases = [c.reduced_formula for c in rxn.reactants]
        non_gasses = [p for p in phases if p not in self.phases.gas_phases]
        mps = [self.phases.get_melting_point(p) for p in non_gasses]
        min_mp = min(mps)

        # Softplus adjustment
        delta_g_adjustment = erf(rxn.energy_per_atom)
        return tamman_score_softplus(self.temp / min_mp) * delta_g_adjustment  

class ConstantScore(BasicScore):

    def score(self, _):
        return 1.0

class GibbsErfScore(BasicScore):
        
    def score(self, rxn: ComputedReaction):
        return erf(rxn.energy_per_atom)
    

class TammanTightLinear(BasicScore):

    def score(self, rxn: ComputedReaction):
        phases = [c.reduced_formula for c in rxn.reactants]
        non_gasses = [p for p in phases if p not in self.phases.gas_phases]
        mps = [self.phases.get_melting_point(p) for p in non_gasses]
        min_mp = min(mps)

        def _score(x):
            return 1/2*(1 + math.erf(20*(x -0.6))) * (1/0.6*x)

        # Softplus adjustment
        delta_g_adjustment = erf(rxn.energy_per_atom)
        return _score(self.temp / min_mp) * delta_g_adjustment  

class DiffusionScorer(BasicScore):
    """
    A scorer that uses the diffusion coefficients of the species in the reaction to score the reaction.
    Args:
        phase_set (SolidPhaseSet): The phase set to use for phase information.
        chem_pot_diagram (ChemicalPotentialDiagram): The chemical potential diagram to use for chemical potential information.
        precursor_size (float): The size of the precursor in meters. Default is 1e-7 (0.1 micro-m).
        scale_factor (float): The scale factor to use for the scoring. Default is 5e13. This is used to scale the score to a reasonable range.
                            Modify this to ensure no math range errors. Might have to adjust based on the system.
        chemsys (list[str] | str): The chemical system to use for the diffusion coefficients. Default is None, which uses the chemical system of the chemical potential diagram.
        temp (int): The temperature to use for the scoring.
        self_diffusion (bool): Whether to use self-diffusion coefficients instead of the full transport tensor. Default is False.
    """
    # https://arxiv.org/abs/2501.08560 (Section 4.5)

    def __init__(self, 
                 phase_set: SolidPhaseSet, 
                 chem_pot_diagram: ChemicalPotentialDiagram, 
                 temp : int,
                 precursor_size : float = 1e-7, 
                 scale_factor : float = 5e13, 
                 chemsys : list[str] | str = None, 
                 self_diffusion=False):
        
        self.phase_set = phase_set
        self.chem_pot_diagram = chem_pot_diagram
        self.max_flux = 1e-3
        self.min_flux = 1e-10
        self.chemsys = chemsys if chemsys else chem_pot_diagram.chemical_system
        self.precursor_size = precursor_size
        self.scale_factor = scale_factor
        self.temp = temp
        self.diff_df = load_diffusivities(chemsys = self.chemsys, self_diffusion=self_diffusion)   # Note: temps in Kelvin
    
    def score(self, reaction: ComputedReaction):
        # This scorer assumes binary reactions and no gasses (A + B -> C + D + E..., not A -> B + C)
        if len(reaction.reactants) == 1:
            reactant_1 = reaction.reactants[0]
            reactant_2 = reaction.reactants[0]
        else:
            reactant_1 = reaction.reactants[0]
            reactant_2 = reaction.reactants[1]
        
        all_fluxes = []
        for p in reaction.products:
            interface = [c.reduced_formula for c in [reactant_1, p, reactant_2]]
            el_fraction = get_el_ratios(p)
            if len(el_fraction) < 2:
                el_fraction = [el_fraction[0], 0] # Edge case for single element products
            
            vol = self.phase_set.get_vol(p)

            if p not in DEFAULT_GASES:
                melt_pt = self.phase_set.get_melting_point(p)
                if self.temp > melt_pt:
                    fluxes = [self.max_flux, self.max_flux, self.max_flux] # Set fluxes to max if above melting point
                    all_fluxes.append(fluxes[0]*vol)
            
            fluxes = get_fluxes_across_interface(
                interface,
                self.diff_df,
                self.temp,
                self.chem_pot_diagram,
                chemsys=self.chemsys
            )
            all_fluxes.append(np.sum([np.abs(fluxes[i]/el_fraction[i]*vol) for i in range(len(fluxes))])) # TODO: check if indices match between Lij data and pmg structure species
        
        selected_flux = max([np.abs(f) for f in all_fluxes])
        dG = reaction.energy_per_atom
        phases = [c.reduced_formula for c in reaction.reactants]
        non_gasses = [p for p in phases if p not in DEFAULT_GASES]
        mps = [self.phase_set.get_melting_point(p) for p in non_gasses]
        min_mp = np.min(mps)        
        
        system_factor = 1/self.precursor_size**2/self.scale_factor/Na/KB/self.temp
        onset_factor = tamman_score_softplus(self.temp/min_mp)
        selectivity_factor = softplus(selected_flux*system_factor*erf(dG))
        score = onset_factor*selectivity_factor
        return score if score < 20 else 20 # Cap score at 20, anything higher doesnt really help with selectivity/rate

def mu_distance(phases : list, chempot : ChemicalPotentialDiagram, mode : str ='min'):
    domains = []
    for phase in phases:
        if phase in list(chempot.domains.keys()):
            domains.append(chempot.domains[phase])
        else:
            domains.append(chempot.metastable_domains[phase])

    if mode == 'max':
        max_distance = 0
        ind_1, ind_2 = 0, 0
        for node1 in range(len(domains[0])):
            for node2 in range(len(domains[1])):
                distance = np.linalg.norm(domains[0][node1] - domains[1][node2])
                if distance > max_distance:
                    max_distance = distance
                    ind_1, ind_2 = node1, node2
        return domains[0][ind_1] - domains[1][ind_2] 
    
    if mode == 'min':
        min_distance = 1000
        ind_1, ind_2 = 0, 0
        for node1 in range(len(domains[0])):
            for node2 in range(len(domains[1])):
                distance = np.linalg.norm(domains[0][node1] - domains[1][node2])
                if distance < min_distance:
                    min_distance = distance
                    ind_1, ind_2 = node1, node2
        return domains[0][ind_1] - domains[1][ind_2]
    
    if mode == 'mean':
        return np.mean(domains[0], axis=0) - np.mean(domains[1], axis=0)

def get_fluxes_across_interface(interface : list, L_data : pd.DataFrame, temperature : int, chempot : ChemicalPotentialDiagram, mode : str ='min', chemsys : str = None):
    
    if len(chempot.chemical_system.split('-')) > 3:
        species = chempot.chemical_system.split('-')
        if 'C' in species:
            del_index = species.index('C')
        mu = np.delete(mu_distance([interface[0], interface[2]], chempot, mode), del_index)
    else:
        mu = mu_distance([interface[0], interface[2]], chempot, mode)
    
    if chempot.chemical_system != chemsys:  # Reorder mu to match the chemical system, i.e., order in which the L_ij values are stored (need better way of handling this)
        returned_order = chempot.chemical_system.split('-')
        indices = [returned_order.index(label) for label in chemsys.split('-')]
        mu = mu[indices]
        
    temps = [unstringify_temp(t) for t in list(set(L_data['Temperature']))]
    diffs = np.abs(np.array(temps) - temperature)
    closest_temp_idx = np.argmin(diffs)
    closest_temp = stringify_temp(temps[closest_temp_idx])
    formula_data =  L_data[(L_data['Formula'] == interface[1]) & (L_data['Temperature'] == closest_temp)]
    # Extract L_ij values
    
    L_ij_values = np.eye(len(mu))*1e15
    
    for i in range(len(mu)):
        for j in range(i, len(mu)):
            try:
                L_ij_values[i][j] = formula_data[f'L{i}{j}'].values[0]
                L_ij_values[j][i] = formula_data[f'L{i}{j}'].values[0]
            except IndexError:
                warnings.warn(f'No L_ij data for {interface[1]} at {temperature}K, setting to 1e15')
                
    return np.dot(L_ij_values, mu)

def load_diffusivities(chemsys : list = None, self_diffusion : bool = False):
    if chemsys in ['Ba-Ti-O', 'Ba-Ti-O-C']:
        if self_diffusion:
            try:
                diff_csv_path = pkg_resources.resource_filename("rxn_ca.phases", "batio_self_coeffs.csv")
            except FileNotFoundError:
                raise FileNotFoundError("Transport data is not available yet. Please contact the authors if you would like to use this feature!")
        else:
            try:
                diff_csv_path = pkg_resources.resource_filename("rxn_ca.phases", "batio_transport_coeffs.csv")
            except FileNotFoundError:
                raise FileNotFoundError("Transport data is not available yet. Please contact the authors if you would like to use this feature!")

    return pd.read_csv(diff_csv_path)

def get_el_ratios(formula_string):
    return [Composition(formula_string).get_atomic_fraction(el) for el in Composition(formula_string).elements]    

def stringify_temp(temp):
    return f"{temp}.0K"

def unstringify_temp(temp_str):
    return int(temp_str.split('.')[0])

def score_rxns(reactions: ReactionSet, scorer: BasicScore, phase_set: SolidPhaseSet = None):
    scored_reactions = []

    for rxn in tqdm(reactions.get_rxns(), desc=f"Scoring reactions... at temp {scorer.temp}"):
        reactants = [r.reduced_formula for r in rxn.reactants]
        non_gases = [r for r in reactants if r not in phase_set.gas_phases]
        if len(non_gases) > 0:
            scored_rxn = ScoredReaction.from_rxn_network(scorer.score(rxn), rxn, phase_set.volumes)
            scored_reactions.append(scored_rxn)

    return scored_reactions