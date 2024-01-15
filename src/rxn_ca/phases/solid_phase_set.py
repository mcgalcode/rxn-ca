from pylattica.discrete.phase_set import PhaseSet
from rxn_network.reactions.reaction_set import ReactionSet
from mp_api.client import MPRester

from typing import List, Dict

import math
import pkg_resources
import pandas as pd
from monty.serialization import loadfn
from .gasses import DEFAULT_GASES

from ..utilities.helpers import normalize_dict, add_values_to_dict_by_addition

from pymatgen.core.composition import Composition

class SolidPhaseSet(PhaseSet):

    FREE_SPACE = "Free Space"

    @classmethod
    def from_dict(cls, set_dict):
        return cls(
            set_dict["phases"],
            set_dict["volumes"],
            set_dict["gas_phases"],
            set_dict["melting_points"],
            set_dict["experimentally_observed"],
        )
    
    @classmethod
    def from_rxn_set(cls, rxn_set: ReactionSet, gas_phases: List[str] = DEFAULT_GASES):
        all_phases = list(set([e.composition.reduced_formula for e in rxn_set.entries]))
        return cls.from_phase_list(all_phases, gas_phases=gas_phases)
    
    @classmethod
    def from_phase_list(cls, solid_phases: List[str], gas_phases: List[str] = DEFAULT_GASES):
        vols = get_phase_vols(solid_phases)
        mps = get_melting_points(solid_phases)
        obs = get_experimentally_observed(solid_phases)
        return cls(solid_phases, gas_phases = gas_phases, volumes=vols, melting_points=mps, experimentally_observed=obs)

    def __init__(self, phases: List[str],
                       volumes: Dict[str, float],
                       gas_phases: List[str] = DEFAULT_GASES,
                       melting_points: Dict[str, float] = None,
                       experimentally_observed: Dict[str, bool] = None):
        phases = list(set(phases + [SolidPhaseSet.FREE_SPACE]))
        self.gas_phases: List[str] = gas_phases
        self.volumes: Dict[str, float] = volumes
        self.melting_points: Dict[str, float] = melting_points
        self.experimentally_observed: Dict[str, bool] = experimentally_observed
        super().__init__(phases)

    def get_vol(self, phase: str) -> float:
        """Returns the molar volume associated with the supplied phase.

        Args:
            phase (str): The formula of the phase of interest

        Returns:
            float: The molar volume
        """
        return self.volumes.get(phase)
    
    def get_melting_point(self, phase: str) -> float:
        """Returns the machine learning estimated melting point of the supplied phase.

        Args:
            phase (str): The formula of the phase of interest

        Returns:
            float: The melting point
        """
        return self.melting_points.get(phase)
    
    def is_theoretical(self, phase: str) -> bool:
        """Indicates whether or not the phase is marked as theoretical in MP

        Args:
            phase (str): The phase of interest

        Returns:
            bool:
        """
        return not self.experimentally_observed[phase]
    
    def is_gas(self, phase: str) -> bool:
        """Indicates whether or not the supplied phase is a gas

        Args:
            phase (str): The formula of the phase

        Returns:
            bool: Whether or not it is a gas
        """
        return phase in self.gas_phases

    def get_melted_phases(self, temp: int) -> List[str]:
        """Return the subset of phases which are melted at the supplied temp

        Args:
            temp (int): The temperature

        Returns:
            List[str]: The melted phases
        """
        return [p for p in self.phases if self.is_melted(p, temp)]
    
    def is_melted(self, phase: str, temp: int) -> bool:
        """Indicates whether or not the supplied phase is melted at the supplied temp

        Args:
            phase (str): The phase
            temp (int): The temperature

        Returns:
            bool: Is it melted?
        """
        return temp > self.get_melting_point(phase)

    def get_theoretical_phases(self) -> List[str]:
        """Returns the phases inside this SolidPhaseSet that are marked
        as theoretical in MP

        Returns:
            List[str]:
        """
        return [phase for phase in self.phases if phase != SolidPhaseSet.FREE_SPACE and self.is_theoretical(phase)]
    
    def mole_amts_to_vols(self, mol_amts: Dict[str, float]) -> Dict[str, float]:
        """Converts amounts expressed in moles to their equivalent volumes
        using the molar volumes of each phase supplied by MP

        Args:
            mol_amts (Dict[str, float]): A map of phase to # of moles

        Returns:
            Dict[str, float]: A map of phase to # of vols
        """
        return { phase: self.moles_to_vol(amt, phase) for phase, amt in mol_amts.items() }

    def moles_to_vol(self, moles: float, phase: str) -> float:
        """Converts a number of moles of a phase to the equivalent
        volume of that phase

        Args:
            moles (float): The amount to convert
            phase (str): The formula of the phase

        Returns:
            float: The equivalent volume
        """
        return moles * self.get_vol(phase)

    def vol_to_moles(self, vol: float, phase: str) -> float:
        """Converts a volume to a number of moles

        Args:
            vol (float): The volume to
            phase (str): The formula of the phase

        Returns:
            float: The number of moles equivalent to the supplied volume
        """
        return vol / self.get_vol(phase)

    def vol_amts_to_moles(self, vol_amts: Dict[str, float]) -> Dict[str, float]:
        """Converts a mapping of phase to volume to a mapping of phase to mole amount

        Args:
            vol_amts (Dict[str, float]): A mapping of phase formula to volume

        Returns:
            Dict[str, float]: A mapping of phase formula to mole amount
        """
        return { phase: self.vol_to_moles(amt, phase) for phase, amt in vol_amts.items() }

    def mole_amts_to_el_amts(self, mole_amts: Dict[str, float]) -> Dict[str, float]:
        """Calculates the moles of each element present in a collection of phases

        Args:
            mole_amts (Dict[str, float]): A mapping from phase formula to molar amount

        Returns:
            Dict[str, float]: A mapping from element to molar amount
        """
        elemental_amounts = {}

        for phase, moles in mole_amts.items():
            comp = Composition(phase)
            comp_dict = comp.as_dict()
            # multiply the moles of the current _phase_ by the stoich coeff. of the each element in the composition
            # i.e. 1 mole of TiO2 has 2 moles of O and 1 mole of Ti
            scaled_comp_dict = { el: comp_amt * moles for el, comp_amt in comp_dict.items()}
            add_values_to_dict_by_addition(elemental_amounts, scaled_comp_dict)
        
        return elemental_amounts
    
    def mole_amts_to_el_fracs(self, mole_amts: Dict[str, float]) -> Dict[str, float]:
        """Gives fractional elemental composition from molar amts

        Args:
            mole_amts (Dict[str, float]): A map of phase to number of

        Returns:
            Dict[str, float]: A map of element to fractional amount
        """
        return normalize_dict(self.mole_amts_to_el_amts(mole_amts))

    def vol_amts_to_el_amts(self, vol_amts: Dict[str, float]) -> Dict[str, float]:
        """Converts a set of phases with associated volumes to the moles of each
        element present in the set

        Args:
            vol_amts (Dict[str, float]): A map from phase formula to volume

        Returns:
            Dict[str, float]: A map from element to number of moles
        """
        mol_amts = self.vol_amts_to_moles(vol_amts)
        return self.mole_amts_to_el_amts(mol_amts)

    def as_dict(self):
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "phases": self.phases,
            "volumes": self.volumes,
            "gas_phases": self.gas_phases,
            "melting_points": self.melting_points,
            "experimentally_observed": self.experimentally_observed
        }

def get_phase_vols(phases: List[str]) -> Dict[str, float]:
    """Retrieves molar volumes for a list of phases from MP.

    Args:
        phases (List[str]): The phases of interest

    Returns:
        Dict[str, float]: A map from formula to molar volume
    """

    volumes = {}

    with MPRester() as mpr:
        res = mpr.summary.search(
            formula = phases,
            fields=["structure", "composition", "task_id", "energy_above_hull"]
        )
        min_energies = {}
        for item in res:
            struct = item.structure
            comp = item.composition
            red_form = comp.reduced_formula
            
            if item.energy_above_hull < min_energies.get(red_form, math.inf):
                volumes[comp.reduced_formula] = struct.volume / comp.get_reduced_composition_and_factor()[1]
                min_energies[red_form] = item.energy_above_hull

    return volumes

def get_experimentally_observed(phases) -> Dict[str, bool]:
    """Retrieves the "theoretical" flag for each of the supplied phases from MP.

    Args:
        phases (List[str]): The phases of interest

    Returns:
        Dict[str, bool]: A map from formula to whether or not the formula has been experimentally observed
    """
    experimentally_observed = {}

    with MPRester() as mpr:
        
        res = mpr.summary.search(
            formula = phases,
            fields=["theoretical", "composition"]
        )
        for item in res:
            form = item.composition.reduced_formula
            prev = experimentally_observed.get(form)
            if not prev:
                experimentally_observed[form] = not item.theoretical
    
    return experimentally_observed

def get_melting_points(phases: List[str]) -> Dict[str, float]:
    """Retrieves the machine learning melting point for the supplied phase

    Args:
        phases (List[str]): The formulas of the phases of interest

    Returns:
        Dict[str, float]: A map of formula to melting point
    """
    mp_json = pkg_resources.resource_filename("rxn_ca.reactions", "melting_points_df_08_08_23.json")
    melting_pt_data = pd.DataFrame(loadfn(mp_json))   # Note: temps in Kelvin
    mps = {}
    for p in phases:
        mps[p] = int(melting_pt_data[melting_pt_data["reduced_formula"] == p].iloc[0]["melting_point"])
    return mps