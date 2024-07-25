from pylattica.discrete.phase_set import PhaseSet
from rxn_network.entries.entry_set import GibbsEntrySet
from mp_api.client import MPRester
from pymatgen.core.structure import Structure
from rxn_network.entries.experimental import ExperimentalReferenceEntry
from rxn_network.entries.gibbs import GibbsComputedEntry
from rxn_network.entries.utils import process_entries
from pymatgen.core.composition import Composition
from typing import List, Dict, Any

import executing
import itertools
import pkg_resources
import pandas as pd
from monty.serialization import loadfn
from .gasses import DEFAULT_GASES

import copy
from enum import Enum
import requests

from ..utilities.helpers import normalize_dict, add_values_to_dict_by_addition


from pymatgen.core.composition import Composition

def process_composition(comp_str):
    return Composition(comp_str).reduced_formula

def process_composition_list(comp_list):
    return [process_composition(c) for c in comp_list]

def process_composition_dict(comp_dict):
    return { process_composition(c): v for c, v in comp_dict.items() }

class MatterPhase(Enum):

    SOLID  = "SOLID"
    LIQUID = "LIQUID"
    GAS    = "GAS"

class SolidPhaseSet(PhaseSet):

    FREE_SPACE = "Free Space"

    @classmethod
    def from_dict(cls, set_dict):
        return cls(
            phases=set_dict["phases"],
            volumes=set_dict["volumes"],
            gas_phases=set_dict["gas_phases"],
            densities=set_dict["densities"],
            melting_points=set_dict["melting_points"],
            experimentally_observed=set_dict["experimentally_observed"],
            phase_metadata=set_dict.get("phase_metadata"),
        )
   
    @classmethod
    def from_entry_set(cls,
                       entry_set: GibbsEntrySet,
                       entry_metadata: Dict = {},
                       gas_phases: List[str] = DEFAULT_GASES,
                       phase_metadata: Dict = None):
        """Constructs a SolidPhaseSet from an EntrySet object. Expects the entry
        set itself and allows direct specification of the metadata for those entries
        along side it (metadata being whether or not that entry is experimentally observed,
        its melting point, its volume, and its density)

        Parameters
        ----------
        entry_set : GibbsEntrySet
            _description_
        entry_metadata : Dict
            _description_

        Returns
        -------
        _type_
            _description_
        """
        phase_names = [e.composition.reduced_formula for e in entry_set.entries_list]
        print("Retrieving entry volumes...")
        volumes, densities = get_densities_and_vols_from_entry_set(entry_set)
        print("Retrieving entry obs....")
        exp_obs = get_exp_observ(entry_set)
        print("Retrieving entry mps...")

        known_mps = process_composition_dict(entry_metadata.get("melting_points", {}))
        remaining = list(set(phase_names) - set(known_mps.keys()))
        melting_points = get_melting_points(remaining)
        
        volumes = { **volumes, **process_composition_dict(entry_metadata.get("volumes", {})) }
        densities = { **densities, **process_composition_dict(entry_metadata.get("densities", {})) }
        melting_points = { **melting_points, **known_mps }
        exp_obs = { **exp_obs, **process_composition_dict(entry_metadata.get("experimentally_observed", {})) }

        return cls(
            phase_names,
            gas_phases = gas_phases,
            volumes=volumes,
            melting_points=melting_points,
            experimentally_observed=exp_obs,
            densities=densities,
            phase_metadata=phase_metadata
        )
    
    @classmethod
    def from_phase_list(cls, solid_phases: List[str], entry_metadata: Dict = {}, gas_phases: List[str] = DEFAULT_GASES, phase_metadata: Dict = None):
        entry_set = get_entry_set_from_phase_list(solid_phases)
        return cls.from_entry_set(entry_set, entry_metadata=entry_metadata, gas_phases=gas_phases, phase_metadata=phase_metadata)

    def __init__(self, phases: List[str],
                       volumes: Dict[str, float],
                       gas_phases: List[str] = DEFAULT_GASES,
                       densities: Dict[str, float] = None,
                       melting_points: Dict[str, float] = None,
                       experimentally_observed: Dict[str, bool] = None,
                       phase_metadata: Dict = None):
        phases = process_composition_list(list(set(phases)))
        self.gas_phases: List[str] = process_composition_list(gas_phases)
        self.volumes: Dict[str, float] = process_composition_dict(volumes)
        self.melting_points: Dict[str, float] = process_composition_dict(melting_points)
        self.experimentally_observed: Dict[str, bool] = process_composition_dict(experimentally_observed)
        self.densities: Dict[str, float] = process_composition_dict(densities)
        self.phase_metadata = phase_metadata
        super().__init__(phases)

    def get_vol(self, phase: str) -> float:
        """Returns the molar volume associated with the supplied phase.

        Args:
            phase (str): The formula of the phase of interest

        Returns:
            float: The molar volume
        """
        return self.volumes.get(process_composition(phase))
    
    def get_melting_point(self, phase: str) -> float:
        """Returns the machine learning estimated melting point of the supplied phase.

        Args:
            phase (str): The formula of the phase of interest

        Returns:
            float: The melting point
        """
        return self.melting_points.get(process_composition(phase))
    
    def get_density(self, phase: str) -> float:
        """Returns the density of the supplied phase.

        Args:
            phase (str): The formula of the phase of interest

        Returns:
            float: The density
        """
        return self.densities.get(process_composition(phase))
    
    def is_theoretical(self, phase: str) -> bool:
        """Indicates whether or not the phase is marked as theoretical in MP

        Args:
            phase (str): The phase of interest

        Returns:
            bool:
        """
        return not self.experimentally_observed.get(process_composition(phase), False)
    
    def is_gas(self, phase: str) -> bool:
        """Indicates whether or not the supplied phase is a gas

        Args:
            phase (str): The formula of the phase

        Returns:
            bool: Whether or not it is a gas
        """
        return process_composition(phase) in self.gas_phases

    def get_matter_phase(self, phase: str, temp: int = None):
        if self.is_gas(phase):
            return MatterPhase.GAS
        if temp is None:
            return MatterPhase.SOLID
        elif self.is_melted(phase, temp):
            return MatterPhase.LIQUID
        else:
            return MatterPhase.SOLID

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
        return temp > self.get_melting_point(process_composition(phase))
    
    def is_non_gaseous_el(self, phase: str) -> bool:
        c = Composition(phase)

        if len(c.elements) > 1:
            return False
        elif c.reduced_formula not in self.gas_phases:
            return True


    def get_theoretical_phases(self) -> List[str]:
        """Returns the phases inside this SolidPhaseSet that are marked
        as theoretical in MP

        Returns:
            List[str]:
        """
        return [phase for phase in self.phases if self.is_theoretical(phase)]
    
    def get_experimentally_observed_phases(self) -> List[str]:
        """Returns the phases inside this SolidPhaseSet that are marked
        as theoretical in MP

        Returns:
            List[str]:
        """
        return [phase for phase in self.phases if not self.is_theoretical(phase)]
    
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

    def vol_to_moles(self, vol: float, phase: str, should_round: int=False) -> float:
        """Converts a volume to a number of moles

        Args:
            vol (float): The volume to
            phase (str): The formula of the phase

        Returns:
            float: The number of moles equivalent to the supplied volume
        """
        if should_round:
            return round(vol / self.get_vol(phase), should_round)
        else:
            return vol / self.get_vol(phase)

    def vol_amts_to_moles(self, vol_amts: Dict[str, float], should_round=False) -> Dict[str, float]:
        """Converts a mapping of phase to volume to a mapping of phase to mole amount

        Args:
            vol_amts (Dict[str, float]): A mapping of phase formula to volume

        Returns:
            Dict[str, float]: A mapping of phase formula to mole amount
        """
        return { phase: self.vol_to_moles(amt, phase, should_round=should_round) for phase, amt in vol_amts.items() }

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
            "densities": self.densities,
            "melting_points": self.melting_points,
            "experimentally_observed": self.experimentally_observed,
            "phase_metadata": self.phase_metadata
        }
    
    def __iter__(self):
        return self.phases

    def __len__(self):
        return len(self.phases)


def get_entry_set_from_phase_list(phases: List[str]) -> List:
    els = set()
    for p in phases:
        comp = Composition(p)
        for el in comp.elements:
            els.add(el)
    
    search_phases = copy.copy(phases)
    search_phases.extend([str(e) for e in els])

    with MPRester() as mpr:
        entries = mpr.get_entries(
            search_phases,
            additional_criteria={"thermo_types": ["GGA_GGA+U"]},
        )

    entry_set = process_entries(entries, 300, 0.1, formulas_to_include=phases)

    for sp in search_phases:
        if sp not in phases:
            reduced_sp_form = Composition(sp).reduced_formula
            matching_entries = [e for e in entry_set.entries_list if e.composition.reduced_formula == reduced_sp_form]
            for e in matching_entries:
                entry_set.discard(e)

    return entry_set

def get_density_from_struct(struct: Structure):
    return struct.volume / struct.composition.weight

def get_molar_volume_from_struct(struct: Structure):
    return struct.volume / struct.composition.get_reduced_composition_and_factor()[1]

def get_molar_volumes_from_structures(structs: List[Structure]) -> Dict[str, float]:
    return { s.composition.reduced_formula: get_molar_volume_from_struct(s) for s in structs}

def get_densities_from_structures(structs: List[Structure]) -> Dict[str, float]:
    return { s.composition.reduced_formula: get_density_from_struct(s) for s in structs}

def get_densities_and_vols_from_entry_set(eset):
    phases = [e.composition.reduced_formula for e in eset.entries]
    densities = {}
    vols = {}
    
    phases_without_vol = []
    
    # First, we calculate vol/dens using volume information attached to the
    # entries themselves
    for p in phases:
        e = eset.get_min_entry_by_formula(p)
        comp = e.composition
        reduced_composition = comp.get_reduced_composition_and_factor()[0] 
        if isinstance(e, GibbsComputedEntry):    
            vols[p] = e.volume_per_atom * reduced_composition.num_atoms
            densities[p] = float(reduced_composition.weight) / vols[p]
        else:
            phases_without_vol.append(p)
    
    # For InterpolatedEntry and ExperimentalReferenceEntry items, we see if MP
    # Can help us supply any volumes
    if len(phases_without_vol) > 0:
        with MPRester() as mpr:
            res = mpr.summary.search(formula=phases_without_vol, fields=["structure", "composition", "formation_energy_per_atom"])


        min_e_structs = []
        for _, group in itertools.groupby(res, key=lambda i: i.composition.reduced_formula):
            min_entry = min(group, key=lambda i: i.formation_energy_per_atom)
            min_e_structs.append(min_entry.structure)

        vols = {**vols, **get_molar_volumes_from_structures(min_e_structs)}
        densities = {**densities, **get_densities_from_structures(min_e_structs)}
        
    phases_without_vol = [p for p in phases if p not in vols]
    if len(phases_without_vol) > 0:
        # For remaining items, should any be left without volume/density information, we
        # use the average molar volume
        avg_vol_per_atom = sum([v / Composition(p).num_atoms for p, v in vols.items()]) / len(vols)
        for p in phases_without_vol:
            comp = Composition(p)
            vol = comp.num_atoms * avg_vol_per_atom
            vols[p] = vol
            densities[p] = comp.weight / vol

    return vols, densities

def get_exp_observ(eset):
    phases = [e.composition.reduced_formula for e in eset.entries]
    with MPRester() as mpr:
        res = mpr.summary.search(
            formula = phases,
            fields=["theoretical", "composition"]
        )

    experimentally_observed = {}
    for comp, group in itertools.groupby(res, key=lambda i: i.composition.reduced_formula):
        experimentally_observed[comp] = any([not i.theoretical for i in group])
        
    for phase in phases:
        if phase not in experimentally_observed:
            experimentally_observed[phase] = False
    
    for phase in phases:
        entry = eset.get_min_entry_by_formula(phase)
        if isinstance(entry, ExperimentalReferenceEntry):
            experimentally_observed[phase] = True
    
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
    unknown_phases = []
    for p in phases:
        try:
            mps[p] = int(melting_pt_data[melting_pt_data["reduced_formula"] == p].iloc[0]["melting_point"])
        except (executing.executing.NotOneValueFound, IndexError):
            unknown_phases.append(p)
    
    if len(unknown_phases) > 0:
        print(f"Couldn't find {len(unknown_phases)} in mem... Using API for retrieval.")
        predicted = predict_melting_points_api(unknown_phases)
        mps = { **mps, **predicted }

    return mps

def predict_melting_points_api(phases: List[str]) -> Dict[str, float]:
    data = [{"9": p} for p in phases]
    url = 'http://206.207.50.58:5007/MT_ML_Qijun_Hong_Predict_noNN'
    res = requests.post(url, json=data)
    result = {}
    for phase, item in zip(phases, res.json()):
        result[phase] = item['melting temperature']
    return result