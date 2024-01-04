from pylattica.discrete.phase_set import PhaseSet
from rxn_network.reactions.reaction_set import ReactionSet
from mp_api.client import MPRester
import math
import pkg_resources
import pandas as pd
from monty.serialization import loadfn
from .gasses import DEFAULT_GASES

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
    def from_rxn_set(cls, rxn_set: ReactionSet, gas_phases = DEFAULT_GASES):
        all_phases = list(set([e.composition.reduced_formula for e in rxn_set.entries]))
        vols = get_phase_vols(all_phases)
        mps = get_melting_points(all_phases)
        obs = get_experimentally_observed(all_phases)
        return cls(all_phases, gas_phases = gas_phases, volumes=vols, melting_points=mps, experimentally_observed=obs)

    def __init__(self, phases,
                       volumes,
                       gas_phases = DEFAULT_GASES,
                       melting_points = None,
                       experimentally_observed = None):
        phases = phases + [SolidPhaseSet.FREE_SPACE]
        self.gas_phases = gas_phases
        self.volumes = volumes
        self.melting_points = melting_points
        self.experimentally_observed = experimentally_observed
        super().__init__(phases)

    def get_vol(self, phase):
        return self.volumes.get(phase)
    
    def get_melting_point(self, phase):
        return self.melting_points.get(phase)
    
    def is_theoretical(self, phase):
        return not self.experimentally_observed[phase]

    def get_theoretical_phases(self):
        return [phase for phase in self.phases if phase != SolidPhaseSet.FREE_SPACE and self.is_theoretical(phase)]
    
    def mole_amts_to_vols(self, mol_amts):
        return { phase: amt * self.get_vol(phase) for phase, amt in mol_amts.items() }

    def vol_amts_to_moles(self, vol_amts):
        { phase: amt / self.get_vol(phase) for phase, amt in vol_amts.items() }

    def mole_amts_to_el_amts(self, mol_amts):
        el_amts = {}

        for phase, p_amt in mol_amts.items():
            comp = Composition(phase)
            for el, el_amt in comp.get_el_amt_dict().items():
                if el in el_amts:
                    el_amts[el] += el_amt * p_amt
                else:
                    el_amts[el] = el_amt * p_amt
        
        return el_amts

    def vol_amts_to_el_amts(self, vol_amts):
        mol_amts = self.vol_amts_to_moles(vol_amts)
        return self.mole_amts_to_el_amts(mol_amts)

    def el_amts_to_el_fracs(self, el_amts):
        total = sum(el_amts.values())
        return {
            el: amt / total for el, amt in el_amts.items()
        }
        

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

def get_phase_vols(phases):

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

def get_experimentally_observed(phases):
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

def get_melting_points(phases):
    mp_json = pkg_resources.resource_filename("rxn_ca", "reactions/melting_points_df_08_08_23.json")
    melting_pt_data = pd.DataFrame(loadfn(mp_json))   # Note: temps in Kelvin
    mps = {}
    for p in phases:
        mps[p] = int(melting_pt_data[melting_pt_data["reduced_formula"] == p].iloc[0]["melting_point"])
    return mps