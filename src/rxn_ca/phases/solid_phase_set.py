from pylattica.discrete.phase_set import PhaseSet
from rxn_network.reactions.reaction_set import ReactionSet
from mp_api.client import MPRester
import math
import pkg_resources
import pandas as pd
from monty.serialization import loadfn

class SolidPhaseSet(PhaseSet):

    FREE_SPACE = "Free Space"

    @classmethod
    def from_dict(cls, set_dict):
        return cls(
            set_dict["phases"],
            set_dict["volumes"],
            set_dict["melting_points"],
            set_dict["experimentally_observed"]
        )
    
    @classmethod
    def from_rxn_set(cls, rxn_set: ReactionSet):
        all_phases = list(set([e.composition.reduced_formula for e in rxn_set.entries]))
        vols = get_phase_vols(all_phases)
        mps = get_melting_points(all_phases)
        obs = get_experimentally_observed(all_phases)
        return cls(all_phases, volumes=vols, melting_points=mps, experimentally_observed=obs)

    def __init__(self, phases, volumes, melting_points = None, experimentally_observed = None):
        phases = phases + [SolidPhaseSet.FREE_SPACE]
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

    def as_dict(self):
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "phases": self.phases,
            "volumes": self.volumes,
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
        print(len(res))
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