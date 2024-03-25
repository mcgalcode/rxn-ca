from pylattica.core import SimulationState
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY
from pymatgen.core.composition import Composition

from ..phases.solid_phase_set import SolidPhaseSet, MatterPhase
from ..core.constants import VOLUME, VOL_MULTIPLIER, GASES_EVOLVED
from ..utilities.helpers import normalize_dict

from typing import Union, List, Dict

from enum import Enum

class AnalysisQuantity(Enum):

    MOLES    = "MOLES"
    ATOMS    = "ATOMS"
    VOLUME   = "VOLUME"
    MASS     = "MASS"
    ELEMENTS = "ELEMENT"

class AnalysisMode(Enum):
    
    ABSOLUTE   = "ABSOLUTE"
    FRACTIONAL = "FRACTIONAL"

class ReactionStepAnalyzer():

    def __init__(self, phase_set: SolidPhaseSet) -> None:
        self.phase_set: SolidPhaseSet = phase_set

    def set_step_group(self, step_group: Union[List[SimulationState], SimulationState]):
        if not isinstance(step_group, list):
            step_group = [step_group]
        self.steps = step_group
        return self
    
    def get_value_general(self,
                          a_type: AnalysisQuantity,
                          mode: AnalysisMode = AnalysisMode.ABSOLUTE,
                          include_matter_phases: List[MatterPhase] = None,
                          temperature: int = None,
                          phase: str = None,
                          phases: List[str] = None):
        values = None
        match a_type.value:
            case AnalysisQuantity.MOLES.value:
                values = self.get_all_absolute_molar_amounts()
            case AnalysisQuantity.ATOMS.value:
                values = self.get_all_absolute_atomic_molar_amts()
            case AnalysisQuantity.VOLUME.value:
                values = self.get_all_absolute_phase_volumes()
            case AnalysisQuantity.MASS.value:
                values = self.get_all_absolute_phase_masses()
            case AnalysisQuantity.ELEMENTS.value:
                values = self.get_molar_elemental_composition()
    
        if mode.value == AnalysisMode.FRACTIONAL.value:
            values = normalize_dict(values)

        if include_matter_phases is not None:

            to_exclude = []
            for p in values.keys():
                if self.phase_set.get_matter_phase(p, temperature) not in include_matter_phases:
                    to_exclude.append(p)

            for e in to_exclude:
                del values[e]

        
        if phases is not None:
            return { k: v for k, v in values.items() if k in phases}
        elif phase is not None:
            return values.get(phase, 0)
        else:
            return values


    def get_all_absolute_phase_volumes(self):
        phase_amts = {}
        for step in self.steps:
            vol_multiplier = step.get_general_state().get(VOL_MULTIPLIER, 1.0)
            for site in step.all_site_states():
                phase = site[DISCRETE_OCCUPANCY]
                if phase != SolidPhaseSet.FREE_SPACE:
                    vol = site[VOLUME] * vol_multiplier
                    if phase in phase_amts:
                        phase_amts[phase] += vol
                    else:
                        phase_amts[phase] = vol

            gaseous = step.get_general_state().get(GASES_EVOLVED, {})

            for phase, vol in gaseous.items():
                if phase in phase_amts:
                    phase_amts[phase] += vol
                else:
                    phase_amts[phase] = vol
                    

        return phase_amts

    def get_all_absolute_phase_masses(self):
        vols = self.get_all_absolute_phase_volumes()
        masses = {}
        for p, v in vols.items():
            masses[p] = v * self.phase_set.get_density(p)
        return masses
    
    def get_all_mass_fractions(self):
        phase_masses = self.get_all_absolute_phase_masses()
        return normalize_dict(phase_masses)    

    def phases_present(self):
        return list(self.get_all_absolute_phase_volumes().keys())

    def get_absolute_phase_volume(self, phase: str):
        return self.get_all_absolute_phase_volumes().get(phase)
    
    def get_total_volume(self):
        return sum(self.get_all_absolute_phase_volumes().values())

    def get_simulation_side_length(self) -> int:
        num_sites = len(self.steps[0].all_site_states())
        return round(num_sites ** (1/3))
    
    def get_simulation_size(self) -> int:
        num_sites = len(self.steps[0].all_site_states())
        return num_sites

    def get_all_volume_fractions(self):
        vols = self.get_all_absolute_phase_volumes()
        return normalize_dict(vols)
    
    def get_phase_volume_fraction(self,
                                  phase: str,):
        return self.get_all_volume_fractions().get(phase)

    def get_all_absolute_molar_amounts(self):
        phase_abs_vols = self.get_all_absolute_phase_volumes()
        return self.phase_set.vol_amts_to_moles(phase_abs_vols)
    
    def get_all_absolute_atomic_molar_amts(self):
        phase_abs_vols = self.get_all_absolute_phase_volumes()
        mole_amts = self.phase_set.vol_amts_to_moles(phase_abs_vols)

        res = {}
        for p, amt in mole_amts.items():    
            comp = Composition(p)
            res[p] = amt * comp.num_atoms
        
        return res

    def get_absolute_molar_amt(self, phase: str):
        return self.get_all_absolute_molar_amounts().get(phase)

    def get_mole_fraction(self, phase: str):
        return self.get_all_mole_fractions().get(phase)

    def get_all_mole_fractions(self):
        phase_moles = self.get_all_absolute_molar_amounts()
        return normalize_dict(phase_moles)

    def get_molar_elemental_composition(self):
        molar_abs_amts = self.get_all_absolute_molar_amounts()
        return self.phase_set.mole_amts_to_el_amts(molar_abs_amts)
    
    def get_fractional_elemental_composition(self):
        ecomp = self.get_molar_elemental_composition()
        return normalize_dict(ecomp)