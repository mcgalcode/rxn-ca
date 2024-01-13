from pylattica.core import SimulationState
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY

from ..phases.solid_phase_set import SolidPhaseSet
from ..core.constants import VOLUME, MELTED_AMTS, VOL_MULTIPLIER, GASES_EVOLVED
from ..utilities.helpers import normalize_dict

from typing import Union, List, Dict

class ReactionStepAnalyzer():

    def __init__(self, phase_set: SolidPhaseSet) -> None:
        super().__init__()
        self.phase_set: SolidPhaseSet = phase_set

    def get_all_absolute_phase_volumes(self, step_group: Union[List[SimulationState], SimulationState],
                                       include_melted: bool = True,
                                       include_evolved_gases: bool = False):
        if not isinstance(step_group, list):
            step_group = [step_group]

        phase_amts = {}
        for step in step_group:
            vol_multiplier = step.get_general_state().get(VOL_MULTIPLIER, 1.0)
            for site in step.all_site_states():
                phase = site[DISCRETE_OCCUPANCY]
                if phase != SolidPhaseSet.FREE_SPACE:
                    vol = site[VOLUME] * vol_multiplier
                    if phase in phase_amts:
                        phase_amts[phase] += vol
                    else:
                        phase_amts[phase] = vol
            
            if include_melted:
                melted = step.get_general_state().get(MELTED_AMTS, {})

                for phase, vol in melted.items():
                    if phase in phase_amts:
                        phase_amts[phase] += vol
                    else:
                        phase_amts[phase] = vol

            if include_evolved_gases:
                melted = step.get_general_state().get(GASES_EVOLVED, {})

                for phase, vol in melted.items():
                    if phase in phase_amts:
                        phase_amts[phase] += vol
                    else:
                        phase_amts[phase] = vol
                    

        return phase_amts

    def phases_present(self,
                       step_group: Union[List[SimulationState], SimulationState],
                       include_melted: bool = True,
                       include_evolved_gases: bool = False):
        return list(self.get_all_absolute_phase_volumes(step_group, include_melted=include_melted, include_evolved_gases=include_evolved_gases).keys())

    def get_absolute_melted_volumes(self, step_group: Union[List[SimulationState], SimulationState], temp: int) -> Dict[str, float]:
        all_volume = self.get_all_absolute_phase_volumes(step_group, include_melted=True)
        filtered = { p: v for p, v in all_volume.items() if self.phase_set.is_melted(p, temp) }
        return filtered
    
    def get_absolute_solid_volumes(self, step_group: Union[List[SimulationState], SimulationState], temp: int) -> Dict[str, float]:
        all_volume = self.get_all_absolute_phase_volumes(step_group, include_melted=True)
        filtered = { p: v for p, v in all_volume.items() if not self.phase_set.is_melted(p, temp) }
        return filtered

    def get_absolute_phase_volume(self,
                                  step_group: Union[List[SimulationState], SimulationState],
                                  phase: str,
                                  include_melted = True,
                                  include_evolved_gases: bool = False):
        return self.get_all_absolute_phase_volumes(step_group, include_melted=include_melted, include_evolved_gases=include_evolved_gases).get(phase)
    
    def get_total_volume(self,
                         step_group: Union[List[SimulationState], SimulationState],
                         include_melted = True,
                         include_evolved_gases: bool = False):
        return sum(self.get_all_absolute_phase_volumes(step_group, include_melted=include_melted, include_evolved_gases=include_evolved_gases).values())

    def get_simulation_size(self, step: SimulationState) -> int:
        num_sites = len(step.all_site_states())
        return round(num_sites ** (1/3))

    def get_ideal_step_volume(self, step: SimulationState) -> float:
        return len(step.all_site_states()) * step.get_general_state().get(VOL_MULTIPLIER, 1.0)
    
    def get_total_solid_volume(self, step_group: Union[List[SimulationState], SimulationState], temp: int) -> float:
        solid_vols = self.get_absolute_solid_volumes(step_group, temp)
        return sum(solid_vols.values())

    def get_all_volume_fractions(self,
                                 step_group: Union[List[SimulationState], SimulationState],
                                 include_melted = True,
                                 include_evolved_gases: bool = False):
        vols = self.get_all_absolute_phase_volumes(step_group, include_melted=include_melted, include_evolved_gases=include_evolved_gases)
        return normalize_dict(vols)
    
    def get_phase_volume_fraction(self,
                                  step_group: Union[List[SimulationState], SimulationState],
                                  phase: str,
                                  include_melted = True,
                                  include_evolved_gases: bool = False):
        return self.get_all_volume_fractions(step_group, include_melted=include_melted, include_evolved_gases=include_evolved_gases).get(phase)

    def get_all_absolute_molar_amounts(self,
                                       step_group: Union[List[SimulationState], SimulationState],
                                       include_melted = True,
                                       include_evolved_gases: bool = False):
        phase_abs_vols = self.get_all_absolute_phase_volumes(step_group, include_melted=include_melted, include_evolved_gases=include_evolved_gases)
        return self.phase_set.vol_amts_to_moles(phase_abs_vols)

    def get_absolute_molar_amt(self,
                               step_group: Union[List[SimulationState], SimulationState],
                               phase: str,
                               include_melted = True,
                               include_evolved_gases: bool = False):
        return self.get_all_absolute_molar_amounts(step_group, include_melted=include_melted, include_evolved_gases=include_evolved_gases).get(phase)

    def get_mole_fraction(self,
                          step_group: Union[List[SimulationState], SimulationState],
                          phase: str,
                          include_melted = True,
                          include_evolved_gases: bool = False):
        return self.get_all_mole_fractions(step_group, include_melted=include_melted, include_evolved_gases=include_evolved_gases).get(phase)

    def get_all_mole_fractions(self,
                               step_group: Union[List[SimulationState], SimulationState],
                               include_melted = True,
                               include_evolved_gases: bool = False):
        phase_moles = self.get_all_absolute_molar_amounts(step_group, include_melted=include_melted, include_evolved_gases=include_evolved_gases)
        return normalize_dict(phase_moles)

    def get_molar_elemental_composition(self,
                                        step_group: Union[List[SimulationState], SimulationState],
                                        include_melted = True,
                                        include_evolved_gases: bool = True):
        molar_abs_amts = self.get_all_absolute_molar_amounts(step_group, include_melted=include_melted, include_evolved_gases=include_evolved_gases)
        return self.phase_set.mole_amts_to_el_amts(molar_abs_amts)
    
    def get_fractional_elemental_composition(self,
                                             step_group: Union[List[SimulationState], SimulationState],
                                             include_melted = True,
                                             include_evolved_gases: bool = True):
        ecomp = self.get_molar_elemental_composition(step_group, include_melted=include_melted, include_evolved_gases=include_evolved_gases)
        return normalize_dict(ecomp)