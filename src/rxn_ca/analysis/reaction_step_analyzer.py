from pylattica.core import SimulationState
from pylattica.discrete import DiscreteStepAnalyzer
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY

from ..phases.solid_phase_set import SolidPhaseSet
from ..core.constants import VOLUME, MELTED_AMTS, VOL_MULTIPLIER
from ..utilities.helpers import normalize_dict

from typing import Union, List

class ReactionStepAnalyzer(DiscreteStepAnalyzer):

    def __init__(self, phase_set: SolidPhaseSet) -> None:
        super().__init__()
        self.phase_set: SolidPhaseSet = phase_set

    def get_all_absolute_phase_volumes(self, step_group: Union[List[SimulationState], SimulationState], include_melted: bool = True):
        if not isinstance(step_group, list):
            step_group = [step_group]

        phase_amts = {}
        for step in step_group:
            vol_multiplier = step.get_general_state().get(VOL_MULTIPLIER, 1.0)
            for site in step.all_site_states():
                phase = site[DISCRETE_OCCUPANCY]
                if phase is not SolidPhaseSet.FREE_SPACE:
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

        return phase_amts
    
    def get_total_volume(self, step_group: Union[List[SimulationState], SimulationState], include_melted = True):
        return sum(self.get_all_absolute_phase_volumes(step_group, include_melted=include_melted).values())

    def get_all_volume_fractions(self, step_group: Union[List[SimulationState], SimulationState], include_melted = True):
        vols = self.get_all_absolute_phase_volumes(step_group, include_melted=include_melted)
        return normalize_dict(vols)

    def get_all_absolute_molar_amounts(self, step_group: Union[List[SimulationState], SimulationState], include_melted = True):
        phase_abs_vols = self.get_all_absolute_phase_volumes(step_group, include_melted=include_melted)
        return self.phase_set.vol_amts_to_moles(phase_abs_vols)

    def get_mole_fraction(self, step_group: Union[List[SimulationState], SimulationState], phase, include_melted = True):
        return self.get_mole_fractions(step_group, include_melted=include_melted).get(phase)

    def get_all_mole_fractions(self, step_group: Union[List[SimulationState], SimulationState], include_melted = True):
        phase_moles = self.get_all_absolute_molar_amounts(step_group, include_melted=include_melted)
        return normalize_dict(phase_moles)

    def get_molar_elemental_composition(self, step_group: Union[List[SimulationState], SimulationState], include_melted = True):
        molar_abs_amts = self.get_all_absolute_molar_amounts(step_group, include_melted=include_melted)
        return self.phase_set.mole_amts_to_el_amts(molar_abs_amts)
    
    def get_fractional_elemental_composition(self, step_group: Union[List[SimulationState], SimulationState], include_melted = True):
        ecomp = self.get_molar_elemental_composition(step_group, include_melted=include_melted)
        return normalize_dict(ecomp)