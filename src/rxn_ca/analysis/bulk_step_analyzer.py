
from pymatgen.core.composition import Composition

from pylattica.core import SimulationState
from pylattica.discrete import DiscreteStepAnalyzer
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY

from ..phases.solid_phase_set import SolidPhaseSet
from ..core.constants import VOLUME, MELTED_AMTS, VOL_MULTIPLIER

from typing import List

class BulkReactionStepAnalyzer(DiscreteStepAnalyzer):

    def __init__(self, phase_set: SolidPhaseSet) -> None:
        super().__init__()
        self.phase_set: SolidPhaseSet = phase_set

    def summary(self, step_group, phases = None):
        if phases is None:
            phases = self.phases_present(step_group)
        
        moles = self.molar_breakdown(step_group)
        vol = self.phase_volumes(step_group)

        for p in phases:
            print(f'{p} moles: ', moles[p])
            print(f'{p} vol: {vol[p]}')

        denom = min([moles[p] for p in phases])
        for p in phases:
            print(f'mole ratio of {p}: ', moles[p] / denom)

        for el, amt in self.elemental_composition(step_group).items():
            print(f'{el} moles: ', amt)


    def phase_volumes(self, step_group: List[SimulationState], include_melted=True):
        phase_amts = {}
        for step in step_group:
            vol_scale = step.get_general_state().get(VOL_MULTIPLIER, 1.0)
            for site in step.all_site_states():
                phase = site[DISCRETE_OCCUPANCY]
                if phase is not SolidPhaseSet.FREE_SPACE:
                    vol = site[VOLUME] * vol_scale
                    if phase in phase_amts:
                        phase_amts[phase] += vol
                    else:
                        phase_amts[phase] = vol
            
            if include_melted:
                for phase, vol in step.get_general_state().get(MELTED_AMTS, {}).items():
                    if phase in phase_amts:
                        phase_amts[phase] += vol
                    else:
                        phase_amts[phase] = vol

        
        return phase_amts
    
    def total_volume(self, step: SimulationState, include_melted=True):
        return sum(self.phase_volumes(step, include_melted=include_melted).values())

    def phase_volume_fractions(self, step_group: List[SimulationState], include_melted=True):
        total = self.total_volume(step_group[0], include_melted=include_melted)
        ratios = { phase: vol_abs / total for phase, vol_abs in self.phase_volumes(step_group, include_melted=include_melted).items()}
        return ratios

    def molar_breakdown(self, step_group: List[SimulationState], include_melted = True):
        phase_moles = {}
        for phase, vol in self.phase_volumes(step_group, include_melted=include_melted).items():
            if phase != SolidPhaseSet.FREE_SPACE:
                moles = vol / self.phase_set.volumes[phase]
                phase_moles[phase] = moles

        return phase_moles

    def molar_fractional_breakdown(self, step_group: List[SimulationState], include_melted = True):
        phase_moles = self.molar_breakdown(step_group, include_melted=include_melted)
        frac_moles = {}
        total = sum(phase_moles.values())
        for phase, amt in phase_moles.items():
            frac_moles[phase] = amt / total
        
        return frac_moles

    def elemental_composition(self, step_group: List[SimulationState], include_melted=True):
        molar_breakdown = self.molar_breakdown(step_group, include_melted=include_melted)
        elemental_amounts = {}
        total = 0
        for phase, moles in molar_breakdown.items():
            comp = Composition(phase)
            for el, am in comp.as_dict().items():
                num_moles = moles * am
                if el in elemental_amounts:
                    elemental_amounts[el] += num_moles
                else:
                    elemental_amounts[el] = num_moles
                total += num_moles

        for el, am in elemental_amounts.items():
            elemental_amounts[el] = am


        return elemental_amounts
    
    def elemental_composition_fractional(self, step_group: List[SimulationState], include_melted=True):
        molar_breakdown = self.molar_breakdown(step_group, include_melted=include_melted)
        elemental_amounts = {}
        total = 0
        for phase, moles in molar_breakdown.items():
            comp = Composition(phase)
            for el, am in comp.as_dict().items():
                num_moles = moles * am
                if el in elemental_amounts:
                    elemental_amounts[el] += num_moles
                else:
                    elemental_amounts[el] = num_moles
                total += num_moles

        for el, am in elemental_amounts.items():
            elemental_amounts[el] = am / total


        return elemental_amounts