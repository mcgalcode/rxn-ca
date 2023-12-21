
from pymatgen.core.composition import Composition

from pylattica.core import SimulationState
from pylattica.discrete import DiscreteStepAnalyzer
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY

from ..phases.solid_phase_set import SolidPhaseSet
from ..core.constants import VOLUME, MELTED_AMTS, VOL_MULTIPLIER


class ReactionStepAnalyzer(DiscreteStepAnalyzer):

    def __init__(self, phase_set: SolidPhaseSet) -> None:
        super().__init__()
        self.phase_set: SolidPhaseSet = phase_set

    def phase_volumes(self, step: SimulationState, include_melted = True):
        phase_amts = {}
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
    
    def total_volume(self, step: SimulationState, include_melted = True):
        return sum(self.phase_volumes(step, include_melted=include_melted).values())

    def phase_volume_fractions(self, step: SimulationState, include_melted = True):
        total = self.total_volume(step, include_melted=include_melted)
        ratios = { phase: vol_abs / total for phase, vol_abs in self.phase_volumes(step, include_melted=include_melted).items()}
        return ratios

    def molar_breakdown(self, step, include_melted = True):
        phase_moles = {}
        for phase, vol in self.phase_volumes(step, include_melted=include_melted).items():
            if phase != SolidPhaseSet.FREE_SPACE:
                moles = vol / self.phase_set.volumes[phase]
                phase_moles[phase] = moles

        return phase_moles

    def molar_fractional_breakdown(self, step, include_melted = True):
        phase_moles = self.molar_breakdown(step, include_melted=include_melted)
        frac_moles = {}
        total = sum(phase_moles.values())
        for phase, amt in phase_moles.items():
            frac_moles[phase] = amt / total
        
        return frac_moles

    def elemental_composition(self, step, include_melted = True):
        molar_breakdown = self.molar_breakdown(step, include_melted=include_melted)
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
    
    def elemental_composition_fractional(self, step, include_melted = True):
        ecomp = self.elemental_composition(step, include_melted=include_melted)
        total = sum(ecomp.values())
        return { phase: val / total for phase, val in ecomp.items() }