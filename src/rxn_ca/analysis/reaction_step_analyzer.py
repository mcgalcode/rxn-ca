
from pymatgen.core.composition import Composition

from pylattica.core import SimulationState
from pylattica.discrete import DiscreteStepAnalyzer
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY

from ..core.solid_phase_set import SolidPhaseSet
from ..core.constants import VOLUME


class ReactionStepAnalyzer(DiscreteStepAnalyzer):

    def __init__(self, phase_set: SolidPhaseSet) -> None:
        super().__init__()
        self.phase_set: SolidPhaseSet = phase_set

    def summary(self, step, phases = None):
        if phases is None:
            phases = self.phases_present(step)
        
        moles = self.molar_fractional_breakdown(step)

        for p in phases:
            print(f'{p} moles: ', moles[p])

        denom = min([moles[p] for p in phases])
        for p in phases:
            print(f'mole ratio of {p}: ', moles[p] / denom)

        for el, amt in self.elemental_composition(step).items():
            print(f'{el} moles: ', amt)

    def get_reaction_choices(self, step: SimulationState):
        rxns = {}
        for site in step.all_site_states():
            rxn = site['rxn']
            rxn_str = str(rxn)
            if rxn_str in rxns:
                rxns[rxn_str] += 1
            else:
                rxns[rxn_str] = 0
        return rxns

    def phase_volumes(self, step: SimulationState):
        phase_amts = {}
        for site in step.all_site_states():
            phase = site[DISCRETE_OCCUPANCY]
            if phase is not SolidPhaseSet.FREE_SPACE:
                vol = site[VOLUME]
                if phase in phase_amts:
                    phase_amts[phase] += vol
                else:
                    phase_amts[phase] = vol
        
        return phase_amts
    
    def total_volume(self, step: SimulationState):
        return sum(self.phase_volumes(step).values())

    def phase_volume_fractions(self, step: SimulationState):
        total = self.total_volume(step)
        ratios = { phase: vol_abs / total for phase, vol_abs in self.phase_volumes(step).items()}
        return ratios

    def molar_breakdown(self, step):
        phase_moles = {}
        for phase, vol in self.phase_volumes(step).items():
            if phase != SolidPhaseSet.FREE_SPACE:
                moles = vol / self.phase_set.volumes[phase]
                phase_moles[phase] = moles

        return phase_moles

    def molar_fractional_breakdown(self, step):
        phase_moles = self.molar_breakdown(step)
        frac_moles = {}
        total = sum(phase_moles.values())
        for phase, amt in phase_moles.items():
            frac_moles[phase] = amt / total
        
        return frac_moles

    def elemental_composition(self, step):
        molar_breakdown = self.molar_breakdown(step)
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