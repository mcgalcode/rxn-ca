
from pymatgen.core.composition import Composition

from pylattica.core import SimulationState
from pylattica.discrete import DiscreteStepAnalyzer

from ..core.solid_phase_set import SolidPhaseSet
from ..reactions.scored_reaction_set import ScoredReactionSet


class ReactionStepAnalyzer(DiscreteStepAnalyzer):

    def __init__(self, reaction_set: ScoredReactionSet) -> None:
        super().__init__()
        self.rxn_set: ScoredReactionSet = reaction_set

    def summary(self, step, phases = None):
        if phases is None:
            phases = self.phases_present(step)

        for p in phases:
            print(f'{p} moles: ', self.cell_count(step, p))

        denom = min([self.cell_count(step, p) for p in phases])
        for p in phases:
            print(f'mole ratio of {p}: ', self.cell_count(step, p) / denom)

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

    def mole_fraction(self, step, phase):
        total_moles = 0
        for p in self.phases_present(step):
            if p is not SolidPhaseSet.FREE_SPACE:
                total_moles += self.moles_of(step, p)

        return self.moles_of(step, phase) 

    def molar_breakdown(self, step):
        mole_fractions = {}
        for phase in self.phases_present(step):
            mole_fractions[phase] = self.mole_fraction(step, phase)
        return mole_fractions

    def normalized_molar_breakdown(self, step):
        moles = self.molar_breakdown(step)
        min_moles = 1
        min_phase = None
        for phase, mole_count in moles.items():
            if mole_count <= min_moles:
                min_phase = phase
                min_moles = mole_count

        return { phase: mole_count / min_moles for phase, mole_count in moles.items() }


    def moles_of(self, step, phase):
        return float(self.cell_count(step, phase) / self.rxn_set.volumes[phase])

    def elemental_composition(self, step):
        phases = self.phases_present(step)
        elemental_amounts = {}
        total = 0
        for p in phases:
            if p is not SolidPhaseSet.FREE_SPACE:
                comp = Composition(p)
                moles = self.moles_of(step, p)
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