import pytest

from rxn_ca.phases import SolidPhaseSet
from rxn_ca.utilities.setup_reaction import setup_reaction
from rxn_ca.analysis.reaction_step_analyzer import ReactionStepAnalyzer

def test_volume_scale():

    SIZE = 15

    expected_mole_ratios = {
        "BaTiO3": 1,
        "TiO2": 1
    }

    phases = SolidPhaseSet.from_phase_list(list(expected_mole_ratios.keys()))
    analyzer = ReactionStepAnalyzer(phases)

    rxn_one = setup_reaction(phases, size=SIZE, phase_mole_ratios=expected_mole_ratios)
    rxn_two = setup_reaction(phases, size=SIZE, phase_mole_ratios=expected_mole_ratios, vol_multiplier=0.5)

    assert analyzer.mole_amt(rxn_one.state, "BaTiO3") == analyzer.mole_amt(rxn_two.state, "BaTiO3") * 2

