import pytest

from rxn_ca.phases import SolidPhaseSet
from rxn_ca.utilities.setup_reaction import setup_reaction
from rxn_ca.analysis.reaction_step_analyzer import ReactionStepAnalyzer

from pylattica.core import Simulation

def test_volume_scale():
    expected_mole_ratios = {
        "BaTiO3": 1,
        "TiO2": 1
    }

    phases = SolidPhaseSet.from_phase_list(list(expected_mole_ratios.keys()))
    analyzer = ReactionStepAnalyzer(phases)

    rxn_one: Simulation = setup_reaction(phases, precursor_mole_ratios=expected_mole_ratios)
    rxn_two: Simulation = setup_reaction(phases, precursor_mole_ratios=expected_mole_ratios, vol_multiplier=0.5)

    analyzer.set_step_group(rxn_one.state)
    amt_one = analyzer.get_absolute_molar_amt("BaTiO3")
    analyzer.set_step_group(rxn_two.state)
    amt_two = analyzer.get_absolute_molar_amt("BaTiO3")
    assert amt_one == amt_two * 2

