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

    assert analyzer.get_absolute_molar_amt(rxn_one.state, "BaTiO3") == analyzer.get_absolute_molar_amt(rxn_two.state, "BaTiO3") * 2

