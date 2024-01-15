import pytest

from rxn_ca.phases import SolidPhaseSet
from rxn_ca.utilities.setup_reaction import setup_reaction
from rxn_ca.analysis.reaction_step_analyzer import ReactionStepAnalyzer
from rxn_ca.core.constants import GASES_CONSUMED, GASES_EVOLVED, MELTED_AMTS, VOL_MULTIPLIER

from pylattica.core import Simulation

NA_CL = "NaCl"
LI2_O = "Li2O"

@pytest.fixture
def basic_phase_set():
    return SolidPhaseSet(
        [NA_CL, LI2_O],
        volumes={
            NA_CL: 2.0,
            LI2_O: 0.5
        },
        melting_points={
            NA_CL: 800,
            LI2_O: 1000
        },
        experimentally_observed={
            NA_CL: True,
            LI2_O: False
        })

def test_general_state(basic_phase_set):
    expected_mole_ratios = {
        NA_CL: 1,
        LI2_O: 1
    }

    analyzer = ReactionStepAnalyzer(basic_phase_set)

    rxn_one: Simulation = setup_reaction(basic_phase_set, precursor_mole_ratios=expected_mole_ratios)
    
    state = rxn_one.state.get_general_state()

    assert GASES_EVOLVED in state
    assert GASES_CONSUMED in state
    assert MELTED_AMTS in state
    assert VOL_MULTIPLIER in state
    
