import pytest

from rxn_ca.analysis import ReactionStepAnalyzer
from rxn_ca.phases import SolidPhaseSet

from rxn_ca.core.constants import VOL_MULTIPLIER, VOLUME

from pylattica.core import SimulationState
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY

NA_CL = "NaCl"
LI2_O = "Li2O"

NUM_CELLS = 20

@pytest.fixture
def half_and_half_sim_step():
    cell_vol = 1.0

    state = SimulationState()
    state.set_general_state({
        VOL_MULTIPLIER: 1.0
    })

    for i in range(0, 10):
        state.set_site_state(i, {
            DISCRETE_OCCUPANCY: NA_CL,
            VOLUME: cell_vol
        })
    
    for i in range(10, 20):
        state.set_site_state(i, {
            DISCRETE_OCCUPANCY: LI2_O,
            VOLUME: cell_vol
        })

    return state

@pytest.fixture
def phase_set():
    return SolidPhaseSet(
        [NA_CL, LI2_O],
        volumes={
            NA_CL: 1.0,
            LI2_O: 2.0
        },
        melting_points = {
            NA_CL: 100,
            LI2_O: 100
        },
        experimentally_observed={
            NA_CL: True,
            LI2_O: True
        }
    )


def test_total_volume(phase_set, half_and_half_sim_step):
    analyzer = ReactionStepAnalyzer(phase_set)
    total_vol = analyzer.get_total_volume(half_and_half_sim_step)

    assert total_vol == NUM_CELLS * 1.0