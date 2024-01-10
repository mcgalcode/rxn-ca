import pytest

from rxn_ca.analysis import ReactionStepAnalyzer
from rxn_ca.phases import SolidPhaseSet

from rxn_ca.core.constants import VOL_MULTIPLIER, VOLUME

from pylattica.core import SimulationState
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY

import numpy as np

NA_CL = "NaCl"
LI2_O = "Li2O"

NUM_CELLS = 20
VOL_SCALE_1 = 0.5

@pytest.fixture
def sim_step_with_vol_scale():
    cell_vol = 1.0

    state = SimulationState()
    state.set_general_state({
        VOL_MULTIPLIER: VOL_SCALE_1
    })

    for i in range(0, 20):
        state.set_site_state(i, {
            DISCRETE_OCCUPANCY: NA_CL,
            VOLUME: cell_vol
        })

    return state

@pytest.fixture
def sim_step_with_different_cell_vol():
    cell_vol = 0.8

    state = SimulationState()
    state.set_general_state({
        VOL_MULTIPLIER: 1.0
    })

    for i in range(0, 20):
        state.set_site_state(i, {
            DISCRETE_OCCUPANCY: NA_CL,
            VOLUME: cell_vol
        })
        
    return state

@pytest.fixture
def sim_step_with_vol_scale_and_different_cell_vol():
    cell_vol = 0.8

    state = SimulationState()
    state.set_general_state({
        VOL_MULTIPLIER: VOL_SCALE_1
    })

    for i in range(0, 20):
        state.set_site_state(i, {
            DISCRETE_OCCUPANCY: NA_CL,
            VOLUME: cell_vol
        })
        
    return state 

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
def complex_half_and_half_sim_step():
    cell_vol = 0.8

    state = SimulationState()
    state.set_general_state({
        VOL_MULTIPLIER: VOL_SCALE_1
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


def test_total_volume_basic(phase_set, half_and_half_sim_step):
    analyzer = ReactionStepAnalyzer(phase_set)
    total_vol = analyzer.get_total_volume(half_and_half_sim_step)

    assert total_vol == 20

def test_total_volume_with_vol_scale(phase_set, sim_step_with_vol_scale):
    analyzer = ReactionStepAnalyzer(phase_set)
    total_vol = analyzer.get_total_volume(sim_step_with_vol_scale)

    assert total_vol == 20 * 0.5 # (10)    

def test_total_volume_with_different_cell_size(phase_set, sim_step_with_different_cell_vol):
    analyzer = ReactionStepAnalyzer(phase_set)
    total_vol = analyzer.get_total_volume(sim_step_with_different_cell_vol)

    assert np.isclose(total_vol, 20 * 0.8) # (16)

def test_total_volume_with_vol_scale_and_different_cell_size(phase_set, sim_step_with_vol_scale_and_different_cell_vol):
    analyzer = ReactionStepAnalyzer(phase_set)
    total_vol = analyzer.get_total_volume(sim_step_with_vol_scale_and_different_cell_vol)

    assert np.isclose(total_vol, 20 * 0.8 * 0.5) # (10)    

def test_phase_volumes(phase_set, half_and_half_sim_step):
    analyzer = ReactionStepAnalyzer(phase_set)
    phase_vols = analyzer.get_all_absolute_phase_volumes(half_and_half_sim_step)

    assert phase_vols[NA_CL] == 10
    assert phase_vols[LI2_O] == 10

def test_phase_volumes_complex(phase_set, complex_half_and_half_sim_step):
    analyzer = ReactionStepAnalyzer(phase_set)
    phase_vols = analyzer.get_all_absolute_phase_volumes(complex_half_and_half_sim_step)

    assert np.isclose(phase_vols[NA_CL], 4)
    assert np.isclose(phase_vols[LI2_O], 4)

def test_phase_volume_fractions(phase_set, half_and_half_sim_step):
    analyzer = ReactionStepAnalyzer(phase_set)
    phase_vol_fracs = analyzer.get_all_volume_fractions(half_and_half_sim_step)

    assert phase_vol_fracs[NA_CL] == 0.5
    assert phase_vol_fracs[LI2_O] == 0.5

def test_phase_volume_fractions_complex(phase_set, complex_half_and_half_sim_step):
    analyzer = ReactionStepAnalyzer(phase_set)
    phase_vol_fracs = analyzer.get_all_volume_fractions(complex_half_and_half_sim_step)

    assert phase_vol_fracs[NA_CL] == 0.5
    assert phase_vol_fracs[LI2_O] == 0.5

def test_molar_amounts(phase_set, half_and_half_sim_step):
    analyzer = ReactionStepAnalyzer(phase_set)
    molar_amts = analyzer.get_all_absolute_molar_amounts(half_and_half_sim_step)

    assert molar_amts[NA_CL] == 10
    assert molar_amts[LI2_O] == 5

def test_molar_amounts_complex(phase_set, complex_half_and_half_sim_step):
    analyzer = ReactionStepAnalyzer(phase_set)
    molar_amts = analyzer.get_all_absolute_molar_amounts(complex_half_and_half_sim_step)

    assert np.isclose(molar_amts[NA_CL], 4)
    assert np.isclose(molar_amts[LI2_O], 2)

def test_mole_fractions(phase_set, half_and_half_sim_step):
    analyzer = ReactionStepAnalyzer(phase_set)
    mole_fractions = analyzer.get_all_mole_fractions(half_and_half_sim_step)

    assert mole_fractions[NA_CL] == 2 / 3
    assert mole_fractions[LI2_O] == 1 / 3

def test_mole_fractions_complex(phase_set, complex_half_and_half_sim_step):
    analyzer = ReactionStepAnalyzer(phase_set)
    mole_fractions = analyzer.get_all_mole_fractions(complex_half_and_half_sim_step)

    assert np.isclose(mole_fractions[NA_CL], 2 / 3)
    assert np.isclose(mole_fractions[LI2_O], 1 / 3)

def test_molar_elemental_composition_complex(phase_set, complex_half_and_half_sim_step):
    analyzer = ReactionStepAnalyzer(phase_set)
    elemental_amts = analyzer.get_molar_elemental_composition(complex_half_and_half_sim_step)

    assert np.isclose(elemental_amts["Na"], 4)
    assert np.isclose(elemental_amts["Cl"], 4)
    assert np.isclose(elemental_amts["Li"], 4)
    assert np.isclose(elemental_amts["O"], 2)

def test_molar_elemental_fractions(phase_set, half_and_half_sim_step):
    analyzer = ReactionStepAnalyzer(phase_set)
    elemental_amts = analyzer.get_fractional_elemental_composition(half_and_half_sim_step)

    assert elemental_amts["Na"] == 2 / 7
    assert elemental_amts["Cl"] == 2 / 7
    assert elemental_amts["Li"] == 2 / 7
    assert elemental_amts["O"] == 1 / 7

def test_molar_elemental_fractions_complex(phase_set, complex_half_and_half_sim_step):
    analyzer = ReactionStepAnalyzer(phase_set)
    elemental_amts = analyzer.get_fractional_elemental_composition(complex_half_and_half_sim_step)

    assert elemental_amts["Na"] == 2 / 7
    assert elemental_amts["Cl"] == 2 / 7
    assert elemental_amts["Li"] == 2 / 7
    assert elemental_amts["O"] == 1 / 7