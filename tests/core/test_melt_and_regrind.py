import pytest

from rxn_ca.phases import SolidPhaseSet
from rxn_ca.core.melt_and_regrind import separate_solid_and_melt, calculate_melted_fraction, calculate_solid_ratio
from rxn_ca.core.constants import VOL_MULTIPLIER, VOLUME, MELTED_AMTS
from rxn_ca.analysis import ReactionStepAnalyzer
from rxn_ca.utilities.setup_reaction import setup_reaction

import numpy as np
from pylattica.core import SimulationState
from pylattica.core.constants import GENERAL
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY

@pytest.fixture
def phases():
    return SolidPhaseSet(
        ["NaCl", "Li2O", "YMnO3"],
        volumes={
            "NaCl": 1.0,
            "Li2O": 1.0,
            "YMnO3": 1.0,
        },
        melting_points={
            "NaCl": 600,
            "Li2O": 700,
            "YMnO3": 800,
        },
    )


def test_calculate_melted_fraction(phases):
    sim = setup_reaction(phases, {
        "NaCl": 1.0,
        "Li2O": 1.0,
        "YMnO3": 1.0,
    })

    assert np.isclose(calculate_melted_fraction(sim.state, phases, 500), 0, atol=0.001)
    assert np.isclose(calculate_melted_fraction(sim.state, phases, 650), 1 / 3, atol=0.001)
    assert np.isclose(calculate_melted_fraction(sim.state, phases, 750), 2 / 3, atol=0.001)
    assert np.isclose(calculate_melted_fraction(sim.state, phases, 850), 1.0, atol=0.001)

def test_calculate_melted_fraction_with_melted_phase(phases):
    sim = setup_reaction(phases, {
        "NaCl": 1.0,
        "Li2O": 1.0,
        "YMnO3": 1.0,
    })

    analyzer = ReactionStepAnalyzer(phases)
    initial_abs_vols = analyzer.get_all_absolute_phase_volumes(sim.state)

    # Add in an equivalent amount of melted NaCl
    sim.state.batch_update({
        GENERAL: {
            MELTED_AMTS: {
                "NaCl": initial_abs_vols.get("NaCl")
            }
        }
    })

    assert np.isclose(calculate_melted_fraction(sim.state, phases, 500), 0, atol=0.001)
    assert np.isclose(calculate_melted_fraction(sim.state, phases, 650), 1 / 3, atol=0.001)
    assert np.isclose(calculate_melted_fraction(sim.state, phases, 750), 2 / 3, atol=0.001)
    assert np.isclose(calculate_melted_fraction(sim.state, phases, 850), 1.0, atol=0.001)

def test_calculate_solid_ratio(phases):
    sim = setup_reaction(phases, {
        "NaCl": 1.0,
        "Li2O": 1.0,
        "YMnO3": 1.0,
    })

    analyzer = ReactionStepAnalyzer(phases)
    initial_abs_vols = analyzer.get_all_absolute_phase_volumes(sim.state)

    assert np.isclose(calculate_solid_ratio(sim.state, phases, 500), 1.0, atol=0.001)
    assert np.isclose(calculate_solid_ratio(sim.state, phases, 650), 2 / 3, atol=0.001)
    assert np.isclose(calculate_solid_ratio(sim.state, phases, 750), 1 / 3, atol=0.001)
    assert np.isclose(calculate_solid_ratio(sim.state, phases, 850), 0, atol=0.001)

    # Add in an equivalent amount of melted NaCl
    sim.state.batch_update({
        GENERAL: {
            MELTED_AMTS: {
                "NaCl": initial_abs_vols.get("NaCl")
            }
        }
    })

    assert np.isclose(calculate_solid_ratio(sim.state, phases, 500), 4 / 3, atol=0.001)
    assert np.isclose(calculate_solid_ratio(sim.state, phases, 650), 2 / 3, atol=0.001)
    assert np.isclose(calculate_solid_ratio(sim.state, phases, 750), 1 / 3, atol=0.001)
    assert np.isclose(calculate_solid_ratio(sim.state, phases, 850), 0, atol=0.001)
    

def test_separate_solid_and_melt(phases):
    sim = setup_reaction(phases, {
        "NaCl": 1.0,
        "Li2O": 1.0,
        "YMnO3": 1.0,
    })

    analyzer = ReactionStepAnalyzer(phases)
    initial_abs_vols = analyzer.get_all_absolute_phase_volumes(sim.state)

    step_one = separate_solid_and_melt(sim.state, phases, 650)
    assert np.isclose(step_one.get_general_state().get(VOL_MULTIPLIER), 2 /3, atol=0.001)
    vol_fracs_one = analyzer.get_all_volume_fractions(step_one, include_melted=False)
    assert vol_fracs_one.get("NaCl") is None
    assert np.isclose(vol_fracs_one.get("Li2O"), 0.5, atol=0.001)
    assert np.isclose(vol_fracs_one.get("YMnO3"), 0.5, atol=0.001)

    step_two = separate_solid_and_melt(step_one, phases, 750)
    assert np.isclose(step_two.get_general_state().get(VOL_MULTIPLIER), 1 /3, atol=0.001)
    vol_fracs_two = analyzer.get_all_volume_fractions(step_two, include_melted=False)
    assert vol_fracs_two.get("NaCl") is None
    assert vol_fracs_two.get("Li2O") is None
    assert np.isclose(vol_fracs_two.get("YMnO3"), 1.0)

    final_abs_vols = analyzer.get_all_absolute_phase_volumes(step_two)
    assert np.isclose(initial_abs_vols.get("NaCl"), final_abs_vols.get("NaCl"), rtol=0.001)
    assert np.isclose(initial_abs_vols.get("Li2O"), final_abs_vols.get("Li2O"), rtol=0.001)
    assert np.isclose(initial_abs_vols.get("YMnO3"), final_abs_vols.get("YMnO3"), rtol=0.001)

def test_resolidfy(phases):
    sim = setup_reaction(phases, {
        "NaCl": 1.0,
        "Li2O": 1.0,
        "YMnO3": 1.0,
    })

    analyzer = ReactionStepAnalyzer(phases)
    initial_abs_vols = analyzer.get_all_absolute_phase_volumes(sim.state)

    melted = separate_solid_and_melt(sim.state, phases, 750)
    assert np.isclose(melted.get_general_state().get(VOL_MULTIPLIER), 1 /3, atol=0.001)
    vol_fracs_melted = analyzer.get_all_volume_fractions(melted, include_melted=False)
    assert vol_fracs_melted.get("NaCl") is None
    assert vol_fracs_melted.get("Li2O") is None
    assert np.isclose(vol_fracs_melted.get("YMnO3"), 1.0)

    resolidified = separate_solid_and_melt(melted, phases, 300)
    assert np.isclose(resolidified.get_general_state().get(VOL_MULTIPLIER), 1.0, atol=0.001)
    vol_fracs_resolidifed = analyzer.get_all_volume_fractions(resolidified, include_melted=False)
    assert np.isclose(vol_fracs_resolidifed.get("NaCl"), 1/3, atol=0.001)
    assert np.isclose(vol_fracs_resolidifed.get("Li2O"), 1/3, atol=0.001)
    assert np.isclose(vol_fracs_resolidifed.get("YMnO3"), 1/3, atol=0.001)

    final_abs_vols = analyzer.get_all_absolute_phase_volumes(resolidified)
    assert np.isclose(initial_abs_vols.get("NaCl"), final_abs_vols.get("NaCl"), rtol=0.001)
    assert np.isclose(initial_abs_vols.get("Li2O"), final_abs_vols.get("Li2O"), rtol=0.001)
    assert np.isclose(initial_abs_vols.get("YMnO3"), final_abs_vols.get("YMnO3"), rtol=0.001)