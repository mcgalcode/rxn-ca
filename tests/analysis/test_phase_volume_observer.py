import numpy as np
import pytest

from rxn_ca.analysis import PhaseVolumeObserver, ReactionStepAnalyzer
from rxn_ca.analysis.reaction_step_analyzer import phase_volumes_for_state
from rxn_ca.phases import SolidPhaseSet
from rxn_ca.core.constants import GASES_EVOLVED, VOL_MULTIPLIER, VOLUME

from pylattica.core import SimulationState
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY

NA_CL = "NaCl"
LI2_O = "Li2O"
CO2 = "CO2"


@pytest.fixture
def phase_set():
    return SolidPhaseSet(
        [NA_CL, LI2_O],
        volumes={NA_CL: 1.0, LI2_O: 2.0},
        densities={NA_CL: 2.16, LI2_O: 2.01},
        melting_points={NA_CL: 100, LI2_O: 100},
        experimentally_observed={NA_CL: True, LI2_O: True},
    )


@pytest.fixture
def half_and_half_sim_step():
    state = SimulationState()
    state.set_general_state({VOL_MULTIPLIER: 1.0})
    for i in range(0, 10):
        state.set_site_state(i, {DISCRETE_OCCUPANCY: NA_CL, VOLUME: 1.0})
    for i in range(10, 20):
        state.set_site_state(i, {DISCRETE_OCCUPANCY: LI2_O, VOLUME: 1.0})
    return state


@pytest.fixture
def sim_step_with_gas():
    state = SimulationState()
    state.set_general_state({VOL_MULTIPLIER: 1.0, GASES_EVOLVED: {CO2: 3.0}})
    for i in range(0, 10):
        state.set_site_state(i, {DISCRETE_OCCUPANCY: NA_CL, VOLUME: 1.0})
    return state


def test_observer_matches_analyzer_single_state(phase_set, half_and_half_sim_step):
    observer = PhaseVolumeObserver()
    observed = observer.observe(half_and_half_sim_step, 0)

    analyzer = ReactionStepAnalyzer(phase_set)
    analyzer.set_step_group(half_and_half_sim_step)
    expected = analyzer.get_all_absolute_phase_volumes()

    assert observed == expected
    assert observed[NA_CL] == 10
    assert observed[LI2_O] == 10


def test_observer_includes_evolved_gases(sim_step_with_gas):
    observer = PhaseVolumeObserver()
    observed = observer.observe(sim_step_with_gas, 0)

    assert observed[NA_CL] == 10
    assert observed[CO2] == 3.0


def test_observer_delegates_to_shared_helper(half_and_half_sim_step):
    # The observer must be the same single source of truth as the analyzer.
    observer = PhaseVolumeObserver()
    assert observer.observe(half_and_half_sim_step, 0) == phase_volumes_for_state(
        half_and_half_sim_step
    )


def test_precomputed_volumes_derive_same_quantities(phase_set, half_and_half_sim_step):
    """Derived quantities from precomputed volumes match the full-state path."""
    full = ReactionStepAnalyzer(phase_set).set_step_group(half_and_half_sim_step)

    volumes = PhaseVolumeObserver().observe(half_and_half_sim_step, 0)
    reduced = ReactionStepAnalyzer(phase_set).set_precomputed_volumes(volumes)

    assert reduced.get_all_absolute_phase_volumes() == full.get_all_absolute_phase_volumes()
    assert reduced.get_all_absolute_molar_amounts() == full.get_all_absolute_molar_amounts()
    assert reduced.get_all_mass_fractions() == full.get_all_mass_fractions()
    assert reduced.get_all_mole_fractions() == full.get_all_mole_fractions()
    assert reduced.get_molar_elemental_composition() == full.get_molar_elemental_composition()
    assert np.isclose(reduced.get_total_volume(), full.get_total_volume())


def test_precomputed_volumes_sum_group(phase_set):
    """A list of volume dicts is summed, mirroring set_step_group."""
    analyzer = ReactionStepAnalyzer(phase_set)
    analyzer.set_precomputed_volumes([{NA_CL: 4}, {NA_CL: 6, LI2_O: 10}])

    vols = analyzer.get_all_absolute_phase_volumes()
    assert vols[NA_CL] == 10
    assert vols[LI2_O] == 10


def test_set_step_group_clears_precomputed(phase_set, half_and_half_sim_step):
    analyzer = ReactionStepAnalyzer(phase_set)
    analyzer.set_precomputed_volumes({NA_CL: 99})
    # Switching back to a full-state group must drop the precomputed volumes.
    analyzer.set_step_group(half_and_half_sim_step)
    assert analyzer._precomputed_volumes is None
    assert analyzer.get_all_absolute_phase_volumes()[NA_CL] == 10


def test_returned_volumes_are_copy(phase_set):
    analyzer = ReactionStepAnalyzer(phase_set)
    analyzer.set_precomputed_volumes({NA_CL: 5})
    vols = analyzer.get_all_absolute_phase_volumes()
    vols[NA_CL] = 0
    # Mutating the returned dict must not corrupt the stored volumes.
    assert analyzer.get_all_absolute_phase_volumes()[NA_CL] == 5
