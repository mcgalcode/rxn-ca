import pytest

from rxn_ca.analysis import BulkReactionAnalyzer, PhaseVolumeObserver
from rxn_ca.phases import SolidPhaseSet
from rxn_ca.core.constants import VOL_MULTIPLIER, VOLUME

from pylattica.core import ObservedResult, SimulationState
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY

NA_CL = "NaCl"
LI2_O = "Li2O"


@pytest.fixture
def phase_set():
    return SolidPhaseSet(
        [NA_CL, LI2_O],
        volumes={NA_CL: 1.0, LI2_O: 2.0},
        densities={NA_CL: 2.16, LI2_O: 2.01},
        melting_points={NA_CL: 100, LI2_O: 100},
        experimentally_observed={NA_CL: True, LI2_O: True},
    )


def _state(phase, n=10, vol=1.0):
    state = SimulationState()
    state.set_general_state({VOL_MULTIPLIER: 1.0})
    for i in range(n):
        state.set_site_state(i, {DISCRETE_OCCUPANCY: phase, VOLUME: vol})
    return state


def _observed_result(phase, num_steps=4):
    """An ObservedResult of per-phase volumes recorded every step."""
    obs = ObservedResult(observer=PhaseVolumeObserver(), compress_freq=1)
    state = _state(phase)
    obs.record(state, 0)
    for step in range(1, num_steps + 1):
        obs.record(state, step)
    obs.finalize(state, num_steps)
    return obs


def test_analyzer_builds_from_observed_results(phase_set):
    results = [_observed_result(NA_CL), _observed_result(LI2_O)]
    analyzer = BulkReactionAnalyzer(results, phase_set, None)

    assert analyzer._analysis_only is True
    # Observed every step from 0..4.
    assert analyzer.loaded_step_idxs == [0, 1, 2, 3, 4]
    assert len(analyzer.loaded_step_groups) == 5


def test_analyze_step_sums_across_results(phase_set):
    results = [_observed_result(NA_CL), _observed_result(LI2_O)]
    analyzer = BulkReactionAnalyzer(results, phase_set, None)

    vols = analyzer.analyze_step(2).get_all_absolute_phase_volumes()
    # Result A: 10 NaCl @ vol 1; Result B: 10 Li2O @ vol 1 -> summed.
    assert vols[NA_CL] == 10
    assert vols[LI2_O] == 10


def test_volume_trace_length(phase_set):
    results = [_observed_result(NA_CL), _observed_result(LI2_O)]
    analyzer = BulkReactionAnalyzer(results, phase_set, None)

    trace = analyzer.get_volume_trace()
    assert len(trace) == len(analyzer.loaded_step_idxs)
    # Total volume is constant here (20 each step).
    assert all(v == 20 for v in trace)


def test_final_molar_breakdown(phase_set):
    results = [_observed_result(NA_CL), _observed_result(LI2_O)]
    analyzer = BulkReactionAnalyzer(results, phase_set, None)

    breakdown = analyzer.get_final_molar_breakdown()
    # 10 vol NaCl -> 10 mol; 10 vol Li2O / 2.0 vol-per-mol -> 5 mol; fractions 2/3, 1/3.
    assert breakdown[NA_CL] == pytest.approx(2 / 3)
    assert breakdown[LI2_O] == pytest.approx(1 / 3)


def test_result_length_uses_total_steps(phase_set):
    results = [_observed_result(NA_CL, num_steps=4)]
    analyzer = BulkReactionAnalyzer(results, phase_set, None)
    # total_steps (4) + 1
    assert analyzer.result_length == 5


def test_get_step_size_unavailable(phase_set):
    results = [_observed_result(NA_CL)]
    analyzer = BulkReactionAnalyzer(results, phase_set, None)
    with pytest.raises(RuntimeError, match="analysis-only"):
        analyzer.get_step_size()
