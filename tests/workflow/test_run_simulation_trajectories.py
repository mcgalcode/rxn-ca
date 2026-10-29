"""Tests for rxn_ca.workflow.jobs.run_simulation.

Focus: the trajectory-alignment fix (commit 6041bc8). run_simulation builds a
per-phase molar-amount trajectory. The earlier implementation appended
get_all_absolute_molar_amounts(step) per step, which omits phases absent at a
given step -- so each phase's trajectory was as long as "the number of steps it
happened to exist for", and the per-phase lists fell out of alignment on the
common step axis. A consumed precursor (present early, gone late) and a
late-forming product then looked identical, producing nonsense mass fractions
downstream.

The fix delegates to PhaseTraceCalculator, which walks every step and 0-fills
absent phases, so every phase trajectory has length == number of steps and is
aligned on the step axis. These tests lock that invariant in.

run_simulation is wrapped in jobflow's @job, so we call run_simulation.original
to execute the function body directly, mocking only the heavy simulation /
analysis machinery. The real PhaseTraceCalculator runs unmocked -- it is the
component under test.
"""

from types import SimpleNamespace
from unittest.mock import patch

import pytest

from rxn_ca.workflow.jobs import run_simulation


class _FakeStepAnalyzer:
    """Stands in for ReactionStepAnalyzer.

    PhaseTraceCalculator calls, per step group:
        step_analyzer.set_step_group(sg).get_value_general(quantity, mode,
                                                           include_matter_phases=...)
    We return a pre-canned analysis dict for each step group, ignoring the
    quantity/mode arguments (their handling is ReactionStepAnalyzer's concern,
    not the trajectory builder's).
    """

    def __init__(self, analyses_by_group):
        self._analyses_by_group = analyses_by_group
        self._current = None

    def set_step_group(self, step_group):
        self._current = step_group
        return self

    def get_value_general(self, *args, **kwargs):
        return self._analyses_by_group[self._current]


def _make_fake_analyzer(analyses_by_group, step_groups):
    """Build a fake BulkReactionAnalyzer exposing only what run_simulation uses.

    get_all_absolute_molar_amounts(idx) returns that step's ragged dict (only the
    phases present at that step) -- exactly what the buggy implementation appended
    per step. Wiring it this way makes these tests a real regression guard: revert
    to per-step appends and the per-phase lists come out unequal and unaligned,
    failing the assertions below.
    """
    step_analyzer = _FakeStepAnalyzer(analyses_by_group)
    return SimpleNamespace(
        loaded_step_idxs=list(range(len(step_groups))),
        loaded_step_groups=step_groups,
        step_analyzer=step_analyzer,
        # PhaseTraceCalculator's set_group_fn (BulkReactionAnalyzer.get_analyzer)
        get_analyzer=step_analyzer.set_step_group,
        last_loaded_step_idx=len(step_groups) - 1,
        get_all_absolute_molar_amounts=lambda idx: analyses_by_group[step_groups[idx]],
        # the recorded temperature run_simulation reads for each step
        get_temperature_at=lambda idx: 300.0 + 100.0 * idx,
    )


def _make_recipe(n_steps):
    """Minimal ReactionRecipe stand-in: as_dict, num_realizations, heating schedule."""
    heating_schedule = SimpleNamespace(
        all_temps=[300.0 + 100.0 * i for i in range(n_steps)],
        temp_at=lambda idx: 300.0 + 100.0 * idx,
    )
    return SimpleNamespace(
        as_dict=lambda: {},
        num_realizations=1,
        heating_schedule=heating_schedule,
    )


# Ragged per-step analyses: each step group reports only the phases present at
# that step. NaCl is consumed (present sg0, sg1; gone sg2); Na2O and Cl2 form.
RAGGED_ANALYSES = {
    "sg0": {"NaCl": 2.0},
    "sg1": {"NaCl": 1.0, "Na2O": 0.5},
    "sg2": {"Na2O": 1.0, "Cl2": 0.5},
}
STEP_GROUPS = ["sg0", "sg1", "sg2"]

FINAL_AMOUNTS = {"Na2O": 1.0, "Cl2": 0.5}


def _run(analyses_by_group, step_groups):
    """Invoke run_simulation's body with the heavy pieces mocked out."""
    n_steps = len(step_groups)
    fake_phase_set = SimpleNamespace(as_dict=lambda: {})
    fake_analyzer = _make_fake_analyzer(analyses_by_group, step_groups)

    with patch(
        "rxn_ca.workflow.jobs._build_reaction_library",
        return_value=(fake_phase_set, object(), {}),
    ), patch(
        "rxn_ca.utilities.single_sim.run_single_sim", return_value=object()
    ), patch("rxn_ca.analysis.BulkReactionAnalyzer") as mock_analyzer_cls:
        mock_analyzer_cls.from_result_doc.return_value = fake_analyzer
        return run_simulation.original(
            recipe=_make_recipe(n_steps),
            chemical_system="Na-Cl",
            save_to_file=False,
        )


def test_all_phase_trajectories_have_equal_length():
    """Every per-phase trajectory must be exactly len(step_indices) long."""
    output = _run(RAGGED_ANALYSES, STEP_GROUPS)

    n_steps = len(STEP_GROUPS)
    assert len(output.step_indices) == n_steps

    lengths = {
        phase: len(traj) for phase, traj in output.molar_amounts_trajectory.items()
    }
    assert lengths, "expected at least one phase trajectory"
    assert all(n == n_steps for n in lengths.values()), (
        f"trajectories not aligned to {n_steps} steps: {lengths}"
    )


@pytest.mark.parametrize(
    "phase, expected",
    [
        ("NaCl", [2.0, 1.0, 0.0]),
        ("Na2O", [0.0, 0.5, 1.0]),
        ("Cl2", [0.0, 0.0, 0.5]),
    ],
)
def test_absent_phases_are_zero_filled_and_aligned(phase, expected):
    """Absent phases are 0-filled in place, not dropped or truncated."""
    output = _run(RAGGED_ANALYSES, STEP_GROUPS)

    traj = output.molar_amounts_trajectory
    assert set(traj) == {"NaCl", "Na2O", "Cl2"}
    assert traj[phase] == pytest.approx(expected)


def test_temperature_trajectory_length_matches_steps():
    """Temperature trajectory shares the common step axis."""
    output = _run(RAGGED_ANALYSES, STEP_GROUPS)
    assert len(output.temperature_trajectory) == len(output.step_indices)


def test_final_molar_amounts_passed_through():
    """Sanity check that the mocked plumbing reaches the output unchanged."""
    output = _run(RAGGED_ANALYSES, STEP_GROUPS)
    assert output.final_molar_amounts == FINAL_AMOUNTS


def test_trajectories_aligned_when_every_step_differs():
    """Alignment must hold even when no two steps share a phase set."""
    analyses = {
        "a": {"A": 1.0},
        "b": {"B": 2.0},
        "c": {"C": 3.0},
        "d": {"A": 0.5, "C": 0.5},
    }
    step_groups = ["a", "b", "c", "d"]

    output = _run(analyses, step_groups)
    traj = output.molar_amounts_trajectory

    assert all(len(t) == 4 for t in traj.values())
    assert traj["A"] == pytest.approx([1.0, 0.0, 0.0, 0.5])
    assert traj["B"] == pytest.approx([0.0, 2.0, 0.0, 0.0])
    assert traj["C"] == pytest.approx([0.0, 0.0, 3.0, 0.5])
