"""Tests for the get_entries utility function."""

import pytest
from pymatgen.entries.computed_entries import ComputedEntry

import rxn_ca.utilities.get_entries as get_entries_module
from rxn_ca.utilities.get_entries import get_entries


@pytest.mark.integration
def test_get_entries_basic():
    """Test basic get_entries functionality."""
    entries = get_entries(
        chem_sys="Li-O",
        metastability_cutoff=0.1,
        exclude_theoretical_phases=False,
    )

    formulas = {e.composition.reduced_formula for e in entries}
    assert "Li2O" in formulas
    assert "O2" in formulas


@pytest.mark.integration
def test_get_entries_ensure_phases_overrides_theoretical_filter():
    """Test that ensure_phases keeps phases even if marked theoretical in MP.

    This is a regression test for a bug where ensure_phases would not
    actually ensure inclusion of phases that were filtered out by
    exclude_theoretical_phases=True.

    Ca3Co4O9 is a real experimentally-synthesized thermoelectric material,
    but it's marked as 'theoretical' in Materials Project. When we have
    experimental evidence that a phase exists (e.g., from synthesis records),
    ensure_phases should override the theoretical filter.
    """
    # Ca3Co4O9 is marked as theoretical in MP but is experimentally real
    ensure_phases = ["Ca3Co4O9", "CaCO3", "Co3O4"]

    entries = get_entries(
        chem_sys="C-Ca-Co-O",
        metastability_cutoff=0.1,
        ensure_phases=ensure_phases,
        exclude_theoretical_phases=True,
    )

    formulas = {e.composition.reduced_formula for e in entries}

    # All ensure_phases should be present, even the "theoretical" one
    for phase in ensure_phases:
        assert phase in formulas, (
            f"{phase} should be included via ensure_phases "
            f"even if marked theoretical in MP"
        )


@pytest.mark.integration
def test_get_entries_ensure_phases_without_theoretical_filter():
    """Test ensure_phases works when theoretical filter is off."""
    ensure_phases = ["Ca3Co4O9"]

    entries = get_entries(
        chem_sys="C-Ca-Co-O",
        metastability_cutoff=0.1,
        ensure_phases=ensure_phases,
        exclude_theoretical_phases=False,
    )

    formulas = {e.composition.reduced_formula for e in entries}
    assert "Ca3Co4O9" in formulas


@pytest.mark.integration
def test_get_entries_theoretical_phases_excluded_by_default():
    """Test that theoretical phases are excluded when not in ensure_phases."""
    # Don't include Ca3Co4O9 in ensure_phases
    entries = get_entries(
        chem_sys="C-Ca-Co-O",
        metastability_cutoff=0.1,
        ensure_phases=[],
        exclude_theoretical_phases=True,
    )

    formulas = {e.composition.reduced_formula for e in entries}

    # Ca3Co4O9 should NOT be present since it's theoretical and not ensured
    assert "Ca3Co4O9" not in formulas, (
        "Ca3Co4O9 should be excluded as theoretical when not in ensure_phases"
    )


# --- exclude_phases: dropped from the Materials Project entries before the
# hull is built. MPRester and process_entries are stubbed.

class _FakeMPRester:
    def __init__(self, entries):
        self.entries = entries

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def get_entries_in_chemsys(self, elements, additional_criteria):
        return self.entries


class _Processed:
    def filter_by_stability(self, cutoff):
        return "stable entries"


def _run(monkeypatch, entries, **kwargs):
    seen = {}

    def fake_process(all_entries, *args, **kw):
        seen.update(entries=all_entries, **kw)
        return _Processed()

    monkeypatch.setattr(get_entries_module, "MPRester", lambda: _FakeMPRester(entries))
    monkeypatch.setattr(get_entries_module, "process_entries", fake_process)
    out = get_entries_module.get_entries("Li-O", **kwargs)
    return out, seen


def test_excluded_phases_never_reach_the_hull(monkeypatch):
    entries = [ComputedEntry("Li2O", -14.0, entry_id="a"), ComputedEntry("Li2O2", -18.0, entry_id="b"),
               ComputedEntry("Li2O2", -17.9, entry_id="c"), ComputedEntry("Li", -1.9, entry_id="d")]

    out, seen = _run(monkeypatch, entries, exclude_phases=["LiO"])  # LiO reduces like Li2O2

    assert out == "stable entries"
    assert sorted(e.entry_id for e in seen["entries"]) == ["a", "d"]  # every Li2O2 entry gone


def test_exclusion_wins_over_ensure_phases(monkeypatch):
    entries = [ComputedEntry("Li2O", -14.0, entry_id="a"), ComputedEntry("Li2O2", -18.0, entry_id="b")]

    _, seen = _run(monkeypatch, entries, exclude_phases=["Li2O2"], ensure_phases=["Li2O2", "Li2O"])

    assert seen["formulas_to_include"] == ["Li2O"]


def test_no_exclusion_keeps_every_entry(monkeypatch):
    entries = [ComputedEntry("Li2O", -14.0, entry_id="a"), ComputedEntry("Li2O2", -18.0, entry_id="b")]

    _, seen = _run(monkeypatch, entries)

    assert sorted(e.entry_id for e in seen["entries"]) == ["a", "b"]


def test_custom_entries_without_structures_fail_before_the_download(monkeypatch):
    monkeypatch.setattr(get_entries_module, "MPRester", lambda: pytest.fail("MP was queried"))

    with pytest.raises(TypeError, match="ComputedStructureEntry"):
        get_entries("Li-O", custom_entries=[ComputedEntry("Li2O", -14.0)])
