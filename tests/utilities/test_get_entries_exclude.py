"""get_entries' exclude_phases: the phases are dropped from the Materials Project
entries before the hull is built. MPRester and process_entries are stubbed."""

from pymatgen.entries.computed_entries import ComputedEntry

import rxn_ca.utilities.get_entries as get_entries_module


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
