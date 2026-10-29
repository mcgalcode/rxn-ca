"""Tests for the setup_reaction_library job's options and output.

_build_reaction_library is stubbed with the bundled BaTiO3 library, so no
Materials Project access or enumeration is needed; the job's own logic (where
the library is saved and what the output carries) runs as is.
"""

import os

import pytest
from pymatgen.core import Composition
from rxn_network.entries.entry_set import GibbsEntrySet
from rxn_network.entries.gibbs import GibbsComputedEntry

from rxn_ca.core.recipe import ReactionRecipe, ScoreTypes, _SCORE_TYPE_MAP
from rxn_ca.reactions import ReactionLibrary
from rxn_ca.workflow import jobs
from rxn_ca.workflow.jobs import run_simulation, setup_reaction_library


@pytest.fixture
def built(get_test_file_path, monkeypatch):
    """Stub the MP download, enumeration and scoring; record what they were asked for."""
    lib = ReactionLibrary.from_file(get_test_file_path("integration/batio3_library.json"))
    calls = []

    def fake_build(*args, **kwargs):
        calls.append((args, kwargs))
        return lib.phases, lib, {"n_entries": 7, "n_reactions": 42}

    monkeypatch.setattr(jobs, "_build_reaction_library", fake_build)
    monkeypatch.setattr(jobs, "_mp_database_version", lambda: "test-db")
    return lib, calls


def _in(tmp_path, fn, **kwargs):
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        return fn(**kwargs)
    finally:
        os.chdir(cwd)


def _setup(tmp_path, **kwargs):
    return _in(tmp_path, setup_reaction_library.original,
               chemical_system="Ba-Ti-O", temperatures=[300, 400], **kwargs)


def test_saved_library_is_referenced_not_embedded(built, tmp_path):
    target = tmp_path / "libs" / "batio3_lib.json"
    target.parent.mkdir()
    out = _setup(tmp_path, library_path=str(target))

    assert out.reaction_library_path == str(target.resolve())
    # the output stays small enough for a MongoDB job store: no library inline
    assert out.reaction_library_dict == {}
    assert out.phase_set_dict and out.phases_available
    assert ReactionLibrary.from_file(out.reaction_library_path).temps == built[0].temps
    assert not list(target.parent.glob("*.partial"))


def test_default_path_is_the_working_directory(built, tmp_path):
    out = _setup(tmp_path)

    assert out.reaction_library_path == str((tmp_path / "reaction_library_Ba_Ti_O.json").resolve())


def test_unsaved_library_is_embedded(built, tmp_path):
    out = _setup(tmp_path, save_to_file=False)

    assert out.reaction_library_path is None
    assert ReactionLibrary.from_dict(out.reaction_library_dict).temps == built[0].temps


def test_options_reach_the_build_and_the_metadata(built, tmp_path):
    out = _setup(tmp_path, thermo_types=["R2SCAN"], exclude_phases=["BaO2"],
                 scorer="GIBBS_ERF_SCORE", metastability_cutoff=0.03)

    args, kwargs = built[1][0]
    assert args[3] == 0.03
    assert kwargs == {"thermo_types": ["R2SCAN"], "exclude_phases": ["BaO2"],
                      "scorer": "GIBBS_ERF_SCORE"}
    assert out.metadata == {
        "n_entries": 7, "n_reactions": 42, "ensure_phases": [], "exclude_phases": ["BaO2"],
        "metastability_cutoff": 0.03, "exclude_theoretical": True, "thermo_types": ["R2SCAN"],
        "scorer": "GIBBS_ERF_SCORE", "mp_database_version": "test-db",
    }


def test_run_simulation_loads_the_path_only_output(built, tmp_path, get_test_file_path):
    lib_data = _setup(tmp_path, library_path=str(tmp_path / "lib.json"))
    recipe = ReactionRecipe.from_file(get_test_file_path("integration/batio3_recipe.json"))

    out = _in(tmp_path, run_simulation.original, recipe=recipe,
              reaction_library_data=lib_data, analysis_only=True, num_frames=5)

    assert out.reaction_library_path == lib_data.reaction_library_path
    assert len(out.final_molar_amounts) > 0


def _entry(formula, energy, entry_id):
    return GibbsComputedEntry(Composition(formula), formation_energy_per_atom=energy,
                              volume_per_atom=10.0, temperature=300, entry_id=entry_id)


def test_unknown_scorer_fails_before_any_download(monkeypatch):
    import rxn_ca.utilities.get_entries as get_entries_module
    monkeypatch.setattr(get_entries_module, "get_entries",
                        lambda **kw: pytest.fail("get_entries was called"))

    with pytest.raises(ValueError):
        jobs._build_reaction_library("Ba-Ti-O", [300], scorer="NOT_A_SCORER")


class _FakeReactionSet(list):
    def compute_at_temperatures(self, temps):
        return {t: self for t in temps}


def test_scorer_and_serial_fallback(monkeypatch):
    # stubbed around the scoring: the real enumerators start Ray and the phase
    # set queries Materials Project
    import rxn_ca.utilities.get_entries as get_entries_module
    import rxn_ca.utilities.get_scored_rxns as get_scored_rxns_module
    import rxn_network.enumerators.utils as enumerator_utils
    from rxn_ca.phases import SolidPhaseSet

    entries = GibbsEntrySet([_entry("Li", 0.0, "a"), _entry("O2", 0.0, "b"),
                             _entry("Li2O", -2.07, "c"), _entry("Li2O2", -1.6, "d")])
    seen = {}
    monkeypatch.setattr(get_entries_module, "get_entries",
                        lambda **kw: seen.update(kw) or entries)
    monkeypatch.setattr(enumerator_utils, "run_enumerators",
                        lambda enumerators, entries: _FakeReactionSet(["rxn1", "rxn2"]))
    monkeypatch.setattr(SolidPhaseSet, "from_entry_set", lambda *a, **kw: "phase set")
    calls = []

    def fake_score(rxn_set, parallel, **kw):
        calls.append((parallel, kw.get("scorer_class")))
        if parallel:
            raise RuntimeError("pool worker killed")
        return "scored library"

    monkeypatch.setattr(get_scored_rxns_module, "get_scored_rxns", fake_score)

    _, lib, stats = jobs._build_reaction_library(
        "Li-O", [300], thermo_types=["R2SCAN"], exclude_phases=["Li2O2"],
        scorer="GIBBS_ERF_SCORE")

    scorer_class = _SCORE_TYPE_MAP[ScoreTypes.GIBBS_ERF_SCORE]
    assert lib == "scored library"
    assert calls == [(True, scorer_class), (False, scorer_class)]
    assert seen["thermo_types"] == ["R2SCAN"] and seen["exclude_phases"] == ["Li2O2"]
    assert stats == {"n_entries": 4, "n_reactions": 2}
