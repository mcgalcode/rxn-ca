import pytest

from rxn_ca.core.recipe import ReactionRecipe
from rxn_ca.reactions import ReactionLibrary
from rxn_ca.utilities.parallel_sim import run_sim_parallel, run_single_sim
from rxn_ca.analysis import BulkReactionAnalyzer
from rxn_ca.computing.schemas.ca_result_schema import RxnCAResultDoc

from pylattica.core import ObservedResult

EXPECTED_PRODUCTS = {"BaTiO3", "Ba2TiO4", "BaTi2O5", "Ba4Ti13O30"}


@pytest.fixture
def recipe(get_test_file_path):
    return ReactionRecipe.from_file(get_test_file_path("integration/batio3_recipe.json"))


@pytest.fixture
def rxn_lib(get_test_file_path):
    return ReactionLibrary.from_file(get_test_file_path("integration/batio3_library.json"))


def test_single_sim_analysis_only_routes_to_observed_results(recipe, rxn_lib):
    result_doc = run_single_sim(
        recipe, reaction_lib=rxn_lib, analysis_only=True, num_frames=5
    )

    # Analysis-only mode discards full ReactionResults and keeps only the
    # lightweight per-phase-volume observations.
    assert result_doc.results == []
    assert len(result_doc.observed_results) == 1
    assert isinstance(result_doc.observed_results[0], ObservedResult)
    # The final live state is still available for chaining / inspection.
    assert result_doc.final_simulation is not None


def test_single_sim_analysis_only_analyzer_recovers_products(recipe, rxn_lib):
    result_doc = run_single_sim(
        recipe, reaction_lib=rxn_lib, analysis_only=True, num_frames=5
    )

    analyzer = BulkReactionAnalyzer.from_result_doc(result_doc)
    assert analyzer._analysis_only is True

    breakdown = analyzer.get_final_molar_breakdown()
    present = {phase for phase, amt in breakdown.items() if amt > 0}
    assert EXPECTED_PRODUCTS & present


def test_analysis_only_serialization_roundtrip(recipe, rxn_lib, tmp_path):
    result_doc = run_single_sim(
        recipe, reaction_lib=rxn_lib, analysis_only=True, num_frames=5
    )

    fname = str(tmp_path / "analysis_only_doc.json")
    result_doc.to_file(fname)

    reloaded = RxnCAResultDoc.from_file(fname)
    assert reloaded.results == []
    assert len(reloaded.observed_results) == 1

    # The reloaded doc drives a working analyzer.
    analyzer = BulkReactionAnalyzer.from_result_doc(reloaded)
    assert analyzer._analysis_only is True
    assert len(analyzer.get_volume_trace()) == len(analyzer.loaded_step_idxs)


def test_analysis_only_retains_temperature_from_state(recipe, rxn_lib):
    result_doc = run_single_sim(
        recipe, reaction_lib=rxn_lib, analysis_only=True, num_frames=5
    )

    analyzer = BulkReactionAnalyzer.from_result_doc(result_doc)

    schedule_temps = set(recipe.heating_schedule.all_temps)
    temps = [analyzer.get_temperature_at(i) for i in analyzer.loaded_step_idxs]

    # Every observed frame carries the temperature captured from the state, and
    # each value is one the heating schedule actually visited.
    assert all(t is not None for t in temps)
    assert set(temps) <= schedule_temps


def test_parallel_sim_analysis_only(recipe, rxn_lib):
    result_doc = run_sim_parallel(
        recipe, reaction_lib=rxn_lib, analysis_only=True, num_frames=5
    )

    assert result_doc.results == []
    assert len(result_doc.observed_results) == recipe.num_realizations
    assert all(isinstance(r, ObservedResult) for r in result_doc.observed_results)

    analyzer = BulkReactionAnalyzer.from_result_doc(result_doc)
    breakdown = analyzer.get_final_molar_breakdown()
    present = {phase for phase, amt in breakdown.items() if amt > 0}
    assert EXPECTED_PRODUCTS & present
