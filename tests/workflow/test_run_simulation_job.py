"""Tests for the run_simulation job's analysis-only and library-path behavior.

These exercise the underlying function (run_simulation.original) directly so no
jobflow execution engine or Materials Project access is required -- a
ReactionLibraryData is built from the bundled test fixtures.
"""

import json
from pathlib import Path

import pytest

from rxn_ca.core.recipe import ReactionRecipe
from rxn_ca.reactions import ReactionLibrary
from rxn_ca.workflow.jobs import run_simulation
from rxn_ca.workflow.schemas import ReactionLibraryData, SimulationOutput


@pytest.fixture
def library_path(get_test_file_path):
    return get_test_file_path("integration/batio3_library.json")


@pytest.fixture
def recipe(get_test_file_path):
    return ReactionRecipe.from_file(get_test_file_path("integration/batio3_recipe.json"))


@pytest.fixture
def library_data(library_path):
    lib = ReactionLibrary.from_file(library_path)
    return ReactionLibraryData(
        phase_set_dict=lib.phases.as_dict(),
        reaction_library_dict=lib.as_dict(),
        chemical_system="Ba-Ti-O",
        temperatures=list(lib.temps),
        phases_available=list(lib.phases.phases),
        reaction_library_path=library_path,
    )


def _run(recipe, library_data, tmp_path, **kwargs):
    # Run inside tmp_path so the job's saved result_doc lands there.
    import os
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        return run_simulation.original(
            recipe=recipe,
            reaction_library_data=library_data,
            **kwargs,
        )
    finally:
        os.chdir(cwd)


def test_run_simulation_analysis_only(recipe, library_data, tmp_path):
    out = _run(recipe, library_data, tmp_path, analysis_only=True, num_frames=5)

    assert isinstance(out, SimulationOutput)
    assert out.reaction_library_path == library_data.reaction_library_path
    assert out.result_doc_path is not None
    # Analyzed trajectories are still produced from the observed volumes.
    assert len(out.final_molar_amounts) > 0
    assert len(out.step_indices) > 0
    # Temperatures come from the recorded state (no None gaps) and align with
    # the step indices.
    assert len(out.temperature_trajectory) == len(out.step_indices)
    assert all(t is not None for t in out.temperature_trajectory)


def test_run_simulation_result_doc_references_library_not_embeds(
    recipe, library_data, tmp_path
):
    out = _run(recipe, library_data, tmp_path, analysis_only=True, num_frames=5)

    # The persisted result doc references the library by path; the heavy library
    # dict is not embedded.
    with open(out.result_doc_path) as f:
        doc_dict = json.load(f)

    assert doc_dict["reaction_library"] is None
    assert doc_dict["reaction_library_path"] == library_data.reaction_library_path


def test_run_simulation_full_mode_also_references_library(
    recipe, library_data, tmp_path
):
    out = _run(recipe, library_data, tmp_path, num_frames=5)

    with open(out.result_doc_path) as f:
        doc_dict = json.load(f)

    # Even in full (non-analysis) mode the library is referenced, not embedded,
    # because setup provided a reaction_library_path.
    assert doc_dict["reaction_library"] is None
    assert doc_dict["reaction_library_path"] == library_data.reaction_library_path
    assert doc_dict["results"]  # full frames retained
