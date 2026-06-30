import pytest

from rxn_ca.core.recipe import ReactionRecipe
from rxn_ca.reactions import ReactionLibrary
from rxn_ca.utilities.parallel_sim import run_single_sim
from rxn_ca.computing.schemas.ca_result_schema import RxnCAResultDoc
from rxn_ca.analysis import BulkReactionAnalyzer


@pytest.fixture
def recipe(get_test_file_path):
    return ReactionRecipe.from_file(get_test_file_path("integration/batio3_recipe.json"))


@pytest.fixture
def library_path(get_test_file_path):
    return get_test_file_path("integration/batio3_library.json")


def test_run_single_sim_references_library_path(recipe, library_path, tmp_path):
    result_doc = run_single_sim(recipe, reaction_lib_path=library_path, num_frames=5)

    # The path is recorded on the doc, and the library was loaded to run the sim.
    assert result_doc.reaction_library_path == library_path
    assert isinstance(result_doc.reaction_library, ReactionLibrary)

    # Serialized form omits the heavy embedded library...
    full_path = str(tmp_path / "with_lib.json")
    ref_path = str(tmp_path / "lib_by_path.json")
    result_doc.to_file(ref_path)

    # ...whereas embedding it (no path) is substantially larger.
    embedded = RxnCAResultDoc(
        recipe=result_doc.recipe,
        results=result_doc.results,
        reaction_library=result_doc.reaction_library,
        phases=result_doc.phases,
        final_simulation=result_doc.final_simulation,
    )
    embedded.to_file(full_path)

    import os
    assert os.path.getsize(ref_path) < os.path.getsize(full_path)

    # Reloaded path-ref doc still drives analysis (phases) and lazily resolves
    # the library when needed.
    reloaded = RxnCAResultDoc.from_file(ref_path)
    assert reloaded._reaction_library is None
    analyzer = BulkReactionAnalyzer.from_result_doc(reloaded)
    assert len(analyzer.loaded_step_idxs) > 0
    assert isinstance(reloaded.reaction_library, ReactionLibrary)
