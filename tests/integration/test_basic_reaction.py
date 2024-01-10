import pytest

from rxn_ca.core.recipe import ReactionRecipe
from rxn_ca.reactions import ReactionLibrary
from rxn_ca.utilities.parallel_sim import run_sim_parallel, run_single_sim
from rxn_ca.analysis import BulkReactionAnalyzer, ReactionStepAnalyzer

def test_basic_reaction(get_test_file_path):
    recipe = ReactionRecipe.from_file(get_test_file_path("integration/batio3_recipe.json"))
    rxn_lib = ReactionLibrary.from_file(get_test_file_path("integration/batio3_library.json"))

    result_doc = run_single_sim(recipe, reaction_lib=rxn_lib)

    result_analyzer = BulkReactionAnalyzer(results=result_doc.results,
                                           phase_set=result_doc.reaction_library.phases,
                                           heating_sched=result_doc.recipe.heating_schedule)

    step_analyzer = ReactionStepAnalyzer(rxn_lib.phases)

    product_phases = step_analyzer.phases_present(result_analyzer.get_final_steps())

    assert "BaTiO3" in product_phases
    assert "Ba2TiO4" in product_phases
    assert "BaTi2O5" in product_phases
    assert "Ba4Ti13O30" in product_phases

def test_parallel_rxn(get_test_file_path):
    recipe = ReactionRecipe.from_file(get_test_file_path("integration/batio3_recipe.json"))
    rxn_lib = ReactionLibrary.from_file(get_test_file_path("integration/batio3_library.json"))

    result_doc = run_sim_parallel(recipe, reaction_lib=rxn_lib)

    result_analyzer = BulkReactionAnalyzer(results=result_doc.results,
                                           phase_set=result_doc.reaction_library.phases,
                                           heating_sched=result_doc.recipe.heating_schedule)

    step_analyzer = ReactionStepAnalyzer(rxn_lib.phases)

    product_phases = step_analyzer.phases_present(result_analyzer.get_final_steps())

    assert "BaTiO3" in product_phases
    assert "Ba2TiO4" in product_phases
    assert "BaTi2O5" in product_phases
    assert "Ba4Ti13O30" in product_phases

    
