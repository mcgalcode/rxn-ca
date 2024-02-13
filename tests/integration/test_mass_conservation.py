import pytest

from rxn_ca.core.recipe import ReactionRecipe
from rxn_ca.reactions import ReactionLibrary
from rxn_ca.utilities.parallel_sim import run_single_sim, run_sim_parallel
from rxn_ca.analysis import BulkReactionAnalyzer, ReactionStepAnalyzer

import numpy as np

def test_basic_mass_conservation(get_test_file_path):
    recipe = ReactionRecipe.from_file(get_test_file_path("integration/batio3_recipe.json"))
    rxn_lib = ReactionLibrary.from_file(get_test_file_path("integration/batio3_library.json"))

    result_doc = run_single_sim(recipe, reaction_lib=rxn_lib)

    result_analyzer = BulkReactionAnalyzer(results=result_doc.results,
                                           phase_set=result_doc.reaction_library.phases,
                                           heating_sched=result_doc.recipe.heating_schedule)

    step_analyzer = ReactionStepAnalyzer(rxn_lib.phases)

    initial_el_comp = step_analyzer.get_molar_elemental_composition(result_analyzer.get_first_steps())

    for i in range(0,len(result_doc.results[0]), 100):
        steps = result_analyzer.get_steps(i)
        elemental_composition = step_analyzer.get_molar_elemental_composition(steps)
        for el, amt in elemental_composition.items():
            initial_amt = initial_el_comp[el]
            fractional_deviation = np.abs(initial_amt - amt) / initial_amt
            assert fractional_deviation < 0.01, f'{el} has deviated by {fractional_deviation} by step {i}'

def test_parallel_mass_conservation(get_test_file_path):
    recipe = ReactionRecipe.from_file(get_test_file_path("integration/batio3_recipe.json"))
    rxn_lib = ReactionLibrary.from_file(get_test_file_path("integration/batio3_library.json"))

    result_doc = run_sim_parallel(recipe, reaction_lib=rxn_lib)

    result_analyzer = BulkReactionAnalyzer(results=result_doc.results,
                                           phase_set=result_doc.reaction_library.phases,
                                           heating_sched=result_doc.recipe.heating_schedule)

    step_analyzer = ReactionStepAnalyzer(rxn_lib.phases)

    initial_el_comp = step_analyzer.get_molar_elemental_composition(result_analyzer.get_first_steps())

    for i in range(0,len(result_doc.results[0]), 100):
        steps = result_analyzer.get_steps(i)
        elemental_composition = step_analyzer.get_molar_elemental_composition(steps)
        for el, amt in elemental_composition.items():
            initial_amt = initial_el_comp[el]
            fractional_deviation = np.abs(initial_amt - amt) / initial_amt
            assert fractional_deviation < 0.01, f'{el} has deviated by {fractional_deviation} by step {i}'