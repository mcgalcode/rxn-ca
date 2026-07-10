import pytest

from rxn_ca.core.recipe import ReactionRecipe
from rxn_ca.reactions import ReactionLibrary
from rxn_ca.utilities.parallel_sim import run_single_sim, run_sim_parallel
from rxn_ca.analysis import BulkReactionAnalyzer, ReactionStepAnalyzer

import numpy as np

# The mean-field solid product/reactant pairing conserves mass only in
# expectation, so single-realization runs drift stochastically. Empirically
# (size 20, 14 trials): mean max-deviation 0.81%, observed max 1.39%; at size
# 25 an outlier of 1.84% was observed. The 2.5% threshold sits above that
# tail while remaining far below the drift produced by genuine conservation
# bugs.
DEVIATION_THRESHOLD = 0.025
SIMULATION_SIZE = 20

def test_basic_mass_conservation(get_test_file_path):
    recipe = ReactionRecipe.from_file(get_test_file_path("integration/batio3_recipe.json"))
    recipe.simulation_size = SIMULATION_SIZE
    rxn_lib = ReactionLibrary.from_file(get_test_file_path("integration/batio3_library.json"))

    result_doc = run_single_sim(recipe, reaction_lib=rxn_lib)

    result_analyzer = BulkReactionAnalyzer(results=result_doc.results,
                                           phase_set=result_doc.reaction_library.phases,
                                           heating_sched=result_doc.recipe.heating_schedule)

    step_analyzer = ReactionStepAnalyzer(rxn_lib.phases)

    step_analyzer.set_step_group(result_analyzer.get_first_steps())
    initial_el_comp = step_analyzer.get_molar_elemental_composition()

    for i in range(0,len(result_doc.results[0]), 100):
        steps = result_analyzer.get_steps(i)
        step_analyzer.set_step_group(steps)
        elemental_composition = step_analyzer.get_molar_elemental_composition()
        for el, amt in elemental_composition.items():
            initial_amt = initial_el_comp[el]
            fractional_deviation = np.abs(initial_amt - amt) / initial_amt
            assert fractional_deviation < DEVIATION_THRESHOLD, f'{el} has deviated by {fractional_deviation} by step {i}'

def test_parallel_mass_conservation(get_test_file_path):
    recipe = ReactionRecipe.from_file(get_test_file_path("integration/batio3_recipe.json"))
    rxn_lib = ReactionLibrary.from_file(get_test_file_path("integration/batio3_library.json"))

    result_doc = run_sim_parallel(recipe, reaction_lib=rxn_lib)

    result_analyzer = BulkReactionAnalyzer(results=result_doc.results,
                                           phase_set=result_doc.reaction_library.phases,
                                           heating_sched=result_doc.recipe.heating_schedule)

    step_analyzer = ReactionStepAnalyzer(rxn_lib.phases)

    step_analyzer.set_step_group(result_analyzer.get_first_steps())
    initial_el_comp = step_analyzer.get_molar_elemental_composition()

    for i in range(0,len(result_doc.results[0]), 100):
        steps = result_analyzer.get_steps(i)
        step_analyzer.set_step_group(steps)
        elemental_composition = step_analyzer.get_molar_elemental_composition()
        for el, amt in elemental_composition.items():
            initial_amt = initial_el_comp[el]
            fractional_deviation = np.abs(initial_amt - amt) / initial_amt
            assert fractional_deviation < DEVIATION_THRESHOLD, f'{el} has deviated by {fractional_deviation} by step {i}'