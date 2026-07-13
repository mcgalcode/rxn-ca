import numpy as np
import pytest

from rxn_ca.core.heating import HeatingSchedule, HeatingStep
from rxn_ca.core.recipe import ReactionRecipe
from rxn_ca.analysis import BulkReactionAnalyzer, ReactionStepAnalyzer
from rxn_ca.phases import SolidPhaseSet
from rxn_ca.reactions import ReactionLibrary, ScoredReaction, ScoredReactionSet
from rxn_ca.utilities.parallel_sim import run_single_sim

TEMP = 1100


@pytest.fixture(scope="module")
def carbonate_rxn_lib():
    # Molar volumes chosen so that one mole of BaCO3 (50 vol units) decomposes
    # into one mole of BaO (25 vol units) and one mole of CO2 (25 vol units)
    phases = SolidPhaseSet(
        ["BaCO3", "BaO", "CO2"],
        volumes={"BaCO3": 50.0, "BaO": 25.0, "CO2": 25.0},
        densities={"BaCO3": 4.29, "BaO": 5.72, "CO2": 1.98},
        melting_points={"BaCO3": 1600, "BaO": 2196, "CO2": 216},
        experimentally_observed={"BaCO3": True, "BaO": True, "CO2": True},
    )

    decomposition_rxn = ScoredReaction(
        reactants={"BaCO3": 50.0},
        products={"BaO": 25.0, "CO2": 25.0},
        competitiveness=1.0,
    )

    lib = ReactionLibrary(phases=phases)
    lib.add_rxns_at_temp(ScoredReactionSet([decomposition_rxn], phases), TEMP)
    return lib


@pytest.mark.parametrize("update_scheme", ["independent", "pairwise"])
def test_gas_evolution_mass_conservation(carbonate_rxn_lib, update_scheme):
    recipe = ReactionRecipe(
        heating_schedule=HeatingSchedule.build(
            HeatingStep.hold(TEMP, duration=2),
        ),
        reactant_amounts={"BaCO3": 1.0},
        simulation_size=8,
        num_realizations=1,
        update_scheme=update_scheme,
    )

    result_doc = run_single_sim(recipe, reaction_lib=carbonate_rxn_lib)

    result_analyzer = BulkReactionAnalyzer(results=result_doc.results,
                                           phase_set=result_doc.reaction_library.phases,
                                           heating_sched=result_doc.recipe.heating_schedule)

    step_analyzer = ReactionStepAnalyzer(carbonate_rxn_lib.phases)

    step_analyzer.set_step_group(result_analyzer.get_first_steps())
    initial_el_comp = step_analyzer.get_molar_elemental_composition()

    # Elemental amounts (grid solids + the evolved gas ledger) are conserved
    # at every checkpoint, including the final state
    num_steps = len(result_doc.results[0])
    checkpoints = list(range(0, num_steps, 100)) + [num_steps - 1]
    for i in checkpoints:
        step_analyzer.set_step_group(result_analyzer.get_steps(i))
        elemental_composition = step_analyzer.get_molar_elemental_composition()
        for el, initial_amt in initial_el_comp.items():
            amt = elemental_composition.get(el, 0)
            fractional_deviation = np.abs(initial_amt - amt) / initial_amt
            assert fractional_deviation < 1e-6, \
                f'{el} has deviated by {fractional_deviation} by step {i}'

    # The test is only meaningful if decomposition actually evolved gas, and
    # the evolved CO2 must match the produced BaO exactly (1:1 molar)
    final_moles = step_analyzer.get_all_absolute_molar_amounts()
    assert final_moles.get("CO2", 0) > 0, "No gas was evolved during the simulation"
    assert final_moles["CO2"] == pytest.approx(final_moles.get("BaO", 0), rel=1e-6)
