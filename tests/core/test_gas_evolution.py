import pytest

from pylattica.core import AsynchronousRunner
from pylattica.core.constants import GENERAL, SITES
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY

from rxn_ca.analysis.reaction_step_analyzer import ReactionStepAnalyzer
from rxn_ca.core.constants import GASES_EVOLVED
from rxn_ca.core.reaction_calculator import ReactionCalculator
from rxn_ca.core.reaction_controller import ReactionController
from rxn_ca.phases import SolidPhaseSet
from rxn_ca.reactions import ScoredReaction, ScoredReactionSet
from rxn_ca.utilities.setup_reaction import setup_reaction


@pytest.fixture(scope="module")
def carbonate_phases():
    # Molar volumes chosen so that one mole of BaCO3 (50 vol units) decomposes
    # into one mole of BaO (25 vol units) and one mole of CO2 (25 vol units)
    return SolidPhaseSet(
        ["BaCO3", "BaO", "CO2"],
        volumes={"BaCO3": 50.0, "BaO": 25.0, "CO2": 25.0},
        densities={"BaCO3": 4.29, "BaO": 5.72, "CO2": 1.98},
        melting_points={"BaCO3": 1600, "BaO": 2196, "CO2": 216},
        experimentally_observed={"BaCO3": True, "BaO": True, "CO2": True},
    )


@pytest.fixture(scope="module")
def decomposition_rxn():
    # BaCO3 -> BaO + CO2 expressed in volume stoichiometry
    return ScoredReaction(
        reactants={"BaCO3": 50.0},
        products={"BaO": 25.0, "CO2": 25.0},
        competitiveness=1.0,
    )


@pytest.fixture(scope="module")
def base_simulation(carbonate_phases):
    return setup_reaction(carbonate_phases, {"BaCO3": 1.0}, size=5)


def test_scored_reaction_gas_accounting(decomposition_rxn):
    assert decomposition_rxn.solid_products == frozenset(["BaO"])
    assert decomposition_rxn.gas_products == frozenset(["CO2"])

    # Half the product volume is solid
    assert decomposition_rxn.solid_product_fraction == pytest.approx(0.5)

    # Consuming one mole of BaCO3 (50 vol units) evolves one mole of CO2 (25 vol units)
    assert decomposition_rxn.gas_amt_from_reactant_vol("CO2", 50.0) == pytest.approx(25.0)
    # Gas evolution scales linearly with the consumed reactant volume
    assert decomposition_rxn.gas_amt_from_reactant_vol("CO2", 1.0) == pytest.approx(0.5)


def test_gas_is_never_chosen_as_cell_replacement(decomposition_rxn, carbonate_phases):
    rxn_set = ScoredReactionSet([decomposition_rxn], carbonate_phases)
    calc = ReactionCalculator(None, scored_rxns=rxn_set)

    for _ in range(50):
        assert calc.get_product_from_reaction(decomposition_rxn) == "BaO"

    # A reaction with no solid products yields no replacement phase
    all_gas_rxn = ScoredReaction(
        reactants={"BaCO3": 50.0},
        products={"CO2": 75.0},
        competitiveness=1.0,
    )
    assert calc.get_product_from_reaction(all_gas_rxn) is None


def test_decomposition_conserves_elements(carbonate_phases, decomposition_rxn, base_simulation):
    rxn_set = ScoredReactionSet([decomposition_rxn], carbonate_phases)

    nb_graph = ReactionController.get_neighborhood_from_structure(base_simulation.structure)
    calc = ReactionCalculator(nb_graph, scored_rxns=rxn_set)
    controller = ReactionController(base_simulation.structure, calc)

    analyzer = ReactionStepAnalyzer(carbonate_phases)
    initial_state = base_simulation.state.copy()
    initial_el_comp = analyzer.set_step_group(initial_state).get_molar_elemental_composition()

    runner = AsynchronousRunner()
    result = runner.run(initial_state, controller, num_steps=1000)
    final_state = result.last_step

    # The test is only meaningful if decomposition actually occurred
    evolved = final_state.get_general_state().get(GASES_EVOLVED)
    assert evolved.get("CO2", 0) > 0

    # Elemental amounts (including the evolved gas ledger) are conserved exactly,
    # not just in expectation
    final_el_comp = analyzer.set_step_group(final_state).get_molar_elemental_composition()
    for el, initial_amt in initial_el_comp.items():
        assert final_el_comp.get(el, 0) == pytest.approx(initial_amt, rel=1e-6), \
            f"{el} was not conserved: {initial_amt} -> {final_el_comp.get(el, 0)}"

    # The evolved CO2 exactly matches the amount of BaO produced (1:1 molar)
    moles = analyzer.set_step_group(final_state).get_all_absolute_molar_amounts()
    assert moles.get("CO2", 0) == pytest.approx(moles.get("BaO", 0), rel=1e-6)


def test_full_volatilization_leaves_free_space(carbonate_phases, base_simulation):
    # An (artificial) reaction whose products are entirely gaseous
    all_gas_rxn = ScoredReaction(
        reactants={"BaCO3": 50.0},
        products={"CO2": 75.0},
        competitiveness=1.0,
    )
    rxn_set = ScoredReactionSet([all_gas_rxn], carbonate_phases)

    nb_graph = ReactionController.get_neighborhood_from_structure(base_simulation.structure)
    calc = ReactionCalculator(nb_graph, scored_rxns=rxn_set)

    state = base_simulation.state.copy()

    update = None
    for _ in range(20):
        for site_id in state.site_ids():
            candidate = calc.get_state_update(site_id, state)
            if candidate.get(SITES):
                update = candidate
                break
        if update is not None:
            break

    assert update is not None, "No volatilization event occurred"

    # Every consumed cell becomes free space...
    consumed_sites = update[SITES]
    for site_update in consumed_sites.values():
        assert site_update[DISCRETE_OCCUPANCY] == SolidPhaseSet.FREE_SPACE

    # ...and each credits its stoichiometric share of gas: a cell of volume 1.0
    # evolves 1.0 * 75 / 50 = 1.5 vol units of CO2
    evolved = update[GENERAL][GASES_EVOLVED]
    assert evolved["CO2"] == pytest.approx(1.5 * len(consumed_sites))
