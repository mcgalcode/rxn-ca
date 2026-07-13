import pytest

from pylattica.core import AsynchronousRunner
from pylattica.core.constants import GENERAL, SITES
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY

from rxn_ca.analysis.reaction_step_analyzer import ReactionStepAnalyzer
from rxn_ca.core.constants import VOLUME, GASES_EVOLVED, RXN_PRODUCT_LEDGER
from rxn_ca.core.pairwise_reaction_calculator import PairwiseReactionCalculator
from rxn_ca.core.reaction_controller import ReactionController
from rxn_ca.phases import SolidPhaseSet
from rxn_ca.reactions import ScoredReaction, ScoredReactionSet
from rxn_ca.utilities.setup_reaction import setup_reaction


def _make_calculator(simulation, rxn_set, **kwargs):
    nb_graph = ReactionController.get_neighborhood_from_structure(simulation.structure)
    return PairwiseReactionCalculator(nb_graph, scored_rxns=rxn_set, **kwargs)


def _drive_until_event(calc, state, max_tries=50):
    """Repeatedly calls get_state_update on every site until an update that
    consumes at least one cell is produced (the proceed draw is stochastic)."""
    for _ in range(max_tries):
        for site_id in state.site_ids():
            update = calc.get_state_update(site_id, state)
            if update.get(SITES):
                return update
    raise AssertionError("No reaction event occurred")


# ---------------------------------------------------------------------------
# Exact pair consumption
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def titanate_phases():
    # Molar volumes chosen so BaO (25) + TiO2 (20) -> BaTiO3 (45) is volume
    # balanced
    return SolidPhaseSet(
        ["BaO", "TiO2", "BaTiO3"],
        volumes={"BaO": 25.0, "TiO2": 20.0, "BaTiO3": 45.0},
        densities={"BaO": 5.72, "TiO2": 4.23, "BaTiO3": 6.02},
        melting_points={"BaO": 2196, "TiO2": 2116, "BaTiO3": 1898},
        experimentally_observed={"BaO": True, "TiO2": True, "BaTiO3": True},
    )


@pytest.fixture(scope="module")
def titanate_rxn():
    # BaO + TiO2 -> BaTiO3 in volume stoichiometry
    return ScoredReaction(
        reactants={"BaO": 25.0, "TiO2": 20.0},
        products={"BaTiO3": 45.0},
        competitiveness=1.0,
    )


def test_exact_pair_consumption(titanate_phases, titanate_rxn):
    rxn_set = ScoredReactionSet([titanate_rxn], titanate_phases)
    simulation = setup_reaction(titanate_phases, {"BaO": 1.0, "TiO2": 1.0}, size=5)
    calc = _make_calculator(simulation, rxn_set)

    # All cells have volume 1.0, so every event involves one BaO cell and one
    # TiO2 cell: delta_v = min(1.0 * 45/25, 1.0 * 45/20) = 1.8, consuming
    # 1.8 * 25/45 = 1.0 of BaO (the limiting cell, freed) and
    # 1.8 * 20/45 = 0.8 of TiO2 (partial), producing 1.8 * 45/45 = 1.8 of
    # BaTiO3 - the consumed volumes are in the exact 25:20 stoichiometric
    # ratio.
    update = _drive_until_event(calc, simulation.state)

    site_updates = list(update[SITES].values())
    assert len(site_updates) == 2

    product_cells = [u for u in site_updates if u[DISCRETE_OCCUPANCY] == "BaTiO3"]
    partner_cells = [u for u in site_updates if u[DISCRETE_OCCUPANCY] == "TiO2"]
    assert len(product_cells) == 1
    assert len(partner_cells) == 1

    # Freed cell receives V_out = delta_v * S / R
    assert product_cells[0][VOLUME] == pytest.approx(1.8)
    # Partner keeps its phase with the stoichiometric remainder
    assert partner_cells[0][VOLUME] == pytest.approx(0.2)

    # Consumed volumes are exactly stoichiometric: BaO 1.0 : TiO2 0.8 = 25:20
    consumed_bao = 1.0
    consumed_tio2 = 1.0 - partner_cells[0][VOLUME]
    assert consumed_bao / consumed_tio2 == pytest.approx(25.0 / 20.0)

    # No gas is involved
    assert GASES_EVOLVED not in update[GENERAL]
    # The product ledger stays calculator-local unless store_ledger=True
    assert RXN_PRODUCT_LEDGER not in update[GENERAL]


# ---------------------------------------------------------------------------
# Deficit-based product selection
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def two_product_phases():
    return SolidPhaseSet(
        ["Ba2TiO4", "BaTiO3", "BaO"],
        volumes={"Ba2TiO4": 50.0, "BaTiO3": 30.0, "BaO": 20.0},
        densities={"Ba2TiO4": 5.5, "BaTiO3": 6.02, "BaO": 5.72},
        melting_points={"Ba2TiO4": 1860, "BaTiO3": 1898, "BaO": 2196},
        experimentally_observed={"Ba2TiO4": True, "BaTiO3": True, "BaO": True},
    )


@pytest.fixture(scope="module")
def two_product_rxn():
    # One volume unit of reactant yields 60% BaTiO3, 40% BaO
    return ScoredReaction(
        reactants={"Ba2TiO4": 1.0},
        products={"BaTiO3": 0.6, "BaO": 0.4},
        competitiveness=1.0,
    )


def test_deficit_product_sequence_is_largest_remainder(two_product_rxn):
    calc = PairwiseReactionCalculator(None, scored_rxns=None)

    # Unit-volume events against a 60:40 product split must alternate
    # (largest remainder): 0.6 vs 0.4 -> BaTiO3; 1.2-1 = 0.2 vs 0.8 -> BaO; ...
    ledger = {}
    sequence = [
        calc.choose_product_by_deficit(two_product_rxn, 1.0, ledger)
        for _ in range(4)
    ]
    assert sequence == ["BaTiO3", "BaO", "BaTiO3", "BaO"]

    # 2:2 over 4 unit events - never more than one event volume from exact
    rxn_ledger = ledger[two_product_rxn._as_str]
    assert rxn_ledger["BaTiO3"] == pytest.approx(2.0)
    assert rxn_ledger["BaO"] == pytest.approx(2.0)


def test_ledger_bounds_product_mix_error(two_product_phases, two_product_rxn):
    rxn_set = ScoredReactionSet([two_product_rxn], two_product_phases)
    simulation = setup_reaction(two_product_phases, {"Ba2TiO4": 1.0}, size=5)
    # store_ledger=True mirrors the calculator-local ledger into the general
    # state for inspection
    calc = _make_calculator(simulation, rxn_set, store_ledger=True)

    state = simulation.state.copy()

    num_events = 0
    while num_events < 30:
        update = _drive_until_event(calc, state)
        state.batch_update(update)
        num_events += 1

        # The mirrored ledger rides along in the general state updates, and
        # at every point the realized product mix is within one event volume
        # (1.0, a whole consumed cell) of exact 60:40 stoichiometry
        ledger = state.get_general_state().get(RXN_PRODUCT_LEDGER)
        assert ledger == calc._ledger
        rxn_ledger = ledger[two_product_rxn._as_str]
        produced_total = sum(rxn_ledger.values())
        for product, frac in [("BaTiO3", 0.6), ("BaO", 0.4)]:
            deviation = abs(rxn_ledger.get(product, 0) - produced_total * frac)
            assert deviation <= 1.0 + 1e-9, \
                f"{product} mix error {deviation} exceeds one event volume"


# ---------------------------------------------------------------------------
# Decomposition and volatilization
# ---------------------------------------------------------------------------

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
def carbonate_simulation(carbonate_phases):
    return setup_reaction(carbonate_phases, {"BaCO3": 1.0}, size=5)


def test_decomposition_frees_cells_and_credits_gas(carbonate_phases, decomposition_rxn, carbonate_simulation):
    rxn_set = ScoredReactionSet([decomposition_rxn], carbonate_phases)
    calc = _make_calculator(carbonate_simulation, rxn_set)

    update = _drive_until_event(calc, carbonate_simulation.state)

    # Single-reactant reactions fully drain the pooled cells (one cell for a
    # decomposition interaction, two for a same-phase pair interaction), so
    # every consumed cell is freed and becomes solid product
    consumed_vol = 0.0
    for site_update in update[SITES].values():
        assert site_update[DISCRETE_OCCUPANCY] == "BaO"
        # Each freed unit cell yields 1.0 * 25/50 = 0.5 of BaO
        assert site_update[VOLUME] == pytest.approx(0.5)
        consumed_vol += 1.0

    # Gas credit is exactly stoichiometric with the consumed volume
    evolved = update[GENERAL][GASES_EVOLVED]
    assert evolved["CO2"] == pytest.approx(0.5 * consumed_vol)


def test_full_volatilization_leaves_free_space(carbonate_phases, carbonate_simulation):
    # An (artificial) reaction whose products are entirely gaseous
    all_gas_rxn = ScoredReaction(
        reactants={"BaCO3": 50.0},
        products={"CO2": 75.0},
        competitiveness=1.0,
    )
    rxn_set = ScoredReactionSet([all_gas_rxn], carbonate_phases)
    calc = _make_calculator(carbonate_simulation, rxn_set)

    update = _drive_until_event(calc, carbonate_simulation.state)

    # Every consumed cell becomes free space...
    consumed_sites = update[SITES]
    for site_update in consumed_sites.values():
        assert site_update[DISCRETE_OCCUPANCY] == SolidPhaseSet.FREE_SPACE
        assert site_update[VOLUME] == pytest.approx(1.0)

    # ...and the event credits its stoichiometric share of gas: a cell of
    # volume 1.0 evolves 1.0 * 75 / 50 = 1.5 vol units of CO2
    evolved = update[GENERAL][GASES_EVOLVED]
    assert evolved["CO2"] == pytest.approx(1.5 * len(consumed_sites))


# ---------------------------------------------------------------------------
# Elemental conservation in a mini simulation
# ---------------------------------------------------------------------------

def test_decomposition_conserves_elements(carbonate_phases, decomposition_rxn, carbonate_simulation):
    rxn_set = ScoredReactionSet([decomposition_rxn], carbonate_phases)
    calc = _make_calculator(carbonate_simulation, rxn_set)
    controller = ReactionController(carbonate_simulation.structure, calc)

    analyzer = ReactionStepAnalyzer(carbonate_phases)
    initial_state = carbonate_simulation.state.copy()
    initial_el_comp = analyzer.set_step_group(initial_state).get_molar_elemental_composition()

    runner = AsynchronousRunner()
    result = runner.run(initial_state, controller, num_steps=1000)
    final_state = result.last_step

    # The test is only meaningful if decomposition actually occurred
    evolved = final_state.get_general_state().get(GASES_EVOLVED)
    assert evolved.get("CO2", 0) > 0

    # Elemental amounts (including the evolved gas ledger) are conserved
    # exactly, not just in expectation
    final_el_comp = analyzer.set_step_group(final_state).get_molar_elemental_composition()
    for el, initial_amt in initial_el_comp.items():
        assert final_el_comp.get(el, 0) == pytest.approx(initial_amt, rel=1e-6), \
            f"{el} was not conserved: {initial_amt} -> {final_el_comp.get(el, 0)}"

    # The evolved CO2 exactly matches the amount of BaO produced (1:1 molar)
    moles = analyzer.set_step_group(final_state).get_all_absolute_molar_amounts()
    assert moles.get("CO2", 0) == pytest.approx(moles.get("BaO", 0), rel=1e-6)
