from .heating_schedule_runner import MeltAndRegrindMultiRunner, HeatingScheduleRunner
from ..core.recipe import ReactionRecipe
from ..reactions import ReactionLibrary

from ..computing.schemas.ca_result_schema import RxnCAResultDoc

from pylattica.core import Simulation
from rxn_network.reactions.reaction_set import ReactionSet

from ..phases import SolidPhaseSet
from ..core.reaction_controller import ReactionController
from ..core.liquid_swap_controller import LiquidSwapController
from ..core.reaction_calculator import ReactionCalculator

from .get_scored_rxns import get_scored_rxns
from .setup_reaction import setup_reaction, setup_noise_reaction


def run_single_sim(recipe: ReactionRecipe,
                   base_reactions: ReactionSet = None,
                   reaction_lib: ReactionLibrary = None,
                   initial_simulation: Simulation = None,
                   phase_set: SolidPhaseSet = None) -> RxnCAResultDoc:

    if base_reactions is None and reaction_lib is None:
        raise ValueError("Must provide either base_reactions or reaction_lib")

    if reaction_lib is None:

        print("================= RETRIEVING AND SCORING REACTIONS =================")

        reaction_lib: ReactionLibrary = get_scored_rxns(
            base_reactions,
            heating_sched=recipe.heating_schedule,
            exclude_phases=recipe.exclude_phases,
            exclude_theoretical=recipe.exclude_theoretical,
            scorer_class=recipe.get_score_class(),
            phase_set=phase_set
        )

    print()
    print()
    print()

    if len(recipe.exclude_phases) > 0:
        reaction_lib = reaction_lib.exclude_phases(recipe.exclude_phases)

    if recipe.exact_phase_set is not None:
        reaction_lib = reaction_lib.limit_phase_set(recipe.exact_phase_set)

    if initial_simulation is None:

        print("================= SETTING UP SIMULATION =================")

        initial_simulation = setup_noise_reaction(
            reaction_lib.phases,
            precursor_mole_ratios = recipe.reactant_amounts,
            size = recipe.simulation_size,
        )

    print(f'================= RUNNING SIMULATION =================')

    atmosphere = {}
    if recipe.atmospheric_phases is not None:
        for p in recipe.atmospheric_phases:
            atmosphere[p] = 1

    rxn_calculator = ReactionCalculator(
        LiquidSwapController.get_neighborhood_from_structure(initial_simulation.structure),
        open_species=atmosphere
    )

    controller = LiquidSwapController(
        initial_simulation.structure,
        rxn_calculator=rxn_calculator,
    )

    runner = HeatingScheduleRunner()

    result = runner.run_multi(
        initial_simulation,
        reaction_lib,
        recipe.heating_schedule,
        controller=controller
    )

    result_doc = RxnCAResultDoc(
        recipe=recipe,
        results=[result],
        reaction_library=reaction_lib,
    )

    return result_doc

