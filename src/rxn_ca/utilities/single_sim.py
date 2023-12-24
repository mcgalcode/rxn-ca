from ..core.multi_stage_runner import run_multi
from ..core.recipe import ReactionRecipe
from ..reactions import ReactionLibrary

from ..computing.schemas.ca_result_schema import RxnCAResultDoc

from pylattica.core import Simulation
from rxn_network.reactions.reaction_set import ReactionSet

from ..phases import SolidPhaseSet

from .get_scored_rxns import get_scored_rxns
from .setup_reaction import setup_reaction


def run_single_sim(recipe: ReactionRecipe,
                   base_reactions: ReactionSet = None,
                   reaction_lib: ReactionLibrary = None,
                   initial_simulation: Simulation = None,
                   phase_set: SolidPhaseSet = None):

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

    if initial_simulation is None:

        print("================= SETTING UP SIMULATION =================")

        initial_simulation = setup_reaction(
            reaction_lib.phases,
            recipe.simulation_size,
            phase_mole_ratios = recipe.reactant_amounts,
        )

    print(f'================= RUNNING SIMULATION =================')


    result = run_multi(
        initial_simulation,
        reaction_lib,
        recipe.heating_schedule
    )

    result_doc = RxnCAResultDoc(
        recipe=recipe,
        results=[result],
        reaction_library=reaction_lib,
    )

    return result_doc

