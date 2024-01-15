from ..core.recipe import ReactionRecipe

from ..reactions import ReactionLibrary
from ..phases import SolidPhaseSet
from ..computing.schemas.ca_result_schema import RxnCAResultDoc

from rxn_network.reactions.reaction_set import ReactionSet
from pylattica.core import Simulation

import multiprocessing as mp

from .single_sim import run_single_sim
from .get_scored_rxns import get_scored_rxns

_reaction_lib = "reaction_lib"
_recipe = "recipe"
_initial_simulation = "initial_simulation"

def _get_result(_):

    result: RxnCAResultDoc = run_single_sim(
        mp_globals[_recipe],
        reaction_lib=mp_globals.get(_reaction_lib),
        initial_simulation=mp_globals.get(_initial_simulation)
    )
    return result.results[0]


def run_sim_parallel(recipe: ReactionRecipe,
                     base_reactions: ReactionSet = None,
                     reaction_lib: ReactionLibrary = None,
                     initial_simulation: Simulation = None,
                     phase_set: SolidPhaseSet = None):

    print("================= RETRIEVING AND SCORING REACTIONS =================")

    if base_reactions is None and reaction_lib is None:
        raise ValueError("Must provide either base_reactions or reaction_lib")

    if reaction_lib is None:
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

    print(f'================= RUNNING SIMULATION w/ {recipe.num_realizations} REALIZATIONS =================')


    global mp_globals

    mp_globals = {
        _reaction_lib: reaction_lib,
        _recipe: recipe,
        _initial_simulation: initial_simulation
    }

    with mp.get_context("fork").Pool(recipe.num_realizations) as pool:
        results = pool.map(_get_result, [_ for _ in range(recipe.num_realizations)])

    good_results = [res for res in results if res is not None]
    print(f'{len(good_results)} results achieved out of {len(results)}')

    result_doc = RxnCAResultDoc(
        recipe=recipe,
        results=good_results,
        reaction_library=reaction_lib,
    )

    return result_doc

