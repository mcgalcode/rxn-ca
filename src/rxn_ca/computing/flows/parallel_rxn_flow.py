from jobflow import Maker, job

from rxn_ca import setup, get_scored_rxns
from ...core.multi_stage_runner import run_multi
from ...core.recipe import ReactionRecipe

from ..schemas.ca_result_schema import RxnCAResultDoc

import multiprocessing as mp
from dataclasses import dataclass

_heating = "_heating"
_rxns = "rxns"
_recipe = "recipe"

def get_result(_):
    sim = setup(
        mp_globals[_rxns].phases,
        mp_globals[_recipe].reactant_amounts,
        mp_globals[_recipe].simulation_size,
        particle_size=mp_globals[_recipe].particle_size
    )
    try:
        return run_multi(
            sim,
            mp_globals[_rxns],
            mp_globals[_heating]
        )
    except Exception as e:
        print(e)
        print("Run failed!")
        return None

@dataclass
class ParallelRxnMaker(Maker):

    name: str = "Rxn Sim. Parallel"

    @job(data="results")
    def make(self,
             recipe: ReactionRecipe,
             output_fname: str = None):


        print("================= RETRIEVING AND SCORING REACTIONS =================")
        rxn_lib = get_scored_rxns(
            recipe.chem_sys,
            heating_sched=recipe.heating_schedule,
            exclude_phases=recipe.exclude_phases
        )

        print()
        print()
        print()

        # print("================= SETTING UP SIMULATION =================")





        print(f'================= RUNNING SIMULATION w/ {recipe.num_realizations} REALIZATIONS =================')


        global mp_globals

        mp_globals = {
            _heating: recipe.heating_schedule,
            _rxns: rxn_lib,
            _recipe: recipe
        }


        with mp.get_context("fork").Pool(recipe.num_realizations) as pool:
            results = pool.map(get_result, [_ for _ in range(recipe.num_realizations)])

        good_results = [res for res in results if res is not None]
        print(f'{len(good_results)} results achieved out of {len(results)}')
  
        result_doc = RxnCAResultDoc(
            recipe=recipe,
            results=good_results,
            reaction_library=rxn_lib
        )

        if output_fname is not None:
            print(f'================= SAVING RESULTS to {output_fname} =================')
            result_doc.to_file(output_fname)

        return result_doc

