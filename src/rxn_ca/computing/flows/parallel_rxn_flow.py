from jobflow import Maker, job

from rxn_ca import setup, get_scored_rxns
from ...core.multi_stage_runner import run_multi
from ...core.recipe import ReactionRecipe

from ..schemas.ca_result_schema import RxnCAResultDoc

import multiprocessing as mp
from dataclasses import dataclass

_heating = "_heating"
_rxns = "rxns"

def get_result(start):
    try:
        return run_multi(
            start,
            mp_globals[_rxns],
            mp_globals[_heating]
        )
    except:
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
        )

        print()
        print()
        print()

        print("================= SETTING UP SIMULATION =================")


        sim = setup(
            rxn_lib.phases,
            recipe.reactant_amounts,
            recipe.simulation_size,
            particle_size=recipe.particle_size
        )


        print(f'================= RUNNING SIMULATION w/ {recipe.num_realizations} REALIZATIONS =================')


        global mp_globals

        mp_globals = {
            _heating: recipe.heating_schedule,
            _rxns: rxn_lib
        }


        with mp.get_context("fork").Pool(recipe.num_realizations) as pool:
            results = pool.map(get_result, [sim for _ in range(recipe.num_realizations)])

        good_results = [res for res in results if res is not None]
        print(f'{len(good_results)} results achieved out of {len(results)}')

        jserializable = [res.as_dict() for res in results if res is not None]
        
        result_doc = RxnCAResultDoc(
            recipe=recipe,
            results=jserializable
        )

        if output_fname is not None:
            print(f'================= SAVING RESULTS to {output_fname} =================')
            result_doc.to_file(output_fname)

        return result_doc

