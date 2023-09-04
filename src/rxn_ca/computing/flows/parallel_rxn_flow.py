from jobflow import Maker, job

from rxn_ca import setup, get_scored_rxns
from ...core.multi_stage_runner import run_multi
from ...core.heating import HeatingSchedule

from ..schemas.ca_result_schema import RxnCAResultDoc

import multiprocessing as mp
from typing import List, Tuple, Union, Dict
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
             chem_sys: Union[str, List[str]],
             heating_schedule: HeatingSchedule,
             reactant_ratios: Dict[str, float],
             simulation_size: int = 21,
             num_realizations: int = 8):


        print("================= RETRIEVING AND SCORING REACTIONS =================")
        rxn_lib = get_scored_rxns(
            chem_sys,
            heating_sched=heating_schedule,
        )
        print()
        print()
        print()

        print("================= SETTING UP SIMULATION =================")


        sim = setup(
            rxn_lib.phases,
            reactant_ratios,
            simulation_size,
            particle_size=0.4
        )


        print(f'================= RUNNING SIMULATION w/ {num_realizations} REALIZATIONS =================')


        global mp_globals

        mp_globals = {
            _heating: heating_schedule,
            _rxns: rxn_lib
        }


        with mp.get_context("fork").Pool(num_realizations) as pool:
            results = pool.map(get_result, [sim for _ in range(num_realizations)])

        good_results = [res for res in results if res is not None]
        print(f'{len(good_results)} results achieved out of {len(results)}')

        jserializable = [res.as_dict() for res in results if res is not None]
        
        return RxnCAResultDoc(
            chem_sys=chem_sys,
            heating_schedule=heating_schedule,
            reactant_amts=reactant_ratios,
            results=jserializable
        )

