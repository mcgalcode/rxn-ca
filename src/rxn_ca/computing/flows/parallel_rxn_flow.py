from jobflow import Maker, job

from rxn_ca import setup, get_raw_rxns
from ...core.multi_stage_runner import run_multi
from ...core.heating import HeatingSchedule

from ..schemas.ca_result_schema import RxnCAResultDoc

import multiprocessing as mp
from typing import List, Tuple, Union, Dict
from dataclasses import dataclass

_heating = "_heating"
_rxns = "rxns"

def get_result(start):
    return run_multi(
        start,
        mp_globals[_rxns],
        mp_globals[_heating]
    )

@dataclass
class ParallelRxnMaker(Maker):

    name: str = "Rxn Sim. Parallel"

    @job(data="results")
    def make(self,
             chem_sys: Union[str, List[str]],
             heating_schedule: List[Tuple[int]],
             reactant_ratios: Dict[str, float],
             simulation_size: int = 21,
             simulation_dimension: int = 3):


        print("================= RETRIEVING AND SCORING REACTIONS =================")
        rxns = get_raw_rxns(chem_sys)
        heating_sched = HeatingSchedule(heating_schedule)
        library = heating_sched.calculate_rxns(rxns)
        print()
        print()
        print()

        print("================= SETTING UP SIMULATION =================")


        sim = setup(library.phases, reactant_ratios, simulation_size, dim=simulation_dimension)
        PROCESSES = mp.cpu_count()

        print(f'================= RUNNING SIMULATION w/ {PROCESSES} REALIZATIONS =================')


        global mp_globals

        mp_globals = {
            _heating: heating_sched,
            _rxns: library
        }


        with mp.get_context("fork").Pool(PROCESSES) as pool:
            results = pool.map(get_result, [sim for _ in range(PROCESSES)])

        jserializable = [res.as_dict() for res in results]
        
        return RxnCAResultDoc(
            chem_sys=chem_sys,
            heating_schedule=heating_schedule,
            reactant_amts=reactant_ratios,
            results=jserializable
        )

