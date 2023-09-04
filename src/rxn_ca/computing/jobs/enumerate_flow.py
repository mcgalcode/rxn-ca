from .enumerate_rxns_maker import EnumerateRxnsMaker
from ..schemas.enumerated_rxns_schema import EnumeratedRxnsModel
from jobflow import Flow, job

from ..utils.automaton_store import AutomatonStore
from rxn_network.reactions.reaction_set import ReactionSet

import typing

ENUMERATE_FLOW = "ENUMERATE_FLOW"

def enumerate_flow(chem_sys,
                   base_temp,
                   stability_cutoff=0.1,
                   open_element=None,
                   chempot=None,
                   other_temps = [],
                   formulas_to_include: typing.List = []
    ):
    enumerate_maker = EnumerateRxnsMaker()

    store = AutomatonStore()
    
    # Retrieve a reaction set if one already exists
    base_rxns = store.get_raw_rxns(
        chem_sys=chem_sys,
        cutoff=stability_cutoff,
        temp=base_temp
    )

    jobs = []
    # If it doesn't, queue up a job to create one
    if base_rxns is None:
        initial_job = enumerate_maker.make(
            chem_sys=chem_sys,
            temp=base_temp,
            stability_cutoff=stability_cutoff,
            open_el=open_element,
            chempot=chempot,
            formulas_to_include = formulas_to_include,
        )

        jobs.append(initial_job)
        base_rxns: ReactionSet = initial_job.output.rxn_set

    for t in other_temps:
        # Check if reactions exist for this temperature
        existing = store.get_raw_rxns(
            chem_sys,
            cutoff=stability_cutoff,
            temp = t
        )
        # If they don't exist, queue up a job using the reactiosn
        # we found at the base temp.
        if existing is None:
            new_rxns_job = get_rxns_new_temp(base_rxns,
                                             chem_sys,
                                             t,
                                             stability_cutoff,
                                             open_element,
                                             chempot)
            jobs.append(new_rxns_job)


    if len(jobs) == 0:
        print(f"No work needed for {chem_sys}")
        # return 
    else:
        return Flow(jobs, name = ENUMERATE_FLOW, output=jobs[0].output)

@job
def get_rxns_new_temp(base_rxns,
                      chem_sys,
                      new_temp,
                      stability_cutoff,
                      open_el,
                      chem_pot):
    print(f"Enumerating {chem_sys} at {new_temp} w cutoff {stability_cutoff} using base reactions")    
    new_rxns = base_rxns.set_new_temperature(new_temp)
    return EnumeratedRxnsModel.from_obj(
        new_rxns,
        chem_sys,
        stability_cutoff=stability_cutoff,
        open_el = open_el,
        chem_pot = chem_pot,
        temperature=new_temp
    )

