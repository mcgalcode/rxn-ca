from .enumerate_rxns_maker import EnumerateRxnsMaker
from ..schemas.enumerated_rxns_schema import EnumeratedRxnsModel
from jobflow import Flow, job

from ..utils.automaton_store import AutomatonStore
from rxn_network.reactions.reaction_set import ReactionSet

import typing

ENUMERATE_FLOW = "ENUMERATE_FLOW"

def enumerate_flow(chem_sys,
                   base_temp,
                   stability_cutoff=0.01,
                   open_element=None,
                   chempot=None,
                   formulas_to_include: typing.List = []
    ):
    enumerate_maker = EnumerateRxnsMaker()

    # If it doesn't, queue up a job to create one
    initial_job = enumerate_maker.make(
        chem_sys=chem_sys,
        temp=base_temp,
        stability_cutoff=stability_cutoff,
        open_el=open_element,
        chempot=chempot,
        formulas_to_include = formulas_to_include,
    )

    base_rxns: ReactionSet = initial_job.output.rxn_set

