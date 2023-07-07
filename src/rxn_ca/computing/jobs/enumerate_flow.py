from .enumerate_rxns_maker import EnumerateRxnsMaker
from jobflow import Flow

import typing

ENUMERATE_FLOW = "ENUMERATE_FLOW"

def enumerate_flow(chem_sys,
                   temp,
                   stability_cutoff=0.1,
                   open_element=None,
                   chempot=None,
                   formulas_to_include: typing.List = []
    ):
    enumerate_maker = EnumerateRxnsMaker()

    enumerate_job = enumerate_maker.make(
        chem_sys=chem_sys,
        temp=temp,
        stability_cutoff=stability_cutoff,
        open_el=open_element,
        chempot=chempot,
        formulas_to_include = formulas_to_include,
    )

    return Flow([enumerate_job], name = ENUMERATE_FLOW, output=enumerate_job.output)
