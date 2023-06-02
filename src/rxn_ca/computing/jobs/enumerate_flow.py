from .enumerate_rxns_maker import EnumerateRxnsMaker
from .score_rxns_in_store import ScoreRxnsMaker

from jobflow import Flow

ENUMERATE_FLOW = "ENUMERATE_FLOW"

# Defunct
def enumerate_flow(chem_sys,
                    temp,
                    stability_cutoff=0.1,
                    open_element=None,
                    chempot=None,
    ):
    enumerate_maker = EnumerateRxnsMaker()

    enumerate_job = enumerate_maker.make(
        chem_sys=chem_sys,
        temp=temp,
        stability_cutoff=stability_cutoff,
        open_el=open_element,
        chempot=chempot
    )

    return Flow([enumerate_job], name = ENUMERATE_FLOW, output=enumerate_job.output)
