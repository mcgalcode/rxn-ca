from .enumerate_rxns_maker import EnumerateRxnsMaker
from .score_rxns_in_store import ScoreRxnsMaker

from jobflow import Flow

ENUMERATE_AND_SCORE_FLOW = "ENUMERATE_AND_SCORE_FLOW"

def enumerate_and_score_flow(chem_sys,
                             temp,
                             stability_cutoff=0.1,
                             open_element=None,
                             chempot=None,
    ):
    enumerate_maker = EnumerateRxnsMaker()
    score_maker = ScoreRxnsMaker()
    enumerate_job = enumerate_maker.make(
        chem_sys=chem_sys,
        temp=temp,
        stability_cutoff=stability_cutoff,
        open_el=open_element,
        chempot=chempot
    )
    score_job = score_maker.make(
        chem_sys=chem_sys,
        temp=temp,
        rxns=enumerate_job.output,
    )
    return Flow([enumerate_job, score_job], name = ENUMERATE_AND_SCORE_FLOW, output=score_job.output)
