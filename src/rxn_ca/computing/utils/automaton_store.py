from maggma.stores.mongolike import MongoStore
from rxn_network.reactions.reaction_set import ReactionSet
from jobflow.settings import JobflowSettings
from ..schemas.job_types import JobTypes
from .functions import format_chem_sys
from ...core.reaction_result import ReactionResult
from ...reactions.scored_reaction_set import ScoredReactionSet
from ...reactions import score_rxns


class AutomatonStore():

    def __init__(self, store: MongoStore = None):
        if store is None:
            settings = JobflowSettings()
            self.store = settings.JOB_STORE
            self.store.connect()
        else:
            self.store = store

    def list_available_sets(self, ):
        result = self.store.query({
            "output.job_type": JobTypes.ENUMERATE_RXNS.value
        }, {
            "output.chem_sys": True,
            "output.stability_cutoff": True,
            "output.temperature": True
        })
        systems = []

        for res in result:
            systems.append((
                res["output"]["chem_sys"],
                res["output"]["stability_cutoff"],
                res["output"]["temperature"]
            ))
        return systems
    
    def delete_rxn_sets(self, chem_sys):
        result = self.store.remove_docs({
            "output.job_type": JobTypes.ENUMERATE_RXNS.value,
            "output.chem_sys": chem_sys,
        })
        return result

    def get_raw_rxns(self, chem_sys, cutoff = None, temp = 300) -> ReactionSet:
        chem_sys = format_chem_sys(chem_sys)
        criteria = {
            "output.job_type": JobTypes.ENUMERATE_RXNS.value,
            "output.chem_sys": chem_sys,
            "output.temperature": temp,
        }
        if cutoff is not None:
            criteria["output.stability_cutoff"] = cutoff
        
        result = self.store.query_one(criteria)

        if result is not None:
            return ReactionSet.from_dict(result["output"]["rxn_set"])
        else:
            return None


    def get_automaton_result_by_task(self, task_id: str):
        result = self.store.query_one({
            "output.task_id": task_id,
            "output.job_type": JobTypes.RUN_RXN_AUTOMATON.value
        })
        if result is not None:
            return ReactionResult.from_dict(result["output"]["result"])
        else:
            return None

