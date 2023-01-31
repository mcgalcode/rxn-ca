from maggma.stores.mongolike import MongoStore
from rxn_network.reactions.reaction_set import ReactionSet
from jobflow.settings import JobflowSettings
from ..schemas.job_types import JobTypes
from .functions import format_chem_sys
from ...core.reaction_result import ReactionResult
from ...reactions.scored_reaction_set import ScoredReactionSet

class AutomatonStore():

    def __init__(self, store = None):
        if store is None:
            settings = JobflowSettings()
            self.store = settings.JOB_STORE
            self.store.connect()
        else:
            self.store = store


    def get_rxn_enumeration_by_task_id(self, task_id: str):
        result = self.store.query_one({
                "output.task_id": task_id,
                "output.job_type": JobTypes.ENUMERATE_RXNS.value
            })
        if result is not None:
            return ReactionSet.from_dict(result["output"]["rxn_set"])
        else:
            return None

    def get_scored_rxns_by_task_id(self, task_id: str):
        result = self.store.query_one({
            "output.task_id": task_id,
            "output.job_type": JobTypes.SCORE_RXNS.value
        })
        if result is not None:
            return ScoredReactionSet.from_dict(result["output"]["scored_rxn_set"])
        else:
            return None

    def list_available_sets(self, ):
        result = self.store.query({
            "output.job_type": JobTypes.SCORE_RXNS.value
        }, {
            "output.chem_sys": True,
            "output.temperature": True
        })
        systems = []

        for res in result:
            systems.append((res["output"]["chem_sys"], res["output"]["temperature"]))
        return systems

    def get_scored_rxns(self, chem_sys, temperature):
        chem_sys = format_chem_sys(chem_sys)
        result = self.store.query_one({
            "output.job_type": JobTypes.SCORE_RXNS.value,
            "output.chem_sys": chem_sys,
            "output.temperature": temperature
        })
        if result is not None:
            return ScoredReactionSet.from_dict(result["output"]["scored_rxn_set"])
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

