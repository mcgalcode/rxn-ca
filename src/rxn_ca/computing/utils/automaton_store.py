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

    def get_rxn_enumeration_by_task_id(self, task_id: str):
        result = self.store.query_one({
                "output.task_id": task_id,
                "output.job_type": JobTypes.ENUMERATE_RXNS.value
            })
        if result is not None:
            return ReactionSet.from_dict(result["output"]["rxn_set"])
        else:
            return None

    def list_available_sets(self, ):
        result = self.store.query({
            "output.job_type": JobTypes.ENUMERATE_RXNS.value
        }, {
            "output.chem_sys": True,
        })
        systems = []

        for res in result:
            systems.append((res["output"]["chem_sys"]))
        return systems
    
    def delete_rxn_sets(self, chem_sys):
        result = self.store.remove_docs({
            "output.job_type": JobTypes.ENUMERATE_RXNS.value,
            "output.chem_sys": chem_sys,
        })
        return result

    def get_raw_rxns(self, chem_sys, cutoff = None) -> ReactionSet:
        chem_sys = format_chem_sys(chem_sys)
        criteria = {
            "output.job_type": JobTypes.ENUMERATE_RXNS.value,
            "output.chem_sys": chem_sys,
        }
        if cutoff is not None:
            criteria["output.stability_cutoff"] = cutoff

        result = self.store.query_one(criteria)
        # print(result["output"]["rxn_set"])

        if result is not None:
            return ReactionSet.from_dict(result["output"]["rxn_set"])
        else:
            return None

    def get_scored_rxns(self, chem_sys, temperature, scorer, **kwargs):
        raw: ReactionSet = self.get_raw_rxns(chem_sys, **kwargs)
        if raw is None:
            raise RuntimeError(f'No reactions found for chem sys {chem_sys}')
        
        print(f'Calculating reactions at new temperature {temperature}')
        raw.set_new_temperature(temperature)
        rxns, phases = score_rxns(raw, scorer)
        return ScoredReactionSet(rxns, phases)


    def get_automaton_result_by_task(self, task_id: str):
        result = self.store.query_one({
            "output.task_id": task_id,
            "output.job_type": JobTypes.RUN_RXN_AUTOMATON.value
        })
        if result is not None:
            return ReactionResult.from_dict(result["output"]["result"])
        else:
            return None

