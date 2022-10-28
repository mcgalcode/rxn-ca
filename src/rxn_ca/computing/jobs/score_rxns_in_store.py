from re import M
from jobflow.core.maker import Maker
from jobflow.core.job import job

from dataclasses import dataclass

from ..schemas.enumerated_rxns_schema import EnumeratedRxnsModel
from ..schemas.scored_rxns_schema import ScoredRxnsModel
from ..schemas.job_types import JobTypes

from ..utils.automaton_store import AutomatonStore

from ...reactions import ScoredReactionSet, ArrheniusScore, score_rxns

from rxn_network.reactions.reaction_set import ReactionSet


@dataclass
class ScoreRxnsMaker(Maker):

    name: str = "Score Rxns In Store"

    @job
    def make(self,
             chem_sys: str,
             temp: int,
             rxns: EnumeratedRxnsModel = None,
        ):
        rxn_set = ReactionSet.from_dict(rxns.rxn_set)
        scorer = ArrheniusScore(temp)
        scored_rxns = score_rxns(rxn_set, scorer)
        scored_rxn_set = ScoredReactionSet(scored_rxns)
        result_model = ScoredRxnsModel.from_obj(
            scored_rxn_set,
            chem_sys = chem_sys,
            temp=temp
        )

        return result_model

