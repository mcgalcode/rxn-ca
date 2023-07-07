from typing import Optional
from jobflow.core.maker import Maker
from jobflow.core.job import job
from jobflow import Response, Flow
from pylattica.core import Runner

from dataclasses import dataclass

from ..schemas.ca_result_schema import RxnCAResultModel
from ..schemas.scored_rxns_schema import ScoredRxnsModel
from ..utils.automaton_store import AutomatonStore
from ...core.reaction_result import ReactionResult

from ...core import SolidPhaseSet, ReactionSetup, ReactionController
from ...reactions import ScoredReactionSet

from uuid import uuid4


@dataclass
class RunRxnAutomatonMaker(Maker):
    name: str = "Run Rxn Automaton"

    @job
    def make(self,
            chem_sys: str,
            temp: int,
            setup_style: str,
            setup_args: dict,
            num_steps: int,
            scored_rxns: ScoredRxnsModel = None,
            scored_rxns_task_id: str = None,
            db_connection_params: dict = {},
            dimensionality: int = 2,
            inertia: float = 1,
            open_species: dict = {},
            free_species: list = [],
            parallel: bool = True,
            db_file: Optional[str] = None,
        ):

        scored_rxn_set = None

        if scored_rxns is not None:
            scored_rxn_set: ScoredReactionSet = ScoredReactionSet.from_dict(scored_rxns.scored_rxn_set)

        if scored_rxn_set is None:

            if db_file is not None:
                store = AutomatonStore.from_db_file(db_file)
            else:
                store = AutomatonStore(**db_connection_params)

            store.connect()

            if scored_rxns_task_id is not None:
                scored_rxn_set = store.get_scored_rxns_by_task_id(scored_rxns_task_id)

            if scored_rxn_set is None:
                scored_rxn_set = store.get_scored_rxns(chem_sys, temperature=temp)

        assert scored_rxn_set is not None, "ScoredRxnSet not found!"
        phase_map: SolidPhaseSet = SolidPhaseSet(scored_rxn_set.phases)

        setup = ReactionSetup(phase_map, volumes = scored_rxn_set.volumes)


        if setup_style == "random_mixture":
            initial_state = setup.setup_random_mixture(**setup_args)
        elif setup_style == "interface":
            initial_state = setup.setup_interface(**setup_args)

        nb = ReactionController.get_neighborhood_from_step(initial_state)
        controller = ReactionController(
            phase_map=phase_map,
            neighborhood=nb,
            scored_rxns=scored_rxn_set,
            inertia=inertia,
            open_species = open_species,
            free_species = free_species,
            temperature = temp
        )

        runner = Runner(parallel=parallel)
        result: ReactionResult = runner.run(initial_state, controller, num_steps)

        return RxnCAResultModel(
            task_id=str(uuid4()),
            chem_sys=chem_sys,
            temperature=temp,
            result=result.as_dict()
        )