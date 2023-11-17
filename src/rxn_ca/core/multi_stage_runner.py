from .reaction_result import ReactionResult
from .reaction_controller import ReactionController
from .heating import HeatingSchedule
from ..reactions.reaction_library import ReactionLibrary

from pylattica.core import Runner, Simulation
from pylattica.core.constants import GENERAL, SITES

from .constants import TEMPERATURE

from typing import List, Dict


def run_multi(simulation: Simulation,
              reaction_lib: ReactionLibrary,
              heating_schedule: HeatingSchedule,
              free_species=None,
              verbose=True,
              inertia=1,
              open_gas=None):
    runner = Runner(is_async=True)

    open_species = {}
    if open_gas is not None:
        open_species = {
            open_gas: 1.0
        }

    controller = ReactionController(
        simulation.structure,
        phases=reaction_lib.phases,
        free_species=free_species,
        inertia=inertia,
        open_species=open_species,
        heating_schedule=heating_schedule
    )
    
    results = []

    starting_state = simulation.state
    
    for step in heating_schedule.steps:
        print(f'Setting new temperature: {step.temp}')
        controller.set_rxn_set(reaction_lib.get_rxns_at_temp(step.temp))

        if len(results) > 0:
            new_updates = {
                GENERAL: {
                    TEMPERATURE: step.temp
                },
                SITES: {}
            }
            starting_state = get_new_starting_state(result, new_updates)

        result = runner.run(
            starting_state,
            controller,
            step.duration,
            structure=simulation.structure,
            verbose=verbose
        )

        results.append(result)


    result = concatenate_results(results)
    return result

def get_new_starting_state(prev_result: ReactionResult, new_updates: Dict):
    new_starting_state = prev_result.output
    new_starting_state.batch_update(new_updates)
    return new_starting_state


def concatenate_results(results: List[ReactionResult]):
    starting_state = results[0].initial_state

    new_result = ReactionResult(
        starting_state
    )

    for res in results:
        for d in res._diffs:
            new_result.add_step(d)
    
    return new_result