from .reaction_result import ReactionResult
from .reaction_simulation import ReactionSimulation
from .reaction_controller import ReactionController
from .heating import HeatingSchedule
from ..reactions.reaction_library import ReactionLibrary

from pylattica.core import Runner, Simulation
from pylattica.core.constants import GENERAL, SITES

from .constants import TEMPERATURE

from typing import List, Dict


def run(simulation: Simulation, reaction_lib: ReactionLibrary, heating_schedule: HeatingSchedule, free_species=None, verbose=True, inertia=0):
    runner = Runner(is_async=True)
    print(f'Running simulation with inertia {inertia}')

    controller = ReactionController(
        simulation.structure,
        phases=reaction_lib.phases,
        free_species=free_species,
        inertia=inertia
    )
    
    results = []

    starting_state = simulation.state
    
    for step in heating_schedule.steps:
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
    prev_result.load_steps()
    new_starting_state = prev_result.last_step
    new_starting_state.batch_update(new_updates)
    return new_starting_state


def concatenate_results(results: List[ReactionResult]):
    starting_state = results[0].initial_state
    rxn_set = results[0].rxn_set

    new_result = ReactionResult(starting_state, rxn_set)

    for res in results:
        for d in res._diffs:
            new_result.add_step(d)
    
    return new_result