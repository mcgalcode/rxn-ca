from ..core.reaction_result import ReactionResult
from ..core.reaction_controller import ReactionController
from ..core.reaction_calculator import ReactionCalculator
from ..core.heating import HeatingSchedule
from ..core.constants import GASES_EVOLVED, GASES_CONSUMED, MELTED_AMTS
from ..reactions.reaction_library import ReactionLibrary
from ..core.melt_and_regrind import melt_and_regrind

from pylattica.core import AsynchronousRunner, Simulation, BasicController

from typing import List, Callable
import numpy as np

class HeatingScheduleRunner():

    def __init__(self, middlewares: List[Callable] = []) -> None:
        self._middlewares = middlewares
        
    def run_multi(self,
                simulation: Simulation,
                reaction_lib: ReactionLibrary,
                heating_schedule: HeatingSchedule,
                controller: BasicController,
                verbose=True):
        runner = AsynchronousRunner()       
        results = []

        starting_state = simulation.state

        # One "step" is visiting every site once
        step_size = len(simulation.structure.site_ids)
        total_steps = len(heating_schedule)

        prev_temp = None

        assert GASES_EVOLVED in starting_state.get_general_state()
        assert MELTED_AMTS in starting_state.get_general_state()
        assert GASES_CONSUMED in starting_state.get_general_state()

        for step_no, step in enumerate(heating_schedule.steps):
            print(f'Running step {step_no + 1} of {total_steps}.')
            if step.temp != prev_temp:
                print(f'Setting new temperature: {step.temp}')
            
            prev_temp = step.temp
            controller.set_temperature(step.temp)
            controller.set_rxn_set(reaction_lib.get_rxns_at_temp(step.temp))

            num_simulation_steps = int(step_size * step.duration)

            if len(results) > 0:
                starting_state = results[-1].output

                for middleware in self._middlewares:
                    starting_state = middleware(starting_state, reaction_lib.phases, step.temp)

            result = runner.run(
                starting_state,
                controller,
                num_simulation_steps,
                verbose=verbose
            )

            results.append(result)


        result = concatenate_results(results)
        return result
    
class MeltAndRegrindMultiRunner(HeatingScheduleRunner):

    def __init__(self) -> None:
        super().__init__([melt_and_regrind])

def concatenate_results(results: List[ReactionResult]):
    starting_state = results[0].initial_state

    new_result = ReactionResult(
        starting_state
    )

    for idx, res in enumerate(results):
        if idx > 0:
            new_result.add_step(res.first_step._state)
        for d in res._diffs:
            new_result.add_step(d)
    
    return new_result