from ..core.reaction_result import ReactionResult
from ..core.reaction_controller import ReactionController
from ..core.reaction_calculator import ReactionCalculator
from ..core.heating import HeatingSchedule, RegrindStep, HeatingStep
from ..core.constants import GASES_EVOLVED, GASES_CONSUMED, MELTED_AMTS
from ..reactions.reaction_library import ReactionLibrary
from ..core.melt_and_regrind import melt_and_regrind
from ..analysis.reaction_step_analyzer import ReactionStepAnalyzer
from .setup_reaction import setup_noise_reaction

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
        sim_size = int(step_size ** (1 / 3))
        total_steps = len(heating_schedule)

        prev_temp = None

        assert GASES_EVOLVED in starting_state.get_general_state()
        assert MELTED_AMTS in starting_state.get_general_state()
        assert GASES_CONSUMED in starting_state.get_general_state()

        reground_state = None

        for step_no, step in enumerate(heating_schedule.steps):
            if isinstance(step, HeatingStep):
                print(f'Running step {step_no + 1} of {total_steps}.')
                if step.temperature != prev_temp:
                    print(f'Setting new temperature: {step.temperature}')
                
                prev_temp = step.temperature
                controller.set_temperature(step.temperature)
                controller.set_rxn_set(reaction_lib.get_rxns_at_temp(step.temperature))

                num_simulation_steps = int(step_size * step.duration)

                if reground_state is not None:
                    starting_state = reground_state
                    reground_state = None
                if len(results) > 0:
                    starting_state = results[-1].output

                    for middleware in self._middlewares:
                        starting_state = middleware(starting_state, reaction_lib.phases, step.temperature)

                result = runner.run(
                    starting_state,
                    controller,
                    num_simulation_steps,
                    verbose=verbose
                )

                results.append(result)
            elif isinstance(step, RegrindStep):
                analyzer = ReactionStepAnalyzer(reaction_lib.phases)
                last_step = results[-1].output
                amts = analyzer.get_all_mole_fractions(last_step)
                new_amts = { p: amt for p, amt in amts.items() if amt > 0.01}

                reground_state = setup_noise_reaction(
                    reaction_lib.phases,
                    precursor_mole_ratios = new_amts,
                    size = sim_size,
                )

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