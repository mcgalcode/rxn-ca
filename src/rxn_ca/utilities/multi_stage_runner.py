from ..core.reaction_result import ReactionResult
from ..core.reaction_controller import ReactionController
from ..core.heating import HeatingSchedule
from ..core.constants import GASES_EVOLVED, GASES_CONSUMED, MELTED_AMTS
from ..reactions.reaction_library import ReactionLibrary
from ..core.melt_and_regrind import melt_and_regrind
from ..analysis import ReactionStepAnalyzer

from pylattica.core import AsynchronousRunner, Simulation

from typing import List
import numpy as np

def run_multi(simulation: Simulation,
              reaction_lib: ReactionLibrary,
              heating_schedule: HeatingSchedule,
              verbose=True,
              inertia=1,
              open_gas=None):
    runner = AsynchronousRunner()

    open_species = {}
    if open_gas is not None:
        open_species = {
            open_gas: 1.0
        }

    controller = ReactionController(
        simulation.structure,
        inertia=inertia,
        open_species=open_species,
    )
    
    results = []

    starting_state = simulation.state
    analyzer = ReactionStepAnalyzer(reaction_lib.phases)

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
        controller.set_rxn_set(reaction_lib.get_rxns_at_temp(step.temp))

        num_simulation_steps = step_size * step.duration

        if len(results) > 0:
            starting_state = melt_and_regrind(
                result.output,
                reaction_lib.phases,
                step.temp,
            )

            old_vols = analyzer.get_all_absolute_phase_volumes(result.output)
            new_vols = analyzer.get_all_absolute_phase_volumes(starting_state)
            for phase, old_vol in old_vols.items():
                new_vol = new_vols.get(phase)
                if not np.isclose(old_vol, new_vol, rtol=0.011):
                    print(old_vols)
                    print(new_vols)
                    raise RuntimeError(f"After melting, volume was not conserved: {phase} {old_vol} -> {new_vol}")
            

        result = runner.run(
            starting_state,
            controller,
            num_simulation_steps,
            verbose=verbose
        )

        results.append(result)


    result = concatenate_results(results)
    return result

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