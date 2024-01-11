from ..core.reaction_result import ReactionResult
from ..core.reaction_controller import ReactionController
from ..core.heating import HeatingSchedule
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
    
    for step in heating_schedule.steps:
        print(f'Setting new temperature: {step.temp}')
        controller.set_rxn_set(reaction_lib.get_rxns_at_temp(step.temp))

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
                if not np.isclose(old_vol, new_vol, rtol=0.01):
                    print(old_vols)
                    print(new_vols)
                    raise RuntimeError(f"After melting, volume was not conserved: {phase} {old_vol} -> {new_vol}")
            

        result = runner.run(
            starting_state,
            controller,
            step.duration,
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

    for res in results:
        new_result.add_step(res.first_step._state)
        for d in res._diffs:
            new_result.add_step(d)
    
    return new_result