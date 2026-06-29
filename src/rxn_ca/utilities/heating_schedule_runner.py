from ..core.reaction_result import ReactionResult
from ..core.reaction_controller import ReactionController
from ..core.reaction_calculator import ReactionCalculator
from ..core.heating import HeatingSchedule, RegrindStep, HeatingStep
from ..core.constants import GASES_EVOLVED, GASES_CONSUMED, MELTED_AMTS, TEMPERATURE
from ..reactions.reaction_library import ReactionLibrary
from ..core.melt_and_regrind import melt_and_regrind
from ..analysis.reaction_step_analyzer import ReactionStepAnalyzer
from ..analysis.phase_volume_observer import PhaseVolumeObserver
from .setup_reaction import setup_noise_reaction

from pylattica.core import AsynchronousRunner, ObservedResult, Simulation, BasicController

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
                verbose=True,
                compress_freq: int = None,
                num_frames: int = None,
                analysis_only: bool = False):
        runner = AsynchronousRunner()
        results: List[ReactionResult] = []

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

                print("Setting temperature state")
                starting_state.set_general_state({TEMPERATURE: step.temperature })

                if analysis_only:
                    # Record only the reduced per-phase-volume view (no full
                    # states). compress_freq / num_frames set the observation
                    # cadence here instead of a frame-retention cadence.
                    result = runner.run(
                        starting_state,
                        controller,
                        num_simulation_steps,
                        verbose=verbose,
                        retain_history=False,
                        observers=[PhaseVolumeObserver()],
                        observe_freq=compress_freq,
                        observe_num_frames=num_frames,
                    )
                else:
                    result = runner.run(
                        starting_state,
                        controller,
                        num_simulation_steps,
                        verbose=verbose,
                        compress_freq=compress_freq,
                        num_frames=num_frames,
                    )

                results.append(result)
            elif isinstance(step, RegrindStep):
                analyzer = ReactionStepAnalyzer(reaction_lib.phases)
                analyzer.set_step_group(results[-1].output)
                amts = analyzer.get_all_mole_fractions()
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
    first_result = results[0]

    # Analysis-only results retain no frames/diffs -- only the merged per-frame
    # observations and the final live state need to be carried forward.
    if not first_result.retain_history:
        return _concatenate_analysis_only(results)

    starting_state = first_result.initial_state

    # Preserve the compression strategy from the first result. Compression is
    # now configured on the result rather than passed to the constructor.
    new_result = ReactionResult(starting_state)
    if first_result.live_compress:
        new_result.configure_compression(compress_freq=first_result.compress_freq)

    # Track cumulative step offset for frame indices
    step_offset = 0

    for idx, res in enumerate(results):
        if idx > 0:
            new_result.add_step(res.first_step._state)
            step_offset += 1

        # If live_compress mode, copy frames with adjusted indices
        if res.live_compress and res._frames:
            for frame_step, frame_state in res._frames.items():
                if frame_step == 0 and idx > 0:
                    # Skip frame 0 for subsequent results (already have that state)
                    continue
                new_result._frames[step_offset + frame_step] = frame_state
            # Update offset by the max frame step
            step_offset += max(res._frames.keys()) if res._frames else 0
            new_result._total_steps = step_offset
        else:
            # Normal mode: copy diffs
            for d in res._diffs:
                new_result.add_step(d)

    return new_result


def _concatenate_analysis_only(results: List[ReactionResult]) -> ReactionResult:
    """Concatenate analysis-only (retain_history=False) heating-step results.

    Each input carries a single ObservedResult of per-phase volumes. The
    observations are merged into one ObservedResult with step numbers offset by
    the cumulative step count, mirroring how concatenate_results stitches frames
    together. The boundary observation (step 0) of each subsequent segment is
    dropped because it duplicates the previous segment's final state.
    """
    new_result = ReactionResult(results[0].initial_state, retain_history=False)

    merged = ObservedResult(
        observer=None, compress_freq=results[0].observed_results[0].compress_freq
    )

    step_offset = 0
    for idx, res in enumerate(results):
        obs = res.observed_results[0]
        steps = obs.observation_steps
        for step in steps:
            if step == 0 and idx > 0:
                continue
            merged._observations[step_offset + step] = obs.observation_at(step)
        step_offset += max(steps) if steps else 0

    merged._total_steps = step_offset

    new_result._observed_results = [merged]
    new_result._live_state = results[-1].output.copy()
    new_result._total_steps = step_offset
    return new_result