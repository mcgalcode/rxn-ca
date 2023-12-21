from .reaction_result import ReactionResult
from .reaction_controller import ReactionController
from ..phases import SolidPhaseSet
from .heating import HeatingSchedule
from ..reactions.reaction_library import ReactionLibrary

from .. import setup_reaction

from ..analysis import ReactionStepAnalyzer

from pylattica.core import AsynchronousRunner, Simulation, SimulationState
from pylattica.core.constants import GENERAL, SITES

from .constants import TEMPERATURE, MELTED_AMTS, VOL_MULTIPLIER

from typing import List


def run_multi(simulation: Simulation,
              reaction_lib: ReactionLibrary,
              heating_schedule: HeatingSchedule,
              free_species=None,
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
        free_species=free_species,
        inertia=inertia,
        open_species=open_species,
    )
    
    results = []

    starting_state = simulation.state
    
    for step in heating_schedule.steps:
        print(f'Setting new temperature: {step.temp}')
        controller.set_rxn_set(reaction_lib.get_rxns_at_temp(step.temp))

        if len(results) > 0:
            starting_state = melt_and_regrind(
                result.output,
                reaction_lib.phases,
                step.temp,
                simulation
            )

        result = runner.run(
            starting_state,
            controller,
            step.duration,
            verbose=verbose
        )

        results.append(result)


    result = concatenate_results(results)
    return result

def melt_and_regrind(previous_final_step: SimulationState, phases: SolidPhaseSet, new_temp: int, sim: Simulation):
    melted_phases = [p for p in phases.phases if phases.get_melting_point(p) is not None and phases.get_melting_point(p) < new_temp]

    analyzer = ReactionStepAnalyzer(phases)

    solid_vol_breakdown = analyzer.phase_volume_fractions(previous_final_step, include_melted=False)

    total_melted_vol_frac = 0
    for phase, vol_frac in solid_vol_breakdown.items():
        if phase in melted_phases:
            total_melted_vol_frac += vol_frac

    should_recalc_melted = total_melted_vol_frac > 0.25

    if should_recalc_melted:
        volume_breakdown = analyzer.phase_volumes(previous_final_step)
        
        melted_vols = { p: amt for p, amt in volume_breakdown.items() if p in melted_phases }
        print("Melted vols from prev state: ", melted_vols)

        solid_vols = { p: amt for p, amt in volume_breakdown.items() if p not in melted_phases }
        print("Solid vols from prev state: ", solid_vols)

        total_vol = analyzer.total_volume(previous_final_step, include_melted=False)

        total_solid_vol = sum(solid_vols.values())

        new_vol_scale = total_solid_vol / total_vol
            
        # Convert current volume composition to moles
        current_total_molar_comp = analyzer.molar_breakdown(previous_final_step, include_melted=False)
        current_total_molar_comp = { p: amt for p, amt in current_total_molar_comp.items() if p not in melted_phases }

        print("Building new sim w", current_total_molar_comp)
        num_sites = len(sim.structure.site_ids)
        sim_size = int(round(num_sites ** (1 / 3)))
        prev_vol_multiplier = previous_final_step.get_general_state().get(VOL_MULTIPLIER, 1)
        new_vol_multiplier = prev_vol_multiplier * new_vol_scale
        print("Using volume multipler: ", new_vol_multiplier)

        new_solid_state = setup_reaction(
            phases,
            sim_size,
            phase_mole_amts=current_total_molar_comp,
            particle_size = 0.5,
            vol_multiplier=new_vol_multiplier
        ).state


        new_updates = {
            GENERAL: {
                TEMPERATURE: new_temp,
                MELTED_AMTS: melted_vols,
                VOL_MULTIPLIER: new_vol_multiplier
            },
            SITES: {}
        }

        new_solid_state.batch_update(new_updates)

        return new_solid_state
    else:
        # current_total_molar_comp = analyzer.molar_breakdown(previous_final_step, include_melted=False)

        # num_sites = len(sim.structure.site_ids)
        # sim_size = int(round(num_sites ** (1 / 3)))

        # new_solid_state = setup_reaction(
        #     phases,
        #     sim_size,
        #     phase_mole_amts=current_total_molar_comp,
        #     particle_size = 0.5,
        #     vol_scale = 1.0,
        # ).state

        new_solid_state = previous_final_step.copy()

        new_updates = {
            GENERAL: {
                TEMPERATURE: new_temp,
            },
            SITES: {}
        }

        new_solid_state.batch_update(new_updates)
        return new_solid_state

def get_new_temp_state(previous_final_step: SimulationState, phases: SolidPhaseSet, new_temp: int, sim: Simulation):

    new_state = previous_final_step.copy()

    new_updates = {
        GENERAL: {
            TEMPERATURE: new_temp,
        },
        SITES: {}
    }

    new_state.batch_update(new_updates)
    return new_state


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