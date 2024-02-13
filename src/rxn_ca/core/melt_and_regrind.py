from ..phases import SolidPhaseSet
from ..utilities.setup_reaction import setup_reaction
from ..analysis import ReactionStepAnalyzer
from .constants import TEMPERATURE, MELTED_AMTS, VOL_MULTIPLIER, GASES_CONSUMED, GASES_EVOLVED

from pylattica.core import SimulationState
from pylattica.core.constants import GENERAL, SITES
from typing import Dict
from tabulate import tabulate
import numpy as np

REGRIND_CUTOFF = 0.15

def _print_dict(d: Dict, key_header: str, val_header: str):
    table = []
    for k, v in d.items():
        table.append([k, v])
    print(tabulate(table, headers=[key_header, val_header]))

def calculate_melted_fraction(step: SimulationState, phases: SolidPhaseSet, temp: int):
    analyzer = ReactionStepAnalyzer(phases)

    unmelted_vol_fractions = analyzer.get_all_volume_fractions(step, include_melted=False)

    total_melted_vol_frac = 0
    for phase, vol_frac in unmelted_vol_fractions.items():
        if phases.is_melted(phase, temp):
            total_melted_vol_frac += vol_frac

    return total_melted_vol_frac

def calculate_solid_ratio(step: SimulationState, phases: SolidPhaseSet, temp: int):
    analyzer = ReactionStepAnalyzer(phases)

    ideal_grid_vol = analyzer.get_ideal_step_volume(step)
    solid_vol = analyzer.get_total_solid_volume(step, temp)

    return solid_vol / ideal_grid_vol

def separate_solid_and_melt(step: SimulationState, phases: SolidPhaseSet, temp: int):
    analyzer = ReactionStepAnalyzer(phases)

    print(f"Total volume of melted material greater than {100 * REGRIND_CUTOFF}%, calculating new state...")
    
    melted_vols = analyzer.get_absolute_melted_volumes(step, temp)
    print("\nMelted vols from prev state: ")
    _print_dict(melted_vols, "Phase", "Volume")

    solid_vols = analyzer.get_absolute_solid_volumes(step, temp)
    print("\nSolid vols from prev state: ")
    _print_dict(solid_vols, "Phase", "Volume")
    print("Total: ", analyzer.get_total_solid_volume(step, temp))

    solid_moles = phases.vol_amts_to_moles(solid_vols)

    print(f"Previous real volume: {analyzer.get_total_volume(step)}")
    print(f"Previous ideal volume: {analyzer.get_ideal_step_volume(step)}")
    sim_size = analyzer.get_simulation_size(step)
    print(f'Deduced sim size to be {sim_size}')

    print("\nConstructing new state with molar solid composition: ", solid_moles)
    prev_vol_multiplier = step.get_general_state().get(VOL_MULTIPLIER, 1)
    solid_ratio = calculate_solid_ratio(step, phases, temp)
    new_vol_multiplier = prev_vol_multiplier * solid_ratio
    print(f"Previous volume multiplier: {prev_vol_multiplier}, current solid / total ratio: {solid_ratio}")
    print("Calculated new volume multipler: ", new_vol_multiplier)

    new_solid_state = setup_reaction(
        phases,
        precursor_mole_ratios=solid_moles,
        vol_multiplier=new_vol_multiplier,
        size=sim_size
    ).state


    new_updates = {
        GENERAL: {
            TEMPERATURE: temp,
            MELTED_AMTS: melted_vols,
            VOL_MULTIPLIER: new_vol_multiplier
        },
        SITES: {}
    }

    new_solid_state.batch_update(new_updates)

    return new_solid_state    

def melt_and_regrind(step: SimulationState, phases: SolidPhaseSet, temp: int):
    total_melted_vol_frac = calculate_melted_fraction(step, phases, temp)
    should_recalc_melted = total_melted_vol_frac > REGRIND_CUTOFF
    analyzer = ReactionStepAnalyzer(phases)

    if should_recalc_melted:
        separated = separate_solid_and_melt(step, phases, temp)

        old_vols = analyzer.get_all_absolute_phase_volumes(step)
        new_vols = analyzer.get_all_absolute_phase_volumes(separated)
        for phase, old_vol in old_vols.items():
            new_vol = new_vols.get(phase)
            if not np.isclose(old_vol, new_vol, rtol=0.011):
                print(old_vols)
                print(new_vols)
                raise RuntimeError(f"After melting, volume was not conserved: {phase} {old_vol} -> {new_vol}")
        
        return separated
    else:
        new_solid_state = step.copy()

        new_updates = {
            GENERAL: {
                TEMPERATURE: temp,
            },
            SITES: {}
        }

        new_solid_state.batch_update(new_updates)
        return new_solid_state
