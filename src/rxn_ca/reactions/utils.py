from mp_api.client import MPRester
import math

def get_phase_vols(phases):

    volumes = {}

    with MPRester() as mpr:
        res = mpr.summary.search(
            formula = phases,
            fields=["structure", "composition", "task_id", "energy_above_hull"]
        )
        min_energies = {}
        for item in res:
            struct = item.structure
            comp = item.composition
            red_form = comp.reduced_formula
            
            if item.energy_above_hull < min_energies.get(red_form, math.inf):
                volumes[comp.reduced_formula] = struct.volume / comp.get_reduced_composition_and_factor()[1]
                min_energies[red_form] = item.energy_above_hull

    return volumes