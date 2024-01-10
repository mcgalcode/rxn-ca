from ..phases.solid_phase_set import SolidPhaseSet
from ..setup import ReactionPreparer

from typing import Dict

def setup_reaction(
        phases: SolidPhaseSet,
        precursor_mole_ratios: Dict,
        size: int = 15,
        vol_multiplier = 1.0,
    ):

    num_sites = 500 # min(num_desired_particles, num_allowed_particles)
    print(f'Nucleating grains at {num_sites} sites')

    preparer = ReactionPreparer(phases, dim=3)
    sim = preparer.prepare_reaction(
        size,
        num_sites,
        phase_mol_ratios=precursor_mole_ratios,
        volume_multiplier=vol_multiplier
    )
    
    return sim