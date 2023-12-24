

from ..core.reaction_controller import NB_HOOD_RADIUS
from ..phases.solid_phase_set import SolidPhaseSet
from ..core.reaction_setup import ReactionSetup

from typing import Dict

def setup_reaction(
        phases: SolidPhaseSet,
        size: int,
        phase_mole_amts: Dict = None,
        phase_mole_ratios: Dict = None,
        dim: int = 3,
        vol_scale = 1.0,
        vol_multiplier = 1.0,
    ):
    buffer_len = 1

    volume = size ** dim
    interaction_vol = 4 / 3 * NB_HOOD_RADIUS ** 3

    print(f"Simulation volume is {volume}")
    print(f"Interaction volume is {interaction_vol}")
    # num_desired_particles = volume /  (interaction_vol * particle_size)
    # num_allowed_particles = int((dim - 0.5) * volume / buffer_vol)

    num_sites = 500 # min(num_desired_particles, num_allowed_particles)
    print(f'Nucleating grains at {num_sites} sites')

    setup = ReactionSetup(phases, dim=dim)
    sim = setup.setup_growth(
        size,
        num_sites,
        phase_mol_amts=phase_mole_amts,
        phase_mol_ratios=phase_mole_ratios,
        volume_scale=vol_scale,
        volume_multiplier=vol_multiplier
    )
    
    return sim