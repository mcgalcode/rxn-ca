
from rxn_ca.reactions.scored_reaction import ScoredReaction
from rxn_ca.reactions.scored_reaction_set import ScoredReactionSet
from .core.reaction_controller import ReactionController, NB_HOOD_RADIUS
from .phases.solid_phase_set import SolidPhaseSet
from .core.reaction_setup import ReactionSetup
from .core.heating import HeatingSchedule
from .reactions.scorers import BasicScore, TammanHuttigScoreErf, score_rxns
from .reactions.reaction_library import ReactionLibrary

from typing import Dict

from .computing import AutomatonStore, enumerate_flow

from jobflow.managers.local import run_locally
from pylattica.core import AsynchronousRunner, Simulation

from typing import List
from rxn_network.reactions.reaction_set import ReactionSet
from tqdm import tqdm

def setup_reaction(
        phases: SolidPhaseSet,
        size: int,
        phase_mole_amts: Dict = None,
        phase_mole_ratios: Dict = None,
        dim: int = 3,
        particle_size = 1,
        vol_scale = 1.0,
        vol_multiplier = 1.0,
    ):
    buffer_len = 1

    volume = size ** dim
    buffer_vol = buffer_len ** dim
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

def run(simulation: Simulation,
        num_steps: int,
        free_species=None,
        verbose=True,
        open_gas=None,
        inertia=0):
    runner = AsynchronousRunner(is_async=True)
    print(f'Running simulation with inertia {inertia}')

    open_species = None
    if open_gas is not None:
        open_species = {
            open_gas: 1.0
        }

    controller = ReactionController(
        simulation.structure,
        simulation.reactions,
        free_species=free_species,
        inertia=inertia,
        open_species=open_species
    )

    result = runner.run(
        simulation.state,
        controller,
        num_steps,
        structure=simulation.structure,
        verbose=verbose
    )
    return result

def get_raw_rxns(chemsys: str):
    store = AutomatonStore()
    rxn_set = store.get_raw_rxns(chemsys)
    return ReactionLibrary(rxn_set)
