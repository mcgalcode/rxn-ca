
from .core.reaction_controller import ReactionController, NB_HOOD_RADIUS
from .core.reaction_result import ReactionResult
from .core.solid_phase_set import SolidPhaseSet
from .core.reaction_setup import ReactionSetup
from .core.reaction_simulation import ReactionSimulation
from .core.heating import HeatingSchedule

from .reactions.scored_reaction import ScoredReaction
from .reactions.scored_reaction_set import ScoredReactionSet
from .reactions.reaction_library import ReactionLibrary

from typing import Dict

from .computing import AutomatonStore, enumerate_flow

from jobflow.managers.local import run_locally
from pylattica.core import Runner, Simulation

import typing

def setup(phases: SolidPhaseSet, phase_mixture: Dict, size: int, dim: int = 2, particle_size = 1):
    phase_names = []
    phase_amts = []
    buffer_len = 5

    volume = size ** dim
    buffer_vol = buffer_len ** dim
    interaction_vol = 4 / 3 * NB_HOOD_RADIUS ** 3

    print(f"Simulation volume is {volume}")
    print(f"Interaction volume is {interaction_vol}")
    num_desired_particles = volume /  (interaction_vol * particle_size)
    num_allowed_particles = int((dim - 0.5) * volume / buffer_vol)

    num_sites = min(num_desired_particles, num_allowed_particles)
    print(f'Nucleating grains at {num_sites} sites')
    for name, amt in phase_mixture.items():
        phase_names.append(name)
        phase_amts.append(amt)

    setup = ReactionSetup(phases, dim=dim)
    sim = setup.setup_growth(size, num_sites, phase_names, phase_mol_ratios=phase_amts)
    return sim

def run(simulation: Simulation, num_steps: int, free_species=None, verbose=True, inertia=0):
    runner = Runner(is_async=True)
    print(f'Running simulation with inertia {inertia}')
    controller = ReactionController(
        simulation.structure,
        simulation.reactions,
        free_species=free_species,
        inertia=inertia
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

def enumerate_rxns(chem_sys,
                   temp = 300,
                   stability_cutoff=0.1,
                   open_element=None,
                   chempot=None,
                   formulas_to_include: typing.List = []
    ):
    store = AutomatonStore()
    flow = enumerate_flow(
        chem_sys,
        temp,
        stability_cutoff=stability_cutoff,
        open_element=open_element,
        chempot=chempot,
        formulas_to_include=formulas_to_include
    )
    run_locally(flow, store=store.store)    