from typing import Dict
from ..core import ReactionSetup
from ..core import ReactionSimulation
from ..core import ReactionController
from ..reactions import ScoredReactionSet

from ..computing import AutomatonStore, enumerate_flow

from jobflow.managers.local import run_locally
from pylattica.core import AsynchronousRunner

import typing

def setup(rxns: ScoredReactionSet, phases: Dict, size: int, num_particles: int, dim: int = 2):
    phase_names = []
    phase_amts = []

    for name, amt in phases.items():
        phase_names.append(name)
        phase_amts.append(amt)

    setup = ReactionSetup(rxns.phases, dim=dim)
    sim = setup.setup_growth(size, num_particles, phase_names, phase_mol_ratios=phase_amts)
    return ReactionSimulation(rxns, sim)

def run(simulation: ReactionSimulation, num_steps: int, free_species=None, verbose=True, inertia=0):
    runner = AsynchronousRunner()
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
        verbose=verbose
    )
    return result

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