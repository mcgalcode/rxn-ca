from pylattica.core import Simulation, Runner

from .reaction_controller import ReactionController
from .solid_phase_set import SolidPhaseSet
from .reaction_result import ReactionResult
from .reaction_setup import ReactionSetup
from .reaction_simulation import ReactionSimulation

from typing import Dict
from ..reactions import ScoredReactionSet, ArrheniusScore
from ..computing import AutomatonStore, enumerate_flow
from jobflow.managers.local import run_locally

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
    runner = Runner(is_async=True)
    controller = ReactionController(
        simulation.structure,
        simulation.reactions,
        free_species=free_species,
        inertia=inertia
    )
    result = runner.run(simulation.state, controller, num_steps, structure=simulation.structure, verbose=verbose)
    return result

def get_arrhenius_rxns(chemsys, temp, **kwargs):
    store = AutomatonStore()
    scorer = ArrheniusScore(temp)
    rxns = store.get_scored_rxns(chemsys, temp, scorer, **kwargs)
    return rxns

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