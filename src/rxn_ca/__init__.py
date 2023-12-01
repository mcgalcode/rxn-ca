
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

def setup(phases: SolidPhaseSet, phase_mixture: Dict, size: int, dim: int = 3, particle_size = 1):
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


def get_scored_rxns(chemsys: str,
                    heating_sched: HeatingSchedule,
                    scorer_class: BasicScore = TammanHuttigScoreErf,
                    exclude_pure_elements: bool = False,
                    exclude_theoretical: bool = True,
                    exclude_phases: List[str]= []):
    store = AutomatonStore()
    base_rxn_set = store.get_raw_rxns(
        chemsys,
        temp=300
    )

    lib = ReactionLibrary(base_rxn_set)

    sched_temps = heating_sched.all_temps

    new_temps = list(set(sched_temps) - set(lib.temps))
    if len(new_temps) == 0:
        return
    
    scorers = [scorer_class(temp=t, phase_set=lib.phases) for t in new_temps]

    rsets: List[ReactionSet] = []
    for t in tqdm(new_temps, desc="Getting reaction energies at temperatures..."):
        existing_rxns = store.get_raw_rxns(
            chemsys,
            temp=t
        )
        if existing_rxns is not None:
            rsets.append(existing_rxns)
        else:
            rsets.append(base_rxn_set.set_new_temperature(t))

    for t, scorer, rset in zip(new_temps, scorers, rsets):
        scored_rxns: List[ScoredReaction] = score_rxns(rset, scorer, phase_set=lib.phases)
        rxn_set = ScoredReactionSet(scored_rxns, lib.phases, identity_score=3)
        if exclude_pure_elements:
            print("Excluding reactions involving pure elements...")
            rxn_set = rxn_set.exclude_pure_els()
        
        if exclude_theoretical:
            print("Excluding reactions involving theoretical compounds...")
            rxn_set = rxn_set.exclude_theoretical(lib.phases)

        if len(exclude_phases) > 0:
            print(f"Excluding reactions including {exclude_phases}")
            rxn_set = rxn_set.exclude_phases(exclude_phases)

        lib.add_rxns_at_temp(rxn_set, t)

    return lib



def enumerate_rxns(chem_sys,
                   base_temp = 300,
                   stability_cutoff=0.1,
                   open_element=None,
                   chempot=None,
                   other_temps = [],
                   formulas_to_include: List = []
    ):
    store = AutomatonStore()
    flow = enumerate_flow(
        chem_sys,
        base_temp=base_temp,
        stability_cutoff=stability_cutoff,
        open_element=open_element,
        chempot=chempot,
        other_temps = other_temps,
        formulas_to_include=formulas_to_include
    )
    if flow is not None:
        run_locally(flow, store=store.store)    