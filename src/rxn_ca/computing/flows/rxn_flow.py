from jobflow import Maker, job

from rxn_ca import setup_reaction, get_scored_rxns
from ...core.multi_stage_runner import run_multi
from ...core.recipe import ReactionRecipe
from ...reactions import ReactionLibrary

from ..schemas.ca_result_schema import RxnCAResultDoc

from pylattica.core import Simulation

from dataclasses import dataclass


@dataclass
class SingleRxnMaker(Maker):

    name: str = "Rxn Sim. Single"

    @job(data="results")
    def make(self,
             recipe: ReactionRecipe,
             output_fname: str = None,
             reactions_filename: str = None,
             initial_simulation_filename: str = None):
        
        print('\n\n\n')
        print(
        "██╗░░░░░███████╗████████╗  ██╗████████╗  ░█████╗░░█████╗░░█████╗░██╗░░██\n"
        "██║░░░░░██╔════╝╚══██╔══╝  ██║╚══██╔══╝  ██╔══██╗██╔══██╗██╔══██╗██║░██╔╝\n"
        "██║░░░░░█████╗░░░░░██║░░░  ██║░░░██║░░░  ██║░░╚═╝██║░░██║██║░░██║█████═╝░\n"
        "██║░░░░░██╔══╝░░░░░██║░░░  ██║░░░██║░░░  ██║░░██╗██║░░██║██║░░██║██╔═██╗░\n"
        "███████╗███████╗░░░██║░░░  ██║░░░██║░░░  ╚█████╔╝╚█████╔╝╚█████╔╝██║░╚██╗\n"
        "╚══════╝╚══════╝░░░╚═╝░░░  ╚═╝░░░╚═╝░░░  ░╚════╝░░╚════╝░░╚════╝░╚═╝░░╚═╝\n")

        print('\n\n\n')

        print("================= RETRIEVING AND SCORING REACTIONS =================")

        if reactions_filename is not None:
            print(f'Reading reactions from {reactions_filename}...')
            rxn_lib = ReactionLibrary.from_file(reactions_filename)
        else:
            rxn_lib = get_scored_rxns(
                recipe.chem_sys,
                heating_sched=recipe.heating_schedule,
                exclude_phases=recipe.exclude_phases,
                exclude_theoretical=recipe.exclude_theoretical,
                scorer_class=recipe.get_score_class()
            )

        print()
        print()
        print()

        print("================= SETTING UP SIMULATION =================")

        if initial_simulation_filename is None:
            initial_simulation = setup_reaction(
                rxn_lib.phases,
                recipe.simulation_size,
                phase_mole_ratios = recipe.reactant_amounts,
                particle_size=recipe.particle_size
            )
        else:
            initial_simulation = Simulation.from_file(initial_simulation_filename)

        print(f'================= RUNNING SIMULATION w/ 1 REALIZATION =================')


        result = run_multi(
            initial_simulation,
            rxn_lib,
            recipe.heating_schedule
        )

        result_doc = RxnCAResultDoc(
            recipe=recipe,
            results=[result],
            reaction_library=rxn_lib,
        )

        if output_fname is not None:
            print(f'================= SAVING RESULTS to {output_fname} =================')
            result_doc.to_file(output_fname)

        return result_doc

