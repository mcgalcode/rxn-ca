from typing import Dict
import numpy as np

from pylattica.core.simulation import Simulation

from ..phases.solid_phase_set import SolidPhaseSet
from ..analysis.reaction_step_analyzer import ReactionStepAnalyzer
from ..core.constants import VOLUME, VOL_MULTIPLIER, GASES_CONSUMED, GASES_EVOLVED, MELTED_AMTS
from .constants import VOLUME_TOLERANCE_FRAC, VOLUME_TOLERANCE_ABS
from pylattica.structures.square_grid import DiscreteGridSetup, PseudoHexagonalNeighborhoodBuilder2D, PseudoHexagonalNeighborhoodBuilder3D
from pylattica.core import AsynchronousRunner, Simulation
from pylattica.discrete.discrete_step_analyzer import DiscreteStepAnalyzer

from .volume_tuning_controller import VolumeTuningController
from .phase_growth_controller import PhaseGrowthController

from tabulate import tabulate

class ReactionPreparer():
    """Sets up SimulationStates for running the reaction automaton.
    The main purpose of this class is to handle converting phase ratios
    (which are interpreted as molar quantities) to volume ratios
    """    

    def __init__(self, phase_set: SolidPhaseSet, dim: int = 3):
        self.dim = dim
        self.volumes = phase_set.volumes
        self.phase_set: SolidPhaseSet = phase_set

    def prepare_reaction(self, 
                         phase_mol_ratios: Dict[str, float],
                         size: int = 15,
                         volume_multiplier: float = 1.0
        ) -> Simulation:

        total_vol = size ** self.dim
        vol_per_nuc_site = 5
        num_sites = int(total_vol / vol_per_nuc_site) / 4

        print(f'Nucleating grains at {num_sites} sites')

        total_vol_available = size ** self.dim * volume_multiplier
        volume_ratios = self.phase_set.mole_amts_to_vols(phase_mol_ratios)
        
        total_vol_ratio = sum(volume_ratios.values())
        normalized_vol_ratios = { p: vol / total_vol_ratio for p, vol in volume_ratios.items() }
        desired_phase_vols = { p: vol * total_vol_available for p, vol in normalized_vol_ratios.items() }

        ## 1. SET UP NUCLEATION SITES BASED ON VOLUME RATIOS
        print("Reactant Phases: ", desired_phase_vols)
        print("Volume Ratios: ", normalized_vol_ratios)
        print("Using volume multiplier: ", volume_multiplier)
        if self.dim == 2:
            nb_spec = PseudoHexagonalNeighborhoodBuilder2D()
        else:
            nb_spec = PseudoHexagonalNeighborhoodBuilder3D()

        setup = DiscreteGridSetup(self.phase_set, dim=self.dim)
        simulation = setup.setup_random_sites(
            size,
            num_sites_desired=num_sites,
            background_spec=self.phase_set.FREE_SPACE,
            nuc_amts=normalized_vol_ratios,
            buffer = 1,
        )

        discrete_analzyer = DiscreteStepAnalyzer()
        rxn_step_analyzer = ReactionStepAnalyzer(self.phase_set)            

        simulation.state.set_general_state({
            MELTED_AMTS: {},
            VOL_MULTIPLIER: volume_multiplier,
            GASES_EVOLVED: {},
            GASES_CONSUMED: {}
        })
        
        for site_id in simulation.state.site_ids():
            simulation.state.set_site_state(site_id, { VOLUME: 1.0 })
        
        ## 2. Grow particles away from nucleation sites according to desired amounts

        controller = PhaseGrowthController(
            self.phase_set,
            simulation.structure,
            desired_phase_vols=desired_phase_vols,
            nb_builder=nb_spec,
            background_phase=self.phase_set.FREE_SPACE,
        )

        empty_count = discrete_analzyer.cell_count(simulation.state, self.phase_set.FREE_SPACE)

        runner = AsynchronousRunner()

        print("\nGrowing grains to fill empty space:")
        while empty_count > 0:
            print(f'Filling remaining {empty_count} vacant cells...')
            res = runner.run(simulation.state, controller, num_steps=3 * int(size**3 / 2))
            simulation = Simulation(res.last_step, simulation.structure)
            empty_count = discrete_analzyer.cell_count(simulation.state, self.phase_set.FREE_SPACE)

        assert simulation.state.get_general_state().get(VOL_MULTIPLIER) == volume_multiplier
        print(f"Simulation is filled! Total volume is {rxn_step_analyzer.get_total_volume(simulation.state)}")

        # At this stage, we check to see how close we are
        # if we are not close, we randomly adjust cells to get closer to the 
        # desired state

        def _check_closeness(ideal_amts, actual_amts):
            close = True

            table = []
            for phase, ideal_amt in ideal_amts.items():
                actual = actual_amts.get(phase, 0)
                abs_diff = np.abs(ideal_amt - actual)
                frac_diff = abs_diff / ideal_amt
                if abs_diff > VOLUME_TOLERANCE_ABS or frac_diff > VOLUME_TOLERANCE_FRAC:
                    close = False
                
                table.append([phase, ideal_amt, actual, abs_diff, frac_diff])

            print(tabulate(table, headers=["Phase", "Target Amt.", "Current Amt.", "Abs. Diff.", "Frac. Diff."])) 
                        
            return close
    
        def _assess_convergence(analyzer, simulation):
    
            print("Assessing convergence using volumes...")
            actual = analyzer.get_all_absolute_phase_volumes(simulation.state, include_melted=False)
            ideal = desired_phase_vols

            close = _check_closeness(ideal, actual)

            return close
        

        done = _assess_convergence(rxn_step_analyzer, simulation)

        if done:
            print("Converged!")
            
            return simulation

        tuner_controller = VolumeTuningController(
            self.phase_set,
            desired_phase_vols,
        )

        print("\n")
        print(f"Tweaking volumes to tune amounts using volume multiplier {volume_multiplier}...")
        print("Targeting desired phase volumes:")
        print(desired_phase_vols)

        close = False
        tries_remaining = 15
        while not close and tries_remaining > 0:
            print("\nAmounts not converged, continuing tuning process...")
            res = runner.run(simulation.state, tuner_controller, num_steps=5 * int(size**3))
            simulation = Simulation(res.last_step, simulation.structure)

            close = _assess_convergence(rxn_step_analyzer, simulation)
            tries_remaining = tries_remaining - 1

        if tries_remaining == 0 and not close:
            raise RuntimeError("Could not generate starting state with correct molar amounts!")
        
        print("Converged!")

        return simulation