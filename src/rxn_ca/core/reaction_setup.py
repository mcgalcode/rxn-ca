from typing import Dict, List
import numpy as np
import random

from pylattica.core.simulation import Simulation

from ..phases.solid_phase_set import SolidPhaseSet
from ..analysis.reaction_step_analyzer import ReactionStepAnalyzer
from .constants import VOLUME, VOL_MULTIPLIER
from pylattica.structures.square_grid import DiscreteGridSetup, PseudoHexagonalNeighborhoodBuilder2D, PseudoHexagonalNeighborhoodBuilder3D, GrowthSetup
from pylattica.core import AsynchronousRunner, SimulationState, Simulation
from pylattica.core import BasicController
from pylattica.core.neighborhood_builders import NeighborhoodBuilder
from pylattica.core.periodic_structure import PeriodicStructure
from pylattica.core.simulation_state import SimulationState
from pylattica.discrete import PhaseSet
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY, VACANT
from pylattica.structures.square_grid.neighborhoods import MooreNbHoodBuilder

class ReactionSetupController(BasicController):

    def __init__(
        self,
        phase_set: PhaseSet,
        periodic_struct: PeriodicStructure,
        desired_phase_vols: Dict,
        background_phase: str = VACANT,
        nb_builder: NeighborhoodBuilder = None,
    ) -> None:
        
        self.background_phase = background_phase
        self.phase_set: PhaseSet = phase_set
        self.analyzer = ReactionStepAnalyzer(self.phase_set)
        self.desired_phase_vols_abs = desired_phase_vols

        self.known_empty_ids = periodic_struct.site_ids.copy()

        if nb_builder is None:
            self.nb_builder = MooreNbHoodBuilder(1, dim=periodic_struct.dim)
        else:
            self.nb_builder = nb_builder

        self.nb_graph = self.nb_builder.get(periodic_struct)
    
    def get_random_site(self, _):
        if len(self.known_empty_ids) > 0:
            return random.choice(self.known_empty_ids)
        else:
            return 1

    def get_state_update(self, site_id: int, prev_state: SimulationState):
        
        curr_state = prev_state.get_site_state(site_id)
        if curr_state[DISCRETE_OCCUPANCY] == self.background_phase:
            nb_phases = set()
            for nb_id in self.nb_graph.neighbors_of(site_id):
                nb_phase = prev_state.get_site_state(nb_id)[DISCRETE_OCCUPANCY]
                if nb_phase != self.background_phase:
                    nb_phases.add(nb_phase)

            if len(nb_phases) > 0:
                chosen_phase = None
                volume_breakdown = self.analyzer.phase_volumes(prev_state, include_melted=False)

                deficient_candidates = []

                for phase in nb_phases:
                    diff = volume_breakdown[phase] - self.desired_phase_vols_abs[phase]
                    deficient_candidates.append((phase, diff))
                
                most_deficient = min(deficient_candidates, key = lambda c: c[1])[0]

                if len(deficient_candidates) > 0:
                    chosen_phase = most_deficient
                    self.known_empty_ids.remove(site_id)
                    return {DISCRETE_OCCUPANCY: chosen_phase}
                else:
                    return {}
            else:
                return {}
        else:
            if site_id in self.known_empty_ids:
                self.known_empty_ids.remove(site_id)
            return {}


class ReactionTunerController(BasicController):

    INC_UP = 1.01
    INC_DOWN = 0.99
    MAX_ADJUSTMENTS = 50

    TOL_FRAC = 0.01

    def __init__(
        self,
        phase_set: PhaseSet,
        ideal_vol_amts: Dict,
        vol_scale: float
    ) -> None:
        
        self.phase_set = phase_set
        self.analyzer = ReactionStepAnalyzer(self.phase_set)
        self.ideal_vol_amts = ideal_vol_amts

        self.vol_scale = vol_scale
        self.lower_vol_lim = vol_scale * (self.INC_DOWN ** self.MAX_ADJUSTMENTS)
        self.upper_vol_lim = vol_scale * (self.INC_UP ** self.MAX_ADJUSTMENTS)

    def get_random_site(self, prev_state: SimulationState):
        
        curr_amt = self.analyzer.phase_volumes(prev_state, include_melted=False)

        deficient_phases = []

        for phase, ideal_vol in self.ideal_vol_amts.items():
            curr_vol = curr_amt.get(phase)
            if np.abs(curr_vol - ideal_vol) / ideal_vol > self.TOL_FRAC:
                deficient_phases.append(phase)
        
        criteria = [
            lambda s: s[DISCRETE_OCCUPANCY] in deficient_phases,
        ]

        valid_sites = self.analyzer.get_sites(
            prev_state,
            state_criteria = criteria
        )

        random.shuffle(valid_sites)
        return valid_sites

    def get_state_update(self, site_id: int, prev_state: SimulationState):
        
        curr_state = prev_state.get_site_state(site_id)
        curr_phase = curr_state[DISCRETE_OCCUPANCY]
        curr_vol = curr_state[VOLUME]

        ideal_amt = self.ideal_vol_amts.get(curr_phase)

        curr_amt = self.analyzer.phase_volumes(prev_state, include_melted=False).get(curr_phase)
        if curr_amt < ideal_amt:
            new_vol = curr_vol * self.INC_UP
        elif curr_amt > ideal_amt:
            new_vol = curr_vol * self.INC_DOWN
        else:
            new_vol = curr_vol
        
        return {VOLUME: new_vol}

class ReactionSetup(DiscreteGridSetup):
    """Sets up SimulationStates for running the reaction automaton.
    The main purpose of this class is to handle converting phase ratios
    (which are interpreted as molar quantities) to volume ratios
    """    

    def __init__(self, phase_set: SolidPhaseSet, dim: int = 2):
        super().__init__(phase_set, dim)
        self.dim = dim
        self.volumes = phase_set.volumes
        self.phase_set: SolidPhaseSet = phase_set

    def setup_growth(self, 
            size: int,
            num_seeds: int,
            phase_mol_ratios: Dict = None,
            phase_mol_amts: Dict = None,
            buffer: int = None,
            volume_scale: float = 1.0,
            volume_multiplier: float = 1.0
        ) -> SimulationState:

        # handle mole ratios
        if phase_mol_ratios is not None:
            total_vol_available = size ** self.dim * volume_scale
            volume_scaled_ratios = self.phase_set.mole_amts_to_vols(phase_mol_ratios)
            
            total_vol_ratio = sum(volume_scaled_ratios.values())
            normalized_vol_ratios = { p: vol / total_vol_ratio for p, vol in volume_scaled_ratios.items() }
            desired_phase_vols = { p: vol * total_vol_available for p, vol in normalized_vol_ratios.items() }
        # handle absolute precursor mole amts
        elif phase_mol_amts is not None:
            desired_phase_vols = self.phase_set.mole_amts_to_vols(phase_mol_amts)
            total_vol = sum(desired_phase_vols.values())
            normalized_vol_ratios = { p: vol / total_vol for p, vol in desired_phase_vols.items() }

        ## 1. SET UP NUCLEATION SITES BASED ON VOLUME RATIOS
        print("Reactant Phases: ", desired_phase_vols)
        print("Volume Ratios: ", normalized_vol_ratios)
        print(f"Volume Scale: {volume_scale}")
        if self.dim == 2:
            nb_spec = PseudoHexagonalNeighborhoodBuilder2D()
        else:
            nb_spec = PseudoHexagonalNeighborhoodBuilder3D()

        setup = DiscreteGridSetup(self.phase_set, dim=self.dim)
        simulation = setup.setup_random_sites(
            size,
            num_sites_desired=num_seeds,
            background_spec=self.phase_set.FREE_SPACE,
            nuc_amts=normalized_vol_ratios,
            buffer=buffer,
        )

        analyzer = ReactionStepAnalyzer(self.phase_set)            

        simulation.state.set_general_state({ VOL_MULTIPLIER: volume_multiplier })
        
        for site_id in simulation.state.site_ids():
            simulation.state.set_site_state(site_id, { VOLUME: volume_scale })
        
        ## 2. Grow particles away from nucleation sites according to desired amounts

        controller = ReactionSetupController(
            self.phase_set,
            simulation.structure,
            desired_phase_vols=desired_phase_vols,
            nb_builder=nb_spec,
            background_phase=self.phase_set.FREE_SPACE,
        )

        empty_count = analyzer.cell_count(simulation.state, self.phase_set.FREE_SPACE)

        runner = AsynchronousRunner()

        while empty_count > 0:
            print(f'Filling remaining {empty_count} vacant cells...')
            res = runner.run(simulation.state, controller, num_steps=3 * int(size**3 / 2))
            simulation = Simulation(res.last_step, simulation.structure)
            empty_count = analyzer.cell_count(simulation.state, self.phase_set.FREE_SPACE)


        # At this stage, we check to see how close we are
        # if we are not close, we randomly adjust cells to get closer to the 
        # desired state

        def _check_closeness(ideal_amts, actual_amts):
            close = True

            for phase, ideal_amt in ideal_amts.items():
                actual = actual_amts.get(phase, 0)
                if np.abs(ideal_amt - actual) / ideal_amt > 0.01:
                    print(f"{phase}: {ideal_amt}, {actual}, {np.abs(ideal_amt - actual) / ideal_amt}")
                    close = False
            
            return close

        tuner_controller = ReactionTunerController(
            self.phase_set,
            desired_phase_vols,
            vol_scale=volume_scale
        )

        close = False
        tries_remaining = 15
        while not close and tries_remaining > 0:
            print("Tweaking volumes to tune amounts...")
            res = runner.run(simulation.state, tuner_controller, num_steps=5 * int(size**3))
            simulation = Simulation(res.last_step, simulation.structure)

            if phase_mol_amts is not None:
                actual = analyzer.molar_breakdown(simulation.state, include_melted=False)
                ideal = phase_mol_amts
                print("Assessing convergence using molar amounts...")
            elif phase_mol_ratios is not None:
                print("Assessing convergence using volume ratios...")
                actual = analyzer.phase_volume_fractions(simulation.state, include_melted=False)
                ideal = normalized_vol_ratios

            close = _check_closeness(ideal, actual)
            tries_remaining = tries_remaining - 1

        
        return simulation

    def setup_interface(self, size: int, p1: str, p2: str) -> Simulation:
        simulation = super().setup_interface(size, p1, p2)
        for site_id in simulation.state.site_ids():
            simulation.state.set_site_state(site_id, { VOLUME: 1.0 })
        return simulation
    
    def setup_random_sites(self,
        size: int,
        num_sites_desired: int,
        background_spec: str,
        nuc_species: List[str],
        nuc_ratios: List[float] = None,
        buffer = 2
    ) -> SimulationState:
        volume_scaled_ratios = [self.volumes[phase] * nuc_ratios[idx] for idx, phase in enumerate(nuc_species)]
        return super().setup_random_sites(
            size,
            num_sites_desired=num_sites_desired,
            background_spec=background_spec,
            nuc_species=nuc_species,
            nuc_ratios=volume_scaled_ratios,
            buffer=buffer
        )