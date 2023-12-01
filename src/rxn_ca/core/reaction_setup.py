from typing import Dict, List
import numpy as np
import random

from pylattica.core.simulation import Simulation

from ..phases.solid_phase_set import SolidPhaseSet
from ..analysis.reaction_step_analyzer import ReactionStepAnalyzer
from .constants import VOLUME
from pylattica.structures.square_grid import DiscreteGridSetup, PseudoHexagonalNeighborhoodBuilder2D, PseudoHexagonalNeighborhoodBuilder3D, GrowthSetup
from pylattica.core import AsynchronousRunner, SimulationState, Simulation
from pylattica.core import BasicController
from pylattica.core.neighborhood_builders import NeighborhoodBuilder
from pylattica.core.periodic_structure import PeriodicStructure
from pylattica.core.simulation_state import SimulationState
from pylattica.discrete import PhaseSet
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY, VACANT
from pylattica.structures.square_grid.neighborhoods import MooreNbHoodBuilder

default_ratios = [1, 1, 1]


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
        self.desired_phase_vols = desired_phase_vols

        if nb_builder is None:
            self.nb_builder = MooreNbHoodBuilder(1, dim=periodic_struct.dim)
        else:
            self.nb_builder = nb_builder

        self.nb_graph = self.nb_builder.get(periodic_struct)

    def get_state_update(self, site_id: int, prev_state: SimulationState):
        curr_state = prev_state.get_site_state(site_id)
        if curr_state[DISCRETE_OCCUPANCY] == self.background_phase:
            counts = {}
            for nb_id in self.nb_graph.neighbors_of(site_id):
                nb_phase = prev_state.get_site_state(nb_id)[DISCRETE_OCCUPANCY]
                if nb_phase != self.background_phase:
                    if nb_phase not in counts:
                        counts[nb_phase] = 1
                    else:
                        counts[nb_phase] += 1

            if len(counts) > 0:
                max_count = 0
                max_spec = None
                molar_bdown = self.analyzer.phase_volume_fractions(prev_state)
                for phase, count in counts.items():
                    if count > max_count and (molar_bdown[phase] <= self.desired_phase_vols[phase] or random.random() < 0.8):
                        max_spec = phase
                        max_count = count

                if max_spec is not None:
                    return {DISCRETE_OCCUPANCY: max_spec}
                else:
                    return {}
            else:
                return {}
        else:
            return {}


class ReactionSetup(DiscreteGridSetup):
    """Sets up SimulationStates for running the reaction automaton.
    The main purpose of this class is to handle converting phase ratios
    (which are interpreted as molar quantities) to volume ratios
    """    

    def __init__(self, phase_set: SolidPhaseSet, dim: int = 2):
        super().__init__(phase_set, dim)
        self.dim = dim
        self.volumes = phase_set.volumes

    def setup_growth(self, 
            size: int,
            num_seeds: int,
            reactant_phases: List[str],
            phase_vol_ratios: List[float] = default_ratios,
            phase_mol_ratios: List[float] = None,
            buffer: int = 2
        ) -> SimulationState:
        if phase_mol_ratios is not None:
            volume_scaled_ratios = [self.volumes[phase] * phase_mol_ratios[idx] for idx, phase in enumerate(reactant_phases)]
        else:
            volume_scaled_ratios = phase_vol_ratios
        
        volume_scaled_ratios = np.array(volume_scaled_ratios) / np.sum(volume_scaled_ratios)
        print("Reactant Phases: ", reactant_phases)
        print("Volume Ratios: ", volume_scaled_ratios)
        if self.dim == 2:
            nb_spec = PseudoHexagonalNeighborhoodBuilder2D()
        else:
            nb_spec = PseudoHexagonalNeighborhoodBuilder3D()

        setup = DiscreteGridSetup(self.phase_set, dim=self.dim)
        simulation = setup.setup_random_sites(
            size,
            num_sites_desired=num_seeds,
            background_spec=self.phase_set.FREE_SPACE,
            nuc_species=reactant_phases,
            nuc_ratios=volume_scaled_ratios,
            buffer=buffer,
        )

        analyzer = ReactionStepAnalyzer(self.phase_set)            

        for site_id in simulation.state.site_ids():
            simulation.state.set_site_state(site_id, { VOLUME: 1.0 })

        desired_phase_vols = {}
        for phase, amt in zip(reactant_phases, volume_scaled_ratios):
            desired_phase_vols[phase] = amt

        controller = ReactionSetupController(
            self.phase_set,
            simulation.structure,
            desired_phase_vols=desired_phase_vols,
            nb_builder=nb_spec,
            background_phase=self.phase_set.FREE_SPACE,
        )

        empty_count = analyzer.cell_count(simulation.state, self.phase_set.FREE_SPACE)
        while empty_count > 0:
            print(f'Filling remaining {empty_count} vacant cells...')
            runner = AsynchronousRunner()
            res = runner.run(simulation.state, controller, num_steps=int(size**3 / 2))
            simulation = Simulation(res.last_step, simulation.structure)
            empty_count = analyzer.cell_count(simulation.state, self.phase_set.FREE_SPACE)


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

    def setup_random_mixture(self,
        side_length: int,
        grain_size: int,
        phases: List[str],
        weights: List[float] = None
    ) -> SimulationState:
        volume_scaled_ratios = [self.volumes[phase] * weights[idx] for idx, phase in enumerate(phases)]
        return super().setup_random_mixture(side_length, grain_size, phases, volume_scaled_ratios)