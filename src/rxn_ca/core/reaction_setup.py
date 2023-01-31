from typing import Dict, List

from .solid_phase_set import SolidPhaseSet
from pylattica.square_grid import DiscreteGridSetup, PseudoHexagonalNeighborhoodBuilder2D, PseudoHexagonalNeighborhoodBuilder3D
from pylattica.models.growth import GrowthController
from pylattica.core import Runner, SimulationState, Simulation

class ReactionSetup(DiscreteGridSetup):
    """Sets up SimulationStates for running the reaction automaton.
    The main purpose of this class is to handle converting phase ratios
    (which are interpreted as molar quantities) to volume ratios
    """    

    def __init__(self, phase_set: SolidPhaseSet, dim: int = 2, volumes: Dict = None):
        super().__init__(phase_set, dim)
        self.dim = dim
        self.volumes = volumes

    def setup_growth(self, size: int, num_seeds: int, reactant_phases: List[str], phase_ratios: List[float]) -> SimulationState:
        volume_scaled_ratios = [self.volumes[phase] * phase_ratios[idx] for idx, phase in enumerate(reactant_phases)]
        
        if self.dim == 2:
            nb_spec = PseudoHexagonalNeighborhoodBuilder2D()
        else:
            nb_spec = PseudoHexagonalNeighborhoodBuilder3D()

        setup = DiscreteGridSetup(reactant_phases, self.dim)
        periodic_state = setup.setup_random_sites(
            size,
            num_seeds,
            self.phase_set.FREE_SPACE,
            reactant_phases,
            volume_scaled_ratios
        )

        controller = GrowthController(
            reactant_phases,
            periodic_state.structure,
            neighborhood_spec = nb_spec,
            background_phase = self.phase_set.FREE_SPACE
        )
        
        runner = Runner(parallel=True)
        res = runner.run(periodic_state.state, controller, num_steps = size)
        return Simulation(res.last_step, periodic_state.structure)

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