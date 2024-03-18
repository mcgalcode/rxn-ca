from ..phases.solid_phase_set import SolidPhaseSet
from ..setup import ReactionPreparer
from ..setup.noise_setup import SetupRandomNoise

from pylattica.core import Simulation
from typing import Dict

def setup_reaction(
        phases: SolidPhaseSet,
        precursor_mole_ratios: Dict,
        size: int = 15,
        vol_multiplier = 1.0,
    ) -> Simulation:

    preparer = ReactionPreparer(phases, dim=3)
    sim = preparer.prepare_reaction(
        phase_mol_ratios=precursor_mole_ratios,
        size=size,
        volume_multiplier=vol_multiplier
    )
    
    return sim

def setup_noise_reaction(
        phases: SolidPhaseSet,
        precursor_mole_ratios: Dict,
        size: int = 15,
        packing_fraction = 1.0
    ):
    return SetupRandomNoise(phases).setup(precursor_mole_ratios, size, packing_efficiency=packing_fraction)