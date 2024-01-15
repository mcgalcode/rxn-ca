from ..phases.solid_phase_set import SolidPhaseSet
from ..setup import ReactionPreparer

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