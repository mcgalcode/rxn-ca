from typing import Dict

from pylattica.structures.square_grid import DiscreteGridSetup
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY
from pylattica.core import Simulation

from ..core.constants import VOLUME, MELTED_AMTS, VOL_MULTIPLIER, GASES_CONSUMED, GASES_EVOLVED
from ..phases import SolidPhaseSet
from random import shuffle

import copy
class SetupRandomNoise():

    def __init__(self, phases: SolidPhaseSet, dim: int = 3):
        self.phase_set = phases
        self.dim = dim

    def setup(self,
            phase_mol_ratios: Dict[str, float],
            size: int = 15,
            packing_efficiency = 0.97
    ):
        total_vol = size ** self.dim * packing_efficiency
        volume_ratios = self.phase_set.mole_amts_to_vols(phase_mol_ratios)
        
        total_vol_ratio = sum(volume_ratios.values())
        normalized_vol_ratios = { p: vol / total_vol_ratio for p, vol in volume_ratios.items() }
        desired_phase_vols = { p: round(vol * total_vol) for p, vol in normalized_vol_ratios.items() }
        setup = DiscreteGridSetup(self.phase_set, dim=self.dim)

        struct = setup.build_structure(size)
        state = setup.setup_solid_phase(struct, self.phase_set.FREE_SPACE)

        cell_occs = [k for k, v in desired_phase_vols.items() for _ in range(v)]
        shuffle(cell_occs)

        site_ids = struct.site_ids
        shuffle(site_ids)

        for i, occ in enumerate(cell_occs):
            state.set_site_state(site_ids[i], {
                DISCRETE_OCCUPANCY: occ,
                VOLUME: 1.0
            })

        for sid in site_ids:
            state.set_site_state(sid, {
                VOLUME: 1.0
            })

        simulation = Simulation(state, struct)

        simulation.state.set_general_state({
            MELTED_AMTS: {},
            VOL_MULTIPLIER: 1.0,
            GASES_EVOLVED: {},
            GASES_CONSUMED: {}
        })

        return simulation
