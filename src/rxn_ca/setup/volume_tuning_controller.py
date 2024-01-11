from typing import Dict
import numpy as np
import random

from ..analysis.reaction_step_analyzer import ReactionStepAnalyzer
from ..core.constants import VOLUME
from pylattica.core import SimulationState
from pylattica.core import BasicController
from pylattica.core.simulation_state import SimulationState
from pylattica.discrete import PhaseSet, DiscreteStepAnalyzer
from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY

class VolumeTuningController(BasicController):

    INC_UP = 1.01
    INC_DOWN = 0.99
    MAX_ADJUSTMENTS = 50

    TOL_FRAC = 0.01

    def __init__(
        self,
        phase_set: PhaseSet,
        ideal_vol_amts: Dict,
    ) -> None:
        
        self.phase_set = phase_set
        self.analyzer = ReactionStepAnalyzer(self.phase_set)
        self.discrete_analyzer = DiscreteStepAnalyzer()
        self.ideal_vol_amts = ideal_vol_amts

    def get_random_site(self, prev_state: SimulationState):
        
        curr_amt = self.analyzer.get_all_absolute_phase_volumes(prev_state, include_melted=False)

        deficient_phases = []

        for phase, ideal_vol in self.ideal_vol_amts.items():
            curr_vol = curr_amt.get(phase)
            if np.abs(curr_vol - ideal_vol) / ideal_vol > self.TOL_FRAC:
                deficient_phases.append(phase)
        
        criteria = [
            lambda s: s[DISCRETE_OCCUPANCY] in deficient_phases,
        ]

        valid_sites = self.discrete_analyzer.get_sites(
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

        curr_amt = self.analyzer.get_all_absolute_phase_volumes(prev_state, include_melted=False).get(curr_phase)
        if curr_amt < ideal_amt:
            new_vol = curr_vol * self.INC_UP
        elif curr_amt > ideal_amt:
            new_vol = curr_vol * self.INC_DOWN
        else:
            new_vol = curr_vol
        
        return {VOLUME: new_vol}