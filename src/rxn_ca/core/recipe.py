from monty.json import MSONable
from dataclasses import dataclass, field
from .heating import HeatingSchedule
from typing import Dict, List
from ..utilities.helpers import format_chem_sys
from ..reactions.scorers import TammanHuttigScoreErf, TammanHuttigScoreExponential, TammanHuttigScoreSoftplus, VirScorer
import json

from enum import Enum

class ScoreTypes(str, Enum):

    TAMMAN_HUTTIG_SOFTPLUS_GIBBS_ERF = "TAMMAN_HUTTIG_SOFTPLUS_GIBBS_ERF"
    TAMMAN_HUTTIG_SOFTPLUS_GIBBS_SOFTPLUS = "TAMMAN_HUTTIG_SOFTPLUS_GIBBS_SOFTPLUS"
    TAMMAN_HUTTIG_EXP_GIBBS_SOFTPLUS = "TAMMAN_HUTTIG_EXP_GIBBS_SOFTPLUS"
    VIR_SCORER = "VIR_SCORER"


_SCORE_TYPE_MAP = {
    ScoreTypes.TAMMAN_HUTTIG_EXP_GIBBS_SOFTPLUS: TammanHuttigScoreExponential,
    ScoreTypes.TAMMAN_HUTTIG_SOFTPLUS_GIBBS_ERF: TammanHuttigScoreErf,
    ScoreTypes.TAMMAN_HUTTIG_SOFTPLUS_GIBBS_SOFTPLUS: TammanHuttigScoreSoftplus,
    ScoreTypes.VIR_SCORER: VirScorer
}


@dataclass
class ReactionRecipe(MSONable):

    heating_schedule: HeatingSchedule
    reactant_amounts: Dict[str, float]
    chem_sys: str
    simulation_size: int = 15
    num_realizations: int = 3
    exclude_phases: List[str] = field(default_factory=list)
    exclude_theoretical: bool = True
    exclude_pure_elements: bool = False
    score_type: str = ScoreTypes.TAMMAN_HUTTIG_SOFTPLUS_GIBBS_ERF
    additional_gas_phases: List[str] = field(default_factory=list)
    exact_phase_set: List[str] = None
    
    def __post_init__(self):
        self.chem_sys = format_chem_sys(self.chem_sys)

    def to_file(self, fname: str):
        with open(fname, 'w+') as f:
            f.write(self.to_json())

    def get_score_class(self):
        return _SCORE_TYPE_MAP[self.score_type]

    @classmethod
    def from_file(self, fname: str):
        with open(fname, 'rb') as f:
            return self.from_dict(json.loads(f.read()))