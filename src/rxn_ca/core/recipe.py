from monty.json import MSONable
from dataclasses import dataclass, field
from .heating import HeatingSchedule
from typing import Dict, List
from ..computing.utils.functions import format_chem_sys
import json

@dataclass
class ReactionRecipe(MSONable):

    heating_schedule: HeatingSchedule
    reactant_amounts: Dict[str, float]
    particle_size: float
    chem_sys: str
    simulation_size: int
    num_realizations: int = 3
    exclude_phases: List[str] = field(default_factory=list)

    def __post_init__(self):
        self.chem_sys = format_chem_sys(self.chem_sys)


    def to_file(self, fname: str):
        with open(fname, 'w+') as f:
            f.write(self.to_json())

    @classmethod
    def from_file(self, fname: str):
        with open(fname, 'rb') as f:
            return self.from_dict(json.loads(f.read()))