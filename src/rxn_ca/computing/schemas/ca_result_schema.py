from typing import Any, List, Tuple, Dict
from pydantic import BaseModel, Field

from ..utils.functions import format_chem_sys
from .job_types import JobTypes

import json

class RxnCAResultDoc(BaseModel):

    chem_sys: str = Field(description="The chemical system used")
    heating_schedule: List[Tuple[int, int]] = Field(description="The temperature of the simulation")
    reactant_amts: Dict[str, float] = Field(description="Initial reactant amounts")
    results: List[Dict] = Field(description="The serialized result object")
    job_type: str = Field(default=JobTypes.RUN_RXN_AUTOMATON.value)

    @classmethod
    def from_file(cls, fname):
        with open(fname, "rb") as f:
            return cls(**json.loads(f.read()))

    def __init__(self, **data: Any):
        super().__init__(**data)
        self.chem_sys = format_chem_sys(self.chem_sys)

    
    def to_file(self, fname):
        with open(fname, "w+") as f:
            f.write(json.dumps(dict(self)))