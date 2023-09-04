from typing import Any, List, Tuple, Dict
from pydantic import BaseModel, Field

from ...core.heating import HeatingSchedule

from ..utils.functions import format_chem_sys
from .job_types import JobTypes

import json

class RxnCAResultDoc(BaseModel):

    chem_sys: str = Field(description="The chemical system used")
    heating_schedule: HeatingSchedule = Field(description="The temperature of the simulation")
    reactant_amts: Dict[str, float] = Field(description="Initial reactant amounts")
    results: List[Dict] = Field(description="The serialized result object")
    job_type: str = Field(default=JobTypes.RUN_RXN_AUTOMATON.value)

    @classmethod
    def from_file(cls, fname):
        with open(fname, "rb") as f:
            return cls.from_dict(json.loads(f.read()))
    
    @classmethod
    def from_dict(cls, d):
        return cls(
            chem_sys = d["chem_sys"],
            heating_schedule = HeatingSchedule.from_dict(d["heating_schedule"]),
            reactant_amts = d["reactant_amts"],
            results = d["results"],
        )

    def __init__(self, **data: Any):
        super().__init__(**data)
        self.chem_sys = format_chem_sys(self.chem_sys)

    def as_dict(self):
        return {
            "chem_sys": self.chem_sys,
            "heating_schedule": self.heating_schedule.as_dict(),
            "reactant_amts": self.reactant_amts,
            "results": self.results,
        }

    def to_file(self, fname):
        with open(fname, "w+") as f:
            f.write(json.dumps(self.as_dict()))