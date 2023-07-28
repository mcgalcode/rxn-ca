from typing import Any, List, Tuple, Dict
from pydantic import BaseModel, Field

from ..utils.functions import format_chem_sys
from .job_types import JobTypes

class RxnCAResultDoc(BaseModel):

    chem_sys: str = Field(description="The chemical system used")
    heating_schedule: List[Tuple[int, int]] = Field(description="The temperature of the simulation")
    reactant_amts: Dict[str, float] = Field(description="Initial reactant amounts")
    results: List[Dict] = Field(description="The serialized result object")
    job_type: str = Field(default=JobTypes.RUN_RXN_AUTOMATON.value)

    def __init__(self, **data: Any):
        super().__init__(**data)
        self.chem_sys = format_chem_sys(self.chem_sys)