from pydantic import BaseModel, Field
from typing import Any
from uuid import uuid4

from .job_types import JobTypes
from ..utils.functions import format_chem_sys
from ...reactions import ScoredReactionSet

class ScoredRxnsModel(BaseModel):

    task_id: str = Field(description="The ID of this result")
    scored_rxn_set: dict = Field(description="The scored reactions")
    chem_sys: str = Field(description="The chemical system containing these reactions")
    job_type: str = Field(default=JobTypes.SCORE_RXNS.value)
    temperature: int = Field(description="The temperature used to score these reactions")

    def __init__(self, **data: Any):
        super().__init__(**data)
        self.chem_sys = format_chem_sys(self.chem_sys)

    @classmethod
    def from_obj(cls, scored_rxn_set: ScoredReactionSet, chem_sys: str, temp: int):
        return cls(
            task_id = str(uuid4()),
            scored_rxn_set = scored_rxn_set.as_dict(),
            chem_sys = chem_sys,
            temperature = temp
        )