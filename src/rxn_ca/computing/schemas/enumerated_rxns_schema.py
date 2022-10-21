from uuid import uuid4
from pydantic import BaseModel, Field

from typing import Optional, Any

from rxn_network.reactions.reaction_set import ReactionSet

from ..utils.functions import format_chem_sys

from .job_types import JobTypes

class EnumeratedRxnsModel(BaseModel):

    task_id: str = Field(description="The ID of this result")
    rxn_set: dict = Field(description="The enumerated reactions")
    chem_sys: str = Field(description="The chemical system containing these reactions")
    stability_cutoff: float = Field(description="The energy tolerance for considering a phase stable")
    open_el: Optional[str] = Field(description="An open element")
    chem_pot: Optional[float] = Field(description="The chemical potential of the open element")
    temperature: int = Field(description="The temperature at which the energy of the reactions is calculated")
    job_type: str = Field(default=JobTypes.ENUMERATE_RXNS.value)

    def __init__(self, **data: Any):
        super().__init__(**data)
        self.chem_sys = format_chem_sys(self.chem_sys)

    @classmethod
    def from_obj(cls,
                 rxn_set: ReactionSet,
                 chem_sys: str,
                 stability_cutoff: float,
                 open_el: str,
                 chem_pot: float,
                 temperature: int):
        return cls(
            task_id = str(uuid4()),
            rxn_set = rxn_set.as_dict(),
            chem_sys = chem_sys,
            stability_cutoff = stability_cutoff,
            open_el = open_el,
            chem_pot = chem_pot,
            temperature = temperature
        )