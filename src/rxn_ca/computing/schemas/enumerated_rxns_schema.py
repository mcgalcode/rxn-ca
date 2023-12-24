from __future__ import annotations

from pydantic import Field

from typing import Optional

from rxn_network.reactions.reaction_set import ReactionSet

from ..utils.functions import format_chem_sys
from ...phases import SolidPhaseSet

from .base_schema import BaseSchema
from dataclasses import dataclass

@dataclass
class EnumeratedRxnsModel(BaseSchema):

    rxn_set: dict = Field(description="The enumerated reactions")
    phases: SolidPhaseSet = Field(description="The phases present in the enumerated reactions")
    chem_sys: str = Field(description="The chemical system containing these reactions")
    stability_cutoff: float = Field(description="The energy tolerance for considering a phase stable")
    open_el: Optional[str] = Field(description="An open element")
    chem_pot: Optional[float] = Field(description="The chemical potential of the open element")

    @classmethod
    def from_obj(cls,
                 rxn_set: ReactionSet,
                 chem_sys: str,
                 stability_cutoff: float,
                 open_el: str,
                 chem_pot: float):
        return cls(
            rxn_set = rxn_set,
            chem_sys = chem_sys,
            stability_cutoff = stability_cutoff,
            open_el = open_el,
            chem_pot = chem_pot,
            phases = SolidPhaseSet.from_rxn_set(rxn_set)
        )

    @classmethod
    def from_dict(cls: EnumeratedRxnsModel, d):

        return cls(
            rxn_set = ReactionSet.from_dict(d['rxn_set']),
            phases = SolidPhaseSet.from_dict(d['phases']),
            chem_sys = d['chem_sys'],
            stability_cutoff = d['stability_cutoff'],
            open_el = d['open_el'],
            chem_pot = d['chem_pot'],
        )

    def as_dict(self):
        d = super().as_dict()
        return { **d, **{
            "rxn_set": self.rxn_set.as_dict(),
            "chem_sys": self.chem_sys,
            "stability_cutoff": self.stability_cutoff,
            "open_el": self.open_el,
            "chem_pot": self.chem_pot,
            "phases": self.phases.as_dict()
        }}