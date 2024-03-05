from __future__ import annotations

from pydantic import Field

from typing import Optional

from rxn_network.reactions.reaction_set import ReactionSet
from rxn_network.entries.entry_set import GibbsEntrySet
from rxn_ca.phases import SolidPhaseSet


from .base_schema import BaseSchema
from dataclasses import dataclass

@dataclass
class BaseReactionInputs(BaseSchema):

    rxn_set: ReactionSet = Field(description="The enumerated reactions")
    solid_phase_set: SolidPhaseSet = Field(description="The phase set for use in the simulation")
    entry_set: GibbsEntrySet = Field(description="The origianl set of entries used to produce these inputs")
    chem_sys: str = Field(description="The chemical system containing these reactions")
    stability_cutoff: float = Field(description="The energy tolerance for considering a phase stable")
    open_el: Optional[str] = Field(description="An open element")
    chem_pot: Optional[float] = Field(description="The chemical potential of the open element")

    @classmethod
    def from_obj(cls,
                 rxn_set: ReactionSet,
                 solid_phase_set: SolidPhaseSet,
                 entry_set: GibbsEntrySet,
                 chem_sys: str,
                 stability_cutoff: float,
                 open_el: str,
                 chem_pot: float):
        return cls(
            rxn_set = rxn_set,
            solid_phase_set = solid_phase_set,
            entry_set = entry_set,
            chem_sys = chem_sys,
            stability_cutoff = stability_cutoff,
            open_el = open_el,
            chem_pot = chem_pot,
        )

    @classmethod
    def from_dict(cls: BaseReactionInputs, d):

        return cls(
            rxn_set = ReactionSet.from_dict(d['rxn_set']),
            solid_phase_set = SolidPhaseSet.from_dict(d['solid_phase_set']),
            entry_set = GibbsEntrySet.from_dict(d['entry_set']),
            chem_sys = d['chem_sys'],
            stability_cutoff = d['stability_cutoff'],
            open_el = d['open_el'],
            chem_pot = d['chem_pot'],
        )

    def as_dict(self):
        d = super().as_dict()
        return { **d, **{
            "rxn_set": self.rxn_set.as_dict(),
            "solid_phase_set": self.solid_phase_set.as_dict(),
            "entry_set": self.entry_set.as_dict(),
            "chem_sys": self.chem_sys,
            "stability_cutoff": self.stability_cutoff,
            "open_el": self.open_el,
            "chem_pot": self.chem_pot,
        }}