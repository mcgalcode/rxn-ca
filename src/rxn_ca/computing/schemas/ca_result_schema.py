from typing import Any, List, Dict

from ...core.recipe import ReactionRecipe
from ...core.reaction_result import ReactionResult
from ...reactions.reaction_library import ReactionLibrary

from pylattica.core import PeriodicStructure

from .job_types import JobTypes
from monty.json import MSONable
from dataclasses import dataclass

import json

@dataclass
class RxnCAResultDoc(MSONable):

    recipe: ReactionRecipe
    results: List[ReactionResult]
    reaction_library: ReactionLibrary
    job_type: str = JobTypes.RUN_RXN_AUTOMATON.value

    @classmethod
    def from_file(cls, fname):
        with open(fname, "rb") as f:
            return cls.from_dict(json.load(f))
    
    @classmethod
    def from_dict(cls, d):

        return cls(
            recipe = ReactionRecipe.from_dict(d['recipe']),
            results = [ReactionResult.from_dict(d) for d in  d["results"]],
            reaction_library = ReactionLibrary.from_dict(d['reaction_library']),
        )

    def as_dict(self):
        d = super().as_dict()
        return { **d, **{
            "recipe": self.recipe.as_dict(),
            "results": [r.as_dict() for r in self.results],
            "reaction_library": self.reaction_library.as_dict(),
        }}

    def to_file(self, fname):
        with open(fname, "w+") as f:
            json.dump(self.as_dict(), f)