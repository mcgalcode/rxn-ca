from typing import Any, List, Dict

from ...core.recipe import ReactionRecipe

from .job_types import JobTypes
from monty.json import MSONable
from dataclasses import dataclass

import json

@dataclass
class RxnCAResultDoc(MSONable):

    recipe: ReactionRecipe
    results: List[Dict]
    job_type: str = JobTypes.RUN_RXN_AUTOMATON.value

    @classmethod
    def from_file(cls, fname):
        with open(fname, "rb") as f:
            return cls.from_dict(json.loads(f.read()))
    
    @classmethod
    def from_dict(cls, d):
        return cls(
            recipe = ReactionRecipe.from_dict(d['recipe']),
            results = d["results"],
        )

    def as_dict(self):
        d = super().as_dict()
        return { **d, **{
            "recipe": self.recipe.as_dict(),
            "results": self.results,
        }}

    def to_file(self, fname):
        with open(fname, "w+") as f:
            f.write(json.dumps(self.as_dict()))