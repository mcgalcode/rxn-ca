from typing import List

from ...core.recipe import ReactionRecipe
from ...core.reaction_result import ReactionResult
from ...reactions.reaction_library import ReactionLibrary

from .base_schema import BaseSchema

from dataclasses import dataclass

@dataclass
class RxnCAResultDoc(BaseSchema):

    recipe: ReactionRecipe
    results: List[ReactionResult]
    reaction_library: ReactionLibrary
    
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
