from typing import List

from ...core.recipe import ReactionRecipe
from ...core.reaction_result import ReactionResult
from ...reactions.reaction_library import ReactionLibrary
from ...phases.solid_phase_set import SolidPhaseSet
from pylattica.core.simulation_result import compress_result

from .base_schema import BaseSchema
from dataclasses import dataclass

@dataclass
class RxnCAResultDoc(BaseSchema):

    recipe: ReactionRecipe
    results: List[ReactionResult]
    reaction_library: ReactionLibrary = None
    phases: SolidPhaseSet = None
    
    # @classmethod
    # def from_dict(cls, d):
    #     if d.get('reaction_library') is not None:
    #         rlib = ReactionLibrary.from_dict(d['reaction_library'])
    #     else:
    #         rlib = None

    #     return cls(
    #         recipe = ReactionRecipe.from_dict(d['recipe']),
    #         results = [ReactionResult.from_dict(d) for d in  d["results"]],
    #         reaction_library = rlib,
    #     )

    # def as_dict(self):
    #     d = super().as_dict()
    #     if self.reaction_library is not None:
    #         rlib = self.reaction_library.as_dict()
    #     else:
    #         rlib = None
    #     return { **d, **{
    #         "recipe": self.recipe.as_dict(),
    #         "results": [r.as_dict() for r in self.results],
    #         "reaction_library": rlib,
    #     }}

def compress_doc(result_doc: RxnCAResultDoc, num_steps=100):
    results = result_doc.results
    compressed = [compress_result(r, num_steps) for r in results]
    return RxnCAResultDoc(recipe=result_doc.recipe, results=compressed, phases=result_doc.phases, reaction_library=result_doc.reaction_library)