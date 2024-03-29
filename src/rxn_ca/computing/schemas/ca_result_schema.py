from typing import List

from ...core.recipe import ReactionRecipe
from ...core.constants import REACTION_CHOSEN
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
    metadata: dict = None

def compress_doc(result_doc: RxnCAResultDoc, num_steps=100):
    results = result_doc.results
    compressed = [compress_result(r, num_steps) for r in results]
    return RxnCAResultDoc(recipe=result_doc.recipe,
                          results=compressed,
                          phases=result_doc.phases,
                          reaction_library=result_doc.reaction_library,
                          metadata=result_doc.metadata)

def get_metadata_from_results(results: List[ReactionResult]):
    return {
        "rxn_choices": [assemble_rxn_choices(r) for r in results]
    }

def assemble_rxn_choices(result: ReactionResult):
    choices = []
    for step in result.steps():
        rxn_choice = step.get_general_state(REACTION_CHOSEN)
        if rxn_choice != None:
            choices.append(rxn_choice)
    return choices