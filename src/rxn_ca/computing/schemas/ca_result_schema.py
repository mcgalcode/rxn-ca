import warnings
from typing import List

from ...core.recipe import ReactionRecipe
from ...core.constants import REACTION_CHOSEN
from ...core.reaction_result import ReactionResult
from ...reactions.reaction_library import ReactionLibrary
from ...phases.solid_phase_set import SolidPhaseSet
from pylattica.core.simulation_result import compress_result
from pylattica.core import ObservedResult, Simulation

from .base_schema import BaseSchema
from dataclasses import dataclass, field

@dataclass
class RxnCAResultDoc(BaseSchema):

    recipe: ReactionRecipe
    results: List[ReactionResult] = field(default_factory=list)
    reaction_library: ReactionLibrary = None
    final_simulation: Simulation = None
    phases: SolidPhaseSet = None
    metadata: dict = None
    # Populated in analysis-only mode instead of `results`: a lightweight,
    # per-frame reduced view (per-phase volumes) from a PhaseVolumeObserver.
    observed_results: List[ObservedResult] = field(default_factory=list)

def compress_doc(result_doc: RxnCAResultDoc, num_steps=100):
    # Analysis-only docs hold already-reduced ObservedResults (no diffs to
    # resample), so compression does not apply -- return the doc unchanged.
    if result_doc.observed_results:
        return result_doc

    results = result_doc.results
    compressed = [compress_result(r, num_steps) for r in results]
    return RxnCAResultDoc(recipe=result_doc.recipe,
                          results=compressed,
                          phases=result_doc.phases,
                          reaction_library=result_doc.reaction_library,
                          final_simulation=result_doc.final_simulation,
                          metadata=result_doc.metadata)

def get_metadata_from_results(results: List[ReactionResult]):
    return {
        "rxn_choices": [assemble_rxn_choices(r) for r in results]
    }

def assemble_rxn_choices(result: ReactionResult):
    # Reaction choices live on the full per-step states. In analysis-only mode
    # those states are not retained (the PhaseVolumeObserver does not capture
    # them), so there is nothing to assemble.
    if not result.retain_history:
        warnings.warn(
            "assemble_rxn_choices is unavailable in analysis-only mode "
            "(retain_history=False); returning an empty list of choices."
        )
        return []

    choices = []
    for step in result.steps():
        rxn_choice = step.get_general_state(REACTION_CHOSEN)
        choices.append(rxn_choice)
    return choices