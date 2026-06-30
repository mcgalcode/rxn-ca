import warnings
from typing import List, Optional

from ...core.recipe import ReactionRecipe
from ...core.constants import REACTION_CHOSEN
from ...core.reaction_result import ReactionResult
from ...reactions.reaction_library import ReactionLibrary
from ...phases.solid_phase_set import SolidPhaseSet
from pylattica.core.simulation_result import compress_result
from pylattica.core import ObservedResult, Simulation

from .base_schema import BaseSchema
from dataclasses import dataclass, field, InitVar

@dataclass
class RxnCAResultDoc(BaseSchema):

    recipe: ReactionRecipe
    results: List[ReactionResult] = field(default_factory=list)
    final_simulation: Simulation = None
    phases: SolidPhaseSet = None
    metadata: dict = None
    # Populated in analysis-only mode instead of `results`: a lightweight,
    # per-frame reduced view (per-phase volumes) from a PhaseVolumeObserver.
    observed_results: List[ObservedResult] = field(default_factory=list)
    # If set, the reaction library lives in this file rather than being embedded
    # in the doc. The (large) library is then resolved lazily on first access
    # and is omitted from as_dict() so it is never duplicated to disk.
    reaction_library_path: Optional[str] = None
    # Passed at construction but stored on a private backing attribute so the
    # public `reaction_library` property can lazy-load from reaction_library_path.
    reaction_library: InitVar[Optional[ReactionLibrary]] = None

    def __post_init__(self, reaction_library):
        # When no value is supplied the dataclass passes the property object
        # itself (the class-level default); normalize that to None.
        if isinstance(reaction_library, property):
            reaction_library = None
        self._reaction_library = reaction_library

    @property
    def reaction_library(self) -> Optional[ReactionLibrary]:
        if self._reaction_library is None and self.reaction_library_path is not None:
            self._reaction_library = ReactionLibrary.from_file(self.reaction_library_path)
        return self._reaction_library

    @reaction_library.setter
    def reaction_library(self, value):
        self._reaction_library = value

    def as_dict(self) -> dict:
        d = super().as_dict()
        # reaction_library is an InitVar (not a dataclass field), so add it here.
        # Read the backing attribute directly so serialization never triggers a
        # lazy load. When a path is available, persist only the path and skip the
        # heavy library entirely.
        if self.reaction_library_path is not None:
            d["reaction_library"] = None
        elif self._reaction_library is not None:
            d["reaction_library"] = self._reaction_library.as_dict()
        else:
            d["reaction_library"] = None
        return d

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
                          reaction_library=result_doc._reaction_library,
                          reaction_library_path=result_doc.reaction_library_path,
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