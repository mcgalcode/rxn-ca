"""Output schemas for rxn-ca workflow jobs."""

from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional


@dataclass
class ReactionLibraryData:
    """Output from the setup_reaction_library job.

    Contains the serialized phase set, plus metadata. A library saved to a file
    is referenced by reaction_library_path and reaction_library_dict is empty;
    otherwise the library is embedded in reaction_library_dict.
    """

    phase_set_dict: Dict[str, Any]
    reaction_library_dict: Dict[str, Any]
    chemical_system: str
    temperatures: List[float]
    phases_available: List[str]

    # If the library was saved to a file, this is the path
    reaction_library_path: Optional[str] = None

    # How the library was built: n_entries, n_reactions, the settings used and
    # the Materials Project database version
    metadata: Optional[Dict[str, Any]] = None

    def as_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "ReactionLibraryData":
        d = {k: v for k, v in d.items() if not k.startswith("@")}
        return cls(**d)


@dataclass
class SimulationOutput:
    """Output from an rxn-ca simulation job.

    Contains analyzed results for quick access, plus file references
    for retrieving full simulation data when needed.

    The phase_set enables conversions between molar/mass/volume amounts:
        phase_set = SolidPhaseSet.from_dict(output.phase_set_dict)
        volume = moles * phase_set.get_vol(phase)
        mass = moles * Composition(phase).weight
    """

    # Analyzed results for quick access
    final_molar_amounts: Dict[str, float]
    molar_amounts_trajectory: Dict[str, List[float]]
    temperature_trajectory: List[float]
    step_indices: List[int]

    # Phase set for unit conversions (relatively small, safe to serialize)
    phase_set_dict: Dict[str, Any]

    # File references for full data (avoids serializing large objects)
    reaction_library_path: Optional[str] = None
    result_doc_path: Optional[str] = None

    # Recipe and metadata
    recipe_dict: Dict[str, Any] = field(default_factory=dict)
    chemical_system: Optional[str] = None
    num_realizations: int = 1
    metadata: Optional[Dict[str, Any]] = None

    def as_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "SimulationOutput":
        d = {k: v for k, v in d.items() if not k.startswith("@")}
        return cls(**d)
