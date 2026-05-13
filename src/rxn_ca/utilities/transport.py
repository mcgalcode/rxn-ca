from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np
import pkg_resources
from monty.json import MontyDecoder
from pymatgen.core.composition import Composition
from pydantic import BaseModel, Field

KB = 8.6173303e-5  # Boltzmann constant in eV/K

L_MIN_VALID = 1e8
L_MAX_VALID = 1e24


def _canonical_formula(formula: str) -> str:
    return Composition(formula).reduced_formula


def _chemsys_from_filename(stem: str) -> str:
    """Extract chemsys string from a JSON filename stem, e.g. 'Ba-Ti-O-transport' -> 'Ba-Ti-O'."""
    for suffix in ("-transport", "-self"):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return stem


def _mean_tensors(
    tensors: list[np.ndarray],
    min_valid: float = L_MIN_VALID,
    max_valid: float = L_MAX_VALID,
) -> np.ndarray:
    """
    Elementwise mean of tensors with out-of-range masking.
    Components outside [min_valid, max_valid] are excluded; if all are excluded
    for a given element, fall back to the plain mean.
    """
    if not tensors:
        raise ValueError("No tensors provided for aggregation.")
    arr = np.stack([np.array(t, dtype=float) for t in tensors], axis=0)
    mask = (arr >= float(min_valid)) & (arr <= float(max_valid))
    masked = np.where(mask, arr, np.nan)
    mean = np.nanmean(masked, axis=0)
    fallback = arr.mean(axis=0)
    return np.where(np.isnan(mean), fallback, mean)


def _aggregate_docs(docs: list[TransportDoc]) -> TransportDoc:
    """Average multiple TransportDocs that share the same (formula, temperature)."""
    if not docs:
        raise ValueError("No docs to aggregate.")
    rf = docs[0].reduced_formula
    T = float(docs[0].temperature)
    if any(d.reduced_formula != rf for d in docs):
        raise ValueError("Cannot aggregate docs with different formulas.")
    if any(float(d.temperature) != T for d in docs):
        raise ValueError("Cannot aggregate docs with different temperatures.")
    if any(d.species != docs[0].species for d in docs):
        raise ValueError("Cannot aggregate docs with different species ordering.")
    if any(d.mapping != docs[0].mapping for d in docs):
        raise ValueError("Cannot aggregate docs with different mapping.")

    L = _mean_tensors([d.L_tensor for d in docs])

    self_tensors = [d.L_tensor_self for d in docs if d.L_tensor_self is not None]
    L_self = _mean_tensors(self_tensors) if self_tensors else None
    L_dis = (L - L_self) if L_self is not None else None

    volumes = [d.volume for d in docs if d.volume is not None]
    volume = float(np.mean(volumes)) if volumes else None

    return TransportDoc(
        reduced_formula=rf,
        temperature=T,
        species=list(docs[0].species),
        mapping=dict(docs[0].mapping),
        L_tensor=L,
        L_tensor_self=L_self,
        L_tensor_dis=L_dis,
        volume=volume,
    )


def _doc_from_dict(d: dict) -> TransportDoc:
    return TransportDoc(
        reduced_formula=d["reduced_formula"],
        temperature=float(d["temperature"]),
        species=list(d["species"]),
        mapping={str(k): int(v) for k, v in d["mapping"].items()},
        L_tensor=np.array(d["L_tensor"], dtype=float),
        L_tensor_self=np.array(d["L_tensor_self"], dtype=float)
        if d.get("L_tensor_self") is not None
        else None,
        L_tensor_dis=np.array(d["L_tensor_dis"], dtype=float)
        if d.get("L_tensor_dis") is not None
        else None,
        volume=float(d["volume"]) if d.get("volume") is not None else None,
    )


def _doc_to_dict(doc: TransportDoc) -> dict:
    return {
        "reduced_formula": doc.reduced_formula,
        "temperature": float(doc.temperature),
        "species": list(doc.species),
        "mapping": dict(doc.mapping),
        "L_tensor": doc.L_tensor.tolist(),
        "L_tensor_self": doc.L_tensor_self.tolist()
        if doc.L_tensor_self is not None
        else None,
        "L_tensor_dis": doc.L_tensor_dis.tolist()
        if doc.L_tensor_dis is not None
        else None,
        "volume": float(doc.volume) if doc.volume is not None else None,
    }


class TransportDoc(BaseModel):
    model_config = {"arbitrary_types_allowed": True}

    reduced_formula: str = Field(..., description="Reduced formula of the phase.")
    temperature: float = Field(..., description="Temperature in K.")
    species: List[str] = Field(..., description="Species ordering for tensor indices.")
    mapping: Dict[str, int] = Field(..., description="Species label -> tensor index.")
    L_tensor: np.ndarray = Field(..., description="Onsager tensor (3x3) in 1/cm/s/eV.")
    L_tensor_self: Optional[np.ndarray] = Field(
        None, description="Self Onsager tensor (3x3) in 1/cm/s/eV."
    )
    L_tensor_dis: Optional[np.ndarray] = Field(
        None, description="Distinct Onsager tensor (3x3) in 1/cm/s/eV."
    )
    volume: Optional[float] = Field(None, description="Cell volume (Å^3 if provided).")

    def get_transport_coefficient(
        self, species: List[str] | str, self_transport: bool = False
    ) -> float:
        if isinstance(species, list) and len(species) == 2:
            i = self.mapping[species[0]]
            j = self.mapping[species[1]]
            return float(self.L_tensor[i, j])
        if isinstance(species, str):
            i = self.mapping[species]
            if self_transport:
                if self.L_tensor_self is None:
                    raise ValueError("L_tensor_self is not available for this doc.")
                return float(self.L_tensor_self[i, i])
            return float(self.L_tensor[i, i])
        raise ValueError("species must be a single str or a list of two str.")

    @classmethod
    def from_py_oats_dict(cls, d: dict) -> "TransportDoc":
        """
        Accept a serialized TransportDoc from py-oats (Monty/pydantic dict),
        and only keep fields needed by rxn_ca scoring.
        """
        decoded = MontyDecoder().process_decoded(d)
        reduced_formula = getattr(decoded, "reduced_formula", None) or decoded.get("reduced_formula")
        if not reduced_formula:
            full_formula = getattr(decoded, "composition", None) or decoded.get("composition")
            reduced_formula = _canonical_formula(str(full_formula))
        temperature = getattr(decoded, "temperature", None) or decoded.get("temperature")
        species = getattr(decoded, "species", None) or decoded.get("species")
        mapping = getattr(decoded, "mapping", None) or decoded.get("mapping")
        L_tensor = getattr(decoded, "L_tensor", None) or decoded.get("L_tensor")
        L_tensor_self = getattr(decoded, "L_tensor_self", None) or decoded.get("L_tensor_self")
        volume = getattr(decoded, "volume", None) or decoded.get("volume")

        if any(v is None for v in [reduced_formula, temperature, species, mapping, L_tensor]):
            raise ValueError("py-oats transport doc is missing required fields.")

        L_tensor = np.array(L_tensor, dtype=float)
        L_tensor_self_arr = (
            np.array(L_tensor_self, dtype=float) if L_tensor_self is not None else None
        )
        L_tensor_dis = (
            L_tensor - L_tensor_self_arr if L_tensor_self_arr is not None else None
        )

        return cls(
            reduced_formula=str(reduced_formula),
            temperature=float(temperature),
            species=[str(s) for s in species],
            mapping={str(k): int(v) for k, v in mapping.items()},
            L_tensor=L_tensor,
            L_tensor_self=L_tensor_self_arr,
            L_tensor_dis=L_tensor_dis,
            volume=float(volume) if volume is not None else None,
        )

    def fluxes(self, mu: np.ndarray, self_transport: bool = False) -> np.ndarray:
        """Compute flux vector J = L @ mu."""
        if self_transport:
            if self.L_tensor_self is None:
                raise ValueError("L_tensor_self is not available for this doc.")
            return np.dot(self.L_tensor_self, mu)
        return np.dot(self.L_tensor, mu)


class TransportDatabase:
    """
    In-memory database of TransportDoc objects for a chemical system.

    Supports ingestion of docs from JSON files, nearest-temperature lookup,
    aggregation of duplicate (formula, T) entries, validation, and flux computation.
    """

    def __init__(self) -> None:
        self._docs: list[TransportDoc] = []
        self._available_formulas: set[str] = set()

    @property
    def docs(self) -> list[TransportDoc]:
        return list(self._docs)

    @property
    def available_formulas(self) -> set[str]:
        return set(self._available_formulas)

    def ingest(self, docs: TransportDoc | Iterable[TransportDoc]) -> None:
        """Add one or more TransportDoc objects to the database."""
        if isinstance(docs, TransportDoc):
            self._docs.append(docs)
            self._available_formulas.add(docs.reduced_formula)
        else:
            docs = list(docs)
            self._docs.extend(docs)
            self._available_formulas.update(d.reduced_formula for d in docs)

    def get_doc(self, formula: str, temperature: float) -> TransportDoc:
        """
        Return the TransportDoc for (formula, temperature).
        Uses nearest-temperature lookup if an exact match is not found.
        Aggregates by averaging if multiple docs match the same (formula, T).
        """
        rf = _canonical_formula(formula)
        t = float(temperature)
        candidates = [d for d in self._docs if d.reduced_formula == rf]
        if not candidates:
            raise KeyError(f"No transport data for formula={rf!r}")
        exact = [d for d in candidates if float(d.temperature) == t]
        if exact:
            return _aggregate_docs(exact) if len(exact) > 1 else exact[0]
        closest_T = min(
            set(float(d.temperature) for d in candidates), key=lambda x: abs(x - t)
        )
        closest = [d for d in candidates if float(d.temperature) == closest_T]
        return _aggregate_docs(closest) if len(closest) > 1 else closest[0]

    def validate(self) -> list[str]:
        """
        Check diagonal L_tensor components against physically expected bounds.
        Returns a list of warning strings; does not raise.
        """
        issues = []
        for doc in self._docs:
            for sp in doc.species:
                i = doc.mapping[sp]
                diag = float(doc.L_tensor[i, i])
                if not (L_MIN_VALID <= abs(diag) <= L_MAX_VALID):
                    issues.append(
                        f"{doc.reduced_formula} @ {doc.temperature}K: "
                        f"L[{sp},{sp}]={diag:.2e} outside "
                        f"[{L_MIN_VALID:.0e}, {L_MAX_VALID:.0e}]"
                    )
        return issues

    def compute_fluxes(
        self,
        formula: str,
        temperature: float,
        mu: dict[str, float],
        self_transport: bool = False,
    ) -> dict[str, float]:
        """
        Compute species fluxes J = L @ mu_vec for the given (formula, temperature).

        mu is a dict mapping element symbol -> chemical potential distance (eV).
        Returns a dict mapping element symbol -> flux value.

        Falls back to identity-scaled flux if formula is not in the database.
        """
        try:
            doc = self.get_doc(formula, temperature)
        except KeyError:
            scale = KB * float(temperature) * 1e14
            return {sp: scale * mu.get(sp, 0.0) for sp in mu}

        mu_vec = np.array([mu.get(sp, 0.0) for sp in doc.species], dtype=float)
        flux_vec = doc.fluxes(mu_vec, self_transport=self_transport)
        return {sp: float(flux_vec[i]) for i, sp in enumerate(doc.species)}

    @classmethod
    def from_json(cls, path: str | Path) -> "TransportDatabase":
        """
        Load a TransportDatabase from a JSON file.
        Accepts both a bare list of doc dicts and the legacy
        {'chemsys': ..., 'docs': [...]} format.
        """
        data = json.loads(Path(path).read_text())
        db = cls()
        docs_raw = data if isinstance(data, list) else data.get("docs", [])
        for doc_d in docs_raw:
            db.ingest(_doc_from_dict(doc_d))
        return db

    @classmethod
    def from_chemsys(cls, chemsys: str | list[str]) -> "TransportDatabase":
        """
        Build a TransportDatabase by loading all JSON files in the bundled
        transport_db/ directory whose chemsys element set matches the requested one.

        Both *-transport.json and *-self.json files are loaded and merged,
        so the returned database contains all available transport data.
        Validation warnings are emitted via the warnings module.
        """
        if isinstance(chemsys, list):
            target_els = frozenset(chemsys)
        else:
            target_els = frozenset(chemsys.split("-"))

        transport_dir = Path(
            pkg_resources.resource_filename("rxn_ca.phases", "transport_db")
        )
        db = cls()
        for json_path in sorted(transport_dir.glob("*.json")):
            file_chemsys = _chemsys_from_filename(json_path.stem)
            if frozenset(file_chemsys.split("-")) <= target_els:
                db.ingest(cls.from_json(json_path).docs)

        for issue in db.validate():
            warnings.warn(issue, stacklevel=2)

        return db

    def to_json(self, path: str | Path) -> None:
        """Serialize all docs to a JSON file as a list of doc dicts."""
        Path(path).write_text(
            json.dumps([_doc_to_dict(d) for d in self._docs], indent=2, sort_keys=True)
        )
