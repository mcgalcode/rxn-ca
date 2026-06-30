import pytest

from rxn_ca.core.reaction_result import ReactionResult
from rxn_ca.reactions import ReactionLibrary
from rxn_ca.computing.schemas.ca_result_schema import (
    RxnCAResultDoc,
    compress_doc,
    assemble_rxn_choices,
)

from pylattica.core import ObservedResult, SimulationState


@pytest.fixture
def library(get_test_file_path):
    return ReactionLibrary.from_file(get_test_file_path("integration/batio3_library.json"))


@pytest.fixture
def library_path(get_test_file_path):
    return get_test_file_path("integration/batio3_library.json")


def test_no_path_embeds_library(library):
    doc = RxnCAResultDoc(recipe=None, reaction_library=library)
    d = doc.as_dict()
    # Backwards-compatible default: the library is embedded when no path is set.
    assert d["reaction_library"] is not None
    assert d["reaction_library_path"] is None


def test_path_omits_library_from_serialization(library, library_path):
    doc = RxnCAResultDoc(
        recipe=None, reaction_library=library, reaction_library_path=library_path
    )
    d = doc.as_dict()
    # The heavy library is not embedded; only the path is persisted.
    assert d["reaction_library"] is None
    assert d["reaction_library_path"] == library_path
    # In-memory access still returns the object we passed (no disk reload).
    assert doc.reaction_library is library


def test_lazy_load_from_path(library_path):
    doc = RxnCAResultDoc(recipe=None, reaction_library_path=library_path)
    # Nothing loaded until accessed.
    assert doc._reaction_library is None
    loaded = doc.reaction_library
    assert isinstance(loaded, ReactionLibrary)
    # Cached after first access.
    assert doc._reaction_library is loaded


def test_path_ref_roundtrip_is_tiny(library, library_path, tmp_path):
    doc = RxnCAResultDoc(
        recipe=None, reaction_library=library, reaction_library_path=library_path
    )
    fname = str(tmp_path / "pathref.json")
    doc.to_file(fname)

    reloaded = RxnCAResultDoc.from_file(fname)
    assert reloaded.reaction_library_path == library_path
    assert reloaded._reaction_library is None
    # Library resolves lazily from the recorded path.
    assert isinstance(reloaded.reaction_library, ReactionLibrary)


def test_compress_doc_preserves_path(library, library_path):
    doc = RxnCAResultDoc(
        recipe=None, reaction_library=library, reaction_library_path=library_path
    )
    compressed = compress_doc(doc)
    assert compressed.reaction_library_path == library_path
    # The backing library is carried without forcing a lazy reload.
    assert compressed._reaction_library is library


def test_compress_doc_passthrough_for_analysis_only():
    obs = ObservedResult(observer=None, compress_freq=1)
    doc = RxnCAResultDoc(recipe=None, observed_results=[obs])

    # Analysis-only docs hold already-reduced observations, so compression is a
    # no-op and the doc is returned unchanged.
    result = compress_doc(doc)
    assert result is doc


def test_assemble_rxn_choices_warns_and_returns_empty_in_analysis_only():
    result = ReactionResult(SimulationState(), retain_history=False)

    with pytest.warns(UserWarning, match="analysis-only"):
        choices = assemble_rxn_choices(result)

    assert choices == []
