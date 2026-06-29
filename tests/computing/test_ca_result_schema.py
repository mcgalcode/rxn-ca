import pytest

from rxn_ca.core.reaction_result import ReactionResult
from rxn_ca.computing.schemas.ca_result_schema import (
    RxnCAResultDoc,
    compress_doc,
    assemble_rxn_choices,
)

from pylattica.core import ObservedResult, SimulationState


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
