import pytest

from rxn_ca.reactions import ScoredReactionSet


def test_from_file(get_test_file_path):
    sr_set = ScoredReactionSet.from_file(get_test_file_path("reactions/scored_rxn_set.json"))
    assert isinstance(sr_set, ScoredReactionSet)

def test_limit_phases(get_test_file_path):
    sr_set = ScoredReactionSet.from_file(get_test_file_path("reactions/scored_rxn_set.json"))
    before_len = len(sr_set)
    limited = sr_set.limit_phases(["BaO", "TiO2"])
    after_len = len(limited)
    assert after_len < before_len