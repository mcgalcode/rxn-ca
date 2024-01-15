from pathlib import Path
import pytest
import json

from rxn_ca.phases import SolidPhaseSet

@pytest.fixture()
def get_test_file_path():
    def _(file_path: str):
        return str(Path(__file__).parent / file_path)

    return _

@pytest.fixture
def ymno3_phases(get_test_file_path):
    fpath = get_test_file_path("core/ymno3_phases.json")
    with open(fpath, 'r+') as f:
        d = json.load(f)
        return SolidPhaseSet.from_dict(d)