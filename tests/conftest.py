from pathlib import Path
import pytest

@pytest.fixture()
def get_test_file_path():
    def _(file_path: str):
        return str(Path(__file__).parent / file_path)

    return _