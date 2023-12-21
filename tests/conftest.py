import pytest

from pylattica.discrete import PhaseSet
from pylattica.structures.square_grid import DiscreteGridSetup

from pylattica.structures.square_grid import SimpleSquare2DStructureBuilder, DiscreteGridSetup
from pylattica.discrete import PhaseSet

@pytest.fixture(scope="module")
def square_grid_2D_2x2():
    return SimpleSquare2DStructureBuilder().build(2)

@pytest.fixture(scope="module")
def square_grid_2D_4x4():
    return SimpleSquare2DStructureBuilder().build(4)

@pytest.fixture(scope="module")
def simple_phase_set():
    return PhaseSet(["A", "B", "C", "D"])

@pytest.fixture(scope="module")
def grid_setup(simple_phase_set):
    return DiscreteGridSetup(simple_phase_set)