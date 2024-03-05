import pytest

from rxn_ca.phases.solid_phase_set import get_melting_points, SolidPhaseSet

import numpy as np

NA_CL = "NaCl"
LI2_O = "Li2O"

@pytest.fixture
def basic_phase_set():
    return SolidPhaseSet(
        [NA_CL, LI2_O],
        volumes={
            NA_CL: 2.0,
            LI2_O: 0.5
        },
        melting_points={
            NA_CL: 800,
            LI2_O: 1000
        },
        experimentally_observed={
            NA_CL: True,
            LI2_O: False
        })

def test_get_melting_points_from_data():
    formula = "BaTiO3"

    mps = get_melting_points([formula])

    assert mps.get(formula) is not None

def test_get_melting_point(basic_phase_set: SolidPhaseSet):
    assert basic_phase_set.get_melting_point(NA_CL) == 800
    assert basic_phase_set.get_melting_point(LI2_O) == 1000

def test_get_vol(basic_phase_set: SolidPhaseSet):
    assert basic_phase_set.get_vol(NA_CL) == 2.0
    assert basic_phase_set.get_vol(LI2_O) == 0.5

def test_is_theoretical(basic_phase_set: SolidPhaseSet):
    assert not basic_phase_set.is_theoretical(NA_CL)
    assert basic_phase_set.is_theoretical(LI2_O)

def test_get_theoretical_phases(basic_phase_set: SolidPhaseSet):
    phases = basic_phase_set.get_theoretical_phases()

    assert len(phases) == 1
    assert phases[0] == LI2_O

def test_volume_conversion(basic_phase_set: SolidPhaseSet):
    nacl_vol = 3.0
    nacl_moles = basic_phase_set.vol_to_moles(nacl_vol, NA_CL)
    assert np.isclose(nacl_moles, 1.5)

    back_to_vol = basic_phase_set.moles_to_vol(nacl_moles, NA_CL)
    assert np.isclose(back_to_vol, nacl_vol)

def test_volume_conversion_multiple(basic_phase_set: SolidPhaseSet):
    vols = {
        NA_CL: 3.0,
        LI2_O: 4.0
    }

    moles = basic_phase_set.vol_amts_to_moles(vols)
    assert np.isclose(moles.get(NA_CL), 1.5)
    assert np.isclose(moles.get(LI2_O), 8.0)

    back_to_vols = basic_phase_set.mole_amts_to_vols(moles)
    assert np.isclose(back_to_vols.get(NA_CL), 3.0)
    assert np.isclose(back_to_vols.get(LI2_O), 4.0)

def test_mole_to_el_conversions(basic_phase_set: SolidPhaseSet):
    moles = {
        NA_CL: 2,
        LI2_O: 3
    }

    els = basic_phase_set.mole_amts_to_el_amts(moles)

    assert np.isclose(els.get("Na"), 2)
    assert np.isclose(els.get("Cl"), 2)
    assert np.isclose(els.get("Li"), 6)
    assert np.isclose(els.get("O"), 3)

def test_vol_to_el_conversions(basic_phase_set: SolidPhaseSet):
    vols = {
        NA_CL: 4,
        LI2_O: 3
    }

    els = basic_phase_set.vol_amts_to_el_amts(vols)

    assert np.isclose(els.get("Na"), 2)
    assert np.isclose(els.get("Cl"), 2)
    assert np.isclose(els.get("Li"), 12)
    assert np.isclose(els.get("O"), 6)

def test_make_phase_set_from_mp():
    phases = [NA_CL, LI2_O]
    pset = SolidPhaseSet.from_phase_list(phases)
    # +1 for the FREE_SPACE phase
    assert len(pset) == len(phases) + 1
