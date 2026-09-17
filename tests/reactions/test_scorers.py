import pytest
from pymatgen.core.composition import Composition

from rxn_ca.phases import SolidPhaseSet
from rxn_ca.reactions.scorers import (
    TammanHuttigStrict,
    TammanStrict,
    erf_tight,
    huttig_score_softplus,
    tamman_score_softplus,
)

TEMP = 1200


class _Rxn:
    """Minimal stand-in exposing the two attributes scorers read from a
    ComputedReaction."""

    def __init__(self, reactant_formulas, energy_per_atom):
        self.reactants = [Composition(f) for f in reactant_formulas]
        self.energy_per_atom = energy_per_atom


@pytest.fixture
def phases():
    return SolidPhaseSet(
        ["BaCO3", "BaO", "TiO2", "BaTiO3", "CO2"],
        volumes={"BaCO3": 50.0, "BaO": 25.0, "TiO2": 30.0, "BaTiO3": 55.0, "CO2": 25.0},
        densities={"BaCO3": 4.29, "BaO": 5.72, "TiO2": 4.23, "BaTiO3": 6.02, "CO2": 1.98},
        melting_points={"BaCO3": 1600, "BaO": 2196, "TiO2": 2116, "BaTiO3": 1898, "CO2": 216},
        experimentally_observed={p: True for p in ["BaCO3", "BaO", "TiO2", "BaTiO3", "CO2"]},
    )


def test_erf_tight_gate_shape():
    # Positive dG is strictly closed, unlike the default erf gate (~0.07 at 0)
    assert erf_tight(0.0) < 1e-6
    assert erf_tight(0.05) < 1e-12

    # Midpoint at -0.10 eV/atom, saturating around -0.16
    assert erf_tight(-0.10) == pytest.approx(0.5)
    assert erf_tight(-0.16) > 0.99

    # Gradual on the negative side: weakly exothermic reactions stay slow
    assert erf_tight(-0.05) < 0.01


def test_tamman_strict(phases):
    scorer = TammanStrict(phases, temp=TEMP)

    # Downhill: strict gate on dG times the Tammann factor at T / min(Tm)
    downhill = _Rxn(["BaO", "TiO2"], -0.2)
    expected = tamman_score_softplus(TEMP / 2116) * erf_tight(-0.2)
    assert scorer.score(downhill) == pytest.approx(expected)

    # Uphill reactions are shut off
    uphill = _Rxn(["BaO", "TiO2"], 0.05)
    assert scorer.score(uphill) < 1e-12


def test_tamman_huttig_strict_branches(phases):
    scorer = TammanHuttigStrict(phases, temp=TEMP)

    # Two solid reactants: Tammann factor
    two_solid = _Rxn(["BaO", "TiO2"], -0.2)
    expected = tamman_score_softplus(TEMP / 2116) * erf_tight(-0.2)
    assert scorer.score(two_solid) == pytest.approx(expected)

    # Single solid reactant (decomposition): Huttig factor, as in the
    # default TammanHuttigScoreErf
    one_solid = _Rxn(["BaCO3"], -0.12)
    expected = huttig_score_softplus(TEMP / 1600) * erf_tight(-0.12)
    assert scorer.score(one_solid) == pytest.approx(expected)

    # Uphill reactions are shut off
    assert scorer.score(_Rxn(["BaO", "TiO2"], 0.05)) < 1e-12
