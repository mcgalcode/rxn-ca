import pytest

from rxn_ca.utilities.get_scored_rxns import get_scored_rxns
from rxn_ca.utilities.setup_reaction import setup_reaction
from rxn_ca.core.heating import HeatingSchedule, HeatingStep
from rxn_ca.phases import SolidPhaseSet
from rxn_ca.analysis import ReactionStepAnalyzer

from pymatgen.core.composition import Composition

import numpy as np

def test_stoich():
    expected_molar_ratios = {
        'BaO': 42.57417234996495,
        'TiO2': 42.364520461311265,
        'Ba2TiO4': 0.18899195946011765,
        'Ba4Ti13O30': 0.01127376098794665,
        'BaTiO3': 0.2825715048280896,
        'BaTi2O5': 0.14207506768818176
    }

    phases = SolidPhaseSet.from_phase_list(list(expected_molar_ratios.keys()))

    expected_elemental_bdown = phases.mole_amts_to_el_amts(expected_molar_ratios)

    sim = setup_reaction(phases, precursor_mole_ratios=expected_molar_ratios)

    analyzer = ReactionStepAnalyzer(phases)
    true_el_bdown = analyzer.get_molar_elemental_composition(sim.state)
    true_phase_bdown = analyzer.get_all_absolute_molar_amounts(sim.state)

    EL_TOL = 0.02
    for phase, amt in expected_elemental_bdown.items():
        assert (true_el_bdown.get(phase) - amt) / true_el_bdown.get(phase) < EL_TOL

    PHASE_ABS_TOL = 0.5
    for phase, amt in expected_molar_ratios.items():
        assert (true_phase_bdown.get(phase) - amt) / true_phase_bdown.get(phase) < PHASE_ABS_TOL

def test_stoich_two():
    expected_molar_ratios = {'YMn2O5': 0.4064785252678766,
        'LiCl': 60.94572336147908,
        'Mn2O3': 10.81542281053917,
        'Y2Mn2O7': 0.023635039329169646,
        'LiMn2O4': 0.013504240936119735,
        'Mn3O4': 0.4494577785146918,
        'Li2O': 0.19844145043221514,
        'Y2O3': 1.2534252647363726,
        'YClO': 8.018652052558396,
        'YCl3': 0.2504141211998018,
        'MnO': 0.32232120436269973,
        'LiYO2': 0.03866374211455529,
        'YMnO3': 0.9250325636410563,
        'Li2MnO3': 0.08037032761918489,
        'LiMnO2': 0.09949602764689137
    }

    phases = SolidPhaseSet.from_phase_list(list(expected_molar_ratios.keys()))

    expected_elemental_bdown = phases.mole_amts_to_el_fracs(expected_molar_ratios)
    
    sim = setup_reaction(
        phases,
        precursor_mole_ratios=expected_molar_ratios,
    )

    analyzer = ReactionStepAnalyzer(phases)
    true_el_bdown = analyzer.get_fractional_elemental_composition(sim.state)
    true_phase_bdown = analyzer.get_all_absolute_molar_amounts(sim.state)

    print(expected_elemental_bdown, true_el_bdown)

    EL_TOL = 0.02
    for phase, amt in expected_elemental_bdown.items():
        dev = np.abs(true_el_bdown.get(phase) - amt) / true_el_bdown.get(phase)
        assert dev < EL_TOL

    PHASE_ABS_TOL = 0.015
    print(expected_molar_ratios)
    print(true_phase_bdown)
    for phase, amt in expected_molar_ratios.items():
        deviation = np.abs(true_phase_bdown.get(phase) - amt) / true_phase_bdown.get(phase)
        assert deviation < PHASE_ABS_TOL