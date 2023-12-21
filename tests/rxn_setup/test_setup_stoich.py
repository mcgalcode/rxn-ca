import pytest

from rxn_ca import get_scored_rxns, setup_reaction
from rxn_ca.core.heating import HeatingSchedule, HeatingStep
from rxn_ca.analysis import ReactionStepAnalyzer

from pymatgen.core.composition import Composition

import numpy as np

def test_stoich():
    expected_mbdown = {
        'BaO': 42.57417234996495,
        'TiO2': 42.364520461311265,
        'Ba2TiO4': 0.18899195946011765,
        'Ba4Ti13O30': 0.01127376098794665,
        'BaTiO3': 0.2825715048280896,
        'BaTi2O5': 0.14207506768818176
    }

    expected_elemental_bdown = {}

    for phase, p_amt in expected_mbdown.items():
        comp = Composition(phase)
        for el, el_amt in comp.get_el_amt_dict().items():
            if el in expected_elemental_bdown:
                expected_elemental_bdown[el] += el_amt * p_amt
            else:
                expected_elemental_bdown[el] = el_amt * p_amt

    sched = HeatingSchedule(HeatingStep.hold(1200, 10000))

    rxn_lib = get_scored_rxns(
        "Ba-Ti-O",
        heating_sched=sched,
    )

    sim = setup_reaction(rxn_lib.phases, expected_mbdown, 30)

    analyzer = ReactionStepAnalyzer(rxn_lib.phases)
    true_el_bdown = analyzer.elemental_composition(sim.state)
    true_phase_bdown = analyzer.molar_breakdown(sim.state)

    EL_TOL = 0.02
    for phase, amt in expected_elemental_bdown.items():
        assert (true_el_bdown.get(phase) - amt) / true_el_bdown.get(phase) < EL_TOL

    PHASE_ABS_TOL = 0.5
    for phase, amt in expected_mbdown.items():
        assert (true_phase_bdown.get(phase) - amt) / true_phase_bdown.get(phase) < PHASE_ABS_TOL

def test_stoich_two():
    expected_molar_bdown = {'YMn2O5': 0.4064785252678766,
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
    
    expected_elemental_bdown = {}

    for phase, p_amt in expected_molar_bdown.items():
        comp = Composition(phase)
        for el, el_amt in comp.get_el_amt_dict().items():
            if el in expected_elemental_bdown:
                expected_elemental_bdown[el] += el_amt * p_amt
            else:
                expected_elemental_bdown[el] = el_amt * p_amt
    
    total_el_mols = sum(expected_elemental_bdown.values())
    expected_elemental_bdown = { el: amt / total_el_mols for el, amt in expected_elemental_bdown.items() }

    sched = HeatingSchedule(HeatingStep.hold(1100, 10000))

    rxn_lib = get_scored_rxns(
        "Y-Mn-Cl-Li-O",
        heating_sched=sched,
    )

    sim = setup_reaction(
        rxn_lib.phases,
        15,
        phase_mole_amts=expected_molar_bdown,
    )

    analyzer = ReactionStepAnalyzer(rxn_lib.phases)
    true_el_bdown = analyzer.elemental_composition_fractional(sim.state)
    true_phase_bdown = analyzer.molar_breakdown(sim.state)

    print(expected_elemental_bdown, true_el_bdown)

    EL_TOL = 0.02
    for phase, amt in expected_elemental_bdown.items():
        dev = np.abs(true_el_bdown.get(phase) - amt) / true_el_bdown.get(phase)
        assert dev < EL_TOL

    PHASE_ABS_TOL = 0.015
    print(expected_molar_bdown)
    print(true_phase_bdown)
    for phase, amt in expected_molar_bdown.items():
        deviation = np.abs(true_phase_bdown.get(phase) - amt) / true_phase_bdown.get(phase)
        assert deviation < PHASE_ABS_TOL