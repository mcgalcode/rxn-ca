"""get_scored_rxns: scorer_kwargs reach the scorer, serially and in parallel, and
a parallel call scores its own reactions, not ones left from an earlier call."""

import pytest
from pymatgen.core import Composition
from rxn_network.entries.gibbs import GibbsComputedEntry
from rxn_network.reactions.computed import ComputedReaction
from rxn_network.reactions.reaction_set import ReactionSet

from rxn_ca.phases import SolidPhaseSet
from rxn_ca.reactions.scorers.core import BasicScore, ConstantScore
from rxn_ca.utilities.get_scored_rxns import get_scored_rxns


class _FixedScore(BasicScore):
    def __init__(self, phase_set, temp=None, value=1.0):
        super().__init__(phase_set, temp)
        self.value = value

    def score(self, rxn):
        return self.value


def _rxn_set(metal):
    """One reaction: metal + peroxide -> oxide, e.g. Li + Li2O2 -> Li2O."""
    entries = [GibbsComputedEntry(Composition(f), formation_energy_per_atom=e,
                                  volume_per_atom=10.0, temperature=300, entry_id=f)
               for f, e in [(metal, 0.0), (f"{metal}2O2", -1.6), (f"{metal}2O", -2.07)]]
    return ReactionSet.from_rxns([ComputedReaction.balance(entries[:2], entries[2:])])


def _phases(*metals):
    phases = [f for m in metals for f in (m, f"{m}2O2", f"{m}2O")]
    return SolidPhaseSet(phases, volumes={p: 10.0 for p in phases},
                         densities={}, melting_points={}, experimentally_observed={})


@pytest.mark.parametrize("parallel", [False, True])
def test_scorer_kwargs_reach_the_scorer(parallel):
    lib = get_scored_rxns(_rxn_set("Li"), temps=[300], scorer_class=_FixedScore,
                          phase_set=_phases("Li"), scorer_kwargs={"value": 0.25}, parallel=parallel)

    assert [r.competitiveness for r in lib.get_rxns_at_temp(300).reactions] == [0.25]


def test_parallel_call_does_not_reuse_an_earlier_calls_reactions():
    li, na = _rxn_set("Li"), _rxn_set("Na")
    phases = _phases("Li", "Na")
    get_scored_rxns(li, temps=[300], scorer_class=ConstantScore, phase_set=phases,
                    rxns_at_temps=li.compute_at_temperatures([300]), parallel=True)

    lib = get_scored_rxns(na, temps=[300], scorer_class=ConstantScore, phase_set=phases,
                          parallel=True)

    assert [sorted(r.reactants) for r in lib.get_rxns_at_temp(300).reactions] == [["Na", "Na2O2"]]
