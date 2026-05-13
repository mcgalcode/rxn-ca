from rxn_network.reactions.reaction_set import ReactionSet

from ..core import HeatingSchedule
from ..phases import SolidPhaseSet

from ..reactions import ReactionLibrary, ScoredReaction, ScoredReactionSet, score_rxns
from ..reactions.scorers.core import BasicScore, TammanScore

from typing import List

import multiprocessing as mp

_scoring_globals = {}

def fn(temp):
    score_class = _scoring_globals.get('score_class')
    phase_set = _scoring_globals.get('phase_set')
    rxn_set  = _scoring_globals.get('base_rxns')

    scorer = score_class(temp=temp, phase_set=phase_set)
    if _scoring_globals.get("rxns_at_tmps") is None:
        rset = rxn_set.set_new_temperature(temp)
    else:
        rset = _scoring_globals.get("rxns_at_tmps").get(temp)

    scored_rxns: List[ScoredReaction] = score_rxns(rset, scorer, phase_set=phase_set)
    scored_rset = ScoredReactionSet(scored_rxns, phase_set)
    return scored_rset

def get_scored_rxns(rxn_set: ReactionSet,
                    heating_sched: HeatingSchedule = None,
                    temps: List = None,
                    scorer_class: BasicScore = TammanScore,
                    phase_set: SolidPhaseSet = None,
                    rxns_at_temps = None,
                    scorer_kwargs: dict = {},
                    parallel=True,
                    existing_lib: ReactionLibrary = None):
    """Score reactions at specified temperatures.

    Args:
        rxn_set: Base reaction set from rxn_network
        heating_sched: Heating schedule to extract temperatures from
        temps: List of temperatures (alternative to heating_sched)
        scorer_class: Scorer class to use (default: TammanScore)
        phase_set: Phase set with volume and melting point data
        rxns_at_temps: Pre-computed temperature-specific reactions (optional)
        scorer_kwargs: Additional kwargs for scorer
        parallel: Whether to score in parallel (default: True)
        existing_lib: Existing ReactionLibrary to reuse scored temps from.
            Only temperatures not in existing_lib will be scored.

    Returns:
        ReactionLibrary with scored reactions at all temperatures
    """
    # Start with existing library or create new one
    if existing_lib is not None:
        lib = existing_lib
    else:
        lib = ReactionLibrary(phases=phase_set)

    if heating_sched is not None:
        temps = heating_sched.all_temps

    # Filter to only temps we need to score
    temps_to_score = lib.get_missing_temps(temps)

    if len(temps_to_score) == 0:
        # All temps already scored, nothing to do
        return lib

    if rxns_at_temps is not None:
        rxns_at_temps = {int(t): r for t, r in rxns_at_temps.items() }

    if parallel:
        global _scoring_globals

        _scoring_globals['score_class'] = scorer_class
        _scoring_globals['phase_set'] = phase_set
        _scoring_globals['base_rxns'] = rxn_set

        if rxns_at_temps is not None:
            _scoring_globals['rxns_at_tmps'] = rxns_at_temps

        with mp.get_context('fork').Pool(mp.cpu_count()) as pool:
            results = pool.map(fn, temps_to_score)
            for t, r in zip(temps_to_score, results):
                lib.add_rxns_at_temp(r, t)
    else:
        if rxns_at_temps is None:
            rxns_at_temps = rxn_set.compute_at_temperatures(temps_to_score)

        for t in temps_to_score:
            scorer = scorer_class(temp=t, phase_set=phase_set, **scorer_kwargs)
            rset = rxns_at_temps.get(t)

            scored_rxns: List[ScoredReaction] = score_rxns(rset, scorer, phase_set=phase_set)
            scored_rset = ScoredReactionSet(scored_rxns, lib.phases)
            lib.add_rxns_at_temp(scored_rset, t)

    return lib