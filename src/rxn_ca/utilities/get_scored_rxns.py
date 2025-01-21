from rxn_network.reactions.reaction_set import ReactionSet

from ..core import HeatingSchedule
from ..phases import SolidPhaseSet

from ..reactions import ReactionLibrary, ScoredReaction, ScoredReactionSet, score_rxns
from ..reactions.scorers import BasicScore, TammanScore

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
                    parallel=True):

    lib = ReactionLibrary(phases=phase_set)

    if heating_sched is not None:
        temps = heating_sched.all_temps

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

            results = pool.map(fn, temps)
            for t, r in zip(temps, results):
                lib.add_rxns_at_temp(r, t)
    else:
        if rxns_at_temps is None:
            rxns_at_temps = rxn_set.compute_at_temperatures(temps)
        
        for t in temps:
            scorer = scorer_class(temp=t, phase_set=phase_set, **scorer_kwargs)
            rset = rxns_at_temps.get(t)

            scored_rxns: List[ScoredReaction] = score_rxns(rset, scorer, phase_set=phase_set)
            scored_rset = ScoredReactionSet(scored_rxns, lib.phases)
            lib.add_rxns_at_temp(scored_rset, t)

    return lib