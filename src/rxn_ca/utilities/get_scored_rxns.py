from rxn_network.reactions.reaction_set import ReactionSet

from ..core import HeatingSchedule
from ..phases import SolidPhaseSet

from ..reactions import ReactionLibrary, ScoredReaction, ScoredReactionSet, score_rxns
from ..reactions.scorers import BasicScore, TammanHuttigScoreErf
from tqdm import tqdm

from typing import List

def get_scored_rxns(rxn_set: ReactionSet,
                    heating_sched: HeatingSchedule,
                    scorer_class: BasicScore = TammanHuttigScoreErf,
                    exclude_pure_elements: bool = False,
                    exclude_theoretical: bool = True,
                    exclude_phases: List[str]= [],
                    phase_set: SolidPhaseSet = None):

    lib = ReactionLibrary(phases=phase_set)
    sched_temps = heating_sched.all_temps
    
    for t in tqdm(sched_temps, desc="Getting reaction energies at temperatures..."):
        scorer = scorer_class(temp=t, phase_set=phase_set)
        rset = rxn_set.set_new_temperature(t)

        scored_rxns: List[ScoredReaction] = score_rxns(rset, scorer, phase_set=phase_set)
        scored_rset = ScoredReactionSet(scored_rxns, lib.phases, identity_score=3)
        if exclude_pure_elements:
            print("Excluding reactions involving pure elements...")
            scored_rset = scored_rset.exclude_pure_els()
        
        if exclude_theoretical:
            print("Excluding reactions involving theoretical compounds...")
            scored_rset = scored_rset.exclude_theoretical()

        if len(exclude_phases) > 0:
            print(f"Excluding reactions including {exclude_phases}")
            scored_rset = scored_rset.exclude_phases(exclude_phases)

        lib.add_rxns_at_temp(scored_rset, t)

    return lib