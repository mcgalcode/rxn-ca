from __future__ import annotations

from typing import Optional

from tqdm import tqdm

from rxn_network.reactions.reaction_set import ReactionSet

from ...phases.solid_phase_set import SolidPhaseSet
from ..scored_reaction import ScoredReaction
from .core import (
    KB,
    Na,
    BasicScore,
    ConstantScore,
    GibbsErfScore,
    TammanHuttigScoreErf,
    TammanHuttigScoreExponential,
    TammanHuttigScoreSoftplus,
    TammanScore,
    TammanTightLinear,
    erf,
    huttig_score_exp,
    huttig_score_softplus,
    softplus,
    tamman_score_exp,
    tamman_score_softplus,
)
from .diffusion import DiffusionScorer, mu_distance


def score_rxns(
    reactions: ReactionSet, scorer: BasicScore, phase_set: Optional[SolidPhaseSet] = None
):
    scored_reactions = []

    for rxn in tqdm(
        reactions.get_rxns(), desc=f"Scoring reactions... at temp {scorer.temp}"
    ):
        reactants = [r.reduced_formula for r in rxn.reactants]
        non_gases = [r for r in reactants if r not in phase_set.gas_phases]
        if len(non_gases) > 0:
            scored_rxn = ScoredReaction.from_rxn_network(
                scorer.score(rxn), rxn, phase_set.volumes
            )
            scored_reactions.append(scored_rxn)

    return scored_reactions


__all__ = [
    "KB",
    "Na",
    "BasicScore",
    "ConstantScore",
    "DiffusionScorer",
    "GibbsErfScore",
    "TammanHuttigScoreErf",
    "TammanHuttigScoreExponential",
    "TammanHuttigScoreSoftplus",
    "TammanScore",
    "TammanTightLinear",
    "erf",
    "huttig_score_exp",
    "huttig_score_softplus",
    "mu_distance",
    "score_rxns",
    "softplus",
    "tamman_score_exp",
    "tamman_score_softplus",
]
