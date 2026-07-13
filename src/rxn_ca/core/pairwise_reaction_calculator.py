import random
from copy import deepcopy
from typing import Dict, List, Tuple

from pylattica.discrete.state_constants import DISCRETE_OCCUPANCY
from pylattica.core.constants import GENERAL, SITE_ID, SITES
from pylattica.core.simulation_state import SimulationState

from ..phases.solid_phase_set import SolidPhaseSet
from ..reactions import ScoredReaction
from .constants import VOLUME, GASES_EVOLVED, REACTION_CHOSEN, RXN_PRODUCT_LEDGER
from .reaction_calculator import ReactionCalculator, choose_from_list

# Relative tolerance below which a cell's leftover volume is treated as fully
# drained (floating point residue from the limiting-reagent arithmetic).
_DRAIN_EPS = 1e-9


class PairwiseReactionCalculator(ReactionCalculator):
    """Stoichiometric pair-consumption update scheme.

    The legacy ReactionCalculator gives each cell in an interaction an
    independent should_reaction_proceed draw, so for A + B -> products the
    realized consumed A:B ratio random-walks away from the reaction
    stoichiometry, while the product credit assumes it is exact. This subclass
    replaces the per-site draws with a single draw per interaction event:

    1. All cells in the interaction are pooled by phase, and a
       limiting-reagent extent is computed so that every proceeding event
       consumes the participating phases in *exact* stoichiometric ratio and
       fully drains at least one cell (the limiting pool).
    2. Partially consumed cells keep their phase with a decremented volume.
    3. Each freed cell is replaced by the solid product currently furthest
       below its stoichiometric share, tracked by a per-reaction production
       ledger (largest-deficit / largest-remainder selection). This bounds
       the product-mix error at roughly one cell's volume per product per
       reaction, instead of a random walk.

    The ledger is kept on the calculator instance rather than in the
    simulation state: routing it through the general state costs a deepcopy
    of the (growing) ledger per event plus a full copy in every stored step
    diff, roughly tripling simulation time. Set store_ledger=True to mirror
    it into the general state under RXN_PRODUCT_LEDGER for debugging or
    analysis - the mirrored copies are snapshots, at debug-only cost.

    Interaction enumeration, interaction/reaction selection and scoring are
    inherited unchanged, so results are directly comparable with the legacy
    scheme.
    """

    def __init__(self,
        neighborhood_graph,
        scored_rxns = None,
        inertia = 0.05,
        atmospheric_species = [],
        store_ledger = False,
    ) -> None:
        super().__init__(
            neighborhood_graph,
            scored_rxns=scored_rxns,
            inertia=inertia,
            atmospheric_species=atmospheric_species,
        )
        self.store_ledger = store_ledger
        self._ledger: Dict = {}

    def get_state_update(self, site_id: int, prev_state: SimulationState):
        updates = {}

        possible_interactions = self.possible_interactions_at_site(site_id, prev_state)
        selected_interaction = self.choose_interaction(possible_interactions)

        if selected_interaction.is_no_op:
            return updates

        updates[GENERAL] = {}
        updates[SITES] = {}

        # Reaction selection is identical to the legacy scheme
        rxns: List[ScoredReaction] = selected_interaction.reactions
        selected_reaction: ScoredReaction = choose_from_list(rxns, [rxn.competitiveness for rxn in rxns])
        selected_reaction_id: int = self.rxn_set.get_rxn_id(selected_reaction)
        updates[GENERAL][REACTION_CHOSEN] = selected_reaction_id

        # Pool the interaction's cells by phase. Two cells may share a phase
        # (e.g. a same-phase neighbor pair matching a single-reactant
        # reaction), in which case their volumes are pooled.
        cells_by_phase: Dict[str, List[Dict]] = {}
        for site_state in selected_interaction.site_states:
            phase = site_state[DISCRETE_OCCUPANCY]
            cells_by_phase.setdefault(phase, []).append(site_state)

        # Limiting-reagent extent: the largest total solid reactant volume
        # consumable while keeping every phase within its available pool.
        # Consuming delta_v * r_i / R from each phase i keeps the consumed
        # ratio exactly stoichiometric, and the limiting pool is fully
        # drained, so every proceeding event frees at least one cell.
        total_r = selected_reaction.total_solid_reactant_stoich
        delta_v = min(
            sum(cell[VOLUME] for cell in cells) * total_r / selected_reaction.reactant_stoich(phase)
            for phase, cells in cells_by_phase.items()
        )

        # One proceed draw per event. Kinetics match the legacy scheme: the
        # legacy per-site draw succeeds with p_i = spf * (r_i / R) / vol_i and
        # consumes vol_i, so E[consumed_i per event] = spf * r_i / R. Here a
        # single draw with q = spf / delta_v consumes delta_v * r_i / R from
        # phase i, giving the same expectation q * delta_v * r_i / R
        # = spf * r_i / R. The legacy 1/vol factor is subsumed by 1/delta_v.
        # Reactions with no solid products are left unscaled by spf, matching
        # should_reaction_proceed.
        if selected_reaction.total_solid_product_stoich > 0:
            proceed_prob = selected_reaction.solid_product_fraction / delta_v
        else:
            proceed_prob = 1.0 / delta_v

        if random.random() >= min(proceed_prob, 1.0):
            return updates

        # Gas evolution is deterministic and atomic with consumption: the
        # event credits the stoichiometric share of each gaseous product for
        # the full consumed solid volume delta_v.
        if len(selected_reaction.gas_products) > 0:
            # NOTE: get_general_state returns a deepcopy, so this dict is
            # safe to mutate
            gas_amts = prev_state.get_general_state().get(GASES_EVOLVED, {})

            for gas_phase in selected_reaction.gas_products:
                evolved = selected_reaction.gas_amt_from_reactant_vol(gas_phase, delta_v)
                gas_amts[gas_phase] = gas_amts.get(gas_phase, 0) + evolved

            updates[GENERAL][GASES_EVOLVED] = gas_amts

        # Consume delta_v * r_i / R from each phase pool, draining cells in
        # order. Fully drained cells are freed (candidates for product
        # placement); a partially drained cell keeps its phase with a
        # decremented volume.
        freed_cells: List[Tuple[int, float]] = []  # (site_id, consumed volume)
        for phase, cells in cells_by_phase.items():
            to_consume = delta_v * selected_reaction.reactant_stoich(phase) / total_r
            for cell in cells:
                cell_vol = cell[VOLUME]
                consumed = min(cell_vol, to_consume)
                to_consume -= consumed
                remaining = cell_vol - consumed

                if remaining <= _DRAIN_EPS * cell_vol:
                    # Any float residue is folded into the freed cell so the
                    # totals stay exact.
                    freed_cells.append((cell[SITE_ID], cell_vol))
                else:
                    updates[SITES][cell[SITE_ID]] = {
                        DISCRETE_OCCUPANCY: phase,
                        VOLUME: remaining,
                    }

        total_s = selected_reaction.total_solid_product_stoich

        if total_s == 0:
            # Full volatilization: every freed cell becomes empty space
            for freed_site_id, _ in freed_cells:
                updates[SITES][freed_site_id] = {
                    DISCRETE_OCCUPANCY: SolidPhaseSet.FREE_SPACE,
                    VOLUME: 1.0,
                }
            return updates

        # Distribute the total solid product volume across the freed cells in
        # proportion to each cell's consumed volume. Product identities are
        # chosen by largest deficit against the calculator-local per-reaction
        # production ledger, updated sequentially within the event.
        v_out_total = delta_v * total_s / total_r
        freed_total = sum(consumed for _, consumed in freed_cells)

        for freed_site_id, consumed in freed_cells:
            v_out = v_out_total * consumed / freed_total
            product_phase = self.choose_product_by_deficit(selected_reaction, v_out, self._ledger)
            updates[SITES][freed_site_id] = {
                DISCRETE_OCCUPANCY: product_phase,
                VOLUME: v_out,
            }

        if self.store_ledger:
            # Snapshot (not the live dict): set_general_state stores the
            # reference, so writing self._ledger directly would let later
            # events mutate previously stored step diffs in place.
            updates[GENERAL][RXN_PRODUCT_LEDGER] = deepcopy(self._ledger)

        return updates

    def choose_product_by_deficit(self, rxn: ScoredReaction, v_out: float, ledger: Dict) -> str:
        """Chooses the solid product currently furthest below its
        stoichiometric share of this reaction's cumulative production, then
        records v_out against it in the ledger (mutated in place).

        The deficit of product j is d_j = (P + v_out) * s_j / S - P_j, where
        P_j is the cumulative volume produced as j so far and P their sum.
        Largest-remainder selection keeps every |P_j - P * s_j / S| bounded
        by about one event volume, so the realized product mix tracks the
        reaction stoichiometry instead of random-walking.
        """
        rxn_ledger = ledger.setdefault(rxn._as_str, {})
        produced_total = sum(rxn_ledger.values())
        total_s = rxn.total_solid_product_stoich

        best_product = None
        best_deficit = None
        # Sorted for deterministic tie-breaking
        for product in sorted(rxn.solid_products):
            target = (produced_total + v_out) * rxn.product_stoich(product) / total_s
            deficit = target - rxn_ledger.get(product, 0)
            if best_deficit is None or deficit > best_deficit:
                best_deficit = deficit
                best_product = product

        rxn_ledger[best_product] = rxn_ledger.get(best_product, 0) + v_out
        return best_product
