from pylattica.core import StateObserver
from pylattica.core.simulation_state import SimulationState

from .reaction_step_analyzer import phase_volumes_for_state


class PhaseVolumeObserver(StateObserver):
    """Reduces a SimulationState to per-phase absolute volumes (incl. gases).

    This records the minimal sufficient statistic for ReactionStepAnalyzer:
    every derived quantity (moles, mass, fractions, elemental composition) is
    computed from per-phase volumes. Capturing only this dict per frame is
    dramatically lighter than retaining full simulation states, and is what the
    analysis-only run mode stores instead of a full SimulationResult.

    Note: this observer does not capture per-step reaction choices; analyses
    that need them (e.g. assemble_rxn_choices) are unavailable in this mode.
    """

    def observe(self, state: SimulationState, step_no: int) -> dict:
        return phase_volumes_for_state(state)
