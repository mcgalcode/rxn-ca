from pylattica.core import StateObserver
from pylattica.core.simulation_state import SimulationState

from ..core.constants import TEMPERATURE
from .reaction_step_analyzer import phase_volumes_for_state

# Keys of the structured observation produced by PhaseVolumeObserver.
VOLUMES = "volumes"


class PhaseVolumeObserver(StateObserver):
    """Reduces a SimulationState to per-phase absolute volumes (incl. gases).

    Each observation is a structured dict ``{VOLUMES: {...}, TEMPERATURE: T}``:

    - ``VOLUMES`` is the per-phase absolute volume mapping -- the minimal
      sufficient statistic for ReactionStepAnalyzer, from which every derived
      quantity (moles, mass, fractions, elemental composition) is computed.
    - ``TEMPERATURE`` is the temperature active in the state at observation
      time, captured so analysis-only runs retain the true per-frame
      temperature (rather than recomputing it from the heating schedule).

    Recording this per frame is dramatically lighter than retaining full
    simulation states, and is what the analysis-only run mode stores instead
    of a full SimulationResult.

    Note: this observer does not capture per-step reaction choices; analyses
    that need them (e.g. assemble_rxn_choices) are unavailable in this mode.
    """

    def observe(self, state: SimulationState, step_no: int) -> dict:
        return {
            VOLUMES: phase_volumes_for_state(state),
            TEMPERATURE: state.get_general_state().get(TEMPERATURE),
        }
