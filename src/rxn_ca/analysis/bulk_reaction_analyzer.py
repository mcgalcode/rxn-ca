from __future__ import annotations

from ..phases.solid_phase_set import SolidPhaseSet
from ..core.reaction_result import ReactionResult
from ..core.heating import HeatingSchedule
from ..phases.gasses import DEFAULT_GASES
from .reaction_step_analyzer import ReactionStepAnalyzer

from ..computing.schemas.ca_result_schema import RxnCAResultDoc

from pylattica.core import ObservedResult

from typing import Tuple, List, Union
import numpy as np
    
def color(color, text):
    return f"<span style='color:{str(color)}'> {str(text)} </span>"

OTHERS = "others"
class BulkReactionAnalyzer():
    """A class that stores the result of running a simulation. Keeps track of all
    the steps that the simulation proceeded through, and the set of reactions that
    was used in the simulation.
    """

    @classmethod
    def from_result_doc_file(cls, fname: str) -> BulkReactionAnalyzer:
        doc: RxnCAResultDoc = RxnCAResultDoc.from_file(fname)
        return cls.from_result_doc(doc)

    @classmethod
    def from_result_doc(cls, doc: RxnCAResultDoc) -> BulkReactionAnalyzer:
        # In analysis-only mode the doc carries lightweight ObservedResults
        # instead of full ReactionResults.
        results = doc.observed_results if doc.observed_results else doc.results
        return cls(results, doc.phases, doc.recipe.heating_schedule)
    
    def __init__(self, results: Union[List[ReactionResult], List[ObservedResult]], phase_set: SolidPhaseSet, heating_sched: HeatingSchedule):
        """Initializes a ReactionResult with the reaction set used in the simulation

        Args:
            results: Either full ReactionResults, or, in analysis-only mode, the
                ObservedResults produced by a PhaseVolumeObserver. In analysis-only
                mode only volume-derived quantities are available (no per-site or
                reaction-choice analyses), and steps can only be analyzed at the
                observed cadence.
            phase_set (SolidPhaseSet):
            heating_sched (HeatingSchedule):
        """
        self.step_analyzer = ReactionStepAnalyzer(phase_set)
        self.heating_schedule = heating_sched
        self.phase_set = phase_set

        # Analysis-only mode: results carry precomputed per-phase volumes
        # (ObservedResult) rather than reconstructable full states.
        self._analysis_only = isinstance(results[0], ObservedResult)

        if self._analysis_only:
            # ObservedResult.__len__ counts observations, not steps; use the
            # tracked total step count to keep result_length consistent.
            self.result_length = results[0]._total_steps + 1
        else:
            self.result_length = len(results[0])
        self.results = results
        self._results_loaded = False
        self._step_idxs = None
        self._step_groups = None


    @property
    def last_loaded_step_idx(self):
        self._get_step_groups()
        return self._step_idxs[-1]
    
    @property
    def loaded_step_idxs(self):
        self._get_step_groups()
        return self._step_idxs
    
    @property
    def loaded_step_groups(self):
        self._get_step_groups()
        return self._step_groups
      
    def get_analyzer(self, step_group):
        # In analysis-only mode a "step group" is precomputed per-phase volumes
        # (one dict per result); otherwise it is full simulation states.
        if self._analysis_only:
            return self.step_analyzer.set_precomputed_volumes(step_group)
        return self.step_analyzer.set_step_group(step_group)

    def analyze_step(self, step_number):
        return self.get_analyzer(self.get_steps(step_number))

    def get_step_size(self):
        if self._analysis_only:
            raise RuntimeError(
                "get_step_size requires full simulation states; it is "
                "unavailable in analysis-only mode."
            )
        return self.get_analyzer(self.get_first_steps()[0]).get_simulation_size()

    def get_elemental_amounts_at(self, step_no):
        return self.analyze_step(step_no).get_molar_elemental_composition()

    def molar_fractional_breakdown(self, step_no):
        return self.analyze_step(step_no).get_all_absolute_molar_amounts()
    
    def mass_fraction_breakdown(self, step_no):
        return self.analyze_step(step_no).get_all_mass_fractions()
    
    def get_final_molar_breakdown(self):
        return self.get_analyzer(self.get_final_steps()).get_all_mole_fractions()

    def get_all_absolute_molar_amounts(self, step_no: int):
        return self.analyze_step(step_no).get_all_absolute_molar_amounts()
    
    def get_steps(self, step_no):
        if self._analysis_only:
            # Only observed steps are available; surface a clear error otherwise.
            return [r.observation_at(step_no) for r in self.results]
        return [r.get_step(step_no) for r in self.results]

    def get_final_steps(self):
        if self._analysis_only:
            return [r.observation_at(r.observation_steps[-1]) for r in self.results]
        return [r.last_step for r in self.results]

    def get_first_steps(self):
        if self._analysis_only:
            return [r.observation_at(r.observation_steps[0]) for r in self.results]
        return [r.first_step for r in self.results]

    def all_phases_present(self, min_mass_fraction_prevalence=0.0):
        analyses = [self.get_analyzer(sg).get_all_mass_fractions() for sg in self.loaded_step_groups]
        all_phases = set()

        for analysis in analyses:
            for phase,amt in analysis.items():
                if amt > min_mass_fraction_prevalence:
                    all_phases.add(phase)
        return all_phases


    def num_phases_present(self, min_mass_fraction_prevalence=0.0):
        return len(self.all_phases_present(min_mass_fraction_prevalence))
    
    def phases_with_prevalence(self, min_prevalence, max_prevalence=1.0):
        analyses = [self.get_analyzer(sg).get_all_mass_fractions() for sg in self.loaded_step_groups]
        
        phases_over_min = set()
        phases_with_excess = set()

        for a in analyses:
            for p, v in a.items():
                if v > min_prevalence:
                    phases_over_min.add(p)
                if v > max_prevalence:
                    phases_with_excess.add(p)
        return list(phases_over_min - phases_with_excess)

    def get_volume_trace(self):
        return [self.get_analyzer(sg).get_total_volume() for sg in self.loaded_step_groups]

    def get_condensed_mass_trace(self):
        return [self.get_analyzer(sg).get_total_mass() for sg in self.loaded_step_groups]

    def _get_step_groups(self) -> Tuple[List[int], List]:
        if self._step_idxs is None:
            first_result = self.results[0]

            if self._analysis_only:
                # Observations are already reduced and stored at their own
                # cadence; group them across results by observed step number.
                self._step_idxs = first_result.observation_steps
                self._results_loaded = True
                self._step_groups = [
                    [r.observation_at(step_idx) for r in self.results]
                    for step_idx in self._step_idxs
                ]
                return self._step_idxs, self._step_groups

            # Check if results use live_compress (frames already stored)
            if hasattr(first_result, '_frames') and first_result._frames:
                # Use the stored frames directly - no reconstruction needed
                self._step_idxs = sorted(first_result._frames.keys())
                self._results_loaded = True
            else:
                # Traditional mode: load steps at computed interval.
                # NOTE: The previous formula (num_points = result_length / 2;
                # step_size = round(result_length / num_points)) always evaluates to 2,
                # so this has always sampled every other step. Intent is unverified, so we
                # preserve that exact behavior here rather than risk changing all analysis
                # outputs. Revisit if a fixed target point count is actually desired.
                step_size = 2
                if not self._results_loaded:
                    [r.load_steps(step_size) for r in self.results]
                    self._results_loaded = True
                self._step_idxs = list(range(0, self.result_length, step_size))

            self._step_groups = [[r.get_step(step_idx) for r in self.results] for step_idx in self._step_idxs]

        return self._step_idxs, self._step_groups
    
    def has_simulation_converged(self, convergence_criteria : float = 0.001) -> bool:
        """Check if a simulation has converged based on the phase data"""

        converged = True
        if self._analysis_only:
            # Only observed steps exist; use the final window of observations.
            final_idxs = self.loaded_step_idxs[-10:]
            phase_amounts = [self.get_all_absolute_molar_amounts(i) for i in final_idxs]
        else:
            phase_amounts = [self.get_all_absolute_molar_amounts(i) for i in range(self.result_length-10, self.result_length)]
        for phase in phase_amounts[0].keys():
            if phase in [*DEFAULT_GASES, SolidPhaseSet.FREE_SPACE]:
                continue
            if np.std([p[phase] for p in phase_amounts]) / np.mean([p[phase] for p in phase_amounts]) > convergence_criteria:
                converged = False
                break
        return converged