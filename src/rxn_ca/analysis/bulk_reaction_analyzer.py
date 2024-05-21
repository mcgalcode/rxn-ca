from __future__ import annotations

from ..phases.solid_phase_set import SolidPhaseSet
from ..core.reaction_result import ReactionResult
from ..core.heating import HeatingSchedule
from .reaction_step_analyzer import ReactionStepAnalyzer

from ..computing.schemas.ca_result_schema import RxnCAResultDoc

from typing import Tuple, List

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
        return cls(doc.results, doc.phases, doc.recipe.heating_schedule)
    
    @classmethod
    def from_result_doc(cls, doc: RxnCAResultDoc) -> BulkReactionAnalyzer:
        return cls(doc.results, doc.phases, doc.recipe.heating_schedule)
    
    def __init__(self, results: List[ReactionResult], phase_set: SolidPhaseSet, heating_sched: HeatingSchedule):
        """Initializes a ReactionResult with the reaction set used in the simulation

        Args:
            rxn_set (ScoredReactionSet):
        """
        self.step_analyzer = ReactionStepAnalyzer(phase_set)
        self.heating_schedule = heating_sched

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
        return self.step_analyzer.set_step_group(step_group)
    
    def analyze_step(self, step_number):
        return self.get_analyzer(self.get_steps(step_number))
    
    def get_step_size(self):
        return self.get_analyzer(self.get_first_steps()[0]).get_simulation_size()

    def get_elemental_amounts_at(self, step_no):
        return self.analyze_step(step_no).get_molar_elemental_composition()

    def molar_fractional_breakdown(self, step_no):
        return self.analyze_step(step_no).get_all_absolute_molar_amounts()
    
    def get_final_molar_breakdown(self):
        return self.get_analyzer(self.get_final_steps()).get_all_mole_fractions()

    def get_all_absolute_molar_amounts(self, step_no: int):
        return self.analyze_step(step_no).get_all_absolute_molar_amounts()
    
    def get_steps(self, step_no):
        return [r.get_step(step_no) for r in self.results]
    
    def get_final_steps(self):
        return [r.last_step for r in self.results]

    def get_first_steps(self):
        return [r.first_step for r in self.results]

    def all_phases_present(self):
        sgs = [self.get_analyzer(sg).phases_present() for sg in self.loaded_step_groups]

        all_phases = set()
        for sg in sgs:
            all_phases.update(sg)
        return list(all_phases)
    
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
            num_points = self.result_length / 2
            step_size = max(1, round(self.result_length / num_points))
            if not self._results_loaded:
                [r.load_steps(step_size) for r in self.results]
                self._results_loaded = True
            self._step_idxs = list(range(0, self.result_length, step_size))
            self._step_groups = [[r.get_step(step_idx) for r in self.results] for step_idx in self._step_idxs]

        return self._step_idxs, self._step_groups
        