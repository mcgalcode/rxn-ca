import matplotlib.pyplot as plt
import plotly.graph_objects as go

from ..phases.solid_phase_set import SolidPhaseSet
from ..core.reaction_result import ReactionResult
from ..reactions.scored_reaction_set import ScoredReactionSet

from .bulk_step_analyzer import BulkReactionStepAnalyzer
from .reaction_step_analyzer import ReactionStepAnalyzer

from typing import Tuple, List


class BulkReactionAnalyzer():
    """A class that stores the result of running a simulation. Keeps track of all
    the steps that the simulation proceeded through, and the set of reactions that
    was used in the simulation.
    """


    def __init__(self, results: List[ReactionResult]):
        """Initializes a ReactionResult with the reaction set used in the simulation

        Args:
            rxn_set (ScoredReactionSet):
        """
        self.rxn_set: ScoredReactionSet = results[0].rxn_set
        self.bulk_step_analyzer = BulkReactionStepAnalyzer(self.rxn_set.phases)
        self.single_step_analyzer = ReactionStepAnalyzer(self.rxn_set.phases)

        self.result_length = len(results[0])
        self.results = results

    def plot_elemental_amounts(self) -> None:
        fig = go.Figure()
        fig.update_layout(width=800, height=800, title="Molar Elemental Amount vs time step")
        print("updated")
        fig.update_layout(yaxis_range=[0,None])

        fig.update_yaxes(title="# of Moles")
        fig.update_xaxes(title="Simulation Step")

        elements = list(self.single_step_analyzer.elemental_composition(self.results[0].first_step).keys())
        traces = []
        
        step_idxs, step_groups = self._get_step_groups()
        amounts = [self.bulk_step_analyzer.elemental_composition(sg) for sg in step_groups]
        for el in elements:
            ys = [a.get(el, 0) for a in amounts]
            traces.append((step_idxs, ys, el))


        for t in traces:
            fig.add_trace(go.Scatter(name=t[2], x=t[0], y=t[1], mode='lines'))

        fig.show()


    def plot_mole_fractions(self, min_prevalence=0.01) -> None:
        """In a Jupyter Notebook environment, plots the phase prevalence traces for the simulation.

        Returns:
            None:
        """

        fig = go.Figure()
        fig.update_layout(width=800, height=800, title="Prevalence by Simulation Step")
        fig.update_yaxes(title="Prevalence", range=(0, None))

        traces = []
        step_idxs, step_groups = self._get_step_groups()
        fig.update_xaxes(range=[0, step_idxs[-1]], title="Simulation Step")
        molar_breakdowns = [self.bulk_step_analyzer.molar_fractional_breakdown(sg) for sg in step_groups]

        phases = set()
        for bd in molar_breakdowns:
            phases = phases.union(set(bd.keys()))
            
        for phase in phases:
            if phase is not SolidPhaseSet.FREE_SPACE:
                ys = [mb.get(phase, 0) for mb in molar_breakdowns]
                traces.append((step_idxs, ys, phase))

        filtered_traces = [t for t in traces if max(t[1]) > min_prevalence]

        for t in filtered_traces:
            fig.add_trace(go.Scatter(name=t[2], x=t[0], y=t[1], mode='lines'))

        fig.show()

    def plot_molar_phase_amounts(self, min_prevalence=0.01) -> None:
        """In a Jupyter Notebook environment, plots the phase prevalence traces for the simulation.

        Returns:
            None:
        """

        fig = go.Figure()
        fig.update_layout(width=800, height=800, title="Absolute Molar Prevalence by Simulation Step")
        fig.update_yaxes(title="# of Moles", range=[0, None])

        traces = []
        step_idxs, step_groups = self._get_step_groups()
        fig.update_xaxes(range=[0, step_idxs[-1]], title="Simulation Step")
        molar_breakdowns = [self.bulk_step_analyzer.molar_breakdown(step_group) for step_group in step_groups]

        phases = set()
        for bd in molar_breakdowns:
            phases = phases.union(set(bd.keys()))
            
        for phase in phases:
            if phase is not SolidPhaseSet.FREE_SPACE:
                ys = [mb.get(phase, 0) for mb in molar_breakdowns]
                traces.append((step_idxs, ys, phase))

        filtered_traces = [t for t in traces if max(t[1]) > min_prevalence]

        for t in filtered_traces:
            fig.add_trace(go.Scatter(name=t[2], x=t[0], y=t[1], mode='lines'))

        fig.show()

    def as_dict(self):
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "steps": [s.as_dict() for s in self.steps],
            "rxn_set": self.rxn_set.as_dict(),
            "phase_set": self.phase_set.as_dict()
        }

    def _get_step_groups(self) -> Tuple[List[int], List]:
        num_points = min(100, self.result_length)
        step_size = max(1, round(self.result_length / num_points))
        [r.load_steps(step_size) for r in self.results]
        step_idxs = list(range(0, self.result_length, step_size))
        return step_idxs, [[r.get_step(step_idx) for r in self.results] for step_idx in step_idxs]