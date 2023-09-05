import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.graph_objs.layout import YAxis,XAxis,Margin

from ..phases.solid_phase_set import SolidPhaseSet
from ..core.reaction_result import ReactionResult
from ..core.heating import HeatingSchedule
from .bulk_step_analyzer import BulkReactionStepAnalyzer
from .reaction_step_analyzer import ReactionStepAnalyzer

from ..computing.schemas.ca_result_schema import RxnCAResultDoc

from typing import Tuple, List

from pymatgen.core.composition import Composition
import numpy as np
class BulkReactionAnalyzer():
    """A class that stores the result of running a simulation. Keeps track of all
    the steps that the simulation proceeded through, and the set of reactions that
    was used in the simulation.
    """

    @classmethod
    def from_result_doc_file(cls, fname):
        doc = RxnCAResultDoc.from_file(fname)
        return cls(doc.results, doc.reaction_library.phases, doc.recipe.heating_schedule)

    def __init__(self, results: List[ReactionResult], phase_set: SolidPhaseSet, heating_sched: HeatingSchedule):
        """Initializes a ReactionResult with the reaction set used in the simulation

        Args:
            rxn_set (ScoredReactionSet):
        """
        self.bulk_step_analyzer = BulkReactionStepAnalyzer(phase_set)
        self.single_step_analyzer = ReactionStepAnalyzer(phase_set)
        self.heating_schedule = heating_sched

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

    def _get_plotly_fig(self, x_label, y_label, title, max_x):
        layout = go.Layout(
            title=title,
            width=1000, height=800,
            xaxis=XAxis(
                title=x_label,
                range=[0, max_x]
            ),
            yaxis2 = YAxis(
                title="Temp (K)",
                overlaying= 'y', 
                side= 'right',
            ),
            yaxis=dict(
                title=y_label,
                range=(0, None),
                side='left'
            ),
            legend=dict(
                x=1.1,
                y=0.6
            )
        )

        fig = go.Figure(layout=layout)        


        heating_xs, heating_ys = self.heating_schedule.get_xy_for_plot()
        fig.add_trace(go.Scatter(
            name="Temperature",
            x=heating_xs,
            y=heating_ys,
            mode='lines',
            yaxis='y2',
            line = dict(color='royalblue', width=4, dash='dash')))
        
        return fig

    def plot_mole_fractions(self, min_prevalence=0.01) -> None:
        """In a Jupyter Notebook environment, plots the phase prevalence traces for the simulation.

        Returns:
            None:
        """

        step_idxs, step_groups = self._get_step_groups()

        fig = self._get_plotly_fig(
            "Simulation Step",
            "Prevalence",
            "Prevalence by Simulation Step",
            step_idxs[-1]
        )

        traces = []
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

    def plot_molar_phase_amounts(self, min_prevalence=0.01, xrd_adjust=True) -> None:
        """In a Jupyter Notebook environment, plots the phase prevalence traces for the simulation.

        Returns:
            None:
        """
        step_idxs, step_groups = self._get_step_groups()

        fig = self._get_plotly_fig(
            "Simulation Step",
            "# of Moles",
            "Absolute Molar Prevalence by Simulation Step",
            None
        )
        
        traces = []
        molar_breakdowns = [self.bulk_step_analyzer.molar_breakdown(step_group) for step_group in step_groups]

        phases = set()
        for bd in molar_breakdowns:
            phases = phases.union(set(bd.keys()))
            
        for phase in phases:
            if phase is not SolidPhaseSet.FREE_SPACE:
                ys = np.array([mb.get(phase, 0) for mb in molar_breakdowns])
                
                if xrd_adjust:
                    comp = Composition(phase)
                    ys = ys * comp.num_atoms
                    
                traces.append((step_idxs, ys, phase))

        filtered_traces = [t for t in traces if max(t[1]) > min_prevalence]

        for t in filtered_traces:
            fig.add_trace(go.Scatter(name=t[2], x=t[0], y=t[1], mode='lines'))

        fig.show()

    def plot_phase_volumes(self):
        """In a Jupyter Notebook environment, plots the phase prevalence traces for the simulation.

        Returns:
            None:
        """

        step_idxs, step_groups = self._get_step_groups()

        fig = self._get_plotly_fig(
            "Simulation Step",
            "Volume",
            "Phase Volume by Simulation Step",
            step_idxs[-1]
        )
        

        traces = []
        vol_breakdowns = [self.bulk_step_analyzer.phase_volumes(step_group) for step_group in step_groups]

        phases = set()
        for bd in vol_breakdowns:
            phases = phases.union(set(bd.keys()))
            
        for phase in phases:
            if phase is not SolidPhaseSet.FREE_SPACE:
                ys = [mb.get(phase, 0) for mb in vol_breakdowns]
                traces.append((step_idxs, ys, phase))

        filtered_traces = [t for t in traces if max(t[1]) > 0.1]

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