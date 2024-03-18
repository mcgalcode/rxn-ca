from __future__ import annotations

import plotly.graph_objects as go
from plotly.graph_objs.layout import YAxis,XAxis,Margin

from ..phases.solid_phase_set import SolidPhaseSet
from ..core.reaction_result import ReactionResult
from ..core.heating import HeatingSchedule
from .reaction_step_analyzer import ReactionStepAnalyzer

from ..computing.schemas.ca_result_schema import RxnCAResultDoc

from typing import Tuple, List

from pymatgen.core.composition import Composition
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
        return cls(doc.results, doc.reaction_library.phases, doc.recipe.heating_schedule)
    
    @classmethod
    def from_result_doc(cls, doc: RxnCAResultDoc) -> BulkReactionAnalyzer:
        return cls(doc.results, doc.reaction_library.phases, doc.recipe.heating_schedule)
    
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

    def _get_trace(self, name, xs, ys, legend="legend1"):
        return go.Scatter(name=name, x=xs, y=ys, mode='lines', line=dict(width=4), legend=legend)

    def plot_elemental_amounts(self, num_points=None, **layout_kwargs) -> None:
        elements = list(self.step_analyzer.get_molar_elemental_composition(self.results[0].first_step).keys())
        traces = []
        
        step_idxs, step_groups = self._get_step_groups(num_points=num_points)
        amounts = [self.step_analyzer.get_molar_elemental_composition(sg) for sg in step_groups]
        for el in elements:
            ys = [a.get(el, 0) for a in amounts]
            traces.append(self._get_trace(el, step_idxs, ys))

        max_x = max(set().union([t.x for t in traces]))
        fig = self._get_plotly_fig(
            "Moles of Element",
            "Molar Elemental Amts. vs time step",
            max_x,
            **layout_kwargs
        )

        for t in traces:
            fig.add_trace(t)

        fig.show()
    
    def plot_elemental_fractions(self, **layout_kwargs) -> None:
        elements = list(self.step_analyzer.get_fractional_elemental_composition(self.results[0].first_step).keys())
        traces = []
        
        step_idxs, step_groups = self._get_step_groups()
        amounts = [self.step_analyzer.get_fractional_elemental_composition(sg) for sg in step_groups]
        for el in elements:
            ys = [a.get(el, 0) for a in amounts]
            traces.append(self._get_trace(el, step_idxs, ys))

        max_x = max(set().union([t.x for t in traces]))
        fig = self._get_plotly_fig(
            "El. Fraction",
            "Elemental Fractions vs time step",
            max_x,
            **layout_kwargs
        )

        for t in traces:
            fig.add_trace(t)

        fig.show()
    
    def get_elemental_amounts_at(self, step_no):
        step_group = self.get_steps(step_no)
        return self.step_analyzer.get_molar_elemental_composition(step_group)

    def get_layout(self, x_label, y_label, title, max_x, max_y=None, **layout_kwargs):
        default_kwargs = dict(title={
                'text' : title,
                'x': 0.42,
                'y': 0.92,
            },
            paper_bgcolor='rgba(0,0,0,0)', 
            plot_bgcolor='rgba(0,0,0,0)',
            width=900, height=800,
            xaxis=XAxis(
                range=[0, max_x],
                ticks='',
                mirror=True,
                showticklabels=False,
                showline=True,
                linecolor='black',
            ),        
            yaxis2 = YAxis(
                title=dict(
                    text="Temp (K)",
                    font=dict(
                        color='crimson'
                    )       
                ),
                gridcolor='red',
                tickfont=dict(
                    color='crimson'
                ),
                overlaying= 'y', 
                side= 'right',
                ticks='outside',
                showline=True,  
                linecolor='black',      
            ),
            yaxis=dict(
                title=dict(
                    text=y_label,
                ),
                showgrid=False,
                range=(0, max_y),
                side='left',
                showticklabels=False,
                ticks='',
                showline=True,    
                linecolor='black',
                mirror=True   
            ),
            legend=dict(
                borderwidth=2
            ),
            font=dict(
                family="Lato",
                size=18,
            )
        )

        updated = {**default_kwargs, **layout_kwargs}

        return go.Layout(
            **updated
        )

    def _get_plotly_fig(self,
                        y_label,
                        title,
                        max_x,
                        include_heating_trace=True,
                        use_heating_xaxis=False,
                        show_y_ticks_and_grids=False,
                        max_y=None,
                        **layout_kwargs):
        if use_heating_xaxis:
            x_label = "Temperature (K)"
        else:
            x_label = "Reaction Coordinate (arb. units)"
        layout = self.get_layout(x_label, y_label, title, max_x, max_y=max_y, **layout_kwargs)

        fig = go.Figure(layout=layout)

        if show_y_ticks_and_grids:
            fig.update_layout(
                yaxis=dict(
                    title=dict(
                        text=y_label,
                    ),
                    ticks='outside',
                    gridcolor='black',
                    showgrid=True,
                    showticklabels=True,
                    range=(0, None),
                    side='left',
                    showline=True,    
                    linecolor='black',
                ),
            )

        if include_heating_trace:
            fig.add_trace(self.get_heating_trace())
            fig.update_layout(
                width=1000,
                legend=dict(
                    x=1.2,
                    y=0.6,
                    borderwidth=2
                ),
            )

        if use_heating_xaxis:
            step_size = len(self.results[0].first_step.all_site_states()) / 100
            xs = []
            xlabels = []
            curr_x = 0
            last_temp = None
            for step in self.heating_schedule.steps:
                if step.temp != last_temp:
                    xs.append(curr_x)
                    xlabels.append(step.temp)
                    last_temp = step.temp
                curr_x = curr_x + step.duration * step_size


            fig.update_layout(
                xaxis = dict(
                    tickmode = 'array',
                    tickvals = xs,
                    ticktext = xlabels,
                    showticklabels = True,
                    showline = True,
                    showgrid=False,
                    title=x_label
                )
            )
        else:
            fig.add_annotation(
                dict(
                    font=dict(
                        color="black",
                        family="Lato",
                        size=22
                    ),
                    showarrow=False,
                    x=0.49,
                    y=-0.08,
                    text=x_label,
                    textangle=0,
                    xref="paper",
                    yref="paper"
                )
            )
        
        return fig
    
    def get_heating_trace(self):
        # NOTE: This 100 is hardcoded to match the `record_every` parameter in the async runner in pylattica
        step_size = len(self.results[0].first_step.all_site_states()) / 100
        heating_xs, heating_ys = self.heating_schedule.get_xy_for_plot(step_size=step_size)
        return go.Scatter(
            name="Temperature",
            x=heating_xs,
            y=heating_ys,
            mode='lines',
            yaxis='y2',
            line = dict(color='crimson', width=4, dash='dash')
        )
    
    def get_mole_fraction_traces(self, min_prevalence=0.01, phases=None):
        """In a Jupyter Notebook environment, plots the phase prevalence traces for the simulation.

        Returns:
            None:
        """

        step_idxs, step_groups = self._get_step_groups()

        traces = []
        mole_fraction_breakdowns = [self.step_analyzer.get_all_mole_fractions(sg) for sg in step_groups]

        if phases is None:
            phases = set()
            for bd in mole_fraction_breakdowns:
                phases = phases.union(set(bd.keys()))
        
        for phase in phases:
            ys = [mb.get(phase, 0) for mb in mole_fraction_breakdowns]
            traces.append((step_idxs, ys, phase))

        filtered_traces = [t for t in traces if max(t[1]) > min_prevalence]

        go_traces = [self._get_trace(t[2], t[0], t[1]) for t in filtered_traces]
        return go_traces
    
    def molar_fractional_breakdown(self, step_no, include_melted: bool = True):
        return self.step_analyzer.get_all_absolute_molar_amounts(self.get_steps(step_no), include_melted=include_melted)
    
    def get_final_molar_breakdown(self):
        _, steps = self._get_step_groups()
        return self.step_analyzer.get_all_mole_fractions(steps[-1])

    def get_all_absolute_molar_amounts(self, step_no: int, include_melted: bool = True):
        return self.step_analyzer.get_all_absolute_molar_amounts(self.get_steps(step_no), include_melted=include_melted)

    def plot_mole_fractions(self,
                            min_prevalence=0.01,
                            **layout_kwargs) -> None:
        """In a Jupyter Notebook environment, plots the phase prevalence traces for the simulation.

        Returns:
            None:
        """

        traces = self.get_mole_fraction_traces(min_prevalence)
        max_x = max(set().union([t.x for t in traces]))

        fig = self._get_plotly_fig(
            "Mole Fraction",
            "Molar Phase Fractions",
            max_x,
            show_y_ticks_and_grids=True,
            **layout_kwargs
        )

        for t in traces:
            fig.add_trace(t)

        return fig
    
    def get_steps(self, step_no):
        return [r.get_step(step_no) for r in self.results]
    
    def get_final_steps(self):
        return [r.last_step for r in self.results]

    def get_first_steps(self):
        return [r.first_step for r in self.results]

    def all_phases_present(self):
        step_idxs, step_groups = self._get_step_groups()
        sgs = [self.step_analyzer.phases_present(sg) for sg in step_groups]
        all_phases = set()
        for sg in sgs:
            all_phases.update(sg)
        return list(all_phases)

    def get_molar_phase_traces(self, min_prevalence=0.01, xrd_adjust=True, phases=None, legend="legend1") -> None:
        """In a Jupyter Notebook environment, plots the phase prevalence traces for the simulation.

        Returns:
            None:
        """
        step_idxs, step_groups = self._get_step_groups()
        
        traces = []
        molar_breakdowns = [self.step_analyzer.get_all_absolute_molar_amounts(step_group) for step_group in step_groups]

        if phases is None:
            phases = set()
            for bd in molar_breakdowns:
                phases = phases.union(set(bd.keys()))
            
        for phase in phases:
            if phase is not SolidPhaseSet.FREE_SPACE:
                ys = np.array([mb.get(phase, 0) for mb in molar_breakdowns])
            
                if xrd_adjust:
                    comp = Composition(phase)
                    ys = ys * comp.num_atoms

                if max(ys) > min_prevalence:
                    t = self._get_trace(phase, step_idxs, ys, legend=legend)
                    traces.append(t)

        return traces
    

    def plot_molar_phase_amounts(self,
                                 min_prevalence=0.01,
                                 xrd_adjust=False,
                                 phases=None,
                                 **layout_kwargs) -> None:
        """In a Jupyter Notebook environment, plots the phase prevalence traces for the simulation.

        Returns:
            None:
        """
        traces = self.get_molar_phase_traces(min_prevalence, xrd_adjust=xrd_adjust, phases=phases)

        if xrd_adjust:
            ylabel = "# of Moles (weighted by # atoms)"
        else:
            ylabel = "# of Moles"

        fig = self._get_plotly_fig(
            ylabel,
            "Absolute Molar Prevalence by Simulation Step",
            None,
            **layout_kwargs
        )

        for t in traces:
            fig.add_trace(t)

        return fig
    
    def get_mole_trace(self, phase_name, xrd_adjust=False):
        step_idxs, step_groups = self._get_step_groups()
        molar_breakdowns = [self.step_analyzer.get_all_absolute_molar_amounts(step_group) for step_group in step_groups]

        ys = np.array([mb.get(phase_name, 0) for mb in molar_breakdowns])
            
        if xrd_adjust:
            comp = Composition(phase_name)
            ys = ys * comp.num_atoms

        return ys
    
    def generate_rip_plots(self,
                           reactants,
                           products,
                           **layout_kwargs):
        all_phases = self.all_phases_present()
        impurities = set(all_phases) - set(reactants) - set(products)
        groups = {
            "Precursor": reactants,
            "Impurities": impurities,
            "Products": products
        }

        details = {
            "Precursor": {
                "fillcolor": "rgb(227, 228, 229)",
                "line": {
                    "color": "white",
                    "width": 4
                }
            },
            "Impurities": {
                "fillcolor": "rgb(248,227,237)",
                "line": {
                    "color": "white",
                    "width": 4
                }
            },
            "Products": {
                "fillcolor": "rgb(201, 225, 215)",
                "line": {
                    "color": "white",
                    "width": 4
                }
            }
        }
        self.plot_groups(groups, group_details=details, **layout_kwargs)
    
    def plot_groups(self,
                    group_defs,
                    group_details = {},
                    **layout_kwargs):
        step_idxs, step_groups = self._get_step_groups()

        traces = []
        molar_breakdowns = [self.step_analyzer.get_all_mass_fractions(sg) for sg in step_groups]


        phases = set()
        for bd in molar_breakdowns:
            phases = phases.union(set(bd.keys()))

        traces = {}

        for group_name, group_members in group_defs.items():
            group_trace = np.array([0 for _ in molar_breakdowns])

            for phase in phases:
                if phase in group_members and phase is not SolidPhaseSet.FREE_SPACE:
                    ys = [mb.get(phase, 0) for mb in molar_breakdowns]
                    group_trace = group_trace + ys
            
            traces[group_name] = group_trace

        go_traces = []

        for trace_name, trace_data in traces.items():
            tr = go.Scatter(
                name=trace_name, 
                x=step_idxs,
                y=trace_data,
                mode='lines',
                stackgroup='one',
                **group_details.get(trace_name, {})
            )

            go_traces.append(tr)

        fig = self._get_plotly_fig(
            "Mass Fraction",
            "Mass Fraction Prevalence by Simulation Step",
            None,
            max_y=1.0,
            **layout_kwargs
        )

        fig.add_traces(go_traces)
        for t in self.get_mass_fraction_traces():
            fig.add_trace(self._get_trace(t[2], t[0], t[1]))
            

        fig.show()
        
        

    def plot_phase_volumes(self, **layout_kwargs):
        """In a Jupyter Notebook environment, plots the phase prevalence traces for the simulation.

        Returns:
            None:
        """

        step_idxs, step_groups = self._get_step_groups()

        fig = self._get_plotly_fig(
            "Volume (arb. units)",
            "Phase Volume by Simulation Step",
            step_idxs[-1],
            **layout_kwargs
        )
        
        traces = []
        vol_breakdowns = [self.step_analyzer.get_all_absolute_phase_volumes(step_group) for step_group in step_groups]

        phases = set()
        for bd in vol_breakdowns:
            phases = phases.union(set(bd.keys()))
            
        for phase in phases:
            if phase != SolidPhaseSet.FREE_SPACE:
                ys = [mb.get(phase, 0) for mb in vol_breakdowns]
                traces.append((step_idxs, ys, phase))

        filtered_traces = [t for t in traces if max(t[1]) > 0.1]

        for t in filtered_traces:
            fig.add_trace(self._get_trace(t[2], t[0], t[1]))

        fig.show()

    def get_mass_traces(self):
        step_idxs, step_groups = self._get_step_groups()
        
        traces = []
        mass_breakdowns = [self.step_analyzer.get_absolute_phase_masses(step_group) for step_group in step_groups]

        phases = set()
        for bd in mass_breakdowns:
            phases = phases.union(set(bd.keys()))
            
        for phase in phases:
            if phase != SolidPhaseSet.FREE_SPACE:
                ys = [mb.get(phase, 0) for mb in mass_breakdowns]
                traces.append((step_idxs, ys, phase))

        filtered_traces = [t for t in traces if max(t[1]) > 0.1]

        return filtered_traces
    
    def plot_phase_masses(self, **layout_kwargs):
        """In a Jupyter Notebook environment, plots the phase prevalence traces for the simulation.

        Returns:
            None:
        """

        step_idxs, _ = self._get_step_groups()

        fig = self._get_plotly_fig(
            "Mass (arb. units)",
            "Phase Mass by Simulation Step",
            step_idxs[-1],
            **layout_kwargs
        )

        for t in self.get_mass_traces():
            fig.add_trace(self._get_trace(t[2], t[0], t[1]))

        fig.show()

    def get_mass_fraction_traces(self):
        step_idxs, step_groups = self._get_step_groups()

        traces = []
        mass_breakdowns = [self.step_analyzer.get_all_mass_fractions(step_group) for step_group in step_groups]

        phases = set()
        for bd in mass_breakdowns:
            phases = phases.union(set(bd.keys()))
            
        for phase in phases:
            if phase != SolidPhaseSet.FREE_SPACE:
                ys = [mb.get(phase, 0) for mb in mass_breakdowns]
                traces.append((step_idxs, ys, phase))

        filtered_traces = [t for t in traces if max(t[1]) > 0.01] 
        return filtered_traces      

    def plot_mass_fractions(self, **layout_kwargs):
        """In a Jupyter Notebook environment, plots the phase prevalence traces for the simulation.

        Returns:
            None:
        """

        step_idxs, _ = self._get_step_groups()

        fig = self._get_plotly_fig(
            "Mass Fraction",
            "Phase Mass Fraction by Simulation Step",
            step_idxs[-1],
            show_y_ticks_and_grids=True,
            **layout_kwargs
        )

        for t in self.get_mass_fraction_traces():
            fig.add_trace(self._get_trace(t[2], t[0], t[1]))

        fig.show()        

    def as_dict(self):
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "steps": [s.as_dict() for s in self.steps],
            "rxn_set": self.rxn_set.as_dict(),
            "phase_set": self.phase_set.as_dict()
        }

    def _get_step_groups(self, num_points=None) -> Tuple[List[int], List]:
        if self._step_idxs is None or num_points is not None:
            if num_points is None:
                num_points = min(50, self.result_length)
            step_size = max(1, round(self.result_length / num_points))
            if not self._results_loaded:
                [r.load_steps(step_size) for r in self.results]
                self._results_loaded = True
            self._step_idxs = list(range(0, self.result_length, step_size))
            self._step_groups = [[r.get_step(step_idx) for r in self.results] for step_idx in self._step_idxs]

        return self._step_idxs, self._step_groups
        