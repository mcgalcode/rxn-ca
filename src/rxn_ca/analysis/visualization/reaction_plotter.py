from __future__ import annotations

from ..bulk_reaction_analyzer import BulkReactionAnalyzer
from ..reaction_step_analyzer import AnalysisMode, AnalysisQuantity
from .phase_trace_calculator import PhaseTraceCalculator, PhaseTraceConfig, PhaseTrace
from .layout import RxnCALayout
from .rip_plotter import RIPPlotter
import plotly.graph_objects as go

from typing import List, Dict

class ReactionPlotter():
    """A class that stores the result of running a simulation. Keeps track of all
    the steps that the simulation proceeded through, and the set of reactions that
    was used in the simulation.
    """
    
    def __init__(self,
                 bulk_analyzer: BulkReactionAnalyzer,
                 trace_config: PhaseTraceConfig = PhaseTraceConfig(),
                 include_heating_trace: bool = False,
                 rip_config: Dict = None,
                 phase_colors: Dict = None):
        """Initializes a ReactionResult with the reaction set used in the simulation

        Args:
            rxn_set (ScoredReactionSet):
        """
        self.bulk_analyzer = bulk_analyzer
        self.trace_config = trace_config
        self.trace_calculator = PhaseTraceCalculator(
            bulk_analyzer.loaded_step_groups,
            bulk_analyzer.step_analyzer,
        )
        self.include_heating_trace = include_heating_trace
        self.layout = RxnCALayout(self.bulk_analyzer.get_step_size(), self.bulk_analyzer.heating_schedule)
        self.rip_config = rip_config
        self.phase_colors = phase_colors

    def get_heating_trace(self):
        step_size = self.bulk_analyzer.get_step_size()
        heating_xs, heating_ys = self.bulk_analyzer.heating_schedule.get_xy_for_plot(step_size=step_size)
        return go.Scatter(
            name="Temperature",
            x=heating_xs,
            y=heating_ys,
            mode='lines',
            yaxis='y2',
            line = dict(color='crimson', width=4, dash='dash')
        )
    
    def _get_plotly_trace(self):
        return go.Scatter(
            mode='lines',
            line=dict(width=4)
        )

    def _get_plotly_phase_trace(self, t: PhaseTrace):
        default_trace = self._get_plotly_trace()
        default_trace.update(
                name=t.name,
                x=self.bulk_analyzer.loaded_step_idxs,
                y=t.ys,            
        )
        if self.phase_colors is not None:
            default_trace.line.update(color=self.phase_colors.get(t.name))

        return default_trace
    
    def _get_rip_trace(self, pt: PhaseTrace, plotly_trace: go.Scatter):
        freq = 5
        xs = self.bulk_analyzer.loaded_step_idxs[::freq]
        ys = pt.ys[::freq]
        if pt.name in self.rip_config.get("reactants"):
            mdict = dict(symbol= "circle", size=12)
        elif pt.name in self.rip_config.get("products"):
            mdict = dict(symbol= "diamond", size=12)
        else:
            mdict = dict(symbol= "x", size=12)        
        
        mdict.update(color=plotly_trace.line.color)
        rip_trace = go.Scatter(
            mode="markers",
            x=xs,
            y=ys,
            marker=mdict,
            showlegend=False,
            name=pt.name,
            hoverinfo='skip'
        )
        return rip_trace
    
    def _get_basic_phase_trace_fig(self, title, y_axis, phase_traces: List[PhaseTrace]):
        fig = self.layout.get_plotly_fig(
            y_axis,
            title,
        )

        fig.layout.xaxis.update(autorange=False, range=(0, self.bulk_analyzer.last_loaded_step_idx))

        for t in phase_traces:
            plotly_trace = self._get_plotly_phase_trace(t)
            fig.add_trace(plotly_trace)
            if self.rip_config is not None:
                rip_trace = self._get_rip_trace(t, fig.data[-1])
                fig.add_trace(rip_trace)
        
        if self.include_heating_trace:
            fig.add_trace(self.get_heating_trace())

        return fig

    def plot_value(self, 
                   quantity: AnalysisQuantity,
                   mode: AnalysisMode,
                   title: str,
                   ylabel: str) -> None:
        phase_traces = self.trace_calculator.get_general_traces(self.trace_config, quantity, mode)

        fig = self._get_basic_phase_trace_fig(
            title,
            ylabel,
            phase_traces
        )

        if mode == AnalysisMode.FRACTIONAL:
            fig.layout.yaxis.update(range=(0, 1.0))
        
        if quantity != AnalysisQuantity.ELEMENTS and mode == AnalysisMode.FRACTIONAL and self.rip_config is not None:
            all_phases = [t.name for t in phase_traces]
            reactants = self.rip_config.get("reactants")
            products = self.rip_config.get("products")
            impurities = set(all_phases) - set(reactants) - set(products)
            rip_generator = RIPPlotter()
            rip_traces = rip_generator.get_rip_traces(reactants, impurities, products, self.bulk_analyzer.loaded_step_idxs, phase_traces)
            for rt in rip_traces[::-1]:
                fig.add_trace(rt)

        fig.data = fig.data[::-1]
        
        return fig


    def plot_elemental_amounts(self) -> None:
        return self.plot_value(
            AnalysisQuantity.ELEMENTS,
            AnalysisMode.ABSOLUTE,
            "Moles of Element",
            "Molar Elemental Amts. vs time step"
        )
    
    def plot_elemental_fractions(self) -> None:
        return self.plot_value(
            AnalysisQuantity.ELEMENTS,
            AnalysisMode.FRACTIONAL,
            "Fraction",
            "Elemental Fraction"
        )

    def plot_molar_phase_fractions(self) -> None:
        return self.plot_value(
            AnalysisQuantity.MOLES,
            AnalysisMode.FRACTIONAL,
            "Fraction",
            "Molar Fraction"
        )

    def plot_molar_phase_amounts(self) -> None:
        return self.plot_value(
            AnalysisQuantity.MOLES,
            AnalysisMode.ABSOLUTE,
            "Amt.",
            "Molar Amts."
        )
        
    def plot_phase_volumes(self):
        return self.plot_value(
            AnalysisQuantity.VOLUME,
            AnalysisMode.ABSOLUTE,
            "Amt.",
            "Volume"
        )
    
    def plot_phase_masses(self):
        return self.plot_value(
            AnalysisQuantity.MASS,
            AnalysisMode.ABSOLUTE,
            "Amt.",
            "Mass"
        ) 

    def plot_mass_fractions(self):
        return self.plot_value(
            AnalysisQuantity.MASS,
            AnalysisMode.FRACTIONAL,
            "Fraction",
            "Mass"
        )


        