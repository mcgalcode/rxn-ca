from __future__ import annotations

from ...core.heating import HeatingSchedule
from ...phases.solid_phase_set import MatterPhase
from ..reaction_step_analyzer import ReactionStepAnalyzer, AnalysisMode, AnalysisQuantity
from typing import List, Dict
from dataclasses import dataclass, field

import plotly.graph_objects as go
import numpy as np

@dataclass
class PhaseTraceConfig():

    exclude_phases: List = field(default_factory=list)
    exact_phase_set: List = field(default_factory=list)
    minimum_required_prevalence: float = field(default=0.05)

@dataclass
class PhaseTrace():
    
    name: str
    ys: List = field(default_factory=List)


class PhaseTraceCalculator():
    """A class that stores the result of running a simulation. Keeps track of all
    the steps that the simulation proceeded through, and the set of reactions that
    was used in the simulation.
    """
        
    def __init__(self, step_groups, analyzer: ReactionStepAnalyzer):
        """Initializes a ReactionResult with the reaction set used in the simulation

        Args:
            rxn_set (ScoredReactionSet):
        """
        self.step_analyzer = analyzer
        self._step_groups = step_groups
    
    def get_traces(self, step_analyses: Dict[str, float], config: PhaseTraceConfig) -> List[PhaseTrace]:
        if len(config.exact_phase_set) == 0:
            phases_to_plot = set()
            for analysis in step_analyses:
                phases_to_plot = phases_to_plot.union(set(analysis.keys()))
        else:
            phases_to_plot = set(config.exact_phase_set)
        
        traces = []


        for phase in phases_to_plot:
            ys = np.array([a.get(phase, 0) for a in step_analyses])

            if max(ys) > config.minimum_required_prevalence:
                t = PhaseTrace(name=phase, ys=ys)
                traces.append(t)

        return traces
    
    def get_general_traces(self,
                           trace_config: PhaseTraceConfig,
                           quantity: AnalysisQuantity,
                           mode: AnalysisMode,
                           matter_phases: List[MatterPhase] = [MatterPhase.SOLID, MatterPhase.LIQUID]) -> List[PhaseTrace]:
        analyses = [self.step_analyzer.set_step_group(sg).get_value_general(quantity, mode, include_matter_phases=matter_phases) for sg in self._step_groups]
        return self.get_traces(analyses, trace_config)    

    def get_absolute_elemental_mole_traces(self, trace_config: PhaseTraceConfig) -> List[PhaseTrace]:
        return self.get_general_traces(trace_config, AnalysisQuantity.ELEMENTS, AnalysisMode.ABSOLUTE, matter_phases=None)

    def get_fractional_elemental_traces(self, trace_config: PhaseTraceConfig) -> List[PhaseTrace]:
        return self.get_general_traces(trace_config, AnalysisQuantity.ELEMENTS, AnalysisMode.FRACTIONAL, matter_phases=None)

    def get_absolute_molar_amount_traces(self, trace_config: PhaseTraceConfig) -> List[PhaseTrace]:
        return self.get_general_traces(trace_config, AnalysisQuantity.MOLES, AnalysisMode.ABSOLUTE)
    
    def get_mole_fraction_traces(self, trace_config: PhaseTraceConfig) -> List[PhaseTrace]:
        return self.get_general_traces(trace_config, AnalysisQuantity.MOLES, AnalysisMode.FRACTIONAL)

    def get_absolute_atomic_molar_amount_traces(self, trace_config: PhaseTraceConfig) -> List[PhaseTrace]:
        return self.get_general_traces(trace_config, AnalysisQuantity.ATOMS, AnalysisMode.ABSOLUTE)
    
    def get_fractional_atomic_molar_amount_traces(self, trace_config: PhaseTraceConfig) -> List[PhaseTrace]:
        return self.get_general_traces(trace_config, AnalysisQuantity.ATOMS, AnalysisMode.FRACTIONAL)

    def get_absolute_phase_volume_traces(self, trace_config: PhaseTraceConfig) -> List[PhaseTrace]:
        return self.get_general_traces(trace_config, AnalysisQuantity.VOLUME, AnalysisMode.ABSOLUTE)

    def get_absolute_mass_traces(self, trace_config: PhaseTraceConfig) -> List[PhaseTrace]:
        return self.get_general_traces(trace_config, AnalysisQuantity.MASS, AnalysisMode.ABSOLUTE)
    
    def get_mass_fraction_traces(self, trace_config: PhaseTraceConfig) -> List[PhaseTrace]:
        return self.get_general_traces(trace_config, AnalysisQuantity.MASS, AnalysisMode.FRACTIONAL)
        