from __future__ import annotations

import plotly.graph_objects as go
import numpy as np

from typing import List

from .phase_trace_calculator import PhaseTrace

class RIPPlotter():
    """A class that stores the result of running a simulation. Keeps track of all
    the steps that the simulation proceeded through, and the set of reactions that
    was used in the simulation.
    """

    def get_rip_traces(self,
                       reactants,
                       impurities,
                       products,
                       xs,
                       value_traces: List[PhaseTrace]):
        groups = {
            "Precursor": reactants,
            "Impurities": impurities,
            "Products": products
        }

        details = {
            "Precursor": {
                "fillcolor": "rgb(240, 240, 240)",
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

        rip_traces = {}

        names_used = set()

        normalization_factors = np.array([sum(g) for g in zip(*[v.ys for v in value_traces])])

        for group_name, group_members in groups.items():
            group_trace = np.array([0 for _ in value_traces[0].ys])

            for v_trace in value_traces:
                if v_trace.name in group_members:
                    names_used.add(v_trace.name)
                    group_trace = group_trace + np.array(v_trace.ys) / normalization_factors
            
            rip_traces[group_name] = group_trace

        go_rip_traces = []

        for trace_name, trace_data in rip_traces.items():
            tr = go.Scatter(
                name=trace_name, 
                x=xs,
                y=trace_data,
                mode='lines',
                stackgroup='one',
                **details.get(trace_name, {})
            )

            go_rip_traces.append(tr)

        return go_rip_traces
        
        
       
        