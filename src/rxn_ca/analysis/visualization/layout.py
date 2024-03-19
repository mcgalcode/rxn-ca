from __future__ import annotations

import plotly.graph_objects as go
from plotly.graph_objs.layout import YAxis,XAxis

from ...core.heating import HeatingSchedule

class RxnCALayout():
    """A class that stores the result of running a simulation. Keeps track of all
    the steps that the simulation proceeded through, and the set of reactions that
    was used in the simulation.
    """
    
    def __init__(self, step_size: int, heating_sched: HeatingSchedule):
        """Initializes a ReactionResult with the reaction set used in the simulation

        Args:
            rxn_set (ScoredReactionSet):
        """
        self.heating_schedule = heating_sched
        self.step_size = step_size
    
    def get_layout(self, y_label, title, **layout_kwargs):
        default_kwargs = dict(title={
                'text' : title,
                'x': 0.42,
                'y': 0.92,
            },
            paper_bgcolor='rgba(0,0,0,0)', 
            plot_bgcolor='rgba(0,0,0,0)',
            width=900, height=800,
            xaxis=XAxis(
                mirror=True,
                showticklabels=False,
                showline=True,
                linecolor='black',
                zeroline=True,
            ),        
            yaxis2 = YAxis(
                title=dict(
                    text="Temp (K)",
                    font=dict(
                        color='crimson'
                    )       
                ),
                # gridcolor='red',
                tickfont=dict(
                    color='crimson'
                ),
                overlaying= 'y', 
                side= 'right',
                ticks='outside',
                showline=True,  
                linecolor='black',      
            ),
            yaxis=YAxis(
                title=dict(
                    text=y_label,
                ),
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

    def get_plotly_fig(self,
                       y_label,
                       title,
                       include_heating_trace=True,
                       use_heating_xaxis=False,
                       **layout_kwargs):
        if use_heating_xaxis:
            x_label = "Temperature (K)"
        else:
            x_label = "Reaction Coordinate (arb. units)"

        
        layout = self.get_layout(y_label, title, **layout_kwargs)

        fig = go.Figure(layout=layout)

        fig.update_layout(
            yaxis=dict(
                title=dict(
                    text=y_label,
                ),
                ticks='outside',
                gridcolor='black',
                showgrid=True,
                showticklabels=True,
                side='left',
                showline=True,    
                linecolor='black',
            ),
        )

        if include_heating_trace:
            fig.update_layout(
                width=1000,
                legend=dict(
                    x=1.2,
                    y=0.6,
                    borderwidth=2
                ),
            )

        # if use_heating_xaxis:
        #     xs = []
        #     xlabels = []
        #     curr_x = 0
        #     last_temp = None
        #     for step in self.heating_schedule.steps:
        #         if step.temperature != last_temp:
        #             xs.append(curr_x)
        #             xlabels.append(step.temperature)
        #             last_temp = step.temperature
        #         curr_x = curr_x + step.duration * self.step_size


        #     fig.update_layout(
        #         xaxis = dict(
        #             tickmode = 'array',
        #             tickvals = xs,
        #             ticktext = xlabels,
        #             showticklabels = True,
        #             showline = True,
        #             showgrid=False,
        #             title=x_label
        #         )
        #     )
        # else:
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