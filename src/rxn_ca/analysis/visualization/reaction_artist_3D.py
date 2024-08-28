from pylattica.visualization.structure_artist import StructureArtist
from pylattica.visualization import DiscreteCellArtist
from pylattica.core import SimulationState, PeriodicStructure
from pylattica.core.periodic_structure import LOCATION, SITE_ID

import numpy as np
import io
import matplotlib.pyplot as plt
from PIL import Image
from matplotlib.lines import Line2D
from ..reaction_step_analyzer import ReactionStepAnalyzer

def formula_to_latex(formula):
    """Convert a chemical formula to a LaTeX formatted string."""
    result = ""
    i = 0
    while i < len(formula):
        if formula[i].isdigit():
            result += '$_{' + formula[i]
            while i + 1 < len(formula) and formula[i + 1].isdigit():
                result += formula[i + 1]
                i += 1
            result += '}$'
        else:
            result += formula[i]
        i += 1
    return result

class ReactionArtist3D(StructureArtist):
    """A helper StructureArtist class for rendering 3D square grids."""

    def __init__(self, structure: PeriodicStructure, cell_artist: DiscreteCellArtist, step_analyzer: ReactionStepAnalyzer):
        """Instantiates the StructureArtist with a structure and the artist used to
        render it's sites.

        Parameters
        ----------
        structure : PeriodicStructure
            The structure to render.
        cell_artist : CellArtist
            The artist for rendering each of the structure's sites.
        """
        self.structure = structure
        self.cell_artist = cell_artist
        self._step_analyzer = step_analyzer

    def _draw_image(self, state: SimulationState, **kwargs):
        shell_only = kwargs.get("shell_only", False)

        size = round(state.size ** (1 / 3))

        shape = [size for _ in range(self.structure.dim)]
        dataset = {}

        dataset["empty"] = np.ones(shape)
        color_cache = {}

        for site in self.structure.sites():
            loc = site[LOCATION]
            if not shell_only or (loc[1] == 0 or loc[0] == size or loc[2] == size):
                site_id = site[SITE_ID]
                site_state = state.get_site_state(site_id)
                color = self.cell_artist.get_color_from_cell_state(site_state)
                color_str = str(color)
                if color_str not in color_cache:
                    color_cache[color_str] = color

                if color_str not in dataset:
                    dataset[color_str] = np.zeros(shape)

                shifted_loc = tuple(int(i) for i in loc)
                dataset[color_str][shifted_loc] = 1
                dataset["empty"][shifted_loc] = 0

        fig = plt.figure(figsize=(18, 10))

        gs = fig.add_gridspec(2, 2, width_ratios=(1,1), height_ratios=[1, 0.12])

        bar_chart = fig.add_subplot(gs[0, 0])

        ax_3d = fig.add_subplot(gs[0, 1], projection="3d")

        self._step_analyzer.set_step_group(state)
        bar_data = self._step_analyzer.get_all_mass_fractions()

        legend = self.cell_artist.get_legend(state)
        colors = self.cell_artist.color_map
        
        formatted_bar_data = { phase: bar_data.get(phase, 0) for phase, _ in self.cell_artist.color_map.items() }
        formatted_bar_data = {formula_to_latex(p): v for p, v in formatted_bar_data.items() }
        phases = list(formatted_bar_data.keys())
        print([legend.get(p) for p in self.cell_artist.color_map.keys() if legend.get(p) is not None])
        bar_colors = [np.array(legend.get(p)) / 255 for p in self.cell_artist.color_map.keys() if legend.get(p) is not None]
        values = list(formatted_bar_data.values())
        bar_chart.bar(phases, values, color=bar_colors)
        bar_chart.set_ylim(0, 0.65)
        # bar_chart.yaxis.set_visible(False)
        bar_chart.tick_params(axis='y', which='both', length=0)
        bar_chart.set_yticklabels([])
        bar_chart.set_ylabel("Phase Prevalence", labelpad=20, fontdict={"family": "Lato", "size": 28})
    
        # Remove top, left, and right spines
        bar_chart.spines['top'].set_visible(False)
        bar_chart.spines['right'].set_visible(False)
        bar_chart.spines['left'].set_visible(False)
        
        # Only keep the bottom spine
        bar_chart.spines['bottom'].set_position('zero')
        bar_chart.tick_params(axis='x', which='major', labelsize=18)

        for color, data in dataset.items():
            if color == "empty":
                colors = [0.8, 0.8, 0.8, 0.2]
                ax_3d.voxels(data, facecolors=colors, edgecolor="k", linewidth=0)
            else:
                # colors = [*list(np.array(color_cache[color]) / 255), 0]
                colors = list(np.array(color_cache[color]) / 255)
                ax_3d.voxels(data, facecolors=colors, edgecolor="k", linewidth=0.25)
        
        ax_3d.axis("off")
        RADIUS = 30.0  # Control this value.
        ax_3d.set_xlim3d(0, RADIUS / 2)
        ax_3d.set_zlim3d(0, RADIUS / 2)
        ax_3d.set_ylim3d(0, RADIUS / 2)   

        legend_ax = fig.add_subplot(gs[1, :])
        legend_ax.axis('off')  # Turn off axis

        # if kwargs.get("show_legend") == True:
        legend_handles = []
        for phase, color in legend.items():
            legend_handles.append(
                Line2D(
                    [0],
                    [0],
                    marker="s",
                    color="w",
                    markerfacecolor=list(np.array(color) / 255),
                    markersize=30,
                    label=phase,
                )
            )

        # Add custom legend to the plot
        legend_font_props = {"family": "Lato", "size": 24}
        legend_ax.legend(
            handles=legend_handles,
            loc="center",
            prop=legend_font_props,
            ncols=5,
            frameon=False,
            bbox_to_anchor=(0.5, 0.2),
        )
        
        plt.tight_layout()

        if kwargs.get("label") is not None:
            x_text, y_text, z_text = -6, -5, 23

            # Add the text
            annotation_font = {
                "size": 24,
                "family": "Lato",
                "color": np.array([194, 29, 63]) / 255,
                "weight": "bold",
            }
            ax_3d.text(
                x_text, y_text, z_text, kwargs.get("label"), fontdict=annotation_font
            )

        buf = io.BytesIO()
        fig.savefig(buf)
        plt.close()
        buf.seek(0)
        img = Image.open(buf)
        return img
