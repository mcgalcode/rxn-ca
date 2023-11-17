import matplotlib.pyplot as plt
import plotly.graph_objects as go
from pylattica.discrete import DiscreteResultAnalyzer

from ..phases.solid_phase_set import SolidPhaseSet
from ..core.reaction_result import ReactionResult
from ..reactions.scored_reaction_set import ScoredReactionSet

from .reaction_step_analyzer import ReactionStepAnalyzer



class ReactionAnalyzer(DiscreteResultAnalyzer):
    """A class that stores the result of running a simulation. Keeps track of all
    the steps that the simulation proceeded through, and the set of reactions that
    was used in the simulation.
    """


    def __init__(self, result: ReactionResult):
        """Initializes a ReactionResult with the reaction set used in the simulation

        Args:
            rxn_set (ScoredReactionSet):
        """
        self.rxn_set: ScoredReactionSet = result.rxn_set
        self.phases = self.rxn_set.phases
        self.heating_schedule = result.heating_schedule
        self.step_analyzer = ReactionStepAnalyzer(self.rxn_set.phases)
        super().__init__(result)

    def get_choices_at(self, step_no: int, top: int = None, exclude_ids = True) -> None:
        data = self.step_analyzer.get_reaction_choices(self.result.get_step(step_no))
        names = list(data.keys())
        values = list(data.values())
        zipped = list(zip(names, values))
        zipped.sort(key = lambda x: -x[1])

        if exclude_ids:
            filtered = list(filter(lambda item: not self.rxn_set.get_rxn_by_str(item[0]).is_identity, zipped))
        else:
            filtered = zipped

        if top is None:
            top = len(filtered)

        return filtered[0:top]

    def show_choices_at(self, step_no: int, top: int = None, exclude_ids = True) -> None:
        choices = self.get_choices_at(step_no, top, exclude_ids)
        for choice in choices:
            print(choice[0])
            print(f'Competitiveness: {self.rxn_set.get_rxn_by_str(choice[0]).competitiveness}')
            if self.rxn_set.get_rxn_by_str(choice[0]).original_rxn is not None:
                print(f'eV / atom: {self.rxn_set.get_rxn_by_str(choice[0]).original_rxn.energy_per_atom}')
            print(f'Count: {choice[1]}')
            print('--------------')

    def plot_choices_at(self, step_no: int, top: int = None, exclude_ids = True) -> None:
        choices = self.get_choices_at(step_no, top, exclude_ids)
        sorted_names = list(map(lambda x: x[1], choices[0:top]))
        sorted_values = list(map(lambda x: x[0], choices[0:top]))
        fig, axs = plt.subplots(1, 1, figsize=(10, 10))
        axs.bar(sorted_names, sorted_values)
        fig.suptitle(f'Chosen Reactions @ Step {step_no}')

        axs.set_xticklabels(sorted_names, rotation = 60)

        fig.show()


    def plot_elemental_amounts(self) -> None:
        fig = go.Figure()
        fig.update_layout(width=800, height=800, title="Molar Elemental Amount vs time step")
        fig.update_yaxes(title="# of Moles")
        fig.update_xaxes(title="Simulation Step")

        elements = list(self.step_analyzer.elemental_composition(self._result.first_step).keys())
        traces = []
        
        step_idxs, steps = self._get_steps_to_plot()
        amounts = [self.step_analyzer.elemental_composition(s) for s in steps]
        for el in elements:
            ys = [a.get(el, 0) for a in amounts]
            traces.append((step_idxs, ys, el))


        for t in traces:
            fig.add_trace(go.Scatter(name=t[2], x=t[0], y=t[1], mode='lines'))
        
        fig.update_yaxes(range=[0, 100])
        fig.update_layout(showlegend=True)
        fig.show()


    def plot_mole_fractions(self, min_prevalence=0.01) -> None:
        """In a Jupyter Notebook environment, plots the phase prevalence traces for the simulation.

        Returns:
            None:
        """

        fig = go.Figure()
        fig.update_layout(width=800, height=800, title="Prevalence by Simulation Step")
        fig.update_yaxes(title="Prevalence")

        traces = []
        step_idxs, steps = self._get_steps_to_plot()
        fig.update_xaxes(range=[0, step_idxs[-1]], title="Simulation Step")
        molar_breakdowns = [self.step_analyzer.molar_fractional_breakdown(step) for step in steps]

        phases = set()
        for bd in molar_breakdowns:
            phases = phases.union(set(bd.keys()))
            
        for phase in phases:
            if phase is not SolidPhaseSet.FREE_SPACE:
                ys = []
                xs = []
                for s, mb in zip(step_idxs, molar_breakdowns):
                    if self.heating_schedule.temp_at(s) <= self.phases.get_melting_point(phase):
                        ys.append(mb.get(phase, 0))
                        xs.append(s)
                traces.append((xs, ys, phase))

        filtered_traces = [t for t in traces if max(t[1]) > min_prevalence]

        for t in filtered_traces:
            fig.add_trace(go.Scatter(name=t[2], x=t[0], y=t[1], mode='lines'))
        fig.update_layout(showlegend=True)

        fig.show()

    def plot_molar_phase_amounts(self, min_prevalence=0.01) -> None:
        """In a Jupyter Notebook environment, plots the phase prevalence traces for the simulation.

        Returns:
            None:
        """

        fig = go.Figure()
        fig.update_layout(width=800, height=800, title="Absolute Molar Prevalence by Simulation Step")
        fig.update_yaxes(title="# of Moles")

        traces = []
        step_idxs, steps = self._get_steps_to_plot()
        fig.update_xaxes(range=[0, step_idxs[-1]], title="Simulation Step")
        molar_breakdowns = [self.step_analyzer.molar_breakdown(step) for step in steps]

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

        fig.update_layout(showlegend=True)

        fig.show()

    def plot_phase_volumes(self, min_prevalence=0.01) -> None:
        """In a Jupyter Notebook environment, plots the phase prevalence traces for the simulation.

        Returns:
            None:
        """

        fig = go.Figure()
        fig.update_layout(width=800, height=800, title="Absolute Phase Volumes by Simulation Step")
        fig.update_yaxes(title="Volume (cc)")

        traces = []
        step_idxs, steps = self._get_steps_to_plot()
        fig.update_xaxes(range=[0, step_idxs[-1]], title="Simulation Step")
        volume_breakdowns = [self.step_analyzer.phase_volumes(step) for step in steps]

        phases = set()
        for bd in volume_breakdowns:
            phases = phases.union(set(bd.keys()))
            
        for phase in phases:
            if phase is not SolidPhaseSet.FREE_SPACE:
                ys = [vb.get(phase, 0) for vb in volume_breakdowns]
                traces.append((step_idxs, ys, phase))

        filtered_traces = [t for t in traces if max(t[1]) > min_prevalence]

        for t in filtered_traces:
            fig.add_trace(go.Scatter(name=t[2], x=t[0], y=t[1], mode='lines'))
        fig.update_layout(showlegend=True)

        fig.show()
    
    def plot_total_volume(self):
        fig = go.Figure()

        step_idxs, steps = self._get_steps_to_plot()
        volumes = [self.step_analyzer.total_volume(step) for step in steps]

        fig.update_xaxes(range=[0, step_idxs[-1]], title="Simulation Step")
        fig.update_layout(width=800, height=800, title="Total Simulation Volume")
        fig.update_yaxes(range=[0, max(volumes) + max(volumes) * 0.05], title="Volume (cc)")
        fig.add_trace(go.Scatter(name='Volume', x=step_idxs, y=volumes, mode='lines'))

        fig.show()      

    def final_mol_fractions(self):
        return self.all_mole_fractions_at(len(self._result))

    def all_mole_fractions_at(self, step):
        fracs = {}
        for phase in self.all_phases():
            fracs[phase] = self.step_analyzer.mole_fraction(self._result.get_step(step), phase)

        return fracs

    def mole_fraction_at(self, step, phase):
        return self.step_analyzer.mole_fraction(self.steps[step - 1], phase)

    def as_dict(self):
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "steps": [s.as_dict() for s in self.steps],
            "rxn_set": self.rxn_set.as_dict(),
            "phase_set": self.phase_set.as_dict()
        }