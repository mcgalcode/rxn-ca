# Step 1
from math import comb
import math
import random
import tqdm
import matplotlib.pyplot as plt
import plotly.graph_objects as go

# Note, try  considering atoms as individual units e.g. Iron(in Fe2O3) so you can count discrete concepts

import numpy as np

from pymatgen.core.composition import Composition

from ..reactions import ScoredReactionSet


class CompletionError(BaseException):

    pass


class KMCCalculator():

    def get_h(self, reaction, current_populations):
        reactants = reaction.reactants
        h = 1
        for reactant in reactants:
            population = current_populations.get(reactant, 0)
            combinations = comb(math.floor(population), max(int(reaction.reactant_stoich(reactant)), 1))
            h = h * combinations

        return h

    def a(self, reaction, current_populations):
        h = self.get_h(reaction, current_populations)
        c = reaction.competitiveness
        return h * c

    def get_tau(self, r1, a_0):
        return (1 / a_0) * math.log(1 / r1)

    def get_mu(self, a_values, r, a_0):
        r2a0 = r * a_0
        total_a = 0
        mu = 0
        for a_val in a_values:
            total_a = total_a + a_val
            mu = mu + 1
            if total_a >= r2a0:
                break

        return mu

class KMCSimulation():

    def __init__(self, rxn_set: ScoredReactionSet):
        self.rxn_set = rxn_set
        self.rxns = rxn_set.reactions
        self.calc = KMCCalculator()

        self._index_map = {}
        for idx, rxn in enumerate(self.rxns):
            self._index_map[str(rxn)] = idx

        self.affected_rxns = {}
        for p in self.rxn_set.phases:
            self.affected_rxns[p] = []
            for r in rxn_set.reactions:
                if r.any_reactants([p]):
                    self.affected_rxns[p].append(str(r))

    def update_populations(self, current_populations, rxn):
        reactants = rxn.reactants
        products = rxn.products
        new_populations = current_populations.copy()

        for reactant in reactants:
            stoich = rxn.reactant_stoich(reactant)
            new_populations[reactant] = max(new_populations.get(reactant, 0) - stoich, 0)

        for product in products:
            stoich = rxn.product_stoich(product)
            new_populations[product] = max(new_populations.get(product, 0) + stoich, 0)

        return new_populations

    def a_value_for(self, rxn):
        return self.a_values[self._index_map[str(rxn)]]

    def step_simulation(self, current_populations, curr_time, rxn_counter):
        a_0 = sum(self.a_values)
        if a_0 == 0:
            raise CompletionError

        r1 = random.random()
        r2 = random.random()

        tau = self.calc.get_tau(r1, a_0)
        mu = self.calc.get_mu(self.a_values, r2, a_0)

        new_time = curr_time + tau
        chosen_rxn = self.rxns[mu - 1]

        new_populations = self.update_populations(current_populations, chosen_rxn)
        new_rxn_counter = rxn_counter + 1

        a_indices_to_update = []
        for r in chosen_rxn.reactants:
            for rxn in self.affected_rxns[r]:
                a_indices_to_update.append(self._index_map[str(rxn)])

        for p in chosen_rxn.products:
            for rxn in self.affected_rxns[p]:
                a_indices_to_update.append(self._index_map[str(rxn)])

        a_indices_to_update = list(set(a_indices_to_update))
        for idx in a_indices_to_update:
            self.a_values[idx] = self.calc.a(self.rxns[idx], new_populations)

        # print(f'Must update {len(a_indices_to_update) / len(self.rxns)} a values')

        # phases_present = sum([1 for _, val in current_populations.items() if val > 0])
        # print(f'{phases_present} phases present')
        # self.a_values = [self.calc.a(rxn, new_populations) for rxn in self.rxns]

        return new_populations, chosen_rxn, new_time, new_rxn_counter

    def run_simulation(self, starting_population, number_of_steps):
        t = 0
        n = 0
        current_pop = starting_population
        populations = [starting_population]
        rxn_choices = []
        self.a_values = [self.calc.a(rxn, starting_population) for rxn in self.rxns]

        for i in tqdm.tqdm(range(number_of_steps)):
            try:
                current_pop, chosen_rxn, t, n = self.step_simulation(current_pop, t, n)
                rxn_choices.append(chosen_rxn)
                populations.append(current_pop)
            except CompletionError:
                print("a_0 was zero")
                break


        return KMCResult(populations, rxn_choices, t, n)


class KMCResult():

    def __init__(self, pops, rxns, time, rxn_coordinate):
        self.populations = pops
        self.rxn_choices = rxns

    def pop_change(self, initial_step, final_step):
        state1 = self.populations[initial_step]
        state2 = self.populations[final_step]
        changes = {}
        for p in state1:
            change = state2[p] - state1[p]
            if change != 0:
                changes[p] = change

        return changes

    def change_at(self, step):
        return self.rxn_choices[step], self.pop_change(step, step + 1)

    def plot_molecule(self, mol, bounds = None):
        if bounds is None:
            bounds = [0, len(self.populations)]

        if bounds[1] > len(self.populations) - 1:
            bounds[1] = len(self.populations) - 1

        if (bounds[1] - bounds[0]) > 1000:
            step_size = int((bounds[1] - bounds[0]) / 1000)
        else:
            step_size = 1
        points = []
        xlabels = []
        for i in range(bounds[0], bounds[1], step_size):
            points.append(self.populations[i].get(mol, 0))
            xlabels.append(i)
        plt.plot(xlabels, points, label=mol)

    def get_molecule_trace(self, mol, bounds = None):
        if bounds is None:
            bounds = [0, len(self.populations)]

        if bounds[1] > len(self.populations) - 1:
            bounds[1] = len(self.populations) - 1

        if (bounds[1] - bounds[0]) > 1000:
            step_size = int((bounds[1] - bounds[0]) / 1000)
        else:
            step_size = 1
        points = []
        xlabels = []
        for i in range(bounds[0], bounds[1], step_size):
            points.append(self.populations[i].get(mol, 0))
            xlabels.append(i)
        return (xlabels, points, mol)

    def plot_molecules(self, mols, bounds = None):
        [self.plot_molecule(mol, bounds=bounds) for mol in mols]
        plt.legend()

    def plot_all(self, bounds=None):
        mols_present = []
        for step in self.populations:
            for mol, amt in step.items():
                if amt > 0:
                    mols_present.append(mol)
        mols_present = list(set(mols_present))

        fig = go.Figure()
        fig.update_layout(width=800, height=800)
        fig.update_yaxes(title="Molecule Count")
        fig.update_xaxes(range=[0, len(self.populations) - 1], title="Simulation Step")

        traces = []
        for mol in mols_present:
            traces.append(self.get_molecule_trace(mol, bounds=bounds))

        for t in traces:
            fig.add_trace(go.Scatter(name=t[2], x=t[0], y=t[1], mode='lines'))

        fig.show()

    def phases_present(self, step_no):
        present = []
        for p, amt in self.populations[step_no].items():
            if amt > 0:
                present.append(p)

        return present

    def elemental_composition(self, step):
        phases = self.phases_present(step)
        elemental_amounts = {}
        total = 0
        for p in phases:
            comp = Composition(p)
            moles = self.populations[step].get(p, 0)
            for el, am in comp.as_dict().items():
                num_moles = moles * am
                if el in elemental_amounts:
                    elemental_amounts[el] += num_moles
                else:
                    elemental_amounts[el] = num_moles
                total += num_moles

        for el, am in elemental_amounts.items():
            elemental_amounts[el] = am / total

        return elemental_amounts

    def population(self, step, show_zero = False):
        pops = self.populations[step]
        to_show = {}
        for phase, amt in pops.items():
            if show_zero:
                to_show[phase] = amt
            elif amt > 0:
                to_show[phase] = amt

        return to_show

    def plot_elemental_amounts(self) -> None:
        fig = go.Figure()
        fig.update_layout(width=800, height=800)
        fig.update_yaxes(title="Relative Prevalence")
        fig.update_xaxes(title="Simulation Step")

        elements = list(self.elemental_composition(0).keys())
        traces = []
        amounts = [self.elemental_composition(s) for s in range(len(self.populations))]
        for el in elements:
            xs = np.arange(len(self.populations))
            ys = [a[el] for a in amounts]
            traces.append((xs, ys, el))

        # filtered_traces = [t for t in traces if max(t[1]) > min_prevalence]

        for t in traces:
            fig.add_trace(go.Scatter(name=t[2], x=t[0], y=t[1], mode='lines'))

        fig.show()

