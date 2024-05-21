from rustworkx import PyDiGraph
from rustworkx.visualization import graphviz_draw
from ...computing.schemas import RxnCAResultDoc
from ...reactions import ReactionLibrary

from rxn_ca.utilities.helpers import is_in_chemsys
from collections import Counter
import plotly.graph_objects as go


class ReactionGraph():
    
    def __init__(self, rxn_lib: ReactionLibrary, result_doc: RxnCAResultDoc):
        self._original_rxn_lib = rxn_lib
        self.phase_set = rxn_lib.phases
        
        r_chosen = []
        for res_choices in result_doc.metadata["rxn_choices"]:
            for choice in res_choices:
                if choice is not None:
                    r_chosen.append(choice)
                    

        self.choices = result_doc.metadata["rxn_choices"]

        self.counter = Counter(r_chosen)
        self._pruned_lib = rxn_lib.get_lib_from_ids(r_chosen)
        self.rxns = self._pruned_lib.get_rxns_at_temp(self._pruned_lib.temps[0])

    def count_rxns(self,
                   products=[],
                   reactants=[]):
        total = 0
        for rid, count in self.counter.items():
            rxn = self.rxns.get_rxn_by_id(rid)
            if not rxn.reactants.issuperset(reactants):
                continue
            if not rxn.products.issuperset(products):
                continue
            
            total += count
        
        return total
    
    def get_rxns(self,
                 products=[],
                 reactants=[],
                 t=300,
                 min_count=0):
        rxns = []
        for rid, count in self.counter.items():
            if count > min_count:
                rxn = self._pruned_lib.get_rxns_at_temp(t).get_rxn_by_id(rid)
                if not rxn.reactants.issuperset(reactants):
                    continue
                if not rxn.products.issuperset(products):
                    continue
                
                rxns.append((rxn, count, rid))
        
        rxns = sorted(rxns, key=lambda p: p[1], reverse=True)
        return rxns
    
    def plot_occurrences(self, rxn_id_groups):
        
        idx_groups = [[] for _ in rxn_id_groups]
        for res in self.choices:
            for idx, choice in enumerate(res):
                for idxs, rxn_ids in zip(idx_groups, rxn_id_groups):
                    if choice in rxn_ids:
                        idxs.append(idx)
        
        fig = go.Figure()
        fig.update_layout(
            paper_bgcolor='rgba(0,0,0,0)', 
            plot_bgcolor='rgba(0,0,0,0)',
            margin={'t':0,'l':0,'b':0,'r':0},
            width=500,
            height=500,
            showlegend=False,
            title=None,
            font=dict(
                family="Lato",
                size=18,
            ),
            xaxis=dict(
                mirror=True,
                showticklabels=False,
                showline=True,
                linecolor='black',
                zeroline=True,
                title=None,
                range=(0, None)
            ),
            yaxis=dict(
                side='left',
                # showticklabels=False,
                ticks='',
                showline=True,    
                linecolor='black',
                mirror=True   
            ),    
        )

        for g in idx_groups:
            fig.add_trace(go.Histogram(x=g))
        fig.update_layout(barmode='overlay')
        # Reduce opacity to see both histograms
        fig.update_traces(opacity=0.75)
        fig.update_layout(xaxis=dict(range=(0,max([max(idxs) for idxs in idx_groups]))))
        return fig


    def build(self,
              reactant_phase_set = None,
              product_phase_set = None,
              product_chemsys = None,
              reactant_chemsys = None,
              phases_to_exclude = None,
              phases_to_include = None,
              required_phase = None,
              exact_phase_set = None,
              min_count = 100):
        
        edges = {}
        nodes = set()


        for rxn_id, count in self.counter.items():
            if count < min_count:
                continue
            rxn = self.rxns.get_rxn_by_id(rxn_id)
            nodes = nodes.union(rxn.reactants)
            nodes = nodes.union(rxn.products)
            for r in rxn.reactants:
                for p in rxn.products:
                    curr = edges.get(r, {})
                    stoich_ratio = rxn.reactant_stoich_fraction(r) * rxn.product_stoich_fraction(p)
                    density_ratio = self.phase_set.get_density(p) / self.phase_set.get_density(r)
                    value = count * stoich_ratio #* density_ratio 
                    if p in curr:
                        curr[p] += value
                    else:
                        curr[p] = value
                    edges[r] = curr

        g = PyDiGraph()

        for reactant, products in edges.items():
            for product, count in products.items():

                if reactant_phase_set is not None and reactant not in reactant_phase_set:
                    continue

                if product_phase_set is not None and product in product_phase_set:
                    continue

                if product_chemsys and not is_in_chemsys(product, product_chemsys):
                    continue

                if reactant_chemsys and not is_in_chemsys(reactant, reactant_chemsys):
                    continue

                if phases_to_exclude is not None and (reactant in phases_to_exclude or product in phases_to_exclude):
                    continue

                if phases_to_include is not None and (reactant not in phases_to_include or product not in phases_to_include):
                    continue

                if required_phase is not None and product != required_phase and reactant != required_phase:
                    continue

                if exact_phase_set is not None and (product not in exact_phase_set or reactant not in exact_phase_set):
                    continue

                p_node = g.find_node_by_weight(reactant)
                if p_node is None:
                    p_node = g.add_node(reactant)

                c_node = g.find_node_by_weight(product)
                if c_node is None:
                    c_node = g.add_node(product)

                g.add_edge(p_node, c_node, count)

        self.graph: PyDiGraph = g
        return self
    
    def show(self, equal_weights=False, max_edge_size = 5):
        def node_fn(n):
            return {
                "label": n
            }
        
        max_edge_weight = max(self.graph.edges())

        def edgefn(e):
            if equal_weights:
                return {}
            else:
                return {
                    "penwidth": str(max(0.25, e / max_edge_weight * max_edge_size))
                }
        
        graph_attrs = {
            "ratio": "compress",
        }
        
        return graphviz_draw(self.graph, node_attr_fn=node_fn, edge_attr_fn=edgefn, graph_attr=graph_attrs)        
        