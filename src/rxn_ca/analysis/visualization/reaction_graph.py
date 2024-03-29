from rustworkx import PyDiGraph
from rustworkx.visualization import graphviz_draw
from ...computing.schemas import RxnCAResultDoc
from ...reactions import ReactionLibrary

from rxn_ca.utilities.helpers import is_in_chemsys
from collections import Counter


class ReactionGraph():
    
    def __init__(self, rxn_lib: ReactionLibrary, result_doc: RxnCAResultDoc):
        self._original_rxn_lib = rxn_lib
        self.phase_set = rxn_lib.phases
        
        r_chosen = []
        for res_choices in result_doc.metadata["rxn_choices"]:
            for choice in res_choices:
                r_chosen.append(choice)

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

    def build(self,
              reactant_phase_set = None,
              product_phase_set = None,
              product_chemsys = None,
              reactant_chemsys = None,
              phases_to_exclude = None,
              phases_to_include = None,
              required_phase = None,
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
                    stoich_ratio = rxn.product_stoich(p) / rxn.reactant_stoich(r)
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

                p_node = g.find_node_by_weight(reactant)
                if p_node is None:
                    p_node = g.add_node(reactant)

                c_node = g.find_node_by_weight(product)
                if c_node is None:
                    c_node = g.add_node(product)

                g.add_edge(p_node, c_node, count)

        self.graph: PyDiGraph = g
        return self
    
    def show(self, equal_weights=False):
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
                    "penwidth": str(max(0.25, e / max_edge_weight * 5))
                }
        
        
        return graphviz_draw(self.graph, node_attr_fn=node_fn, edge_attr_fn=edgefn)        
        