from typing import List, Dict

import json

from monty.json import MontyDecoder, MontyEncoder, MSONable

from .scored_reaction import ScoredReaction
from ..phases.solid_phase_set import SolidPhaseSet
from ..phases.gasses import DEFAULT_GASES

from pymatgen.core.composition import Composition

import matplotlib.pyplot as plt

class ScoredReactionSet():
    """A set of ScoredReactions that capture the events that can occur during a simulation. Typically
    includes every reaction possible in the chemical system defined by the precursors and open
    elements
    """

    @classmethod
    def from_file(cls, fpath):
        with open(fpath, 'r+') as f:
            return cls.from_dict(json.load(f,))

    @classmethod
    def from_dict(cls, rxn_set_dict: Dict):
        return cls(
            [ScoredReaction.from_dict(r) for r in rxn_set_dict["reactions"]],
            SolidPhaseSet.from_dict(rxn_set_dict["phase_set"])
        )

    IDENTITY = "IDENTITY"

    def __init__(self, reactions: list[ScoredReaction], phase_set: SolidPhaseSet):
        """Initializes a SolidReactionSet object. Requires a list of possible reactions
        and the elements which should be considered available in the atmosphere of the
        simulation.

        Args:
            reactions (list[Reaction]):
        """
        if phase_set is None:
            raise ValueError("phase_set is required when instantiating a ScoredReactionSet")
        self.phases = phase_set
        self.reactant_map = {}
        self.reactions: List[ScoredReaction] = []
        self.rxn_map = {}
        
        # Replace strength of identity reaction with the depth of the hull its in
        for r in reactions:
            self.add_rxn(r)  

    def rescore(self, scorer):
        rescored = [rxn.rescore(scorer) for rxn in self.reactions]
        return ScoredReactionSet(rescored, self.phases)

    def add_rxn(self, rxn: ScoredReaction) -> None:
        reactant_set = frozenset(rxn.reactants)
        if self.reactant_map.get(reactant_set) is None:
            self.reactant_map[reactant_set] = [rxn]
        else:
            self.reactant_map[reactant_set].append(rxn)
            self.reactant_map[reactant_set] = sorted(self.reactant_map[reactant_set], key = lambda rxn: rxn.competitiveness, reverse = True)
    
        self.rxn_map[str(rxn)] = rxn
        self.reactions.append(rxn)

    def exclude_pure_els(self):
        filtered = ScoredReactionSet([])
        for r in self.reactions:
            contains_element = False
            for p in r.all_phases:
                if len(Composition(p).elements) == 1:
                    contains_element = True
            if not contains_element:
                filtered.add_rxn(r)
        
        return filtered

    def exclude_theoretical(self, ensure_phases: List[str] = []):
        filtered = ScoredReactionSet([], self.phases)
        theoretical_phases = self.phases.get_theoretical_phases()
        for r in self.reactions:
            contains_theoretical = False
            for p in r.all_phases:
                if p in theoretical_phases and p not in ensure_phases:
                    contains_theoretical = True

            if not contains_theoretical:
                filtered.add_rxn(r)
        
        return filtered
    
    def exclude_metastable(self, metastability_cutoff: float, ensure_phases: List[str] = []):
        filtered = ScoredReactionSet([], self.phases)
        for r in self.reactions:
            contains_metastable = False
            for p in r.all_phases:
                if self.phases.get_e_above_hull(p) > metastability_cutoff and p not in ensure_phases:
                    contains_metastable = True

            if not contains_metastable:
                filtered.add_rxn(r)
        
        return filtered        

    def exclude_phases(self, phase_list: List[str]):
        filtered = ScoredReactionSet([], self.phases)

        for r in self.reactions:
            contains_exclude = False
            for p in r.all_phases:
                if p in phase_list:
                    contains_exclude = True

            if not contains_exclude:
                filtered.add_rxn(r)
        
        return filtered
    
    def limit_phases(self, phase_list: List[str]):
        comps = [Composition(p) for p in phase_list]
        filtered_rxn_set = ScoredReactionSet([], self.phases)
        for r in self.reactions:
            if all([Composition(p) in comps for p in r.all_phases]):
                filtered_rxn_set.add_rxn(r)
        return filtered_rxn_set

    def get_reactions(self, reactants: list[str]) -> List[ScoredReaction]:
        """Given a list of reactants, returns the list of reactions which
        consume exactly that set of precursors

        Args:
            reactants (list[str]): The list of reactants to match with

        Returns:
            Reaction: The matching reaction, if it exists, otherwise None.
        """
        return self.reactant_map.get(frozenset(reactants), [])

    def get_rxn_by_str(self, rxn_str: str) -> ScoredReaction:
        """Retrieves a reaction from this set by it's serialized string form

        Args:
            rxn_str (str):

        Returns:
            Reaction:
        """
        return self.rxn_map.get(rxn_str)

    def search_products(self, products: list[str]) -> list[ScoredReaction]:
        """Returns all the reactions in this SolidReactionSet that produce all of the
        product phases specified.

        Args:
            products (list[str]): The products which matching reactions will produce.

        Returns:
            list[Reaction]: The matching reactions.
        """
        return [rxn for rxn in self.reactions if set(rxn.products).issuperset(products)]

    def search_all(self, products: list[str], reactants: list[str]) -> list[ScoredReaction]:
        return [rxn for rxn in self.reactions if set(rxn.products).issuperset(products) and set(rxn.reactants).issuperset(reactants)]

    def search_reactants(self, reactants: list[str], exact = False) -> list[ScoredReaction]:
        """Returns all the reactions in this SolidReactionSet that produce all of the
        reactant phases specified.

        Args:
            reactants (list[str]): The reactants which matching reactions will produce.

        Returns:
            list[Reaction]: The matching reactions.
        """
        if not exact:
            return [rxn for rxn in self.reactions if set(rxn.reactants).issuperset(reactants)]
        else:
            return [rxn for rxn in self.reactions if set(rxn.reactants) == set(reactants)]
    
    def plot_energies(self, bins=300):
        es = [r.energy_per_atom for r in self.reactions]
        plt.title("Reaction Energies")
        plt.xlabel("Delta G (eV/atom)")
        plt.ylabel("Count")
        plt.hist(es, bins=bins)

    def plot_scores(self, bins=300):
        scores = [r.competitiveness for r in self.reactions]
        plt.title("Reaction Scores")
        plt.xlabel("Score")
        plt.ylabel("Count")        
        plt.hist(scores, bins=bins)

    def as_dict(self):
        return {
            "reactions": [r.as_dict() for r in self.reactions],
            "phase_set": self.phases.as_dict(),
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
        }
    
    def __len__(self):
        return len(self.rxn_map)