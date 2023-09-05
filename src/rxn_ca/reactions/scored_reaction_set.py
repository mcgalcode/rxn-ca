from typing import List

import json

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
            return cls.from_dict(json.loads(f.read()))

    @classmethod
    def from_dict(cls, rxn_set_dict):
        return cls(
            [ScoredReaction.from_dict(r) for r in rxn_set_dict["reactions"]],
        )

    def __init__(self, reactions: list[ScoredReaction], phase_set: SolidPhaseSet = None, identity_score = 1):
        """Initializes a SolidReactionSet object. Requires a list of possible reactions
        and the elements which should be considered available in the atmosphere of the
        simulation.

        Args:
            reactions (list[Reaction]):
        """
        self.reactant_map = {}
        self.reactions: List[ScoredReaction] = []
        self.rxn_map = {}
        # Replace strength of identity reaction with the depth of the hull its in

        for r in reactions:
            self.add_rxn(r)

        if phase_set is not None:
            for phase in phase_set.phases:
                if phase not in DEFAULT_GASES and phase is not SolidPhaseSet.FREE_SPACE:
                    self_rxn = ScoredReaction.self_reaction(phase, strength = identity_score)
                    existing = self.get_reactions([phase])
                    if len(existing) > 0 and not any(map(lambda rxn: rxn.is_identity, existing)):
                        self.add_rxn(self_rxn)
                    elif len(existing) == 0:
                        self.add_rxn(self_rxn)

    def rescore(self, scorer):
        rescored = [rxn.rescore(scorer) for rxn in self.reactions if not rxn.is_identity]
        skip_vols = bool(self.volumes)
        return ScoredReactionSet(rescored, skip_vols)

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

    def exclude_theoretical(self, phase_set: SolidPhaseSet):
        filtered = ScoredReactionSet([])
        theoretical_phases = phase_set.get_theoretical_phases()
        for r in self.reactions:
            containts_theoretical = False
            for p in r.all_phases:
                if p in theoretical_phases:
                    containts_theoretical = True

            if not containts_theoretical:
                filtered.add_rxn(r)
        
        return filtered

    def exclude_phases(self, phase_list: List[str]):
        filtered = ScoredReactionSet([])

        for r in self.reactions:
            contains_exclude = False
            for p in r.all_phases:
                if p in phase_list:
                    contains_exclude = True

            if not contains_exclude:
                filtered.add_rxn(r)
        
        return filtered

    def get_reactions(self, reactants: list[str]) -> ScoredReaction:
        """Given a list of string reaction names, returns a reaction that uses exactly those
        reactants as precursors.

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

    def search_reactants(self, reactants: list[str]) -> list[ScoredReaction]:
        """Returns all the reactions in this SolidReactionSet that produce all of the
        reactant phases specified.

        Args:
            reactants (list[str]): The reactants which matching reactions will produce.

        Returns:
            list[Reaction]: The matching reactions.
        """
        return [rxn for rxn in self.reactions if set(rxn.reactants).issuperset(reactants)]
    
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
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
        }