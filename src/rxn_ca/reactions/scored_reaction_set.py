from typing import Dict
from pymatgen.ext.matproj import MPRester

import json

from .scored_reaction import ScoredReaction, phases_to_str

def get_phase_vols(phases):

    volumes = {}

    with MPRester() as mpr:
        res = mpr.summary.search(
            formula = phases,
            energy_above_hull = (None, 0.1),
            fields=["structure", "composition", "task_id"]
        )
        for item in res:
            struct = item.structure
            comp = item.composition
            volumes[comp.reduced_formula] = struct.volume / comp.get_reduced_composition_and_factor()[1]

    return volumes

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
            volumes = rxn_set_dict.get('volumes'),
        )

    def __init__(self, reactions: list[ScoredReaction], skip_vols: bool = False, volumes: Dict[str, float] = None):
        """Initializes a SolidReactionSet object. Requires a list of possible reactions
        and the elements which should be considered available in the atmosphere of the
        simulation.

        Args:
            reactions (list[Reaction]):
        """
        self.reactant_map = {}
        self.reactions = []
        self.rxn_map = {}
        self.phases = []
        # Replace strength of identity reaction with the depth of the hull its in

        for r in reactions:
            self.add_rxn(r)

        # for phase in self.phases:
        #     self_rxn = ScoredReaction.self_reaction(phase, strength = 0.1)
        #     existing = self.get_reaction([phase])
        #     if existing is not None and not existing.is_identity:
        #         self.add_rxn(self_rxn)
        #     elif existing is None:
        #         self.add_rxn(self_rxn)

        if volumes is not None:
            self.volumes = volumes
        elif not skip_vols:
            self.volumes = get_phase_vols(self.phases)

    def rescore(self, scorer):
        rescored = [rxn.rescore(scorer) for rxn in self.reactions if not rxn.is_identity]
        skip_vols = bool(self.volumes)
        return ScoredReactionSet(rescored, skip_vols)

    def add_rxn(self, rxn: ScoredReaction) -> None:
        self.reactant_map[rxn.reactant_str()] = rxn
        self.rxn_map[str(rxn)] = rxn
        self.reactions.append(rxn)
        for phase in rxn.all_phases:
            self.phases = list(set(self.phases + [phase]))

    def get_reaction(self, reactants: list[str]) -> ScoredReaction:
        """Given a list of string reaction names, returns a reaction that uses exactly those
        reactants as precursors.

        Args:
            reactants (list[str]): The list of reactants to match with

        Returns:
            Reaction: The matching reaction, if it exists, otherwise None.
        """
        return self.reactant_map.get(phases_to_str(reactants), None)

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

    def as_dict(self):
        return {
            "reactions": [r.as_dict() for r in self.reactions],
            "volumes": self.volumes,
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
        }