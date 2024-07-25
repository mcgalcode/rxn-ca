from __future__ import annotations
import typing
from numbers import Number

from rxn_network.reactions.basic import BasicReaction
from ..phases.gasses import DEFAULT_GASES
from ..phases import SolidPhaseSet

def stoich_map_to_str(stoich_map: typing.Dict[str, Number]) -> str:
    """Generates a string that encapsulates a stoichiometry map. For example
    the map { "Na": 1, "Cl": 1 } will become 1Na + 1Cl. Useful for serializing
    reactions as strings.

    Args:
        stoich_map (dict): The stoichiometry map as exemplified in the description.

    Returns:
        string:
    """
    result = ""

    for phase, stoich in stoich_map.items():
        result = result + f"{stoich}{phase}+"

    result = result[:-1]
    return result


def phases_to_str(phases: list[str]) -> str:
    """Generates a string of phases concatenated with a plus sign.

    Args:
        phases (list[str]): _description_

    Returns:
        _type_: _description_
    """
    phases = sorted(list(set(phases)))
    return "+".join(phases)

class ScoredReaction:

    NO_RXN = "NO_RXN"

    @classmethod
    def from_dict(cls, rxn_dict):
        return cls(
            rxn_dict["reactants"],
            rxn_dict["products"],
            rxn_dict["competitiveness"],
            energy_per_atom = rxn_dict.get("energy_per_atom")
        )

    @classmethod
    def from_rxn_network(cls, score, original_rxn: BasicReaction, volumes: typing.Dict) -> ScoredReaction:
        react_dict = { comp.reduced_formula: round(-coeff * volumes.get(comp.reduced_formula), 2) for comp, coeff in original_rxn.reactant_coeffs.items() }
        product_dict = { comp.reduced_formula: round(coeff * volumes.get(comp.reduced_formula), 2) for comp, coeff in original_rxn.product_coeffs.items() }
        return ScoredReaction(react_dict, product_dict, score, energy_per_atom=original_rxn.energy_per_atom)

    def __init__(self, reactants, products, competitiveness, energy_per_atom = None):
        """Instantiate a reaction object by providing stoichiometry maps describing the
        reactant and product stoichiometry, and the relative competitiveness of this
        reaction.

        NOTE: Stoichiometry should be provided in terms of _volume_ for the purpose of this model.

        Args:
            reactants (typing.Dict[str, Number]): A map representing the stoichiometry of the
            reactants, e.g. { "Na": 1, "Cl": 1 }
            products (typing.Dict[str, Number]): A map representing the stoichiometry of the products.
            competitiveness (Number): A competitiveness score for the reaction.
        """
        self._reactants: typing.Dict[str, Number] = reactants
        self._products: typing.Dict[str, Number] = products

        self.reactants = frozenset(self._reactants.keys())
        self.products = frozenset(self._products.keys())

        self.solid_reactants = frozenset([ r for r in self._reactants.keys() if r not in DEFAULT_GASES])
        self.solid_products = frozenset([ r for r in self._products.keys() if r not in DEFAULT_GASES])

        self.is_identity = self.reactants == self.products

        self.total_reactant_stoich = sum(reactants.values())
        self.total_product_stoich = sum(products.values())

        self.total_solid_reactant_stoich = sum([reactants[r] for r in self.solid_reactants])
        self.total_solid_product_stoich = sum([products[r] for r in self.solid_products])

        self.product_reactant_stoich_ratio = self.total_product_stoich / self.total_reactant_stoich
        self.solid_product_reactant_stoich_ratio = self.total_solid_product_stoich / self.total_solid_reactant_stoich

        self.competitiveness: Number = competitiveness
        self._as_str = f"{stoich_map_to_str(self._reactants)}->{stoich_map_to_str(self._products)}"
        self.energy_per_atom = energy_per_atom

    def rescore(self, scorer) -> ScoredReaction:
        new_score = scorer.score(self)
        return ScoredReaction(self._reactants, self._products, new_score)

    def can_proceed_with(self, reactants: list[str]) -> bool:
        """Helper method that, given a list of reactants, returns true if it is the same
        as the list of reactants for this reaction. Note that this is an exact match.

        Args:
            reactants (list[str]): A list of reactant phase names

        Returns:
            bool: True if the reactants match, otherwise False
        """
        return set(reactants) == self.reactants

    def reactant_str(self) -> str:
        """Returns a string representing the reactant side of this reaction. For instance,
        1Na+1Cl

        Returns:
            str:
        """
        return phases_to_str(self.reactants)

    @property
    def all_phases(self):
        return list(set(list(self.reactants) + list(self.products)))

    def stoich_ratio(self, phase1, phase2) -> Number:
        all_phases = {**self._reactants, **self._products}
        return all_phases[phase1] / all_phases[phase2]

    def product_stoich(self, phase: str) -> Number:
        """Returns the stoichiometry in this reaction for the desired product phase.

        Args:
            phase (str): The phase whose stoichiometry is desired

        Returns:
            Number:
        """
        return self._products[phase]

    def reactant_stoich(self, phase: str) -> Number:
        """Returns the stoichiometry in this reaction for the desired product phase.

        Args:
            phase (str): The phase whose stoichiometry is desired

        Returns:
            Number:
        """
        return self._reactants[phase]

    def reactant_stoich_fraction(self, phase: str) -> Number:
        """Returns the stoichiometry in this reaction for the desired reactant phase.

        Args:
            phase (str): The phase whose stoichiometry is desired

        Returns:
            Number:
        """
        try:
            return self._reactants[phase] / self.total_reactant_stoich
        except:
            print(phase, str(self))

    def product_stoich_fraction(self, phase: str) -> Number:
        """Returns the stoichiometry in this reaction for the desired product phase.

        Args:
            phase (str): The phase whose stoichiometry is desired

        Returns:
            Number:
        """
        try:
            return self._products[phase] / self.total_product_stoich
        except:
            print(phase, str(self))

    def solid_reactant_stoich_fraction(self, phase: str) -> Number:
        """Returns the stoichiometry in this reaction for the desired product phase.

        Args:
            phase (str): The phase whose stoichiometry is desired

        Returns:
            Number:
        """
        try:
            return self._reactants[phase] / self.total_solid_reactant_stoich
        except:
            print(phase, str(self))
    
    def convert_reactant_amt_to_product_amt(self, reactant: str, reactant_vol: float, product: str) -> float:
        """Converts a volume of reactant to a volume of product by ratio

        Args:
            reactant (str): The formula of the reactant
            reactant_vol (float): The volume to convert,
            product (str): The product to produce

        Returns:
            float: The volume of product produced
        """
        ratio = self._products[product] / self._reactants[reactant]
        return reactant_vol * self.solid_product_reactant_stoich_ratio
    
    def convert_to_moles(self, phase_set: SolidPhaseSet):
        reactants_moles = phase_set.vol_amts_to_moles(self._reactants, should_round=3)
        products = phase_set.vol_amts_to_moles(self._products, should_round=3)

        return ScoredReaction(reactants_moles, products, competitiveness=self.competitiveness, energy_per_atom=self.energy_per_atom)

    def any_reactants(self, phases):
        return len(self.reactants.intersection(phases)) > 0

    def __str__(self):
        return f"{self._as_str}, Score: {self.competitiveness}, E/atom: {self.energy_per_atom}"

    def as_dict(self):
        return {
            "reactants": self._reactants,
            "products": self._products,
            "competitiveness": self.competitiveness,
            "energy_per_atom": self.energy_per_atom,
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
        }
