from jobflow.core.maker import Maker
from jobflow.core.job import job

from pymatgen.ext.matproj import MPRester

from rxn_network.enumerators.minimize import MinimizeGibbsEnumerator
from rxn_network.entries.entry_set import GibbsEntrySet
from rxn_network.reactions.reaction_set import ReactionSet

from dataclasses import dataclass

from ..schemas.enumerated_rxns_schema import EnumeratedRxnsModel


@dataclass
class EnumerateRxnsMaker(Maker):
    """
    Downloads reaction given a chemical system to a filesystem store.
    """

    name: str = "Enumerate Rxns"

    @job
    def make(self,
             chem_sys: str,
             temp: int,
             stability_cutoff: float,
             open_el: str = None,
             chempot: float = None,
             mp_api_key: str = None) -> EnumeratedRxnsModel:

        if mp_api_key is None:
            with MPRester() as mpr:  # insert your Materials Project API key here if it's not stored in .pmgrc.yaml
                entries = mpr.get_entries_in_chemsys(chem_sys, inc_structure="final")
        else:
            with MPRester(api_key=mp_api_key) as mpr:  # insert your Materials Project API key here if it's not stored in .pmgrc.yaml
                entries = mpr.get_entries_in_chemsys(chem_sys, inc_structure="final")

        gibbs_entries = GibbsEntrySet.from_entries(entries, temp)
        entry_set = gibbs_entries.filter_by_stability(stability_cutoff)

        gibbs_enumerator = MinimizeGibbsEnumerator()
        gibbs_rxns = gibbs_enumerator.enumerate(entry_set)
        rxn_set = ReactionSet.from_rxns(
                    gibbs_rxns, entry_set, open_elem=open_el, chempot=chempot
                )

        result_model = EnumeratedRxnsModel.from_obj(
            rxn_set,
            chem_sys,
            temperature=temp,
            stability_cutoff=stability_cutoff,
            open_el=open_el,
            chem_pot=chempot,
        )

        return result_model