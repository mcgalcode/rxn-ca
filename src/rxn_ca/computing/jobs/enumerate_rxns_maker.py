from jobflow.core.maker import Maker
from jobflow.core.job import job

from pymatgen.ext.matproj import MPRester

from rxn_network.enumerators.minimize import MinimizeGibbsEnumerator
from rxn_network.entries.entry_set import GibbsEntrySet
from rxn_network.reactions.reaction_set import ReactionSet
from rxn_network.jobs.core import GetEntrySetMaker, ReactionEnumerationMaker

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

        entry_set_maker = GetEntrySetMaker(
            temperature=1000,
            e_above_hull=0.1,
            MP_API_KEY=mp_api_key
        )

        eset = entry_set_maker.make.original(entry_set_maker, chem_sys)
        enumerator = MinimizeGibbsEnumerator()
        enumeration_maker = ReactionEnumerationMaker()
        enumerators = [enumerator]
        entries = eset.entries
        rxns = enumeration_maker.make.original(enumeration_maker, enumerators, entries)
        

        result_model = EnumeratedRxnsModel.from_obj(
            rxns.rxns,
            chem_sys,
            temperature=temp,
            stability_cutoff=stability_cutoff,
            open_el=open_el,
            chem_pot=chempot,
        )

        return result_model