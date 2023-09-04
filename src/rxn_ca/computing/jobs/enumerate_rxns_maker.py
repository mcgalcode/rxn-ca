from jobflow.core.maker import Maker
from jobflow.core.job import job

from rxn_network.enumerators.minimize import MinimizeGibbsEnumerator
from rxn_network.jobs.core import GetEntrySetMaker, ReactionEnumerationMaker
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
             temp: int = 300,
             stability_cutoff: float = 0.1,
             open_el: str = None,
             chempot: float = None,
             mp_api_key: str = None,
             formulas_to_include: dict = []) -> EnumeratedRxnsModel:

        # First we enumerate entries
        entry_set_maker = GetEntrySetMaker(
            temperature=temp,
            e_above_hull=stability_cutoff,
            MP_API_KEY=mp_api_key,
            formulas_to_include=formulas_to_include
        )

        eset: GibbsEntrySet = entry_set_maker.make.original(entry_set_maker, chem_sys)

        print(f'Using {len(eset.entries)} entries')
        gibbs_enumerator = MinimizeGibbsEnumerator()
        enumerators = [gibbs_enumerator]
            
        enumeration_maker = ReactionEnumerationMaker()
        entries = eset.entries
        
        rxns: ReactionSet = enumeration_maker.make.original(enumeration_maker, enumerators, entries)
        result_model = EnumeratedRxnsModel.from_obj(
            rxns.rxns,
            chem_sys,
            temperature=temp,
            stability_cutoff=stability_cutoff,
            open_el=open_el,
            chem_pot=chempot,
        )

        return result_model