from jobflow.core.maker import Maker
from jobflow.core.job import job

from rxn_network.enumerators.minimize import MinimizeGibbsEnumerator
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
             stability_cutoff: float = 0.1,
             open_el: str = None,
             chempot: float = None,
             mp_api_key: str = None) -> EnumeratedRxnsModel:

        entry_set_maker = GetEntrySetMaker(
            temperature=temp,
            e_above_hull=stability_cutoff,
            MP_API_KEY=mp_api_key
        )


        # The difference between the open enumeration and the basic enumeration is just that another element
        # is included - you can actually reset the chemical potential _after_ enumeration
        #
        # To set a NEW CHEMICAL POTENTIAL - reinitialize a new ReactionSet object with the same underlying
        # data but with different chempot and open_el values
        #
        # To set a new temperature - iterate over every reaction in the set and call get_new_temperature
        # we will ned up with many more reactions than we had previously, so it is probably good to
        # filter out reactions with positive delta Gs, say more than 0.1eV per atom
        #
        # potential way to use the grand potential enumerator is to run it at a variety of chemical potentials,
        # then add all the resulting reaction sets to gether, and set the chempot to the desired value for all of them
        # can calculate oxygen chempots with ktlnPO2 / 2

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