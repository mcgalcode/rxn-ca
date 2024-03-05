from rxn_network.enumerators.minimize import MinimizeGibbsEnumerator
from rxn_network.entries.entry_set import GibbsEntrySet
from rxn_network.jobs.core import GetEntrySetMaker, ReactionEnumerationMaker

from ..computing.schemas.enumerated_rxns_schema import EnumeratedRxnsModel

def enumerate_rxns(entries: GibbsEntrySet,
                   open_el: str = None,
                   chempot: float = None) -> EnumeratedRxnsModel:


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


    print(f'Using {len(entries)} entries')

    enumerator = MinimizeGibbsEnumerator()
    enumeration_maker = ReactionEnumerationMaker()
    enumerators = [enumerator]
    rxns = enumeration_maker.make.original(enumeration_maker, enumerators, entries)

    result_model = EnumeratedRxnsModel.from_obj(
        rxns.rxns,
        chem_sys,
        stability_cutoff=stability_cutoff,
        open_el=open_el,
        chem_pot=chempot,
    )

    return result_model