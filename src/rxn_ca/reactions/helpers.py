from mp_api.client import MPRester

def enumerate_rxns(self,
            chem_sys: str,
            temp: int,
            stability_cutoff: float,
            open_el: str = None,
            chempot: float = None,
            mp_api_key: str = None):

    entry_set_maker = GetEntrySetMaker(
        temperature=1000,
        e_above_hull=0.1,
        MP_API_KEY=mp_api_key
    )

    with MPRester() as mpr:
        mpr.get_entries_in_chemsys(chem_sys)

    eset = entry_set_maker.make.original(entry_set_maker, chem_sys)
    enumerator = MinimizeGibbsEnumerator()
    enumeration_maker = ReactionEnumerationMaker()
    enumerators = [enumerator]
    entries = eset.entries
    rxns = enumeration_maker.make.original(enumeration_maker, enumerators, entries)
        