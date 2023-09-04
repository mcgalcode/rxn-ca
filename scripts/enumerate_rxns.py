from rxn_ca import AutomatonStore, enumerate_rxns

chem_syses = [
    "Ba-Cl-Na-Ti-O",
    "Ba-S-Na-Ti-O",
    "Ba-Ti-O",
    "Ba-Ti-O-H",
    "Ba-Ti-O-S",
    "Fe-Si-S",
    "Li-Mn-O-H",
    "Li-Mn-Ti-F-C-O",
    "Y-Ba-Cu-O",
    "Y-Mn-C-Cl-O-Li",
    "Cl-Li-Mn-O-Y",
    "Yb-Ru-Sn",
    "Y-Mn-O",
    "Bi-Fe-O",
    "Co-Cl-H-O",
    "Ca-S-O-H",
    "Cu-O"
]

special_cutoffs = {
    "Li-Mn-Ti-F-C-O": [0.01],
    "Y-Mn-C-Cl-O-Li": [0.01],
    "Y-Mn-Cl-O-Li": [0.01]    
}

full_extras = list(range(400,1500,100))

extra_temps = {
    "Ba-Cl-Na-Ti-O": full_extras,
    "Ba-S-Na-Ti-O": full_extras,
    "Ba-Ti-O": full_extras,
    "Ba-Ti-O-H": full_extras,
    "Ba-Ti-O-S": full_extras,
    "Fe-Si-S": full_extras,
    "Li-Mn-O-H": full_extras,
    "Li-Mn-Ti-F-C-O": full_extras,
    "Y-Ba-Cu-O": full_extras,
    "Y-Mn-C-Cl-O-Li": full_extras,
    "Cl-Li-Mn-O-Y": full_extras,
    "Yb-Ru-Sn": full_extras,
    "Bi-Fe-O": full_extras,
    "Co-Cl-H-O": full_extras,
    "Ca-S-O-H": full_extras,
    "Cu-O": full_extras
}

store = AutomatonStore()

store.list_available_sets()

default_cutoffs = [0.01]
extended_cutoffs = [0.01, 0.1, 0.2, 0.5]

for sys in chem_syses:
    if sys in special_cutoffs:
        cutoffs = special_cutoffs.get(sys)
    else:
        cutoffs = default_cutoffs
        
    for cutoff in cutoffs:
        print(f"Enumerating rxns for {sys} at {cutoff} ")          
        enumerate_rxns(sys,
                       stability_cutoff=cutoff,
                       base_temp=300,
                       other_temps=extra_temps.get(sys))