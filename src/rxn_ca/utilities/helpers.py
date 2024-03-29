from typing import Dict, Union, List

from copy import copy
from pymatgen.core.composition import Composition

def normalize_dict(d: Dict):
    total = sum(d.values())
    return { k: v / total for k, v in d.items() }

def add_values_to_dict_by_addition(target, new_vals):
    for k, v in new_vals.items():
        if k in target:
            target[k] = target[k] + v
        else:
            target[k] = v

def format_chem_sys(chem_sys: Union[List, str]):
    if type(chem_sys) is str:
        arr = chem_sys.split("-")
        arr.sort()
        return "-".join(arr)
    else:
        copy_arr = copy(chem_sys)
        copy_arr.sort()
        return copy_arr
    
def is_in_chemsys(phase, chemsys: Union[List, str]):
    els = [str(el) for el in Composition(phase).elements]
    if type(chemsys) is str:
        chemsys = chemsys.split("-")
    return set(chemsys).issuperset(els)