from typing import Union, List
import copy

def format_chem_sys(chem_sys: Union[List, str]):
    if type(chem_sys) is str:
        arr = chem_sys.split("-")
        arr.sort()
        return "-".join(arr)
    else:
        copy_arr = copy(chem_sys)
        copy_arr.sort()
        return copy_arr