from typing import Dict

def normalize_dict(d: Dict):
    total = sum(d.values())
    return { k: v / total for k, v in d.items() }

def add_values_to_dict_by_addition(target, new_vals):
    for k, v in new_vals.items():
        if k in target:
            target[k] = target[k] + v
        else:
            target[k] = v