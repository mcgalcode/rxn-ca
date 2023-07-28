from rxn_ca.core import setup, run, get_arrhenius_rxns

import multiprocessing as mp
import json

print("================= RETRIEVING AND SCORING REACTIONS =================")
temp = 700
rxns = get_arrhenius_rxns("Y-Mn-C-Cl-O-Li", temp)
print()
print()
print()

print("================= SETTING UP SIMULATION =================")
phases = {
    "Mn2O3": 1,
    "YCl3": 2,
    "Li2O": 3
}

sim = setup(rxns, phases, 15, 30, dim=3)

PROCESSES = mp.cpu_count()

print(f'================= RUNNING SIMULATION w/ {PROCESSES} REALIZATIONS =================')


global mp_globals
mp_globals = {}


def get_result(start):
    return run(start, 30000)

with mp.get_context("fork").Pool(PROCESSES) as pool:
    results = pool.map(get_result, [sim for _ in range(PROCESSES)])

jserializable = [res.as_dict() for res in results]

with open("out", "w+") as f:
    f.write(json.dumps(jserializable))