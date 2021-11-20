import itertools
import json
import multiprocessing as mp
from time import perf_counter

import networkx as nx
import numpy as np
from drugable.utils.interleave import roundrobin
from matplotlib import pyplot as plt
from tqdm import tqdm

import chembl.chembl as chembl
import id_graph
import opentargets.opentargets as ot
import pharmgkb.pharmgkb as pharmgkb
import ttdb.ttdb as ttdb
from chemlinks import get_chembl_pc_linker
from DBObj import Drug, ObjContainer, Target, unify


def grouper_it(n, iterable):
    it = iter(iterable)
    while True:
        chunk_it = itertools.islice(it, n)
        try:
            first_el = next(chunk_it)
        except StopIteration:
            return
        yield itertools.chain((first_el,),chunk_it)

def drug_eq(tup):
    if not set(tup[0]._atomic_ids_not_null()).isdisjoint(tup[1]._atomic_ids_not_null()):
        return True 
    return False


def mp_check(group):

    with mp.Manager() as manager:
            with manager.Pool(processes=10) as pool:
                results = pool.map_async(func=drug_eq, iterable=list(group))
                output = [res for res in results.get()]
                return output

if __name__ == '__main__':
    pko = pharmgkb.pharmgkb()
    pk_drugs = pko.load_drugs()
    oto = ot.opentargets()
    ot_drugs = oto.load_drugs()
    ttdbo = ttdb.ttdb()
    ttdb_drugs = ttdbo.load_drugs()
    chemblo = chembl.chembl()
    chembl_drugs = chemblo.load_drugs()
    chembl_do = [Drug.from_dict(x) for x in chembl_drugs]
    ttdb_do = [Drug.from_dict(x) for x in ttdb_drugs]
    ot_do = [Drug.from_dict(x) for x in ot_drugs]
    pk_do = [Drug.from_dict(x) for x in pk_drugs]
    for o in chembl_do: 
        o._src = ['chembl']
    for o in ttdb_do:
        o._src = ['ttdb']
    for o in ot_do:
        o._src = ['opentargets']
    for o in pk_do:
        o._src = ['pharmgkb']
    all_do_i = list(roundrobin(chembl_do, ttdb_do, ot_do, pk_do))
    all_do = chembl_do + ttdb_do + ot_do + pk_do

    test_iter = itertools.combinations(all_do, 2)
    groups = grouper_it(1000, test_iter)
    all_done = []
    for group in tqdm(groups):
        res = mp_check(group)
        all_done.extend(zip(group, res))
