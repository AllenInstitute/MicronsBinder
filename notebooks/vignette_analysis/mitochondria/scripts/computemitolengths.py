"""Computing the path length covered by each mitochondrion"""
import sys
import argparse
from collections import defaultdict

from tqdm import tqdm
import numpy as np
import pandas as pd

sys.path.append(".")
from lib import compartment, skel, u


def main(mitodffilename, idfilename, mitotoskelfilename, outputfilename):

    mitodf = pd.read_csv(mitodffilename, index_col=0)
    mitotoskeldf = compartment.readmitotoskel(mitotoskelfilename)
    ids = u.readids(idfilename)

    lengths = computemitolens(mitodf, mitotoskeldf, ids)

    lengths.to_csv(outputfilename)


def computemitolens(mitodf, mitotoskeldf, ids):
    subdfs = list()
    for cellid in tqdm(ids):
        subdfs.append(computecellmitolengths(mitodf, mitotoskeldf, cellid))

    return pd.concat(subdfs)


def computecellmitolengths(mitodf, mitotoskeldf, cellid):
    cellmitodf = mitodf[mitodf.cellid == cellid].copy()
    cellmitotoskel = mitotoskeldf[mitotoskeldf.cellid == cellid]

    cellskel = skel.read_neuron_skel(cellid)

    assocmitolookup = dict()
    for (nodeid, mitoid) in zip(cellmitotoskel.nodeid, cellmitotoskel.mitoid):
        if nodeid not in assocmitolookup:
            assocmitolookup[nodeid] = {mitoid}
        else:
            assocmitolookup[nodeid].add(mitoid)

    edgelengths = np.linalg.norm(cellskel.vertices[cellskel.edges[:, 1]]
                                 - cellskel.vertices[cellskel.edges[:, 0]],
                                 axis=1)

    emptyset = set()
    mitolengths = defaultdict(float)
    for ((n1, n2), edgelength) in zip(cellskel.edges, edgelengths):
        set1 = assocmitolookup.get(n1, emptyset)
        set2 = assocmitolookup.get(n2, emptyset)
        spanningmitos = set1.intersection(set2)
        for m in spanningmitos:
            mitolengths[m] += edgelength

    cellmitodf["pathlength"] = np.array([mitolengths[m]
                                         for m in cellmitodf.index])

    return cellmitodf


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("mitodffilename")
    parser.add_argument("idfilename")
    parser.add_argument("mitotoskelfilename")
    parser.add_argument("outputfilename")

    main(**vars(parser.parse_args()))
