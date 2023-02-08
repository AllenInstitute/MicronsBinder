"""
Associating mitochondria to a set of skeleton nodes

More precisely, finding the set of skeleton nodes that is the closest to
some mesh vertex of each mitochondrion. Labels each such mitochondrion
with the set of skeleton node indices that is paired with it in this way.
"""
import sys
import time
import operator
import argparse
import itertools
from collections import Counter

import numpy as np
import pandas as pd
from scipy.spatial import KDTree
import h5py
from multiwrapper import multiprocessing_utils as mu

sys.path.append(".")
from lib import skel, mesh, u


def main(mito_filename, id_filename,
         output_prefix, parallel=1, verbose=True):

    mitodf = pd.read_csv(mito_filename, index_col=0)

    segids = u.readids(id_filename)

    multiargs = make_batches(segids, mitodf, output_prefix, parallel)

    if parallel == 1:
        find_assoc_mitos_batch(multiargs[0])

    else:
        mu.multiprocess_func(
            find_assoc_mitos_batch,
            multiargs, n_threads=parallel,
            verbose=verbose)


def make_batches(segids, mitodf, output_prefix, parallel=1):
    n_jobs = parallel * 2 if parallel > 1 else 1
    segid_batches = np.array_split(segids, n_jobs)

    mitoid_batches = list()
    for batch in segid_batches:
        mitoid_batch = [list(mitodf.index[mitodf.cellid == segid])
                        for segid in batch]
        mitoid_batches.append(mitoid_batch)

    assert len(segid_batches) == len(mitoid_batches)

    multiargs = [[segid_batch, mitoid_batch, output_prefix]
                 for segid_batch, mitoid_batch
                 in zip(segid_batches, mitoid_batches)]

    return multiargs


def find_assoc_mitos_batch(args):
    segids, mitoid_sets, output_prefix = args

    for (segid, mitoids) in zip(segids, mitoid_sets):
        assocs = find_assoc_mitos(segid, mitoids)
        write_assocs(assocs, segid, output_prefix)


def find_assoc_mitos(segid, mitoids, parallel=1):
    """Find the skeleton nodes associated with each mitochondrion"""
    start = time.time()
    print(f"Processing seg {segid}")

    if len(mitoids) == 0:
        print(f"No overlapping mitos for seg {segid}")
        return dict()

    skel_kdt = read_skel_kdtree(segid)
    mesh.download_meshes(mitoids, parallel=parallel)

    mito_to_nodes = dict()
    for mitoid in mitoids:
        mitomesh = mesh.read_mito_mesh(mitoid)
        add_associations(skel_kdt, mitomesh, mitoid, mito_to_nodes)

    print(f"Processing seg {segid} finished in {time.time() - start}s")

    return mito_to_nodes


def read_skel_kdtree(segid):
    segskel = skel.read_neuron_skel(segid)
    print(len(segskel.vertices))

    return KDTree(segskel.vertices)


def add_associations(kdtree, mesh, meshid, records):
    records[meshid] = Counter(kdtree.query(v)[1] for v in mesh.vertices)

    return records


def write_assocs(assocs, segid, output_prefix):
    items = assocs.items()

    def getcounteritems(assoc_item):
        # extracting the items from each counter
        return assoc_item[1].items()

    data = np.array(
             list(itertools.chain.from_iterable(
                      map(getcounteritems, items))),
             dtype=np.uint32)

    mitoids = np.array(
                  list(map(operator.itemgetter(0), items)),
                  dtype=np.uint32)

    ptr = 0
    ptrs = list()
    for (mitoid, v) in items:
        ptrs.append((ptr, ptr+len(v)))
        ptr += len(v)

    output_filename = f"{output_prefix}_{segid}.h5"
    with h5py.File(output_filename, 'w') as f:
        f.create_dataset(f"{segid}/data", data=data)
        f.create_dataset(f"{segid}/mitoids", data=mitoids)
        f.create_dataset(f"{segid}/ptrs", data=np.array(ptrs))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("mito_filename")
    parser.add_argument("id_filename")
    parser.add_argument("output_prefix")
    parser.add_argument("--parallel", type=int, default=1)

    args = parser.parse_args()

    main(**vars(args))
