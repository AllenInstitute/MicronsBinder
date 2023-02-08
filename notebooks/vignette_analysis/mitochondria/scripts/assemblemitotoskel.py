"""
Assembling data associating mitochondria to skeleton nodes
"""
import sys
sys.path.append('.')
import numpy as np
import pandas as pd
from scipy.sparse import csgraph
import h5py

from lib import skel, compartment, u


def main(mito_filename, id_filename, output_filename,
         refine=False, no_diam=False):

    mitodf = pd.read_csv(mito_filename, index_col=0)
    segids = u.readids(id_filename)

    seg_dist_data = list()
    for (i, segid) in enumerate(segids):
        print(f"#{i+1} of {len(segids)}: segid {segid}")#, end="\r")
        seg_dist_data.append(assemble_dist_data(segid, mitodf,
                             refine=refine, no_diam=no_diam))

    print("")
    print("Writing results")

    seg_dist_data = consolidate_data(seg_dist_data)
    write_all_data(seg_dist_data, output_filename)


def assemble_dist_data(cellid, mitodf, refine=False, no_diam=False):
    segskel = skel.read_neuron_skel(cellid)
    nodelbls = skel.read_neuron_complbls(cellid)
    if refine:
        nodelbls = compartment.refine_labels(segskel, nodelbls)
    nodesomadist = compute_node_distance_to_soma(segskel, nodelbls)

    mito_to_node = skel.read_skel_assoc(cellid)

    t_mitoids = list()
    t_mitovols = list()
    t_nodeids = list()
    t_counts = list()
    t_nodelbls = list()
    t_nodedists = list()
    t_cellids = list()

    for (mitoid, nodeids) in mito_to_node.items():
        mitovol = int(mitodf.mito_vx.loc[mitoid])
        for (nodeid, count) in nodeids:
            t_cellids.append(cellid)
            t_mitoids.append(mitoid)
            t_mitovols.append(mitovol)
            t_nodeids.append(nodeid)
            t_counts.append(count)
            t_nodelbls.append(nodelbls[nodeid])
            t_nodedists.append(nodesomadist[nodeid])

    data_dict = dict(mitoid=t_mitoids, mitovol=t_mitovols,
                     nodeid=t_nodeids, nodelbl=t_nodelbls,
                     nodedist=t_nodedists, cellid=t_cellids,
                     count=t_counts)

    if not no_diam:
        data_dict["diams"] = t_diams

    return data_dict


def compute_node_distance_to_soma(segskel, nodelbls):
    soma_ids = np.nonzero(nodelbls == 0)

    dists = csgraph.dijkstra(
                segskel.csgraph, directed=False, indices=soma_ids)

    return np.min(np.min(dists, axis=0), axis=0)


def consolidate_data(seg_dist_data):
    assert len(seg_dist_data) > 0, "no data provided"
    all_keys = [d.keys() for d in seg_dist_data]
    assert [all_keys[0] == ks for ks in all_keys], "mismatched keys"

    ks = all_keys[0]
    single_dict = dict()
    for k in ks:
        single_dict[k] = np.concatenate([d[k] for d in seg_dist_data])

    return single_dict


def write_all_data(seg_dist_data, output_filename):
    with h5py.File(output_filename, 'w') as f:
        for (k, v) in seg_dist_data.items():
            f.create_dataset(k, data=v)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("mito_filename")
    parser.add_argument("id_filename")
    parser.add_argument("output_filename")
    parser.add_argument("--refine", action="store_true")
    parser.add_argument("--no_diam", action="store_true")

    args = parser.parse_args()

    main(**vars(args))
