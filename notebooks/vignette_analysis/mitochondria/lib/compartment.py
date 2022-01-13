"""Labeling mitochondria and neurons by compartment (axonal, dendritic, somatic)"""
from collections import Counter

import numpy as np
import pandas as pd
from scipy.spatial import KDTree
from scipy.sparse import csgraph
import h5py

from . import skel
from . import u


ENGLISHLABELS = {
    0: "Somatic", 1: "Axonal", 2: "Basal",
    3: "Apical", 4: "Unknown", 5: "Unknown dendritic"
}
MITOTOSKEL_FILENAME = "data/mitotoskel_v185.h5"


def toskelfields(filename=MITOTOSKEL_FILENAME):
    with h5py.File(filename, 'r') as f:
        return list(f.keys())


def readtoskeldata(dset_name, filename=MITOTOSKEL_FILENAME):
    with h5py.File(filename, 'r') as f:
        return f[dset_name][()]


def readmitotoskel(filename=MITOTOSKEL_FILENAME):
    fields = dict()
    for field in toskelfields(filename=filename):
        fields[field] = readtoskeldata(field, filename=filename)

    return pd.DataFrame(fields)


def refine_labels(segskel, complbls, n_edges=1):
    """
    Restricting the soma label to a smaller distance threshold from the root node
    """
    assert n_edges >= 0, "need a nonnegative edge threshold"

    # Find all nodes within n_edges of the root
    soma_edge_dist = skel.path_length_to_root(segskel, segskel.root, binary=True)
    soma_nodes = np.nonzero(soma_edge_dist <= n_edges)[0]

    newlbls = np.zeros((len(complbls),), dtype=np.uint8)
    for path in skel.full_cover_paths(segskel):
        c = Counter(complbls[path])
        newlbls[path] = c.most_common()[0][0]
        newlbls[path[np.isin(path, soma_nodes)]] = 0

    return newlbls


def majority_vote_label_old(mitotoskel, colname="nodelbls"):
    u_mitoids = np.unique(mitotoskel.mitoids)

    countdf = mitotoskel.groupby(["mitoids", colname])["nodedists"].count()\
                  .reset_index()
    c = Counter()
    for (mitoid, nodelbl, count) in zip(countdf.mitoids,
                                        countdf[colname],
                                        countdf["nodedists"]):
        c[(mitoid, nodelbl)] = count

    lbls = np.empty((len(u_mitoids),), dtype=np.uint8)
    for (i, v) in enumerate(u_mitoids):
        maxval = 0
        maxlbl = -1
        for j in range(6):
            if c[(v, j)] > maxval:
                maxval = c[(v, j)]
                maxlbl = j

        lbls[i] = maxlbl

    return u_mitoids, lbls


def majority_vote_label(mitotoskel, colname="nodelbl"):
    u_mitoids = np.unique(mitotoskel.mitoid)

    countdf = mitotoskel.groupby(["mitoid", colname])["count"]\
                  .agg(sum).reset_index()
    c = Counter()
    for (mitoid, nodelbl, count) in zip(countdf.mitoid,
                                        countdf[colname],
                                        countdf["count"]):
        c[(mitoid, nodelbl)] = count

    lbls = np.empty((len(u_mitoids),), dtype=np.uint8)
    for (i, v) in enumerate(u_mitoids):
        maxval = 0
        maxlbl = -1
        for j in range(6):
            if c[(v, j)] > maxval:
                maxval = c[(v, j)]
                maxlbl = j

        lbls[i] = maxlbl

    return u_mitoids, lbls


def compute_node_distance_to_soma(segskel, nodelbls, binary=False):
    soma_ids = np.nonzero(nodelbls == 0)

    graph = segskel.csgraph_binary if binary else segskel.csgraph
    dists = csgraph.dijkstra(graph, directed=False, indices=soma_ids)

    return np.min(np.min(dists, axis=0), axis=0)

    
def find_segment_components(segskel, bin_centers, complbl):
    components = dict()
    component_centers = dict()
    component_labels = dict()
    
    lbls = np.unique(complbl)
    centers = np.unique(bin_centers[~np.isnan(bin_centers)])
    for lbl in lbls:
        compmsk = complbl == lbl
        for center in centers:
            segmentmsk = np.logical_and(compmsk, bin_centers == center)
            segmentids = np.flatnonzero(segmentmsk)
            ncs, cs = graph_components(segskel.csgraph, segmentids)

            for i in range(ncs):
                component = segmentids[cs == i]
                componentid = min(component)
                component_centers[componentid] = center
                component_labels[componentid] = lbl
                for j in component:
                    components[j] = componentid

    return components, component_centers, component_labels


def graph_components(graph, subgraph_inds=None):
    if subgraph_inds is not None:
        graph = graph[subgraph_inds, :][:, subgraph_inds]

    return csgraph.connected_components(graph, directed=False)    
