"""
Skeleton handling

Label conventions:
{0: 'Somatic',
 1: 'Axonal',
 2: 'Basal dendrite',
 3: 'Apical dendrite',
 4: 'Unknown/Ambiguous dendrite',
 5: 'Unknown/Ambiguous'}
"""
import os
import json

import h5py
import numpy as np
import cloudvolume as cv
from meshparty import skeleton_io
from meshparty.skeleton import Skeleton
from scipy.sparse import csgraph

from . import ngl


CVOL = cv.CloudVolume(ngl.MITO_CVPATH, parallel=8, progress=False)
SKELDIR = "../data/smoothed_skeletons_v185"
SKEL_ASSOC_DIR = "data/mitotoskel"


def read_neuron_skel(segid, scale=True):
    filename = f"{SKELDIR}/{segid}_skeleton.h5"

    return skeleton_io.read_skeleton_h5(filename)


def read_neuron_complbls(segid):
    filename = f"{SKELDIR}/{segid}_skeleton_label.npy"

    return np.load(smoothed_filename)


def read_skel_and_labels(segid, refine=False):
    segskel = skel.read_neuron_skel(segid)
    complbl = skel.read_neuron_complbls(segid)

    if refine:
        complbl = compartment.refine_labels(segskel, complbl)

    return segskel, complbl


def read_skel_assoc(segid, include_counts=False):
    filename = f"{SKEL_ASSOC_DIR}/assoc__{segid}.h5"
    with h5py.File(filename, "r") as f:
        data = f[f"{segid}/data"][()]
        mitoids = f[f"{segid}/mitoids"][()]
        ptrs = f[f"{segid}/ptrs"][()]

    assocs = dict()
    for (i, v) in enumerate(mitoids):
        rng = ptrs[i, :]
        if include_counts:
            assocs[v] = data[rng[0]:rng[1]]
        else:
            assocs[v] = list(map(tuple, data[rng[0]:rng[1]]))

    return assocs


def path_length(segskel, r, binary=False):
    graph = segskel.csgraph_binary if binary else segskel.csgraph

    return csgraph.dijkstra(graph, directed=False, indices=r)


def path_length_to_leaves(segskel, binary=False, return_inds=False):
    pls = path_length(segskel, segskel.end_points, binary=binary)

    if return_inds:
        return np.min(pls, axis=0), segskel.end_points[np.argmin(pls, axis=0)]
    else:
        return np.min(pls, axis=0)


def associate_to_node(distdf, summary_fn, colname="mitoids"):
    nodeids = dict()

    for (i, subdf) in distdf.groupby(colname):
        summ_dist = summary_fn(subdf.nodedists)
        j = np.argmin(np.abs(np.array(subdf.nodedists) - summ_dist))
        nodeids[i] = subdf.nodeids.tolist()[j]

    return nodeids


def full_cover_paths(segskel):
    return [segskel.path_to_root(ep) for ep in segskel.end_points]


def edge_weights(segskel):
    xs = segskel.vertices[segskel.edges[:, 0]]
    ys = segskel.vertices[segskel.edges[:, 1]]
    return np.linalg.norm(xs - ys, axis=1)


def branches(skel, compartmentlabels, somadists=None):
    """Splitting a skeleton into branches

    Each branch is defined as the connected components
    of the skeleton with the same compartment label.
    """
    branchids = np.empty((len(skel.vertices),), dtype=np.int)
    branchlbl = dict()

    lbls = np.unique(compartmentlabels)
    for lbl in lbls:
        compartmentids = np.flatnonzero(compartmentlabels == lbl)
        ncs, cs = graphcomponents(skel.csgraph, compartmentids)

        for i in range(ncs):
            component = compartmentids[cs == i]

            if somadists is None:
                branchid = min(component)
            else:
                branchid = component[np.argmin(somadists[component])]

            branchids[component] = branchid
            branchlbl[branchid] = lbl

    return branchids, branchlbl


def graph_components(graph, subgraph_inds=None):
    if subgraph_inds is not None:
        graph = graph[subgraph_inds, :][:, subgraph_inds]

    return csgraph.connected_components(graph, directed=False)    
