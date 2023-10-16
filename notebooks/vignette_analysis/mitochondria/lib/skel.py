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


SKELDIR = "../data/smoothed_skeletons_v185"
SKEL_ASSOC_DIR = "data/mitotoskel"


def read_neuron_skel(segid, scale=True):
    filename = f"{SKELDIR}/{segid}_skeleton.h5"

    return skeleton_io.read_skeleton_h5(filename)


def read_neuron_complbls(segid):
    filename = f"{SKELDIR}/{segid}_skeleton_label.npy"

    return np.load(smoothed_filename)


def read_skel_and_labels(segid, refine=False):
    segskel = read_neuron_skel(segid)
    complbl = read_neuron_complbls(segid)

    if refine:
        complbl = compartment.refine_labels(segskel, complbl)

    return segskel, complbl


def read_skel_assoc(segid, include_counts=False):
    filename = f"{SKEL_ASSOC_DIR}/assoc_{segid}.h5"
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


def path_length_to_root(segskel, r, binary=False):
    graph = segskel.csgraph_binary if binary else segskel.csgraph

    return csgraph.dijkstra(graph, directed=False, indices=r)


def path_length_to_leaves(segskel, binary=False, return_inds=False):
    pls = path_length_to_root(segskel, segskel.end_points, binary=binary)

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


def distalbranches(skel, compartmentlabels, somadists, distthresh=60_000):
    templabels = np.copy(compartmentlabels)
    templabels[np.logical_and(somadists < distthresh,
                              compartmentlabels != 0)] = -1
    templabels[np.isinf(somadists)] = -2

    return branches(skel, templabels, somadists)


def branchsegments(skel, compartmentlabels, somadists=None):
    templabels = np.zeros((len(skel.vertices),), dtype=np.int)
    paths = full_cover_paths(skel)
    bps = skel.branch_points

    for path in paths:
        bps_crossed = np.cumsum(np.isin(path, bps))
        templabels[path] = np.maximum(templabels[path], bps_crossed)

    # Making the soma a separate component
    templabels[compartmentlabels == 0] = -1

    return branches(skel, templabels, somadists)


def distancebins(skel, compartmentlabels, somadists, resolution=10_000):
    templabels = (somadists // resolution).astype(int)
    templabels[compartmentlabels == 0] = -1

    branchids, templbl = branches(skel, templabels, somadists)
    branchlbl = {i: compartmentlabels[i] for i in np.unique(branchids)}

    return branchids, branchlbl


def coveredornot(skel, compartmentlabels, somadists, coverednodes):
    nodecovered = np.isin(np.arange(len(skel.vertices)), coverednodes)

    # use the negative soma distances to label the node furthest from
    # the soma instead of the closest
    branchids_, branchlbl_ = branches(skel, nodecovered, -somadists)

    # Adding 6 to the label for uncovered components
    branchlbl = {k: compartmentlabels[k] + 6*(1-branchlbl_[k])
                 for k in branchlbl_.keys()}

    return branchids_, branchlbl


def graphcomponents(graph, subgraph_inds=None):
    if subgraph_inds is not None:
        graph = graph[subgraph_inds, :][:, subgraph_inds]

    return csgraph.connected_components(graph, directed=False)    


def pathlength(cellskel, grouplbl=None, groupid=None, scale=1000.):
    """Computes the path length of the edges within a group (often a branch)"""
    if grouplbl is None or groupid is None:
        groupedges = np.arange(len(cellskel.edges))

    else:
        vertinds = np.flatnonzero(grouplbl == groupid)

        groupedges = np.flatnonzero(np.isin(cellskel.edges, vertinds).min(1))

    distances = np.linalg.norm(
                    cellskel.vertices[cellskel.edges[groupedges, 1]]
                    - cellskel.vertices[cellskel.edges[groupedges, 0]],
                    axis=1)

    return sum(distances) / scale
