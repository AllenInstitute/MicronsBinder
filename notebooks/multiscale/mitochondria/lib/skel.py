"""
Skeleton handling

Label conventions:
{1: 'axon',
 2: 'basal dendrite',
 3: 'apical dendrite',
 4: 'ambiguous dendrite',
 5: 'ambiguous',
 0: 'soma'}
"""
import os
import json

import h5py
import numpy as np
import cloudvolume as cv
from meshparty import skeleton_io
from meshparty.skeleton import Skeleton
from scipy.sparse import csgraph


SKELDIR = "data/smoothed_skeletons_v185"
SKEL_ASSOC_DIR = "../../mito-analysis/intermeds/mito_to_skel"


def read_neuron_skel(segid, scale=False):
    filename = f"{SKELDIR}/{segid}_skeleton.h5"

    if os.path.exists(filename):
        skel = skeleton_io.read_skeleton_h5(filename)
    else:
        raise Exception(f"{filename} not found. Have you set the SKELDIR?")

    if scale:
        skel.voxel_scaling = [0.895, 0.895, 1]

    return skel


def read_smoothed_skel(segid):
    filename = f"{SKELDIR}/{segid}_skeleton.h5"

    return skeleton_io.read_skeleton_h5(filename)


def read_neuron_complbls(segid):
    smoothed_filename = f"{SKELDIR}/{segid}_skeleton_label.npy"

    if os.path.exists(smoothed_filename):
        return np.load(smoothed_filename)
    else:
        raise Exception(f"{filename} not found. Have you set the SKELDIR?")


def read_smoothed_complbls(segid):
    filename = f"{SKELDIR}/{segid}_skeleton_label.npy"

    return np.load(filename)


def skel_length(segid):
    return CVOL.skeleton.get(segid).cable_length()


def read_skel_assoc(segid):
    filename = f"{SKEL_ASSOC_DIR}/assoc_{segid}.h5"
    with h5py.File(filename, "r") as f:
        data = f[f"{segid}/data"][()]
        mitoids = f[f"{segid}/mitoids"][()]
        ptrs = f[f"{segid}/ptrs"][()]

    assocs = dict()
    for (i, v) in enumerate(mitoids):
        rng = ptrs[i, :]
        assocs[v] = data[rng[0]:rng[1]]

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


def remove_dups(segskel):
    weights = edge_weights(segskel)

    to_merge = np.max(segskel.edges[weights == 0], axis=1)
    targets = np.min(segskel.edges[weights == 0], axis=1)

    vertexmask = np.ones((len(segskel.vertices),), dtype=np.uint8)
    vertexmask[to_merge] = 0
    vertexinds = np.flatnonzero(vertexmask)

    newverts = segskel.vertices[vertexinds]
    newvertprops = {k: np.array(v)[vertexinds] for (k, v)
                    in segskel.vertex_properties.items()}

    remapping = np.zeros((len(segskel.vertices),), dtype=np.uint32)
    remapping[vertexinds] = np.arange(len(newverts))
    remapping[to_merge] = remapping[targets]

    newedges = remapping[segskel.edges[weights != 0]]
    newm2s = remapping[segskel.mesh_to_skel_map]

    return Skeleton(vertices=newverts, edges=newedges,
                    mesh_to_skel_map=newm2s,
                    vertex_properties=newvertprops,
                    root=remapping[segskel.root])


def edge_weights(segskel):
    xs = segskel.vertices[segskel.edges[:, 0]]
    ys = segskel.vertices[segskel.edges[:, 1]]
    return np.linalg.norm(xs - ys, axis=1)
