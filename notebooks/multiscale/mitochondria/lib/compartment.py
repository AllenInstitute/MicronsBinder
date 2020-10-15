"""Labeling mitochondria and neurons by compartment (axonal, dendritic, somatic)"""
from collections import Counter

import numpy as np
import pandas as pd
from scipy.spatial import KDTree
from scipy.sparse import csgraph
import h5py

from . import skel
from . import u


DIST_DATA_FILENAME = "intermeds/200702_pyr_dist_data_refine.h5"


def label_mitos_skelcompdf(mitodf):
    soma_ids = list(set(mitodf.v185_rt_id))

    labels = np.array([""] * len(mitodf))
    label_mapping = {0: "axon", 1: "dend", 2: "soma"}
    for (i, soma_id) in enumerate(soma_ids):
        print(f"#{i+1} of {len(soma_ids)}: {soma_id}", end="\r")
        mito_mask = mitodf.v185_rt_id == soma_id
        mito_pts = u.extract_centroids(mitodf[mito_mask])
        neuron_skel = skel.read_neuron_skel(soma_id)
        neuron_complbl = skel.read_neuron_complbls(soma_id)

        mito_labels = label_mitos_skelcomp(
                          mito_pts, neuron_skel, neuron_complbl)

        labels[mito_mask] = [label_mapping[v] for v in mito_labels]
    print("")

    return labels


def label_mitos_skelcomp(mito_pts, neuron_skel, skel_labels):
    """
    A slightly less rough scheme to label mitochondria as axonal,
    dendritic, or somatic. Mitochondria whose coordinate is closer to
    a point with a particular label adopts that point

    input label conventions:
    {1: 'axon',
     2: 'basal dendrite',
     3: 'apical dendrite',
     4: 'ambiguous dendrite',
     5: 'ambiguous',
     0: 'soma'}

    Args:
        mito_pts (2darray):
        neuron_skel (meshparty.skeleton or cloudvolume.skeleton):
        skel_labels (1darray):

    Returns:
        labels (1darray): 0 => axonal, 1 => dendritic, 2 => somatic
    """
    axon_pts = neuron_skel.vertices[skel_labels == 1]
    dend_pts = neuron_skel.vertices[np.isin(skel_labels, [2, 3, 4])]
    soma_pts = neuron_skel.vertices[skel_labels == 0]

    return closest_pt_set(mito_pts, axon_pts, dend_pts, soma_pts)


def label_mitos_prepost(mito_pts, soma_pt, presyn_pts, postsyn_pts,
                        soma_thresh):
    """
    A rough scheme to label mitochondria as axonal, dendritic, or somatic.
    Mitochondria whose coordinate is within soma_thresh distance of the soma_pt
    is labeled as somatic, those who are closer to a presynaptic point than
    a postsynaptic point are labeled axonal, and others are dendritic.

    Args:
        mito_pts (2darray):
        soma_pt (tuple):
        presyn_pts (2darray):
        postsyn_pts (2darray):
        soma_thresh (float):

    Returns:
        labels (1darray): 0 => axonal, 1 => dendritic, 2 => somatic
    """
    labels = closest_pt_set(mito_pts, presyn_pts, postsyn_pts)
    labels[dist_to_pt(mito_pts, soma_pt) <= soma_thresh] = 2

    return labels


def closest_pt_set(pts, *pt_sets):
    """Label points {pts} by the set which has the nearest neighbor"""
    labels = np.empty((len(pts),), dtype=np.uint8)

    # create label lookup
    label_lookup = dict()
    for (i, pt_set) in enumerate(pt_sets):
        for pt in pt_set:
            label_lookup[tuple(pt)] = i

    all_sets = np.concatenate(pt_sets)
    kdt = KDTree(all_sets)

    for (i, pt) in enumerate(pts):
        closest_pt = all_sets[kdt.query(pt)[1], :]
        labels[i] = label_lookup[tuple(closest_pt)]

    return labels


def dist_to_pt(pts, target_pt, voxel_res=[4, 4, 40]):
    """Compute the distance between a set of points and a single target"""
    return np.linalg.norm((pts - target_pt) * voxel_res, axis=1)


def dist_data_fields(filename=DIST_DATA_FILENAME):
    with h5py.File(filename, 'r') as f:
        return list(f.keys())


def read_dist_data(dset_name, filename=DIST_DATA_FILENAME):
    with h5py.File(filename, 'r') as f:
        return f[dset_name][()]


def read_dist_df(filename=DIST_DATA_FILENAME):
    fields = dict()
    for field in dist_data_fields(filename=filename):
        fields[field] = read_dist_data(field, filename=filename)

    return pd.DataFrame(fields)


def label_skel_comp_prepostpath(
    segskel, presyn_pts, postsyn_pts, soma_thresh=15000):
    """
    Labeling skeleton nodes by compartment using path topology and
    synaptic assignment information.

    * Assign all nodes within soma_thresh path_length distance as somatic
    * Find nodes closest to all pre/post synaptic terminals
    * Classify paths as axonal or dendritic by majority vote of prox to site type
    * Fill in non-somatic labels of each path by this classification

    0 -> somatic, 1 -> axonal, 2 -> dendritic, 3 -> unknown
    """
    skel_kdt = KDTree(segskel.vertices)

    # Find closest skeleton node for each synaptic point
    presyn_nodes = [skel_kdt.query(pt)[1] for pt in presyn_pts]
    postsyn_nodes = [skel_kdt.query(pt)[1] for pt in postsyn_pts]

    soma_dists = csgraph.dijkstra(
                     segskel.csgraph, directed=False, indices=segskel.root)

    rootlbls = np.ones((len(segskel.vertices),), dtype=np.uint8) * 3
    # Label postsyn nodes first to give presyn an advantage - generally fewer presyn
    rootlbls[postsyn_nodes] = 2
    rootlbls[presyn_nodes] = 1
    rootlbls[soma_dists < soma_thresh] = 0

    paths = skel.full_cover_paths(segskel)
    initlbls = vote_by_roots(rootlbls, paths)
    # second pass for "consistency"
    finallbls = vote_by_roots(initlbls, paths)

    return finallbls


def vote_by_roots(rootlbls, paths):
    skellbls = np.ones((len(rootlbls),), dtype=np.uint8) * 3
    for path in paths:
        pathlbls = rootlbls[path]
        c = Counter(pathlbls)
        if c[1] == 0 and c[2] == 0:
            if c[0] == 0:
                continue
            else:
                pathlbls[:] = 0

        if c[2] > c[1]:
            pathlbls[pathlbls != 0] = 2
        else:
            pathlbls[pathlbls != 0] = 1
        # Fancy indexing makes a copy?
        skellbls[path] = pathlbls

    return skellbls


def refine_labels(segskel, complbls, n_edges=1):
    """
    Restricting the soma label to a smaller distance threshold from the root node
    """
    assert n_edges >= 0, "need a nonnegative edge threshold"

    # Find all nodes within n_edges of the root
    soma_edge_dist = skel.path_length(segskel, segskel.root, binary=True)
    soma_nodes = np.nonzero(soma_edge_dist <= n_edges)[0]

    newlbls = np.zeros((len(complbls),), dtype=np.uint8)
    for path in skel.full_cover_paths(segskel):
        c = Counter(complbls[path])
        newlbls[path] = c.most_common()[0][0]
        newlbls[path[np.isin(path, soma_nodes)]] = 0

    return newlbls


def majority_vote_label(distdf):
    u_mitoids = np.unique(distdf.mitoids)

    c = Counter(zip(distdf.mitoids, distdf.nodelbls))

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


def mitochondrial_coverage_by_comp(
    segid, covered_nodeids, binwidth=50, refine=True):

    segskel = skel.read_neuron_skel(segid)
    complbl = skel.read_neuron_complbls(segid)
    if refine:
        complbl = refine_labels(segskel, complbl)
    somadists = compute_node_distance_to_soma(segskel, complbl) / 1000.

    # Can cause inf -> NaN
    prevsettings = np.seterr(divide="ignore", invalid="ignore")
    bin_centers = somadists // binwidth * binwidth + binwidth // 2
    np.seterr(**prevsettings)

    edgedf = construct_edge_frame(segskel, somadists, covered_nodeids)

    # Finding connected components within distance bounds and
    #  adding their info to the dataframe
    cs, cbins, clbls = find_segment_components(segskel, bin_centers, complbl)
    edgedf.loc[:, "component1"] = [cs.get(i, -1) for i in edgedf["nodeid1"]]
    edgedf.loc[:, "component2"] = [cs.get(i, -1) for i in edgedf["nodeid2"]]

    # Computing fractions per segment
    edgedf = edgedf[(edgedf.component1 != -1) &
                    (edgedf.component2 != -1) &
                    (edgedf.component1 == edgedf.component2)]
    covdf = edgedf[(edgedf.covered1 == 1) &
                   (edgedf.covered2 == 1)]

    componentids = np.unique(edgedf.component1)
    covsums = dict(covdf.groupby("component1").sum()["pathlength"])
    allsums = dict(edgedf.groupby("component1").sum()["pathlength"])
    fracs = [covsums.get(k, 0.) / allsums[k] for k in componentids]
    component_centers = [cbins[i] for i in componentids]
    component_labels = [clbls[i] for i in componentids]
    fracdf = pd.DataFrame({"componentid": componentids,
                           "coverage": fracs,
                           "bin_center": component_centers,
                           "nodelbl": component_labels})

    return fracdf, cs


def construct_edge_frame(segskel, somadists, covered_nodeids):

    nodeid1 = segskel.edges[:, 0]
    nodeid2 = segskel.edges[:, 1]

    pathlength = np.linalg.norm(
                     segskel.vertices[nodeid1] - segskel.vertices[nodeid2],
                     axis=1)

    covered_nodeids = list(covered_nodeids)
    covered1 = np.isin(nodeid1, covered_nodeids)
    covered2 = np.isin(nodeid2, covered_nodeids)

    return pd.DataFrame(
               dict(nodeid1=nodeid1, nodeid2=nodeid2,
                    pathlength=pathlength,
                    covered1=covered1, covered2=covered2))

    
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
