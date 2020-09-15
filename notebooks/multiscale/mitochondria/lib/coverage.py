"""Computing mitochondrial coverage"""
import numpy as np
import pandas as pd
from scipy.sparse import csgraph

from . import skel
from . import compartment


def compute_mitocovfactor(
    distdf, cellids=None, binwidth=20, refine=True,
    shrink=True, shrinkthresh=10, verbose=True):

    cellids = np.unique(distdf.cellid) if cellids is None else cellids
    subdfs = list()

    for (i, cellid) in enumerate(cellids):
        if verbose:
            print(f"#{i+1} of {len(cellids)}", end="\r")

        subdistdf = distdf[distdf.cellid == cellid]
        subdf, _ = mitocovfactor(
                       cellid, subdistdf, binwidth=binwidth, refine=refine,
                       shrink=shrink, shrinkthresh=shrinkthresh)

        cellid_series = pd.Series([cellid] * len(subdf), dtype=np.uint64)
        subdf.loc[:, "cellid"] = cellid_series
        subdfs.append(subdf)

    return pd.concat(subdfs, ignore_index=True)


def compute_bulk_mitocovfactor(
    distdf, cellids=None, refine=True,
    shrink=True, shrinkthresh=10, verbose=True):

    cellids = np.unique(distdf.cellid) if cellids is None else cellids
    subdfs = list()

    for (i, cellid) in enumerate(cellids):
        if verbose:
            print(f"#{i+1} of {len(cellids)}", end="\r")

        subdistdf = distdf[distdf.cellid == cellid]
        subdf = bulk_mitocovfactor(
                    cellid, subdistdf, refine=refine,
                    shrink=shrink, shrinkthresh=shrinkthresh)

        cellid_series = pd.Series([cellid] * len(subdf), dtype=np.uint64)
        subdf.loc[:, "cellid"] = cellid_series
        subdfs.append(subdf)

    return pd.concat(subdfs, ignore_index=True)


def path_length_near_somas(
    cellids, distthresh=50, refine=True,
    shrink=True, shrinkthresh=10, verbose=True):

    subdfs = list()

    for (i, cellid) in enumerate(cellids):
        if verbose:
            print(f"#{i+1} of {len(cellids)}", end="\r")

        subdf = path_length_near_soma(
                    cellid, refine=refine, shrink=shrink,
                    shrinkthresh=shrinkthresh)

        cellid_series = pd.Series([cellid] * len(subdf), dtype=np.uint64)
        subdf.loc[:, "cellid"] = cellid_series
        subdfs.append(subdf)

    return pd.concat(subdfs, ignore_index=True)


def compute_bulk_coverage(
    distdf, cellids=None, refine=True,
    shrink=True, shrinkthresh=10, verbose=True):

    cellids = np.unique(distdf.cellid) if cellids is None else cellids
    subdfs = list()

    for (i, cellid) in enumerate(cellids):
        if verbose:
            print(f"#{i+1} of {len(cellids)}", end="\r")

        covered = distdf.nodeids[distdf.cellid == cellid].values
        subdf = bulk_mitochondrial_coverage_by_comp(
                    cellid, covered, refine=refine,
                    shrink=shrink, shrinkthresh=shrinkthresh)

        cellid_series = pd.Series([cellid] * len(subdf), dtype=np.uint64)
        subdf.loc[:, "cellid"] = cellid_series
        subdfs.append(subdf)

    return pd.concat(subdfs, ignore_index=True)


def bulk_mitochondrial_coverage_by_comp(
    cellid, coveredids, refine=True,
    shrink=True, shrinkthresh=10, synapse_nodes=None):

    cellskel, complbl = read_skel_and_labels(cellid, refine=refine)

    disttoleaf = skel.path_length(
                     cellskel, cellskel.end_points).min(axis=0) / 1000.

    edgedf = construct_edge_frame(cellskel, coveredids, disttoleaf=disttoleaf)

    edgedf.loc[:, "nodelbl1"] = [complbl[i] for i in edgedf["nodeid1"]]
    edgedf.loc[:, "nodelbl2"] = [complbl[i] for i in edgedf["nodeid2"]]

    edgedf = edgedf[edgedf.nodelbl1 == edgedf.nodelbl2]

    if shrink:
        edgedf = edgedf[edgedf.disttoleaf > shrinkthresh]

    covdf = edgedf[(edgedf.covered1 == 1) &
                   (edgedf.covered2 == 1)]

    compartmentids = np.unique(edgedf.nodelbl1)
    coveredlen = dict(covdf.groupby("nodelbl1")["pathlength"].sum())
    fulllen = dict(edgedf.groupby("nodelbl2")["pathlength"].sum())

    coverage = [coveredlen.get(i, 0.) / fulllen[i] for i in compartmentids]
    pathlength = [fulllen[i] / 1000. for i in compartmentids]

    data_dict = {"coverage": coverage,
                 "compartment": compartmentids,
                 "pathlength": pathlength}

    if synapse_nodes is not None:
        synapse_counts = {i: 0 for i in compartmentids}
        for (node, assoc_synapses) in synapse_nodes.items():
            nodelbl = complbl[node]
            if nodelbl not in synapse_counts:
                continue
            synapse_counts[nodelbl] += assoc_synapses
        data_dict["synapsecount"] = [synapse_counts[i]
                                     for i in compartmentids]

    res = pd.DataFrame(data_dict)

    if synapse_nodes is not None:
        res["synapse_density"] = res["synapsecount"] / res["pathlength"]

    return res


def path_length_near_soma(
    cellid, distthresh=50, refine=True, shrink=True, shrinkthresh=10):

    cellskel, complbl = read_skel_and_labels(cellid, refine)
    disttosoma = compartment.compute_node_distance_to_soma(
                    cellskel, complbl) / 1000.
    disttoleaf = skel.path_length_to_leaves(cellskel) / 1000.

    edgedf = construct_edge_frame(
                 cellskel, [], disttosoma=disttosoma, disttoleaf=disttoleaf)

    if shrink:
        edgedf = edgedf[edgedf.disttoleaf > shrinkthresh]

    edgedf.loc[:, "nodelbl1"] = [complbl[i] for i in edgedf["nodeid1"]]
    edgedf.loc[:, "nodelbl2"] = [complbl[i] for i in edgedf["nodeid2"]]

    edgedf = edgedf[edgedf.nodelbl1 == edgedf.nodelbl2]
    neardf = edgedf[edgedf.disttosoma < distthresh]

    nearpl = neardf.groupby("nodelbl1").sum()["pathlength"]
    allpl = edgedf.groupby("nodelbl1").sum()["pathlength"]

    return (nearpl / allpl).reset_index()


def mitocovfactor(
    cellid, distdf, binwidth=50, refine=True,
    shrink=False, shrinkthresh=5):

    cellskel, complbl = read_skel_and_labels(cellid, refine=refine)
    disttosoma = compartment.compute_node_distance_to_soma(
                    cellskel, complbl) / 1000.
    disttoleaf = skel.path_length_to_leaves(cellskel) / 1000.

    # Can cause inf -> NaN
    prevsettings = np.seterr(divide="ignore", invalid="ignore")
    bin_centers = disttosoma // binwidth * binwidth + binwidth // 2
    np.seterr(**prevsettings)

    edgedf = construct_edge_frame(
                 cellskel, [], disttoleaf=disttoleaf)

    # Finding connected components within distance bounds and
    #  adding their info to the dataframe
    cs, cbins, clbls = find_segment_components(cellskel, bin_centers, complbl)
    edgedf.loc[:, "component1"] = [cs.get(i, -1) for i in edgedf["nodeid1"]]
    edgedf.loc[:, "component2"] = [cs.get(i, -1) for i in edgedf["nodeid2"]]

    if shrink:
        edgedf = edgedf[edgedf.disttoleaf > shrinkthresh]

    edgedf.loc[:, "nodelbl1"] = [complbl[i] for i in edgedf["nodeid1"]]
    edgedf.loc[:, "nodelbl2"] = [complbl[i] for i in edgedf["nodeid2"]]

    edgedf = edgedf[edgedf.component1 == edgedf.component2]
    edgedf = edgedf[edgedf.component1 != -1]

    edgedf["factor"] = coverage_factor(edgedf, distdf)
    edgedf["scaled"] = edgedf["pathlength"] * edgedf["factor"]

    pls = edgedf.groupby("component1").sum()[["pathlength", "scaled"]]
    res = (pls["scaled"] / pls["pathlength"]).reset_index()
    res["bin_center"] = [cbins[i] for i in res.component1]
    res["nodelbl"] = [clbls[i] for i in res.component1]

    minframe = edgedf.groupby("component1").min()
    leafdist_lookup = dict(zip(minframe.index, minframe.disttoleaf))
    res["disttoleaf"] = [leafdist_lookup[i] for i in res.component1]

    return res.rename(
               {0: "mitocovfactor", "component1": "componentid"}, axis=1), cs


def bulk_mitocovfactor(cellid, distdf,
                       refine=True, shrink=True, shrinkthresh=10):

    cellskel, complbl = read_skel_and_labels(cellid, refine=refine)
    disttoleaf = skel.path_length_to_leaves(cellskel) / 1000.

    edgedf = construct_edge_frame(
                 cellskel, [], disttoleaf=disttoleaf)

    if shrink:
        edgedf = edgedf[edgedf.disttoleaf > shrinkthresh]

    edgedf.loc[:, "nodelbl1"] = [complbl[i] for i in edgedf["nodeid1"]]
    edgedf.loc[:, "nodelbl2"] = [complbl[i] for i in edgedf["nodeid2"]]

    edgedf = edgedf[edgedf.nodelbl1 == edgedf.nodelbl2]

    edgedf["factor"] = coverage_factor(edgedf, distdf)
    edgedf["scaled"] = edgedf["pathlength"] * edgedf["factor"]

    pls = edgedf.groupby("nodelbl1").sum()[["pathlength", "scaled"]]
    res = (pls["scaled"] / pls["pathlength"]).reset_index()

    return res.rename({0: "mitocovfactor"}, axis=1)


def construct_edge_frame(
    segskel, covered_nodeids, **distances):
    nodeid1 = segskel.edges[:, 0]
    nodeid2 = segskel.edges[:, 1]

    pathlength = np.linalg.norm(
                     segskel.vertices[nodeid1] - segskel.vertices[nodeid2],
                     axis=1)

    covered_nodeids = list(covered_nodeids)
    covered1 = np.isin(nodeid1, covered_nodeids)
    covered2 = np.isin(nodeid2, covered_nodeids)

    data_dict = dict(nodeid1=nodeid1, nodeid2=nodeid2,
                     pathlength=pathlength,
                     covered1=covered1, covered2=covered2)

    for (name, distance) in distances.items():
        data_dict[name] = np.minimum(distance[nodeid1],
                                     distance[nodeid2])

    return pd.DataFrame(data_dict)


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


def coverage_factor(edgedf, distdf):

    mitoid_lookup = dict()
    for (nodeid, mitoid) in zip(distdf.nodeids, distdf.mitoids):
        if nodeid not in mitoid_lookup:
            mitoid_lookup[nodeid] = {mitoid}
        else:
            mitoid_lookup[nodeid].add(mitoid)

    factors = list()
    emptyset = set()
    for (nodeid1, nodeid2) in zip(edgedf.nodeid1, edgedf.nodeid2):
        set1 = mitoid_lookup.get(nodeid1, emptyset)
        set2 = mitoid_lookup.get(nodeid2, emptyset)
        factors.append(len(set1.intersection(set2)))

    return np.array(factors)


def read_skel_and_labels(segid, refine=True):
    segskel = skel.read_neuron_skel(segid)
    complbl = skel.read_neuron_complbls(segid)
    if refine:
        complbl = compartment.refine_labels(segskel, complbl)

    return segskel, complbl
