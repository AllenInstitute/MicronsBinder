"""
Mitochondrial density analysis

For each compartment branch, we count the number of synapses associated
with each branch and compute its total mitochondrial volume. We then
compute synapse and mitochondrial density metrics, and write out
two dataframes as a result - one describing a branch per row, and
another describing a cell per row, where each branch of a compartment
is summed within each cell.

The mitochondrial segmentation was performed in a larger volume than the
cellular segmentation. This means that several mitochondria poke out from
the tips of branches at the ends of the volume. To quantify mitochondrial
density accurately at these points, we only count the volume that overlaps
with a given segmentation ID.
"""
import sys
import argparse
from collections import Counter

import numpy as np
import pandas as pd

sys.path.append(".")
from lib import skel, compartment, u


MITOSTATS_FILENAME = "data/fullmitostats.csv"
ID_FILENAME = "data/pycids.csv"
MITOTOSKEL_FILENAME = "data/mitotoskel.h5"
POSTSYN_FILENAME = "data/pycinputs.csv"
BRANCHSTATS_FILENAME = "data/basicbranchstats.csv"
OVERLAP_FILENAME = "data/mitochondriaoverlap.csv"
NUCLEUSSTATS_FILENAME = "data/pni_nucleus_segments.csv"
NUCLEUSLOOKUP_FILENAME = "data/pni_nucleus_lookup.csv"

VOXELRES = [3.58, 3.58, 40]  # voxel resolution in nm
NUCLEUSVOXELRES = [VOXELRES[0]*16, VOXELRES[1]*16, VOXELRES[2]]
# volume of a voxel in um^3
VOXELVOL = (VOXELRES[0]/1000.
            * VOXELRES[1]/1000.
            * VOXELRES[2]/1000.)
NUCLEUSVOXELVOL = (NUCLEUSVOXELRES[0]/1000.
                   * NUCLEUSVOXELRES[1]/1000.
                   * NUCLEUSVOXELRES[2]/1000.)


def main(mitostats_filename, id_filename, mitotoskel_filename,
         postsyn_filename, basicbranchstats_filename, overlap_filename,
         nucleusstats_filename, nucleuslookup_filename, branchtype,
         fullbranchstats_filename, cellwisestats_filename):

    print("Reading data")
    mitodf = pd.read_csv(mitostats_filename, index_col=0)
    mitotoskeldf = compartment.readmitotoskel(mitotoskel_filename)
    postsyndf = pd.read_csv(postsyn_filename, index_col=0)
    ids = u.readids(id_filename)
    overlapdf = formatoverlap(pd.read_csv(overlap_filename, index_col=0))
    branchdf = pd.read_csv(basicbranchstats_filename, index_col=0)
    nucleusdf = pd.read_csv(nucleusstats_filename, index_col=0)
    nucleuslookup = u.read_lookup(nucleuslookup_filename)

    mitodf = pd.merge(mitodf, overlapdf.drop("cellid", axis=1),
                      left_index=True, right_index=True)

    if branchtype == "full":
        branchfn = skel.branches
    elif branchtype == "branchpointsegments":
        branchfn = skel.branchsegments
    elif branchtype == "distancebins":
        branchfn = skel.distancebins
    elif branchtype == "distalbranches":
        branchfn = skel.distalbranches
    elif branchtype == "coveredornot":
        branchfn = skel.coveredornot
    else:
        raise ValueError(f"unknown branch type: {branchtype}")

    print("")
    print(len(branchdf))
    print("Adding synapse counts")
    branchdf = addsynapsecounts(branchdf, ids, postsyndf,
                                mitotoskeldf, groupfn=branchfn)
    print("")
    print(len(branchdf))
    print("Adding compartment branch mitochondrion volume")
    branchdf = addmitovolumes(branchdf, ids, mitodf,
                              mitotoskeldf, groupfn=branchfn)
    print("")
    print(len(branchdf))
    print("Adding mito coverage factor")
    branchdf = addmitocoverage(branchdf, ids, mitodf,
                               mitotoskeldf, groupfn=branchfn)
    print("")
    print(len(branchdf))
    print("Computing densities")
    branchdf = computedensities(branchdf)
    print("")
    print(len(branchdf))

    print("Summing across compartments")
    cellwisedf = sumcompartments(branchdf, nucleusdf, nucleuslookup)
    print("Adding compartment branch mitochondrion volume")
    cellwisedf = addmitovolumes(cellwisedf, ids, mitodf, mitotoskeldf,
                                dfgroupname="complbl", groupfn=compartments)
    print("Adding mito coverage factor")
    cellwisedf = addmitocoverage(cellwisedf, ids, mitodf, mitotoskeldf,
                                 dfgroupname="complbl", groupfn=compartments)
    print("Computing densities")
    cellwisedf = computedensities(cellwisedf, volumecolname="cytovol")

    print("Writing results")
    branchdf.to_csv(fullbranchstats_filename)
    cellwisedf.to_csv(cellwisestats_filename)


def formatoverlap(overlapdf):
    """Changing some formatting to match the other data frames"""
    overlapdf = overlapdf.rename({
        "row_id": "mito_id",
        "col_id": "cellid",
        "vals": "overlap_vx"
    }, axis=1)

    # scaling to match mitodf["mito_vx"] voxel res
    overlapdf["overlap_vx"] *= 4

    overlapdf = overlapdf.set_index("mito_id")

    return overlapdf


def compute_and_merge(branchdf, dfgroupname, fn, *args, **kwargs):
    """
    Runs a function and merges the result into
    the branch dataframe by cell ID and branch ID
    """
    newbranchstatdf = fn(*args, **kwargs)

    return pd.merge(branchdf, newbranchstatdf,
                    left_on=["cellid", dfgroupname],
                    right_on=["cellid", "groupid"]).drop("groupid", axis=1)


def compute_for_each_id(fn, ids, colnames, *args, **kwargs):
    """
    Runs a function for a list of ids and returns
    a consolidated data frame
    """
    # init - holds the branch-wise dataframe for each cell
    subdfs = list()
    for (i, cellid) in enumerate(ids):
        print(f"Cell #{i+1} of {len(ids)}", end="\r")
        # run the function
        returnvalues = fn(cellid, *args, **kwargs)

        # match return values with column names
        if not isinstance(returnvalues, tuple):  # single return value
            returnvalues = (returnvalues,)
        else:
            assert all(len(v) == len(returnvalues[0])
                       for v in returnvalues)

        if isinstance(colnames, str):
            colnames = [colnames]

        assert len(colnames) == len(returnvalues)

        # construct dataframe of returned values
        groupids = list(returnvalues[0].keys())
        datadict = {
            "cellid": cellid,
            "groupid": groupids,
        }

        for colname, dictbygroup in zip(colnames, returnvalues):
            groupvalues = [dictbygroup[i] for i in groupids]
            datadict[colname] = groupvalues

        subdf = pd.DataFrame(datadict)

        subdfs.append(subdf)

    return pd.concat(subdfs, ignore_index=True)


def addsynapsecounts(branchdf, ids, postsyndf, mitotoskeldf,
                     dfgroupname="branchid", groupfn=skel.branches):
    """
    Counts the synapses associated with each compartment branch
    for a list of ids and merges these counts into the branch data frame.
    """
    return compute_and_merge(branchdf, dfgroupname,
                             compute_for_each_id,
                             synapsesbygroup,
                             ids, "synapsecount", postsyndf, mitotoskeldf,
                             groupfn=groupfn)


def synapsesbygroup(cellid, postsyndf, mitotoskeldf, groupfn=skel.branches):
    """
    Counts the synapses associated with each compartment branch
    of a single cell (or other group as labeled by the groupfn).
    Returns a dictionary mapping group/branch ID to a synapse count.
    """
    centroids = u.extract_coords_dlcsv(
        postsyndf[postsyndf.post_pt_root_id == cellid]) * VOXELRES

    cellskel, complbl = skel.read_skel_and_labels(cellid, refine=False)
    somadists = compartment.compute_node_distance_to_soma(
                    cellskel, complbl, binary=False)
    if groupfn == skel.coveredornot:
        coverednodes = list(set(mitotoskeldf.nodeid[
                                    mitotoskeldf.cellid == cellid]))
        grouplbl, _ = groupfn(cellskel, complbl, somadists, coverednodes)
    else:
        grouplbl, _ = groupfn(cellskel, complbl, somadists)

    kdt = u.KDTree(cellskel.vertices)
    nodecounts = Counter(kdt.query(centroids)[1])

    groupcount = Counter()
    for (node, count) in nodecounts.items():
        groupcount[grouplbl[node]] += count

    return {i: groupcount[i] for i in np.unique(grouplbl)}


def addmitovolumes(branchdf, ids, mitodf, mitotoskeldf,
                   dfgroupname="branchid", groupfn=skel.branches):
    """
    Computes the mitochondrion volume for each compartment branch
    for a list of cell ids and adds these measurements to the branch
    data frame.
    """
    return compute_and_merge(branchdf, dfgroupname,
                             compute_for_each_id,
                             mitovolbygroup,
                             ids, "mitovol",
                             mitodf, mitotoskeldf, groupfn=groupfn)


def mitovolbygroup(cellid, mitodf, mitotoskeldf, groupfn=skel.branches):
    """
    Finds the mitochondria associated with each compartment branch
    of a single cell and computes the sum of their volumes.
    Returns a dictionary mapping branch ID to a total volume measurement.
    """
    cellskel, complbl = skel.read_skel_and_labels(cellid, refine=False)
    somadists = compartment.compute_node_distance_to_soma(
                    cellskel, complbl, binary=False)
    if groupfn == skel.coveredornot:
        coverednodes = list(set(mitotoskeldf.nodeid[
                                    mitotoskeldf.cellid == cellid]))
        grouplbl, _ = groupfn(cellskel, complbl, somadists, coverednodes)
    else:
        grouplbl, _ = groupfn(cellskel, complbl, somadists)

    # Each mitochondrion gives volume to a group proportional
    # to the number of mesh associations it has with the skeleton nodes
    # of that group.
    groupdf = mitotoskeldf[mitotoskeldf.cellid == cellid].copy()
    groupdf["grouplbl"] = grouplbl[groupdf.nodeid]
    countdf = groupdf.groupby(["mitoid", "grouplbl"])["count"].sum()
    weights = countdf.groupby(level=0)\
                  .apply(lambda x: x / float(x.sum()))\
                  .reset_index()\
                  .rename({"count": "weight"}, axis=1)

    mitovols = dict()
    groupids = np.unique(grouplbl)
    for i in groupids:
        subdf = weights[weights.grouplbl == i]
        
        mitovols[i] = sum(mitodf.loc[subdf.mitoid].overlap_vx
                          * np.array(subdf.weight)
                          * VOXELVOL)
        # OLD
        #nodeids = np.flatnonzero(grouplbl == i)
        #groupcomplbl = complbl[nodeids[0]]
        #engcomplbl = compartment.ENGLISHLABELS[groupcomplbl]

        # Collect the set of mitochondria that are associated with some
        # skeleton node in this branch and match the branch compartment
        # label
        #mitoids = set(mitotoskeldf.mitoid[
        #              (mitotoskeldf.cellid == cellid) &
        #              (mitotoskeldf.nodeid.isin(nodeids))]
        #              ).intersection(
        #                  mitodf.index[mitodf.compartment == engcomplbl]
        #              )

        #mitovols[i] = sum(mitodf.loc[mitoids].overlap_vx) * VOXELVOL

    return mitovols


def addmitocoverage(branchdf, ids, mitodf, mitotoskeldf,
                    dfgroupname="branchid", groupfn=skel.branches):
    """
    Computes the mitochondrion volume for each compartment branch
    for a list of cell ids and adds these measurements to the branch
    data frame.
    """
    return compute_and_merge(branchdf, dfgroupname,
                             compute_for_each_id,
                             mitocoveragebygroup,
                             ids,
                             ["coveredpathlength", "uncoveredpathlength",
                              "extracoveragelength"],
                             mitodf, mitotoskeldf, groupfn=groupfn)


def mitocoveragebygroup(cellid, mitodf, mitotoskeldf, groupfn=skel.branches):
    cellskel, complbl = skel.read_skel_and_labels(cellid, refine=False)
    somadists = compartment.compute_node_distance_to_soma(
                    cellskel, complbl, binary=False)
    if groupfn == skel.coveredornot:
        coverednodes = list(set(mitotoskeldf.nodeid[
                                    mitotoskeldf.cellid == cellid]))
        grouplbl, _ = groupfn(cellskel, complbl, somadists, coverednodes)
    else:
        grouplbl, _ = groupfn(cellskel, complbl, somadists)

    cellmitotoskel = mitotoskeldf[mitotoskeldf.cellid == cellid]

    assocmitolookup = dict()
    for (nodeid, mitoid) in zip(cellmitotoskel.nodeid, cellmitotoskel.mitoid):
        if nodeid not in assocmitolookup:
            assocmitolookup[nodeid] = {mitoid}
        else:
            assocmitolookup[nodeid].add(mitoid)

    factors = list()
    emptyset = set()
    for (n1, n2) in cellskel.edges:
        set1 = assocmitolookup.get(n1, emptyset)
        set2 = assocmitolookup.get(n2, emptyset)
        factors.append(len(set1.intersection(set2)))
    factors = np.array(factors)

    coveredpathlength = dict()
    uncoveredpathlength = dict()
    extracoveragelength = dict()
    for i in np.unique(grouplbl):
        nodeids = np.flatnonzero(grouplbl == i)
        edgemask = np.isin(cellskel.edges, nodeids).min(1)
        edges = cellskel.edges[edgemask]

        dists = np.linalg.norm(cellskel.vertices[edges[:, 1]]
                               - cellskel.vertices[edges[:, 0]], axis=1)
        groupedgefactors = factors[edgemask]
        coveredpathlength[i] = sum(dists * groupedgefactors) / 1000.  # to um
        uncoveredpathlength[i] = sum(dists * (groupedgefactors == 0)) / 1000.
        extracoveragelength[i] = sum(dists
                                     * (groupedgefactors > 1)
                                     * (groupedgefactors - 1)) / 1000.

    return coveredpathlength, uncoveredpathlength, extracoveragelength


def computedensities(df, volumecolname="volume"):
    """Computing mitochondrial and synapse density

    Mitochondrial volume is normalized by volume, and synapse density
    is normalized by surface area
    """
    df["mitovoldensity"] = df["mitovol"] / df[volumecolname]
    df["linearmitocoverage"] = df["coveredpathlength"] / df["pathlength"]
    df["surfsynapsedensity"] = df["synapsecount"] / df["surfacearea"]
    df["linearsynapsedensity"] = df["synapsecount"] / df["pathlength"]

    return df


def sumcompartments(branchdf, nucleusdf, nucleuslookup):
    """Agglomerating branch data into cells and computing cytosolic volume"""
    celldf = branchdf.groupby(["cellid", "complbl"]).sum()[
        ["volume", "synapsecount", "surfacearea", "pathlength"]].reset_index()

    # Adding nucleus volume
    nucleusdf["vol"] = nucleusdf["size"] * NUCLEUSVOXELVOL
    nucleusvol = [nucleusdf["vol"].loc[nucleuslookup[i]]
                  for i in celldf.loc[celldf.complbl == 0, "cellid"]]
    celldf.loc[celldf.complbl == 0, "nucleusvol"] = nucleusvol
    celldf.loc[celldf.complbl != 0, "nucleusvol"] = 0

    # Computing cytosolic volume (volume - nucleusvol)
    celldf["cytovol"] = celldf["volume"] - celldf["nucleusvol"]

    return celldf


def compartments(cellskel, complbl, somadists):
    return complbl, None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("mitostats_filename", default=MITOSTATS_FILENAME)
    parser.add_argument("id_filename", default=ID_FILENAME)
    parser.add_argument("mitotoskel_filename", default=MITOTOSKEL_FILENAME)
    parser.add_argument("postsyn_filename", default=POSTSYN_FILENAME)
    parser.add_argument("basicbranchstats_filename",
                        default=BRANCHSTATS_FILENAME)
    parser.add_argument("overlap_filename", default=OVERLAP_FILENAME)
    parser.add_argument("nucleusstats_filename", default=NUCLEUSSTATS_FILENAME)
    parser.add_argument("nucleuslookup_filename", default=NUCLEUSLOOKUP_FILENAME)
    parser.add_argument("fullbranchstats_filename")
    parser.add_argument("cellwisestats_filename")
    parser.add_argument("--branchtype", default="full")

    main(**vars(parser.parse_args()))
