"""Mitochondrial density analysis"""
import sys
import argparse
from collections import Counter

import numpy as np
import pandas as pd

sys.path.append(".")
from lib import skel, compartment, u


MITOSTATS_FILENAME = "intermeds/200811_pinky100_mito_cellswskel.csv"
ID_FILENAME = "intermeds/200712_clean_and_complete.csv"
MITOTOSKEL_FILENAME = "intermeds/210329_dist_data_pyr_refine.h5"
POSTSYN_FILENAME = "intermeds/200527_neuron_postsyn.csv"
BRANCHSTATS_FILENAME = "intermeds/210324_basicbranchstats.csv"
OVERLAP_FILENAME = "intermeds/200509_v185_overlaps.df"
NUCLEUSSTATS_FILENAME = "intermeds/200713_nucleus_info.df"
NUCLEUSLOOKUP_FILENAME = "intermeds/200713_nucleus_lookup.csv"

VOXELRES = [3.58, 3.58, 40]  # voxel resolution in nm
# volume of a voxel in um^3
VOXELVOL = [VOXELRES[0]/1000., VOXELRES[1]/1000., VOXELRES[2]/1000.]


def main(mitostats_filename, id_filename, mitotoskel_filename,
         postsyn_filename, basicbranchstats_filename, overlap_filename,
         nucleusstats_filename, nucleuslookup_filename,
         fullbranchstats_filename, cellwisestats_filename):

    print("Reading data")
    mitodf = pd.read_csv(mitostats_filename, index_col=0)
    mitotoskeldf = compartment.read_mitotoskel(mitotoskel_filename)
    postsyndf = pd.read_csv(postsyn_filename, index_col=0)
    ids = u.readids(id_filename)
    overlapdf = formatoverlap(pd.read_csv(overlap_filename))
    branchdf = pd.read_csv(basicbranchstats_filename, index_col=0)
    nucleusdf = pd.read_csv(nucleusstats_filename, index_col=0)
    nucleuslookup = u.read_lookup(nucleuslookup_filename)

    mitodf = pd.merge(mitodf, overlapdf,
                      left_on=["mitoid", "cellid"],
                      right_on=["mitoid", "cellid"])

    print("Adding synapse counts")
    branchdf = addsynapsecounts(branchdf, ids, postsyndf)
    print("Adding compartment branch mitochondrion volume")
    branchdf = addmitovolumes(branchdf, ids, mitodf, mitotoskeldf)
    print("Computing densities")
    branchdf = computedensities(branchdf)

    print("Grouping data by cell")
    cellwisedf = groupcells(branchdf, nucleusdf, nucleuslookup)

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

    return overlapdf


def compute_and_merge(branchdf, fn, *args, **kwargs):
    """
    Runs a function and merges the result into
    the branch dataframe by cell ID and branch ID
    """
    newbranchstatdf = fn(*args, **kwargs)

    return pd.merge(branchdf, newbranchstatdf,
                    left_on=["cellid", "branchid"],
                    right_on=["cellid", "branchid"])


def compute_for_each_id(branchfn, ids, *args, **kwargs):
    """
    Runs a function for a list of ids and returns
    a consolidated data frame
    """
    # init - holds the branch-wise dataframe for each cell
    subdfs = list()
    for (i, cellid) in enumerate(ids):
        print(f"#{i+1} of {len(ids)}", end="\r")
        bybranch = branchfn(cellid, *args, **kwargs)
        branchids = list(bybranch.keys())
        branchcounts = [bybranch[i] for i in branchids]

        subdf = pd.DataFrame(dict(cellid=cellid,
                                  branchid=branchids,
                                  synapsecount=branchcounts))

        subdfs.append(subdf)

    return pd.concat(subdfs, ignore_index=True)


def addsynapsecounts(branchdf, ids, postsyndf):
    """
    Counts the synapses associated with each compartment branch
    for a list of ids and merges these counts into the branch data frame.
    """
    return compute_and_merge(compute_for_each_id, synapsesbybranch,
                             ids, postsyndf)


def synapsesbybranch(cellid, postsyndf):
    """
    Counts the synapses associated with each compartment branch
    of a single cell. Returns a dictionary mapping branch ID to
    a synapse count.
    """
    centroids = u.extract_coords_dlcsv(
        postsyndf[postsyndf.post_pt_root_id == cellid]) * VOXELRES

    cellskel, complbl = skel.read_skel_and_labels(cellid, refine=False)
    branchlbl = skel.branches(cellskel, complbl)

    kdt = u.KDTree(cellskel.vertices)
    nodecounts = Counter(kdt.query(centroids)[1])

    branchcount = Counter()
    for (node, count) in nodecounts.items():
        branchcount[branchlbl[node]] += count

    return {i: branchcount[i] for i in np.unique(branchlbl)}


def addmitovolumes(branchdf, ids, mitodf, mitotoskeldf):
    """
    Computes the mitochondrion volume for each compartment branch
    for a list of cell ids and adds these measurements to the branch
    data frame.
    """
    return compute_and_merge(compute_for_each_id, mitovolbybranch,
                             ids, mitodf, mitotoskeldf)


def mitovolbybranch(cellid, mitodf, mitotoskeldf):
    """
    Finds the mitochondria associated with each compartment branch
    of a single cell and computes the sum of their volumes.
    Returns a dictionary mapping branch ID to a total volume measurement.
    """
    cellskel, complbl = skel.read_skel_and_labels(cellid, refine=False)
    branchlbl = skel.branches(cellskel, complbl)

    mitovols = dict()
    branchids = np.unique(branchlbl)
    for i in branchids:
        nodeids = np.flatnonzero(branchlbl == i)
        branchcomplbl = complbl[nodeids[0]]
        engcomplbl = compartment.LABELMAP[branchcomplbl]

        # Collect the set of mitochondria that are associated with some
        # skeleton node in this branch and match the branch compartment
        # label
        mitoids = set(mitotoskeldf.mitoid[
                      (mitotoskeldf.cellid == cellid) &
                      (mitotoskeldf.nodeid.isin(nodeids))]
                      ).intersection(
                          mitodf.index[mitodf.compartment == engcomplbl]
                      )

        mitovols[i] = sum(mitodf.loc[mitoids].overlap_vx) * VOXELVOL

    return mitovols


def computedensities(branchdf):
    pass


def groupcells(branchdf, nucleusdf, nucleuslookup):
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("mitostats_filename", default=MITOSTATS_FILENAME)
    parser.add_argument("id_filename", default=ID_FILENAME)
    parser.add_argument("mitotoskel_filename", default=MITOTOSKEL_FILENAME)
    parser.add_argument("postsyn_filename", default=POSTSYN_FILENAME)
    parser.add_argument("basicbranchstats_filename",
                        default=BRANCHSTATS_FILENAME)
    parser.add_argument("overlap_filename", default=OVERLAP_FILENAME)
    parser.add_argument("fullbranchstats_filename")
    parser.add_argument("cellwisestats_filename")

    main(**vars(parser.parse_args()))
