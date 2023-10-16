import argparse

from tqdm import tqdm
import numpy as np
import pandas as pd

from lib import skel, compartment, u


def main(cellidfilename, mitotoskelfilename,
         outputfilename, branchtype="full"):
    ids = u.readids(cellidfilename)
    mitotoskeldf = compartment.readmitotoskel(mitotoskelfilename)

    celldfs = list()
    for i in tqdm(ids):
        celldfs.append(groupstats(i, mitotoskeldf, branchtype))

    fulldf = pd.concat(celldfs, ignore_index=True)
    fulldf.to_csv(outputfilename)


def groupstats(cellid, mitotoskeldf, branchtype="full"):
    skl = skel.read_neuron_skel(cellid)
    lbl = skel.read_neuron_complbls(cellid)

    coverednodes = list(set(mitotoskeldf.nodeid[
                                mitotoskeldf.cellid == cellid]))

    somadists = compartment.compute_node_distance_to_soma(
                    skl, lbl, binary=False)

    with np.errstate(invalid="ignore", divide="ignore"):
        if branchtype == "full":
            grouplbl, _ = skel.branches(skl, lbl, somadists)
        elif branchtype == "distancebins":
            grouplbl, _ = skel.distancebins(skl, lbl, somadists)
        elif branchtype == "distalbranches":
            grouplbl, _ = skel.distalbranches(skl, lbl, somadists)
        elif branchtype == "coveredornot":
            grouplbl, _ = skel.coveredornot(skl, lbl, somadists, coverednodes)
        elif branchtype == "cellwise":
            grouplbl = lbl
        else:
            raise ValueError(f"unknown branch type: {branchtype}")

    groupids = np.unique(grouplbl)
    diameters = np.array(skl.vertex_properties["rs"])
    groupdiams = [diameters[grouplbl == i] for i in groupids]

    mindiams = [np.min(d) for d in groupdiams]
    meandiams = [np.mean(d) for d in groupdiams]
    mediandiams = [np.median(d) for d in groupdiams]
    maxdiams = [np.max(d) for d in groupdiams]

    return pd.DataFrame({
               "cellid": cellid, "branchid": groupids,
               "mindiam": mindiams, "meandiam": meandiams,
               "mediandiam": mediandiams, "maxdiam": maxdiams,
               })


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("cellidfilename")
    parser.add_argument("mitotoskelfilename")
    parser.add_argument("outputfilename")
    parser.add_argument("--branchtype", default="full")

    main(**vars(parser.parse_args()))
