import argparse

from tqdm import tqdm
import numpy as np
import pandas as pd

from lib import skel, compartment, u


def main(cellidfilename, outputfilename, branchtype="full"):
    ids = u.readids(cellidfilename)

    celldfs = list()
    for i in tqdm(ids):
        celldfs.append(branchsegstats(i, branchtype))

    fulldf = pd.concat(celldfs, ignore_index=True)
    fulldf.to_csv(outputfilename)


def branchsegstats(cellid, branchtype):
    skl = skel.read_neuron_skel(cellid)
    lbl = skel.read_neuron_complbls(cellid)

    somadists = compartment.compute_node_distance_to_soma(
                    skl, lbl, binary=False)
    with np.errstate(divide='ignore',invalid='ignore'):
        if branchtype == "full":
            grouplbl, _ = skel.branches(skl, lbl, somadists)
        elif branchtype == "distancebins":
            grouplbl, _ = skel.distancebins(skl, lbl, somadists)
        elif branchtype == "distalbranches":
            grouplbl, _ = skel.distalbranches(skl, lbl, somadists)
        else:
            raise ValueError(f"unknown branch type: {branchtype}")

    groupids = np.unique(grouplbl)
    groupdists = [somadists[grouplbl == i] / 1000. for i in groupids]
    minsomadists = [np.min(d) for d in groupdists]
    meansomadists = [np.mean(d) for d in groupdists]
    maxsomadists = [np.max(d) for d in groupdists]

    leafdists = skel.path_length_to_leaves(skl)
    minleafdists = [np.min(leafdists[grouplbl == i]) / 1000. for i in groupids]

    return pd.DataFrame({
               "cellid": cellid,
               "branchid": groupids,
               "mindisttosoma": minsomadists,
               "meandisttosoma": meansomadists,
               "maxdisttosoma": maxsomadists,
               "mindisttoleaf": minleafdists,
               })


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("cellidfilename")
    parser.add_argument("outputfilename")
    parser.add_argument("--branchtype", default="full")

    main(**vars(parser.parse_args()))
