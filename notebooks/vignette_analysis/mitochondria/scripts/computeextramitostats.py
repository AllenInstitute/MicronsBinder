import sys
import math
import argparse

import numpy as np
import pandas as pd
import trimesh

sys.path.append(".")
from lib import mesh, compartment


VXVOL = (3.58/1000. * 3.58/1000. * 40/1000.)


def main(mitodffilename, mitotoskelfilename, outputfilename):
    print("Reading input data")
    mitodf = pd.read_csv(mitodffilename, index_col=0)
    mitotoskel = compartment.readmitotoskel(mitotoskelfilename)

    print("Computing surface area")
    mitodf["surface_area"] = computesurfaceareas(mitodf)

    print("Computing complexity index")
    mitodf["complexityindex"] = computeMCIs(mitodf)

    print("Labeling compartments")
    mitodf["compartment"] = labelcompartments(mitodf, mitotoskel)

    mitodf.to_csv(outputfilename)


def labelcompartments(mitodf, mitotoskel, thresh=1_000):
    """Propagating skeleton compartment labels to mitochondria

    Most mitochondria are labeled by majority vote across their associated
    skeleton nodes. Some nodes "close" to the soma in path length are
    labeled as unknown as these can bias different measurements.
    """
    ids, majlbls_ = compartment.majority_vote_label(mitotoskel)
    majlbls = [compartment.ENGLISHLABELS[lbl] for lbl in majlbls_]
    labellookup = dict(zip(ids, majlbls))

    unknownlabel = compartment.ENGLISHLABELS[4]
    #mindist = mitotoskel.groupby("mitoid")["nodelbl", "nodedist"].min()
    #idsnearsoma = set(mindist.index[(mindist.nodedist < thresh) &
    #                                (mindist.nodelbl != 0)].tolist())
    #for i in idsnearsoma:
    #    labellookup[i] = unknownlabel

    return [labellookup.get(i, unknownlabel) for i in mitodf.index]


def computesurfaceareas(mitodf):
    """Computing the surface area of each mitochondrion mesh

    Largely just relies on trimesh.Trimesh.area
    """
    surfaceareas = list()
    for (i, mitoid) in enumerate(mitodf.index):
        print(f"#{i}/{len(mitodf)} - {mitoid}           ", end="\r")
        try:
            msh = mesh.read_mito_mesh(mitoid)
            temp = trimesh.Trimesh(msh.vertices/1000., msh.faces)
            surfaceareas.append(temp.area)

        except AssertionError:
            surfaceareas.append(np.nan)

        except Exception as e:
            print("unexpected exception on #{i}-{mitoid}")
            print(e)
            surfaceareas.append(np.nan)
    print("")

    return surfaceareas


def computeMCIs(mitodf):
    """Computes the mitochondria complexity index"""
    numerator = mitodf["surface_area"] ** 3
    vol = mitodf["mito_vx"] * VXVOL
    denominator = 16 * math.pi**2 * vol**2

    return numerator / denominator


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("mitodffilename")
    parser.add_argument("mitotoskelfilename")
    parser.add_argument("outputfilename")

    main(**vars(parser.parse_args()))
