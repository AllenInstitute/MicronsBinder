"""
'Basic' branch stats computation

Here, a branch refers to a connected component of a skeleton that possesses
the same compartment label. This means that branches can also refer to the
soma.

This script computes a data frame with an ID, volume, and surface area
measurement for each branch. The branch ID is the minimum skeleton node index
that is contained within the branch.
"""
import os
import argparse
import sys

import numpy as np
import pandas as pd
from meshparty import trimesh_io

sys.path.append(".")
from lib import compartment, skel, mesh, u


def main(outputfilename, idfilename,
         mitotoskelfilename=None, branchtype="full"):

    ids = u.readids(idfilename)
    mitotoskel = (compartment.readmitotoskel(mitotoskelfilename)
                  if mitotoskelfilename is not None else None)

    branchdf, errors = branchstats(ids, branchtype, mitotoskel)

    branchdf.to_csv(outputfilename)

    if len(errors) > 0:
        d, base = os.path.split(outputfilename)
        errfilename = os.path.join(d, f"errors_{base}")
        writeerrors(errors, errfilename)


def branchstats(cellids, branchtype="full", mitotoskel=None):
    """
    Computes the basic branch stats (volume, surface area) for a
    set of cell ids. Returns a pandas DataFrame, and a list of
    errors specified as a tuple (cellid, branchid).
    """
    errors = []
    # Specifying a dataframe by rows (i.e. records)
    records = list()
    for (i, cellid) in enumerate(cellids):
        try:
            cellskel, complbl = skel.read_skel_and_labels(cellid, refine=False)
            cellmesh = mesh.read_cell_mesh(cellid)
            cellmitotoskel = (mitotoskel[mitotoskel.cellid == cellid]
                              if mitotoskel is not None else None)

            # Label the skeleton by branch, and propagate those labels
            # to the mesh
            branchlbl, branchclbl = skelbranches(cellskel, complbl,
                                                 branchtype, cellmitotoskel)
            meshbranchlbl = meshbranches(cellmesh, cellskel, branchlbl)

            branchids = np.unique(branchlbl)
            for (j, branchid) in enumerate(branchids):
                print(f"Cell #{i+1}/{len(cellids)}"
                      f" - Branch #{j+1}/{len(branchids)}"
                      "                                  ",
                      end="\r")

                branchmesh = mesh.maskmesh(
                    cellmesh, meshbranchlbl == branchid)

                pl = skel.pathlength(cellskel, branchlbl, branchid)
                vol = meshvolume(branchmesh)
                surf = meshsurfacearea(branchmesh)

                records.append({
                    "cellid": cellid,
                    "branchid": branchid,
                    "complbl": branchclbl[branchid],
                    "pathlength": pl,
                    "volume": vol,
                    "surfacearea": surf
                })

        # Mostly for debugging, the final run never raised an exception
        except Exception as e:
            print("")
            print(f"ERROR on cell {(i, cellid)} branch {(j, branchid)}")
            print(e)
            errors.append((i, j))
            continue

    return pd.DataFrame.from_records(records), errors


def skelbranches(cellskel, complbl, branchtype="full", cellmitotoskel=None):
    """Labels the branches of a skeleton given the compartment labels"""
    somadists = compartment.compute_node_distance_to_soma(
                    cellskel, complbl, binary=False)

    if branchtype == "full":
        branchlbl, branchclbl = skel.branches(cellskel, complbl, somadists)
    elif branchtype == "branchpointsegments":
        branchlbl, branchclbl = skel.branchsegments(cellskel, complbl, somadists)
    elif branchtype == "distancebins":
        branchlbl, branchclbl = skel.distancebins(cellskel, complbl, somadists)
    elif branchtype == "distalbranches":
        branchlbl, branchclbl = skel.distalbranches(cellskel, complbl, somadists)
    elif branchtype == "coveredornot":
        coverednodes = list(set(cellmitotoskel.nodeid))
        branchlbl, branchclbl = skel.coveredornot(cellskel, complbl,
                                                   somadists, coverednodes)
    else:
        raise Exception(f"unknown branch type: {branchtype}")

    return branchlbl, branchclbl


def meshbranches(cellmesh, cellskel, branchlbl):
    """Propagates the branch labels from a skeleton to a mesh by proximity"""
    meshtoskel = mesh.node_associations(cellmesh, cellskel.vertices)
    matched = meshtoskel != -9999

    meshbranchlbl = np.zeros((len(cellmesh.vertices),), dtype=int)
    meshbranchlbl[matched] = branchlbl[meshtoskel[matched]]
    meshbranchlbl[~matched] = -9999

    return meshbranchlbl


def meshvolume(m):
    """Computes the volume of a mesh after some preprocessing"""
    # Seals any holes from cutting it from the rest of the mesh
    sealed = mesh.simple_mesh_seal(m)
    centered = mesh.center_mesh(sealed)

    return mesh.mesh_volume(centered).item()


def meshsurfacearea(m):
    """Computes the surface area of a mesh"""
    temp = trimesh_io.Mesh(m.vertices / 1000., m.faces)

    return temp.area


def writeerrors(errors, filename):
    with open(filename, 'w+') as f:
        for err in errors:
            f.write(f"{err}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("idfilename")
    parser.add_argument("outputfilename")
    parser.add_argument("--mitotoskelfilename", default=None)
    parser.add_argument("--branchtype", default="full")

    main(**vars(parser.parse_args()))
