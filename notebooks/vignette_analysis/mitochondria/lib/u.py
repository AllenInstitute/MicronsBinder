import re

import numpy as np
import pandas as pd
from scipy.spatial import KDTree  # aliasing this
import matplotlib.pyplot as plt


POS_REGEXP = re.compile(r"\[([0-9]*) +([0-9]*) +([0-9]*)\]")


def extract_centroids(df):
    # old df format
    # return np.array(list(zip(df.centroid_x, df.centroid_y, df.centroid_z)))
    return np.array(list(zip(df.ctr_pos_x_nm,
                             df.ctr_pos_y_nm,
                             df.ctr_pos_z_nm)))


def joinbybranchid(basedf, *otherdfs):
    for otherdf in otherdfs:
        basedf = pd.merge(basedf, otherdf,
                          left_on=["cellid", "branchid"],
                          right_on=["cellid", "branchid"])

    return basedf


def extract_centroids_dl(df, colname="ctr_pt_position"):
    if len(df) == 0:
        return np.empty((0, 3), dtype=np.int)

    return np.array(df[colname].tolist())


def extract_coords_dlcsv(df, colname="ctr_pt_position"):
    if len(df) == 0:
        return np.empty((0, 3), dtype=np.int)

    return np.array(list(map(parse_pos, df[colname])))


def parse_pos(pos, regexp=POS_REGEXP):
    """Extracts a tuple from the coordinates from DAL-derived csvs"""
    m = regexp.match(pos)
    return tuple(map(int, m.groups()))


def scale_to_nm(coord, voxel_res=[4, 4, 40]):
    return (coord[0]*voxel_res[0],
            coord[1]*voxel_res[1],
            coord[2]*voxel_res[2])


def scale_to_vx(coord, voxel_res=[4, 4, 40], asint=True):
    vx_coord = (coord[0]/voxel_res[0],
                coord[1]/voxel_res[1],
                coord[2]/voxel_res[2])

    if asint:
        vx_coord = tuple(map(int, vx_coord))

    return vx_coord


def logspace_bins(arr, n, eps=1e-10):
    return np.logspace(np.log10(arr.min())-eps,
                       np.log10(arr.max())+eps,
                       num=n)


def bins(arr, n, eps=1e-10):
    return np.linspace(arr.min()-eps,
                       arr.max()+eps,
                       num=n)


def read_node_lookup(filename, sep=';'):
    return read_lookup(filename, sep=sep)


def read_lookup(filename, sep=';'):
    lookup = dict()
    with open(filename) as f:
        for l in f:
            fields = l.strip().split(sep)
            lookup[eval(fields[0])] = eval(fields[1])

    return lookup


def readids(filename):
    with open(filename) as f:
        return list(map(int, f))


def writeids(ids, filename):
    with open(filename, "w+") as f:
        for i in ids:
            f.write(f"{i}\n")
