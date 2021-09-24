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
import numpy as np
from meshparty import skeleton_io


SKELDIR = "../data/smoothed_skeletons_v185"
COMPLABELDIR = "../data/smoothed_skeletons_v185"

def read_neuron_skel(segid):
    filename = f"{SKELDIR}/{segid}_skeleton.h5"
    return skeleton_io.read_skeleton_h5(filename)


def read_neuron_complbls(segid):
    filename = f"{COMPLABELDIR}/{segid}_skeleton_label.npy"
    return np.load(filename)


def compartment_pathlength(skl, lbl, compartmentid, scale=1000.):
    """The path length covered by a compartment id"""
    pathlengths = np.linalg.norm(skl.vertices[skl.edges[:, 0]]
                                 - skl.vertices[skl.edges[:, 1]], axis=1)
    edgemask = (lbl[skl.edges] == compartmentid).min(axis=1)

    return pathlengths[edgemask].sum() / scale
