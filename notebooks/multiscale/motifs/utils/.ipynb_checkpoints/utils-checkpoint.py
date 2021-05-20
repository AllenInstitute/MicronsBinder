import ast
import numpy as np
import cloudvolume as cv
import json
from meshparty import skeleton_io

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


SMOOTHED_SKELDIR = "data/smoothed_skels"
SMOOTHED_LABELDIR = "data/smoothed_skels"

def read_smoothed_neuron_skel(segid):
    filename = f"{SMOOTHED_SKELDIR}/{segid}_skeleton.h5"
    return skeleton_io.read_skeleton_h5(filename)


def read_smoothed_neuron_complbls(segid):
    filename = f"{SMOOTHED_LABELDIR}/{segid}_skeleton_label.npy"
    return np.load(filename)

def skel_length(segid):
    return CVOL.skeleton.get(segid).cable_length()

def string2array(array_string):
    array_string = ','.join(array_string.replace('[ ', '[').split())
    return np.array(ast.literal_eval(array_string))