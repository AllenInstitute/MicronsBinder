import itertools

import numpy as np

from . import skel, mesh, u


def somavol_by_nearby_meshvol(cellid, distthresh=15_000):
    """
    Measuring soma volume by measuring the volume of the mesh within
    a threshold distance of the soma pt
    """

    segskel = skel.read_neuron_skel(cellid)
    segmesh = mesh.read_neuron_mesh(cellid)

    soma_pt = segskel.vertices[segskel.root]

    within_dist = np.linalg.norm(
                      segmesh.vertices - soma_pt, axis=1) < distthresh

    somamesh = mesh.mesh_by_inds(segmesh, np.nonzero(within_dist)[0])
    sealed = mesh.simple_mesh_seal(somamesh)
    centered = mesh.center_mesh(sealed)
    scaled = mesh.scale_mesh(centered, [4, 4, 40], [3.58, 3.58, 40])

    return mesh.mesh_volume(scaled), scaled


def somatic_ids(mitodf, cellid, distthresh=15_000):
    subdf = mitodf[mitodf.v185_rt_id == cellid]

    centroids = u.extract_centroids(subdf)
    segskel = skel.read_neuron_skel(cellid)
    soma_pt = segskel.vertices[segskel.root]

    within_thresh = np.linalg.norm(centroids - soma_pt, axis=1) < distthresh

    return subdf.index[within_thresh]


def all_somatic_rows(mitodf, ids=None):

    if ids is None:
        ids = set(mitodf.v185_rt_id)

    somatic_inds = list()
    for (i, cellid) in enumerate(ids):
        print(i, end="\r")
        somatic_inds.append(somatic_ids(mitodf, cellid))
    print("")

    all_inds = itertools.chain.from_iterable(somatic_inds)

    return mitodf.loc[all_inds]
