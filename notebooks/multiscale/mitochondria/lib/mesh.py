import os
import operator
from collections import Counter

import numpy as np
from scipy import sparse as sp
from scipy.sparse import csgraph
from meshparty import trimesh_io
import networkx as nx
import h5py


MESHMETA = trimesh_io.MeshMeta()
MITOMESHDIR = ""
CELLMESHDIR = "data/neuronmeshes"


def download_meshes(segids, overwrite=False, parallel=1,
                    meshdir=MITOMESHDIR, cvpath=None, **kwargs):
    trimesh_io.download_meshes(
        segids, meshdir, cvpath,
        overwrite=overwrite, n_threads=parallel, **kwargs)


def read_mito_mesh(segid):
    return read_mesh(f"{MITOMESHDIR}/{segid}.h5")


def read_cell_mesh(segid):
    filename = f"{CELLMESHDIR}/{segid}_mesh.h5"

    return read_mesh(filename)


def read_neuron_mesh(segid):
    return read_cell_mesh(segid)


def read_mesh_id(segid, meshdir=MITOMESHDIR):
    return read_mesh(f"{meshdir}/{segid}.h5")


def read_mesh(filename, scale=True):
    mesh = MESHMETA.mesh(filename)

    if scale:
        mesh.voxel_scaling = [0.895, 0.895, 1]

    return mesh


def read_mesh_h5(filename):
    with h5py.File(filename) as f:
        vertices = f["vertices"][()]
        faces = f["faces"][()].reshape(-1, 3)
        link_edges = f["link_edges"][()]

    return vertices, faces


def write_cv_mesh(mesh, filename):
    vertices = mesh.vertices
    faces = mesh.faces
    link_edges = np.empty((0, 2), dtype=np.uint8)

    with h5py.File(filename, "w") as f:
        f.create_dataset("vertices", data=vertices)
        f.create_dataset("faces", data=faces)
        f.create_dataset("link_edges", data=link_edges)


def mesh_length(graph):
    first_node = next(iter(graph.nodes))
    initial_dists = nx.single_source_dijkstra_path_length(graph, first_node)
    furthest_node = max(initial_dists.items(), key=operator.itemgetter(1))[0]

    lengths = nx.single_source_dijkstra_path_length(graph, furthest_node)
    other_end, length = max(lengths.items(), key=operator.itemgetter(1))
    return length, (furthest_node, other_end)


def simple_mesh_seal(mesh):
    """Seals a mesh by adding one vertex at the centroid of each boundary"""
    components = find_mesh_boundary_components(mesh)
    for comp in components:
        new_pt = np.mean(mesh.vertices[comp], axis=0)
        mesh = attach_new_pt(mesh, new_pt, boundary=comp)

    return mesh


def find_mesh_boundary_vertices(mesh):
    boundary_edges = find_mesh_boundary_edges(mesh)
    return list(set(itertools.chain.from_iterable(boundary_edges)))


def find_mesh_boundary_components(mesh):
    boundary_edges = find_mesh_boundary_edges(mesh)

    return edge_components(boundary_edges)


def edge_components(edges):
    assert len(edges) > 0, "no edges passed"
    rows, cols = zip(*edges)

    ids = list(set(rows + cols))
    idmap = {v: i for (i, v) in enumerate(ids)}
    rows = [idmap[v] for v in rows]
    cols = [idmap[v] for v in cols]

    vals = np.ones((len(rows),), dtype=np.uint8)
    num_ids = max(max(rows), max(cols)) + 1
    g = sp.coo_matrix((vals, (rows, cols)), shape=(num_ids, num_ids)).tocsr()

    num_comps, labels = csgraph.connected_components(
                            g, directed=False, return_labels=True)
    graph_comps = [np.nonzero(labels == i)[0] for i in range(num_comps)]

    return [[ids[v] for v in comp] for comp in graph_comps]


def find_mesh_boundary_edges(mesh):
    faces = mesh.faces.reshape((-1, 3))
    edges = [tuple(sorted(edge)) for edge in
             np.vstack((faces[:, :2],
                        faces[:, 1:],
                        faces[:, [2, 0]]))]

    counter = Counter(edges)
    return [edge for (edge, count) in counter.items() if count == 1]


def dfs_order_vertices(edges):
    rows, cols = zip(*edges)

    ids = list(set(rows + cols))
    idmap = {v: i for (i, v) in enumerate(ids)}
    rows = [idmap[v] for v in rows]
    cols = [idmap[v] for v in cols]

    vals = np.ones((len(rows),), dtype=np.uint8)
    num_ids = max(max(rows), max(cols)) + 1
    g = sp.coo_matrix((vals, (rows, cols)), shape=(num_ids, num_ids)).tocsr()

    return [ids[i] for i in csgraph.depth_first_order(g, 0, directed=False)[0]]


def flip_if_mostly_inwards(faces, vertices, new_pt):
    bad_faces = find_bad_face_normals(faces, vertices, new_pt)
    if sum(bad_faces) > len(faces) // 2:
        faces = faces[:, [1, 0, 2]]
    return faces


def fix_bad_face_normals(faces, vertices, new_pt):
    bad_faces = find_bad_face_normals(faces, vertices, new_pt)
    faces = faces.copy()
    faces[bad_faces] = faces[bad_faces][:, [1, 0, 2]]

    return faces


def find_bad_face_normals(faces, vertices, new_pt):
    centroid = np.mean(vertices, axis=0)
    normals = compute_normals(faces, vertices)

    return (normals @ (new_pt - centroid)) < 0


def compute_normals(faces, vertices):
    v0 = vertices[faces[:, 0]]
    v1 = vertices[faces[:, 1]]
    v2 = vertices[faces[:, 2]]

    return np.cross(v1-v0, v2-v0)


def mesh_volume(mesh, scale=1000):

    vs = mesh.vertices / scale
    fs = mesh.faces.reshape((-1, 3))

    v0 = vs[fs[:, 0]]
    v1 = vs[fs[:, 1]]
    v2 = vs[fs[:, 2]]

    return np.abs(np.sum(np.cross(v0, v1) * v2 / 6.))


def center_mesh(mesh):
    centroid = np.mean(mesh.vertices, axis=0)
    return trimesh_io.Mesh(vertices=mesh.vertices-centroid,
                           faces=mesh.faces)


def scale_mesh(mesh, old_res, new_res):
    verts = (mesh.vertices / old_res) * new_res
    return trimesh_io.Mesh(vertices=verts, faces=mesh.faces)


def mesh_by_inds(mesh, inds):
    assert len(inds) > 0, "empty inds"
    new_verts = mesh.vertices[inds]

    ind_map = np.empty((max(inds)+1,), dtype=inds.dtype)
    ind_map[inds] = np.arange(len(inds))
    explicit_faces = mesh.faces.reshape((-1, 3))
    face_inds = np.all(np.isin(explicit_faces, inds), axis=1)
    new_faces = ind_map[explicit_faces[face_inds]].ravel()

    return trimesh_io.Mesh(vertices=new_verts, faces=new_faces)

def attach_new_pt(mesh, new_pt, boundary=None):
    new_verts = np.vstack((mesh.vertices, new_pt))
    new_i = len(mesh.vertices)

    boundary_edges = find_mesh_boundary_edges(mesh)
    if boundary is not None:
        bset = set(boundary)
        boundary_edges = [edge for edge in boundary_edges
                          if (edge[0] in bset) and (edge[1] in bset)]

    new_faces = np.array(assemble_consistent_faces(boundary_edges, new_i))
    new_faces = flip_if_mostly_inwards(new_faces, new_verts, new_pt)
    all_new_faces = np.hstack((mesh.faces, new_faces.ravel()))

    return trimesh_io.Mesh(vertices=new_verts, faces=all_new_faces)


def assemble_consistent_faces(edges, new_i):
    vertex_order = dfs_order_vertices(edges)
    return [[vertex_order[i-1], vertex_order[i], new_i]
            for i in range(len(vertex_order))]


def dfs_order_vertices(edges):
    rows, cols = zip(*edges)

    ids = list(set(rows + cols))
    idmap = {v: i for (i, v) in enumerate(ids)}
    rows = [idmap[v] for v in rows]
    cols = [idmap[v] for v in cols]

    vals = np.ones((len(rows),), dtype=np.uint8)
    num_ids = max(max(rows), max(cols)) + 1
    g = sp.coo_matrix((vals, (rows, cols)), shape=(num_ids, num_ids)).tocsr()

    return [ids[i] for i in csgraph.depth_first_order(g, 0, directed=False)[0]]


def flip_if_mostly_inwards(faces, vertices, new_pt):
    bad_faces = find_bad_face_normals(faces, vertices, new_pt)
    if sum(bad_faces) > len(faces) // 2:
        faces = faces[:, [1, 0, 2]]
    return faces


def fix_bad_face_normals(faces, vertices, new_pt):
    bad_faces = find_bad_face_normals(faces, vertices, new_pt)
    faces = faces.copy()
    faces[bad_faces] = faces[bad_faces][:, [1, 0, 2]]

    return faces


def find_bad_face_normals(faces, vertices, new_pt):
    centroid = np.mean(vertices, axis=0)
    normals = compute_normals(faces, vertices)

    return (normals @ (new_pt - centroid)) < 0


def compute_normals(faces, vertices):
    v0 = vertices[faces[:, 0]]
    v1 = vertices[faces[:, 1]]
    v2 = vertices[faces[:, 2]]

    return np.cross(v1-v0, v2-v0)
