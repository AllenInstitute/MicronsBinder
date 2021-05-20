import numpy as np
import networkx as nx
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


SMOOTHED_SKELDIR = "data/smoothed_skeletons_v185"
SMOOTHED_LABELDIR = "data/smoothed_skeletons_v185"

def read_smoothed_neuron_skel(segid):
    filename = f"{SMOOTHED_SKELDIR}/{segid}_skeleton.h5"
    print(filename)
    return skeleton_io.read_skeleton_h5(filename)


def read_smoothed_neuron_complbls(segid):
    filename = f"{SMOOTHED_LABELDIR}/{segid}_skeleton_label.npy"
    return np.load(filename)


def skel_length(segid):
    return CVOL.skeleton.get(segid).cable_length()


def get_Euclidean_dist(p0,p1):
    distvec = [(p0[q]-p1[q])**2 for q in range(len(p0))]
    dist = np.sqrt(np.sum(distvec))
    return dist

def get_total_lengths(segid,centroid):
    
    # With smoothed skeletons you should be able to access
    # both pyramidal and inhibitory data in the same way
    skelcurr = skel.read_smoothed_neuron_skel(segid)
    skelbls = skel.read_smoothed_neuron_complbls(segid)
    
        
    verts = skelcurr.vertices
    edges = skelcurr.edges

    
    # Somatic nodes are classified so that nodes that are both
    # <15 um from the soma centroid and 
    # <100 nodes from the soma centroid
    # So we classify all nodes directly from their labels.
    # classify remaining nodes as axons or dendrites
    
    nodes = [q for q in range(len(skelbls))]
    somanodes = [q for q in nodes if skelbls[q] == 0]
    axnodes = [q for q in nodes if skelbls[q] == 1]
    dendnodes = [q for q in nodes if skelbls[q] in [2,3,4]]
    
    if centroid is None:
        centroid_x = np.mean([verts[q][0] for q in somanodes])
        centroid_y = np.mean([verts[q][1] for q in somanodes])
        centroid_z = np.mean([verts[q][2] for q in somanodes])
        centroid = [centroid_x,centroid_y,centroid_z]
        print(centroid)
        
    dists_centroid = [get_Euclidean_dist(verts[q],centroid) for q in nodes]
    somaroot = nodes[np.argmin(dists_centroid)]
    
    
    # print('numbers of nodes',len(axnodes),len(dendnodes))
    axedges = []
    dendedges = []
    
    for e in edges:
        # explicitly remove any edges connecting nodes outside the soma threshold radius
        # to nodes inside the radius
        isinax = [q for q in e if q in axnodes]
        isindend = [q for q in e if q in dendnodes]
        if len(isinax) > 0:
            axedges.append(e)
        if len(isindend) > 0:
            dendedges.append(e)
            
    # print('number of edges',len(axedges),len(dendedges))

    # Make separate axonal and dendritic graphs
    axgraph = nx.Graph()
    dendgraph = nx.Graph()
    
    axgraph.add_nodes_from(axnodes)
    ax_weighted_edges = []
    for e in axedges:
        wcurr = get_Euclidean_dist(verts[e[0]],verts[e[1]])
        ax_weighted_edges.append((e[0],e[1],wcurr))
    axgraph.add_weighted_edges_from(ax_weighted_edges)
    ax_total_len = axgraph.size(weight="weight") / 1000.0
    
    dendgraph.add_nodes_from(dendnodes)
    dend_weighted_edges = []
    for e in dendedges:
        wcurr = get_Euclidean_dist(verts[e[0]],verts[e[1]])
        dend_weighted_edges.append((e[0],e[1],wcurr))
    dendgraph.add_weighted_edges_from(dend_weighted_edges)
    den_total_len = dendgraph.size(weight="weight") / 1000.0
            
    return ax_total_len, den_total_len