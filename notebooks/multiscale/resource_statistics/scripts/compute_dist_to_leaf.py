import sys
import numpy as np
import pandas as pd
import networkx as nx
import time
from tqdm import tqdm

sys.path.append(".")
from lib import skel
from lib import string2list as s2l


# Define volume properties
xres = yres = 3.58 # in nm/px
zres = 40


# Get all pyramidal and interneuron segment IDs
# These are the IDs from the soma subgraph and the consensus inhibitory cell table
pyrs = pd.read_csv('data/soma_ids/p100_pyr_soma_IDs_v185.csv',index_col=0)
inhs = pd.read_csv('data/soma_ids/p100_inh_soma_IDs_v185.csv',index_col=0)
inh_subtypes = pd.read_csv('data/inh_subtypes.csv', index_col=0)


def get_Euclidean_dist(p0,p1):
    distvec = [(p0[q]-p1[q])**2 for q in range(len(p0))]
    dist = np.sqrt(np.sum(distvec))
    return dist


def get_process_path_lengths(segid,centroid,is_pyr=False):
    
    skelcurr = skel.read_neuron_skel(segid)
    skelbls = skel.read_neuron_complbls(segid)

    verts = skelcurr.vertices
    edges = skelcurr.edges

    
    # Somatic nodes are classified so that nodes that are both
    # <15 um from the soma centroid and 
    # <100 nodes from the soma centroid
    # So we classify all nodes directly from their labels.
    # classify remaining nodes as axons or dendrites
    
    nodes = np.arange(len(skelbls))
    somanodes = np.flatnonzero(skelbls == 0)
    axnodes = np.flatnonzero(skelbls == 1)
    dendnodes = np.flatnonzero(np.isin(skelbls, [2, 3, 4]))
    
    if centroid == None:
        centroid = np.mean(verts[somanodes], axis=0)
        
    dists_centroid = [get_Euclidean_dist(verts[q],centroid) for q in nodes]
    somaroot = nodes[np.argmin(dists_centroid)]
    
    
    # print('numbers of nodes',len(axnodes),len(dendnodes))
    axedges = []
    dendedges = []
    
    axedges = edges[(skelbls[edges] == 1).max(1)]
    dendedges = edges[(np.isin(skelbls[edges], [2, 3, 4])).max(1)]

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
    
    dendgraph.add_nodes_from(dendnodes)
    dend_weighted_edges = []
    for e in dendedges:
        wcurr = get_Euclidean_dist(verts[e[0]],verts[e[1]])
        dend_weighted_edges.append((e[0],e[1],wcurr))
    dendgraph.add_weighted_edges_from(dend_weighted_edges)
    
    # Break the axonal and dendritic graphs into connected subgraph components
    asgs = list(axgraph.subgraph(c)
                for c in nx.connected_components(axgraph))
    dsgs = list(dendgraph.subgraph(c)
                for c in nx.connected_components(dendgraph))



    # Define containers for axon and dendritic path lengths
    axdists = []
    denddists = []
    if is_pyr:
        apical_denddists = []
        basal_denddists = []
        ambig_denddists = []
    
    # For each subgraph of a pyramidal cell (where skeletons are broken),
    # find the node closest to the soma root node and add an edge between that
    # node and the soma root in the full graph
        # in this case, just find the leaf nodes in each connected subgraph
        # and use dijkstra to find path lengths
        for asg in asgs:
            aleaves = [q for q in asg.nodes if asg.degree(q)==1]
            # the source node is the one among these closest to the somaroot
            asourceid = np.argmin([get_Euclidean_dist(verts[somaroot],verts[q]) for q in aleaves])
            asource = aleaves[asourceid]
            atargs = [aleaves[q] for q in range(len(aleaves)) if q != asourceid]
            # print(asource,atargs) # debugging
            for at in atargs:
                acurr = nx.shortest_path_length(asg,source=asource,target=at,weight='weight')/1000.0 # convert to um
                axdists.append(acurr)
            
        for dsg in dsgs:
            dleaves = [q for q in dsg.nodes if dsg.degree(q)==1]
            dsourceid = np.argmin([get_Euclidean_dist(verts[somaroot],verts[q]) for q in dleaves])
            dsource = dleaves[dsourceid]
            dtargs = [dleaves[q] for q in range(len(dleaves)) if q != dsourceid]
            #print(dsource,dtargs) # debugging
            for dt in dtargs:
                dcurr = nx.shortest_path_length(dsg,source=dsource,target=dt,weight='weight')/1000.0 # convert to um
                denddists.append(dcurr)
                if skelbls[dt] == 2:
                    basal_denddists.append(dcurr)
                elif skelbls[dt] == 3:
                    apical_denddists.append(dcurr)
                elif skelbls[dt] == 4:
                    ambig_denddists.append(dcurr)
            
    else:       
        # Find axonal and dendritic leaf nodes
        for asg in asgs:
            aleaves = [q for q in asg.nodes if asg.degree(q)==1]
            asourceid = np.argmin([get_Euclidean_dist(verts[somaroot],verts[q]) for q in aleaves])
            asource = aleaves[asourceid]
            atargs = [aleaves[q] for q in range(len(aleaves)) if q != asourceid]
            for at in atargs:
                acurr = nx.shortest_path_length(axgraph,source=asource,target=at,weight='weight')/1000.0 # convert to um
                axdists.append(acurr)
                
        for dsg in dsgs:
            dleaves = [q for q in dsg.nodes if dsg.degree(q)==1]
            dsourceid = np.argmin([get_Euclidean_dist(verts[somaroot],verts[q]) for q in dleaves])
            dsource = dleaves[dsourceid]
            dtargs = [dleaves[q] for q in range(len(dleaves)) if q != dsourceid]
            for dt in dtargs:
                dcurr = nx.shortest_path_length(dendgraph,source=dsource,target=dt,weight='weight')/1000.0 # um
                denddists.append(dcurr)
                
    if is_pyr:
        return axdists,denddists,basal_denddists,apical_denddists,ambig_denddists
    else:
        return axdists,denddists


if __name__ == "__main__":
    # Process pyramidal neurons
    start_time = time.time()
    pyr_neurite_info = pd.DataFrame(index=pyrs.index,columns=['axon_lengths_um','apical_dend_lengths_um','basal_dend_lengths_um','ambiguous_dend_lengths_um'])
    for idx in tqdm(pyrs.index):
        rc = pyrs.loc[idx]
        rtid = rc['pt_root_id']
        centroidraw = s2l.string2list(rc['pt_position'])
        centroid = [centroidraw[0]*xres, centroidraw[1]*yres, centroidraw[2]*zres] # convert to nm
        neuritelengths = get_process_path_lengths(rtid,centroid,is_pyr=True)
        pyr_neurite_info.at[idx,'axon_lengths_um'] = neuritelengths[0]
        pyr_neurite_info.at[idx,'basal_dend_lengths_um'] = neuritelengths[2]
        pyr_neurite_info.at[idx,'apical_dend_lengths_um'] = neuritelengths[3]
        pyr_neurite_info.at[idx,'ambiguous_dend_lengths_um'] = neuritelengths[4]
    end_time = time.time()

    pyrs_w_neurites = pyrs.join(pyr_neurite_info)
    pyrs_w_neurites.to_csv('data/pyr_dist_to_leaves.csv',index=True)

    print('Time for excitatory neurons: {0:.02} seconds.'.format(end_time - start_time))


    # Process inhibitory interneurons
    start_time_2 = time.time()
    inh_neurite_info = pd.DataFrame(index=inhs.index,columns=['axon_lengths_um','dendrite_lengths_um'])
    for idx in tqdm(inhs.index):
        rc = inhs.loc[idx]
        rtid = rc['pt_root_id']
        axlengths,dendlengths = get_process_path_lengths(rtid,centroid=None,is_pyr=False)
        inh_neurite_info.at[idx,'axon_lengths_um'] = axlengths
        inh_neurite_info.at[idx,'dendrite_lengths_um'] = dendlengths
    end_time_2 = time.time()

    inhs_w_neurites = inhs.join(inh_neurite_info)
    inhs_w_neurites = pd.merge(inhs_w_neurites, inh_subtypes,
                               left_on="pt_root_id", right_on="pt_root_id")
    inhs_w_neurites.to_csv('data/inh_dist_to_leaves.csv',index=True)

    print('Time for inhibitory neurons: {0:.02} seconds.'.format(end_time_2 - start_time_2))
