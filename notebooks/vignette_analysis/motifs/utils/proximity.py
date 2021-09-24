import pandas as pd
import numpy as np
import networkx as nx

def create_proximity_model(g, locs, avg_th=10):
    # get connection vs dist stats
    conn_p_dist = []
    for source_node in g.nodes():
        for target_node in g.nodes():
            if source_node == target_node: continue
            distance = np.sqrt(((locs[target_node] / 1000 - locs[source_node] / 1000)**2).sum())
            conn = 1 if (source_node, target_node) in g.edges() else 0
            conn_p_dist.append([distance, conn])
    conn_p_dist = np.array(conn_p_dist)

    edge_order = conn_p_dist[:,0].argsort()
    conn_p_dist = conn_p_dist[edge_order]
    edge_order_rindex = edge_order.argsort()

    # estimate connection probablity from moving average
    # window size = avg_th
    avg_conn_p_dist = []
    lind = 0
    rind = 0
    max_ind = len(conn_p_dist)
    for i in range(len(conn_p_dist)):
        while conn_p_dist[lind, 0] < conn_p_dist[i, 0] - avg_th:
            lind += 1
        while rind < max_ind-1 and conn_p_dist[rind, 0] < conn_p_dist[i, 0] + avg_th:
            rind += 1
        inds = np.arange(lind, rind)
        dist = conn_p_dist[:,0][inds].mean()
        conn = conn_p_dist[:,1][inds].mean()
        counts = len(inds)
        avg_conn_p_dist.append([dist, conn, counts])

    avg_conn_p_dist = np.array(avg_conn_p_dist)


    prox = pd.DataFrame({
        "soma_dist":avg_conn_p_dist[:,0],
        "p_connect":avg_conn_p_dist[:,1],
        "p_connect_std":avg_conn_p_dist[:,2],
        "edge_order":edge_order,
        })
    
    return prox


def prox_sample_degs(g, ps, edge_order, selfloop=True):
    cp_indegs = []
    cp_outdegs = []
    for i in range(100):
        if selfloop:
            # ER model with selfloop
            edges = g.edges()
            inds = np.random.rand(len(edges)) < ps
            edges = np.array(edges)[edge_order][inds]
            g_sampled = nx.DiGraph()
            g_sampled.add_nodes_from(g.nodes())
            g_sampled.add_edges_from(edges)
        else:
            g_sampled = nx.erdos_renyi_graph(len(g.nodes()), p, directed=True)
            
        cp_indegs.append([g_sampled.in_degree(n) for n in g_sampled.nodes()])
        cp_outdegs.append([g_sampled.out_degree(n) for n in g_sampled.nodes()])
    cp_indegs = np.array(cp_indegs).reshape(-1)
    cp_outdegs = np.array(cp_outdegs).reshape(-1)
    return cp_indegs, cp_outdegs