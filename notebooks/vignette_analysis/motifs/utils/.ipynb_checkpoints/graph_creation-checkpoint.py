import networkx as nx
import pandas as pd
import numpy as np

def edges_to_graph(df, multi=False):
    if multi:
        G = nx.MultiDiGraph()
    else:
        G = nx.DiGraph()
    for k, row in df.iterrows():
        G.add_node(row['pre'])
        G.add_node(row['post'])
        G.add_edge(row['pre'], row['post'])
    return G

def graph_to_sets(g):
    return set(g.nodes).copy(), set(g.edges).copy()

def set_to_graph(e):
    return nx.DiGraph(list(e))

def synapses_to_connections(df):
    return df[['pre','post']].drop_duplicates(subset=['pre','post'])

def remove_autapses(g):
    rg = g.copy()
    for e in g.edges:
        if e[0] == e[1]:
            rg.remove_edge(*e)
    return rg    

def format_synapses(df):
    df = df[['pre_pt_root_id','post_pt_root_id']]
    df2 = df.rename(columns={'pre_pt_root_id': 'pre', 'post_pt_root_id': 'post'})
    df2.reset_index(drop=True, inplace=True)
    return df2    

def get_recurrent_graph(g):
    rg = g.copy()
    while (np.array([g.out_degree(n) for n in g.nodes()]) == 0).sum() > 0:
        for n in g.nodes:
            if g.out_degree(n) == 0:
                rg.remove_node(n)
        g = rg.copy()
    return rg

def get_thresholded_graph(g, axls, th):
    rg = g.copy()
    for n in g.nodes:
        if axls[n] <= th:
            rg.remove_node(n)
    return rg

def get_central_graph(g):
    rg = g.copy()
    for n in g.nodes:
        if g.out_degree(n) == 0:
            rg.remove_node(n)
    g = rg.copy()
    return rg