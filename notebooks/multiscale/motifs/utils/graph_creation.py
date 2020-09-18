import networkx as nx
import pandas as pd

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
    df = df[['pre_pt_root_id','post_pt_root_id','size']]
    df2 = df.rename(columns={'pre_pt_root_id': 'pre', 'post_pt_root_id': 'post', 'ctr_pt_position':'position'})
    df2.reset_index(drop=True, inplace=True)
    return df2    