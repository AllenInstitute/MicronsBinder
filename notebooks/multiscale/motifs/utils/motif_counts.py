import numpy as np
import random
import networkx as nx
from time import time
import pandas as pd
from .graph_creation import *
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from itertools import repeat
from functools import partial
from tqdm.auto import tqdm, trange

### swtich-and-hold algorithms for sampling from configuration model
def switch_and_hold(edges, niters=1000, verbose=False):
    """Iterative implementation of the configuration model
    
    Args:
        edges: set of vertex pairs
        niters: int for number of iterations
    """
    holds = 0
    
    def get_random_edge_pair(edges):
        return random.sample(edges, 2)

    def switch_allowed(edges, edge_pair):
        u, v = edge_pair
        if (u[0] == v[1]) or (v[0] == u[1]):
            return False
        x, y = (u[0], v[1]), (v[0], u[1])
        if (x in edges) or (y in edges):
            return False
        return True

    def make_switch(edges, edge_pair):
        u, v = edge_pair
        x, y = (u[0], v[1]), (v[0], u[1])
        edges.remove(u)
        edges.remove(v)
        edges.add(x)
        edges.add(y)    
    
    start = time()
    for i in range(niters):
        ep = get_random_edge_pair(edges)
        if switch_allowed(edges, ep):
            make_switch(edges, ep)
        else:
            holds += 1
    end = time()
    if verbose:
        print('holds: {}, {:0.2f}'.format(holds, holds / niters))
        print('time: {:0.2f} s'.format(end-start))
    return edges



def switch_and_hold_GE(edges, niters=1000, verbose=False):
    """configuration model while hoding 2 neuron motifs
    
    Args:
        edges: set of vertex pairs
        niters: int for number of iterations
    """
    holds = 0
    
    def get_random_edge_pair(edges):
        return random.sample(edges, 2)

    def switch_allowed(edges, edge_pair):
        u, v = edge_pair
        if (u[0] == v[1]) or (v[0] == u[1]):
            return False
        x, y = (u[0], v[1]), (v[0], u[1])
        if (x in edges) or (y in edges):
            return False
        return True

    def make_switch(edges, edge_pair):
        
        # calculate 2-node motif invariant
        u, v = edge_pair
        ur, vr = (u[1], u[0]), (v[1], v[0])
        xr, yr = (u[1], v[0]), (v[1], u[0])
        ext = lambda e: 1 if e in edges else 0
        inv = ext(ur) + ext(vr) - ext(xr) - ext(yr)
        
        u, v = edge_pair
        x, y = (u[0], v[1]), (v[0], u[1])
        edges.remove(u)
        edges.remove(v)
        edges.add(x)
        edges.add(y)
        
        return inv
    
    start = time()
    inv = 0
    for i in range(niters):
        ep = get_random_edge_pair(edges)
        if switch_allowed(edges, ep):
            inv += make_switch(edges, ep)
        else:
            holds += 1
    
    extra_inters = 0
    while inv != 0:
        extra_inters += 1
        ep = get_random_edge_pair(edges)
        if switch_allowed(edges, ep):
            inv += make_switch(edges, ep)
        else:
            holds += 1
        if extra_inters > 5 * niters:
            print("constraints can not be satisfied..")
            break
            
    end = time()
    if verbose:
        print('holds: {}, {:0.2f}'.format(holds, holds / (niters+extra_inters)))
        print('time: {:0.2f} s'.format(end-start))
    return edges


### count two-neuron motifs

def get_unidirectional_only(E):
    unidirectional = set()
    for e in E:
        if (e[::-1] not in E) and (e[0] != e[1]):
            unidirectional.add(e)
    return unidirectional

def get_bidirectional(E):
    bidirectional = set()
    for e in E:
        if (e[::-1] in E) and (e[::-1] not in bidirectional) and (e[0] != e[1]):
            bidirectional.add(e)
    return bidirectional

def get_autapses(E):
    autapses = set()
    for e in E:
        if e[0] == e[1]:
            autapses.add(e)
    return autapses

def split_edge_set_to_two_patterns(E):
    return get_unidirectional_only(E).copy(), get_bidirectional(E).copy()

def merge_two_patterns_to_edge_set(uni_edges, bi_edges):
    return uni_edges.union(bi_edges).union({e[::-1] for e in bi_edges})

def graph_to_edges_sets(g):
    V, E = graph_to_sets(g)
    return split_edge_set_to_two_patterns(E)

def edge_sets_to_graph(uni_edges, bi_edges):
    return nx.DiGraph(merge_two_patterns_to_edge_set(uni_edges, bi_edges))

def count_two_neuron_motifs_graph(g):
    return count_two_neuron_motifs(*graph_to_sets(g))

def count_two_neuron_motifs(V, E):
    d = {}
    d['neurons'] = len(V)
    d['autapses'] = len(get_autapses(E))
    d['actual_edges'] = len(E) - d['autapses']
    d['uni'] = len(get_unidirectional_only(E))
    d['bi'] = len(get_bidirectional(E))
    d['potential_edges'] = (d['neurons']*(d['neurons']-1))
    d['null'] = (d['potential_edges'] // 2) - d['uni'] - d['bi']
    return d    

def compute_ER_two_neuron_motifs(g):
    tn = count_two_neuron_motifs_graph(g)
    p = tn['actual_edges'] / tn['potential_edges']
    n = tn['potential_edges'] // 2
    obs = [tn['null'], tn['uni'], tn['bi']]    
    return {'null': n*(1-p)**2, 
            'uni': 2*n*p*(1-p),
            'bi': n*p**2}

### count three-neuron motifs

class Triplet():
    """Class representation of a triplet for easy comparison 
    """
    def __init__(self, vertices, edges=set()):
        """Args:
            vertices: list of vertices
            edges: set of pairs of vertices as tuples
        """
        vertices.sort()
        self.vertices = {k: v for k, v in enumerate(vertices)}
        self.edges = set()
        for e in edges: self.add_edge(*e)
    
    @property
    def ids(self):
        return {v: k for k, v in self.vertices.items()}
    
    def get_id(self, v):
        return self.ids[v]
    
    def get_edges(self):
        return [(self.vertices[e[0]], self.vertices[e[1]]) for e in self.edges]
    
    def add_edge(self, u, v):
        e = (self.get_id(u), self.get_id(v))
        self.edges.add(e)
    
    def __len__(self):
        return len(self.edges)
    
    def __repr__(self):
        return str(self.vertices) + ' ' + str(sorted(self.edges))

    def __hash__(self):
        return hash(self.__repr__())
    
    def __eq__(self, X):
        """Is this triplet equivalent (vertices & edges) without any permutation?
        """
        return self.vertices == X.vertices and self.edges == X.edges

def reflect(edges, identity):
    """Reflect triplet edges across axis of reflection_vertex, return new set of edges
    
    Args:
        edges: Triplet.edges
        identity: int in [0,1,2] for vertex on axis of reflection
    """
    reflect = [x for x in range(3) if x != identity]
    remap = {identity: identity, reflect[0]: reflect[1], reflect[1]: reflect[0]} # old_id: new_id
    return set([(remap[u], remap[v]) for u, v in edges])

def match_edges(A, B, reflections=range(3)):
    """Are these two edge lists the same motif? (equivalent under any reflection symmetry)
    Explores tree of reflections with recursion, e.g. A.reflect(0).reflect(1).reflect(2), etc
    
    Args:
        A: Triplet.edges
        B: Triplet.edges
    """
    if A == B:
        return True
    if len(reflections) > 1:
        for k in reflections:
            if match_edges(reflect(A, k), B, reflections=[r for r in reflections if r != k]):
                return True
    return False

def match(A, B):
    """Are these two triplets the same motif?
    
    Args:
        A: Triplet
        B: Triplet
    """
    return match_edges(A.edges, B.edges)

def collect_triplets_graph(G):
    return collect_triplets(G.nodes(), G.edges())

def collect_triplets(V, E):
    tri = set([])
    for a, b in E:
        for c in V:
            if c != a and c != b:
                t = Triplet([a, b, c])
                t.add_edge(a,b)
                if (a, c) in E:
                    t.add_edge(a,c)
                if (c, a) in E:
                    t.add_edge(c,a)
                if (b, c) in E:
                    t.add_edge(b,c)
                if (c, b) in E:
                    t.add_edge(c,b)
                if (b, a) in E:
                    t.add_edge(b,a)
                tri.add(t)
    return tri

t = [0,1,2]
motifs = {
    1: Triplet(t, edges=set([])),
    2: Triplet(t, edges=set([(2,1)])),
    3: Triplet(t, edges=set([(2,1),(1,2)])),
    4: Triplet(t, edges=set([(1,0),(1,2)])),
    5: Triplet(t, edges=set([(0,1),(2,1)])),
    6: Triplet(t, edges=set([(0,1),(1,2)])),
    7: Triplet(t, edges=set([(0,1),(1,2),(2,1)])),
    8: Triplet(t, edges=set([(1,0),(1,2),(2,1)])),
    9: Triplet(t, edges=set([(0,1),(1,0),(1,2),(2,1)])),
   10: Triplet(t, edges=set([(0,2),(1,0),(1,2)])),
   11: Triplet(t, edges=set([(0,2),(2,1),(1,0)])),
   12: Triplet(t, edges=set([(0,1),(0,2),(1,2),(2,1)])),
   13: Triplet(t, edges=set([(0,2),(1,0),(1,2),(2,1)])),
   14: Triplet(t, edges=set([(1,0),(2,0),(1,2),(2,1)])),
   15: Triplet(t, edges=set([(0,2),(0,1),(1,0),(1,2),(2,1)])),
   16: Triplet(t, edges=set([(0,2),(2,0),(0,1),(1,0),(1,2),(2,1)])),
}
# count all permutations of triplets
all_triplets = []
for i in [[None], [(0,1)], [(1,0)], [(0,1),(1,0)]]:
    for j in [[None], [(1,2)], [(2,1)], [(1,2),(2,1)]]:
        for k in [[None], [(0,2)], [(2,0)], [(0,2),(2,0)]]:
            edges = set([])
            for e in i + j + k:
                if e is not None:
                    edges.add(e)
            all_triplets.append(Triplet(t, edges=edges))            
            
# count how many triplet permutations correspond to each motif
motif_factor = {}
for k, m in motifs.items():
    motif_factor[k] = sum([match(m, t) for t in all_triplets])

def collect_three_neuron_motifs(V, E, motifs):
    tri_g = collect_triplets(V, E)
    matches = {k: [] for k in motifs.keys()}
    n = len(tri_g)
    for i, t in enumerate(tri_g):
        # print('{} / {}'.format(n, i), end='\r', flush=True)
        t_match = False
        for k, m in motifs.items():
            if match(t, m):
                t_match = True
                matches[k].append(t)
    return matches		    

def count_three_neuron_motifs(V, E, motifs):
    n = len(V)
    N = n*(n-1)*(n-2) // 6
    tri_motifs = collect_three_neuron_motifs(V, E, motifs)
    tri_counts = {k: len(v) for k, v in tri_motifs.items()}
    # compute unconnected triplets (motif #1)
    tri_counts[1] = N - sum([len(v) for v in tri_motifs.values()])
    return tri_counts, tri_motifs


def compute_three_neuron_motif_probabilities(g):
    # generate base expectations
    tn = count_two_neuron_motifs_graph(g)
    n = tn['potential_edges']
    p = (tn['uni'] + tn['bi']*2) / n
    pr_uni = 2 * p * (1-p) / 2
    pr_bi = p**2
    pr_null = (1-p)**2
    triplet_pr = {
        1: pr_null**3,
        2: pr_null**2 * pr_uni,
        3: pr_null**2 * pr_bi,
        4: pr_null * pr_uni**2,
        5: pr_null * pr_uni**2,
        6: pr_null * pr_uni**2,
        7: pr_null * pr_uni * pr_bi,
        8: pr_null * pr_uni * pr_bi,
        9: pr_null * pr_bi**2,
       10: pr_uni**3,
       11: pr_uni**3,
       12: pr_uni**2 * pr_bi,
       13: pr_uni**2 * pr_bi,
       14: pr_uni**2 * pr_bi,
       15: pr_uni * pr_bi**2,
       16: pr_bi**3,
    }

    # adjust motif probability with motif factor
    triplet_pr = {k: v*motif_factor[k] for k, v in triplet_pr.items()}
    pair_pr = {'pr_uni': pr_uni, 'pr_bi': pr_bi, 'pr_null': pr_null}
    return triplet_pr, pair_pr


def compute_three_neuron_motif_probabilities_GE(g):
    # generate base expectations w/ preseved biedge stats
    tn = count_two_neuron_motifs_graph(g)
    n = tn['potential_edges'] // 2
    pr_uni = tn['uni'] / n / 2
    pr_bi = tn['bi'] / n
    pr_null = 1 - pr_uni*2 - pr_bi
    triplet_pr = {
        1: pr_null**3,
        2: pr_null**2 * pr_uni,
        3: pr_null**2 * pr_bi,
        4: pr_null * pr_uni**2,
        5: pr_null * pr_uni**2,
        6: pr_null * pr_uni**2,
        7: pr_null * pr_uni * pr_bi,
        8: pr_null * pr_uni * pr_bi,
        9: pr_null * pr_bi**2,
       10: pr_uni**3,
       11: pr_uni**3,
       12: pr_uni**2 * pr_bi,
       13: pr_uni**2 * pr_bi,
       14: pr_uni**2 * pr_bi,
       15: pr_uni * pr_bi**2,
       16: pr_bi**3,
    }

    # adjust motif probability with motif factor
    triplet_pr = {k: v*motif_factor[k] for k, v in triplet_pr.items()}
    pair_pr = {'pr_uni': pr_uni, 'pr_bi': pr_bi, 'pr_null': pr_null}
    return triplet_pr, pair_pr


def compute_expected_three_neuron_motifs(g, prob_dict):
    n = len(g.nodes())
    N = n*(n-1)*(n-2) // 6
    return {k: p*N for k, p in prob_dict.items()}


### sample and count
def sample_config_two_neuron_motifs(G, samples=100, niters=int(1e5)):
    df = pd.DataFrame()
    V, E = graph_to_sets(G)
    for i in range(samples):
        print('{} / {} samples'.format(i+1, samples), end='\r', flush=True)
        rE = switch_and_hold(E, niters=niters)
        tn = count_two_neuron_motifs(V, rE)
        df = df.append(tn, ignore_index=True)
    return df


def sample_config_three_neuron_motifs(G, samples, niters=int(1e5)):
    df = pd.DataFrame()
    V, E = graph_to_sets(G)
    for i in range(samples):
#         print('{} / {} samples'.format(i+1, samples), end='\r', flush=True)
        rE = switch_and_hold(E, niters=niters)
        tn, _ = count_three_neuron_motifs(V, rE, motifs)
        df = three_df.append(tn, ignore_index=True)
    return df


def single_sample_two_neuron_motifs(V, E, samples, niters, i):
    print('{} / {} samples'.format(i+1, samples), end='\r', flush=True)
    rE = switch_and_hold(E, niters=niters)
    df = count_two_neuron_motifs(V, rE)
    return df


def continuous_sample_two_neuron_motifs(V, E, samples, niters, thread, i):
    # print('{} / {} threads'.format(i+1, thread))
    df = pd.DataFrame()
    for s in range(samples):
        # print('{}: {} / {} samples'.format(i+1, s+1, samples))
        E = switch_and_hold(E, niters=niters)
        tdf = count_two_neuron_motifs(V, E)
        df = pd.concat((df, pd.DataFrame([tdf])))
    return df

def continuous_sample_two_neuron_motifs_permitted(V, E, pE, samples, niters, thread, i):
    # print('{} / {} threads'.format(i+1, thread))
    df = pd.DataFrame()
    for s in range(samples):
        # print('{}: {} / {} samples'.format(i+1, s+1, samples))
        E = switch_and_hold_permitted(E, pE, niters=niters)
        tdf = count_two_neuron_motifs(V, E)
        df = pd.concat((df, pd.DataFrame([tdf])))
    return df

def continuous_sample_three_neuron_motifs(V, E, samples, niters, thread, i):
    # print('{} / {} threads'.format(i+1, thread))
    df = pd.DataFrame()
    for s in range(samples):
#         print('{}: {} / {} samples'.format(i+1, s+1, samples))
        E = switch_and_hold(E, niters=niters)
        tdf, _ = count_three_neuron_motifs(V, E, motifs)
        df = pd.concat((df, pd.DataFrame([tdf])), ignore_index=True)
    return df


def continuous_sample_three_neuron_motifs_GE(V, E, samples, niters, thread, i):
#     print('{} / {} threads'.format(i+1, thread))
    df = pd.DataFrame()
    for s in range(samples):
#         print('{}: {} / {} samples'.format(i+1, s+1, samples))
        E = switch_and_hold_GE(E, niters=niters)
        tdf, _ = count_three_neuron_motifs(V, E, motifs)
        df = pd.concat((df, pd.DataFrame([tdf])), ignore_index=True)
    return df

def sample_motifs_parallel_continuous(fn, G, samples, niters=int(1e5), threads=8):
    V, E = graph_to_sets(G)
    f = partial(fn, V, E, samples, niters, threads)
    with ProcessPoolExecutor(max_workers=threads) as executor:
        dfs = executor.map(f, range(threads))
    return pd.concat(dfs, ignore_index=True)

### Clustering Coefficient
def clustering_coef(counts):
    return 3 * counts[np.arange(10,17)].sum(1) / (counts[np.arange(4,10)].sum(1) + 3 * counts[np.arange(10,17)].sum(1))
