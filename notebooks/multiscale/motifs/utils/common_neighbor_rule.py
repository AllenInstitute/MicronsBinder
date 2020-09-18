import itertools
from scipy import stats
from .motif_counts import *

def get_cnr_stats(g_):
    g = g_.to_undirected()
    cn_counts = np.zeros(10)
    scn_counts = np.zeros(10)
    pcn_counts = np.zeros(10)
    ccn_counts = np.zeros(10)
    
    cn_connected_counts = np.zeros(10)
    scn_connected_counts = np.zeros(10)
    scn_uni_connected_counts = np.zeros(10)
    scn_bi_connected_counts = np.zeros(10)
    pcn_connected_counts = np.zeros(10)
    pcn_uni_connected_counts = np.zeros(10)
    pcn_bi_connected_counts = np.zeros(10)
    ccn_connected_counts = np.zeros(10)
    ccn_uni_connected_counts = np.zeros(10)
    ccn_bi_connected_counts = np.zeros(10)
    
    precent_conn= np.zeros(10)
    sprecent_conn= np.zeros(10)
    sprecent_uni_conn= np.zeros(10)
    sprecent_bi_conn= np.zeros(10)
    pprecent_conn= np.zeros(10)
    pprecent_uni_conn= np.zeros(10)
    pprecent_bi_conn = np.zeros(10)
    cprecent_conn= np.zeros(10)
    cprecent_uni_conn= np.zeros(10)
    cprecent_bi_conn= np.zeros(10)
    
    node_pairs = list(itertools.combinations(list(g_.nodes), 2))
    for ns in node_pairs:
        
        # common neighborhood
        cns = len(list(nx.common_neighbors(g, ns[0], ns[1])))
        if cns > 9: 
            cns = 9
        cn_counts[cns] += 1
        if g.has_edge(ns[0], ns[1]):
            cn_connected_counts[cns] += 1
            
        # common successor - unidirection
        cns_1 = list(g_.successors(ns[0]))
        cns_2 = list(g_.successors(ns[1]))
        if len(cns_1) > 0 and len(cns_2) > 0:
            intersect = np.intersect1d(cns_1, cns_2)
            cns = len(intersect)
            # remove bidirectional edges..
            for n in intersect:
                n_succ = g_.successors(n)
                if ns[0] in n_succ or ns[1] in n_succ:
                    cns -= 1
        else:
            cns = 0
        
        if cns > 9: 
            cns = 9
    
        scn_counts[cns] += 1
        
        if g_.has_edge(ns[0], ns[1]) or g_.has_edge(ns[1], ns[0]) :
            scn_connected_counts[cns] += 1
        if g_.has_edge(ns[0], ns[1]) and g_.has_edge(ns[1], ns[0]) :
            scn_bi_connected_counts[cns] += 1
            
        # common predecessors - unidirection
        cns_1 = list(g_.predecessors(ns[0]))
        cns_2 = list(g_.predecessors(ns[1]))
        cns = len(np.intersect1d(cns_1, cns_2))
        if len(cns_1) > 0 or len(cns_1) > 0:
            intersect = np.intersect1d(cns_1, cns_2)
            cns = len(intersect)
            # remove bidirectional edges..
            for n in intersect:
                n_pred = g_.predecessors(n)
                if ns[0] in n_pred or ns[1] in n_pred:
                    cns -= 1
        else:
            cns = 0
        
        if cns > 9: 
            cns = 9
    
        pcn_counts[cns] += 1
        
        if g_.has_edge(ns[0], ns[1]) or g_.has_edge(ns[1], ns[0]) :
            pcn_connected_counts[cns] += 1
        if g_.has_edge(ns[0], ns[1]) and g_.has_edge(ns[1], ns[0]) :
            pcn_bi_connected_counts[cns] += 1
    
    scn_uni_connected_counts = scn_connected_counts - scn_bi_connected_counts
    pcn_uni_connected_counts = pcn_connected_counts - pcn_bi_connected_counts
            
    for cnt in range(10):
        if cn_counts[cnt] > 0:
            precent_conn[cnt] = cn_connected_counts[cnt]/cn_counts[cnt]
        if scn_counts[cnt] > 0:
            sprecent_conn[cnt] = scn_connected_counts[cnt]/scn_counts[cnt]
            sprecent_uni_conn[cnt] = scn_uni_connected_counts[cnt]/scn_counts[cnt]
            sprecent_bi_conn[cnt] = scn_bi_connected_counts[cnt]/scn_counts[cnt]
        if pcn_counts[cnt] > 0:
            pprecent_conn[cnt] = pcn_connected_counts[cnt]/pcn_counts[cnt]
            pprecent_uni_conn[cnt] = pcn_uni_connected_counts[cnt]/pcn_counts[cnt]
            pprecent_bi_conn[cnt] = pcn_bi_connected_counts[cnt]/pcn_counts[cnt]
        
    return pd.DataFrame({
    		# number of pairs with [bin] common neighbor
            "undir_pair": cn_counts,
            # number of pairs with [bin] common successor (strict, not a predecessor to any)
            "dir_spair": scn_counts,
            # number of pairs with [bin] common predecessor (strict, not a successor to any)
            "dir_ppair": pcn_counts,
            # number of connected pairs with [bin] common neighbor
            "undir_conn": cn_connected_counts,
            # number of connected pairs with [bin] common successor (strict, not a predecessor to any)
            "dir_sconn": scn_connected_counts,
            # number of connected pairs with [bin] common predecessor (strict, not a successor to any)
            "dir_pconn": pcn_connected_counts,
            # number of unidirectionally connected pairs with [bin] common successor (strict, not a predecessor to any)
            "dir_uni_sconn": scn_uni_connected_counts,
            # number of unidirectionally connected pairs with [bin] common predecessor (strict, not a successor to any)
            "dir_uni_pconn": pcn_uni_connected_counts,
            # number of bidirectionally connected pairs with [bin] common successor (strict, not a predecessor to any)
            "dir_bi_sconn": scn_bi_connected_counts,
            # number of bidirectionally connected pairs with [bin] common predecessor (strict, not a successor to any)
            "dir_bi_pconn": pcn_bi_connected_counts,
            # percentage of pairs with [bin] common neighbor to be connected
            "undir_perc": precent_conn,
            # percentage of pairs with [bin] strict common successor to be connected
            "dir_sperc": sprecent_conn,
            # percentage of pairs with [bin] strict common predcessor to be connected
            "dir_pperc": pprecent_conn}) 

def continuous_sample_cnr(g, samples, niters):
    V, E = graph_to_sets(g)
    ft_config_ucnrs = []
    hrates = []
#     print('{} / {} threads'.format(i+1, thread))
    df = pd.DataFrame()
    for s in trange(samples):
#         print('{}: {} / {} samples'.format(i+1, s+1, samples))
        E = switch_and_hold(E, niters=niters)
        g_sampled = set_to_graph(E)
        ft_config_ucnr = get_cnr_stats(g_sampled)
        ft_config_ucnrs.append(ft_config_ucnr.to_numpy())
    return ft_config_ucnrs


def continuous_sample_cnr_GE(g, samples, niters):
    V, E = graph_to_sets(g)
    ft_config_ucnrs = []
    hrates = []
#     print('{} / {} threads'.format(i+1, thread))
    df = pd.DataFrame()
    for s in trange(samples):
#         print('{}: {} / {} samples'.format(i+1, s+1, samples))
        E = switch_and_hold_GE(E, niters=niters)
        g_sampled = set_to_graph(E)
        ft_config_ucnr = get_cnr_stats(g_sampled)
        ft_config_ucnrs.append(ft_config_ucnr.to_numpy())
    return ft_config_ucnrs

def cnr_mean_std(ft_config_ucnrs):
	ft_config_ucnr = pd.DataFrame(columns=["undir_pair", "dir_spair", "dir_ppair", "undir_conn", "dir_sconn", "dir_pconn", "dir_uni_sconn", 
										   "dir_uni_pconn", "dir_bi_sconn", "dir_bi_pconn", "undir_perc", "dir_sperc", "dir_pperc"],
                          data=np.mean(ft_config_ucnrs,0))
	ft_config_ucnr_std = pd.DataFrame(columns=["undir_pair_std", "dir_spair_std", "dir_ppair_std", "undir_conn_std", 
                                       "dir_sconn_std", "dir_pconn_std", "dir_uni_sconn_std", "dir_uni_pconn_std", 
                                       "dir_bi_sconn_std", "dir_bi_pconn_std", "undir_perc_std", "dir_sperc_std", 
                                       "dir_pperc_std"],
                          data=np.std(ft_config_ucnrs,0))
	return ft_config_ucnr, ft_config_ucnr_std

def pearson_r2(x, y):
    r, p = stats.pearsonr(np.arange(10), y)
    return r ** 2, p	