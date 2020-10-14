import random
import scipy
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import rc, rcParams
rc('text', usetex=False)
rc('axes', linewidth=1.5)
rc('font', weight='regular')

cmap = plt.cm.tab20
cmaplist = [cmap(17*i % 20) for i in range(cmap.N)]

def plot_two_neuron_counts_violin(obs_counts, er_counts, config_counts, subtitle='', ylim=None, figsize=(6,4), fname="temp"):
    """Compare observed, ER, and configuration two-neuron counts
    
    Args:
        obs_counts: pd.DataFrame from count_two_neuron_motifs
        er_counts: pd.DataFrame from compute_ER_two_neuron_motifs
        config_counts: pd.DataFrame from sample_config_two_neuron_motifs
    """
    p = obs_counts['actual_edges'] / obs_counts['potential_edges']
    print("connect probability = ", p)
    f, ax = plt.subplots(1, 1, figsize=figsize)
    x_labels = ['o   o', 'o$\Rightarrow$o', 'o$\Leftrightarrow$o']
    keys = ['null', 'uni', 'bi']
    
    data = [obs_counts[k] / er_counts[k] if er_counts[k] != 0 else 0 for k in keys]
    ax.scatter(np.arange(len(x_labels))+1, data, color='#D8545D', 
               s=150, marker='x', zorder=100, label="OBS/ER")
    
    mean = config_counts.mean()
    std = config_counts.std()    
    rdm = np.array([mean[k] / er_counts[k] if er_counts[k] != 0 else 0 for k in keys])
    err = [std[k] / er_counts[k] for k in keys]
#     ax.errorbar(np.arange(len(x_labels)), rdm, yerr=err, color='r', fmt='o', alpha=1)
    data = [config_counts[k] / er_counts[k] if er_counts[k] != 0 else 0 for k in keys]
    parts = ax.violinplot(data, widths=0.7,
                 showmeans=False, showmedians=False, showextrema=False)
    for pc in parts['bodies']:
#         pc.set_facecolor([0.0,0.5,0.9,0.5])
        pc.set_facecolor('#C0C0C0')
#         pc.set_edgecolor('black')
        pc.set_alpha(0.5)

    quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=1)
    p5, p95 = np.percentile(data, [5, 95], axis=1)
    inds = np.arange(1, len(medians) + 1)
    ax.scatter(inds, medians, marker='o', color='white', edgecolors='k', 
               s=80, zorder=3, label="CFG/ER")
    ax.vlines(inds, quartile1, quartile3, color='k', alpha=0.6, linestyle='-', lw=2.5)
    ax.vlines(inds, p5, p95, color='k', alpha=0.6, linestyle='-.', lw=1.0)
    mid, mlen = (p5+p95)/2, (p95-p5)/2
    eb = plt.errorbar(inds, (p5+p95)/2, yerr = mlen, ls='none', color='k', capsize=2, elinewidth=1.0) 
    eb[-1][0].set_linestyle('-.')
    
    ax.legend(bbox_to_anchor=(0.18, 0.73), loc=8, fontsize=14)
    
    ax.plot(range(0,len(keys)+2), [1]*(len(keys)+2), '--', color='#745e4f', alpha=1.0)
    
    
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')

    # was: axes.spines['bottom'].set_position(('data',1.1*X.min()))
    ax.spines['bottom'].set_position(('axes', -0.05))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('axes', -0.05))
    
    ax.tick_params(length=8, width=1.5, labelsize=15);
    ax.set_xlim(left=0.5,right=3.5)
    ax.set_xticks(np.arange(1, len(x_labels)+1))
    ax.set_ylim(bottom=0)
    ax.set_xticklabels(x_labels, fontsize=20);
    ax.set_ylabel('Counts / ER Mean Counts', fontsize=15);
    
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_title('{}'.format(subtitle), fontsize=20);
    f.savefig("figures/{}.pdf".format(fname), bbox_inches='tight')
    d = pd.DataFrame.from_dict({'obs': {k: obs_counts[k] for k in keys}, 
                                'er_exp': er_counts, 'config_exp': mean, 
                                'config_std': std}, orient='columns')
    d['config_z'] = (d['obs'] - d['config_exp']) / d['config_std']
    # print(d.loc[['null','uni','bi'],:].round({'obs': 0, 'null_exp': 0, 'config_exp': 0, 'config_std': 2, 'config_z': 1}))
    return d

def plot_two_neuron_counts_scatter(obs_two_counts, er_two_counts, config_two_counts, figsize=(5,5), fname="2-scatter"):
    f, ax = plt.subplots(1, 1, figsize=figsize)

    ax.scatter([obs_two_counts['null']], [er_two_counts['null']], 
                marker='^', color='#1f77b4', s=120, zorder=2, alpha=0.9, label='ER   o   o')
    ax.scatter(obs_two_counts['null'], config_two_counts['null'].mean(), 
                marker='o', color='#1f77b4', s=120, zorder=4, alpha=0.9, label='CFG o   o')
    ax.plot([obs_two_counts['null'], obs_two_counts['null']], 
            [er_two_counts['null'], config_two_counts['null'].mean()], 
             '-', color='#1f77b4', linewidth=1, zorder=4, alpha=1.0)

    ax.scatter([obs_two_counts['uni']], [er_two_counts['uni']], 
                marker='^', color='#ff7f0e', s=120, zorder=3, alpha=0.9, label='ER   o$\Rightarrow$o')
    ax.scatter(obs_two_counts['uni'], config_two_counts['uni'].mean(),  
                marker='o', color='#ff7f0e', s=120, zorder=4, alpha=0.9, label='CFG o$\Rightarrow$o')
    ax.plot([obs_two_counts['uni'], obs_two_counts['uni']], 
            [er_two_counts['uni'], config_two_counts['uni'].mean()], 
             '-', color='#ff7f0e', linewidth=1, zorder=4, alpha=1.0)

    ax.scatter([obs_two_counts['bi']], [er_two_counts['bi']], 
                marker='^', color='#2ca02c', s=120, zorder=4, alpha=0.9, label='ER   o$\Leftrightarrow$o')
    ax.scatter(obs_two_counts['bi'], config_two_counts['bi'].mean(), 
                marker='o', color='#2ca02c', s=120, zorder=4, alpha=0.9, label='CFG o$\Leftrightarrow$o')
    ax.plot([obs_two_counts['bi'], obs_two_counts['bi']], 
            [er_two_counts['bi'], config_two_counts['bi'].mean()], 
             '-', color= '#2ca02c', linewidth=1, zorder=4, alpha=1.0)

    x = np.linspace(0, obs_two_counts['null'], 100)
    ax.plot(x,x, '--k', zorder=1)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines['left'].set_linewidth(1.25)
    ax.yaxis.set_tick_params(width=2)
    ax.spines['bottom'].set_linewidth(1.25)
    ax.xaxis.set_tick_params(width=2)

    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')

    # was: axes.spines['bottom'].set_position(('data',1.1*X.min()))
    ax.spines['bottom'].set_position(('axes', -0.05))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('axes', -0.05))

    ax.set_xlim([10,1.2*obs_two_counts['null']])
    ax.set_ylim([10,1.2*obs_two_counts['null']])

    ax.tick_params(length=8, width=1.5, labelsize=15);
    ax.legend(fontsize=15, bbox_to_anchor=(0.20, 0.45), loc=8)

    plt.xlabel('Actual counts', fontsize=20)
    plt.ylabel('Sampled counts', fontsize=20)

    plt.xscale('log')
    plt.yscale('log')
    f.savefig("figures/{}.pdf".format(fname), bbox_inches='tight')

def plot_three_neuron_counts_violin(obs_counts, er_counts, config_counts, subtitle='', ylim=None, fname="temp", figsize=(10,3)):
    f, ax = plt.subplots(1, 1, figsize=figsize)
    x_labels = list(obs_counts.keys()) #['o   o', 'o$\Rightarrow$o', 'o$\Leftrightarrow$o']
    obs_data = [obs_counts[k] / er_counts[k] for k in x_labels]
    # ax.bar(np.arange(len(x_labels))+1, obs_data, bottom=1, color='k', error_kw={'elinewidth': 1, 'capsize': 3})
    ax.scatter(np.arange(len(x_labels))+1, obs_data, color='#D8545D', s=100, marker='x', zorder=150, label="OBS/ER")
    
    mean = config_counts.mean()
    std = config_counts.std()    
    rdm = np.array([mean[k] / er_counts[k] for k in x_labels])
    err = [std[k] / er_counts[k] for k in x_labels]
#     ax.errorbar(np.arange(len(x_labels))+1, rdm, yerr=err, color='r', fmt='o', alpha=1)
    
    
    data = [config_counts[k] / er_counts[k] if er_counts[k] != 0 else 0 for k in x_labels]
    parts = ax.violinplot(data, widths=0.7,
                 showmeans=False, showmedians=False, showextrema=False)
    
    for pc in parts['bodies']:
        pc.set_facecolor('#C0C0C0')
#         pc.set_edgecolor('black')
        pc.set_alpha(0.5)
        
    quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=1)
    p5, p95 = np.percentile(data, [5, 95], axis=1)
    inds = np.arange(1, len(medians) + 1)
    ax.scatter(inds, medians, marker='o', color='w', edgecolors='k', s=80, zorder=3, label="CFG/ER")
    ax.vlines(inds, quartile1, quartile3, color='k', alpha=0.6, linestyle='-', lw=2.5)
#     ax.vlines(inds, p5, p95, color='k', alpha=0.6, linestyle='-.', lw=0.75)
    mid, mlen = (p5+p95)/2, (p95-p5)/2
    eb = plt.errorbar(inds, (p5+p95)/2, yerr = mlen, ls='none', color='k', capsize=2, elinewidth=1.0)
    eb[-1][0].set_linestyle('-.')
    
    ax.plot(range(0,len(x_labels)+2), [1]*(len(x_labels)+2), '--', color='#745e4f', alpha=1.0)
    
    ax.legend(bbox_to_anchor=(0.12, 0.71), loc=8, fontsize=15)
    
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')

    # was: axes.spines['bottom'].set_position(('data',1.1*X.min()))
    ax.spines['bottom'].set_position(('axes', -0.05))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('axes', -0.05))
    
    ax.tick_params(length=8, width=1.5, labelsize=15);
    ax.set_xlim(left=0.5,right=16.5)
    ax.set_xticks(np.arange(len(x_labels))+1)
    # ax.set_xticklabels(x_labels, fontsize=20);
    
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_ylabel('Counts / ER Mean Counts', fontsize=15);
#     ax.set_title('Three-neuron motifs{}' \
#                  '\nbootstrap iterations={}' \
#                  '\np(o   o)={:.3f}' \
#                  '\np(o$\Rightarrow$o)={:.3f}' \
#                  '\np(o$\Leftrightarrow$o)v={:.4f}'.format(subtitle,
#                                                           len(bootstrap),
#                                                           pair_pr['pr_null'], 
#                                                           pair_pr['pr_uni'], 
#                                                           pair_pr['pr_bi']), fontsize=18);
    ax.set_title(''.format(subtitle), fontsize=18);
    f.savefig("figures/{}.pdf".format(fname), bbox_inches='tight')
    d = pd.DataFrame.from_dict({'obs':obs_counts, 'er_exp': er_counts, 'config_exp': mean, 'config_std': std}, orient='columns')
    d['config_z'] = (d['obs'] - d['config_exp']) / d['config_std']
    # print(d.round({'obs': 0, 'null_exp': 0, 'config_exp': 0, 'config_std': 2, 'config_z': 1}))
    return d

def plot_three_neuron_counts_diff_violin(obs_counts, er_counts, config_counts, subtitle='', ylim=None, fname="temp", figsize=(10,3)):
    f, ax = plt.subplots(1, 1, figsize=figsize)
    x_labels = list(obs_counts.keys()) #['o   o', 'o$\Rightarrow$o', 'o$\Leftrightarrow$o']
    obs_data = [obs_counts[k] for k in x_labels]
    # ax.bar(np.arange(len(x_labels))+1, obs_data, bottom=1, color='k', error_kw={'elinewidth': 1, 'capsize': 3})
#     ax.scatter(np.arange(len(x_labels))+1, obs_data, color='g', s=100, marker='x', zorder=100, label="OBS/ER")
    
    ax.scatter(np.arange(len(x_labels))+1, obs_data-config_counts.mean(0), color='#D8545D', s=150, marker='x', 
               zorder=100, label="OBS")
    
    mean = config_counts.mean()
    std = config_counts.std()    
    rdm = np.array([mean[k] / er_counts[k] for k in x_labels])
    err = [std[k] / er_counts[k] for k in x_labels]
#     ax.errorbar(np.arange(len(x_labels))+1, rdm, yerr=err, color='r', fmt='o', alpha=1)
    
    
    data = [(config_counts[k] - config_counts[k].mean()) for k in x_labels]
    parts = ax.violinplot(data, widths=0.7,
                 showmeans=False, showmedians=False, showextrema=False)
    
    for pc in parts['bodies']:
        pc.set_facecolor('#C0C0C0')
#         pc.set_edgecolor('black')
        pc.set_alpha(0.5)
        
    quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=1)
    quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=1)
    p5, p95 = np.percentile(data, [5, 95], axis=1)
    inds = np.arange(1, len(medians) + 1)
#     ax.scatter(inds, medians, marker='o', color='w', edgecolors='k', s=50, zorder=3, label="CFG/ER")
    ax.vlines(inds, quartile1, quartile3, color='k', alpha=0.6, linestyle='-', lw=2.5)
#     ax.vlines(inds, p5, p95, color='k', alpha=0.6, linestyle='-.', lw=0.75)
    mid, mlen = (p5+p95)/2, (p95-p5)/2
    eb = plt.errorbar(inds, (p5+p95)/2, yerr = mlen, ls='none', color='k', capsize=2, elinewidth=1.0) 
    eb[-1][0].set_linestyle('-.')
    
    ax.plot(range(0,len(x_labels)+2), [1]*(len(x_labels)+2), '--', color='#745e4f', alpha=1.0)
    
#     ax.legend(bbox_to_anchor=(0.8, 0.71), loc=8, fontsize=15)
    
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')

    # was: axes.spines['bottom'].set_position(('data',1.1*X.min()))
    ax.spines['bottom'].set_position(('axes', -0.05))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('axes', -0.05))
    
    ax.tick_params(length=8, width=1.5, labelsize=14);
    ax.set_xticks(np.arange(len(x_labels))+1)
    ax.set_xlim(left=0.5,right=16.5)
    # ax.set_xticklabels(x_labels, fontsize=20);
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_ylabel('Counts - CFG Mean Counts', fontsize=15);
#     ax.set_title('Three-neuron motifs{}' \
#                  '\nbootstrap iterations={}' \
#                  '\np(o   o)={:.3f}' \
#                  '\np(o$\Rightarrow$o)={:.3f}' \
#                  '\np(o$\Leftrightarrow$o)v={:.4f}'.format(subtitle,
#                                                           len(bootstrap),
#                                                           pair_pr['pr_null'], 
#                                                           pair_pr['pr_uni'], 
#                                                           pair_pr['pr_bi']), fontsize=18);
    ax.set_title('{}'.format(subtitle), fontsize=18);
    f.savefig("figures/{}.pdf".format(fname), bbox_inches='tight')
    d = pd.DataFrame.from_dict({'obs':obs_counts, 'er_exp': er_counts, 'config_exp': mean, 'config_std': std}, orient='columns')
    d['config_z'] = (d['obs'] - d['config_exp']) / d['config_std']
    # print(d.round({'obs': 0, 'null_exp': 0, 'config_exp': 0, 'config_std': 2, 'config_z': 1}))
    return d

def plot_three_neuron_counts_scatter(obs_three_counts, er_three_counts, config_three_counts, figsize=(12,12), fname="temp"):
    f, ax = plt.subplots(1, 1, figsize=figsize)

    for i in range(1,17):
        ax.scatter([obs_three_counts[i]], [er_three_counts[i]], 
                    marker='^', color=cmaplist[i], s=250, zorder=2*i+1, alpha=0.75, label='ER   {} '.format(i))
        ax.scatter(obs_three_counts[i], config_three_counts[i].mean(), 
                    marker='o', color=cmaplist[i], s=250, zorder=2*i+2, alpha=0.75, label='CFG {} '.format(i))
    #     ax.errorbar(obs_three_counts[i], config_three_counts[i].mean(), yerr=config_three_counts[i].std(), 
    #                 color=cmaplist[i], alpha=1)
        ax.plot([obs_three_counts[i], obs_three_counts[i]], 
            [er_three_counts[i], config_three_counts[i].mean()], 
             color=cmaplist[i], linewidth=1, alpha=0.75)

    x = np.linspace(-1, obs_three_counts[1], 100)
    ax.plot(x,x, '--k', lw=2, zorder=1)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines['left'].set_linewidth(2)
    ax.yaxis.set_tick_params(width=2)
    ax.spines['bottom'].set_linewidth(2)

    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')

    # was: axes.spines['bottom'].set_position(('data',1.1*X.min()))
    ax.spines['bottom'].set_position(('axes', -0.025))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('axes', -0.025))

    ax.set_xlim([-0.5,1.2*obs_three_counts[1]])
    ax.set_ylim([-0.5,1.2*obs_three_counts[1]])

    ax.tick_params(length=8, width=2, labelsize=15)

    # ax.set_xlim([-1.0,1.2*obs_three_counts[15]])
    # ax.set_ylim([-1.0,1.2*obs_three_counts[15]])

    # ax.legend(bbox_to_anchor=(1.4, 0.00), loc=8, fontsize=15)

    plt.title("", fontsize=28)
    # plt.title("Three-neuron motifs\n (linear scale, motifs 15 & 16)", fontsize=18)
    plt.xlabel('Actual counts', fontsize=28)
    plt.ylabel('Sampled counts', fontsize=28)

    plt.xscale('symlog')
    plt.yscale('symlog')
    f.savefig("figures/{}.pdf".format(fname), bbox_inches='tight')    

def plot_cnr(actual_ucnr, config_ucnr, config_ucnr_std, ER_p_null, 
             figsize=(10,4), fname="temp", fontsizes=(20, 22, 14, 18)):
    f, ax = plt.subplots(1, 1, figsize=figsize)
    x = np.arange(10)

    sns.regplot(x=x, y=actual_ucnr["undir_perc"]*100, marker='x',
                label="Observation", color="#D8545D", ci=None, scatter_kws={'s':100, "zorder":10})
    sns.regplot(x=x, y=config_ucnr["undir_perc"]*100, label="Configuration Model", color='gray', ci=None, 
                scatter_kws={'s':100})
    # sns.regplot(x=x, y=erdos_cnr*100, label="Erdős-Rényi")
    ax.errorbar(x, config_ucnr["undir_perc"]*100, yerr=config_ucnr_std["undir_perc_std"]*100, 
                color='gray', fmt='o', alpha=0.8)
    plt.hlines((1-ER_p_null)*100, 0, 9, linestyle="--", label="Erdős-Rényi")

    ax.tick_params(length=8, width=1.5, labelsize=16)
    ax.set_ylabel('Percent connected', fontsize=fontsizes[0])
    ax.set_xlabel('Number of common neighbors to a pair', fontsize=fontsizes[1])
    x_labels = [0, 1,2,3,4,5,6,7,8,"9+"]
    ax.set_xticks(np.arange(len(x_labels)));
    ax.set_xticklabels(x_labels, fontsize=fontsizes[2]);

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines['left'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)

    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('axes', -0.05))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('axes', -0.05))

    ax.set_ylim(bottom=-1.4)

    ax.legend()
    ax.legend(fontsize=fontsizes[3])

    f.savefig("figures/{}.pdf".format(fname), bbox_inches='tight')  


def plot_common_successor(ft_config_ucnrs, actual_ucnr, figsize=(4,4), fname="temp"):
    f, ax = plt.subplots(1, 1, figsize=figsize)
    data = np.array([ft_config_ucnrs[i][1:,4].sum() for i in range(100)])

    z = (actual_ucnr["dir_sconn"][1:].sum() - data.mean()) / data.std()
    pv = scipy.stats.norm.sf(z)*2
    pc = (data>actual_ucnr["dir_sconn"][1:].sum()).sum()/100
    print(z, pv, pc)

    sns.distplot(data, color="gray", label="CFG")
    plt.vlines(actual_ucnr["dir_sconn"][1:].sum(), 0, 0.042, color='r', label="OBS")
    ax.tick_params(length=8, width=1.5, labelsize=15)
    ax.set_xlabel('# connected pairs with at least\n one common successor', fontsize=15)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines['left'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)

    ax.xaxis.set_ticks_position('bottom')
    # ax.spines['bottom'].set_position(('axes', -0.05))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('axes', -0.05))

    ax.legend(bbox_to_anchor=(0.20, 0.73), loc=8, fontsize=14)

    f.savefig("figures/{}.pdf".format(fname), bbox_inches='tight')

def plot_common_predecessor(ft_config_ucnrs, actual_ucnr, figsize=(4,4), fname="temp"):
    f, ax = plt.subplots(1, 1, figsize=figsize)
    data = np.array([ft_config_ucnrs[i][1:,5].sum() for i in range(100)])

    z = (actual_ucnr["dir_pconn"][1:].sum() - data.mean()) / data.std()
    pv = scipy.stats.norm.sf(z)*2
    pc = (data>actual_ucnr["dir_pconn"][1:].sum()).sum()/100
    print(z, pv, pc)

    sns.distplot(data, color="gray", label="CFG")
    plt.vlines(actual_ucnr["dir_pconn"][1:].sum(), 0, 0.02, color='r', label="OBS")
    ax.tick_params(length=8, width=1.5, labelsize=15)
    ax.set_xlabel('# connected pairs with at least\n one common predecessor', fontsize=15)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines['left'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)

    ax.xaxis.set_ticks_position('bottom')
    # ax.spines['bottom'].set_position(('axes', -0.05))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('axes', -0.05))

    ax.legend(bbox_to_anchor=(0.20, 0.73), loc=8, fontsize=14)

    f.savefig("figures/{}.pdf".format(fname), bbox_inches='tight')

    
def plot_two_neuron_counts_violin_adapt(obs_counts, er_counts, config_counts, subtitle='', ylim=None, figsize=(6,4), fname="temp"):
    """Compare observed, ER, and configuration two-neuron counts
    
    Args:
        obs_counts: pd.DataFrame from count_two_neuron_motifs
        er_counts: pd.DataFrame from compute_ER_two_neuron_motifs
        config_counts: pd.DataFrame from sample_config_two_neuron_motifs
    """
    p = obs_counts['actual_edges'] / obs_counts['potential_edges']
    print("connect probability = ", p)
    f, ax = plt.subplots(1, 1, figsize=figsize)
    x_labels = ['o   o', 'o$\Rightarrow$o', 'o$\Leftrightarrow$o']
    keys = ['null', 'uni', 'bi']
    
    data = [obs_counts[k] / er_counts[k] if er_counts[k] != 0 else 0 for k in keys]
    ax.scatter(np.arange(len(x_labels))+1, data, color='#D8545D', 
               s=150, marker='x', zorder=100, label="OBS/ER")
    
    mean = config_counts.mean()
    std = config_counts.std()    
    rdm = np.array([mean[k] / er_counts[k] if er_counts[k] != 0 else 0 for k in keys])
    err = [std[k] / er_counts[k] for k in keys]
#     ax.errorbar(np.arange(len(x_labels)), rdm, yerr=err, color='r', fmt='o', alpha=1)
    data = [config_counts[k] / er_counts[k] if er_counts[k] != 0 else 0 for k in keys]
    parts = ax.violinplot(data, widths=0.7,
                 showmeans=False, showmedians=False, showextrema=False)
    for pc in parts['bodies']:
#         pc.set_facecolor([0.0,0.5,0.9,0.5])
        pc.set_facecolor('#C0C0C0')
#         pc.set_edgecolor('black')
        pc.set_alpha(0.5)

    quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=1)
    p5, p95 = np.percentile(data, [5, 95], axis=1)
    inds = np.arange(1, len(medians) + 1)
    ax.scatter(inds, medians, marker='o', color='white', edgecolors='k', 
               s=80, zorder=3, label="CFG/ER")
    ax.vlines(inds, quartile1, quartile3, color='k', alpha=0.6, linestyle='-', lw=2.5)
    ax.vlines(inds, p5, p95, color='k', alpha=0.6, linestyle='-.', lw=1.0)
    mid, mlen = (p5+p95)/2, (p95-p5)/2
    eb = plt.errorbar(inds, (p5+p95)/2, yerr = mlen, ls='none', color='k', capsize=2, elinewidth=1.0) 
    eb[-1][0].set_linestyle('-.')
    
#     ax.legend(bbox_to_anchor=(0.28, 0.73), loc=8, fontsize=14)
    
    ax.plot(range(0,len(keys)+2), [1]*(len(keys)+2), '--', color='#745e4f', alpha=1.0)
    
    
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')

    # was: axes.spines['bottom'].set_position(('data',1.1*X.min()))
    ax.spines['bottom'].set_position(('axes', -0.05))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('axes', -0.05))
    
    ax.tick_params(length=8, width=1.5, labelsize=15);
    ax.set_xlim(left=0.5,right=3.5)
    ax.set_xticks(np.arange(1, len(x_labels)+1))
    ax.set_ylim(bottom=0)
    ax.set_xticklabels(x_labels, fontsize=20);
    ax.set_ylabel('Counts / ER Mean', fontsize=15);
    
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_title('{}'.format(subtitle), fontsize=20);
    f.savefig("figures/{}.pdf".format(fname), bbox_inches='tight')
    d = pd.DataFrame.from_dict({'obs': {k: obs_counts[k] for k in keys}, 
                                'er_exp': er_counts, 'config_exp': mean, 
                                'config_std': std}, orient='columns')
    d['config_z'] = (d['obs'] - d['config_exp']) / d['config_std']
    # print(d.loc[['null','uni','bi'],:].round({'obs': 0, 'null_exp': 0, 'config_exp': 0, 'config_std': 2, 'config_z': 1}))
    return d

def plot_two_neuron_counts_scatter_adapt(obs_two_counts, er_two_counts, config_two_counts, figsize=(5,5), fname="2-scatter"):
    f, ax = plt.subplots(1, 1, figsize=figsize)

    ax.scatter([obs_two_counts['null']], [er_two_counts['null']], 
                marker='^', color='#1f77b4', s=120, zorder=2, alpha=0.9, label='ER   o   o')
    ax.scatter(obs_two_counts['null'], config_two_counts['null'].mean(), 
                marker='o', color='#1f77b4', s=120, zorder=4, alpha=0.9, label='CFG o   o')
    ax.plot([obs_two_counts['null'], obs_two_counts['null']], 
            [er_two_counts['null'], config_two_counts['null'].mean()], 
             '-', color='#1f77b4', linewidth=1, zorder=4, alpha=1.0)

    ax.scatter([obs_two_counts['uni']], [er_two_counts['uni']], 
                marker='^', color='#ff7f0e', s=120, zorder=3, alpha=0.9, label='ER   o$\Rightarrow$o')
    ax.scatter(obs_two_counts['uni'], config_two_counts['uni'].mean(),  
                marker='o', color='#ff7f0e', s=120, zorder=4, alpha=0.9, label='CFG o$\Rightarrow$o')
    ax.plot([obs_two_counts['uni'], obs_two_counts['uni']], 
            [er_two_counts['uni'], config_two_counts['uni'].mean()], 
             '-', color='#ff7f0e', linewidth=1, zorder=4, alpha=1.0)

    ax.scatter([obs_two_counts['bi']], [er_two_counts['bi']], 
                marker='^', color='#2ca02c', s=120, zorder=4, alpha=0.9, label='ER   o$\Leftrightarrow$o')
    ax.scatter(obs_two_counts['bi'], config_two_counts['bi'].mean(), 
                marker='o', color='#2ca02c', s=120, zorder=4, alpha=0.9, label='CFG o$\Leftrightarrow$o')
    ax.plot([obs_two_counts['bi'], obs_two_counts['bi']], 
            [er_two_counts['bi'], config_two_counts['bi'].mean()], 
             '-', color= '#2ca02c', linewidth=1, zorder=4, alpha=1.0)

    x = np.linspace(0, obs_two_counts['null'], 100)
    ax.plot(x,x, '--k', zorder=1)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines['left'].set_linewidth(1.25)
    ax.yaxis.set_tick_params(width=2)
    ax.spines['bottom'].set_linewidth(1.25)
    ax.xaxis.set_tick_params(width=2)

    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')

    # was: axes.spines['bottom'].set_position(('data',1.1*X.min()))
    ax.spines['bottom'].set_position(('axes', -0.05))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('axes', -0.05))

    ax.set_xlim([10,1.2*obs_two_counts['null']])
    ax.set_ylim([10,1.2*obs_two_counts['null']])

    ax.tick_params(length=8, width=1.5, labelsize=15);
#   ax.legend(fontsize=15, bbox_to_anchor=(0.20, 0.45), loc=8)

    plt.xlabel('Actual counts', fontsize=20)
    plt.ylabel('Sampled counts', fontsize=20)

    plt.xscale('log')
    plt.yscale('log')
    f.savefig("figures/{}.pdf".format(fname), bbox_inches='tight')
    
def plot_three_neuron_counts_violin_adapt(obs_counts, er_counts, config_counts, subtitle='', ylim=None, fname="temp", figsize=(10,3)):
    f, ax = plt.subplots(1, 1, figsize=figsize)
    x_labels = list(obs_counts.keys()) #['o   o', 'o$\Rightarrow$o', 'o$\Leftrightarrow$o']
    obs_data = [obs_counts[k] / er_counts[k] for k in x_labels]
    # ax.bar(np.arange(len(x_labels))+1, obs_data, bottom=1, color='k', error_kw={'elinewidth': 1, 'capsize': 3})
    ax.scatter(np.arange(len(x_labels))+1, obs_data, color='#D8545D', s=100, marker='x', zorder=150, label="OBS/ER")
    
    mean = config_counts.mean()
    std = config_counts.std()    
    rdm = np.array([mean[k] / er_counts[k] for k in x_labels])
    err = [std[k] / er_counts[k] for k in x_labels]
#     ax.errorbar(np.arange(len(x_labels))+1, rdm, yerr=err, color='r', fmt='o', alpha=1)
    
    
    data = [config_counts[k] / er_counts[k] if er_counts[k] != 0 else 0 for k in x_labels]
    parts = ax.violinplot(data, widths=0.7,
                 showmeans=False, showmedians=False, showextrema=False)
    
    for pc in parts['bodies']:
        pc.set_facecolor('#C0C0C0')
#         pc.set_edgecolor('black')
        pc.set_alpha(0.5)
        
    quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=1)
    p5, p95 = np.percentile(data, [5, 95], axis=1)
    inds = np.arange(1, len(medians) + 1)
    ax.scatter(inds, medians, marker='o', color='w', edgecolors='k', s=80, zorder=3, label="CFG/ER")
    ax.vlines(inds, quartile1, quartile3, color='k', alpha=0.6, linestyle='-', lw=2.5)
#     ax.vlines(inds, p5, p95, color='k', alpha=0.6, linestyle='-.', lw=0.75)
    mid, mlen = (p5+p95)/2, (p95-p5)/2
    eb = plt.errorbar(inds, (p5+p95)/2, yerr = mlen, ls='none', color='k', capsize=2, elinewidth=1.0)
    eb[-1][0].set_linestyle('-.')
    
    ax.plot(range(0,len(x_labels)+2), [1]*(len(x_labels)+2), '--', color='#745e4f', alpha=1.0)
    
#     ax.legend(bbox_to_anchor=(0.12, 0.71), loc=8, fontsize=15)
    
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')

    # was: axes.spines['bottom'].set_position(('data',1.1*X.min()))
    ax.spines['bottom'].set_position(('axes', -0.05))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('axes', -0.05))
    
    ax.tick_params(length=8, width=1.5, labelsize=15);
    ax.set_xlim(left=0.5,right=16.5)
    ax.set_xticks(np.arange(len(x_labels))+1)
    # ax.set_xticklabels(x_labels, fontsize=20);
    
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_ylabel('Counts / ER Mean', fontsize=15);
#     ax.set_title('Three-neuron motifs{}' \
#                  '\nbootstrap iterations={}' \
#                  '\np(o   o)={:.3f}' \
#                  '\np(o$\Rightarrow$o)={:.3f}' \
#                  '\np(o$\Leftrightarrow$o)v={:.4f}'.format(subtitle,
#                                                           len(bootstrap),
#                                                           pair_pr['pr_null'], 
#                                                           pair_pr['pr_uni'], 
#                                                           pair_pr['pr_bi']), fontsize=18);
    ax.set_title(''.format(subtitle), fontsize=18);
    f.savefig("figures/{}.pdf".format(fname), bbox_inches='tight')
    d = pd.DataFrame.from_dict({'obs':obs_counts, 'er_exp': er_counts, 'config_exp': mean, 'config_std': std}, orient='columns')
    d['config_z'] = (d['obs'] - d['config_exp']) / d['config_std']
    # print(d.round({'obs': 0, 'null_exp': 0, 'config_exp': 0, 'config_std': 2, 'config_z': 1}))
    return d    



def plot_cnr_adapt(actual_ucnr, config_ucnr, config_ucnr_std, ER_p_null, 
             figsize=(10,4), fname="temp", fontsizes=(20, 22, 14, 18)):
    f, ax = plt.subplots(1, 1, figsize=figsize)
    x = np.arange(10)

    sns.regplot(x=x, y=actual_ucnr["undir_perc"]*100, marker='x',
                label="Observation", color="#D8545D", ci=None, scatter_kws={'s':100, "zorder":10})
    sns.regplot(x=x, y=config_ucnr["undir_perc"]*100, label="Configuration Model", color='gray', ci=None, 
                scatter_kws={'s':100})
    # sns.regplot(x=x, y=erdos_cnr*100, label="Erdős-Rényi")
    ax.errorbar(x, config_ucnr["undir_perc"]*100, yerr=config_ucnr_std["undir_perc_std"]*100, 
                color='gray', fmt='o', alpha=0.8)
    plt.hlines((1-ER_p_null)*100, 0, 9, linestyle="--", label="Erdős-Rényi")

    ax.tick_params(length=8, width=1.5, labelsize=16)
    ax.set_ylabel('Percent connected', fontsize=fontsizes[0])
    ax.set_xlabel('Number of common neighbors', fontsize=fontsizes[1])
    x_labels = [0, 1,2,3,4,5,6,7,8,"9+"]
    ax.set_xticks(np.arange(len(x_labels)));
    ax.set_xticklabels(x_labels, fontsize=fontsizes[2]);

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines['left'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)

    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('axes', -0.05))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('axes', -0.05))

    ax.set_ylim(bottom=-1.4)

#     ax.legend()
#     ax.legend(fontsize=fontsizes[3])

    f.savefig("figures/{}.pdf".format(fname), bbox_inches='tight')  
    
def plot_three_neuron_counts_scatter_adapt(obs_three_counts, er_three_counts, config_three_counts, figsize=(12,12), fname="temp"):
    f, ax = plt.subplots(1, 1, figsize=figsize)

    for i in range(1,17):
        ax.scatter([obs_three_counts[i]], [er_three_counts[i]], 
                    marker='o', color=cmaplist[i], s=250, zorder=2*i+1, alpha=0.75, label='ER   {} '.format(i))
        ax.scatter(obs_three_counts[i], config_three_counts[i], 
                    marker='*', color=cmaplist[i], s=250, zorder=2*i+2, alpha=0.75, label='CFG {} '.format(i))
    #     ax.errorbar(obs_three_counts[i], config_three_counts[i].mean(), yerr=config_three_counts[i].std(), 
    #                 color=cmaplist[i], alpha=1)
        ax.plot([obs_three_counts[i], obs_three_counts[i]], 
            [er_three_counts[i], config_three_counts[i]], 
             color=cmaplist[i], linewidth=1, alpha=0.75)

    x = np.linspace(-1, obs_three_counts[1], 100)
    ax.plot(x,x, '--k', lw=2, zorder=1)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines['left'].set_linewidth(2)
    ax.yaxis.set_tick_params(width=2)
    ax.spines['bottom'].set_linewidth(2)

    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')

    # was: axes.spines['bottom'].set_position(('data',1.1*X.min()))
    ax.spines['bottom'].set_position(('axes', -0.025))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('axes', -0.025))

    ax.set_xlim([-0.5,1.2*obs_three_counts[1]])
    ax.set_ylim([-0.5,1.2*obs_three_counts[1]])

    ax.tick_params(length=8, width=2, labelsize=15)

    # ax.set_xlim([-1.0,1.2*obs_three_counts[15]])
    # ax.set_ylim([-1.0,1.2*obs_three_counts[15]])

    # ax.legend(bbox_to_anchor=(1.4, 0.00), loc=8, fontsize=15)

    plt.title("", fontsize=28)
    # plt.title("Three-neuron motifs\n (linear scale, motifs 15 & 16)", fontsize=18)
    plt.xlabel('Actual counts', fontsize=28)
    plt.ylabel('Sampled counts', fontsize=28)

    plt.xscale('symlog')
    plt.yscale('symlog')
    f.savefig("figures/{}.pdf".format(fname), bbox_inches='tight')