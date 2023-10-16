import random
import scipy
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import rc, rcParams
from mpl_toolkits.mplot3d import Axes3D
rc('text', usetex=False)
rc('axes', linewidth=1.5)
rc('font', weight='regular')

cmap = plt.cm.tab20
cmaplist = [cmap(17*i % 20) for i in range(cmap.N)]

def hist_outline(dataIn, *args, **kwargs):
    (histIn, binsIn) = np.histogram(dataIn, *args, **kwargs)
    histIn = np.insert(histIn, 0, 0)
    histIn = np.append(histIn, 0)
    bins = [t for b in zip(binsIn, binsIn) for t in b]
    data = [t for b in zip(histIn[:-1], histIn[1:]) for t in b]
    return (bins, data)

def overlay_hists_log(x, y, xlabels, ylabels, title, density=False, nbins=8, fname="temp"):
    f, axes = plt.subplots(2, 2, figsize=(12,8)) 
    for axs, d, dlabels in zip(axes, (x,y), (xlabels, ylabels)):
        ex = -0.2
        colors = ["#D8545D", "k"]
        lss = ["-", "--"]
        bins = np.linspace(0, np.max([np.max(di) for di in d]), nbins)
        bins_exp = np.concatenate(([ex], np.e**np.linspace(0, np.log(np.max([np.max(di) for di in d])), nbins)))
        for l, di, c, ls in zip(dlabels, d, colors, lss):
            nzx = di.astype(np.float64)
            nzx[nzx == 0] = np.e**ex
            dx, dy = hist_outline(di, bins=bins, density=density)
            dx_exp, dy_exp = hist_outline(nzx, bins=bins_exp, density=density)
            if (l == "Expected from ER") or (l == "Expected from gER") \
                or (l == "Expected from PM") and not density:
                dy = np.array(dy) / 100
                dy_exp = np.array(dy_exp) / 100
            axs[0].plot(dx, dy, label=l, c=c, linestyle=ls)
            axs[1].plot(dx_exp, np.array(dy_exp), label=l, c=c, linestyle=ls)
        axs[0].legend(loc='best', fontsize=14)
        axs[1].legend(loc='best', fontsize=14)
        axs[1].set_xscale('log')
        axs[1].set_yscale('log')
        axs[0].tick_params(length=8, width=1.5, labelsize=14)
        axs[1].tick_params(length=8, width=1.5, labelsize=14)
    axes[0][0].set_xlabel('In degree', fontsize=20);
    axes[0][1].set_xlabel('In degree', fontsize=20);
    axes[0][0].set_ylabel('Number of cells', fontsize=20);
    axes[0][1].set_ylabel('Number of cells', fontsize=20);
    axes[1][0].set_xlabel('Out degree', fontsize=20)
    axes[1][1].set_xlabel('Out degree', fontsize=20);
    axes[1][0].set_ylabel('Number of cells', fontsize=20);
    axes[1][1].set_ylabel('Number of cells', fontsize=20);
    f.suptitle(title, fontsize=20);
    f.tight_layout(rect=[0,0,1,0.94])
    f.savefig("figures/{}.pdf".format(fname), bbox_inches='tight')
    

def plot_soma_dist(locs, axls, th, figname="tmp"):
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111, projection='3d')

    for n in locs.keys():
        (x, y, z) = locs[n] / 1000
        if axls[n] <= th:
            ax.scatter(x, z, y, c='gray', marker='^', alpha=0.8)
        else:
            ax.scatter(x, z, y, c='#D8545D', marker='o')

    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    ax.set_xlabel('   x (μm)', fontsize=16)
    ax.set_ylabel('z (μm)', fontsize=16)
    ax.set_zlabel('y (μm)', fontsize=16)

    ax.view_init(-160, -20)

    plt.show()
    fig.savefig("figures/{}.pdf".format(figname), bbox_inches='tight')

    
def plot_somadist_p(prox, avg_th, fname="tmp"):
    sns.set_context('talk')
    # calculate 95% confidence interval
    ci_th100 = 1.96 * np.sqrt(prox["p_connect"]*(1-prox["p_connect"])/
                    prox["p_connect_std"])
    
    f, ax = plt.subplots(1, 1, figsize=(6,4))
    ax.plot(prox["soma_dist"], prox["p_connect"], c="#D8545D")
    ax.fill_between(prox["soma_dist"], 
                    prox["p_connect"]-ci_th100, 
                    prox["p_connect"]+ci_th100, 
                    color="#D8545D", alpha = 0.1)
    ax.set_xlabel("Soma distance (μm) ± {}μm".format(avg_th))
    ax.set_ylabel("Connection probability")
    f.savefig("figures/{}.pdf".format(fname), bbox_inches='tight')
    plt.show()    
    
    
def plot_two_neuron_counts_violin(obs_counts, er_counts, config_counts, subtitle='', ylim=None, figsize=(6,4), er_name="ER", cfg_name="CFG", fname="temp"):
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
               s=150, marker='x', zorder=100, label="OBS/{}".format(er_name))
    
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
               s=80, zorder=3, label="{}/{}".format(cfg_name, er_name))
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
    ax.set_ylabel('Counts / {} Mean Counts'.format(er_name), fontsize=15);
    
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_title('{}'.format(subtitle), fontsize=20);
    f.savefig("figures/{}.pdf".format(fname), bbox_inches='tight')
    if cfg_name != "PM":
        d = pd.DataFrame.from_dict({'obs': {k: obs_counts[k] for k in keys}, 
                                    'er_exp': er_counts, 'config_exp': mean, 
                                    'config_std': std}, orient='columns')
        d['proximity_z'] = (d['obs'] - d['config_exp']) / d['config_std']
    elif cfg_name == "PM":
        d = pd.DataFrame.from_dict({'obs': {k: obs_counts[k] for k in keys}, 
                                    'er_exp': er_counts, 'proximity_exp': mean, 
                                    'proximity_std': std}, orient='columns')
        d['proximity_z'] = (d['obs'] - d['proximity_exp']) / d['proximity_std']
    # print(d.loc[['null','uni','bi'],:].round({'obs': 0, 'null_exp': 0, 'config_exp': 0, 'config_std': 2, 'config_z': 1}))
    return d

def plot_two_neuron_counts_scatter(obs_two_counts, er_two_counts, config_two_counts, figsize=(5,5), er_name="ER", cfg_name="CFG", fname="2-scatter"):
    f, ax = plt.subplots(1, 1, figsize=figsize)

    ax.scatter([obs_two_counts['null']], [er_two_counts['null']], 
                marker='^', color='#1f77b4', s=120, zorder=2, alpha=0.9, label='{}   o   o'.format(er_name))
    ax.scatter(obs_two_counts['null'], config_two_counts['null'].mean(), 
                marker='o', color='#1f77b4', s=120, zorder=4, alpha=0.9, label='{} o   o'.format(cfg_name))
    ax.plot([obs_two_counts['null'], obs_two_counts['null']], 
            [er_two_counts['null'], config_two_counts['null'].mean()], 
             '-', color='#1f77b4', linewidth=1, zorder=4, alpha=1.0)

    ax.scatter([obs_two_counts['uni']], [er_two_counts['uni']], 
                marker='^', color='#ff7f0e', s=120, zorder=3, alpha=0.9, label='{}   o$\Rightarrow$o'.format(er_name))
    ax.scatter(obs_two_counts['uni'], config_two_counts['uni'].mean(),  
                marker='o', color='#ff7f0e', s=120, zorder=4, alpha=0.9, label='{} o$\Rightarrow$o'.format(cfg_name))
    ax.plot([obs_two_counts['uni'], obs_two_counts['uni']], 
            [er_two_counts['uni'], config_two_counts['uni'].mean()], 
             '-', color='#ff7f0e', linewidth=1, zorder=4, alpha=1.0)

    ax.scatter([obs_two_counts['bi']], [er_two_counts['bi']], 
                marker='^', color='#2ca02c', s=120, zorder=4, alpha=0.9, label='{}   o$\Leftrightarrow$o'.format(er_name))
    ax.scatter(obs_two_counts['bi'], config_two_counts['bi'].mean(), 
                marker='o', color='#2ca02c', s=120, zorder=4, alpha=0.9, label='{} o$\Leftrightarrow$o'.format(cfg_name))
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

def plot_three_neuron_counts_violin(obs_counts, er_counts, config_counts, subtitle='', ylim=None, er_name="ER", cfg_name="CFG", fname="temp", figsize=(10,3)):
    f, ax = plt.subplots(1, 1, figsize=figsize)
    x_labels = list(obs_counts.keys()) #['o   o', 'o$\Rightarrow$o', 'o$\Leftrightarrow$o']
    obs_data = [obs_counts[k] / er_counts[k] for k in x_labels]
    # ax.bar(np.arange(len(x_labels))+1, obs_data, bottom=1, color='k', error_kw={'elinewidth': 1, 'capsize': 3})
    ax.scatter(np.arange(len(x_labels))+1, obs_data, color='#D8545D', s=100, marker='x', zorder=150, label="OBS/{}".format(er_name))
    
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
    ax.scatter(inds, medians, marker='o', color='w', edgecolors='k', s=80, zorder=3, label="{}/{}".format(cfg_name, er_name))
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
    ax.set_ylabel('Counts / {} Mean Counts'.format(er_name), fontsize=15);
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
    if cfg_name != "PM":
        d = pd.DataFrame.from_dict({'obs':obs_counts, 'er_exp': er_counts, 
                                    'config_exp': mean, 'config_std': std}, orient='columns')
        d['config_z'] = (d['obs'] - d['config_exp']) / d['config_std']
    elif cfg_name == "PM":
        d = pd.DataFrame.from_dict({'obs':obs_counts, 'er_exp': er_counts, 
                                    'proximity_exp': mean, 'proximity_std': std}, orient='columns')
        d['proximity_z'] = (d['obs'] - d['proximity_exp']) / d['proximity_std']
    # print(d.round({'obs': 0, 'null_exp': 0, 'config_exp': 0, 'config_std': 2, 'config_z': 1}))
    return d


def plot_three_neuron_counts_diff_violin(obs_counts, er_counts, config_counts, subtitle='', ylim=None, er_name="ER", cfg_name="CFG", fname="temp", figsize=(10,3)):
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
    ax.set_ylabel('Counts - {} Mean Counts'.format(cfg_name), fontsize=15);
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

def plot_three_neuron_counts_scatter(obs_three_counts, er_three_counts, config_three_counts, figsize=(12,12), er_name="ER", cfg_name="CFG", fname="temp"):
    f, ax = plt.subplots(1, 1, figsize=figsize)

    for i in range(1,17):
        ax.scatter([obs_three_counts[i]], [er_three_counts[i]], 
                    marker='^', color=cmaplist[i], s=250, zorder=2*i+1, alpha=0.75, label='{}   {} '.format(er_name, i))
        ax.scatter(obs_three_counts[i], config_three_counts[i].mean(), 
                    marker='o', color=cmaplist[i], s=250, zorder=2*i+2, alpha=0.75, label='{} {} '.format(cfg_name, i))
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
             figsize=(10,4), er_name="Erdős-Rényi", cfg_name="Configuration Model", fname="temp", fontsizes=(20, 22, 14, 18)):
    f, ax = plt.subplots(1, 1, figsize=figsize)
    x = np.arange(10)

    sns.regplot(x=x, y=actual_ucnr["undir_perc"]*100, marker='x',
                label="Observation", color="#D8545D", ci=None, scatter_kws={'s':100, "zorder":10})
    sns.regplot(x=x, y=config_ucnr["undir_perc"]*100, label=cfg_name, color='gray', ci=None, 
                scatter_kws={'s':100})
    # sns.regplot(x=x, y=erdos_cnr*100, label="Erdős-Rényi")
    ax.errorbar(x, config_ucnr["undir_perc"]*100, yerr=config_ucnr_std["undir_perc_std"]*100, 
                color='gray', fmt='o', alpha=0.8)
    plt.hlines((1-ER_p_null)*100, 0, 9, linestyle="--", label=er_name)

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

    
def plot_two_neuron_counts_violin_adapt(obs_counts, er_counts, config_counts, subtitle='', ylim=None, figsize=(6,4), er_name="ER", cfg_name="CFG", fname="temp"):
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
               s=150, marker='x', zorder=100, label="OBS/{}".format(er_name))
    
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
               s=80, zorder=3, label="{}/{}".format(cfg_name, er_name))
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
    ax.set_ylabel('Counts / {} Mean'.format(er_name), fontsize=15);
    
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

def plot_two_neuron_counts_scatter_adapt(obs_two_counts, er_two_counts, config_two_counts, figsize=(5,5), er_name="ER", cfg_name="CFG", fname="2-scatter"):
    f, ax = plt.subplots(1, 1, figsize=figsize)

    ax.scatter([obs_two_counts['null']], [er_two_counts['null']], 
                marker='^', color='#1f77b4', s=120, zorder=2, alpha=0.9, label='{}   o   o'.format(er_name))
    ax.scatter(obs_two_counts['null'], config_two_counts['null'].mean(), 
                marker='o', color='#1f77b4', s=120, zorder=4, alpha=0.9, label='{} o   o'.format(cfg_name))
    ax.plot([obs_two_counts['null'], obs_two_counts['null']], 
            [er_two_counts['null'], config_two_counts['null'].mean()], 
             '-', color='#1f77b4', linewidth=1, zorder=4, alpha=1.0)

    ax.scatter([obs_two_counts['uni']], [er_two_counts['uni']], 
                marker='^', color='#ff7f0e', s=120, zorder=3, alpha=0.9, label='{}   o$\Rightarrow$o'.format(er_name))
    ax.scatter(obs_two_counts['uni'], config_two_counts['uni'].mean(),  
                marker='o', color='#ff7f0e', s=120, zorder=4, alpha=0.9, label='{} o$\Rightarrow$o'.format(cfg_name))
    ax.plot([obs_two_counts['uni'], obs_two_counts['uni']], 
            [er_two_counts['uni'], config_two_counts['uni'].mean()], 
             '-', color='#ff7f0e', linewidth=1, zorder=4, alpha=1.0)

    ax.scatter([obs_two_counts['bi']], [er_two_counts['bi']], 
                marker='^', color='#2ca02c', s=120, zorder=4, alpha=0.9, label='{}   o$\Leftrightarrow$o'.format(er_name))
    ax.scatter(obs_two_counts['bi'], config_two_counts['bi'].mean(), 
                marker='o', color='#2ca02c', s=120, zorder=4, alpha=0.9, label='{} o$\Leftrightarrow$o'.format(cfg_name))
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
    
def plot_three_neuron_counts_violin_adapt(obs_counts, er_counts, config_counts, subtitle='', ylim=None, fname="temp", er_name="ER", cfg_name="CFG", figsize=(10,3)):
    f, ax = plt.subplots(1, 1, figsize=figsize)
    x_labels = list(obs_counts.keys()) #['o   o', 'o$\Rightarrow$o', 'o$\Leftrightarrow$o']
    obs_data = [obs_counts[k] / er_counts[k] for k in x_labels]
    # ax.bar(np.arange(len(x_labels))+1, obs_data, bottom=1, color='k', error_kw={'elinewidth': 1, 'capsize': 3})
    ax.scatter(np.arange(len(x_labels))+1, obs_data, color='#D8545D', s=100, marker='x', zorder=150, label="OBS/{}".format(er_name))
    
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
    ax.scatter(inds, medians, marker='o', color='w', edgecolors='k', s=80, zorder=3, label="{}/{}".format(cfg_name, er_name))
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
    ax.set_ylabel('Counts / {} Mean'.format(er_name), fontsize=15);
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
             figsize=(10,4), er_name="Erdős-Rényi", cfg_name="Configuration Model",  fname="temp", fontsizes=(20, 22, 14, 18)):
    f, ax = plt.subplots(1, 1, figsize=figsize)
    x = np.arange(10)

    sns.regplot(x=x, y=actual_ucnr["undir_perc"]*100, marker='x',
                label="Observation", color="#D8545D", ci=None, scatter_kws={'s':100, "zorder":10})
    sns.regplot(x=x, y=config_ucnr["undir_perc"]*100, label=cfg_name, color='gray', ci=None, 
                scatter_kws={'s':100})
    # sns.regplot(x=x, y=erdos_cnr*100, label="Erdős-Rényi")
    ax.errorbar(x, config_ucnr["undir_perc"]*100, yerr=config_ucnr_std["undir_perc_std"]*100, 
                color='gray', fmt='o', alpha=0.8)
    plt.hlines((1-ER_p_null)*100, 0, 9, linestyle="--", label=er_name)

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
    
def plot_three_neuron_counts_scatter_adapt(obs_three_counts, er_three_counts, config_three_counts, figsize=(12,12), er_name="ER", cfg_name="CFG", fname="temp"):
    f, ax = plt.subplots(1, 1, figsize=figsize)

    for i in range(1,17):
        ax.scatter([obs_three_counts[i]], [er_three_counts[i]], 
                    marker='o', color=cmaplist[i], s=250, zorder=2*i+1, alpha=0.75, label='{}   {} '.format(er_name, i))
        ax.scatter(obs_three_counts[i], config_three_counts[i], 
                    marker='*', color=cmaplist[i], s=250, zorder=2*i+2, alpha=0.75, label='{} {} '.format(cfg_name, i))
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