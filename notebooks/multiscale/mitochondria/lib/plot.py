import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
import matplotlib.patches as patches
from matplotlib.ticker import ScalarFormatter
from scipy import stats
import seaborn as sns

from . import u
from . import modelling


COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]
PALETTE = {
    "axonal": COLORS[0],
    "somatic": "grey",
    "basal": COLORS[1],
    "apical": COLORS[3]
}

NUM_BINS = 30
DEFAULT_VX_VOL = (0.00716 * 0.00716 * 0.04)
AXIS_FONTSIZE = 25
TICK_FONTSIZE = 13
MU = "$\mathrm{\mu}$"
RAD_AXISLBL = "Mesh Radius ($\mathrm{\mu}$m)"  #noqa
VOL_AXISLBL = "Mitochondrion Volume ({MU}m$^3$)"  # noqa
AXISLBLS = {
    "mitoindex": f"Mitochondrial coverage factor",
    "dendmitoindex": f"Dendritic coverage factor",
    "mitocovfactor": f"Mitochondrial coverage factor",
    "dendmitocovfactor": f"Dendritic coverage factor",
    "basalmitocoverage": f"Basal coverage (%{MU}m)",
    "basalpathlength": f"Basal path length ({MU}m)",
    "basalsyncount": "Basal synapse count",
    "basalsyndensity": f"Basal synapse density (/{MU}m)",
    "basalmitoindex": f"Basal coverage factor",
    "basalmitocovfactor": f"Basal coverage factor",
    "basal %pl within 20um": f"Basal %pl within 20{MU}m",
    "dendsyndensity": f"Dendritic synapse density (/{MU}m)",
    "vol": f"Mitochondrion volume ({MU}m$^3$)",
    "somavol": f"Somatic volume ({MU}m$^3$)",
    "nucvol": f"Nuclear volume ({MU}m$^3$)",
    "somamitodensity": f"Somatic coverage factor",
    "somamitovol": f"Somatic mitochondrial volume (/{MU}m$^3$)",
    "somasyndensity": f"Somatic synapse density (/{MU}m$^2$)",
    "depth": "Cortical depth (a.u. for now)",
    "apicalmitocoverage": f"Apical coverage",
    "apicalpathlength": f"Apical path length ({MU}m)",
    "apicalsyncount": "Apical synapse count",
    "apicalsyndensity": f"Apical synapse density (/{MU}m)",
    "apicalmitocovfactor": "Apical coverage factor",
    "apical %pl within 20um": f"Apical %pl within 20{MU}m"
}
OUTLIERID = 648518346349534079

matplotlib.rc('font', family='sans-serif')
matplotlib.rc('font', serif='Helvetica')

def make_patches(key_to_color):
    return [patches.Patch(facecolor=c, label=key, edgecolor='k')
            for (key, c) in key_to_color.items()]


def default_style():
    global COLORS
    plt.style.use("default")
    COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]


def dark_background():
    global COLORS
    plt.style.use("dark_background")
    COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]


#==========================================================
# Panel D - Volume comparison by compartment

def panelD(df, idsnearsoma=None):
    matplotlib.rc('font', family='sans-serif')
    matplotlib.rc('font', serif='Helvetica')

    df.loc[df.nodelbl == 0, "compartment"] = 's'
    df.loc[df.nodelbl == 1, "compartment"] = 'a'
    df.loc[df.nodelbl == 2, "compartment"] = 'd'
    df.loc[df.nodelbl == 3, "compartment"] = 'p'
    df.loc[df.nodelbl == 4, "compartment"] = 'o'
    if idsnearsoma is not None:
        df.loc[df.index.isin(idsnearsoma), "compartment"] = 'o'

    plt.tick_params(axis='y', labelsize=15)
    plt.tick_params(axis='x', labelsize=20)

    plotdf = plot_vol_by_comp_box(df)

    #plt.xticks(rotation=45)
    format_axes()

    plt.yscale("log")

    return plotdf


def plot_vol_by_comp_box(
    df, varname="mito_vx", vx_vol=DEFAULT_VX_VOL, palette=PALETTE):

    df = df.copy()
    data = df[varname] * vx_vol
    df["vol"] = data

    #Ordering categories and removing "other" label
    df = pd.concat((
             df[df.compartment == 'a'],
             df[df.compartment == 's'],
             df[df.compartment == 'p'],
             df[df.compartment == 'd']))

    df.loc[df.compartment == 's', "compartment"] = "somatic"
    df.loc[df.compartment == 'a', "compartment"] = "axonal"
    if sum(df.compartment == 'p') > 0:  # apical labeled separately
        df.loc[df.compartment == 'p', "compartment"] = "apical"
        df.loc[df.compartment == 'd', "compartment"] = "basal"
    else:
        df.loc[df.compartment == 'd', "compartment"] = "dendritic"

    aux_args = dict(width=0.7, fliersize=0, whis=(5, 95), zorder=10)
    palette = {compartment : "white" for compartment in np.unique(df.compartment)}
    if palette is None:
        ax = sns.boxplot(x="compartment", y="vol", data=df, **aux_args)
    else:
        ax = sns.boxplot(x="compartment", y="vol", **aux_args,
                         palette=palette, data=df)
    make_black_and_white(ax)

    plt.ylabel(AXISLBLS["vol"], fontsize=AXIS_FONTSIZE)
    plt.xlabel("")

    return df


def make_black_and_white(ax):
    # iterate over boxes
    for i,box in enumerate(ax.artists):
        box.set_edgecolor('black')
        box.set_facecolor('white')
    
        # iterate over whiskers and median lines
        for j in range(6*i,6*(i+1)):
             ax.lines[j].set_color('black')


#==========================================================
# Panel E - Dendritic coverage by somatic density

def panelE(covdf, outlierids=[OUTLIERID],
           xcol="somamitodensity",
           ycol="basalmitocovfactor"):
    return plot_scatter_w_outlier(
               covdf, xcol, ycol,
               outlierids=outlierids)


def plot_scatter_w_outlier(covdf, xcolname, ycolname,
                           xlbl=None, ylbl=None, outlierids=None, color='k'):

    if xlbl is None:
        assert xcolname in AXISLBLS, "xcolname not in lookup w/ no default"
        xlbl = AXISLBLS[xcolname]

    if ylbl is None:
        assert ycolname in AXISLBLS, "ycolname not in lookup w/ no default"
        ylbl = AXISLBLS[ycolname]

    r, p, fit = plot_scatter(covdf[xcolname], covdf[ycolname],
                             xlbl, ylbl, color=color)

    if outlierids is not None:
        highlight_outliers(covdf, outlierids, xcolname, ycolname)

    return r, p, fit


def plot_scatter(x, y, xlbl=None, ylbl=None,
                 show_r=False, linear_fit=True,
                 size=None, marker=True, color='k'):
    matplotlib.rc('font', family='sans-serif')
    matplotlib.rc('font', serif='Helvetica')

    #plt.errorbar(x, y, fmt='.', c='k')
    if marker is True:
        sns.scatterplot(x, y, color=color, alpha=0.8)
    else:
        plt.scatter(x, y, color=color, s=size, marker=marker)

    if linear_fit:
        fit = plot_linear_fit(x, y)
    else:
        fit = None

    if len(x) > 1:
        r, p = stats.pearsonr(x, y)
    else:
        r, p = None, None

    format_axes()
    if xlbl is not None:
        plt.xlabel(xlbl, fontsize=AXIS_FONTSIZE)
    if ylbl is not None:
        plt.ylabel(ylbl, fontsize=AXIS_FONTSIZE)

    return r, p, fit


def plot_linear_fit(X, y, lw=3, color=COLORS[3]):
    fit, X = modelling.fit_linear_model(X, y)

    xlbl = [colname for colname in X.columns if colname != "const"]
    assert len(xlbl) == 1, "multiple potential x-axis variables"
    x = X[xlbl[0]]
    dummy_preds = fit.get_prediction(X)
    confint = dummy_preds.summary_frame(alpha=0.05)

    plt.plot(x, fit.fittedvalues, color='r')
    plt.fill_between(
        x, confint["obs_ci_lower"], confint["obs_ci_upper"],
        alpha=0.2, color="gray")

    return fit


def format_axes(ax=None):
    ax = plt.gca() if ax is None else ax

    ax.spines["bottom"].set_linewidth(1.5)
    ax.spines["left"].set_linewidth(1.5)
    ax.spines["right"].set_linewidth(0)
    ax.spines["top"].set_linewidth(0)

    ax.tick_params(length=10, width=1.5, labelsize=18)
    ax.tick_params(length=5, width=1, which="minor")
    

def highlight_outliers(covdf, outlierids, xcolname, ycolname):
    x = covdf[xcolname][covdf.cellid.isin(outlierids)]
    y = covdf[ycolname][covdf.cellid.isin(outlierids)]

    plot_scatter(x, y, linear_fit=False, marker='x', size=120)


#==========================================================
# Panel F - Dendritic mitochondrion coverage vs. synapse density

def panelF(covdf, outlierids=[OUTLIERID],
           xcol="basalsyndensity",
           ycol="basalmitocovfactor"):
    return plot_scatter_w_outlier(
               covdf, xcol, ycol,
               outlierids=outlierids)


#==========================================================
# Panel G - Somatic synapse density vs. mitochondrion density

def panelG(covdf, outlierids=[OUTLIERID],
           xcol="somasyndensity",
           ycol="somamitodensity"):
    return plot_scatter_w_outlier(
               covdf, "somasyndensity", "somamitodensity",
               outlierids=outlierids)


#==========================================================
# Extended Data Fig 4

def extPanel4B(covdf, varname="mitocovfactor"):
    
    fig, subplots = plt.subplots(3, 1, sharex=True,
                                 gridspec_kw = {'wspace':0, 'hspace':0.1},
                                 figsize=(8, 7))

    # apical dendrites
    ax = subplots[0]
    plot_mean_gradient(covdf[covdf.nodelbl == 3], textlbl="apical",
                       color='k', varname=varname, ax=ax)
    ax.set_xticks([])
    ax.tick_params(axis='both', which='major', labelsize=15)

    # basal dendrites
    ax = subplots[1]
    plot_mean_gradient(covdf[covdf.nodelbl == 2], textlbl="basal",
                       color='k', varname=varname, ax=ax)
    ax.set_xticks([])
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_ylabel(AXISLBLS[varname], fontsize=AXIS_FONTSIZE)
    
    # axons
    ax = subplots[2]
    plot_mean_gradient(covdf[covdf.nodelbl == 1], textlbl="axonal",
                       color='k', varname=varname, ax=ax)
    ax.set_yticks([], True)
    ax.set_yticks([0, 0.15, 0.3, 0.45])
    #for axis in [ax.yaxis]:
    #    formatter = ScalarFormatter()
    #    formatter.set_scientific(False)
    #    axis.set_major_formatter(formatter)
    ax.tick_params(axis='both', which='major', labelsize=15)

    plt.xlim(0, 355)
    plt.xlabel("Distance from soma node ($\mathrm{\mu}$m)", fontsize=AXIS_FONTSIZE)


def plot_mean_gradient(df, varname="coverage", textlbl=None,
                       color=None, ax=None):
    ax = plt.gca() if ax is None else ax

    bin_models = modelling.fit_models_to_bins(df, "norm", varname=varname)
    ns = df.groupby("bin_center")[varname].count().values

    plot_gradient(np.unique(df.bin_center), bin_models, ns, "norm", errtype="stderr",
                  tick_spacing=100, color=color, ax=ax)

    ax.text(0.95, 0.8, textlbl,
            horizontalalignment="right", transform=ax.transAxes, fontsize=15)


def plot_gradient(
    centers, models, ns, modeltype, color=None, tick_spacing=50,
    sns_layer=False, ax=None, errtype="stderr", saturation=0.75, **kwargs):

    if sns_layer:
        centers = sns_centers(centers)

    if saturation != 1:
        color = sns.desaturate(color, saturation)

    plot_errorbars(centers, models, ns, color=color, ax=ax,
                   errtype=errtype, **kwargs)

    if tick_spacing is not None:
        binwidth = centers[1] - centers[0]
        upper_bound = centers[-1] + binwidth/2.
        num_ticks = upper_bound // tick_spacing + 1  # +1 for 0
        ticks = (np.arange(num_ticks) * tick_spacing).astype(int)

        plt.xticks(ticks)
        plt.xlim(0, upper_bound)


def plot_errorbars(xs, models, ns, errtype="stderr",
                   color=None, ax=None, **kwargs):
    ax = plt if ax is None else ax
    color = 'k' if color is None else color

    mvs = [model.stats() for model in models]
    means = np.array([mv[0] for mv in mvs])
    stddevs = np.array([np.sqrt(mv[1]) for mv in mvs])

    if errtype == "stderr":
        err = stddevs / np.sqrt(ns)
    elif errtype == "sqrt(n)":
        err = np.sqrt(means * ns) / ns
    elif errtype == "zero":
        err = np.zeros(means.shape)

    mask = ns > 10
    xs = xs[mask]
    means = means[mask]
    err = err[mask]

    ax.errorbar(xs, means, yerr=err,
                fmt="o", lw=3, color=color, capsize=10, capthick=3,
                zorder=10, **kwargs)


def format_distdf(
    distdf, binwidth=50, dist_scale=1/1000., vx_vol=DEFAULT_VX_VOL,
    varname="mitovols"):
    validdf = distdf[((~np.isinf(distdf.nodedists)) &
                      (distdf.nodedists > 0))]

    if dist_scale != 1:
        validdf.loc[:, "nodedists"] = validdf["nodedists"] * dist_scale

    if varname == "mitovols":
        validdf.loc[:, "mitovols"] = validdf["mitovols"] * vx_vol

    validdf.loc[:, "bin_center"] = ((validdf["nodedists"] // binwidth)
                                    * binwidth + binwidth // 2)

    bins = np.unique(validdf["bin_center"])
    counts = np.array(validdf.groupby("bin_center")[varname].count())

    return validdf, bins, counts


def extPanel4C(distdf):
    
    fig, subplots = plt.subplots(3, 1, sharex=True,
                                 gridspec_kw = {'wspace':0, 'hspace':0.1},
                                 figsize=(8, 7))

    # apical dendrites
    ax = subplots[0]
    plot_mean_gradient(distdf[distdf.nodelbl == 3], textlbl="apical",
                       color='k', ax=ax, varname="vol")
    ax.set_xticks([])
    ax.tick_params(axis='both', which='major', labelsize=15)

    # basal dendrites
    ax = subplots[1]
    plot_mean_gradient(distdf[distdf.nodelbl == 2], textlbl="basal",
                       color='k', ax=ax, varname="vol")
    ax.set_xticks([])
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_ylabel(f"Mitochondrion volume ({MU}m$^3$)", fontsize=AXIS_FONTSIZE)
    
    # axons
    ax = subplots[2]
    plot_mean_gradient(distdf[distdf.nodelbl == 1], textlbl="axonal",
                       color='k', ax=ax, varname="vol")
    ax.set_yticks([0.1, 0.05])
    ax.set_yticks([], True)
    for axis in [ax.yaxis]:
        formatter = ScalarFormatter()
        formatter.set_scientific(False)
        axis.set_major_formatter(formatter)
    ax.tick_params(axis='both', which='major', labelsize=15)

    plt.xlim(0, 355)
    plt.xlabel("Distance from soma node ($\mathrm{\mu}$m)", fontsize=AXIS_FONTSIZE)


# copied from StackOverflow
# https://stackoverflow.com/questions/11517986/indicating-the-statistically-significant-difference-in-bar-graph  # noqa
def barplot_annotate_brackets(
    num1, num2, data, center, height, yerr=None,
    dh=.05, barh=.05, fs=None, maxasterix=None):
    """ 
    Annotate barplot with p-values.

    :param num1: number of left bar to put bracket over
    :param num2: number of right bar to put bracket over
    :param data: string to write or number for generating asterixes
    :param center: centers of all bars (like plt.bar() input)
    :param height: heights of all bars (like plt.bar() input)
    :param yerr: yerrs of all bars (like plt.bar() input)
    :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
    :param barh: bar height in axes coordinates (0 to 1)
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    """

    if type(data) is str:
        text = data
    else:
        # * is p < 0.05
        # ** is p < 0.005
        # *** is p < 0.0005
        # etc.
        text = ''
        p = .05

        while data < p:
            text += '*'
            p /= 10.

            if maxasterix and len(text) == maxasterix:
                break

        if len(text) == 0:
            text = 'n. s.'

    lx, ly = center[num1], height[num1]
    rx, ry = center[num2], height[num2]

    if yerr:
        ly += yerr[num1]
        ry += yerr[num2]

    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh *= (ax_y1 - ax_y0)
    barh *= (ax_y1 - ax_y0)

    y = max(ly, ry) + dh

    barx = [lx, lx, rx, rx]
    bary = [y, y+barh, y+barh, y]
    mid = ((lx+rx)/2, y+barh)

    plt.plot(barx, bary, c='black')

    kwargs = dict(ha='center', va='bottom', fontsize=20)
    if fs is not None:
        kwargs['fontsize'] = fs

    plt.text(*mid, text, **kwargs)


def extPanel4D(covdf, outlierids=[OUTLIERID],
               xcol="somasyndensity",
               ycol="basalsyndensity"):
    return plot_scatter_w_outlier(
               covdf, xcol, ycol,
               outlierids=outlierids)


#==========================================================
# Extended Data Fig 5


def extPanel5C(covdf, outlierids=[OUTLIERID],
               xcol="somamitodensity",
               ycol="basalmitocovfactor"):
    plotdf = covdf[~(covdf.cellid.isin(outlierids))]
    return panelE(plotdf, outlierids=None, xcol=xcol, ycol=ycol)


def extPanel5D(covdf, outlierids=[OUTLIERID],
               xcol="basalsyndensity",
               ycol="basalmitocovfactor"):
    plotdf = covdf[~(covdf.cellid.isin(outlierids))]
    return panelF(plotdf, outlierids=None, xcol=xcol, ycol=ycol)


def extPanel5E(covdf, outlierids=[OUTLIERID],
               xcol="somasyndensity",
               ycol="somamitodensity"):
    plotdf = covdf[~(covdf.cellid.isin(outlierids))]
    return panelG(plotdf, outlierids=None, xcol=xcol, ycol=ycol)
