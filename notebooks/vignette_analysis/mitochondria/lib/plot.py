import warnings
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn
import seaborn as sns
from scipy import stats

from . import modelling, u

COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]
# https://coolors.co/f42c04-107e7d-04080f-eca72c-8963ba
COOLORMAP = ["#f42c04", "#107e7d", "#04080f", "#eca72c", "#8963ba"]
DEFAULT_VX_VOL = (3.58/1000. * 3.58/1000. * 40/1000.)
AXIS_FONTSIZE = 25
TICK_FONTSIZE = 13
MU = "$\mathrm{\mu}$"
AXISLBLS = {
    "volume": "Volume (" + MU + "m$^3$)",
    "surfacearea": "Surface area (" + MU + "m$^2$)",
    "pathlength": "Path length (" + MU + "m)",
    "mito_vx": "Mito volume (" + MU + "m$^3$)",
    "mito_vol": "Mito volume (" + MU + "m$^3$)",
    "mitovol": "Mito volume (" + MU + "m$^3$)",
    "nucleusvol": "Nucleus volume (" + MU + "m$^3$)",
    "cytovol": "Cytosolic volume (" + MU + "m$^3$)",
    "complexityindex": "Mitochondrial complexity index",
    "linearsynapsedensity": "Linear synapse density (" + MU + "m$^{-1}$)",
    "surfsynapsedensity": "Areal synapse density (" + MU + "m$^{-2}$)",
    "synapsecount": "Synapse count",
    "mitovoldensity": "Volume mito density",
    "linearmitocoverage": "Mito path length / path length",
    "mediandiam": "Median dend diameter ($\mathrm{\mu}$m)",
    "mitoxsec": "Effective mito cross section ($\mathrm{\mu}$m$^2$)",
    "samcs": "Mito vol / path length ($\mathrm{\mu}$m$^2$)",
    "camcs": "Covered AMCS ($\mathrm{\mu}$m$^2$)",
    "nodedist": "Distance to soma node ($\mathrm{\mu}$m)",
    "mindisttosoma": "Distance to soma node ($\mathrm{\mu}$m)",
    "meandisttosoma": "Distance to soma node ($\mathrm{\mu}$m)",
    "minsynsz": "Minimum synapse size (vx)",
    "meansynsz": "Mean synapse size (vx)",
    "mediansynsz": "Median synapse size (vx)",
    "maxsynsz": "Max synapse size (vx)",
    "synsover300": "# synapses over 300vx",
    "synsover500": "# synapses over 500vx",
    "mitodiam": "Effective mito diameter (" + MU + "m)",
    "cylindersurfarea": "Cylinder surface area (" + MU + "m$^2$)",
    "surf/cylinder": "Surface area / cylinder area ratio",
    "mitoperarea": "Mito vol / surface area (" + MU + "m)",
    "dendareaperpl": "Dend surface area / path length (" + MU + "m)",
    "cylindersynapsedensity": "Cylinder areal synapse density (" + MU + "m$^{-2}$)",
}

matplotlib.rc('font', family='sans-serif')
matplotlib.rc('font', serif='Helvetica')

#==========================================================
# Histograms

def hist_by_comp(df, varname="mito_vx", varscale=1, **kwargs):
    """Makes a histogram separated by compartment"""
    
    df = df.copy()
    df["scaled"] = df[varname] * varscale

    #df.loc[df.nodelbl == 0, "compartment"] = 'Somatic'
    #df.loc[df.nodelbl == 1, "compartment"] = 'Axonal'
    #df.loc[df.nodelbl == 2, "compartment"] = 'Basal'
    #df.loc[df.nodelbl == 3, "compartment"] = 'Apical'
    #df.loc[df.nodelbl == 4, "compartment"] = 'Unknown'

    #Ordering categories and removing "Unknown" label
    df = pd.concat((
             df[df.compartment == 'Axonal'],
             df[df.compartment == 'Somatic'],
             df[df.compartment == 'Apical'],
             df[df.compartment == 'Basal']))

    histogram(df, "compartment", "scaled", **kwargs)

    plt.xlabel(AXISLBLS[varname], fontsize=AXIS_FONTSIZE)
    plt.ylabel("Count", fontsize=AXIS_FONTSIZE)
    plt.legend(fontsize=20)

    return df


def histogram(df, groupcol, xcol, palette=None, zorder=10, n=10, **kwargs):
    """Histogram with formatting"""
    aux_args = dict(zorder=zorder, lw=3, histtype="step")
    aux_args.update(kwargs)

    bins = u.logspace_bins(df[xcol], n=n)
    for v, c in zip(reversed(pd.unique(df[groupcol])), palette):
        plt.hist(df[xcol][df[groupcol] == v], bins=bins, label=v,
                 edgecolor=c, **aux_args)

    format_axes()


#==========================================================
# Box plots


def boxplot_by_comp(df, varname="mito_vx", varscale=1, **kwargs):
    """Makes a boxplot separated by compartment"""

    df = df.copy()
    df["scaled"] = df[varname] * varscale

    #df.loc[df.nodelbl == 0, "compartment"] = 'Somatic'
    #df.loc[df.nodelbl == 1, "compartment"] = 'Axonal'
    #df.loc[df.nodelbl == 2, "compartment"] = 'Basal'
    #df.loc[df.nodelbl == 3, "compartment"] = 'Apical'
    #df.loc[df.nodelbl == 4, "compartment"] = 'Unknown'

    #Ordering categories and removing "Unknown" label
    df = pd.concat((
             df[df.compartment == 'Axonal'],
             df[df.compartment == 'Somatic'],
             df[df.compartment == 'Apical'],
             df[df.compartment == 'Basal']))

    boxplot(df, "compartment", "scaled", **kwargs)

    plt.ylabel(AXISLBLS[varname], fontsize=AXIS_FONTSIZE)
    plt.xlabel("")

    return df


def boxplot(df, xcol, ycol,
            width=0.7, fliersize=0, whis=(5, 95), zorder=10,
            **kwargs):
    """Boxplot with formatting"""
    aux_args = dict(width=width, fliersize=fliersize,
                    whis=whis, zorder=zorder, **kwargs)

    palette = None
    #palette = {compartment : "white" for compartment in np.unique(df.compartment)}
    if palette is None:
        ax = sns.boxplot(x=xcol, y=ycol, data=df, **aux_args)
    else:
        ax = sns.boxplot(x=xcol, y=ycol, **aux_args, palette=palette, data=df)

    make_black_and_white(ax)

    format_axes()


def make_black_and_white(ax):
    """Removes color from the boxes in a boxplot"""
    # iterate over boxes
    for i,box in enumerate(ax.artists):
        box.set_edgecolor('black')
        box.set_facecolor('white')
    
        # iterate over whiskers and median lines
        for j in range(6*i,6*(i+1)):
             ax.lines[j].set_color('black')


#==========================================================
# Scatter plots


def scatterplot(x, y, xlbl=None, ylbl=None,
                show_r=False, linear_fit=True, method="matplotlib",
                size=45, color='k', fontsize=AXIS_FONTSIZE, **kwargs):
    """Plots a scatter plot along with a linear fit and pearson correlation"""

    defaultkwargs = dict(alpha=0.8, edgecolor=None)
    defaultkwargs.update(**kwargs)

    # getting rid of an ugly warning message
    warnings.filterwarnings("ignore", category=FutureWarning)
    if method == "mpl":
        mplscatterdensity(np.array(x), np.array(y), **kwargs)

    elif method == "density_scatter":
        density_scatter(np.array(x), np.array(y), s=size, **kwargs)

    elif method == "seaborn":
        sns.scatterplot(x, y, color=color, s=size, **defaultkwargs)

    elif method == "matplotlib":
        if "s" in defaultkwargs:
            plt.scatter(x, y, color=color, **defaultkwargs)
        else:
            plt.scatter(x, y, color=color, s=size, **defaultkwargs)

    else:
        raise ValueError(f"unrecognized plotting method: {method}")
    
    linecolor = 'k' if method in ["density_scatter", "mpl"] else COLORS[3]
    fit = plot_linear_fit(x, y, color=linecolor) if linear_fit else None

    if len(x) > 1:
        r, p = stats.pearsonr(x, y)
    else:
        r, p = None, None

    format_axes()
    if xlbl is not None:
        plt.xlabel(xlbl, fontsize=fontsize)
    if ylbl is not None:
        plt.ylabel(ylbl, fontsize=fontsize)

    return r, p, fit


def mplscatterdensity(x, y, **kwargs):
    ax = plt.gca()
    ax.scatter_density(x, y, **kwargs)


# Density scatter plot from https://stackoverflow.com/a/53865762
def density_scatter(x, y, ax=None, sort=True, bins=50,
                    cmap=plt.get_cmap("hot"), **kwargs):
    """
    Scatter plot colored by 2d histogram
    """
    fig = plt.gcf()
    ax = plt.gca() if ax is None else ax
    data , x_e, y_e = np.histogram2d(x, y, bins=bins, density=True)

    z = interpn((0.5 * (x_e[1:] + x_e[:-1]), 0.5 * (y_e[1:] + y_e[:-1])),
                data , np.vstack([x,y]).T , method = "splinef2d",
                bounds_error = False)

    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    ax.scatter(x, y, c=z, cmap=cmap, **kwargs)

    return ax


def plot_linear_fit(X, y, lw=3, color=COLORS[3],
                    predconf=0.8, alpha=0.3):
    """Plots a linear fit to a covariate with a prediction interval

    Only supports a single regressor
    """
    fit, X = modelling.fit_linear_model(X, y)

    xlbl = [colname for colname in X.columns if colname != "const"]
    assert len(xlbl) == 1, "multiple potential x-axis variables"
    
    # Extracting prediction interval @ prediction confidence
    dummy_preds = fit.get_prediction(X)
    predint = dummy_preds.summary_frame(alpha=1-predconf)

    x = X[xlbl[0]]
    plt.plot(x, fit.fittedvalues, color=color)
    plt.fill_between(
        x, predint["obs_ci_lower"], predint["obs_ci_upper"],
        alpha=alpha, color="gray")

    return fit


def scatterplot_df(covdf, xcolname, ycolname, *args,
                   xlbl=None, ylbl=None, **kwargs):

    if xlbl is None:
        assert xcolname in AXISLBLS, "xcolname not in lookup w/ no default"
        xlbl = AXISLBLS[xcolname]

    if ylbl is None:
        assert ycolname in AXISLBLS, "ycolname not in lookup w/ no default"
        ylbl = AXISLBLS[ycolname]

    r, p, fit = scatterplot(covdf[xcolname], covdf[ycolname],
                            xlbl, ylbl, *args, **kwargs)

    return r, p, fit


def diff_scatterplot_df(df, distvarname, measurements, increment=10,
                        basegroup=65, endgroup=145, **kwargs):
    tempdf = df.copy()
    tempdf["center"] = df[distvarname] // increment * increment + increment//2

    diff_scatterplot(tempdf, measurements, groupname="center",
                     basegroup=basegroup, endgroup=endgroup, **kwargs)


def diff_scatterplot(df, measurements, groupname,
                     basegroup=65, endgroup=145, **kwargs):
    means = df.groupby(groupname)[measurements].median()

    absolutediffs = []
    relativediffs = []
    absolutediffs = means.loc[endgroup] - means.loc[basegroup]
    relativediffs = absolutediffs / means.loc[basegroup] * 100

    for i, (measurement, absdiff, reldiff) in enumerate(zip(measurements,
                                                            absolutediffs,
                                                            relativediffs)):
       plt.scatter(absdiff, reldiff, color=COLORS[i],
                   label=AXISLBLS[measurement]) 

    plt.xlabel("Absolute Difference", fontsize=AXIS_FONTSIZE)
    plt.ylabel("% Difference", fontsize=AXIS_FONTSIZE)

    format_axes()

    plt.axhline(0, color='gray', lw=1)
    plt.axvline(0, color='gray', lw=1)

    maxx = max(np.abs(absolutediffs)) * 1.1
    maxy = max(np.abs(relativediffs)) * 1.1

    plt.xlim(-maxx, maxx)
    plt.ylim(-maxy, maxy)
    #plt.legend()


def plot_cellpts_df(df, distvarname, yvarname,
                    nearlowerbnd=-1, nearupperbnd=40,
                    farlowerbnd=60, farupperbnd=155,
                    ax=None, seed=403579, axislbl=None,
                    size=10, **kwargs):
    ax = plt if ax is None else ax

    tempdf = df.copy()
    tempdf["near"] = ((nearlowerbnd <= df[distvarname])
                      & (df[distvarname] < nearupperbnd))
    tempdf["far"] = ((farlowerbnd <= df[distvarname])
                     & (df[distvarname] < farupperbnd))

    nearavgs = tempdf.groupby(["cellid", "near"])[yvarname]\
                          .mean().reset_index()
    faravgs = tempdf.groupby(["cellid", "far"])[yvarname]\
                          .mean().reset_index()

    joined = pd.merge(nearavgs[nearavgs.near][["cellid", yvarname]],
                      faravgs[faravgs.far][["cellid", yvarname]],
                      on="cellid", how="inner")

    print(f"# cells: {len(joined)}")

    if axislbl is None:
        axislbl = AXISLBLS[yvarname]

    axislbl = axislbl[0].lower() + axislbl[1:]
    scatterplot_df(joined, f"{yvarname}_x", f"{yvarname}_y",
                   xlbl=f"Proximal {axislbl}", ylbl=f"Distal {axislbl}",
                   size=size, **kwargs)

#==========================================================
# Distance gradient plots

def plot_distance_gradient_df(df, distvarname, yvarname, increment=10,
                              plottype="margin", statstype="IQR", xlbl=None,
                              collapsecol=None,
                              **kwargs):
    tempdf = df.copy()
    tempdf["center"] = df[distvarname] // increment * increment + increment//2

    if collapsecol is not None:
        tempdf = tempdf.groupby(["center", collapsecol])[yvarname]\
                     .first().reset_index()

    if xlbl is None:
        assert distvarname in AXISLBLS, "distvarname not in lookup w/ no arg"
        xlbl = AXISLBLS[distvarname]

    plot_gradient_df(tempdf, varname=yvarname, groupname="center",
                     plottype=plottype, statstype=statstype,
                     xlbl=xlbl, **kwargs)


def plot_gradient_df(df, varname, groupname, statstype,
                     plottype="margin", color=None, ax=None,
                     xlbl=None, ylbl=None, **kwargs):
    ax = plt.gca() if ax is None else ax

    stats = modelling.computestatspergroup(df, varname, groupname, statstype)

    if xlbl is None:
        assert groupname in AXISLBLS, "groupname not in lookup w/ no arg"
        xlbl = AXISLBLS[groupname]

    if ylbl is None:
        assert varname in AXISLBLS, "varname not in lookup w/ no arg"
        ylbl = AXISLBLS[varname]

    plot_gradient(stats,
                  plottype=plottype, color=color, ax=ax,
                  xlbl=xlbl, ylbl=ylbl, **kwargs)


def plot_gradient(stats, color=None, tick_spacing=50,
                  ax=None, errtype="stderr", saturation=0.75,
                  plottype="errorbar", xlbl=None, ylbl=None,
                  reqdn=10, yto0=True, printcenters=False, **kwargs):

    if saturation != 1 and color is not None:
        color = sns.desaturate(color, saturation)

    if printcenters:
        print(stats.centers)

    if plottype == "basic":
        plot_basic(stats, color=color, reqdn=reqdn, ax=ax, **kwargs)
    elif plottype == "errorbar":
        plot_errorbars(stats, color=color, reqdn=reqdn, ax=ax, **kwargs)
    elif plottype == "margin":
        plot_margin(stats, color=color, reqdn=reqdn, ax=ax, **kwargs)

    if tick_spacing is not None:
        mask = stats.ns > reqdn
        centers = stats.xs[mask]
        binwidth = centers[1] - centers[0]
        upper_bound = centers[-1] + binwidth/2.
        num_ticks = upper_bound // tick_spacing + 1  # +1 for 0
        ticks = (np.arange(num_ticks) * tick_spacing).astype(int)

        plt.xticks(ticks)

    format_axes()

    # extending y axis to 0 by default
    if yto0:
        ylim = plt.ylim()
        plt.ylim(0, ylim[1])

    if xlbl is not None:
        plt.xlabel(xlbl, fontsize=AXIS_FONTSIZE)
    if ylbl is not None:
        plt.ylabel(ylbl, fontsize=AXIS_FONTSIZE)


def plot_basic(stats, color=None, ax=None, reqdn=10, baseline=None,
               marker='o', **kwargs):
    ax = plt if ax is None else ax
    color = 'k' if color is None else color

    mask = stats.ns > reqdn
    xs = stats.xs[mask]
    centers = stats.centers[mask]

    if baseline is not None:
        baselinevalue = centers[xs == baseline][0]
        centers = (centers - baselinevalue) / baselinevalue * 100

        centers = centers[xs >= baseline]
        xs = xs[xs >= baseline]

    ax.plot(xs, centers, marker=marker, color=color, **kwargs)


def plot_margin(stats, color=None, ax=None, reqdn=10, **kwargs):
    ax = plt if ax is None else ax
    color = '#aaaaaa' if color is None else color

    mask = stats.ns > reqdn
    xs = stats.xs[mask]
    centers = stats.centers[mask]
    lo = stats.lowerbnds[mask] if stats.lowerbnds is not None else None
    hi = stats.upperbnds[mask] if stats.upperbnds is not None else None

    if lo is not None:
        assert hi is not None, "mismatched upperbnd, lowerbnd"
        ax.fill_between(xs, lo, hi, color=color, **kwargs)
    ax.plot(xs, centers, marker="o", color='k')


def plot_errorbars(stats, color=None, ax=None, reqdn=10, baseline=None, **kwargs):
    ax = plt if ax is None else ax
    color = 'r' if color is None else color

    mask = stats.ns > reqdn
    xs = stats.xs[mask]
    centers = stats.centers[mask]
    lo = (centers - stats.lowerbnds[mask]
          if stats.lowerbnds is not None else None)
    hi = (stats.upperbnds[mask] - centers
          if stats.upperbnds is not None else None)

    if baseline is not None:
        baselinevalue = centers[xs == baseline][0]
        centers = (centers - baselinevalue) / baselinevalue * 100
        lo = lo / baselinevalue * 100
        hi = hi / baselinevalue * 100

        centers = centers[xs >= baseline]
        lo = lo[xs >= baseline]
        hi = hi[xs >= baseline]
        xs = xs[xs >= baseline]

    plotkwargs = dict(lw=1, capsize=10, capthick=1, zorder=10, fmt='o')
    plotkwargs.update(kwargs)

    if lo is not None:
        assert hi is not None, "mismatched upperbnd, lowerbnd"
        ax.errorbar(xs, centers, yerr=(lo, hi), color=color, **plotkwargs)
    else:
        assert hi is None, "mismatched upperbnd, lowerbnd"
        ax.errorbar(xs, centers, color=color, **plotkwargs)


def yto0(ax=None):
    ax = plt.gca() if ax is None else ax
    ylim = ax.get_ylim()
    ax.set_ylim(0, ylim[1])

#==========================================================
# Helper Functions


def format_axes(ax=None):
    ax = plt.gca() if ax is None else ax

    ax.spines["bottom"].set_linewidth(1.5)
    ax.spines["left"].set_linewidth(1.5)
    ax.spines["right"].set_linewidth(0)
    ax.spines["top"].set_linewidth(0)

    ax.tick_params(length=5, width=1, which="minor")
    ax.tick_params(length=10, width=1.5, labelsize=18, which="major")

    updates = {"fontname": "Helvetica"}
    xticklabels = ax.get_xticklabels()
    yticklabels = ax.get_yticklabels()
    for l in xticklabels:
        l.update(updates)
    for l in yticklabels:
        l.update(updates)

    xlabel = ax.get_xlabel()
    ax.set_xlabel(xlabel, fontname="Helvetica")
    ylabel = ax.get_ylabel()
    ax.set_ylabel(ylabel, fontname="Helvetica")


def barplot_annotate_brackets(
    num1, num2, data, center, height, yerr=None,
    dh=.05, barh=.05, fs=None, fontsize=AXIS_FONTSIZE,
    logy=False, maxasterix=None, **kwargs):
    """ 
    Annotate plot with p-values.

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
    # edited from StackOverflow
    # https://stackoverflow.com/questions/11517986/indicating-the-statistically-significant-difference-in-bar-graph  # noqa

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
    barx = [lx, lx, rx, rx]    
    if not logy:
        dh *= (ax_y1 - ax_y0)
        barh *= (ax_y1 - ax_y0)

        y = max(ly, ry) + dh

        bary = [y, y+barh, y+barh, y]
        mid = ((lx+rx)/2, y+barh)
        
    else:
        y_ = max(ly, ry)
        
        h = np.log10(ax_y1) - np.log10(ax_y0)
        dh *= h
        barh *= h
        
        y = 10 ** (np.log10(y_) + dh)
        top = 10 ** (np.log10(y_) + dh + barh)
        
        bary = [y, top, top, y]
        mid = ((lx+rx)/2, top)


    plt.plot(barx, bary, c='black', **kwargs)

    kwargs = dict(ha='center', va='bottom', fontsize=fontsize)
    if fs is not None:
        kwargs['fontsize'] = fs

    plt.text(*mid, text, **kwargs)
