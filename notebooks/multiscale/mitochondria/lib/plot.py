import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import seaborn as sns

from . import modelling


COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]
DEFAULT_VX_VOL = (3.58/1000. * 3.58/1000. * 40/1000.)
AXIS_FONTSIZE = 25
TICK_FONTSIZE = 13
MU = "$\mathrm{\mu}$"
AXISLBLS = {
    "mito_vx": "Mitochondrion volume (" + MU + "m$^3$)",
    "complexityindex": "Mitochondrial complexity index",
    "synapsedensity": "Surface synapse density (" + MU + "m$^{-2}$)",
    "mitovoldensity": "Vol mito density (" + MU + "m$^{-3}$)",
}

matplotlib.rc('font', family='sans-serif')
matplotlib.rc('font', serif='Helvetica')


#==========================================================
# Box plots


def boxplot_by_comp(df, varname="mito_vx", varscale=1):
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

    boxplot(df, "compartment", "scaled")

    plt.ylabel(AXISLBLS[varname], fontsize=AXIS_FONTSIZE)
    plt.xlabel("")

    return df


def boxplot(df, xcol, ycol,
            width=0.7, fliersize=0, whis=(5, 95), zorder=10, **kwargs):
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
                show_r=False, linear_fit=True,
                size=15, marker=True, color='k'):
    """Plots a scatter plot along with a linear fit and pearson correlation"""

    if marker is True:
        sns.scatterplot(x, y, color=color, alpha=0.8, edgecolor=None, s=size)
    else:
        #plt.errorbar(x, y, fmt='.', c='k')
        plt.scatter(x, y, color=color, s=size, marker=marker)

    fit = plot_linear_fit(x, y) if linear_fit else None

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


def plot_linear_fit(X, y, lw=3, color=COLORS[3],
                    predconf=0.8, alpha=0.2):
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
    plt.plot(x, fit.fittedvalues, color='r')
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

#==========================================================
# Helper Functions


def format_axes(ax=None):
    ax = plt.gca() if ax is None else ax

    ax.spines["bottom"].set_linewidth(1.5)
    ax.spines["left"].set_linewidth(1.5)
    ax.spines["right"].set_linewidth(0)
    ax.spines["top"].set_linewidth(0)

    ax.tick_params(length=10, width=1.5, labelsize=18)
    ax.tick_params(length=5, width=1, which="minor")


def barplot_annotate_brackets(
    num1, num2, data, center, height, yerr=None,
    dh=.05, barh=.05, fs=None, maxasterix=None):
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
    # copied from StackOverflow
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
