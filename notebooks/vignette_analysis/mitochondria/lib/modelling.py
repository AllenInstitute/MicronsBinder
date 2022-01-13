import numpy as np
from collections import namedtuple
from scipy import stats
from scipy import optimize
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sklearn.decomposition import FastICA


GroupStats = namedtuple("GroupStats",
                        ["xs", "centers", "ns", "lowerbnds", "upperbnds"])


def computestatspergroup(df, varname, groupname, statstype):
    xs = np.unique(df[groupname])

    cs, ns, lbs, ubs = [], [], [], []
    for x in xs:
        subdf = df[df[groupname] == x]
        n = len(subdf)
        c, lb, ub = computestat(subdf, varname, statstype)

        cs.append(c)
        ns.append(n)

        if lb is not None:
            lbs.append(lb)
        if ub is not None:
            ubs.append(ub)

    cs = np.array(cs)
    ns = np.array(ns)
    lbs = np.array(lbs)
    ubs = np.array(ubs)

    assert len(lbs) == 0 or len(lbs) == len(cs)
    assert len(ubs) == 0 or len(ubs) == len(cs)
    assert len(lbs) == len(ubs)

    if len(lbs) == 0:
        return GroupStats(xs, cs, ns, None, None)
    else:
        return GroupStats(xs, cs, ns, lbs, ubs)


def computestat(subdf, varname, statstype):
    if statstype == "meanstddev":
        return meanstddev(subdf, varname)
    if statstype == "mean":
        return meanstddev(subdf, varname)
    if statstype == "median":
        return meanstddev(subdf, varname)
    elif statstype == "meanstderr":
        return meanstderr(subdf, varname)
    elif statstype == "IQR":
        return iqr(subdf, varname)
    elif statstype == "meanIQR":
        center = np.mean(subdf[varname])
        _, lb, ub = iqr(subdf, varname)
        return center, lb, ub
    elif statstype == "COV":
        return cov(subdf, varname)
    if statstype == "mean":
        return np.mean(subdf[varname]), None, None
    if statstype == "median":
        return np.percentile(subdf[varname], 50), None, None
    else:
        raise ValueError(f"unknown statstype: {statstype}")


def meanstddev(subdf, varname):
    center = np.mean(subdf[varname])
    stddev = np.std(subdf[varname])
    lowerbnd = center - stddev
    upperbnd = center + stddev

    return center, lowerbnd, upperbnd


def meanstderr(subdf, varname):
    center = np.mean(subdf[varname])
    with np.errstate(all="ignore"):
        stderr = stats.sem(subdf[varname])
    lowerbnd = center - stderr
    upperbnd = center + stderr

    return center, lowerbnd, upperbnd


def iqr(subdf, varname):
    return np.percentile(subdf[varname], [50, 25, 75])


def cov(subdf, varname):
    mean = np.mean(subdf[varname])
    stddev = np.std(subdf[varname])

    return stddev/mean, None, None


def regression(df, colnames, ycolname, groupcolname=None, return_model=False):
    if groupcolname is None:
        model = smf.ols(f"{ycolname} ~ {' + '.join(colnames)}", df)
    else:
        model = smf.mixedlm(f"{ycolname} ~ {' + '.join(colnames)}",
                            df, groups=df[groupcolname])

    if return_model:
        return model.fit(), model
    else:
        return model.fit()


def residualregression(df, colnames1, colnames2, ycolname, groupcolname=None):
    initlm = regression(df, colnames1, ycolname, groupcolname=groupcolname)

    tempdf = df.copy()
    tempdf["resid"] = initlm.resid

    return regression(tempdf, colnames2, "resid"), tempdf


def ICAregression(df, colnames, ycolname, groupcolname=None, seed=None):
    df = df.copy()
    Xraw = df[colnames].values

    # standardizing
    X = (Xraw - Xraw.mean(0)) / Xraw.std(0)

    if seed is not None:
        ica = FastICA(n_components=len(colnames), random_state=seed).fit(X)
    else:
        ica = FastICA(n_components=len(colnames)).fit(X)
    T = ica.transform(X)

    Cnames = list()
    for i in range(len(colnames)):
        Cnames.append(f"C{i}")
        df[Cnames[-1]] = T[:, i]

    if groupcolname is None:
        lm = smf.ols(f"{ycolname} ~ {' + '.join(Cnames)}", df).fit()
    else:
        lm = smf.mixedlm(f"{ycolname} ~ {' + '.join(Cnames)}",
                         df, groups=df[groupcolname]).fit()

    return lm, ica, T


SVD = namedtuple("SVD", ["U", "s", "Vt"])


def PCAregression(df, colnames, ycolname, groupcolname=None):
    df = df.copy()
    Xraw = df[colnames].values

    # standardizing
    X = (Xraw - Xraw.mean(0)) / Xraw.std(0)

    U, s, Vt = np.linalg.svd(X.T, full_matrices=False)

    Cnames = list()
    for (i, v) in enumerate(Vt):
        Cnames.append(f"C{i}")
        df[Cnames[-1]] = v.ravel()

    if groupcolname is None:
        lm = smf.ols(f"{ycolname} ~ {' + '.join(Cnames)}", df).fit()
    else:
        lm = smf.mixedlm(f"{ycolname} ~ {' + '.join(Cnames)}",
                         df, groups=df[groupcolname]).fit()

    svd = SVD(U, s, Vt)

    return lm, svd


def fit_linear_model(X, y, sort=True):
    if sort:
        # assume a single regressor for now
        order = np.argsort(X)
        X = X.iloc[order]
        y = y.iloc[order]

    X = sm.add_constant(X)
    return sm.OLS(y, X).fit(), X


def fit_models_to_bins(plotdf, modeltype, varname="mitovols"):
    centers = np.unique(plotdf.bin_center)

    models = list()
    for center in centers:
        subset = plotdf[plotdf.bin_center == center]
        models.append(fit_model(subset[varname], modeltype))

    return models


def fit_models_to_compartments(df, modeltype, varname="mitovols"):
    models = dict()
    for (c, subset) in df.groupby("compartment"):
        models[c] = fit_model(subset[varname], modeltype)

    return models


def fit_model(d, modeltype):
    if modeltype == "norm":
        return Norm(d)
    elif modeltype == "lognorm":
        return LogNorm(d)
    elif modeltype == "skewnorm":
        return SkewNorm(d)
    elif modeltype == "logskewnorm":
        return LogSkewNorm(d)
    else:
        raise Exception(f"unknown model type: {modeltype}")


class Model(object):

    def __init__(self, data, modeltype):
        self.modeltype = modeltype
        self.params = self.fit(data)

    def fit(self, data):
        return self.modeltype.fit(data)

    def pdf(self, samples):
        model = self.modeltype(*self.params)
        return model.pdf(samples)

    def stats(self, moments="mv"):
        model = self.modeltype(*self.params)
        return model.stats(moments)


class Norm(Model):

    def __init__(self, data):
        super(Norm, self).__init__(data, stats.norm)


class LogNorm(Model):

    def __init__(self, data):
        super(LogNorm, self).__init__(data, stats.lognorm)


class SkewNorm(Model):

    def __init__(self, data):
        super(SkewNorm, self).__init__(data, stats.skewnorm)


class LogSkewNorm(Model):

    def __init__(self, data):
        super(LogSkewNorm, self).__init__(
            np.log10(data), stats.skewnorm)

    def pdf(self, samples):
        return self.modeltype(*self.params).pdf(np.log10(samples))


class PiecewiseLinear(object):

    def __init__(self, xs, ys):
        self.params, self.pcov = self.fit(xs, ys)

    @staticmethod
    def pred(x, bx=None, by=None, w1=None, w2=None):
        bx = self.params[0] if bx is None else bx
        by = self.params[1] if by is None else by
        w1 = self.params[2] if w1 is None else w1
        w2 = self.params[3] if w2 is None else w2

        return np.piecewise(x, [x < bx],
                            [lambda x: by + w1*(x-bx),
                             lambda x: by + w2*(x-bx)])

    def fit(self, xs, ys):
        p0 = (np.mean(xs), np.mean(ys), 1, -1)
        return optimize.curve_fit(self.pred, xs, ys, p0=p0, absolute_sigma=True)

    def perr(self):
        return np.sqrt(np.diag(self.pcov))

    def summary(self):
        bx, by, w1, w2 = self.params
        bxerr, byerr, w1err, w2err = self.perr()
        print(f"Piecewise Linear Fit: \n"
              f" x-breakpoint: {bx:3f} +/- {bxerr:3f} \n"
              f" y-breakpoint: {by:3f} +/- {byerr:3f} \n"
              f" slope before bp: {w1:3f} +/- {w1err:3f} \n"
              f" slope after bp: {w2:3f} +/- {w2err:3f}")
