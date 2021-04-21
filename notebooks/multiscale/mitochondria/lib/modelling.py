import numpy as np
from scipy import stats
from scipy import optimize
import statsmodels.api as sm


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
