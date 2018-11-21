""" Compute statistical significance with spatial autoregressive modeling. """

from scipy.optimize import curve_fit
import numpy as np
import pysal as ps


def fit_exp_decay(x, y):
    """ Fit an exponentially decaying function of the form y = exp(-x/x0) for
    dependent variable `y` evaluated at positions `x`.

    Parameters
    ----------
    x : np.ndarray
        independent variable
    y : np.ndarray
        dependent variable

    Returns
    -------
    x0 : float
        OLS estimate of parameter x0

    """
    def _func(z, z0):
        return np.exp(-z/z0)
    popt, pcov = curve_fit(_func, x, y)
    return popt[0]


def fit_empirical_autocorr_scale(corrmat, distmat, dmax=None):
    """ Fit spatial autocorrelation length scale to empirical data.

    Parameters
    ----------
    corrmat : np.ndarray
        correlation matrix whose element [i, j] is the correlation (of a measure
        such as gene expression profiles) between parcels i and j
    distmat : np.ndarray
        matrix whose i,j-th element is the geodesic distance btwn parcels i & j
    dmax : float, optional
        cutoff distance; if not None, only elements of `corrmat` and `distmat`
        for which the element in `distmat` is less than `dmax` are used when
        estimating the empirical spatial autocorrelation scale, d0

    Returns
    -------
    d0 : float
        characteristic scale of spatial autocorrelation

    Notes
    -----
    The characteristic scale of spatial autocorrelation is estimated by fitting
    the parameter d0 in the equation ``r = exp(-d/d0)`` using OLS, where r is
    an element of the correlation matrix `corr`, and d is the corresponding
    element of the geodesic distance matrix. The upper triangular parts of the
    distance and correlation matrices are flattened before estimating d0.

    """
    nr, nc = corrmat.shape
    assert nr == nc and (nr, nc) == distmat.shape
    triu_inds = np.triu_indices(nr, k=1)
    x = distmat[triu_inds]
    y = corrmat[triu_inds]
    if dmax is not None:
        n = x.size
        mask = np.less(x, dmax)
        x = x[mask]
        y = y[mask]
        print("# Fitting d0 on %i of %i elements for which d < (dmax = %f)" % (
            sum(mask), n, dmax))
    d0 = fit_exp_decay(x, y)
    return d0


def pysal_weight_matrix(d, d0, normalize=True):
    """ Construct a weight matrix as a PySAL W weights object.

    Parameters
    ----------
    d : np.ndarray
        matrix whose i,j-th element is the geodesic distance btwn parcels i & j
    d0 : float
        characteristic scale of spatial autocorrelation (e.g., returned from
        fit_empirical_autocorr_scale)
    normalize : bool, optional
        if True, normalize each row of the weight matrix

    Returns
    -------
    ps.W
        PySAL W weights object

    Notes
    -----
    This function returns a PySAL W weights object, for use with the spatial
    autoregressive modeling routine pysal.spreg.ml_lag.ML_Lag().

    """

    nr, nc = d.shape
    assert nr == nc

    # Create a weight matrix using the input spatial autocorrelation scale
    weights = np.exp(-d / d0)

    # Mask diagonal elements to remove self-coupling; optionally normalize
    diagmask = np.eye(nr, dtype=bool)
    weights[diagmask] = 0
    if normalize:
        weights /= weights.sum(axis=1)

    # Construct pysal weight object directly from the np.ndarray
    return ps.weights.full2W(weights)


def slm_ml_beta_pval(x, y, w):
    """ Compute p-value associated with the maximum likelihood (ML) estimate of
    spatial lag model (SLM) parameter \beta.

    Parameters
    ----------
    x : np.ndarray
        independent variable
    y : np.ndarray
        dependent variable
    w : ps.W
        PySAL W weights object

    Returns
    -------
    beta : float
        direct (i.e. local) impact on the dependent variable `y` due to a unit
        change in independent variable `x` assuming spatial structure defined
        by the PySAL weights object `w`

    Notes
    -----
    The SLM has the functional form ``y = \rho W y + x \beta + \nu``, where \rho
    scales the strength of spatial autocorrelation; W is a user-defined weight
    matrix that implicitly specifies the form of spatial structure in the data;
    and \nu is normally distributed. This function fits an SLM using the input
    x, y and w, and returns the statistical significance for parameter \beta.

    """

    assert x.size == y.size
    assert x.ndim == y.ndim == 1

    # Transform 1d numpy arrays to column vectors
    x = np.expand_dims(x, 1)
    y = np.expand_dims(y, 1)

    # Compute maximum likelihood estimation of spatial lag model
    res = ps.spreg.ml_lag.ML_Lag(y, x, w)

    # Return parameter beta, which reflects the direct (i.e. local) impact on
    # dependent variable y due to a unit change in independent variable x
    beta = res.z_stat[1][1]

    return beta
