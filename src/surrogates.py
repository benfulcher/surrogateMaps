""" Generate surrogate maps with a specific autocorrelation structure. """

from os import path
from scipy.optimize import leastsq
from scipy.stats import boxcox
import numpy as np


def construct_weight_matrix(d, d0):
    """ Construct weight matrix with functional form W_{ij} = exp(-d_{ij} / d0).

    Parameters
    ----------
    d : np.ndarray
        distance matrix where element [i, j] is distance between parcels i and j
    d0 : float
        characteristic scale of spatial autocorrelation

    Returns
    -------
    w : np.ndarray
        weight matrix

    """
    nr, nc = d.shape
    assert nr == nc
    w = np.exp(-d / d0)
    identity = np.eye(nr, dtype=bool)
    w[identity] = 0  # Zero diagonal elements to remove self-coupling
    # w /= w.sum(axis=1)  # Normalize rows to account for variation in parcel size
    return w


def fit_parameters(d, neuro_map):
    """ Fit spatial autocorrelation parameters for a neuroimaging map.

    Parameters
    ----------
    d : np.ndarray
        distance matrix where element [i, j] is distance between parcels i and j
    neuro_map : np.ndarray
        parcellated neuroimaging/neurophenotype scalar map, e.g. the myelin map

    Returns
    -------
    np.ndarray
        weight matrix

    """
    if not np.all(neuro_map > 0):
        raise Exception(
            "Input map must be positive definite for Box-Cox transformation")

    # Transform `neuro_map` scalars to be approximately Gaussian using Box-Cox
    y, _ = boxcox(neuro_map)

    # Standardize (i.e. z-score) values
    y -= y.mean()
    y /= y.std()

    # Use OLS to estimate \rho &  d0 using spatial lag model (SLM) of the form
    # y = \rho * W * y, where \rho is a scalar, W is the weight matrix whose
    # elements W_{ij} have the functional form W_{ij} = exp(-d_{ij} / d0), and
    # d_{ij} is the distance between parcels i and j
    def _slm(params):
        rho_, d0_ = params
        w = construct_weight_matrix(d, d0_)
        return y - rho_ * w.dot(y)
    x0 = np.array([0.5, 0.2])  # Initial parameters estimates (i.e., rho, d0)
    x, cov_x, infodict, mesg, ier = leastsq(_slm, x0=x0, full_output=True)
    rho, d0 = x

    return rho, d0


def generate(d, rho, d0):
    """ Generate a surrogate map with spatial autocorrelation structure.

    Parameters
    ----------
    d : np.ndarray
        distance matrix where element [i, j] is distance between parcels i and j
    rho : float
        strength of spatial autocorrelation
    d0 : float
        characteristic spatial autocorrelation length scale

    Returns
    -------
    np.ndarray
        surrogate map with spatial structure determined by `d`, `rho`, and `d0`

    """

    nr, nc = d.shape
    assert nr == nc

    # Construct weight matrix
    w = construct_weight_matrix(d, d0)

    # Generate random vector of normally distributed samples
    u = np.random.randn(nr)

    # Introduce spatial structure into random variates
    identity = np.eye(nr, dtype=bool)
    x = np.linalg.inv(identity - rho * w).dot(u)

    return x


def save_surrogates(d, neuro_map, saveto, n_maps=10000):
    """ Generate and save surrogate maps with spatial autocorrelation structure.

    Parameters
    ----------
    d : np.ndarray
        distance matrix where element [i, j] is distance between parcels i and j
    neuro_map : np.ndarray
        parcellated neuroimaging/neurophenotype scalar map, e.g., the myelin map
    saveto : str
        file to which surrogate maps are saved (with a function call to np.save)
    n_maps : int, optional
        number of surrogate maps to generate

    Returns
    -------
    np.ndarray
        random surrogate maps with shape (n_maps, n_parcels)

    """

    nr, nc = d.shape
    assert nr == nc == neuro_map.size

    # Set random seed generator
    np.random.seed(137)

    if path.exists(saveto):
        print("# Warning! Overwriting surrogates saved to %s." % saveto)

    # Fit spatial autocorrelation parameters for input map
    rho, d0 = fit_parameters(d, neuro_map)
    print("# Parameter estimates: rho = %f; d0 = %f" % (rho, d0))

    print("# Generating %u surrogates... " % n_maps)
    maps = np.empty((n_maps, nr))

    # Generate random surrogates. In the code block below, surrogate maps are
    # generated with spatial autocorrelation structure matched to the empirical
    # map. To match surrogate map value distributions to the distribution of
    # values in the empirical map, rank-ordered random surrogate map values are
    # re-assigned the corresponding rank-ordered values in the empirical map.
    # Note that this approach approximates a spatial autocorrelation-preserving
    # permutation test of the empirical neuroimaging map.
    sorted_neuro_map = np.sort(neuro_map)
    for i in range(n_maps):
        surr = generate(d, rho, d0)
        ii = np.argsort(surr)
        np.put(surr, ii, sorted_neuro_map)
        maps[i] = surr

    if saveto:
        np.save(saveto, maps)
        print("# Surrogate maps saved to %s." % saveto)

    return maps


def load_surrogates(f):
    """Load surrogate maps saved as a .npy file with np.save(). """
    return np.load(f)
