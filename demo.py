"""

Demonstration script which illustrates how to:
    i) compute parcels' pairwise geodesic distances,
    ii) generate 10,000 random surrogate myelin maps,
    iii) use the spatial lag model to compute spatial-autocorrelation-corrected
      statistical significance values, and
    iv) use surrogate maps to construct null distributions for a test statistic.

"""

from src import utils, files, surrogates, spatial
from scipy.stats import pearsonr
import numpy as np

recompute_distMat = False  # recompute geodesic distance matrix (slow!)

# Load mean pairwise geodesic distance matrix. If `recompute_distMat` is True,
# it will be recomputed using Connectome Workbench commands and saved to
# `outputs`. Note that recomputing the matrix takes ~ 1-2 hours to run on my
# machine. Once it's been computed, subsequent calls to utils.geodesic_distance
# with argument `regenerate=False` will load the saved result from memory.
distMat = utils.geodesic_distance(regenerate=recompute_distMat)

# Load the parcellated human myelin map (for the left hemisphere)
parcel_myelin = utils.load_nifti(files.human_pmyelin_left)

# Generate and save 10,000 random surrogate maps with spatial autocorrelation
# structure matched to the empirical myelin map.
myelin_surrogate_maps = surrogates.save_surrogates(
    distMat, parcel_myelin, n_maps=10000, saveto=files.myelin_surrogates)

# Construct a PySAL W weights object using the spatial scale, 25.3mm, estimated
# for empirical gene expression with the command (not included in this demo)
# `spatial.fit_empirical_autocorr_scale(corrmat, distMat, dmax=100)`, where
# `corrmat` was the pairwise gene expression correlation matrix across parcels.
weights = spatial.pysal_weight_matrix(distMat, d0=25.3, normalize=True)

# Compute p-value for the \beta parameter in the spatial lag model fit to the
# empirical myelin map and expression profile for gene GRIN2B, and compare the
# result to the standard Pearson correlation coefficient
grin2b_map = utils.load_nifti(files.grin2b_expression)
grin2b_map = grin2b_map - grin2b_map.min() + 0.1
surrogates.fit_parameters(distMat, grin2b_map)
p_slm = spatial.slm_ml_beta_pval(x=parcel_myelin, y=grin2b_map, w=weights)
r_pearson, p_pearson = pearsonr(parcel_myelin, grin2b_map)
print("# P-value, SLM: %.3e" % p_slm)
print("# P-value, Pearsonr: %.3e" % p_pearson)

# Alternatively, we coould compute a significance value by constructing a null
# distribution using the surrogate maps:
null_dist = np.array([pearsonr(m, grin2b_map)[0] for m in myelin_surrogate_maps])
n_more_extreme = np.greater(np.abs(null_dist), np.abs(r_pearson)).sum()
p_surrogates = n_more_extreme / float(null_dist.size)
print("# P-value, surrogates: %.3e" % p_surrogates)

# This was a trivial example, but the surrogate map method above can be easily
# extended to construct the null distribution for any desired test statistic.
