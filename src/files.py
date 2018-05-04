"""Specification of absolute paths to input and output files. """

from os import path

root = path.dirname(path.dirname(path.abspath(__file__)))  # root directory path
data = path.join(root, "data")  # input data directory
outputs = path.join(root, "outputs")  # output file directory

# Midthickness surface file, used to compute geodesic distances
midthickness_surface = path.join(
    data, "Q1-Q6_RelatedValidation210.L.midthickness_"
          "MSMAll_2_d41_WRN_DeDrift.32k_fs_LR.surf.gii")

# HCP MMP1.0 parcel label file for the left cortical hemisphere
parcel_labels_left = path.join(
    data, "L_Q1-Q6_RelatedValidation210.CorticalAreas_"
          "dil_Final_Final_Areas_Group_Colors.32k_fs_LR.label.gii")

# Parcellated GRIN2B gene expression profile (for demonstration purposes only)
grin2b_expression = path.join(data, "GRIN2B.pscalar.nii")

# HCP human, group-averaged (N=339), parcellated, bias-corrected myelin map
human_pmyelin_left = path.join(
    data, "L_Mean.339.MyelinMap_BC_MSMAll.32k_fs_LR.pscalar.nii")

# Geodesic distance matrix produced with Connectome Workbench commands
distance_matrix = path.join(outputs, "human_geodesic_distance_matrix.npy")

# 10k surrogate maps with spatial structure matched to the empirical myelin map
myelin_surrogates = path.join(outputs, "myelin_surrogate_maps.npy")
