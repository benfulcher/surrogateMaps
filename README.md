# New

Generating correlated spatial maps for mouse and human.
We can run:
```bash
python3 GenerateMapsFixed.py
```
after manually setting either `human` or `mouse` in the script.

---

# From [the original Murray Lab Bitbucket repository](https://bitbucket.org/murraylab/surrogates/src/master/):

The functions defined in `src` operate on the data in `data` to produce the
outputs in `outputs`. All core functionality is demonstrated in `demo.py`.

Descriptions of package contents are provided below.

### data
Input data, treated as read-only.

GRIN2B.pscalar.nii
    The parcellated expression map for gene GRIN2B (for demonstration purposes).

L_Mean.339.MyelinMap_BC_MSMAll.32k_fs_LR.pscalar.nii
    HCP group-averaged myelin map values for parcels in the left cortical
    hemisphere.

L_Q1-Q6_RelatedValidation210.[...]Group_Colors.32k_fs_LR.label.gii
    HCP MMP1.0 cortical parcel labels for surface vertices in the left-hemisphere
    surface mesh.

Q1-Q6_RelatedValidation210.L.midthickness_[...]32k_fs_LR.surf.gii
    The HCP midthickness surface mesh for the left cortical hemisphere.

### docs
The latest version of manuscript, which contains detailed descriptions of the
methods implemented in this package.

### src
Source code, i.e. definitions of functions which operate on files in `data`.

`files.py`
    Absolute paths to files in the `data` and 'outputs' directories.

`spatial.py`
    Functions to compute statistical significance using spatial autoregressive
    modeling.

`surrogates.py`
    Functions to generate and save surrogate maps with a specific
    autocorrelation structure.

`utils.py`
    Utility functions for data I/O and to generate and save the geodesic
    distance matrix.

### outputs
Directory in which outputs produced by the functions defined in `src` are saved.
