""" Utility functions for data I/O (including geodesic distance matrix). """

from os import path, remove, system
import nibabel as nib
import numpy as np
from src import files


def load_gifti(f):
    """ Load data stored in a GIFTI (.gii) neuroimaging file.

    Parameters
    ----------
    f : str
        path to gifti (.gii) file

    Returns
    -------
    np.ndarray

    """
    return nib.load(f).darrays[0].data


def load_nifti(f):
    """ Load data stored in a NIFTI (.nii) neuroimaging file.

    Parameters
    ----------
    f : str
        path to nifti (.nii) file

    Returns
    -------
    np.ndarray

    """
    return np.array(nib.load(f).get_data()).squeeze()


def geodesic_distance(regenerate=False):
    """Compute/load average parcel-wise geodesic distance matrix.

    Parameters
    ----------
    regenerate : bool, optional
        if True, recompute geodesic distance matrix using Connectome Workbench

    Returns
    -------
    distance_matrix : np.ndarray
        element [i, j] is the geodesic distance between parcels i and j

    """

    if not regenerate:
        if path.exists(files.distance_matrix):
            return np.load(files.distance_matrix)
        else:
            raise Exception("Distance file %s doesn't exist! Re-run with "
                            "regenerate=True" % files.distance_matrix)

    print "\n## Computing geodesic distance matrix ##"
    print "# Input surface file: %s" % files.midthickness_surface
    print "# Output distance file: %s" % files.distance_matrix

    # Files produced during run-time execution bv Connectome Workbench commands
    coordinate_metric_file = path.join(
        files.outputs, "vertex_coordinates.func.gii")
    distance_metric_file = path.join(
        files.outputs, "geodesic_distance.func.gii")

    # Create a metric file containing the coordinates of each surface vertex
    system(
        'wb_command -surface-coordinates-to-metric "%s" "%s"' % (
            files.midthickness_surface, coordinate_metric_file))

    # Load left-hemispheric surface vertex parcel labels
    labels = load_gifti(files.parcel_labels_left)

    # Ensure the number of parcel labels equals the number of surface vertices
    vertices = load_gifti(files.midthickness_surface)
    assert labels.size == vertices.shape[0]

    # Skip parcel label 0 if present -- not a parcel
    unique_labels = np.unique(labels)
    nparcels = unique_labels.size
    if 0 in unique_labels:
        unique_labels = unique_labels[unique_labels != 0]
        assert unique_labels.size == (nparcels - 1)
        nparcels -= 1

    # Create vertex-level mask for each unique cortical parcel
    parcel_vertex_mask = {l: labels == l for l in unique_labels}

    # Loop over pairwise parcels at the level of surface vertices
    distance_matrix = np.zeros((nparcels, nparcels))
    for i, li in enumerate(unique_labels[:-1]):

        # Labels of parcels for which to compute mean geodesic distance
        labels_lj = unique_labels[i+1:]

        # Initialize lists in which to store pairwise vertex-level distances
        parcel_distances = {lj: [] for lj in labels_lj}

        # Loop over vertices with parcel label i
        li_vertices = np.where(parcel_vertex_mask[li])[0]
        for vi in li_vertices:

            # Compute and load distance from vertex vi to every other vertex
            system(
                'wb_command -surface-geodesic-distance "%s" %i "%s" ' % (
                    files.midthickness_surface, vi, distance_metric_file))
            distance_from_vi = nib.load(distance_metric_file).darrays[0].data

            # Update lists w/ distances from vertex vi to vertices in parcel j
            for lj in labels_lj:
                vi_lj_distances = distance_from_vi[parcel_vertex_mask[lj]]
                parcel_distances[lj].append(vi_lj_distances)

        # Compute average geodesic distances
        for j, lj in enumerate(labels_lj):
            mean_distance = np.mean(parcel_distances[lj])
            distance_matrix[i, i+j+1] = mean_distance

        print "# Parcel label %s complete." % str(li)

    # Remove intermediate files produced during run-time execution
    remove(coordinate_metric_file)
    remove(distance_metric_file)

    # Make resultant matrix symmetric
    i, j = np.triu_indices(nparcels, k=1)
    distance_matrix[j, i] = distance_matrix[i, j]

    np.save(files.distance_matrix, distance_matrix)
    return distance_matrix
