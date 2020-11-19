import os
from collections import Counter
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


def find_centroids(filename):
    """This (crude) algorithm iterates through the pixels of an image and returns the
    RA and Dec coordinates of the maxima of the image as two 1D numpy arrays.
    """
    with fits.open(filename) as file:
        # load in RA/Dec coordinates
        xs = np.array(file[1].data["X"])
        ys = np.array(file[1].data["Y"])

        # convert to nearest integers
        xs = np.rint(xs)
        ys = np.rint(ys)

        # determine x and y bins
        event_coordinates = np.stack((xs, ys), axis=-1)
        n_unique_events = len(np.unique(event_coordinates))

        print("n_event_coordinates = " + str(len(event_coordinates)))
        print("n_unique_elements = " + str(n_unique_events))

        XMIN, XMAX = np.min(xs), np.max(xs)
        YMIN, YMAX = np.min(ys), np.max(ys)
        SIZE_COORDINATE_MESH = (XMAX-XMIN) * (YMAX-YMIN)

        # xbinedges = np.histogram_bin_edges(xs, bins='fd')
        # ybinedges = np.histogram_bin_edges(ys, bins='fd')

        # iterate through the bins on each axis
        image_maxes = []
        for j in range(len(ybinedges)-1):
            for i in range(len(xbinedges)-1):
                cell_xmin, cell_xmax = xbinedges[i], xbinedges[i+1]
                cell_ymin, cell_ymax = ybinedges[j], ybinedges[j+1]

                # isolate the events within each bin
                cell_constraints = np.where((xs < cell_xmax) & (
                    xs >= cell_xmin) & (ys < cell_ymax) & (ys >= cell_ymin))

                cell_xs = xs[cell_constraints]
                cell_ys = ys[cell_constraints]
                coord_pairs = list(zip(cell_xs, cell_ys))

                # find the maximum point in the bin
                coord_counts = Counter(coord_pairs).most_common(1)
                for result in coord_counts:
                    image_maxes.append(result)

        # sort the identified maxima by count and return their x(RA) and y(Dec) coordinates
        image_maxes = sorted(image_maxes, key=lambda x: x[1], reverse=True)

        centroid_xs = np.array([item[0][0] for item in image_maxes])
        centroid_ys = np.array([item[0][1] for item in image_maxes])
        return centroid_xs, centroid_ys


print(find_centroids("/Users/hunterholland/Documents/Research/Laidlaw/Data/Modified/L1517/Swift/xrt/event/sw00034249004xpcw3po_cl.evt.gz"))
