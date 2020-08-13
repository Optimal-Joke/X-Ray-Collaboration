import os
from ctypes import CDLL, POINTER, c_int32
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np


def halo(filepath, saveimg=False):
    """This method iterates through the input event file and uses an implementation of the Fuzzy Hough Transform
    to determine a series of the most optimal halos of incident photons. It then superimposes the most optimal
    halo onto the image and returns a PNG file which gets saved to the event file's directory.
    """

    # get event data from FITS file
    hdul = fits.open(filepath)
    xs = np.array(hdul[1].data['X'])
    ys = np.array(hdul[1].data['Y'])
    photonCount = len(xs)

    # cpath should be the absolute path to the C file which does the Fuzzy Hough calculation (halo_calculation.c)
    cpath = "/Users/hunterholland/Documents/Research/Laidlaw/Collab. Repository/X-Ray-Collaboration/halo_calculation.c"
    cfilename = os.path.basename(cpath)
    sodir = os.path.dirname(filepath)
    sopath = f"{sodir}/{cfilename}"
    # create C library files in the file's directory
    os.system(
        f'gcc -c -fPIC -lm {cpath}.c -o {sopath}.o; gcc {sopath}.o -shared -o {sopath}.so')
    houghlib = CDLL(f'{sopath}.so')

    c_int_p = POINTER(c_int32)

    # convert numpy arrays from int16 to a C-compatable int32
    xs = xs.astype(np.int32)
    ys = ys.astype(np.int32)

    # cast arrays as ctypes
    cxs = xs.ctypes.data_as(c_int_p)
    cys = ys.ctypes.data_as(c_int_p)

    # set the types for the inputs and outputs to the C hough function
    houghlib.hough.argtypes = [c_int_p, c_int_p, c_int32]
    houghlib.hough.restype = c_int_p
    # call the hough function and store the returned, "most likey" parameters in halo
    haloinfo = houghlib.hough(cxs, cys, photonCount)
    # stores the (x, y) coordinates  of the center
    center = (haloinfo[0], haloinfo[1])
    r = haloinfo[2]  # stores the radius
    houghlib.free(haloinfo)
    print('\nHighest accumulator halo: center = ' +  # note that this only outputs the single highest accumulator-value halo
          str(center) + ', r = ' + str(r))

    # create image with the most likely halo superimposed
    _fig, ax = plt.subplots()
    ax.set_facecolor('black')
    plt.scatter(xs, ys, s=0.1, color='white')
    ax.add_artist(plt.Circle(center, r, color='r', fill=False))

    # save image
    if saveimg == True:  # save PNG
        plt.savefig(
            f"{sodir}/halo.png", dpi=250, format="png")
        print(f"Histogram saved to {sodir}.")

    # show image
    else:
        plt.show()

