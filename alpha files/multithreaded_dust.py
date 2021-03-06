# %%
from random import choices
import concurrent.futures
from os import cpu_count
import numpy as np
import pandas as pd

from test import *

# %%


def generate_weighted_angles(angles):
    """Return a weighted array of possible scattering angles determined from the
    angles passed as input.

    Parameters:\n
    angles: array-like\n
    \tThe angles needing to be weighted based on their probability of being real
    scattering angles.

    Returns:\n
    p: array-like\n
    \tThe possible real scattering angles weighted by their likelihood of occurence.
    """
    # find probablilities of scattering at each angle
    SD = 0.0292
    thetaoffset = angles + SD
    exponent = np.power(thetaoffset, 2)/-0.0017059

    i = np.exp(exponent)
    j = np.sin(thetaoffset)
    k = 22.523 * np.sqrt(2*np.pi) / SD

    almost_p = np.multiply(i, j)
    p = np.multiply(almost_p, k)

    # find total probability of scattering at ANY angle
    p_total = np.sum(p)

    # given that a photon scatters, determine the probabilities of scattering
    # at each angle
    p_list = np.divide(p, p_total)

    # create list of possible scattering angles for the photons weighted
    # by their likelihood (angles can appear more than once in the list, with
    # more likely angles appearing more often)
    scatter_angles = np.round(choices(angles, weights=p_list, k=10**5), 3)
    return scatter_angles


def dust_chunk(chunkphotons: int):
    """Return a Pandas dataframe of the specified number of photons with random
    values of Right Ascension, Declination, and line-of-sight distance associated
    with them.

    Parameters:\n
    nphotons: int\n
    \tThe number of photons (and therefore dust grains) to be generated.

    Returns:\n
    data: DataFrame\n
    \tA 2d array containing location information for each generated photon.
    """
    # create dataframe
    chunkdata = pd.DataFrame()

    # generate specifed photons with random RA, Dec, and distance
    # (d is the distance to the dust grain along the line of sight
    # to the source)
    chunkdata["ra"] = np.random.randint(0, 36000, size=chunkphotons)/100
    chunkdata["dec"] = np.random.randint(0, 9000, size=chunkphotons)/100
    chunkdata["d"] = np.random.randint(0, 100, size=chunkphotons)

    # ra_diff and dec_diff are the differences in position of the
    # dust grains with respect to the source and along the line of
    # sight
    ra_diff, dec_diff = -chunkdata["ra"]+source[0], -chunkdata["dec"]+source[1]
    # r is the degree separations between the source and the dust
    # grains
    chunkdata["r"] = np.sqrt(np.power(ra_diff, 2) + np.power(dec_diff, 2))

    # angle is the angles the photons need to scatter to reach us, given
    # the locations of the dust grains
    chunkdata["angle"] = np.round(np.arctan(chunkdata["r"]/chunkdata["d"]), 3)

    # scatter angle is the angles at which each dust grain would actually
    # scatter a photon
    chunkdata["scatter angle"] = choices(weighted_scatter, k=chunkphotons)
    # print(chunkdata.info())
    # # isolate all photons that scattered to us
    # scattered = chunkdata[chunkdata["scatter angle"] == chunkdata["angle"]]

    # # calculate true distances to the grains that scattered them
    # scattered["d_grain"] = scattered["d"]/np.cos(scattered["scatter angle"])

    # # calculate the light travel time from grains to us
    # speedoflight = 3
    # scattered["t"] = scattered["d_grain"]/speedoflight

    # return populated DataFrame
    # return_list.append(chunkdata)
    # return return_list
    return chunkdata


def batch_dust(nphotons):
    """
    """
    # determine number and size of chunks to be created
    NCHUNKS = cpu_count()
    n = nphotons//NCHUNKS

    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = [executor.submit(dust_chunk, n) for _ in range(NCHUNKS)]
        for f in concurrent.futures.as_completed(results):
            print(f.result())
        return results 


if __name__ == "__main__":
    # create list of angles
    thetas = np.arange(0, np.pi, step=10**-5)

    # generate weighted list of possible scattering angles
    weighted_scatter = generate_weighted_angles(thetas)

    # set source location and number of photons to be generated
    source = (180.0, 45.0)
    NPHOTONS = 10**8

    # generate dust layer/photons
    a = batch_dust(NPHOTONS)
    # print(a)
    # for i in a:
    #     print(i)


# %%
