"""This program simulates the scattering of a specified number of
X-ray photons through a cloud in the interstellar medium.

authors: Susannah Abrams, Hunter Holland
"""

# %%
# ___________________________IMPORTING NECESSARY MODULES___________________________
from random import choice, choices
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# %%
# _______________________________DEFINING FUNCTIONS________________________________


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


def generate_dust(nphotons):
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
    data = pd.DataFrame(index=np.arange(0, nphotons))

    # generate specifed photons with random RA, Dec, and distance
    # (d is the distance to the dust grain along the line of sight
    # to the source)
    data["ra"] = np.random.randint(0, 36000, size=nphotons)/100
    data["dec"] = np.random.randint(0, 9000, size=nphotons)/100
    data["d"] = np.random.randint(0, 100, size=nphotons)

    # ra_diff and dec_diff are the differences in position of the
    # dust grains with respect to the source and along the line of
    # sight
    ra_diff, dec_diff = -data["ra"]+source[0], -data["dec"]+source[1]
    # r is the degree separations between the source and the dust
    # grains
    data["r"] = np.sqrt(np.power(ra_diff, 2) + np.power(dec_diff, 2))

    # angle is the angles the photons need to scatter to reach us, given
    # the locations of the dust grains
    data["angle"] = np.round(np.arctan(data["r"]/data["d"]), 3)

    # scatter angle is the angles at which each dust grain would actually
    # scatter a photon
    data["scatter angle"] = choices(weighted_scatter, k=nphotons)

    # return populated DataFrame
    return data


# %%
# ______________________CREATING PHOTON LOCATION DISTRIBUTION______________________
# create list of angles
thetas = np.arange(0, np.pi, step=10**-5)

# generate weighted list of possible scattering angles
weighted_scatter = generate_weighted_angles(thetas)

# set source location and number of photons to be generated
source = (180.0, 45.0)
NPHOTONS = 10**8

# generate dust layer/photons
dust = generate_dust(NPHOTONS)

# %%
# ___________________________ISOLATING SCATTERED PHOTONS___________________________
# isolate all photons that scattered to us
scattered = dust[dust["scatter angle"] == dust["angle"]]

# calculate true distances to the grains that scattered them
scattered["d_grain"] = scattered["d"]/np.cos(scattered["scatter angle"])

# calculate the light travel time from grains to us
speedoflight = 3
scattered["t"] = scattered["d_grain"]/speedoflight

# sort the lists by time—
# time_sorted = scattered.sort_values(by="t")

# %%
# ________________________________PLOTTING RESULTS_________________________________
# plot scattered photons
plt.scatter(scattered["ra"], scattered["dec"],
            c=scattered["t"], marker=".", linewidths=.25)

# set colormap and colorbar
plt.set_cmap('jet_r')
cbar = plt.colorbar()

# plot source
plt.plot(source[0], source[1], 'k*', markersize=8)

# create legend
plt.legend(labels=["source", "scattered photon"], loc=1, fontsize="small",
           handletextpad=0.1, markerscale=0.9, scatteryoffsets=[0.4])

# label axes and colorbar
plt.xlabel("Right Ascension")
plt.ylabel("Declination")
cbar.set_label("Photon Travel Time")

# display plot
plt.show()

# %%
