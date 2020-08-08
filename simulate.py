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
    data = pd.DataFrame()

    # generate specifed photons with random RA, Dec, and distance
    # (d is the distance to the dust grain along the line of sight
    # to the source)
    data["photons"] = np.arange(0, nphotons)
    data["ra"] = np.random.randint(0, 3600, size=data["photons"].shape)/10
    data["dec"] = np.random.randint(0, 900, size=data["photons"].shape)/10
    data["d"] = np.random.randint(0, 100, size=data["photons"].shape)

    # ra_diff and dec_diff are the differences in position of the
    # dust grains with respect to the source and along the line of
    # sight
    ra_diff, dec_diff = -data["ra"]+source[0], -data["dec"]+source[1]
    # r is the degree separations between the source and the dust
    # grains
    data["r"] = np.sqrt(np.power(ra_diff, 2) + np.power(dec_diff, 2))

    # scatter angle is the angles the photons need to scatter to
    # reach us, given the locations of the dust grains
    data["scatter angle"] = np.round(np.arctan(data["r"]/data["d"]), 3)

    # return populated DataFrame
    return data


def equation(angles):
    """Return an array containing the photon scattering probability at each 
    angle passed as input. 

    Parameters:\n
    angles: array-like\n
    \tThe angles for which scattering probabilities need to be calculated.

    Returns:\n
    p: array-like\n
    \tThe probabilities of scattering at each angle passed as input.
    """
    SD = 0.0292
    thetaoffset = angles + SD
    exponent = np.power(thetaoffset, 2)/-0.0017059

    i = np.exp(exponent)
    j = np.sin(thetaoffset)
    k = 22.523 * np.sqrt(2*np.pi) / SD

    almost_p = np.multiply(i, j)
    p = np.multiply(almost_p, k)
    return p


# %%
# ______________CREATING PHOTON LOCATION AND SCATTERING DISTRIBUTIONS______________
# set source location and number of photons to be generated
source = (180.0, 45.0)
NPHOTONS = 10**8

# generate dust layer/photons
dust = generate_dust(NPHOTONS)

# create distribution of angles
thetas = np.arange(0, np.pi, step=10**-5)

# generate the scattering probability for each angle
probs = equation(thetas)

# create list of potential scattering angles weighted by their likelihood
p_total = np.sum(probs)
p_list = np.divide(probs, p_total)
a = np.round(choices(thetas, weights=p_list, k=10**5), 3)
# np.histogram(a, bins=100)

# randomly choose the optimal scattering angle
sa = choice(a)

# isolate all photons that scattered to us
scattered = dust[dust["scatter angle"] == sa]

# calculate true distances to the grains that scattered them
scattered["d_grain"] = scattered["d"]/np.cos(sa)

# calculate the light travel time from grains to us
speedoflight = 3
scattered["t"] = scattered["d_grain"] / speedoflight

# sort the lists by time—
time_sorted = scattered.sort_values(by="t")

# %%
# ________________________________PLOTTING RESULTS_________________________________
# plot scattered photons
plt.scatter(time_sorted["ra"], time_sorted["dec"],
            c=time_sorted["t"], marker=".", linewidths=.25)

# set colormap and colorbar
plt.set_cmap('jet_r')
cbar = plt.colorbar()

# plot source
plt.plot(180, 45, 'k*', markersize=8)

# label axes and colorbar
plt.xlabel("Right Ascension")
plt.ylabel("Declination")
cbar.set_label("Photon Travel Time")

# display plot
plt.show()

# %%
