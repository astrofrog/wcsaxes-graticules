import os
import glob

import numpy as np

from astropy.wcs import WCS
from astropy.io import fits

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection

from graticules import get_lon_lat_path


class FITSWCSWrapper(object):
    def __init__(self, wcs):
        self.wcs = wcs
    def world2pix(self, world):
        return self.wcs.wcs_world2pix(world, 1)
    def pix2world(self, pixel):
        return self.wcs.wcs_pix2world(pixel, 1)


for filename in glob.glob(os.path.join('data', '*.fits')):

    print(filename)

    # Read in header and zoom out
    header = fits.getheader(filename)
    header['CDELT1'] *= 50
    header['CDELT2'] *= 50
    header['CRPIX1'] = header['NAXIS1'] / 2.
    header['CRPIX2'] = header['NAXIS2'] / 2.

    # Parse WCS transformation
    wcs = FITSWCSWrapper(WCS(header))

    # Set range of coordinates to draw grid
    wmin = [-179.9999, -89.9999]
    wmax = [+179.9999, +89.9999]

    # Set transformation
    trans = wcs

    # Create figure
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
    ax.set_xlim(0.5, header['NAXIS1'] + 0.5)
    ax.set_ylim(0.5, header['NAXIS2'] + 0.5)

    ax.text(0.5, 0.9, os.path.basename(filename), ha='center', transform=ax.transAxes, size=18)

    # Define grid lines
    NG = 18
    N = 1000

    paths = []

    lon = np.linspace(wmin[0], wmax[0], N)
    for latval in np.linspace(wmin[1], wmax[1], NG):
        lat = np.repeat(latval, N)
        lon_lat = np.vstack([lon, lat]).transpose()
        paths.append(get_lon_lat_path(ax, trans, lon_lat))

    ax.add_collection(PathCollection(paths, edgecolors='purple', facecolors='none', alpha=0.4))

    paths = []

    lat = np.linspace(wmin[1], wmax[1], N)
    for lonval in np.linspace(wmin[0], wmax[0], NG)[1:]:
        lon = np.repeat(lonval, N)
        lon_lat = np.vstack([lon, lat]).transpose()
        paths.append(get_lon_lat_path(ax, trans, lon_lat))

    ax.add_collection(PathCollection(paths, edgecolors='b', facecolors='none', alpha=0.4))

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    fig.savefig(filename.replace('.fits', '.png'))
    plt.close(fig)
