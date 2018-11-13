import astropy.io.ascii as ascii
import astropy.table
from astropy.coordinates import SkyCoord, Angle, SkyOffsetFrame, ICRS
from astropy import units as u
from astropy.table import Table
from astropy.io.ascii import FixedWidth
import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
import os
from glob import glob


def build_aspecs_catalog(initial_catalog=None):
    hdu_list = fits.open(initial_catalog)
    initial_catalog = hdu_list[1].data
    ra_dec = SkyCoord(initial_catalog['ra'] * u.deg, initial_catalog['dc'] * u.deg, frame='icrs')

    coords = []
    freqs = []

    with open(os.path.join("data", "ASPECS_lines.txt")) as data_file:
        # Read in the locations and
        for line in data_file:
            no_circle = line.split("(")[1]
            first_split = no_circle.split(",")
            ra = float(first_split[0])
            dec = float(first_split[1])
            coords.append(SkyCoord(ra*u.deg, dec*u.deg, frame='icrs'))
            frequency = first_split[2].split("{")[1]
            frequency = float(frequency.split("}")[0])
            freqs.append(frequency)

    coords = SkyCoord(coords)
    idx, d2d, d3d = coords.match_to_catalog_sky(ra_dec)
    print(idx)
    print(initial_catalog[idx])
    print("Number of Matches: " + str(len(idx)) + "/" + str(len(initial_catalog[idx])))
    print("Distances: ")
    for index, id in enumerate(idx):
        print("\nMatch: " + str(index))
        print("Distance: " + str(coords[index].separation(ra_dec[id]).arcsecond))
    return NotImplementedError

build_aspecs_catalog(os.path.join("data", "jacob_aspecs_catalog_fixed_magphys_jcb3.fits"))