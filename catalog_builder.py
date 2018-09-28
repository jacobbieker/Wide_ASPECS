import astropy.io.ascii as ascii
import astropy.table
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
import os
from glob import glob

photometry_dir = os.path.join("data", "matched_jun18")

photometry_files = glob(os.path.join("data", "matched_jun18", "*.cat"))

full_catalog = None
tmp_catalog = None

for catalog_file in photometry_files:
    try:
        if tmp_catalog is None:
            tmp_catalog = ascii.read(catalog_file)
        else:
            full_catalog = astropy.table.join(tmp_catalog, ascii.read(catalog_file))
    except:
        print("Could not add: ")
        print(catalog_file)
        continue

ra_dec = full_catalog['ra', 'dc']

# The Franco dataset

ra = 53.118815 * u.deg
dec = -27.782889 * u.deg
ags1 = SkyCoord(ra, dec, frame='J2000')

ra = 53.063867 * u.deg
dec = -27.843792 * u.deg
ags2 = SkyCoord(ra, dec, frame='J2000')

ags3 = SkyCoord(53.148839*u.deg, -27.821192*u.deg, frame='J2000')

ags4 = SkyCoord(53.142778-27.827888)
ags5 = SkyCoord(53.158392-27.733607)
ags6 = SkyCoord(53.183458-27.776654)
ags7 = SkyCoord(53.082738-27.866577)
ags8 = SkyCoord(53.020356-27.779905)
ags9 = SkyCoord(53.092844-27.801330)
ags10 = SkyCoord(53.082118-27.767299)
ags11 = SkyCoord(53.108818-27.869055)
ags12 = SkyCoord(53.160634-27.776273)
ags13 = SkyCoord(53.131122-27.773194)
ags14 = SkyCoord(53.223156-27.826771)
ags15 = SkyCoord(53.074847-27.875880)
ags16 = SkyCoord(53.039724-27.784557)
ags17 = SkyCoord(53.079374-27.870770)
ags18 = SkyCoord(53.181355-27.777544)
ags19 = SkyCoord(53.108041-27.813610)
ags20 = SkyCoord(53.092365-27.826829)

ags21 = SkyCoord(53.070274-27.845586)
ags22 = SkyCoord(53.108695-27.848332)
ags23 = SkyCoord(53.086623-27.810272)


# TODO Convert to degress the ra and dec ones, then Skycoordinate

# TODO Make graphs of the redshift clustering, star mass, etc. using mean flux values as the cutoff

