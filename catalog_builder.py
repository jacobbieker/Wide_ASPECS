import astropy.io.ascii as ascii
import astropy.table
from astropy.coordinates import SkyCoord, Angle, SkyOffsetFrame, ICRS
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
ra_dec = SkyCoord(ra_dec['ra'] * u.deg, ra_dec['dc'] * u.deg, frame='icrs')

# The Franco dataset

ra = 53.118815 * u.deg
dec = -27.782889 * u.deg
ags1 = SkyCoord(ra, dec, frame='icrs')

ra = 53.063867 * u.deg
dec = -27.843792 * u.deg
ags2 = SkyCoord(ra, dec, frame='icrs')

ags3 = SkyCoord(53.148839*u.deg, -27.821192*u.deg, frame='icrs')
ags4 = SkyCoord(53.142778*u.deg, -27.827888*u.deg, frame='icrs')
ags5 = SkyCoord(53.158392*u.deg, -27.733607*u.deg, frame='icrs')
ags6 = SkyCoord(53.183458*u.deg, -27.776654*u.deg, frame='icrs')
ags7 = SkyCoord(53.082738*u.deg, -27.866577*u.deg, frame='icrs')
ags8 = SkyCoord(53.020356*u.deg, -27.779905*u.deg, frame='icrs')
ags9 = SkyCoord(53.092844*u.deg, -27.801330*u.deg, frame='icrs')
ags10 = SkyCoord(53.082118*u.deg, -27.767299*u.deg, frame='icrs')
ags11 = SkyCoord(53.108818*u.deg, -27.869055*u.deg, frame='icrs')
ags12 = SkyCoord(53.160634*u.deg, -27.776273*u.deg, frame='icrs')
ags13 = SkyCoord(53.131122*u.deg, -27.773194*u.deg, frame='icrs')
ags14 = SkyCoord(53.223156*u.deg, -27.826771*u.deg, frame='icrs')
ags15 = SkyCoord(53.074847*u.deg, -27.875880*u.deg, frame='icrs')
ags16 = SkyCoord(53.039724*u.deg, -27.784557*u.deg, frame='icrs')
ags17 = SkyCoord(53.079374*u.deg, -27.870770*u.deg, frame='icrs')
ags18 = SkyCoord(53.181355*u.deg, -27.777544*u.deg, frame='icrs')
ags19 = SkyCoord(53.108041*u.deg, -27.813610*u.deg, frame='icrs')
ags20 = SkyCoord(53.092365*u.deg, -27.826829*u.deg, frame='icrs')

ags21 = SkyCoord(53.070274*u.deg, -27.845586*u.deg, frame='icrs')
ags22 = SkyCoord(53.108695*u.deg, -27.848332*u.deg, frame='icrs')
ags23 = SkyCoord(53.086623*u.deg, -27.810272*u.deg, frame='icrs')


# TODO Convert to degress the ra and dec ones, then Skycoordinate

franco_1mm = SkyCoord([ags1, ags2, ags3, ags4, ags5, ags6, ags7, ags8, ags9, ags10, ags11, ags12, ags13, ags14, ags15, ags16, ags17, ags18, ags19, ags20, ags21, ags22, ags23])

idx, d2d, d3d = franco_1mm.match_to_catalog_sky(ra_dec)

print(idx)
print("----------------------------------------------------------------------------------------------------------------")
print(d2d)

# Get the stuff now

i = 0
for galaxy in idx:
    print(full_catalog[galaxy])
    print(ra_dec[galaxy].ra.hms)
    print(d2d[i])
    i += 1
# TODO Make graphs of the redshift clustering, star mass, etc. using mean flux values as the cutoff

