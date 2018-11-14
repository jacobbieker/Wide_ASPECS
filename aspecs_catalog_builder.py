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
from astropy.coordinates import FK5


def build_aspecs_catalog(initial_catalog=None, dec_key='dec', ra_key='ra', frame='fk5'):
    hdu_list = fits.open(initial_catalog)
    initial_catalog = hdu_list[1].data
    ra_dec = SkyCoord(initial_catalog[ra_key] * u.deg, initial_catalog[dec_key] * u.deg, frame=frame)

    coords = []
    freqs = []

    with open(os.path.join("data", "ASPECS_lines.txt")) as data_file:
        # Read in the locations and
        for line in data_file:
            no_circle = line.split("(")[1]
            first_split = no_circle.split(",")
            ra = float(first_split[0])
            dec = float(first_split[1])
            coords.append(SkyCoord(ra*u.deg, dec*u.deg, frame=frame))
            frequency = first_split[2].split("{")[1]
            frequency = float(frequency.split("}")[0])
            freqs.append(frequency)

    coords = SkyCoord(coords)
    idx, d2d, d3d = coords.match_to_catalog_sky(ra_dec)
    print("\n----------------- Number of Matches: " + str(len(idx)) + "/" + str(len(initial_catalog[idx])))
    print("Distances: ")
    num_in_skelton = 0
    for index, id in enumerate(idx):
        print("\nMatch: " + str(index))
        print("Distance: " + str(coords[index].separation(ra_dec[id]).arcsecond))
        print("Location: " + str(coords[index]))
        try:
            print("Catalog ID: " + str(initial_catalog[id]['id']))
            print("In Skelton et al. Catalog: " + str(initial_catalog[id]['flag_3dh']))
        except:
            try:
                print("Catalog ID: " + str(initial_catalog[id]['unique_id']))
                print("Skelton et al. ID: " + str(initial_catalog[id]['skelton_id']))
            except:
                continue

    # Get the IDs of the matched values
    #catalog_ids = initial_catalog[idx]['id']


def compare_catalog_locations(roberto_catalog, initial_catalog, ra_key='ra', dec_key='dec', frame='fk5'):
    hdu_list = fits.open(initial_catalog)
    initial_catalog = hdu_list[1].data
    ra_dec = SkyCoord(initial_catalog[ra_key] * u.deg, initial_catalog[dec_key] * u.deg, frame=frame)
    hdu_list = fits.open(roberto_catalog)
    roberto_catalog = hdu_list[1].data
    roberto_ra_dec = SkyCoord(roberto_catalog['ra'] * u.deg, roberto_catalog['dc'] * u.deg, frame=frame)

    # Now compare the two catalogs to find matches
    idx, d2d, d3d = roberto_ra_dec.match_to_catalog_sky(ra_dec)
    # Matches less than 0.5 arc seconds
    less_than_5 = []
    # Less than 0.25 arc seconds
    less_than_25 = []
    # Less than 1 arc second
    less_than_1 = []
    for index, id in enumerate(idx):
        if roberto_ra_dec[index].separation(ra_dec[id]).arcsecond < 1:
            less_than_1.append([index, id])
        if roberto_ra_dec[index].separation(ra_dec[id]).arcsecond < 0.5:
            less_than_5.append([index, id])
        if roberto_ra_dec[index].separation(ra_dec[id]).arcsecond < 0.25:
            less_than_25.append([index, id])
            #print("\nMatch: " + str(index))
            #print("Location: " + str(roberto_ra_dec[index]))
            #print("Distance: " + str(roberto_ra_dec[index].separation(ra_dec[id]).arcsecond))
            try:
                if not roberto_catalog[index]['flag_3dh']:
                    print("MUSE Catalog Sep: " + str(initial_catalog[id]['skelton_sep']))
                    print("Difference Between MUSE and Roberto: " + str(initial_catalog[id]['skelton_sep'] - roberto_ra_dec[index].separation(ra_dec[id]).arcsecond))
                    print("Catalog ID: " + str(roberto_catalog[index]['id']))
                    print("In Skelton et al. Catalog: " + str(roberto_catalog[index]['flag_3dh']))
                    print("Skelton et al. ID: " + str(initial_catalog[id]['skelton_id']))
                    print("MUSE Catalog ID: " + str(initial_catalog[id]['unique_id']))
            except:
                continue

    print("Less than 1.0 arcseconds: " + str(len(less_than_1)))
    print("Less than 0.5 arcseconds: " + str(len(less_than_5)))
    print("Less than 0.25 arcseconds: " + str(len(less_than_25)))

compare_catalog_locations(os.path.join("data", "jacob_aspecs_catalog_fixed_magphys_jcb3.fits"), os.path.join("data", "MW_44fields_main_table_v1.0.fits"))
exit()
build_aspecs_catalog(os.path.join("data", "jacob_aspecs_catalog_fixed_magphys_jcb3.fits"), dec_key='dc')
build_aspecs_catalog(os.path.join("data", "MW_44fields_main_table_v1.0.fits"))