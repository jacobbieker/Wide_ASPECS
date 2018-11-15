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

def get_aspecs_radec(frame='fk5'):
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
    return coords, freqs

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
    num_in_close = 0
    for index, id in enumerate(idx):
        if coords[index].separation(ra_dec[id]).arcsecond < 1.0:
            num_in_close += 1
            print("\nMatch: " + str(index))
            print("Distance: " + str(coords[index].separation(ra_dec[id]).arcsecond))
            print("Location (RA): " + str(coords[index].ra.hms))
            print("Location (Dec): " + str(coords[index].dec.hms))
            print("Location (Deg): " + str(coords[index]))
            try:
                print("Catalog ID: " + str(initial_catalog[id]['id']))
                print("In Skelton et al. Catalog: " + str(initial_catalog[id]['flag_3dh']))
            except:
                try:
                    print("Catalog ID: " + str(initial_catalog[id]['unique_id']))
                    print("Skelton et al. ID: " + str(initial_catalog[id]['skelton_id']))
                except:
                    continue
    print("Number Close to Catalog: ", num_in_close)

    # Get the IDs of the matched values
    #catalog_ids = initial_catalog[idx]['id']


def compare_catalog_locations(roberto_catalog, initial_catalog, ra_key='ra', dec_key='dec', frame='fk5'):
    hdu_list = fits.open(initial_catalog)
    initial_catalog = hdu_list[1].data
    ra_dec = SkyCoord(initial_catalog[ra_key] * u.deg, initial_catalog[dec_key] * u.deg, frame=frame)
    hdu_list = fits.open(roberto_catalog)
    rob_cat = hdu_list[1].data
    roberto_ra_dec = SkyCoord(rob_cat['ra'] * u.deg, rob_cat['dc'] * u.deg, frame=frame)

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
            print("\nMatch: " + str(index))
            print("Location: " + str(roberto_ra_dec[index]))
            print("Distance: " + str(roberto_ra_dec[index].separation(ra_dec[id]).arcsecond))
            try:
                print("MUSE Catalog Sep: " + str(initial_catalog[id]['skelton_sep']))
                print("Difference Between MUSE and Roberto: " + str(initial_catalog[id]['skelton_sep'] - roberto_ra_dec[index].separation(ra_dec[id]).arcsecond))
                print("Catalog ID: " + str(rob_cat[index]['id']))
                print("In Skelton et al. Catalog: " + str(rob_cat[index]['flag_3dh']))
                print("Skelton et al. ID: " + str(initial_catalog[id]['skelton_id']))
                print("MUSE Catalog ID: " + str(initial_catalog[id]['unique_id']))
            except:
                continue

    print("Less than 1.0 arcseconds: " + str(len(less_than_1)))
    print("Less than 0.5 arcseconds: " + str(len(less_than_5)))
    print("Less than 0.25 arcseconds: " + str(len(less_than_25)))

    roberto_catalog = Table.read(roberto_catalog, format='fits')
    # Now add those matches within a given constraint to the FITS file to get add another Z fit
    roberto_catalog['muse_wide_z'] = np.zeros(len(rob_cat['ra']))
    roberto_catalog['muse_wide_z_err'] = np.zeros(len(rob_cat['ra']))
    print(idx)
    for index, galaxy in enumerate(idx):
        # Change this if want better/worse coverage
        if roberto_ra_dec[index].separation(ra_dec[galaxy]).arcsecond < 1:
            roberto_catalog[galaxy]['muse_wide_z'] = initial_catalog[galaxy]['z']
            roberto_catalog[galaxy]['muse_wide_z_err'] = initial_catalog[galaxy]['z_err']

    return roberto_catalog

"""
build_aspecs_catalog(os.path.join("data", "jacob_aspecs_catalog_fixed_magphys_jcb3.fits"), dec_key='dc')
build_aspecs_catalog(os.path.join("data", "MW_44fields_main_table_v1.0.fits"))

compare_catalog_locations(os.path.join("data", "jacob_aspecs_catalog_fixed_magphys_jcb3.fits"), os.path.join("data", "MW_44fields_main_table_v1.0.fits"))
"""