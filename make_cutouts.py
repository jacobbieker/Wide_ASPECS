import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy.nddata import Cutout2D

import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
import os

from astropy.table import Table, join
from matplotlib.patches import Circle
from aspecs_catalog_builder import get_aspecs_radec


catalog_goodss = fits.open("/home/jacob/Research/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS")
skelton_goodss = catalog_goodss[1].data
print(catalog_goodss[0].header)
print(skelton_goodss.columns)
print(skelton_goodss['id'])

def create_cutout(ax, wcs_header, image, ra, dec, aspecs_ra, aspecs_dec):
    """
    :param ax: Matplotlib ax to use
    :param wcs_header: Image header with WCS info
    :param image: Image
    :param ra: RA coordiantes for center in J2000
    :param dec: DEC coordinates for center in J2000
    :param size: Size, in degrees, of the cutout
    :return:
    """

    w = wcs.WCS(wcs_header)
    center = SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
    aspecs_ra_dec_center = SkyCoord(aspecs_ra * u.deg, aspecs_dec * u.deg, frame='fk5')
    size = np.ceil(center.separation(aspecs_ra_dec_center).arcsecond * 2.)
    if size < 4.0:
        size = 4
    # then make an array cutout
    aspecs_cutout = Cutout2D(image, aspecs_ra_dec_center, size=size * u.arcsec, wcs=w)
    co = Cutout2D(image, center, size=size * u.arcsec, wcs=w)
    ax.imshow(co.data, origin='lower', cmap='gray')
    center_image = Circle((co.center_cutout[0], co.center_cutout[1]), 5, fill=False, color='r')
    ax.add_patch(center_image)

    aspecs_loc_x, aspecs_loc_y = co.to_cutout_position(aspecs_cutout.center_original)
    ax.errorbar(aspecs_loc_x, aspecs_loc_y, yerr=5, xerr=5, color='r')

    return ax


def create_ax_cutout(ax, name, fit_data, aspecs_coordinate, catalog_coordinate):
    ax = create_cutout(ax, fit_data[0].header, fit_data[0].data, catalog_coordinate.ra.deg, catalog_coordinate.dec.deg,
                       aspecs_coordinate.ra.deg, aspecs_coordinate.dec.deg)
    ax.set_title(name)
    ax.tick_params(direction='in', colors='w', bottom=True, top=True, left=True, right=True, labelbottom=True,
                   labeltop=False, labelleft=True, labelright=False)
    return ax


def create_aspecs_cutouts(aspecs_coordinates, fits_files, fits_names, wcs_data, catalog_coordinates, id, aspecs_freqs):
    # get size of file
    shape_file = int(np.ceil(np.sqrt(len(fits_files))))
    f = plt.figure(figsize=(20, 20))
    f.suptitle('ID: ' + str(id) + " ALMA Freq: " + str(aspecs_freqs[id]))
    for index, image in enumerate(fits_files):
        ax = f.add_subplot(shape_file, shape_file, index + 1, projection=wcs_data)
        create_ax_cutout(ax, fits_names[index], image, aspecs_coordinates, catalog_coordinates)

    return f


hdu_list = fits.open("data/jacob_aspecs_catalog_fixed_magphys_jcb3.fits")
print(hdu_list[1].columns)
full_catalog = hdu_list[1].data

from astropy.utils import data
from spectral_cube import SpectralCube

# ASPECS_Data Cubes
aspecs_a1_chn = SpectralCube.read("/home/jacob/Research/gs_A1_2chn.fits")
aspecs_a2_chn = SpectralCube.read("/home/jacob/Research/gs_A2_2chn.fits")

print(aspecs_a1_chn)
velo, dec, ra = aspecs_a1_chn.unitless.world[0,:,:]
print(dec[0])
alma_ra_dec = SkyCoord(ra, dec, frame='fk5')

aspecs_a1_chn_fits = fits.open("/home/jacob/Research/gs_A1_2chn.fits")

f160w_goodss = fits.open("/home/jacob/Research/goodss_3dhst_v4.0_f160w/goodss_3dhst.v4.0.F160W_orig_sci.fits")
w = wcs.WCS(f160w_goodss[0].header)
ax = plt.subplot(projection=w)
#ax.plot(aspecs_a1_chn[0,:,:], origin='lower')
#ax.imshow(image_data, cmap='gray')
aspecs_ra_dec, aspecs_freqs = get_aspecs_radec()
print(aspecs_ra_dec.ra.deg)
plt.scatter(full_catalog['ra'], full_catalog['dc'], transform=ax.get_transform('fk5'), label='All', s=2)
ax.scatter(aspecs_ra_dec.ra.deg, aspecs_ra_dec.dec.deg, transform=ax.get_transform('fk5'), s=50,
           edgecolor='black', facecolor='none')
ax.scatter(ra, dec, transform=ax.get_transform('fk5'), label='All', s=6)
plt.ylabel("DEC")
plt.xlabel("RA")
plt.savefig("test.png", dpi=300)
#plt.show()
exit()

roberto_muse = Table.read("roberto_catalog_muse.fits", format='fits')
test_roberto = Table.read("/home/jacob/Development/Wide_ASPECS/mapghys_in_nov2018_all.fits", format='fits')
# Add in RA and Dec to test_roberto

roberto_muse = join(test_roberto, roberto_muse, keys='id')

muse_catalog = fits.open(os.path.join("data", "MW_44fields_main_table_v1.0.fits"))[1].data
roberto_ra_dec = SkyCoord(roberto_muse['ra'] * u.deg, roberto_muse['dc'] * u.deg, frame='fk5')
muse_ra_dec = SkyCoord(muse_catalog['RA'] * u.deg, muse_catalog['DEC'] * u.deg, frame='fk5')

indicies_to_use = []
muse_ids_to_use = []

# For ones with same ID
from matplotlib.text import OffsetFrom

# For ones with off Z


f125w_goodss = fits.open("/home/jacob/Research/goodss_3dhst_v4.0_f125w/goodss_3dhst.v4.0.F125W_orig_sci.fits")
f140w_goodss = fits.open("/home/jacob/Research/goodss_3dhst_v4.0_f140w/goodss_3dhst.v4.0.F140W_orig_sci.fits")
f160w_goodss = fits.open("/home/jacob/Research/goodss_3dhst_v4.0_f160w/goodss_3dhst.v4.0.F160W_orig_sci.fits")
f435w_goodss = fits.open("/home/jacob/Research/goodss_3dhst_v4.0_f435w/goodss_3dhst.v4.0.F435W_orig_sci.fits")
f606w_goodss = fits.open("/home/jacob/Research/goodss_3dhst.v4.0.f606wcand/goodss_3dhst.v4.0.F606Wcand_orig_sci.fits")
f775w_goodss = fits.open("/home/jacob/Research/goodss_3dhst_v4.0_f775w/goodss_3dhst.v4.0.F775W_orig_sci.fits")
f850lp_goodss = fits.open("/home/jacob/Research/goodss_3dhst_v4.0_f850lp/goodss_3dhst.v4.0.F850LP_orig_sci.fits")
R_goodss = fits.open("/home/jacob/Research/GOODS-S_R/GOODS-S_R_sci.fits")
U38_goodss = fits.open("/home/jacob/Research/GOODS-S_WFI_U38/GOODS-S_WFI_U38_sci.fits")
V_goodss = fits.open("/home/jacob/Research/GOODS-S_WFI_V/GOODS-S_WFI_V_sci.fits")
B_goodss = fits.open("/home/jacob/Research/GOODS-S_WFI_B/GOODS-S_WFI_B_sci.fits")
# U_goodss = fits.open("/home/jacob/Research/GOODS-S_U/GOODS-S_U_sci.fits")
J_goodss = fits.open("/home/jacob/Research/GOODS-S_convJ/GOODS-S_convJ_sci.fits")
H_goodss = fits.open("/home/jacob/Research/GOODS-S_convH/GOODS-S_convH_sci.fits")
f814w_goodss = fits.open("/home/jacob/Research/goodss_3dhst.v4.0.f814wcand/goodss_3dhst.v4.0.F814Wcand_orig_sci.fits")
I_goodss = fits.open("/home/jacob/Research/GOODS-S_WFI_I/GOODS-S_WFI_I_sci.fits")

# Tenis Ones
tKs_goodss = fits.open("/home/jacob/Research/GOODS-S_tenisK/GOODS-S_tenisK_sci.fits")
tJ_goodss = fits.open("/home/jacob/Research/GOODS-S_tenisJ/GOODS-S_tenisJ_sci.fits")

# MUSYC
#ia427_goodss = fits.open("/home/jacob/Research/GOODS-S_IA427/GOODS-S_IA427_sci.fits")
#ia445_goodss = fits.open("/home/jacob/Research/GOODS-S_IA445/GOODS-S_IA445_sci.fits")
#ia505_goodss = fits.open("/home/jacob/Research/GOODS-S_IA505/GOODS-S_IA505_sci.fits")
#ia527_goodss = fits.open("/home/jacob/Research/GOODS-S_IA527/GOODS-S_IA527_sci.fits")
#ia550_goodss = fits.open("/home/jacob/Research/GOODS-S_IA550/GOODS-S_IA550_sci.fits")
#ia574_goodss = fits.open("/home/jacob/Research/GOODS-S_IA574/GOODS-S_IA574_sci.fits")
#ia624_goodss = fits.open("/home/jacob/Research/GOODS-S_IA624/GOODS-S_IA624_sci.fits")
#ia651_goodss = fits.open("/home/jacob/Research/GOODS-S_IA651/GOODS-S_IA651_sci.fits")
#ia679_goodss = fits.open("/home/jacob/Research/GOODS-S_IA679/GOODS-S_IA679_sci.fits")
#ia738_foodss = fits.open("/home/jacob/Research/GOODS-S_IA738/GOODS-S_IA738_sci.fits")
#ia797_foodss = fits.open("/home/jacob/Research/GOODS-S_IA797/GOODS-S_IA797_sci.fits")

# IRAC
irac1_goodss = fits.open("/home/jacob/Research/GOODS-S_SEDS1/GOODS-S_SEDS1_sci_sub.fits")
irac2_goodss = fits.open("/home/jacob/Research/GOODS-S_SEDS2/GOODS-S_SEDS2_sci_sub.fits")
irac3_goodss = fits.open("/home/jacob/Research/GOODS-S_irac3/GOODS-S_irac3_s1_sci.fits")
irac4_goodss = fits.open("/home/jacob/Research/GOODS-S_irac4/GOODS-S_irac4_s1_sci.fits")

# Add to the list with names

#ia_ones = [ia427_goodss, ia445_goodss, ia505_goodss, ia527_goodss, ia550_goodss, ia574_goodss,
#           ia624_goodss, ia651_goodss, ia679_goodss, ia738_foodss, ia797_foodss, ]
fits_files = [f125w_goodss, f140w_goodss, f160w_goodss, f435w_goodss, f606w_goodss, f775w_goodss, f850lp_goodss,
              f814w_goodss, R_goodss, U38_goodss, V_goodss, B_goodss, J_goodss, H_goodss, I_goodss, tKs_goodss,
              tJ_goodss, irac1_goodss, irac2_goodss,
              irac3_goodss, irac4_goodss]
ia_nmes = ["IA427", "IA445", "IA505", "IA527", "IA550", "IA574",
           "IA624", "IA651", "IA679", "IA738", "IA797", ]

fits_names = ["F125W", "F140W", "F160W", "F435W", "F606W", "F775W", "F850LP",
              "F814W", "R", "U38", "V", "B", "J", "H", "I", "tKs",
              "tJ",
              "IRAC1", "IRAC2", "IRAC3", "IRAC4"]
catalog_goodss = fits.open("/home/jacob/Research/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS")
skelton_goodss = catalog_goodss[1].data


def create_overlap_cutout(ax, wcs_header, image, muse, row, other_row, index, other_index):
    """
    :param ax: Matplotlib ax to use
    :param wcs_header: Image header with WCS info
    :param image: Image
    :param ra: RA coordiantes for center in J2000
    :param dec: DEC coordinates for center in J2000
    :param size: Size, in degrees, of the cutout
    :return:
    """

    w = wcs.WCS(wcs_header)
    muse_mask = (muse['UNIQUE_ID'] == int(other_row['muse_id']))
    center = SkyCoord(muse[muse_mask]['RA'] * u.deg, muse[muse_mask]['DEC'] * u.deg, frame='fk5')
    row_center = roberto_ra_dec[index]
    other_row_center = roberto_ra_dec[other_index]
    size = np.ceil(
        np.max([center.separation(row_center).arcsecond * 2., center.separation(other_row_center).arcsecond * 2.]))
    if size < 3.0:
        size = 3
    # then make an array cutout
    row_cutout = Cutout2D(image, row_center, size=size * u.arcsec, wcs=w)
    other_row_cutout = Cutout2D(image, other_row_center, size=size * u.arcsec, wcs=w)
    co = Cutout2D(image, center, size=size * u.arcsec, wcs=w)
    ax.imshow(co.data, origin='lower', cmap='gray')
    center_image = Circle((co.center_cutout[0], co.center_cutout[1]), 3, fill=False, color='r')
    ax.add_patch(center_image)
    ax.annotate(str(int(row['muse_id'])), xy=(co.center_cutout[0], co.center_cutout[1]), textcoords='offset pixels',
                xytext=(2, 1), color='r')

    aspecs_loc_x, aspecs_loc_y = co.to_cutout_position(row_cutout.center_original)
    aspecs_loc_x_other, aspecs_loc_y_other = co.to_cutout_position(other_row_cutout.center_original)
    first_image = Circle((aspecs_loc_x, aspecs_loc_y), 3, fill=False, color='g')
    ax.add_patch(first_image)
    ax.annotate(row['id'], xy=(aspecs_loc_x, aspecs_loc_y), textcoords='offset pixels', xytext=(2, 1), color='g')
    other_image = Circle((aspecs_loc_x_other, aspecs_loc_y_other), 3, fill=False, color='w')
    ax.add_patch(other_image)
    ax.annotate(other_row['id'], xy=(aspecs_loc_x_other, aspecs_loc_y_other), textcoords='offset pixels',
                xytext=(2, 1), color='w')

    return ax


def create_overlap_ax_cutout(ax, name, fit_data, aspecs_coordinate, catalog_coordinate, other_row, index, other_index):
    ax = create_overlap_cutout(ax, fit_data[0].header, fit_data[0].data, catalog_coordinate,
                               aspecs_coordinate, other_row, index, other_index)
    ax.set_title(name)
    ax.tick_params(direction='in', colors='w', bottom=True, top=True, left=True, right=True, labelbottom=True,
                   labeltop=False, labelleft=True, labelright=False)
    return ax

"""
row_mask = (roberto_muse['id'] == 51778) | (roberto_muse['id'] == 57545) |(roberto_muse['id'] == 62887) |(roberto_muse['id'] == 18816)

row_masked = roberto_muse[row_mask]
row_masked[0]['muse_id'] = 125009025
row_masked[1]['muse_id'] = 125009025

row_masked[2]['muse_id'] = 119002002
row_masked[3]['muse_id'] = 119002002

print(row_masked)

shape_file = int(np.ceil(np.sqrt(len(fits_files))))
f = plt.figure(figsize=(20, 20))
f.suptitle(
    'MUSE ID: ' + str(int(row_masked[0]['muse_id'])) + " Roberto IDs: " + str(row_masked[0]['id']) + "/" + str(row_masked[1]['id']))
for third_index, image in enumerate(fits_files):
    ax = f.add_subplot(shape_file, shape_file, third_index + 1, projection=w)
    create_overlap_ax_cutout(ax, fits_names[third_index], image, catalog_coordinate=muse_catalog, aspecs_coordinate=row_masked[0],
                             other_row=row_masked[1], index=46881,
                             other_index=46963)
plt.show()
f.savefig(str("Overlap_2_MUSE_Cutout_" + str(46881) + "MUSE" + str(row_masked[0]['muse_id']) + ".png"), dpi=300)
f.clf()

shape_file = int(np.ceil(np.sqrt(len(fits_files))))
f = plt.figure(figsize=(20, 20))
f.suptitle(
    'MUSE ID: ' + str(int(row_masked[0]['muse_id'])) + " Roberto IDs: " + str(row_masked[2]['id']) + "/" + str(row_masked[3]['id']))
for third_index, image in enumerate(fits_files):
    ax = f.add_subplot(shape_file, shape_file, third_index + 1, projection=w)
    create_overlap_ax_cutout(ax, fits_names[third_index], image, catalog_coordinate=muse_catalog, aspecs_coordinate=row_masked[2],
                             other_row=row_masked[3], index=50461,
                             other_index=50464)
plt.show()
f.savefig(str("Overlap_2_MUSE_Cutout_" + str(50461) + "MUSE" + str(row_masked[3]['muse_id']) + ".png"), dpi=300)
f.clf()
"""

def create_multi_overlap_cutout(ax, wcs_header, image, aspecs, matches, aspecs_index, ra_dec=roberto_ra_dec):
    """
    :param ax: Matplotlib ax to use
    :param wcs_header: Image header with WCS info
    :param image: Image
    :param ra: RA coordiantes for center in J2000
    :param dec: DEC coordinates for center in J2000
    :param size: Size, in degrees, of the cutout
    :return:
    """

    w = wcs.WCS(wcs_header)
    center = aspecs

    other_centers = []
    for coord in matches:
        other_centers.append(ra_dec[coord])
    size = 2
    cutouts = []
    for row_center in other_centers:
        # then make an array cutout
        cutouts.append(Cutout2D(image, row_center, size=size * u.arcsec, wcs=w))
    co = Cutout2D(image, center, size=size * u.arcsec, wcs=w)
    ax.imshow(co.data, origin='lower', cmap='gray')
    center_image = Circle((co.center_cutout[0], co.center_cutout[1]), 3, fill=False, color='r')
    ax.add_patch(center_image)
    #ax.annotate(aspecs_index, xy=(co.center_cutout[0], co.center_cutout[1]), textcoords='offset pixels',
    #            xytext=(2, 1), color='r')

    for idx, cutout in enumerate(cutouts):
        aspecs_loc_x, aspecs_loc_y = co.to_cutout_position(cutout.center_original)
        first_image = Circle((aspecs_loc_x, aspecs_loc_y), 3, fill=False, color='g')
        ax.add_patch(first_image)
        ax.annotate(freqs[matches[idx]], xy=(aspecs_loc_x, aspecs_loc_y), textcoords='offset pixels', xytext=(2, 1), color='g')

    return ax


def create_multi_overlap_ax_cutout(ax, name, fit_data, catalog_coordinate, matches, index, ra_dec=roberto_ra_dec):
    ax = create_multi_overlap_cutout(ax, fit_data[0].header, fit_data[0].data, aspecs=catalog_coordinate,
                               matches=matches, aspecs_index=index, ra_dec=ra_dec)
    ax.set_title(name)
    ax.tick_params(direction='in', colors='w', bottom=True, top=True, left=True, right=True, labelbottom=True,
                   labeltop=False, labelleft=True, labelright=False)
    return ax

coords = []
freqs = []

with open(os.path.join("data", "ASPECS_lines.txt")) as data_file:
    # Read in the locations and
    for line in data_file:
        no_circle = line.split("(")[1]
        first_split = no_circle.split(",")
        ra = float(first_split[0])
        dec = float(first_split[1])
        coords.append(SkyCoord(ra*u.deg, dec*u.deg, frame='fk5'))
        frequency = first_split[2].split("{")[1]
        frequency = float(frequency.split("}")[0])
        freqs.append(frequency)

coords = SkyCoord(coords)

idx, d2d, d3d = coords.match_to_catalog_sky(roberto_ra_dec)

aspecs_matches = [[] for _ in range(92)]
back_match = {}

for index, id in enumerate(idx):
    if coords[index].separation(roberto_ra_dec[id]).arcsecond < 1.0:
        aspecs_matches[index].append(id)
        if id in back_match.keys():
            back_match[id].append(index)
        else:
            back_match[id] = [index]

# Now have the matches, plot them on the sky



"""
for idx, sublist in enumerate(aspecs_matches):
    if len(sublist) > 0:
        # Make the cutouts
        shape_file = int(np.ceil(np.sqrt(len(fits_files))))
        f = plt.figure(figsize=(20, 20))
        f.suptitle(
            'ASPECS ID: ' + str(idx) + " Number Matches: " + str(len(sublist)))
        for third_index, image in enumerate(fits_files):
            ax = f.add_subplot(shape_file, shape_file, third_index + 1, projection=w)
            create_multi_overlap_ax_cutout(ax, fits_names[third_index], image, catalog_coordinate=coords[idx], matches=sublist, index=idx)
        plt.show()
        f.savefig(str("Overlap_ASPECS_Sources_Cutout_" + str(idx) + ".png"), dpi=300)
        f.clf()

for key, values in back_match.items():
    if len(values) > 0:
        # Make the cutouts
        shape_file = int(np.ceil(np.sqrt(len(fits_files))))
        f = plt.figure(figsize=(20, 20))
        test_mask = (test_roberto['id'] == roberto_muse[key]['id'])
        freq_valus = []
        for value in values:
            print("Value: ", value)
            print("Freq: ", freqs[value])
            freq_valus.append(freqs[value])
        f.suptitle(
            'Roberto ID: ' + str(key) + "Z: " + str(test_roberto[test_mask]['z']) + " Matches: " + str(freq_valus))
        for third_index, image in enumerate(fits_files):
            ax = f.add_subplot(shape_file, shape_file, third_index + 1, projection=w)
            create_multi_overlap_ax_cutout(ax, fits_names[third_index], image, catalog_coordinate=roberto_ra_dec[key], matches=values, index=idx, ra_dec=coords)
        # plt.show()
        f.savefig(str("Overlap_ASPECS_From_Roberto_Sources_Cutout_" + str(key) + ".png"), dpi=300)
        f.clf()
        plt.close()
# Now work backwards
exit()

for index, row in enumerate(roberto_muse):
    if row['id'] in [46881, 46963, 50461, 50464]:
        if row['muse_id'] > 0:
            for other_index, other_row in enumerate(roberto_muse):
                # Ones before current index should already have been found
                if int(other_row['muse_id']) == int(row['muse_id']) and int(other_row['id']) != int(row['id']):
                    shape_file = int(np.ceil(np.sqrt(len(fits_files))))
                    f = plt.figure(figsize=(20, 20))
                    f.suptitle(
                        'MUSE ID: ' + str(int(row['muse_id'])) + " Roberto IDs: " + str(row['id']) + "/" + str(other_row['id']))
                    for third_index, image in enumerate(fits_files):
                        ax = f.add_subplot(shape_file, shape_file, third_index + 1, projection=w)
                        create_overlap_ax_cutout(ax, fits_names[third_index], image, catalog_coordinate=muse_catalog, aspecs_coordinate=row, other_row=other_row, index=index,
                                                 other_index=other_index)
                    f.savefig(str("Overlap_2_MUSE_Cutout_" + str(index) + "MUSE" + str(row['muse_id']) + ".png"), dpi=300)
                    f.clf()
                    plt.close()
exit()
"""
aspecs_ra_dec, aspecs_freqs = get_aspecs_radec()
ra_dec = SkyCoord(skelton_goodss['ra'] * u.deg, skelton_goodss['dec'] * u.deg, frame='fk5')

w = wcs.WCS(f160w_goodss[0].header)
ax = plt.subplot(projection=w)
ax.imshow(f160w_goodss[0].data, cmap='gray')
ax.scatter(aspecs_ra_dec.ra.deg, aspecs_ra_dec.dec.deg, transform=ax.get_transform('fk5'), s=100,
           edgecolor='white', facecolor='none')
plt.show()
idx, d2d, d3d = aspecs_ra_dec.match_to_catalog_sky(ra_dec)

for third_index, coord in enumerate(aspecs_ra_dec):
    f = create_aspecs_cutouts(coord, fits_files, fits_names, wcs_data=w, catalog_coordinates=ra_dec[idx[third_index]],
                              id=third_index, aspecs_freqs=aspecs_freqs)
    f.savefig(str("Skelton_ASPECS_Cutout_" + str(third_index) + "_Freq_" + str(
        aspecs_freqs[third_index]) + "_Sep_{:0.3f}" + ".png").format(
        coord.separation(ra_dec[idx[third_index]]).arcsecond), dpi=300)
    f.clf()
    plt.close()

# Now do it with Full Catalog
catalog_goodss = fits.open("data/jacob_aspecs_catalog_fixed_magphys_jcb3.fits")
full_goodss = catalog_goodss[1].data
aspecs_ra_dec, aspecs_freqs = get_aspecs_radec()
ra_dec = SkyCoord(full_goodss['ra'] * u.deg, full_goodss['dc'] * u.deg, frame='fk5')
idx, d2d, d3d = aspecs_ra_dec.match_to_catalog_sky(ra_dec)

for third_index, coord in enumerate(aspecs_ra_dec):
    f = create_aspecs_cutouts(coord, fits_files, fits_names, wcs_data=w, catalog_coordinates=ra_dec[idx[third_index]],
                              id=third_index, aspecs_freqs=aspecs_freqs)
    f.savefig(str("Full_ASPECS_Cutout_" + str(third_index) + "_Freq_" + str(
        aspecs_freqs[third_index]) + "_Sep_{:0.3f}" + ".png").format(
        coord.separation(ra_dec[idx[third_index]]).arcsecond), dpi=300)
    f.clf()
    plt.close()

# Now do it with the Roberto's Catalog, getting the Ra and Dec from the full catalog
roberto_catalog = fits.open("data/magphys_in.fits")
spec_ids = full_goodss['id']
diff_ids = []
roberto_catalog = roberto_catalog[1].data
for id in roberto_catalog['id']:
    if id in spec_ids:
        diff_ids.append(id)
# now create smaller catalog with full catalog info

rows_to_use = []
for third_index, row in enumerate(full_goodss):
    if row['id'] in diff_ids:
        rows_to_use.append(third_index)

smaller_catalog = full_goodss[rows_to_use]
aspecs_ra_dec, aspecs_freqs = get_aspecs_radec()
ra_dec = SkyCoord(smaller_catalog['ra'] * u.deg, smaller_catalog['dc'] * u.deg, frame='fk5')
idx, d2d, d3d = aspecs_ra_dec.match_to_catalog_sky(ra_dec)

for third_index, coord in enumerate(aspecs_ra_dec):
    f = create_aspecs_cutouts(coord, fits_files, fits_names, wcs_data=w, catalog_coordinates=ra_dec[idx[third_index]],
                              id=third_index, aspecs_freqs=aspecs_freqs)
    f.savefig(str("Roberto_ASPECS_Cutout_" + str(third_index) + "_Freq_" + str(
        aspecs_freqs[third_index]) + "_Sep_{:0.3f}" + ".png").format(
        coord.separation(ra_dec[idx[third_index]]).arcsecond), dpi=300)
    f.clf()
    plt.close()
