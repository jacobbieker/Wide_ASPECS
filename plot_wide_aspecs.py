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

catalog_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS")
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
roberto_muse = hdu_list[1].data

from astropy.utils import data
from spectral_cube import SpectralCube

# ASPECS_Data Cubes
# aspecs_a1_chn = SpectralCube.read("/home/jacob/Research/gs_A1_2chn.fits")
# aspecs_a2_chn = SpectralCube.read("/home/jacob/Research/gs_A2_2chn.fits")

# plt.imshow(aspecs_a1_chn.unitless[1,:,:])
# plt.show()

# print(aspecs_a1_chn)
#summed_cub = fits.open("/home/jacob/Research/Wide_ASPECS/Data/gs_C1_2chn.fits")
#summed_cub = summed_cub[0].data
#summed_cub = np.reshape(summed_cub, (480, 2048, 2048))
#print(summed_cub.shape)
#summed_cub = np.sum(summed_cub, axis=0)
#plt.imshow(summed_cub)
#plt.show()

f160w_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/goodss_3dhst_v4.0_f160w/goodss_3dhst.v4.0.F160W_orig_sci.fits")
#w = wcs.WCS(f160w_goodss[0].header)

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


f125w_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/goodss_3dhst_v4.0_f125w/goodss_3dhst.v4.0.F125W_orig_sci.fits")
f140w_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/goodss_3dhst_v4.0_f140w/goodss_3dhst.v4.0.F140W_orig_sci.fits")
f160w_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/goodss_3dhst_v4.0_f160w/goodss_3dhst.v4.0.F160W_orig_sci.fits")
f435w_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/goodss_3dhst_v4.0_f435w/goodss_3dhst.v4.0.F435W_orig_sci.fits")
f606w_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/goodss_3dhst.v4.0.f606wcand/goodss_3dhst.v4.0.F606Wcand_orig_sci.fits")
f775w_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/goodss_3dhst_v4.0_f775w/goodss_3dhst.v4.0.F775W_orig_sci.fits")
f850lp_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/goodss_3dhst_v4.0_f850lp/goodss_3dhst.v4.0.F850LP_orig_sci.fits")
R_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/GOODS-S_R/GOODS-S_R_sci.fits")
U38_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/GOODS-S_WFI_U38/GOODS-S_WFI_U38_sci.fits")
V_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/GOODS-S_WFI_V/GOODS-S_WFI_V_sci.fits")
B_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/GOODS-S_WFI_B/GOODS-S_WFI_B_sci.fits")
# U_goodss = fits.open("/home/jacob/Research/GOODS-S_U/GOODS-S_U_sci.fits")
J_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/GOODS-S_convJ/GOODS-S_convJ_sci.fits")
H_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/GOODS-S_convH/GOODS-S_convH_sci.fits")
f814w_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/goodss_3dhst.v4.0.f814wcand/goodss_3dhst.v4.0.F814Wcand_orig_sci.fits")
I_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/GOODS-S_WFI_I/GOODS-S_WFI_I_sci.fits")

# Tenis Ones
tKs_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/GOODS-S_tenisK/GOODS-S_tenisK_sci.fits")
tJ_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/GOODS-S_tenisJ/GOODS-S_tenisJ_sci.fits")

# MUSYC
# ia427_goodss = fits.open("/home/jacob/Research/GOODS-S_IA427/GOODS-S_IA427_sci.fits")
# ia445_goodss = fits.open("/home/jacob/Research/GOODS-S_IA445/GOODS-S_IA445_sci.fits")
# ia505_goodss = fits.open("/home/jacob/Research/GOODS-S_IA505/GOODS-S_IA505_sci.fits")
# ia527_goodss = fits.open("/home/jacob/Research/GOODS-S_IA527/GOODS-S_IA527_sci.fits")
# ia550_goodss = fits.open("/home/jacob/Research/GOODS-S_IA550/GOODS-S_IA550_sci.fits")
# ia574_goodss = fits.open("/home/jacob/Research/GOODS-S_IA574/GOODS-S_IA574_sci.fits")
# ia624_goodss = fits.open("/home/jacob/Research/GOODS-S_IA624/GOODS-S_IA624_sci.fits")
# ia651_goodss = fits.open("/home/jacob/Research/GOODS-S_IA651/GOODS-S_IA651_sci.fits")
# ia679_goodss = fits.open("/home/jacob/Research/GOODS-S_IA679/GOODS-S_IA679_sci.fits")
# ia738_foodss = fits.open("/home/jacob/Research/GOODS-S_IA738/GOODS-S_IA738_sci.fits")
# ia797_foodss = fits.open("/home/jacob/Research/GOODS-S_IA797/GOODS-S_IA797_sci.fits")

# IRAC
irac1_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/GOODS-S_SEDS1/GOODS-S_SEDS1_sci_sub.fits")
irac2_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/GOODS-S_SEDS2/GOODS-S_SEDS2_sci_sub.fits")
irac3_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/GOODS-S_irac3/GOODS-S_irac3_s1_sci.fits")
irac4_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/GOODS-S_irac4/GOODS-S_irac4_s1_sci.fits")

# Add to the list with names

# ia_ones = [ia427_goodss, ia445_goodss, ia505_goodss, ia527_goodss, ia550_goodss, ia574_goodss,
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
catalog_goodss = fits.open("/home/jacob/Research/Wide_ASPECS/Historical_Data/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS")
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


def create_multi_overlap_cutout(ax, wcs_header, image, aspecs, matches, ra_dec=roberto_ra_dec, rob_z=0):
    """
    :param ax: Matplotlib ax to use
    :param wcs_header: Image header with WCS info
    :param image: Image
    :param ra: RA coordiantes for center in J2000
    :param dec: DEC coordinates for center in J2000
    :param size: Size, in degrees, of the cutout
    :return:
    :return:
    """

    w = wcs.WCS(wcs_header)
    center = aspecs

    other_centers = []
    for coord in matches:
        other_centers.append(ra_dec[coord])
    size = 3
    cutouts = []
    for row_center in other_centers:
        # then make an array cutout
        cutouts.append(Cutout2D(image, row_center, size=size * u.arcsec, wcs=w))
    co = Cutout2D(image, center, size=size * u.arcsec, wcs=w)
    ax.imshow(co.data, origin='lower', cmap='gray_r')
    center_image = Circle((co.center_cutout[0], co.center_cutout[1]), 3, fill=False, color='r')
    ax.add_patch(center_image)
    if rob_z > 0:
        ax.annotate(str(np.round(rob_z,3)), xy=(co.center_cutout[0], co.center_cutout[1]), textcoords='offset pixels',
                    xytext=(2, 1), color='r')

    for idx, cutout in enumerate(cutouts):
        aspecs_loc_x, aspecs_loc_y = co.to_cutout_position(cutout.center_original)
        first_image = Circle((aspecs_loc_x, aspecs_loc_y), 3, fill=False, color='g')
        ax.add_patch(first_image)
        #ax.annotate(freqs[matches[idx]], xy=(aspecs_loc_x, aspecs_loc_y), textcoords='offset pixels', xytext=(5, 20),
        #            color='g')
        #ax.annotate(np.round(z_s[matches[idx]],3), xy=(aspecs_loc_x, aspecs_loc_y), textcoords='offset pixels', xytext=(5, -20),
        #            color='orange')

    return ax



def create_multi_overlap_ax_cutout(ax, name, fit_data, catalog_coordinate, matches, ra_dec=roberto_ra_dec, rob_z=0):
    ax = create_multi_overlap_cutout(ax, fit_data[0].header, fit_data[0].data, aspecs=catalog_coordinate,
                                     matches=matches, ra_dec=ra_dec, rob_z=rob_z)
    ax.set_title(name)
    ax.tick_params(direction='in', colors='w', bottom=True, top=True, left=True, right=True, labelbottom=True,
                   labeltop=False, labelleft=True, labelright=False)
    return ax


def create_multi_matches_ax_cutout(ax, name, fit_data, catalog_coordinates, matches, ra_dec=roberto_ra_dec):
    """
    Makes a plot of all the matches on the same wavelength cutout
    :param ax:
    :param name:
    :param fit_data: The FIT image to use
    :param catalog_coordinates: Coordinates for the catalog matches
    :param matches: The coordinates of the matches
    :param redshifts: The redshifts in tuple format (CO, match)
    :param ra_dec:
    :return:
    """

    ax.set_title(name)
    ax.tick_params(direction='in', colors='w', bottom=True, top=True, left=True, right=True, labelbottom=True,
               labeltop=False, labelleft=True, labelright=False)
    return ax

def convert_to_rest_frame_ghz(z, ghz):
    """
    Take a measured GHz value, and calculates the restframe GHz value based on the given z of the matched galaxy
    :param z:
    :param ghz:
    :return:
    """

    # First step is to convert to nm rom rest rame GHz

    nm = (ghz * u.GHz).to(u.nm, equivalencies=u.spectral())

    # Second step is to convert from nm observed to restframe nm

    # Obseved/(z+1) = emitted

    nm_emitted = nm / (z + 1)
    print("Nm Emitted: {}, Z: {}".format(nm_emitted, z))

    # third step is to convert from restframe nm back to GHz

    final_ghz = (nm_emitted).to(u.GHz, equivalencies=u.spectral())

    return final_ghz

#aspecs_lines = Table.read("ASPECS_Line_Candidates_Z44_Total_Z_Limit.txt", format="ascii", header_start=0, data_start=1)


aspecs_lines = Table.read("/home/jacob/Development/Wide_ASPECS/independent/ASPECS_Line_Candidates_all_closest_Sep_1.5_SN_9.0.ecsv", format='ascii.ecsv')

#sn_cut = 9.5

#aspecs_lines = aspecs_lines[aspecs_lines['rsnrrbin'] > sn_cut]

transitions = {"1-0": [0.0030, 0.3694, 115.271],
               "2-1": [1.0059, 1.7387, 230.538],
               "3-2": [2.0088, 3.1080, 345.796],
               "4-3": [3.0115, 4.4771, 461.041],
               "5-4": [4.0142, 5.8460, 576.268],
               "6-5": [5.0166, 7.2146, 691.473],
               "7-6": [6.0188, 8.5829, 806.652],
               "C1 1-0": [3.2823, 4.8468, 492.161],
               "C1 2-1": [6.0422, 8.6148, 809.342]}

coords = SkyCoord(aspecs_lines['RA (J2000)'] * u.deg, aspecs_lines['DEC (J2000)'] * u.deg, frame='fk5')
freqs = aspecs_lines['Observed CO (GHz)']
z_s = aspecs_lines['Z (CO)']
rob_z = aspecs_lines['Z (Matched)']

# Now plot all Radio Sources and see what is around them for all ones without a match

for index, row in enumerate(aspecs_lines):
    # Make the cutouts
    shape_file = int(np.ceil(np.sqrt(len(fits_files))))
    f = plt.figure(figsize=(20, 20))
    # no counterpart
    distances = [0]
    freq_valus = [np.round(row['Observed CO (GHz)'], 3)]
    rest_frame_ghz = [np.round(row['Restframe CO (GHz)'], 3)]
    f.suptitle(
        " Z: " + str(row['Z (CO)']) + " Delta_Z: " + str(
            row['Delta Z']) + " Roberto ID: " + str(row['Roberto ID']) + " Separation: " + str(row['Separation (Arcsecond)']) + " S/N: " + str(row['S/N']) +
        " Observed Freq: " + str(freq_valus) + "\n Spec Z: " + str(
            row['Spec Z']) + " Transition: " + str(row['Transition']) +
        " RA: " + str(row['RA (J2000)']) + " DEC: " + str(row['DEC (J2000)']) +
        "\n Rest Frame GHz: " + str(rest_frame_ghz) + " MStar: " + str(row['Log(M*)']) + " SFR: " + str(row['Log(SFR)']))
    for third_index, image in enumerate(fits_files):
        wcs_header = image[0].header
        w = wcs.WCS(wcs_header)
        ax = f.add_subplot(shape_file, shape_file, third_index + 1, projection=w)
        create_multi_overlap_ax_cutout(ax, fits_names[third_index], image,
                                       catalog_coordinate=coords[index],
                                       matches=[index], ra_dec=coords)
    f.savefig(str("May_Output/matched/ASPECS_Cutout_NoCounter_Sep1.5_SN9.0_" + str(index) + ".png"), dpi=300)
    f.clf()
    plt.close()



#idx, d2d, d3d = coords.match_to_catalog_sky(roberto_ra_dec)

from astropy.coordinates import match_coordinates_sky, search_around_sky

idxc, idxcatalog, d2d, d3d = search_around_sky(coords, roberto_ra_dec, 1.5 * u.arcsecond)

def make_skycoords(source, ra='ra', dec='dec', distance=None):
    """
    Makes and returns a SkyCoord array from given source
    :param source: Source with information
    :param ra: Key for RA
    :param dec: Key for Dec
    :return: SkyCoord list
    """
    try:
        skycoords = SkyCoord(source[ra] * u.deg, source[dec] * u.deg, frame='icrs')
    except:
        skycoords = SkyCoord(source[ra], source[dec], unit=(u.hour, u.deg), frame='icrs')

    return skycoords

def load_table(ascii_table, header=0, start=1):
    ascii_table_data = Table.read(ascii_table, format="ascii", header_start=header, data_start=start)
    return ascii_table_data

#aspecs_lines = load_table("/home/jacob/Development/Wide_ASPECS/check_SN/line_search_P3_wa_crop.out")

#sn_cut = 9.5
#aspecs_lines = aspecs_lines[aspecs_lines['rsnrrbin'] > sn_cut]

#coords = make_skycoords(aspecs_lines,  ra='rra', dec='rdc')

#idxc, idxcatalog, d2d, d3d = search_around_sky(coords, roberto_ra_dec, 1.5 * u.arcsecond)

aspecs_matches = [[] for _ in range(len(aspecs_lines))]
back_match = {}
z_specs = {}

for index, id in enumerate(idxc):
    if coords[idxc[index]].separation(roberto_ra_dec[idxcatalog[index]]).arcsecond < 1.5:# and np.abs(aspecs_lines[idxc[index]]['Z (CO)'] - roberto_muse[idxcatalog[index]]['z_1']) < 0.3:
        test_mask = (roberto_muse['id'] == roberto_muse[idxcatalog[index]]['id'])
        test_rob = roberto_muse[test_mask]
        spec_z_mask = (test_rob["z_spec_3dh"] > 0.001) | (test_rob["zm_vds"] > 0.001) | (
                test_rob["zm_coeS"] > 0.001) | (test_rob['muse_wide_z'] > 0.0001) \
                      | (test_rob["zs_mor"] > 0.001) | (test_rob["zm_ina"] > 0.001) | (test_rob["zm_her"] > 0.001)
        #if int(aspecs_lines[idxc[index]]["Roberto ID"]) == int(roberto_muse[idxcatalog[index]]['id']):
        aspecs_matches[idxc[index]].append(roberto_muse[idxcatalog[index]]['id'])
        if roberto_muse[idxcatalog[index]]['id'] in back_match.keys():
            back_match[roberto_muse[idxcatalog[index]]['id']].append(idxc[index])
            z_specs[roberto_muse[idxcatalog[index]]['id']].append(len(test_rob[spec_z_mask]))
        else:
            back_match[roberto_muse[idxcatalog[index]]['id']] = [idxc[index]]
            z_specs[roberto_muse[idxcatalog[index]]['id']] = [len(test_rob[spec_z_mask])]

'''
for fits_index in range(len(fits_files)):
    f = plt.figure(figsize=(20, 20))
    f.suptitle(
        'Continuum Lines Matched To Galaxies')
    image_index = 0
    for key, values in back_match.items():
        if len(values) > 0:
            # Make the cutouts
            shape_file = int(np.ceil(np.sqrt(len(fits_files))))
            test_mask = (roberto_muse['id'] == key)
            roberto_ra_dec_index = 1e30
            for index, i in enumerate(roberto_muse):
                if i['id'] == key:
                    roberto_ra_dec_index = index
            for third_index, image in enumerate([fits_files[fits_index]]):
                ax = f.add_subplot(shape_file, shape_file, image_index + 1, projection=w)
                create_multi_overlap_ax_cutout(ax, fits_names[fits_index], image,
                                               catalog_coordinate=roberto_ra_dec[roberto_ra_dec_index],
                                               matches=values, ra_dec=coords, rob_z=0)
                ax.set_title(aspecs_lines[values]['name'])
        image_index += 1
    # plt.show()
    f.savefig(str("Continuum/ASPECS_Continuum_{}.png".format(fits_names[fits_index])), dpi=300)
    f.clf()
    plt.close()
'''
#exit()
exit()
# Now have the matches, plot them on the sky

all_restframe_ghz = {}

# Since all of them have a single match, just check if each one has a value > 0 and go with those
for fits_index in range(len(fits_files)):
    f = plt.figure(figsize=(20, 20))
    f.suptitle(
        'CO Lines Matched To Galaxies')
    image_index = 0
    for key, values in back_match.items():
        if len(values) > 0:
            # Make the cutouts
            shape_file = int(np.ceil(np.sqrt(len(fits_files))))
            test_mask = (roberto_muse['id'] == key)
            roberto_ra_dec_index = 1e30
            for index, i in enumerate(roberto_muse):
                if i['id'] == key:
                    roberto_ra_dec_index = index
            for third_index, image in enumerate([fits_files[fits_index]]):
                ax = f.add_subplot(shape_file, shape_file, image_index + 1, projection=w)
                create_multi_overlap_ax_cutout(ax, fits_names[fits_index], image,
                                               catalog_coordinate=roberto_ra_dec[roberto_ra_dec_index],
                                               matches=values, ra_dec=coords, rob_z=rob_z[values[0]])
        image_index += 1
    # plt.show()
    f.savefig(str("May_Output/ASPECS_Cutout_Matched_Galaxies_Sep1.5_SN6.0_{}.png".format(fits_names[fits_index])), dpi=300)
    f.clf()
    plt.close()


all_restframe_ghz = {}

for key, values in back_match.items():
    if len(values) > 0:
        # Make the cutouts
        shape_file = int(np.ceil(np.sqrt(len(fits_files))))
        f = plt.figure(figsize=(20, 20))
        test_mask = (roberto_muse['id'] == key)
        roberto_ra_dec_index = 1e30
        for index, i in enumerate(roberto_muse):
            if i['id'] == key:
                roberto_ra_dec_index = index

        distances = []
        freq_valus = []
        rest_frame_ghz = []
        for value in values:
            print("Value: ", value)
            print("Freq: ", freqs[value])
            freq_valus.append(freqs[value])
            print("Rest Frame GHz: " + str(convert_to_rest_frame_ghz(roberto_muse[test_mask]['z_1'][0], freqs[value])))
            rest_frame_ghz.append(
                np.round(convert_to_rest_frame_ghz(np.round(roberto_muse[test_mask]['z_1'][0], 3), freqs[value]), 3))
            for index, id in enumerate(idxc):
                if index == value:
                    distances.append(np.round(coords[index].separation(roberto_ra_dec[id]).arcsecond, 4))
                    all_restframe_ghz[value] = (np.round(roberto_muse[test_mask]['z_1'][0], 3), np.round(
                        convert_to_rest_frame_ghz(np.round(roberto_muse[test_mask]['z_1'][0], 3), freqs[value]), 3),
                                                key)
        f.suptitle(
            'Roberto ID: ' + str(key) + " Z_1: " + str(roberto_muse[test_mask]['z_1'][0]) + " Z_2: " + str(
                roberto_muse[test_mask]['z_2'][0]) +
            " Matches: " + str(freq_valus) + " \nDistance: " + str(distances) + "\n Spec Z: " + str(
                z_specs[roberto_muse[test_mask]['id'][0]]) + ""
            + "\n Rest Frame GHz: " + str(rest_frame_ghz))
        for third_index, image in enumerate(fits_files):
            ax = f.add_subplot(shape_file, shape_file, third_index + 1, projection=w)
            create_multi_overlap_ax_cutout(ax, fits_names[third_index], image,
                                           catalog_coordinate=roberto_ra_dec[roberto_ra_dec_index],
                                           matches=values, index=idx, ra_dec=coords)
        # plt.show()
        f.savefig(str("May_Output/ASPECS_Cutout_" + str(key) + ".png"), dpi=300)
        f.clf()
        plt.close()


for key, values in enumerate(aspecs_matches):
    if len(values) > 0:
        # Make the cutouts
        shape_file = int(np.ceil(np.sqrt(len(fits_files))))
        f = plt.figure(figsize=(20, 20))
        test_mask = (roberto_muse['id'] == key)
        roberto_ra_dec_index = 1e30
        for index, i in enumerate(roberto_muse):
            if i['id'] == key:
                roberto_ra_dec_index = index

        distances = []
        freq_valus = []
        rest_frame_ghz = []
        for value in values:
            print("Value: ", value)
            print("Freq: ", freqs[value])
            freq_valus.append(freqs[value])
            print("Rest Frame GHz: " + str(convert_to_rest_frame_ghz(roberto_muse[test_mask]['z_1'][0], freqs[value])))
            rest_frame_ghz.append(
                np.round(convert_to_rest_frame_ghz(np.round(roberto_muse[test_mask]['z_1'][0], 3), freqs[value]), 3))
            for index, id in enumerate(idxc):
                if index == value:
                    distances.append(np.round(coords[index].separation(roberto_ra_dec[id]).arcsecond, 4))
                    all_restframe_ghz[value] = (np.round(roberto_muse[test_mask]['z_1'][0], 3), np.round(
                        convert_to_rest_frame_ghz(np.round(roberto_muse[test_mask]['z_1'][0], 3), freqs[value]), 3),
                                                key)
        f.suptitle(
            'Roberto ID: ' + str(key) + " Z_1: " + str(roberto_muse[test_mask]['z_1'][0]) + " Z_2: " + str(
                roberto_muse[test_mask]['z_2'][0]) +
            " Matches: " + str(freq_valus) + " \nDistance: " + str(distances) + "\n Spec Z: " + str(
                z_specs[roberto_muse[test_mask]['id'][0]])
            + "\n Rest Frame GHz: " + str(rest_frame_ghz))
        f.suptitle(
            " Z: " + str(row['Z (CO)']) + " Delta_Z: " + str(
                row['Delta Z']) + " S/N: " + str(row['S/N']) +
            " Observed Freq: " + str(freq_valus) + "\n Spec Z: " + str(
                row['Spec Z'] + " Transition: " + str(row['Transition'])) +
            " RA: " + str(row['RA (J2000)']) + " DEC: " + str(row['DEC (J2000)']) +
            "\n Rest Frame GHz: " + str(rest_frame_ghz) + " MStar: " + str(row['Mstar']) + " SFR: " + str(row['SFR']))
        for third_index, image in enumerate(fits_files):
            ax = f.add_subplot(shape_file, shape_file, third_index + 1, projection=w)
            create_multi_overlap_ax_cutout(ax, fits_names[third_index], image,
                                           catalog_coordinate=roberto_ra_dec[roberto_ra_dec_index],
                                           matches=values, ra_dec=coords)
        # plt.show()
        f.savefig(str("May_Output/ASPECS_Cutout_" + str(key) + ".png"), dpi=300)
        f.clf()
        plt.close()

exit()
# Now work backwards
print(all_restframe_ghz)
print(aspecs_matches)

results_dict = {1: (1.094, 218.71 * u.GHz, 23419),
                2: (1.675, 279.602 * u.GHz, 21510),
                7: (0.085, 101.745 * u.GHz, 14630),
                8: (0.872, 177.885 * u.GHz, 18282),
                9: (1.685, 249.474 * u.GHz, 17508),
                10: (0.925, 181.568 * u.GHz, 9258),
                20: (0.96, 209.308 * u.GHz, 22705),
                21: (0.96, 209.338 * u.GHz, 22705),
                24: (3.11, 377.384 * u.GHz, 13131),
                32: (0.163, 110.94 * u.GHz, 14957),
                49: (1.191, 234.301 * u.GHz, 28901),
                50: (1.535, 235.043 * u.GHz, 17170),
                51: (1.535, 235.063 * u.GHz, 17170),
                52: (1.019, 211.949 * u.GHz, 28851),
                53: (1.537, 241.314 * u.GHz, 52578),
                54: (2.835, 364.747 * u.GHz, 58258),
                55: (3.496, 420.165 * u.GHz, 52994),
                57: (1.597, 243.105 * u.GHz, 58246),
                58: (1.597, 243.084 * u.GHz, 58246),
                59: (1.597, 243.084 * u.GHz, 58246),
                60: (1.597, 243.126 * u.GHz, 58246),
                61: (1.597, 243.126 * u.GHz, 58246),
                62: (1.599, 246.073 * u.GHz, 62572),
                63: (1.599, 246.053 * u.GHz, 62572),
                69: (0.738, 184.459 * u.GHz, 25004),
                79: (1.269, 211.816 * u.GHz, 52367),
                83: (1.599, 246.359 * u.GHz, 57822),
                86: (0.355, 143.493 * u.GHz, 22160),
                89: (2.32, 311.433 * u.GHz, 18553),
                90: (2.32, 311.406 * u.GHz, 18553),
                91: (1.665, 246.118 * u.GHz, 19491)}

"""
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
plt.savefig("92_Line_Can_Locations.png", dpi=300)
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
    f.savefig(str("FINAL_Full_ASPECS_Cutout_" + str(third_index) + "_Freq_" + str(
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
    f.savefig(str("FINAL_Roberto_ASPECS_Cutout_" + str(third_index) + "_Freq_" + str(
        aspecs_freqs[third_index]) + "_Sep_{:0.3f}" + ".png").format(
        coord.separation(ra_dec[idx[third_index]]).arcsecond), dpi=300)
    f.clf()
    plt.close()
