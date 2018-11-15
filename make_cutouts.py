import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy.nddata import Cutout2D

import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.patches import Circle
from aspecs_catalog_builder import get_aspecs_radec


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
        ax = f.add_subplot(shape_file, shape_file, index+1, projection=wcs_data)
        create_ax_cutout(ax, fits_names[index], image, aspecs_coordinates, catalog_coordinates)

    return f


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
ia427_goodss = fits.open("/home/jacob/Research/GOODS-S_IA427/GOODS-S_IA427_sci.fits")
ia445_goodss = fits.open("/home/jacob/Research/GOODS-S_IA445/GOODS-S_IA445_sci.fits")
ia505_goodss = fits.open("/home/jacob/Research/GOODS-S_IA505/GOODS-S_IA505_sci.fits")
ia527_goodss = fits.open("/home/jacob/Research/GOODS-S_IA527/GOODS-S_IA527_sci.fits")
ia550_goodss = fits.open("/home/jacob/Research/GOODS-S_IA550/GOODS-S_IA550_sci.fits")
ia574_goodss = fits.open("/home/jacob/Research/GOODS-S_IA574/GOODS-S_IA574_sci.fits")
ia624_goodss = fits.open("/home/jacob/Research/GOODS-S_IA624/GOODS-S_IA624_sci.fits")
ia651_goodss = fits.open("/home/jacob/Research/GOODS-S_IA651/GOODS-S_IA651_sci.fits")
ia679_goodss = fits.open("/home/jacob/Research/GOODS-S_IA679/GOODS-S_IA679_sci.fits")
ia738_foodss = fits.open("/home/jacob/Research/GOODS-S_IA738/GOODS-S_IA738_sci.fits")
ia797_foodss = fits.open("/home/jacob/Research/GOODS-S_IA797/GOODS-S_IA797_sci.fits")

# IRAC
irac1_goodss = fits.open("/home/jacob/Research/GOODS-S_SEDS1/GOODS-S_SEDS1_sci_sub.fits")
irac2_goodss = fits.open("/home/jacob/Research/GOODS-S_SEDS2/GOODS-S_SEDS2_sci_sub.fits")
irac3_goodss = fits.open("/home/jacob/Research/GOODS-S_irac3/GOODS-S_irac3_s1_sci.fits")
irac4_goodss = fits.open("/home/jacob/Research/GOODS-S_irac4/GOODS-S_irac4_s1_sci.fits")

# Add to the list with names
fits_files = [f125w_goodss, f140w_goodss, f160w_goodss, f435w_goodss, f606w_goodss, f775w_goodss, f850lp_goodss,
              f814w_goodss, R_goodss, U38_goodss, V_goodss, B_goodss, J_goodss, H_goodss, I_goodss, tKs_goodss,
              tJ_goodss, ia427_goodss, ia445_goodss, ia505_goodss, ia527_goodss, ia550_goodss, ia574_goodss,
              ia624_goodss, ia651_goodss, ia679_goodss, ia738_foodss, ia797_foodss, irac1_goodss, irac2_goodss,
              irac3_goodss, irac4_goodss]

fits_names = ["F125W", "F140W", "F160W", "F435W", "F606W", "F775W", "F850LP",
              "F814W", "R", "U38", "V", "B", "J", "H", "I", "tKs",
              "tJ", "IA427", "IA445", "IA505", "IA527", "IA550", "IA574",
              "IA624", "IA651", "IA679", "IA738", "IA797",
              "IRAC1", "IRAC2", "IRAC3", "IRAC4"]
catalog_goodss = fits.open("/home/jacob/Research/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS")
skelton_goodss = catalog_goodss[1].data
aspecs_ra_dec, aspecs_freqs = get_aspecs_radec()
ra_dec = SkyCoord(skelton_goodss['ra'] * u.deg, skelton_goodss['dec'] * u.deg, frame='fk5')

w = wcs.WCS(f160w_goodss[0].header)
ax = plt.subplot(projection=w)
ax.imshow(f160w_goodss[0].data, cmap='gray')
ax.scatter(aspecs_ra_dec.ra.deg, aspecs_ra_dec.dec.deg, transform=ax.get_transform('fk5'), s=100,
           edgecolor='white', facecolor='none')
plt.show()
idx, d2d, d3d = aspecs_ra_dec.match_to_catalog_sky(ra_dec)

for index, coord in enumerate(aspecs_ra_dec):
    f = create_aspecs_cutouts(coord, fits_files, fits_names, wcs_data=w, catalog_coordinates=ra_dec[idx[index]],
                              id=index, aspecs_freqs=aspecs_freqs)
    f.savefig("Skelton_ASPECS_Cutout_" + str(index) + "_Freq_" + str(aspecs_freqs[index]) + ".png", dpi=300)
    f.clf()
    plt.close()

# Now do it with Full Catalog
catalog_goodss = fits.open("data/jacob_aspecs_catalog_fixed_magphys_jcb3.fits")
full_goodss = catalog_goodss[1].data
aspecs_ra_dec, aspecs_freqs = get_aspecs_radec()
ra_dec = SkyCoord(full_goodss['ra'] * u.deg, full_goodss['dec'] * u.deg, frame='fk5')

for index, coord in enumerate(aspecs_ra_dec):
    f = create_aspecs_cutouts(coord, fits_files, fits_names, wcs_data=w, catalog_coordinates=ra_dec[idx[index]],
                              id=index, aspecs_freqs=aspecs_freqs)
    f.savefig("Full_ASPECS_Cutout_" + str(index) + "_Freq_" + str(aspecs_freqs[index]) + ".png", dpi=300)
    f.clf()
    plt.close()


# Now do it with the Roberto's Catalog, getting the Ra and Dec from the full catalog
roberto_catalog = fits.open("data/magphys_in.fits")
spec_ids = full_goodss['id']
diff_ids = []
for id in roberto_catalog['id']:
    if id in spec_ids:
        diff_ids.append(id)
# now create smaller catalog with full catalog info

rows_to_use = []
for index, row in enumerate(full_goodss):
    if row['id'] in diff_ids:
        rows_to_use.append(index)

smaller_catalog = full_goodss[rows_to_use]
aspecs_ra_dec, aspecs_freqs = get_aspecs_radec()
ra_dec = SkyCoord(smaller_catalog['ra'] * u.deg, smaller_catalog['dec'] * u.deg, frame='fk5')

for index, coord in enumerate(aspecs_ra_dec):
    f = create_aspecs_cutouts(coord, fits_files, fits_names, wcs_data=w, catalog_coordinates=ra_dec[idx[index]],
                              id=index, aspecs_freqs=aspecs_freqs)
    f.savefig("Roberto_ASPECS_Cutout_" + str(index) + "_Freq_" + str(aspecs_freqs[index]) + ".png", dpi=300)
    f.clf()
    plt.close()
