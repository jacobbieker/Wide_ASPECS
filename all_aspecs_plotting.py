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
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import vstack
from astropy.table import Table, hstack
import astropy.units as u
from astropy import constants as const

from scipy import stats
from numpy.polynomial.polynomial import polyfit
from aspecs_catalog_builder import compare_catalog_locations

transitions = {"1-0": [0.0030, 0.3694, 115.271, 0.2801],
               "2-1": [1.0059, 1.7387, 230.538, 1.4277],
               "3-2": [2.0088, 3.1080, 345.796, 2.6129],
               "4-3": [3.0115, 4.4771, 461.041, 3.8030],
               "5-4": [4.0142, 5.8460, 576.268, 4.9933],
               "6-5": [5.0166, 7.2146, 691.473, 6.1843],
               "7-6": [6.0188, 8.5829, 806.652, 7.3750],
               "C1 1-0": [3.2823, 4.8468, 492.161, 4.1242],
               "C1 2-1": [6.0422, 8.6148, 809.342, 7.4031]}


def get_observed_ghz(z, transition):
    """
    Get the observed GHz for given redshift based on transition, only for C3-2 and above, from Wide ASPECS paper
    :param z: Z to calculate for
    :param transition:
    :return:
    """
    nm = (transitions[transition][2] * u.GHz).to(u.nm, equivalencies=u.spectral())
    observed_nm = nm * (z + 1)

    observed_ghz = observed_nm.to(u.GHz, equivalencies=u.spectral())
    return observed_ghz


def get_kms(channels, observed_ghz):
    channel_width = 3.9025 * u.MHz
    channel_width = channel_width.to(u.GHz)
    observed_ghz = observed_ghz * u.GHz
    width = channels * channel_width

    print(width)

    print(width * 299792.458)
    fwhm = (width * 299792.458) / observed_ghz

    fwhm = fwhm * (u.km/u.s)

    return fwhm

def get_estimated_z(ghz):
    """
    Estimate the CO line based on Wide-ASPECS one, (3-2), z > 2 or higher J, calculate possible Z's and find which Z is closest
    to the <z> value from Wide ASPECS
    :param ghz:
    :return: transition, estimated_z
    """
    differences = []
    for key, values in transitions.items():
        if key != "1-0" and key != "2-1":
            observed_ghz = get_observed_ghz(values[3], key)
            differences.append((key, np.abs(ghz - observed_ghz)))

    min_diff = np.min([i[1] for i in differences])

    for index, element in enumerate(differences):
        if np.isclose(element[1], min_diff):
            return element[0], transitions[element[0]][3]


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
    # print("Nm Emitted: {}, Z: {}".format(nm_emitted, z))

    # third step is to convert from restframe nm back to GHz

    final_ghz = (nm_emitted).to(u.GHz, equivalencies=u.spectral())

    return final_ghz


def get_delta_z(z, rest_ghz, ghz):
    """
    Take a measured GHz value, and calculates the restframe GHz value based on the given z of the matched galaxy
    :param z:
    :param ghz:
    :return:
    """

    # First step is to convert to nm rom rest frame GHz
    for key, values in transitions.items():
        if values[0] < z < values[1]:
            print("GHz: ", ghz)
            line_to_observed_ghz = get_observed_ghz(z, key)
            print("CO Line Observed: ", line_to_observed_ghz)
            print("Diff: ", (ghz - line_to_observed_ghz))

            temp_observed = line_to_observed_ghz.to(u.nm, equivalencies=u.spectral())
            temp_measured = ghz.to(u.nm, equivalencies=u.spectral())

            print("Diff Measured: ", ((temp_measured - temp_observed) / temp_observed))

            print("Rest GHz: ", rest_ghz)
            sghz = values[2] * u.GHz
            print("GHz: ", ghz)
            sghz = sghz.to(u.nm, equivalencies=u.spectral())
            rest_ghz = rest_ghz.to(u.nm, equivalencies=u.spectral())
            print("Diff: ", ((rest_ghz - sghz) / sghz))
            print("")
            set_z = ((rest_ghz - sghz) / sghz)
    return set_z

aspecs_lines = Table.read("data/line_search_P3_wa_crop.out", format="ascii", header_start=0, data_start=1)

print(aspecs_lines)

# hdu_list = fits.open(os.path.join("data", "jacob_mapghys_in_nov2018_all_jcb4_magphys_jcb4.fits"))
# initial_catalog = Table.read(os.path.join("data", "jacob_mapghys_in_nov2018_all_jcb4_magphys_jcb4.fits"), format='fits')#hdu_list[1].data
initial_catalog = Table.read("magphys_catalog.fits", format='fits')  # hdu_list[1].data

roberto_catalog = Table.read("roberto_catalog_muse_skelton.fits", format='fits')
# initial_catalog['ra'] = np.zeros(len(initial_catalog['z']))
# initial_catalog['dec'] = np.zeros(len(initial_catalog['z']))

# Now have that, match by skelton_id, then if id not zero, match by Sky
# for index, row in enumerate(initial_catalog):
#    skelton_id = row['id']
#    muse_mask = (np.isclose(roberto_catalog['id'], skelton_id))
#    if len(roberto_catalog[muse_mask]['ra']) == 1:
#        initial_catalog[index]['ra'] = roberto_catalog[muse_mask]['ra']
#        initial_catalog[index]['dec'] = roberto_catalog[muse_mask]['dc']
#    else:
#        print(skelton_id)
#        print(roberto_catalog[muse_mask]['ra'])


# initial_catalog.write("magphys_catalog.fits", format='fits')

# Transitions


ra_dec = SkyCoord(initial_catalog['ra'] * u.deg, initial_catalog['dec'] * u.deg, frame='fk5')

coords = SkyCoord(aspecs_lines['rra'] * u.deg, aspecs_lines['rdc'] * u.deg, frame='fk5')

idx, d2d, d3d = coords.match_to_catalog_sky(ra_dec)
print("\n----------------- Number of Matches: " + str(len(idx)) + "/" + str(len(initial_catalog[idx])))
print("Distances: ")
num_in_close = 0
catalog_ids = []
aspecs_redshifts = []
for index, id in enumerate(idx):
    # print(coords[index].separation(ra_dec[id]).arcsecond)
    if coords[index].separation(ra_dec[id]).arcsecond < 0.25:
        num_in_close += 1
        rest_ghz = convert_to_rest_frame_ghz(initial_catalog[id]['z'], aspecs_lines[index]['rfreq'])
        for key, values in transitions.items():
            if values[0] < initial_catalog[id]['z'] < values[1]:
                diff_freq = values[2] * u.GHz - rest_ghz
                delta_z = get_delta_z(initial_catalog[id]['z'], rest_ghz, aspecs_lines[index]['rfreq'] * u.GHz)
                if aspecs_lines[index]['rsnrrbin'] > 0.:
                    aspecs_redshifts.append((id, index, rest_ghz, diff_freq, key, aspecs_lines[index]['rsnrrbin'],
                                             delta_z, initial_catalog[id]['z']))
        # print("\nMatch: " + str(index))
        # print("Distance: " + str(coords[index].separation(ra_dec[id]).arcsecond))
        # print("Location (RA): " + str(coords[index].ra.hms))
        # print("Location (Dec): " + str(coords[index].dec.hms))
        # print("Location (Deg): " + str(coords[index]))
        try:
            # print("Catalog ID: " + str(initial_catalog[id]['id']))
            catalog_ids.append((initial_catalog[id]['id'], index))
            print("In Skelton et al. Catalog: " + str(initial_catalog[id]['flag_3dh']))
        except:
            continue

print("Number Close to Catalog: ", num_in_close)
print("Catalog IDs")
print(catalog_ids)
print(aspecs_redshifts)

sub_arr = [i[4] for i in aspecs_redshifts]
reds_arr = [i[6] for i in aspecs_redshifts]
match_arr = [i[0] for i in aspecs_redshifts]

ids_matched, counts = np.unique(match_arr, return_counts=True)

for index, value in enumerate(counts):
    if value == 3:
        transition_one = []
        for i, element in enumerate(match_arr):
            if element == ids_matched[index]:
                transition_one.append(aspecs_redshifts[i][4])
                transition_one.append(aspecs_redshifts[i][6])
                transition_one.append(aspecs_redshifts[i][7])
        print(transition_one)

print("Min Delt_Z: ", np.min(reds_arr))
print("Max Delt_Z: ", np.max(reds_arr))
print("Mean Delt_Z: ", np.mean(reds_arr))

print(np.unique(sub_arr, return_counts=True))

catalog_ids_used = [i[1] for i in catalog_ids]

aspecs_catalog = initial_catalog[catalog_ids_used]
full_catalog = initial_catalog

exit()


# Create Graph


def create_points_and_error(column_base_name, full_catalog):
    centerpoints = full_catalog[str(column_base_name + "_50")]
    lower_error = full_catalog[str(column_base_name + "_16")]
    upper_error = full_catalog[str(column_base_name + "_84")]
    centerpoints = np.nan_to_num(centerpoints)
    zero_mask = centerpoints != 0.0
    centerpoints = centerpoints[zero_mask]
    lower_error = centerpoints - lower_error[zero_mask]
    upper_error = upper_error[zero_mask] - centerpoints
    error_bars = [lower_error, upper_error]
    print(error_bars[0].shape)

    return centerpoints, error_bars


def create_points_and_error_by_z(column_base_name, full_catalog, low_z, high_z):
    z_mask = (full_catalog["z"] > low_z) & (full_catalog["z"] < high_z)
    centerpoints = full_catalog[z_mask][str(column_base_name + "_50")]
    lower_error = full_catalog[z_mask][str(column_base_name + "_16")]
    upper_error = full_catalog[z_mask][str(column_base_name + "_84")]
    z_values = full_catalog[z_mask]["z"]
    centerpoints = np.nan_to_num(centerpoints)
    lower_error = np.nan_to_num(lower_error)
    upper_error = np.nan_to_num(upper_error)
    zero_mask = centerpoints != 0.0
    centerpoints = centerpoints[zero_mask]
    lower_error = centerpoints - lower_error[zero_mask]
    upper_error = upper_error[zero_mask] - centerpoints
    z_values = z_values[zero_mask]
    error_bars = [lower_error, upper_error]
    return centerpoints, error_bars, z_values


def perform_cuts(full_catalog):
    """
    Perform quality cuts like expected
    :param full_catalog:
    :return:
    """
    quality_cuts = (full_catalog["Q0"] < 2.0) & (full_catalog["Q2"] < 1.0) & (
            (full_catalog["Mstar_50"] - full_catalog["Mstar_16"]) < 0.5) & \
                   ((full_catalog["Mstar_84"] - full_catalog["Mstar_50"]) < 0.5) & \
                   ((full_catalog["SFR_50"] - full_catalog["SFR_16"]) < 0.5) & \
                   ((full_catalog["SFR_84"] - full_catalog["SFR_16"]) < 0.5)

    return full_catalog[quality_cuts]


aspecs_catalog = perform_cuts(aspecs_catalog)
full_catalog = perform_cuts(full_catalog)

# ASPECS Ones

aspecs_low_z_sfr, low_z_sfr_error, low_z_sfr_z = create_points_and_error_by_z("SFR", aspecs_catalog, 0, 1)
aspecs_mid_z_sfr, mid_z_sfr_error, mid_z_sfr_z = create_points_and_error_by_z("SFR", aspecs_catalog, 1, 2)
aspecs_high_z_sfr, high_z_sfr_error, high_z_sfr_z = create_points_and_error_by_z("SFR", aspecs_catalog, 2, 3)
aspecs_vhigh_z_sfr, vhigh_z_sfr_error, vhigh_z_sfr_z = create_points_and_error_by_z("SFR", aspecs_catalog, 3, 4)

aspecs_low_z_mass, low_z_mass_error, low_z_mass_z = create_points_and_error_by_z("Mstar", aspecs_catalog, 0, 1)
aspecs_mid_z_mass, mid_z_mass_error, mid_z_mass_z = create_points_and_error_by_z("Mstar", aspecs_catalog, 1, 2)
aspecs_high_z_mass, high_z_mass_error, high_z_mass_z = create_points_and_error_by_z("Mstar", aspecs_catalog, 2, 3)
aspecs_vhigh_z_mass, vhigh_z_mass_error, vhigh_z_mass_z = create_points_and_error_by_z("Mstar", aspecs_catalog, 3, 4)

print("Lenghth of ASPECS VHigh: {}".format(len(aspecs_vhigh_z_sfr)))
print("Lenghth of ASPECS High: {}".format(len(aspecs_high_z_sfr)))
print("Lenghth of ASPECS Mid: {}".format(len(aspecs_mid_z_sfr)))
print("Lenghth of ASPECS Low: {}".format(len(aspecs_low_z_sfr)))

low_z_sfr, low_z_sfr_error, low_z_sfr_z = create_points_and_error_by_z("SFR", full_catalog, 0, 1)
mid_z_sfr, mid_z_sfr_error, mid_z_sfr_z = create_points_and_error_by_z("SFR", full_catalog, 1, 2)
high_z_sfr, high_z_sfr_error, high_z_sfr_z = create_points_and_error_by_z("SFR", full_catalog, 2, 3)
vhigh_z_sfr, vhigh_z_sfr_error, vhigh_z_sfr_z = create_points_and_error_by_z("SFR", full_catalog, 3, 4)

low_z_mass, low_z_mass_error, low_z_mass_z = create_points_and_error_by_z("Mstar", full_catalog, 0, 1)
mid_z_mass, mid_z_mass_error, mid_z_mass_z = create_points_and_error_by_z("Mstar", full_catalog, 1, 2)
high_z_mass, high_z_mass_error, high_z_mass_z = create_points_and_error_by_z("Mstar", full_catalog, 2, 3)
vhigh_z_mass, vhigh_z_mass_error, vhigh_z_mass_z = create_points_and_error_by_z("Mstar", full_catalog, 3, 4)

# Just the All

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='all', sharey='all')
ax1.errorbar(low_z_mass, low_z_sfr, yerr=low_z_sfr_error, xerr=low_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='darkgrey', elinewidth=1)
ax1.plot(np.unique(low_z_mass), np.poly1d(np.polyfit(low_z_mass, low_z_sfr, 1))(np.unique(low_z_mass)), label='All fit',
         color='black', zorder=10)
ax1.scatter(aspecs_low_z_mass, aspecs_low_z_sfr, marker='.',
            s=5, c='red', label='ASPECS', zorder=20)
ax1.set_title('0 < Z < 1')
ax2.errorbar(mid_z_mass, mid_z_sfr, yerr=mid_z_sfr_error, xerr=mid_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='darkgrey', elinewidth=1)
ax2.plot(np.unique(mid_z_mass), np.poly1d(np.polyfit(mid_z_mass, mid_z_sfr, 1))(np.unique(mid_z_mass)), color='black',
         zorder=10)
ax2.scatter(aspecs_mid_z_mass, aspecs_mid_z_sfr, marker='.',
            s=5, c='red', zorder=20)
ax2.set_title('1 < Z < 2')
ax3.errorbar(high_z_mass, high_z_sfr, yerr=high_z_sfr_error, xerr=high_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='darkgrey', elinewidth=1)
ax3.plot(np.unique(high_z_mass), np.poly1d(np.polyfit(high_z_mass, high_z_sfr, 1))(np.unique(high_z_mass)),
         color='black',
         zorder=10)
ax3.scatter(aspecs_high_z_mass, aspecs_high_z_sfr, marker='.',
            s=5, c='red', zorder=20)
ax3.set_title('2 < Z < 3')
ax4.errorbar(vhigh_z_mass, vhigh_z_sfr, yerr=vhigh_z_sfr_error, xerr=vhigh_z_mass_error, ecolor='lightgrey', fmt='.',
             ms=1, mec='darkgrey', elinewidth=1)
ax4.plot(np.unique(vhigh_z_mass), np.poly1d(np.polyfit(vhigh_z_mass, vhigh_z_sfr, 1))(np.unique(vhigh_z_mass)),
         color='black', zorder=10)
ax4.scatter(aspecs_vhigh_z_mass, aspecs_vhigh_z_sfr, marker='.',
            s=5, c='red', zorder=20)
# ax4.errorbar(vhigh_z_mass, vhigh_z_sfr, yerr=vhigh_z_sfr_error, xerr=vhigh_z_mass_error)
ax4.set_title('3 < Z < 4')
handles, labels = ax1.get_legend_handles_labels()
f.legend(loc='best', handles=handles, labels=labels, prop={'size': 6})
f.text(0.5, 0.01, 'Log(M*)', ha='center')
f.text(0.01, 0.5, 'Log(SFR)', va='center', rotation='vertical')
f.savefig("AllMStarSFR_ASPECS_MATCHES.png", bbox_inches='tight', dpi=300)
f.show()
exit()
