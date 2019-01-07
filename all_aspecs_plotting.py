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
from astropy.table import Table, hstack, join
import astropy.units as u
from astropy import constants as const

from scipy import stats
from numpy.polynomial.polynomial import polyfit
from aspecs_catalog_builder import compare_catalog_locations

transitions = {"1-0": [0.0030, 0.3694, 115.271, 0.2801, 89],
               "2-1": [1.0059, 1.7387, 230.538, 1.4277, 1920],
               "3-2": [2.0088, 3.1080, 345.796, 2.6129, 3363],
               "4-3": [3.0115, 4.4771, 461.041, 3.8030, 4149],
               "5-4": [4.0142, 5.8460, 576.268, 4.9933, 4571],
               "6-5": [5.0166, 7.2146, 691.473, 6.1843, 4809],
               "7-6": [6.0188, 8.5829, 806.652, 7.3750, 4935], }

C1_transitions = {"C1 1-0": [3.2823, 4.8468, 492.161, 4.1242, 4287],
                  "C1 2-1": [6.0422, 8.6148, 809.342, 7.4031, 4936]}


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

    #print(width)

    #print(width * 299792.458)
    fwhm = (width * 299792.458) / observed_ghz

    fwhm = fwhm * (u.km / u.s)

    return fwhm

kms = get_kms(5, 115.271)
kms /= (70 * (u.km / u.s / u.Mpc))
print((kms * (1 + 0.2801)**3)**3)


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
            print(key)
            delta_z, matched_key = get_delta_z(values[3], ghz * u.GHz)
            print(delta_z)
            delta_z1, matched_key1 = get_delta_z(values[0], ghz * u.GHz)
            print(delta_z1)
            delta_z2, matched_key2 = get_delta_z(values[1], ghz * u.GHz)
            print(delta_z2)

            differences.append((matched_key, delta_z))
            differences.append((matched_key1, delta_z1))
            differences.append((matched_key2, delta_z2))


    min_diff = np.min([np.abs(i[1]) for i in differences])
    print("Differences: ", differences)

    for index, element in enumerate(differences):
        if np.isclose(np.abs(element[1]), min_diff):
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


def get_delta_z(z, rest_ghz, ghz=None):
    """
    Take a measured GHz value, and calculates the restframe GHz value based on the given z of the matched galaxy
    :param z:
    :param ghz:
    :return:
    """

    # First step is to convert to nm rom rest frame GHz
    set_zs = []
    for key, values in transitions.items():
        if values[0] <= z <= values[1]:
            #print("GHz: ", ghz)
            line_to_observed_ghz = get_observed_ghz(z, key)
            #print("CO Line Observed: ", line_to_observed_ghz)
            #print("Diff: ", (ghz - line_to_observed_ghz))

            #print("Rest GHz: ", rest_ghz)
            sghz = values[2] * u.GHz
            #print("GHz: ", ghz)
            sghz = sghz.to(u.nm, equivalencies=u.spectral())
            rest_ghz = rest_ghz.to(u.nm, equivalencies=u.spectral())
            #print("Diff: ", ((rest_ghz - sghz) / sghz))
            #print("")
            set_z = np.round(((rest_ghz - sghz) / sghz), 3)
            set_zs.append((key, set_z))
    set_z = np.min([np.abs(i[1]) for i in set_zs])
    print(set_zs)
    print(set_z)
    for element in set_zs:
        if np.isclose(np.abs(element[1]),set_z):
            return element[1], element[0]


aspecs_lines = Table.read("data/line_search_P3_wa_crop.out", format="ascii", header_start=0, data_start=1)

#print(aspecs_lines)

# hdu_list = fits.open(os.path.join("data", "jacob_mapghys_in_nov2018_all_jcb4_magphys_jcb4.fits"))
# initial_catalog = Table.read(os.path.join("data", "jacob_mapghys_in_nov2018_all_jcb4_magphys_jcb4.fits"), format='fits')#hdu_list[1].data
initial_catalog = Table.read("magphys_catalog.fits", format='fits')  # hdu_list[1].data

roberto_catalog = Table.read("roberto_catalog_muse_skelton_matched_manFix.fits", format='fits')

initial_catalog = join(initial_catalog, roberto_catalog, keys='id')
#print(initial_catalog)
spec_z_mask = (initial_catalog["z_spec_3dh"] > 0.0001) | (initial_catalog["zm_vds"] > 0.0001) | (
        initial_catalog["zm_coeS"] > 0.0001) \
              | (initial_catalog["zs_mor"] > 0.0001) | (initial_catalog["zm_ina"] > 0.0001) | (initial_catalog["zm_her"] > 0.0001) \
              | (initial_catalog['muse_wide_z'] > 0.0001)

def create_spec_z_mask(spec_catalog):
    #print(spec_catalog['z_spec_3dh'])
    spec_z_mask = (spec_catalog["z_spec_3dh"] > 0.0001) | (spec_catalog["zm_vds"] > 0.0001) | (
            spec_catalog["zm_coeS"] > 0.0001) \
                  | (spec_catalog["zs_mor"] > 0.0001) | (spec_catalog["zm_ina"] > 0.0001) | (spec_catalog["zm_her"] > 0.0001) \
                  | (spec_catalog['muse_wide_z'] > 0.0001)
    #print(spec_z_mask)
    return spec_z_mask

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


ra_dec = SkyCoord(initial_catalog['ra_1'] * u.deg, initial_catalog['dec'] * u.deg, frame='fk5')

coords = SkyCoord(aspecs_lines['rra'] * u.deg, aspecs_lines['rdc'] * u.deg, frame='fk5')

idx, d2d, d3d = coords.match_to_catalog_sky(ra_dec)
print("\n----------------- Number of Matches: " + str(len(idx)) + "/" + str(len(initial_catalog[idx])))
print("Distances: ")
num_in_close = 0
num_out = 0
snr_limit = 5.
catalog_ids = []
aspecs_redshifts = []
aspecs_no_fit = []
aspecs_table = Table(names=('RA (J2000)', 'DEC (J2000)', 'Roberto ID', 'Observed CO (GHz)', 'Restframe CO (GHz)', 'Transition', 'Z',
                            'Spec Z', 'Delta Z', 'Km/s', 'Separation (Arcsecond)', 'S/N', 'Flux Density at Peak (Jy/beam)',
                            'Integrated Flux (Jy km/s)', 'Width (Channels)', 'Cosmic Volume (Mpc^3)'),
                     dtype=('f8', 'f8', 'int32', 'f16', 'f16', 'U6', 'f16', 'bool', 'f16', 'f16', 'f16', 'f16', 'f16', 'f16', 'int16', 'f16'))

for index, id in enumerate(idx):
    # print(coords[index].separation(ra_dec[id]).arcsecond)
    if coords[index].separation(ra_dec[id]).arcsecond < 1.0:
        num_in_close += 1
        rest_ghz = convert_to_rest_frame_ghz(initial_catalog[id]['z_1'], aspecs_lines[index]['rfreq'])
        for key, values in transitions.items():
            if values[0] < initial_catalog[id]['z_1'] < values[1]:
                diff_freq = values[2] * u.GHz - rest_ghz
                delta_z, matched_key = get_delta_z(initial_catalog[id]['z_1'], rest_ghz, aspecs_lines[index]['rfreq'] * u.GHz)
                kms = get_kms(aspecs_lines[index]['width'], aspecs_lines[index]['rfreq'])
                has_spec_z = create_spec_z_mask(initial_catalog[id])
                if aspecs_lines[index]['rsnrrbin'] > snr_limit:
                    aspecs_redshifts.append((initial_catalog[id]['id'], index, rest_ghz, diff_freq, matched_key, aspecs_lines[index]['rsnrrbin'],
                                             delta_z, initial_catalog[id]['z_1'], kms, aspecs_lines[index]['rfreq']))
                    aspecs_table.add_row((np.round(aspecs_lines[index]['rra'], 6), np.round(aspecs_lines[index]['rdc'], 6),
                                          np.int(initial_catalog[id]['id']), aspecs_lines[index]['rfreq'],
                                          rest_ghz, matched_key, initial_catalog[id]['z_1'], has_spec_z, delta_z, kms,
                                          np.round(coords[index].separation(ra_dec[id]).arcsecond, 4), aspecs_lines[index]['rsnrrbin'],
                                          aspecs_lines[index]['rpeak'], aspecs_lines[index]['rflux'], aspecs_lines[index]['width'], transitions[matched_key][4]))
        # print("\nMatch: " + str(index))
        # print("Distance: " + str(coords[index].separation(ra_dec[id]).arcsecond))
        # print("Location (RA): " + str(coords[index].ra.hms))
        # print("Location (Dec): " + str(coords[index].dec.hms))
        # print("Location (Deg): " + str(coords[index]))
        try:
            # print("Catalog ID: " + str(initial_catalog[id]['id']))
            catalog_ids.append((initial_catalog[id]['id_1'], index))
            print("In Skelton et al. Catalog: " + str(initial_catalog[id]['flag_3dh']))
        except:
            continue
    else:
        # Outside the limit, so say no prior, make estimated guess
        num_out += 1
        estimated_transition, estimated_z = get_estimated_z(aspecs_lines[index]['rfreq'])
        kms = get_kms(aspecs_lines[index]['width'], aspecs_lines[index]['rfreq'])
        rest_ghz = convert_to_rest_frame_ghz(estimated_z, aspecs_lines[index]['rfreq'])
        delta_z, matched_key = get_delta_z(estimated_z, rest_ghz, aspecs_lines[index]['rfreq'] * u.GHz)
        has_spec_z = False

        if aspecs_lines[index]['rsnrrbin'] > snr_limit:
            aspecs_no_fit.append((initial_catalog[id]['id'], index, rest_ghz, 0, matched_key, aspecs_lines[index]['rsnrrbin'],
                                     delta_z, initial_catalog[id]['z_1'], kms, aspecs_lines[index]['rfreq']))
            aspecs_table.add_row((np.round(aspecs_lines[index]['rra'], 6), np.round(aspecs_lines[index]['rdc'], 6),
                                  -999, aspecs_lines[index]['rfreq'],
                                  rest_ghz, matched_key, estimated_z+delta_z, has_spec_z, delta_z, kms,
                                  -999, aspecs_lines[index]['rsnrrbin'],
                                  aspecs_lines[index]['rpeak'], aspecs_lines[index]['rflux'], aspecs_lines[index]['width'], transitions[matched_key][4]))

print("Number Close to Catalog: ", num_in_close)
print("Catalog IDs")
print(catalog_ids)
print(aspecs_redshifts)
print(aspecs_table)

sub_arr = [i[4] for i in aspecs_redshifts]
reds_arr = [i[6] for i in aspecs_redshifts]
match_arr = [i[0] for i in aspecs_redshifts]
no_match_arr = [i[4] for i in aspecs_no_fit]

ids_matched, counts = np.unique(match_arr, return_counts=True)
ghz_ids_matched = np.unique(no_match_arr, return_counts=True)
print("GHZ Id No Match: ", ghz_ids_matched)

for index, value in enumerate(counts):
    if value == 3:
        transition_one = [[],[]]
        for i, element in enumerate(match_arr):
            if element == ids_matched[index]:
                transition_one[0].append(aspecs_redshifts[i][4])
                transition_one[1].append(aspecs_redshifts[i][2])
                transition_one[1].append(aspecs_redshifts[i][9])
                transition_one[1].append(aspecs_redshifts[i][6])
                transition_one[1].append(aspecs_redshifts[i][7])
        print(transition_one)

print("Min Delt_Z: ", np.min(reds_arr))
print("Max Delt_Z: ", np.max(reds_arr))
print("Mean Delt_Z: ", np.mean(reds_arr))

print(np.unique(sub_arr, return_counts=True))

catalog_ids_used = [i[1] for i in catalog_ids]

aspecs_catalog = initial_catalog[catalog_ids_used]

# Now a Prior vs. Non-Prior detection Z values
prior_z = []
no_prior = []
for row in aspecs_table:
    if row['Roberto ID'] > -1:
        prior_z.append(row['Z'])
    else:
        no_prior.append(row['Z'])

spec_prior_z = []
spec_no_prior = []
for row in aspecs_table:
    if row['Roberto ID'] > -1:
        if row['Spec Z']:
            spec_prior_z.append(row['Z'])
        else:
            spec_no_prior.append(row['Z'])

x = (prior_z, no_prior)

# Now plot the stacked one
transitions = {"1-0": [0.0030, 0.3694, 115.271, 0.2801, 89],
               "2-1": [1.0059, 1.7387, 230.538, 1.4277, 1920],
               "3-2": [2.0088, 3.1080, 345.796, 2.6129, 3363],
               "4-3": [3.0115, 4.4771, 461.041, 3.8030, 4149],
               "5-4": [4.0142, 5.8460, 576.268, 4.9933, 4571],
               "6-5": [5.0166, 7.2146, 691.473, 6.1843, 4809],
               "7-6": [6.0188, 8.5829, 806.652, 7.3750, 4935], }

fig, ax = plt.subplots()

ax.axvspan(0.0030, 0.3694, facecolor='grey', alpha=0.5)
ax.axvspan(1.0059, 1.7387, facecolor='grey', alpha=0.5)
ax.axvspan(2.0088, 3.1080, facecolor='grey', alpha=0.5)
ax.axvspan(3.0115, 4.4771, facecolor='grey', alpha=0.5)
ax.axvspan(4.0142, 5.8460, facecolor='grey', alpha=0.5)
ax.axvspan(5.0166, 7.2146, facecolor='grey', alpha=0.5)
ax.axvspan(6.0188, 7.3750, facecolor='grey', alpha=0.5)

ax.text((0.0030 + 0.3694)/2, 750, "1-0")
ax.text((1.0059 + 1.7387)/2, 750, "2-1")
ax.text((2.0088 + 3.1080)/2, 750, "3-2")
ax.text((3.0115 + 4.4771)/2, 750, "4-3")
ax.text((4.0142 + 5.8460)/2, 750, "5-4")
ax.text((5.0166 + 7.2146)/2, 750, "6-5")
ax.text((6.0188 + 7.3750)/2, 750, "7-6")

ax.hist(x, 50, histtype='bar', label=('Roberto Prior', 'ASPECS No Prior'), color=('blue', 'green'), stacked=True, log=True)


ax.set_xlabel("Redshift (Z)")
ax.set_ylabel("Number of Detected Sources")
ax.set_title("Detected Sources S/N > " + str(np.round(snr_limit, 2)))
fig.legend(loc='best')
fig.savefig("AllCOLines_SN" + str(np.round(snr_limit, 2)) + ".png", bbox_inches='tight', dpi=300)
fig.show()

fig, ax = plt.subplots()

x = (spec_prior_z, spec_no_prior)
ax.axvspan(0.0030, 0.3694, facecolor='grey', alpha=0.5)
ax.axvspan(1.0059, 1.7387, facecolor='grey', alpha=0.5)
ax.axvspan(2.0088, 3.1080, facecolor='grey', alpha=0.5)
ax.axvspan(3.0115, 4.4771, facecolor='grey', alpha=0.5)
ax.axvspan(4.0142, 5.8460, facecolor='grey', alpha=0.5)
ax.axvspan(5.0166, 7.2146, facecolor='grey', alpha=0.5)
ax.axvspan(6.0188, 7.3750, facecolor='grey', alpha=0.5)

ax.text((0.0030 + 0.3694)/2, 100, "1-0")
ax.text((1.0059 + 1.7387)/2, 100, "2-1")
ax.text((2.0088 + 3.1080)/2, 100, "3-2")
ax.text((3.0115 + 4.4771)/2, 100, "4-3")
ax.text((4.0142 + 5.8460)/2, 100, "5-4")
ax.text((5.0166 + 7.2146)/2, 100, "6-5")
ax.text((6.0188 + 7.3750)/2, 100, "7-6")
ax.hist(x, 50, histtype='bar', label=('Spec Prior', 'Photo Prior'), color=('blue', 'green'), stacked=True, log=True)


ax.set_xlabel("Redshift (Z)")
ax.set_ylabel("Number of Detected Sources")
ax.set_title("Detected Sources S/N > " + str(np.round(snr_limit, 2)))
fig.legend(loc='best')
fig.savefig("OnlyMatchedCOLines_SN" + str(np.round(snr_limit, 2)) + ".png", bbox_inches='tight', dpi=300)
fig.show()

#exit()

# Create Graph


def create_points_and_error(column_base_name, initial_catalog):
    centerpoints = initial_catalog[str(column_base_name + "_50_1")]
    lower_error = initial_catalog[str(column_base_name + "_16_1")]
    upper_error = initial_catalog[str(column_base_name + "_84_1")]
    centerpoints = np.nan_to_num(centerpoints)
    zero_mask = centerpoints != 0.0
    centerpoints = centerpoints[zero_mask]
    lower_error = centerpoints - lower_error[zero_mask]
    upper_error = upper_error[zero_mask] - centerpoints
    error_bars = [lower_error, upper_error]
    print(error_bars[0].shape)

    return centerpoints, error_bars


def create_points_and_error_by_z(column_base_name, initial_catalog, low_z, high_z):
    z_mask = (initial_catalog["z_1"] > low_z) & (initial_catalog["z_1"] < high_z)
    centerpoints = initial_catalog[z_mask][str(column_base_name + "_50_1")]
    lower_error = initial_catalog[z_mask][str(column_base_name + "_16_1")]
    upper_error = initial_catalog[z_mask][str(column_base_name + "_84_1")]
    z_values = initial_catalog[z_mask]["z_1"]
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


def perform_cuts(initial_catalog):
    """
    Perform quality cuts like expected
    :param initial_catalog:
    :return:
    """
    quality_cuts = (initial_catalog["Q0_1"] < 2.0) & (initial_catalog["Q2_1"] < 1.0) & (
            (initial_catalog["Mstar_50_1"] - initial_catalog["Mstar_16_1"]) < 0.5) & \
                   ((initial_catalog["Mstar_84_1"] - initial_catalog["Mstar_50_1"]) < 0.5) & \
                   ((initial_catalog["SFR_50_1"] - initial_catalog["SFR_16_1"]) < 0.5) & \
                   ((initial_catalog["SFR_84_1"] - initial_catalog["SFR_16_1"]) < 0.5)

    return initial_catalog[quality_cuts]


aspecs_catalog = perform_cuts(aspecs_catalog)
initial_catalog = perform_cuts(initial_catalog)

spec_z_mask = create_spec_z_mask(initial_catalog)
spec_catalog = initial_catalog[spec_z_mask]

muse_z_mask = (initial_catalog['zm_ina'] > 0.001) | (initial_catalog['zm_her'] > 0.001) | (initial_catalog['muse_wide_z'] > 0.0001)

muse_catalog = initial_catalog[muse_z_mask]
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

low_z_sfr, low_z_sfr_error, low_z_sfr_z = create_points_and_error_by_z("SFR", initial_catalog, 0, 1)
mid_z_sfr, mid_z_sfr_error, mid_z_sfr_z = create_points_and_error_by_z("SFR", initial_catalog, 1, 2)
high_z_sfr, high_z_sfr_error, high_z_sfr_z = create_points_and_error_by_z("SFR", initial_catalog, 2, 3)
vhigh_z_sfr, vhigh_z_sfr_error, vhigh_z_sfr_z = create_points_and_error_by_z("SFR", initial_catalog, 3, 4)

low_z_mass, low_z_mass_error, low_z_mass_z = create_points_and_error_by_z("Mstar", initial_catalog, 0, 1)
mid_z_mass, mid_z_mass_error, mid_z_mass_z = create_points_and_error_by_z("Mstar", initial_catalog, 1, 2)
high_z_mass, high_z_mass_error, high_z_mass_z = create_points_and_error_by_z("Mstar", initial_catalog, 2, 3)
vhigh_z_mass, vhigh_z_mass_error, vhigh_z_mass_z = create_points_and_error_by_z("Mstar", initial_catalog, 3, 4)

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

muse_low_z_sfr, low_z_sfr_error, low_z_sfr_z = create_points_and_error_by_z("SFR", muse_catalog, 0, 1)
muse_mid_z_sfr, mid_z_sfr_error, mid_z_sfr_z = create_points_and_error_by_z("SFR", muse_catalog, 1, 2)
muse_high_z_sfr, high_z_sfr_error, high_z_sfr_z = create_points_and_error_by_z("SFR", muse_catalog, 2, 3)
muse_vhigh_z_sfr, vhigh_z_sfr_error, vhigh_z_sfr_z = create_points_and_error_by_z("SFR", muse_catalog, 3, 4)

muse_low_z_mass, low_z_mass_error, low_z_mass_z = create_points_and_error_by_z("Mstar", muse_catalog, 0, 1)
muse_mid_z_mass, mid_z_mass_error, mid_z_mass_z = create_points_and_error_by_z("Mstar", muse_catalog, 1, 2)
muse_high_z_mass, high_z_mass_error, high_z_mass_z = create_points_and_error_by_z("Mstar", muse_catalog, 2, 3)
muse_vhigh_z_mass, vhigh_z_mass_error, vhigh_z_mass_z = create_points_and_error_by_z("Mstar", muse_catalog, 3, 4)

leinhardt_low_z_sfr, low_z_sfr_error, low_z_sfr_z = create_points_and_error_by_z("SFR", initial_catalog, 0, 1)
leinhardt_mid_z_sfr, mid_z_sfr_error, mid_z_sfr_z = create_points_and_error_by_z("SFR", initial_catalog, 1, 2)
leinhardt_high_z_sfr, high_z_sfr_error, high_z_sfr_z = create_points_and_error_by_z("SFR", initial_catalog, 2, 3)
leinhardt_vhigh_z_sfr, vhigh_z_sfr_error, vhigh_z_sfr_z = create_points_and_error_by_z("SFR", initial_catalog, 3, 4)

leinhardt_low_z_mass, low_z_mass_error, low_z_mass_z = create_points_and_error_by_z("Mstar", initial_catalog, 0, 1)
leinhardt_mid_z_mass, mid_z_mass_error, mid_z_mass_z = create_points_and_error_by_z("Mstar", initial_catalog, 1, 2)
leinhardt_high_z_mass, high_z_mass_error, high_z_mass_z = create_points_and_error_by_z("Mstar", initial_catalog, 2, 3)
leinhardt_vhigh_z_mass, vhigh_z_mass_error, vhigh_z_mass_z = create_points_and_error_by_z("Mstar", initial_catalog, 3, 4)


f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='all')
ax1.hist(low_z_mass, color='lightgrey', histtype='step', bins=10)
ax1.set_yscale('log')
ax1.set_title('0 < Z < 1')
ax1.hist(aspecs_low_z_mass, color='yellow', histtype='step', bins=10)
ax1.hist(muse_low_z_mass, color='green', histtype='step', bins=10)
ax1.hist(leinhardt_low_z_mass, color='orange', histtype='step', bins=10)
ax2.set_title('0 < Z < 1')
ax2.hist(low_z_sfr, color='lightgrey', histtype='step', bins=10)
ax2.hist(aspecs_low_z_sfr, color='yellow', histtype='step', bins=10)
ax2.hist(muse_low_z_sfr, color='green', histtype='step', bins=10)
ax2.hist(leinhardt_low_z_sfr, color='orange', histtype='step', bins=10)
ax3.set_title('1 < Z < 2')
ax3.hist(mid_z_mass, color='lightgrey', histtype='step', bins=10)
ax3.hist(aspecs_mid_z_mass, color='yellow', histtype='step', bins=10)
ax3.hist(muse_mid_z_mass, color='green', histtype='step', bins=10)
ax3.hist(leinhardt_mid_z_mass, color='orange', histtype='step', bins=10)
ax4.set_title('1 < Z < 2')
ax4.hist(mid_z_sfr, color='lightgrey', histtype='step', label='All', bins=10)
ax4.hist(aspecs_mid_z_sfr, color='yellow', histtype='step', label='Spectroscopic', bins=10)
ax4.hist(muse_mid_z_sfr, color='green', histtype='step', label='MUSE', bins=10)
ax4.hist(leinhardt_mid_z_sfr, color='orange', histtype='step', label='Leindert', bins=10)
handles, labels = ax4.get_legend_handles_labels()
f.legend(loc='best', handles=handles, labels=labels, prop={'size': 6})
f.text(0.25, 0.01, 'Log Stellar Mass (Mstar)', ha='center')
f.text(0.75, 0.01, 'Log SFR (Mstar)', ha='center')
f.text(0.01, 0.25, 'Count', va='center', rotation='vertical')
f.text(0.01, 0.75, 'Count', va='center', rotation='vertical')
f.savefig("LeinhardHist02.png", bbox_inches='tight', dpi=300)
f.show()

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='all')
ax1.hist(high_z_mass, color='lightgrey', histtype='step', bins=10)
ax1.set_yscale('log')
ax1.set_title('2 < Z < 3')
ax1.hist(aspecs_high_z_mass, color='yellow', histtype='step', bins=10)
ax1.hist(muse_high_z_mass, color='green', histtype='step', bins=10)
ax1.hist(leinhardt_high_z_mass, color='orange', histtype='step', bins=10)
ax2.set_title('2 < Z < 3')
ax2.hist(high_z_sfr, color='lightgrey', histtype='step', bins=10)
ax2.hist(aspecs_high_z_sfr, color='yellow', histtype='step', bins=10)
ax2.hist(muse_high_z_sfr, color='green', histtype='step', bins=10)
ax2.hist(leinhardt_high_z_sfr, color='orange', histtype='step', bins=10)
ax3.set_title('3 < Z < 4')
ax3.hist(vhigh_z_mass, color='lightgrey', histtype='step', bins=10)
ax3.hist(aspecs_vhigh_z_mass, color='yellow', histtype='step', bins=10)
ax3.hist(muse_vhigh_z_mass, color='green', histtype='step', bins=10)
ax3.hist(leinhardt_vhigh_z_mass, color='orange', histtype='step', bins=10)
ax4.set_title('3 < Z < 4')
ax4.hist(vhigh_z_sfr, color='lightgrey', histtype='step', label='All', bins=10)
ax4.hist(aspecs_vhigh_z_sfr, color='yellow', histtype='step', label='ASPECS', bins=10)
ax4.hist(muse_vhigh_z_sfr, color='green', histtype='step', bins=10, label='MUSE')
ax4.hist(leinhardt_vhigh_z_sfr, color='orange', histtype='step', label='Leindert', bins=10)
handles, labels = ax4.get_legend_handles_labels()
f.legend(loc='best', handles=handles, labels=labels, prop={'size': 6})
f.text(0.25, 0.01, 'Log Stellar Mass (Mstar)', ha='center')
f.text(0.75, 0.01, 'Log SFR (Mstar)', ha='center')
f.text(0.01, 0.25, 'Count', va='center', rotation='vertical')
f.text(0.01, 0.75, 'Count', va='center', rotation='vertical')
f.savefig("LeinhardHist24.png", bbox_inches='tight', dpi=300)
f.show()


exit()
