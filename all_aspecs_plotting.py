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
               "4-3": [3.0115, 4.4771, 461.041, 3.8030, 4149],}

DON_NOT_USE = {               "5-4": [4.0142, 5.8460, 576.268, 4.9933, 4571],
                       "6-5": [5.0166, 7.2146, 691.473, 6.1843, 4809],
                       "7-6": [6.0188, 8.5829, 806.652, 7.3750, 4935],
                       "C1 1-0": [3.2823, 4.8468, 492.161, 4.1242, 4287],
                       "C1 2-1": [6.0422, 8.6148, 809.342, 7.4031, 4936]}
C1_transitions = {"C1 1-0": [3.2823, 4.8468, 492.161, 4.1242, 4287],
                  "C1 2-1": [6.0422, 8.6148, 809.342, 7.4031, 4936]}

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
def get_comoving_volume(start_z, end_z, sq_arcminutes):
    start_vol = cosmo.comoving_volume(start_z)
    end_vol = cosmo.comoving_volume(end_z)
    fraction_sky = sq_arcminutes / 148510800.

    comoving_volume = (end_vol - start_vol) * fraction_sky

    return comoving_volume
print(cosmo.H(0))
print(cosmo.comoving_volume(0.0873))
print(get_comoving_volume(0.0, 0.0873, 1.024))

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

#TODO Check Channel Width: ALTRVAL = 3.021582003656E+06 Hz

def get_kms(channels, observed_ghz):
    channel_width = 3.021582003656E+06 *u.Hz
    #channel_width = 3.9025 * u.MHz
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

delt_z_stuff = []

def get_estimated_z(ghz):
    """
    Estimate the CO line based on Wide-ASPECS one, (3-2), z > 2 or higher J, calculate possible Z's and find which Z is closest
    to the <z> value from Wide ASPECS
    :param ghz:
    :return: transition, estimated_z
    """
    differences = []
    for key, values in transitions.items():
        if True: #key != "1-0" and key != "2-1":
            #print(key)
            # Convert ghz to rest_ghz of the Z value, otherwise always closest to lowest one
            sghz = convert_to_rest_frame_ghz(values[3], ghz)
            delta_z, matched_key = get_delta_z(values[3], sghz)
            #delta_z = values[3] - delta_zget_delta_z
            differences.append((matched_key, delta_z))


    min_diff = np.min([np.abs(i[1]) for i in differences])
    print("Differences: ", differences)

    for index, element in enumerate(differences):
        if np.isclose(np.abs(element[1]), min_diff):
            delt_z_stuff.append(element[1])
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
            sghz = values[2] * u.GHz # Gets the GHz of the CO line
            sghz = sghz.to(u.nm, equivalencies=u.spectral()) # Now Lambda_Emitted
            rest_ghz = rest_ghz.to(u.nm, equivalencies=u.spectral()) # Lambda_Observed in rest frame
            set_z = np.round(((sghz - rest_ghz) / rest_ghz), 3) # (Lambda_Observed - Lambda_Emitted) / Lambda_Emitted
            # (Z+1)Lambda_Emitted = Lambda_Observed
            print("Z: {} Set Z: {}".format(z, set_z))
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
snr_limit = 6.
z_sep = 0.3
catalog_ids = []
aspecs_redshifts = []
aspecs_no_fit = []
aspecs_table = Table(names=('RA (J2000)', 'DEC (J2000)', 'Roberto ID', 'Observed CO (GHz)', 'Restframe CO (GHz)', 'Transition', 'Z',
                            'Spec Z', 'Delta Z', 'Km/s', 'Separation (Arcsecond)', 'S/N', 'Flux Density at Peak (Jy/beam)',
                            'Integrated Flux (Jy km/s)', 'Width (Channels)', 'Cosmic Volume (Mpc^3)'),
                     dtype=('f8', 'f8', 'int32', 'f4', 'f4', 'U6', 'f4', 'bool', 'f4', 'f8', 'f4', 'f4', 'f4', 'f4', 'int8', 'f4'))

all_more_snr_5 = 0

for row in aspecs_lines:
    if row['rsnrrbin'] > 6.:
        all_more_snr_5 += 1

print(all_more_snr_5)
print(len(aspecs_lines))
#exit()

for index, id in enumerate(idx):
    # print(coords[index].separation(ra_dec[id]).arcsecond)
    if coords[index].separation(ra_dec[id]).arcsecond < 1.0:
        num_in_close += 1
        rest_ghz = convert_to_rest_frame_ghz(initial_catalog[id]['z_1'], aspecs_lines[index]['rfreq'])
        matched_z = False
        for key, values in transitions.items():
            if initial_catalog[id]['z_1'] < 0.4 or 1.1 <= initial_catalog[id]['z_1'] <= 1.8 or 2.2 < initial_catalog[id]['z_1'] < 4.4:
                if values[0] < initial_catalog[id]['z_1'] < values[1]:
                    matched_z = True
                    diff_freq = values[2] * u.GHz - rest_ghz
                    delta_z, matched_key = get_delta_z(initial_catalog[id]['z_1'], rest_ghz, aspecs_lines[index]['rfreq'] * u.GHz)
                    if np.abs(delta_z) < z_sep:
                        kms = get_kms(aspecs_lines[index]['width'], aspecs_lines[index]['rfreq'])
                        has_spec_z = create_spec_z_mask(initial_catalog[id])
                        comoving_volume = get_comoving_volume(values[0], values[1], 52.5)
                        if aspecs_lines[index]['rsnrrbin'] > snr_limit:
                            catalog_ids.append((initial_catalog[id]['id'], index))
                            aspecs_redshifts.append((initial_catalog[id]['id'], index, rest_ghz, diff_freq, matched_key, aspecs_lines[index]['rsnrrbin'],
                                                     delta_z, initial_catalog[id]['z_1'], kms, aspecs_lines[index]['rfreq']))
                            aspecs_table.add_row((np.round(aspecs_lines[index]['rra'], 6), np.round(aspecs_lines[index]['rdc'], 6),
                                                  np.int(initial_catalog[id]['id']), aspecs_lines[index]['rfreq'],
                                                  rest_ghz, matched_key, initial_catalog[id]['z_1'], has_spec_z, delta_z, kms,
                                                  np.round(coords[index].separation(ra_dec[id]).arcsecond, 4), aspecs_lines[index]['rsnrrbin'],
                                                  aspecs_lines[index]['rpeak'], aspecs_lines[index]['rflux'], aspecs_lines[index]['width'], np.round(comoving_volume, 3)))
        if not matched_z:
            num_out += 1
            estimated_transition, estimated_z = get_estimated_z(aspecs_lines[index]['rfreq'])
            if estimated_z < 0.4 or 1.1 <= estimated_z <= 1.8 or 2.2 < estimated_z < 4.4:
                kms = get_kms(aspecs_lines[index]['width'], aspecs_lines[index]['rfreq'])
                rest_ghz = convert_to_rest_frame_ghz(estimated_z, aspecs_lines[index]['rfreq'])
                delta_z, matched_key = get_delta_z(estimated_z, rest_ghz, aspecs_lines[index]['rfreq'] * u.GHz)
                if np.abs(delta_z) < z_sep:
                    has_spec_z = False
                    comoving_volume = get_comoving_volume(transitions[estimated_transition][0], transitions[estimated_transition][1], 52.5)

                    if aspecs_lines[index]['rsnrrbin'] > snr_limit:
                        aspecs_no_fit.append((initial_catalog[id]['id'], index, rest_ghz, 0, matched_key, aspecs_lines[index]['rsnrrbin'],
                                              delta_z, initial_catalog[id]['z_1'], kms, aspecs_lines[index]['rfreq']))
                        aspecs_table.add_row((np.round(aspecs_lines[index]['rra'], 6), np.round(aspecs_lines[index]['rdc'], 6),
                                              -999, aspecs_lines[index]['rfreq'],
                                              rest_ghz, matched_key, estimated_z+delta_z, has_spec_z, delta_z, kms,
                                              -999, aspecs_lines[index]['rsnrrbin'],
                                              aspecs_lines[index]['rpeak'], aspecs_lines[index]['rflux'], aspecs_lines[index]['width'], np.round(comoving_volume, 3)))

        # print("\nMatch: " + str(index))
        # print("Distance: " + str(coords[index].separation(ra_dec[id]).arcsecond))
        # print("Location (RA): " + str(coords[index].ra.hms))
        # print("Location (Dec): " + str(coords[index].dec.hms))
        # print("Location (Deg): " + str(coords[index]))
    else:
        # Outside the limit, so say no prior, make estimated guess
        num_out += 1
        estimated_transition, estimated_z = get_estimated_z(aspecs_lines[index]['rfreq'])
        if estimated_z < 0.4 or 1.1 <= estimated_z <= 1.8 or 2.2 < estimated_z < 4.4:
            kms = get_kms(aspecs_lines[index]['width'], aspecs_lines[index]['rfreq'])
            rest_ghz = convert_to_rest_frame_ghz(estimated_z, aspecs_lines[index]['rfreq'])
            delta_z, matched_key = get_delta_z(estimated_z, rest_ghz, aspecs_lines[index]['rfreq'] * u.GHz)
            if np.abs(delta_z) < z_sep:
                has_spec_z = False
                comoving_volume = get_comoving_volume(transitions[estimated_transition][0], transitions[estimated_transition][1], 52.5)

                if aspecs_lines[index]['rsnrrbin'] > snr_limit:
                    aspecs_no_fit.append((initial_catalog[id]['id'], index, rest_ghz, 0, matched_key, aspecs_lines[index]['rsnrrbin'],
                                             delta_z, initial_catalog[id]['z_1'], kms, aspecs_lines[index]['rfreq']))
                    aspecs_table.add_row((np.round(aspecs_lines[index]['rra'], 6), np.round(aspecs_lines[index]['rdc'], 6),
                                          -999, aspecs_lines[index]['rfreq'],
                                          rest_ghz, matched_key, estimated_z+delta_z, has_spec_z, delta_z, kms,
                                          -999, aspecs_lines[index]['rsnrrbin'],
                                          aspecs_lines[index]['rpeak'], aspecs_lines[index]['rflux'], aspecs_lines[index]['width'], np.round(comoving_volume, 3)))

print(aspecs_table['Delta Z'])
print("Number Close to Catalog: ", num_in_close)
print("Catalog IDs")
print(catalog_ids)
print(aspecs_redshifts)
print(aspecs_table)

ascii.write(aspecs_table, "ASPECS_Line_Candidates_ECSV", format='ecsv')

ascii.write(aspecs_table, "ASPECS_Line_Candidates_Z44_Total_Z_Limit.txt", format='fixed_width', bookend=False, delimiter=None, formats={'RA (J2000)': '%2.6f', 'DEC (J2000)': '%2.6f', 'Observed CO (GHz)': '%3.4f', 'Restframe CO (GHz)': '%3.4f', 'Z': '%2.3f',
                                                                            'Delta Z': '%2.3f', 'Km/s': '%4.3f', 'Separation (Arcsecond)': '%2.4f', 'S/N': '%2.3f', 'Flux Density at Peak (Jy/beam)': '%2.4f',
                                                                              'Integrated Flux (Jy km/s)': '%2.4f', 'Cosmic Volume (Mpc^3)': '%8.0f'})
#exit()

def channel_sn_hist(aspecs_table):
    integrated_fluxes = aspecs_table['Integrated Flux (Jy km/s)']
    channels = aspecs_table['Width (Channels)']
    peak_fluxes = aspecs_table["Flux Density at Peak (Jy/beam)"]
    km_s = aspecs_table['Km/s']

    print(integrated_fluxes)

    integrated_vs_channels = integrated_fluxes/ channels
    peak_vs_channels = peak_fluxes / channels
    integrated_vs_kms = integrated_fluxes / km_s

    plt.hist(integrated_vs_channels)
    plt.title("Integrated vs Channels SN>{}".format(snr_limit))
    plt.ylabel("Count")
    plt.xlabel("Integrated Flux / Channels (Jy km/s /channel)")
    plt.show()

    plt.hist(integrated_vs_kms)
    plt.title("Integrated vs Channels SN>{}".format(snr_limit))
    plt.ylabel("Count")
    plt.xlabel("Integrated Flux / Km/s (Jy)")
    plt.show()

    plt.hist(peak_vs_channels)
    plt.title("Peak vs Channels SN>{}".format(snr_limit))
    plt.ylabel("Count")
    plt.xlabel("Peak Flux / Channels (Jy/beam/channel)")
    plt.show()

    return NotImplementedError


channel_sn_hist(aspecs_table)

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

print("Min Delt_Z: ", np.min(np.abs(reds_arr)))
print("Max Delt_Z: ", np.max(np.abs(reds_arr)))
print("Mean Delt_Z: ", np.mean(np.abs(reds_arr)))

plt.hist(reds_arr)
plt.title("Delta Z Distribution S/N > {}".format(snr_limit))
plt.xlabel("Delta Z")
plt.ylabel("Counts")
plt.savefig("Delta_Z_Dist_SN_{}_ZSep_{}.png".format(snr_limit, z_sep), bbox_inches='tight', dpi=300)
plt.show()

#exit()

print(np.unique(sub_arr, return_counts=True))

catalog_ids_used = [i[1] for i in catalog_ids]
print("Len catalog_ids: ", len(catalog_ids_used))
print("Len Catalog Ids: ", len(catalog_ids))

aspecs_catalog = initial_catalog[catalog_ids_used]
print("Len ASPECS Catalog: ", len(aspecs_catalog))


#aspecs_catalog.write("aspecs_SN6_ZSPEC0.3.fits", format='fits')

print("Delt Z Mean: ", np.mean(delt_z_stuff))
print("Delt Z STD: ", np.std(delt_z_stuff))
print("Delt Z Min: ", np.min(delt_z_stuff))
print("Delt Z Max: ", np.max(delt_z_stuff))
'''
#exit()
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

fig, ax = plt.subplots()

ax.axvspan(0.0030, 0.3694, facecolor='grey', alpha=0.5)
ax.axvspan(1.0059, 1.7387, facecolor='grey', alpha=0.5)
ax.axvspan(2.0088, 3.1080, facecolor='grey', alpha=0.5)
ax.axvspan(3.0115, 4.4771, facecolor='grey', alpha=0.5)
#ax.axvspan(4.0142, 5.8460, facecolor='grey', alpha=0.5)
#ax.axvspan(5.0166, 7.2146, facecolor='grey', alpha=0.5)
#ax.axvspan(6.0188, 7.3750, facecolor='grey', alpha=0.5)
#ax.axvspan(3.2823, 4.8468, facecolor='red', alpha=0.5)
#ax.axvspan(6.0422, 8.6148, facecolor='red', alpha=0.5)

ax.text((0.0030 + 0.3694)/2, 750, "1-0")
ax.text((1.0059 + 1.7387)/2, 750, "2-1")
ax.text((2.0088 + 3.1080)/2, 750, "3-2")
ax.text((3.0115 + 4.4771)/2, 750, "4-3")
#ax.text((4.0142 + 5.8460)/2, 750, "5-4")
#ax.text((5.0166 + 7.2146)/2, 750, "6-5")
#ax.text((6.0188 + 7.3750)/2, 750, "7-6")

#ax.text((3.2823 + 4.8468)/2, 900, "C1 1-0")
#ax.text((6.0422 + 8.6148)/2, 900, "C1 2-1")

ax.hist(x, 50, histtype='bar', label=('Roberto Prior', 'ASPECS No Prior'), color=('blue', 'green'), stacked=True, log=True)


ax.set_xlabel("Redshift (Z)")
ax.set_ylabel("Number of Detected Sources")
ax.set_title("Detected Sources S/N > " + str(np.round(snr_limit, 2)))
fig.legend(loc='best')
fig.savefig("AllCOLines_Z44_Limit_SN" + str(np.round(snr_limit, 2)) + ".png", bbox_inches='tight', dpi=300)
fig.show()

fig, ax = plt.subplots()

x = (spec_prior_z, spec_no_prior)
ax.axvspan(0.0030, 0.3694, facecolor='grey', alpha=0.5)
ax.axvspan(1.0059, 1.7387, facecolor='grey', alpha=0.5)
ax.axvspan(2.0088, 3.1080, facecolor='grey', alpha=0.5)
ax.axvspan(3.0115, 4.4771, facecolor='grey', alpha=0.5)
#ax.axvspan(4.0142, 5.8460, facecolor='grey', alpha=0.5)
#ax.axvspan(5.0166, 7.2146, facecolor='grey', alpha=0.5)
#ax.axvspan(6.0188, 7.3750, facecolor='grey', alpha=0.5)
#ax.axvspan(3.2823, 4.8468, facecolor='red', alpha=0.5)
#ax.axvspan(6.0422, 8.6148, facecolor='red', alpha=0.5)

ax.text((0.0030 + 0.3694)/2, 100, "1-0")
ax.text((1.0059 + 1.7387)/2, 100, "2-1")
ax.text((2.0088 + 3.1080)/2, 100, "3-2")
ax.text((3.0115 + 4.4771)/2, 100, "4-3")
#ax.text((4.0142 + 5.8460)/2, 100, "5-4")
#ax.text((5.0166 + 7.2146)/2, 100, "6-5")
#ax.text((6.0188 + 7.3750)/2, 100, "7-6")

#ax.text((3.2823 + 4.8468)/2, 200, "C1 1-0")
#ax.text((6.0422 + 8.6148)/2, 200, "C1 2-1")

ax.hist(x, 50, histtype='bar', label=('Spec Prior', 'Photo Prior'), color=('blue', 'green'), stacked=True, log=True)


ax.set_xlabel("Redshift (Z)")
ax.set_ylabel("Number of Detected Sources")
ax.set_title("Detected Sources S/N > " + str(np.round(snr_limit, 2)))
fig.legend(loc='best')
fig.savefig("OnlyMatchedCOLines_Z44_Limit_SN" + str(np.round(snr_limit, 2)) + ".png", bbox_inches='tight', dpi=300)
fig.show()
'''
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
    z_mask = (initial_catalog["z_1"] >= low_z) & (initial_catalog["z_1"] <= high_z)
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

#TODO: Create SFR M* slope from papers:
"""

Schreiber C. et al. 2015

log_10(SFR) = m - m_0 + a_0*r - a_1*[max(0, m-m_1-a_2*r)]^2

with m_0 = 0.5 +-0.07, a_0 = 1.5 +- 0.15, a_ = 0.3 +- 0.08, m_1 = 0.36 +- 0.3, and a_2 = 2.5 +- 0.6

r = log_10(1+z), and m = log_10(M*/10^9 M_solar)


Second one is Whitaker, KE, et al. 2014

log(SFR) = a + b*log(M*/M_sun) + c*log(M*/M_sun)^2

0.5 < z < 1.0: a = -17.4+-1.91 b = 5.02 +- 0.39 c = -0.22+-0.02
1.0 < z < 1.5 a = -26.03 +- 1.69 b = 4.62 +- 0.34 c = -0.19 +- 0.02
1.5 < z < 2.0 -24.04 +- 2.08 b = 4.17 +- 0.40 c = -0.16 +- 0.02
2.0 < z < 2.5 a = -19.99 +- 1.87 b = 3.44 +- 0.36 c = -0.13 +- 0.02


"""

def whitaker_main_sequence(z, mass_start, mass_end):
    mass_range = np.linspace(mass_start, mass_end, 100)

    sfr = -19.99 + 3.44*mass_range + -0.13*mass_range**2

    return mass_range, sfr


def schrieber_main_sequence(z, mass_start, mass_end):
    """"
    Returns the x and y values for the Schreiber C. et al. 2015 Main Sequence for galaxies

    Because M* and SFR are given in log10 space anyway, need to do 10^value for linspace

    But not for SFR output, because it equals log10(SFR) which is what is given
    """

    r = np.log10(1+z)

    mass_range = np.linspace(10**mass_start, 10**mass_end, 100)

    print(mass_start)
    print(np.min(mass_range))
    print(mass_end)
    print(np.max(mass_range))


    m_0 = 0.5 #+-0.07
    a_0 = 1.5 #+- 0.15
    a_1 = 0.3 #+- 0.08
    m_1 = 0.36 # +- 0.3
    a_2 = 2.5 # +- 0.6

    m_s = np.log10(mass_range)
    print(np.min(m_s))
    print(np.max(m_s))

    sfr = []
    for m in m_s:
        if m-m_1-a_2*r > 0:
            sfr.append((m - m_0 + a_0*r - a_1*((m-m_1-a_2*r)**2)))
        else:
            sfr.append((m - m_0 + a_0*r))

    return m_s, sfr



print("Len ASPECS Catalog: ", len(aspecs_catalog))
print(aspecs_catalog)
#aspecs_catalog = perform_cuts(aspecs_catalog)
initial_catalog = perform_cuts(initial_catalog)

spec_z_mask = create_spec_z_mask(initial_catalog)
spec_catalog = initial_catalog[spec_z_mask]

muse_z_mask = (initial_catalog['zm_ina'] > 0.001) | (initial_catalog['zm_her'] > 0.001) | (initial_catalog['muse_wide_z'] > 0.0001)

muse_catalog = spec_catalog #initial_catalog[muse_z_mask]
# ASPECS Ones

aspecs_low_z_sfr, alow_z_sfr_error, low_z_sfr_z = create_points_and_error_by_z("SFR", aspecs_catalog, 0, 0.4)
aspecs_mid_z_sfr, amid_z_sfr_error, mid_z_sfr_z = create_points_and_error_by_z("SFR", aspecs_catalog, 1.1, 1.8)
aspecs_high_z_sfr, ahigh_z_sfr_error, high_z_sfr_z = create_points_and_error_by_z("SFR", aspecs_catalog, 2.2, 3)
aspecs_vhigh_z_sfr, avhigh_z_sfr_error, vhigh_z_sfr_z = create_points_and_error_by_z("SFR", aspecs_catalog, 3, 4.4)

aspecs_low_z_mass, alow_z_mass_error, low_z_mass_z = create_points_and_error_by_z("Mstar", aspecs_catalog, 0, 0.4)
aspecs_mid_z_mass, amid_z_mass_error, mid_z_mass_z = create_points_and_error_by_z("Mstar", aspecs_catalog, 1.1, 1.8)
aspecs_high_z_mass, ahigh_z_mass_error, high_z_mass_z = create_points_and_error_by_z("Mstar", aspecs_catalog, 2.2, 3)
aspecs_vhigh_z_mass, avhigh_z_mass_error, vhigh_z_mass_z = create_points_and_error_by_z("Mstar", aspecs_catalog, 3, 4.4)

print("Lenghth of ASPECS VHigh: {}".format(len(aspecs_vhigh_z_sfr)))
print("Lenghth of ASPECS High: {}".format(len(aspecs_high_z_sfr)))
print("Lenghth of ASPECS Mid: {}".format(len(aspecs_mid_z_sfr)))
print("Lenghth of ASPECS Low: {}".format(len(aspecs_low_z_sfr)))

low_z_sfr, low_z_sfr_error, low_z_sfr_z = create_points_and_error_by_z("SFR", initial_catalog, 0, 0.4)
mid_z_sfr, mid_z_sfr_error, mid_z_sfr_z = create_points_and_error_by_z("SFR", initial_catalog, 1.1, 1.8)
high_z_sfr, high_z_sfr_error, high_z_sfr_z = create_points_and_error_by_z("SFR", initial_catalog, 2.2, 3)
vhigh_z_sfr, vhigh_z_sfr_error, vhigh_z_sfr_z = create_points_and_error_by_z("SFR", initial_catalog, 3, 4.4)

low_z_mass, low_z_mass_error, low_z_mass_z = create_points_and_error_by_z("Mstar", initial_catalog, 0, 0.4)
mid_z_mass, mid_z_mass_error, mid_z_mass_z = create_points_and_error_by_z("Mstar", initial_catalog, 1.1, 1.8)
high_z_mass, high_z_mass_error, high_z_mass_z = create_points_and_error_by_z("Mstar", initial_catalog, 2.2, 3)
vhigh_z_mass, vhigh_z_mass_error, vhigh_z_mass_z = create_points_and_error_by_z("Mstar", initial_catalog, 3, 4.4)

# Get the Main Sequence for the general population

schreiber_mass_11, schreiber_sfr_11 = whitaker_main_sequence(2.2, np.min(high_z_mass), np.max(high_z_mass))

schreiber_mass_18, schreiber_sfr_18 = whitaker_main_sequence(3.0, np.min(high_z_mass), np.max(high_z_mass))

# Compare Values properties M*, SFR, etc.
full_sfr = initial_catalog['SFR_50_1']
full_mstar = initial_catalog['Mstar_50_1']
#full_sfr = np.append(full_sfr, [high_z_sfr, mid_z_sfr, vhigh_z_sfr])
#full_mstar = np.append(full_mstar, [mid_z_mass, high_z_mass, vhigh_z_mass])

aspecs_mfr = aspecs_catalog["Mstar_50_1"]
aspecs_sfr = aspecs_catalog["SFR_50_1"]

print("Mstar ASPECS Mean: {} Std: {} | All Mean: {} Std: {}".format(
    np.round(np.nanmean(aspecs_mfr), 3), np.round(np.nanstd(aspecs_mfr),3), np.round(np.mean(full_mstar),3), np.round(np.std(full_mstar), 3)))
print("SFR ASPECS Mean: {} Std: {} | All Mean: {} Std: {}".format(
    np.round(np.nanmean(aspecs_sfr), 3), np.round(np.nanstd(aspecs_sfr),3), np.round(np.mean(full_sfr),3), np.round(np.std(full_sfr), 3)))



# Just the All

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='all', sharey='all')
ax1.errorbar(low_z_mass, low_z_sfr, yerr=low_z_sfr_error, xerr=low_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='darkgrey', elinewidth=1)
ax1.plot(np.unique(low_z_mass), np.poly1d(np.polyfit(low_z_mass, low_z_sfr, 1))(np.unique(low_z_mass)), label='All fit',
         color='black', zorder=10)
ax1.errorbar(aspecs_low_z_mass, aspecs_low_z_sfr, yerr=alow_z_sfr_error, xerr=alow_z_mass_error, ecolor='red', fmt='.',
            ms=5, mec='red', label='ASPECS', zorder=20)
ax1.set_title('0 < Z < 0.4')
ax2.errorbar(mid_z_mass, mid_z_sfr, yerr=mid_z_sfr_error, xerr=mid_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='darkgrey', elinewidth=1)
ax2.plot(np.unique(mid_z_mass), np.poly1d(np.polyfit(mid_z_mass, mid_z_sfr, 1))(np.unique(mid_z_mass)), color='black',
         zorder=10)
ax2.errorbar(aspecs_mid_z_mass, aspecs_mid_z_sfr, yerr=amid_z_sfr_error, xerr=amid_z_mass_error, ecolor='red', fmt='.',
             ms=5, mec='red', zorder=20)
ax2.set_title('1.1 < Z < 1.8')
ax3.errorbar(high_z_mass, high_z_sfr, yerr=high_z_sfr_error, xerr=high_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='darkgrey', elinewidth=1)
ax3.plot(np.unique(high_z_mass), np.poly1d(np.polyfit(high_z_mass, high_z_sfr, 1))(np.unique(high_z_mass)),
         color='black',
         zorder=10)
ax3.errorbar(aspecs_high_z_mass, aspecs_high_z_sfr, yerr=ahigh_z_sfr_error, xerr=ahigh_z_mass_error, ecolor='red', fmt='.',
             ms=5, mec='red', zorder=20)

ax3.plot(schreiber_mass_11, schreiber_sfr_11, color='orange', zorder=20)
ax3.plot(schreiber_mass_18, schreiber_sfr_18, color='green', zorder=20)
ax3.set_title('2.2 < Z < 3')
ax4.errorbar(vhigh_z_mass, vhigh_z_sfr, yerr=vhigh_z_sfr_error, xerr=vhigh_z_mass_error, ecolor='lightgrey', fmt='.',
             ms=1, mec='darkgrey', elinewidth=1)
ax4.plot(np.unique(vhigh_z_mass), np.poly1d(np.polyfit(vhigh_z_mass, vhigh_z_sfr, 1))(np.unique(vhigh_z_mass)),
         color='black', zorder=10)
ax4.errorbar(aspecs_vhigh_z_mass, aspecs_vhigh_z_sfr, yerr=avhigh_z_sfr_error, xerr=avhigh_z_mass_error, ecolor='red', fmt='.',
             ms=5, mec='red', zorder=20)
# ax4.errorbar(vhigh_z_mass, vhigh_z_sfr, yerr=vhigh_z_sfr_error, xerr=vhigh_z_mass_error)
ax4.set_title('3 < Z < 4.4')
handles, labels = ax1.get_legend_handles_labels()
f.legend(loc='best', handles=handles, labels=labels, prop={'size': 6})
f.text(0.5, 0.01, 'Log(M*)', ha='center')
f.text(0.01, 0.5, 'Log(SFR)', va='center', rotation='vertical')
f.savefig("Wide_ASPECS_Final_MstarVsSFR_Limit_SN" + str(np.round(snr_limit,3)) + ".png", bbox_inches='tight', dpi=300)
f.show()

exit()

muse_low_z_sfr, low_z_sfr_error, low_z_sfr_z = create_points_and_error_by_z("SFR", muse_catalog, 0, 0.4)
muse_mid_z_sfr, mid_z_sfr_error, mid_z_sfr_z = create_points_and_error_by_z("SFR", muse_catalog, 1.1, 1.8)
muse_high_z_sfr, high_z_sfr_error, high_z_sfr_z = create_points_and_error_by_z("SFR", muse_catalog, 2.2, 3)
muse_vhigh_z_sfr, vhigh_z_sfr_error, vhigh_z_sfr_z = create_points_and_error_by_z("SFR", muse_catalog, 3, 4.4)

muse_low_z_mass, low_z_mass_error, low_z_mass_z = create_points_and_error_by_z("Mstar", muse_catalog, 0, 0.4)
muse_mid_z_mass, mid_z_mass_error, mid_z_mass_z = create_points_and_error_by_z("Mstar", muse_catalog, 1.1, 1.8)
muse_high_z_mass, high_z_mass_error, high_z_mass_z = create_points_and_error_by_z("Mstar", muse_catalog, 2.2, 3)
muse_vhigh_z_mass, vhigh_z_mass_error, vhigh_z_mass_z = create_points_and_error_by_z("Mstar", muse_catalog, 3, 4.4)

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='all')
ax1.hist(low_z_mass, color='black', histtype='step', bins=10)
ax1.set_yscale('log')
ax1.set_title('0 < Z < 0.4')
ax1.hist(aspecs_low_z_mass, color='green', histtype='step', bins=10)
ax1.hist(muse_low_z_mass, color='blue', histtype='step', bins=10)
ax2.set_title('0 < Z < 0.4')
ax2.hist(low_z_sfr, color='black', histtype='step', bins=10)
ax2.hist(aspecs_low_z_sfr, color='green', histtype='step', bins=10)
ax2.hist(muse_low_z_sfr, color='blue', histtype='step', bins=10)
ax3.set_title('1.1 < Z < 1.8')
ax3.hist(mid_z_mass, color='black', histtype='step', bins=10)
ax3.hist(aspecs_mid_z_mass, color='green', histtype='step', bins=10)
ax3.hist(muse_mid_z_mass, color='blue', histtype='step', bins=10)
ax4.set_title('1.1 < Z < 1.8')
ax4.hist(mid_z_sfr, color='black', histtype='step', label='All', bins=10)
ax4.hist(aspecs_mid_z_sfr, color='green', histtype='step', label='ASPECS Prior', bins=10)
ax4.hist(muse_mid_z_sfr, color='blue', histtype='step', label='MUSE', bins=10)
handles, labels = ax4.get_legend_handles_labels()
f.legend(loc='best', handles=handles, labels=labels, prop={'size': 6})
f.text(0.25, 0.01, 'Log Stellar Mass (Mstar)', ha='center')
f.text(0.75, 0.01, 'Log SFR (Mstar)', ha='center')
f.text(0.01, 0.25, 'Count', va='center', rotation='vertical')
f.text(0.01, 0.75, 'Count', va='center', rotation='vertical')
f.savefig("Wide_ASPECS_LeinerdtPlot01_Limit_SN" + str(np.round(snr_limit,3)) + ".png", bbox_inches='tight', dpi=300)
f.show()

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='all')
ax1.hist(high_z_mass, color='black', histtype='step', bins=10)
ax1.set_yscale('log')
ax1.set_title('2.2 < Z < 3')
ax1.hist(aspecs_high_z_mass, color='green', histtype='step', bins=10)
ax1.hist(muse_high_z_mass, color='blue', histtype='step', bins=10)
ax2.set_title('2.2 < Z < 3')
ax2.hist(high_z_sfr, color='black', histtype='step', bins=10)
ax2.hist(aspecs_high_z_sfr, color='green', histtype='step', bins=10)
ax2.hist(muse_high_z_sfr, color='blue', histtype='step', bins=10)
ax3.set_title('3 < Z < 4.4')
ax3.hist(vhigh_z_mass, color='black', histtype='step', bins=10)
ax3.hist(aspecs_vhigh_z_mass, color='green', histtype='step', bins=10)
ax3.hist(muse_vhigh_z_mass, color='blue', histtype='step', bins=10)
ax4.set_title('3 < Z < 4.4 ')
ax4.hist(vhigh_z_sfr, color='black', histtype='step', label='All', bins=10)
ax4.hist(aspecs_vhigh_z_sfr, color='green', histtype='step', label='ASPECS', bins=10)
ax4.hist(muse_vhigh_z_sfr, color='blue', histtype='step', bins=10, label='MUSE')
handles, labels = ax4.get_legend_handles_labels()
f.legend(loc='best', handles=handles, labels=labels, prop={'size': 6})
f.text(0.25, 0.01, 'Log Stellar Mass (Mstar)', ha='center')
f.text(0.75, 0.01, 'Log SFR (Mstar)', ha='center')
f.text(0.01, 0.25, 'Count', va='center', rotation='vertical')
f.text(0.01, 0.75, 'Count', va='center', rotation='vertical')
f.savefig("Wide_ASPECS_LeinerdtPlot02_Limit_SN" + str(np.round(snr_limit,3)) + ".png", bbox_inches='tight', dpi=300)
f.show()


exit()
