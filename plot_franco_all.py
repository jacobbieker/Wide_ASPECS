import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import vstack
from astropy.table import Table, hstack
import astropy.units as u

from scipy import stats
from numpy.polynomial.polynomial import polyfit
from aspecs_catalog_builder import compare_catalog_locations

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


franco_catalog = fits.open("/home/jacob/Research/MAGPHYS/franco/catalog.fits")

franco_catalog = franco_catalog[1].data

hdu_list = fits.open("data/jacob_aspecs_catalog_fixed_magphys_jcb3.fits")
print(hdu_list[1].columns)
full_catalog = hdu_list[1].data

# List of the fits files to open
roberto_magphys_output = ["/home/jacob/Research/MAGPHYS/roberto2/catalog.fits",
                          "/home/jacob/Research/MAGPHYS/roberto10/catalog.fits",
                          "/home/jacob/Research/MAGPHYS/roberto3/catalog.fits",
                          "/home/jacob/Research/MAGPHYS/roberto4/catalog.fits",
                          "/home/jacob/Research/MAGPHYS/roberto5/catalog.fits",
                          "/home/jacob/Research/MAGPHYS/roberto6/catalog.fits",
                          "/home/jacob/Research/MAGPHYS/roberto7/catalog.fits",
                          "/home/jacob/Research/MAGPHYS/roberto8/catalog.fits",
                          "/home/jacob/Research/MAGPHYS/roberto9/catalog.fits"]

roberto_catalog = Table.read("/home/jacob/Research/MAGPHYS/roberto1/catalog.fits", format='fits')

for filename in roberto_magphys_output:
    hdu_list = Table.read(filename, format='fits')
    roberto_catalog = vstack([roberto_catalog, hdu_list])

spec_z_mask = (full_catalog["z_spec_3dh"] > 0.001) | (full_catalog["zm_vds"] > 0.001) | (
        full_catalog["zm_coeS"] > 0.001) \
              | (full_catalog["zs_mor"] > 0.001) | (full_catalog["zm_ina"] > 0.001) | (full_catalog["zm_her"] > 0.001)

muse_z_mask = (full_catalog['zm_ina'] > 0.001) | (full_catalog['zm_her'] > 0.001)

spec_catalog = full_catalog[spec_z_mask]
muse_catalog = full_catalog[muse_z_mask]


#leinhardt_catalog = compare_catalog_locations(os.path.join("data", "jacob_aspecs_catalog_fixed_magphys_jcb3.fits"), os.path.join("data", "MW_44fields_main_table_v1.0.fits"))

#leinhardt_muse_z_mask = (leinhardt_catalog['zm_ina'] > 0.001) | (leinhardt_catalog['zm_her'] > 0.001) | (leinhardt_catalog['muse_wide_z'] > 0.0001)

#leinhardt_catalog = leinhardt_catalog[leinhardt_muse_z_mask]

roberto_muse = Table.read("roberto_catalog_muse.fits", format='fits')
test_roberto = Table.read("/home/jacob/Development/Wide_ASPECS/mapghys_in_nov2018_all.fits", format='fits')
# Add in RA and Dec to test_roberto
from astropy.table import join

roberto_muse = join(test_roberto, roberto_muse, keys='id')

roberto_matched_ids = [[], [23419], [21510], [], [], [], [], [14630], [18282], [17508], [9258], [], [], [], [], [], [], [], [], [], [22705], [22705], [], [], [13131], [], [], [], [], [], [], [], [14957], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [28901], [17170], [17170], [28851], [52578], [58258], [52994], [], [58246], [58246], [58246], [58246], [58246], [62572], [62572], [], [], [], [], [], [25004], [], [], [], [], [], [], [], [], [], [52367], [], [], [], [57822], [], [], [22160], [], [], [18553], [18553], [19491]]

results_dict = {1: (1.094, 218.71 * u.GHz, 230.538 *u.GHz, 23419),
                2: (1.675, 279.602 * u.GHz, 230.538 *u.GHz, 21510),
                7: (0.085, 101.745 * u.GHz, 115.271 *u.GHz, 14630),
                8: (0.872, 177.885 * u.GHz, 115.271 *u.GHz, 18282),
                9: (1.685, 249.474 * u.GHz, 230.538 *u.GHz, 17508),
                10: (0.925, 181.568 * u.GHz, 115.271 *u.GHz, 9258),
                20: (0.96, 209.308 * u.GHz, 115.271 *u.GHz,  22705),
                21: (0.96, 209.338 * u.GHz, 115.271 *u.GHz, 22705),
                24: (3.11, 377.384 * u.GHz, 345.796 * u.GHz, 13131),
                32: (0.163, 110.94 * u.GHz, 115.271 *u.GHz, 14957),
                49: (1.191, 234.301 * u.GHz, 230.538 *u.GHz, 28901),
                50: (1.535, 235.043 * u.GHz, 230.538 *u.GHz, 17170),
                51: (1.535, 235.063 * u.GHz, 230.538 *u.GHz, 17170),
                52: (1.019, 211.949 * u.GHz, 230.538 *u.GHz, 28851),
                53: (1.537, 241.314 * u.GHz, 230.538 *u.GHz, 52578),
                54: (2.835, 364.747 * u.GHz, 345.796 *u.GHz,  58258),
                55: (3.496, 420.165 * u.GHz, 461.041 * u.GHz, 52994),
                57: (1.597, 243.105 * u.GHz, 230.538 *u.GHz, 58246),
                58: (1.597, 243.084 * u.GHz, 230.538 *u.GHz,58246),
                59: (1.597, 243.084 * u.GHz, 230.538 *u.GHz,58246),
                60: (1.597, 243.126 * u.GHz, 230.538 *u.GHz, 58246),
                61: (1.597, 243.126 * u.GHz, 230.538 *u.GHz,58246),
                62: (1.599, 246.073 * u.GHz, 230.538 *u.GHz, 62572),
                63: (1.599, 246.053 * u.GHz, 230.538 *u.GHz, 62572),
                69: (0.738, 184.459 * u.GHz, 115.271 *u.GHz, 25004),
                79: (1.269, 211.816 * u.GHz, 230.538 *u.GHz, 52367),
                83: (1.599, 246.359 * u.GHz, 230.538 *u.GHz, 57822),
                86: (0.355, 143.493 * u.GHz, 115.271 *u.GHz, 22160),
                89: (2.32, 311.433 * u.GHz, 345.796 *u.GHz, 18553),
                90: (2.32, 311.406 * u.GHz, 345.796 *u.GHz, 18553),
                91: (1.665, 246.118 * u.GHz, 230.538 *u.GHz, 19491)}



GHz_diff = {}
for key, value in results_dict.items():
    GHz_diff[key] = value[1] - value[2]

print(GHz_diff)

aspecs_number= []
roberto_id = []
for key, value in results_dict.items():
    roberto_id.append(value[2])
    aspecs_number.append(key)

aspecs_catalog_mask = (full_catalog['id'] == 19491) | (full_catalog['id'] == 18553) | (full_catalog['id'] == 22160) | \
                      (full_catalog['id'] == 57822) | (full_catalog['id'] == 52367) | (full_catalog['id'] == 25004) | \
                      (full_catalog['id'] == 62572) | (full_catalog['id'] == 58246) | (full_catalog['id'] == 52994) | \
                      (full_catalog['id'] == 58258) | (full_catalog['id'] == 52578) | (full_catalog['id'] == 28851) | \
                      (full_catalog['id'] == 17170) | (full_catalog['id'] == 28901) | (full_catalog['id'] == 14957) | \
                      (full_catalog['id'] == 28901) | (full_catalog['id'] == 14957) | (full_catalog['id'] == 13131) | \
                      (full_catalog['id'] == 22705) | (full_catalog['id'] == 9258) | (full_catalog['id'] == 17508) | \
                      (full_catalog['id'] == 18282) | (full_catalog['id'] == 14630) | (full_catalog['id'] == 21510) | \
                      (full_catalog['id'] == 23419)

aspecs_catalog = full_catalog[aspecs_catalog_mask]
#leinhardt_catalog = perform_cuts(leinhardt_catalog)

# Now the quality cuts
#aspecs_catalog = perform_cuts(aspecs_catalog)
print("Lenghth of ASPECS: {}".format(len(aspecs_catalog)))
#exit()
roberto_catalog = perform_cuts(roberto_catalog)
spec_catalog = perform_cuts(spec_catalog)
franco_catalog = perform_cuts(franco_catalog)
full_catalog = perform_cuts(full_catalog)
muse_catalog = perform_cuts(muse_catalog)
count = 0
spec_ids = full_catalog['id']
diff_ids = []
for id in roberto_catalog['id']:
    if id in spec_ids:
        diff_ids.append(id)
# now create smaller catalog with full catalog info

rows_to_use = []
for index, row in enumerate(full_catalog):
    if row['id'] in diff_ids:
        rows_to_use.append(index)

smaller_catalog = full_catalog[rows_to_use]

# Plot on the skymap to see if localized there
plt.scatter(full_catalog['ra'], full_catalog['dc'], label='All', s=2)
plt.scatter(smaller_catalog['ra'], smaller_catalog['dc'], label='Roberto', s=2)
plt.ylabel("Dec")
plt.xlabel("Ra")
plt.xlim(np.max(full_catalog['ra']), np.min(full_catalog['ra']))
plt.legend()
plt.savefig("RobertoCatalog_SkyMap.pdf", dpi=300)
plt.show()
print("Fraction of differing galaxies (Roberto): " + str(len(diff_ids) / len(roberto_catalog)))
print("Fraction of differing galaxies (Spec): " + str(len(diff_ids) / len(spec_catalog)))
av, av_error = create_points_and_error("A_V", full_catalog)
rob_av, rob_av_error = create_points_and_error("A_V", roberto_catalog)
spec_av, spec_av_error = create_points_and_error("A_V", spec_catalog)

plt.hist(av, histtype='step', color='lightgrey', label='All')
plt.hist(rob_av, histtype='step', color='black', label='Roberto')
plt.hist(spec_av, histtype='step', color='yellow', label='Spec')
plt.yscale('log')
plt.xlabel("Av")
plt.ylabel("Count")
plt.title("Av Distribution")
plt.legend(loc='best')
plt.savefig("Av_Distribution.png", bbox_inches='tight', dpi=300)
plt.savefig("Av_Distribution.pdf", bbox_inches='tight', dpi=300)
plt.show()

# Now plotting the three of them in various ways

# Now only the Spectroscopic ones
spec_low_z_sfr, low_z_sfr_error, low_z_sfr_z = create_points_and_error_by_z("SFR", spec_catalog, 0, 1)
spec_mid_z_sfr, mid_z_sfr_error, mid_z_sfr_z = create_points_and_error_by_z("SFR", spec_catalog, 1, 2)
spec_high_z_sfr, high_z_sfr_error, high_z_sfr_z = create_points_and_error_by_z("SFR", spec_catalog, 2, 3)
spec_vhigh_z_sfr, vhigh_z_sfr_error, vhigh_z_sfr_z = create_points_and_error_by_z("SFR", spec_catalog, 3, 4)

spec_low_z_mass, low_z_mass_error, low_z_mass_z = create_points_and_error_by_z("Mstar", spec_catalog, 0, 1)
spec_mid_z_mass, mid_z_mass_error, mid_z_mass_z = create_points_and_error_by_z("Mstar", spec_catalog, 1, 2)
spec_high_z_mass, high_z_mass_error, high_z_mass_z = create_points_and_error_by_z("Mstar", spec_catalog, 2, 3)
spec_vhigh_z_mass, vhigh_z_mass_error, vhigh_z_mass_z = create_points_and_error_by_z("Mstar", spec_catalog, 3, 4)

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

# Now the Roberto Ones

rob_low_z_sfr, low_z_sfr_error, low_z_sfr_z = create_points_and_error_by_z("SFR", roberto_catalog, 0, 1)
rob_mid_z_sfr, mid_z_sfr_error, mid_z_sfr_z = create_points_and_error_by_z("SFR", roberto_catalog, 1, 2)
rob_high_z_sfr, high_z_sfr_error, high_z_sfr_z = create_points_and_error_by_z("SFR", roberto_catalog, 2, 3)
rob_vhigh_z_sfr, vhigh_z_sfr_error, vhigh_z_sfr_z = create_points_and_error_by_z("SFR", roberto_catalog, 3, 4)

rob_low_z_mass, low_z_mass_error, low_z_mass_z = create_points_and_error_by_z("Mstar", roberto_catalog, 0, 1)
rob_mid_z_mass, mid_z_mass_error, mid_z_mass_z = create_points_and_error_by_z("Mstar", roberto_catalog, 1, 2)
rob_high_z_mass, high_z_mass_error, high_z_mass_z = create_points_and_error_by_z("Mstar", roberto_catalog, 2, 3)
rob_vhigh_z_mass, vhigh_z_mass_error, vhigh_z_mass_z = create_points_and_error_by_z("Mstar", roberto_catalog, 3, 4)

# Now the Franco data ones

franco_low_z_sfr, low_z_sfr_error, low_z_sfr_z = create_points_and_error_by_z("SFR", franco_catalog, 0, 1)
franco_mid_z_sfr, mid_z_sfr_error, mid_z_sfr_z = create_points_and_error_by_z("SFR", franco_catalog, 1, 2)
franco_high_z_sfr, high_z_sfr_error, high_z_sfr_z = create_points_and_error_by_z("SFR", franco_catalog, 2, 3)
franco_vhigh_z_sfr, vhigh_z_sfr_error, vhigh_z_sfr_z = create_points_and_error_by_z("SFR", franco_catalog, 3, 4)

franco_low_z_mass, low_z_mass_error, low_z_mass_z = create_points_and_error_by_z("Mstar", franco_catalog, 0, 1)
franco_mid_z_mass, mid_z_mass_error, mid_z_mass_z = create_points_and_error_by_z("Mstar", franco_catalog, 1, 2)
franco_high_z_mass, high_z_mass_error, high_z_mass_z = create_points_and_error_by_z("Mstar", franco_catalog, 2, 3)
franco_vhigh_z_mass, vhigh_z_mass_error, vhigh_z_mass_z = create_points_and_error_by_z("Mstar", franco_catalog, 3, 4)

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
            s=30, c='red',  label='ASPECS', zorder=20)
ax1.set_title('0 < Z < 1')
ax2.errorbar(mid_z_mass, mid_z_sfr, yerr=mid_z_sfr_error, xerr=mid_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='darkgrey', elinewidth=1)
ax2.plot(np.unique(mid_z_mass), np.poly1d(np.polyfit(mid_z_mass, mid_z_sfr, 1))(np.unique(mid_z_mass)), color='black',
         zorder=10)
ax2.scatter(aspecs_mid_z_mass, aspecs_mid_z_sfr, marker='.',
             s=30, c='red',zorder=20)
ax2.set_title('1 < Z < 2')
ax3.errorbar(high_z_mass, high_z_sfr, yerr=high_z_sfr_error, xerr=high_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='darkgrey', elinewidth=1)
ax3.plot(np.unique(high_z_mass), np.poly1d(np.polyfit(high_z_mass, high_z_sfr, 1))(np.unique(high_z_mass)), color='black',
         zorder=10)
ax3.scatter(aspecs_high_z_mass, aspecs_high_z_sfr, marker='.',
             s=30, c='red',zorder=20)
ax3.set_title('2 < Z < 3')
ax4.errorbar(vhigh_z_mass, vhigh_z_sfr, yerr=vhigh_z_sfr_error, xerr=vhigh_z_mass_error, ecolor='lightgrey', fmt='.',
             ms=1, mec='darkgrey', elinewidth=1)
ax4.plot(np.unique(vhigh_z_mass), np.poly1d(np.polyfit(vhigh_z_mass, vhigh_z_sfr, 1))(np.unique(vhigh_z_mass)),
         color='black', zorder=10)
ax4.scatter(aspecs_vhigh_z_mass, aspecs_vhigh_z_sfr, marker='.',
             s=30, c='red',zorder=20)
# ax4.errorbar(vhigh_z_mass, vhigh_z_sfr, yerr=vhigh_z_sfr_error, xerr=vhigh_z_mass_error)
ax4.set_title('3 < Z < 4')
handles, labels = ax1.get_legend_handles_labels()
f.legend(loc='best', handles=handles, labels=labels, prop={'size': 6})
f.text(0.5, 0.01, 'Log(M*)', ha='center')
f.text(0.01, 0.5, 'Log(SFR)', va='center', rotation='vertical')
f.savefig("AllMStarSFR_ASPECS_MATCHES.png", bbox_inches='tight', dpi=300)
f.show()
exit()

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='all', sharey='all')
ax1.errorbar(low_z_mass, low_z_sfr, yerr=low_z_sfr_error, xerr=low_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='lightgrey', elinewidth=1, zorder=-32)
ax1.plot(np.unique(low_z_mass), np.poly1d(np.polyfit(low_z_mass, low_z_sfr, 1))(np.unique(low_z_mass)), label='All fit',
         color='red')
ax1.scatter(spec_low_z_mass, spec_low_z_sfr, color='yellow', label='Spec', s=1)
ax1.plot(np.unique(spec_low_z_mass),
         np.poly1d(np.polyfit(spec_low_z_mass, spec_low_z_sfr, 1))(np.unique(spec_low_z_mass)), label='Spec Fit',
         color='b')
ax1.scatter(rob_low_z_mass, rob_low_z_sfr, color='black', marker='s', label='Roberto', s=1)
ax1.plot(np.unique(rob_low_z_mass), np.poly1d(np.polyfit(rob_low_z_mass, rob_low_z_sfr, 1))(np.unique(rob_low_z_mass)),
         label='Roberto Fit', color='green')
ax1.scatter(franco_low_z_mass, franco_low_z_sfr, color='yellow', marker='*', s=8)
ax1.set_title('0 < Z < 1')
ax2.errorbar(mid_z_mass, mid_z_sfr, yerr=mid_z_sfr_error, xerr=mid_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='lightgrey', elinewidth=1, zorder=-32)
ax2.plot(np.unique(mid_z_mass), np.poly1d(np.polyfit(mid_z_mass, mid_z_sfr, 1))(np.unique(mid_z_mass)), color='red')
ax2.scatter(spec_mid_z_mass, spec_mid_z_sfr, color='yellow', s=1)
ax2.plot(np.unique(spec_mid_z_mass),
         np.poly1d(np.polyfit(spec_mid_z_mass, spec_mid_z_sfr, 1))(np.unique(spec_mid_z_mass)), color='b')
ax2.scatter(rob_mid_z_mass, rob_mid_z_sfr, color='black', marker='s', s=1)
ax2.plot(np.unique(rob_mid_z_mass), np.poly1d(np.polyfit(rob_mid_z_mass, rob_mid_z_sfr, 1))(np.unique(rob_mid_z_mass)),
         color='green')
ax2.scatter(franco_mid_z_mass, franco_mid_z_sfr, color='yellow', marker='*', s=8)
ax2.set_title('1 < Z < 2')
ax3.errorbar(high_z_mass, high_z_sfr, yerr=high_z_sfr_error, xerr=high_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='lightgrey', elinewidth=1, zorder=-32)
ax3.plot(np.unique(high_z_mass), np.poly1d(np.polyfit(high_z_mass, high_z_sfr, 1))(np.unique(high_z_mass)), color='red')
ax3.scatter(spec_high_z_mass, spec_high_z_sfr, color='yellow', s=1)
ax3.plot(np.unique(spec_high_z_mass),
         np.poly1d(np.polyfit(spec_high_z_mass, spec_high_z_sfr, 1))(np.unique(spec_high_z_mass)), color='b')
ax3.scatter(rob_high_z_mass, rob_high_z_sfr, color='black', marker='s', s=1)
ax3.plot(np.unique(rob_high_z_mass),
         np.poly1d(np.polyfit(rob_high_z_mass, rob_high_z_sfr, 1))(np.unique(rob_high_z_mass)), color='green')
ax3.scatter(franco_high_z_mass, franco_high_z_sfr, color='yellow', marker='*', s=8)
ax3.set_title('2 < Z < 3')
ax4.errorbar(vhigh_z_mass, vhigh_z_sfr, yerr=vhigh_z_sfr_error, xerr=vhigh_z_mass_error, ecolor='lightgrey', fmt='.',
             ms=1, mec='lightgrey', elinewidth=1, zorder=-32)
ax4.plot(np.unique(vhigh_z_mass), np.poly1d(np.polyfit(vhigh_z_mass, vhigh_z_sfr, 1))(np.unique(vhigh_z_mass)),
         color='red')
ax4.scatter(spec_vhigh_z_mass, spec_vhigh_z_sfr, color='yellow', s=1)
ax4.plot(np.unique(spec_vhigh_z_mass),
         np.poly1d(np.polyfit(spec_vhigh_z_mass, spec_vhigh_z_sfr, 1))(np.unique(spec_vhigh_z_mass)), color='b')
ax4.scatter(rob_vhigh_z_mass, rob_vhigh_z_sfr, color='black', marker='s', s=1)
ax4.plot(np.unique(rob_vhigh_z_mass),
         np.poly1d(np.polyfit(rob_vhigh_z_mass, rob_vhigh_z_sfr, 1))(np.unique(rob_vhigh_z_mass)), color='green')
ax4.scatter(franco_vhigh_z_mass, franco_vhigh_z_sfr, color='yellow', marker='*', s=8)
# ax4.errorbar(vhigh_z_mass, vhigh_z_sfr, yerr=vhigh_z_sfr_error, xerr=vhigh_z_mass_error)
ax4.set_title('3 < Z < 4')
handles, labels = ax1.get_legend_handles_labels()
f.legend(loc='best', handles=handles, labels=labels, prop={'size': 6})
f.text(0.5, 0.01, 'Log Stellar Mass (Mstar)', ha='center')
f.text(0.01, 0.5, 'Log Star Formation Rate', va='center', rotation='vertical')
f.savefig("MstarVsSFRbyZ04.png", bbox_inches='tight', dpi=300)
f.show()
exit()
# Only Spectroscopic and Roberto

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='all', sharey='all')
ax1.scatter(spec_low_z_mass, spec_low_z_sfr, color='yellow', label='Spec', s=1)
ax1.plot(np.unique(spec_low_z_mass),
         np.poly1d(np.polyfit(spec_low_z_mass, spec_low_z_sfr, 1))(np.unique(spec_low_z_mass)), label='Spec fit',
         color='b')
ax1.scatter(rob_low_z_mass, rob_low_z_sfr, color='black', marker='s', label='Roberto', s=1)
ax1.plot(np.unique(rob_low_z_mass), np.poly1d(np.polyfit(rob_low_z_mass, rob_low_z_sfr, 1))(np.unique(rob_low_z_mass)),
         label='Roberto fit', color='green')
ax1.set_title('0 < Z < 1')
ax2.scatter(spec_mid_z_mass, spec_mid_z_sfr, color='yellow', s=1)
ax2.plot(np.unique(spec_mid_z_mass),
         np.poly1d(np.polyfit(spec_mid_z_mass, spec_mid_z_sfr, 1))(np.unique(spec_mid_z_mass)), color='b')
ax2.scatter(rob_mid_z_mass, rob_mid_z_sfr, color='black', marker='s', s=1)
ax2.plot(np.unique(rob_mid_z_mass), np.poly1d(np.polyfit(rob_mid_z_mass, rob_mid_z_sfr, 1))(np.unique(rob_mid_z_mass)),
         color='green')
ax2.set_title('1 < Z < 2')
ax3.scatter(spec_high_z_mass, spec_high_z_sfr, color='yellow', s=1)
ax3.plot(np.unique(spec_high_z_mass),
         np.poly1d(np.polyfit(spec_high_z_mass, spec_high_z_sfr, 1))(np.unique(spec_high_z_mass)), color='b')
ax3.scatter(rob_high_z_mass, rob_high_z_sfr, color='black', marker='s', s=1)
ax3.plot(np.unique(rob_high_z_mass),
         np.poly1d(np.polyfit(rob_high_z_mass, rob_high_z_sfr, 1))(np.unique(rob_high_z_mass)), color='green')
ax3.set_title('2 < Z < 3')
ax4.scatter(spec_vhigh_z_mass, spec_vhigh_z_sfr, color='yellow', s=1)
ax4.plot(np.unique(spec_vhigh_z_mass),
         np.poly1d(np.polyfit(spec_vhigh_z_mass, spec_vhigh_z_sfr, 1))(np.unique(spec_vhigh_z_mass)), color='b')
ax4.scatter(rob_vhigh_z_mass, rob_vhigh_z_sfr, color='black', marker='s', s=1)
ax4.plot(np.unique(rob_vhigh_z_mass),
         np.poly1d(np.polyfit(rob_vhigh_z_mass, rob_vhigh_z_sfr, 1))(np.unique(rob_vhigh_z_mass)), color='green')
# ax4.errorbar(vhigh_z_mass, vhigh_z_sfr, yerr=vhigh_z_sfr_error, xerr=vhigh_z_mass_error, fmt='none', color='lightgrey')
ax4.set_title('3 < Z < 4')
handles, labels = ax1.get_legend_handles_labels()
f.legend(loc='best', handles=handles, labels=labels, prop={'size': 6})
f.text(0.5, 0.01, 'Log Stellar Mass (Mstar)', ha='center')
f.text(0.01, 0.5, 'Log Star Formation Rate', va='center', rotation='vertical')
f.savefig("MstarVsSFRbyZ_RobertoAndSpec04.png", bbox_inches='tight', dpi=300)
f.show()

# Now plot them all together in one plot, not by Z

sfr, sfr_error = create_points_and_error("SFR", full_catalog)

mass, mass_error = create_points_and_error("Mstar", full_catalog)

# Now only the Spectroscopic ones
spec_sfr, spec_sfr_error = create_points_and_error("SFR", spec_catalog)
spec_mass, spec_mass_error = create_points_and_error("Mstar", spec_catalog)

# Now the Roberto Ones

rob_sfr, rob_sfr_error = create_points_and_error("SFR", roberto_catalog)

rob_mass, rob_mass_error = create_points_and_error("Mstar", roberto_catalog)

# Now the Franco data ones

franco_sfr, france_sfr_error = create_points_and_error("SFR", franco_catalog)
print(np.min(franco_sfr))
print(np.max(franco_sfr))

franco_mass, franco_mass_error = create_points_and_error("Mstar", franco_catalog)
print(np.min(franco_mass))
print(np.max(franco_mass))

plt.errorbar(mass, sfr, yerr=sfr_error, xerr=mass_error, ecolor='lightgrey', fmt='.', mec='lightgrey', ms=1,
             label='All', elinewidth=1, zorder=-32)
plt.scatter(spec_mass, spec_sfr, color='yellow', label='Spectroscopic Z', s=1)
plt.scatter(rob_mass, rob_sfr, color='black', label='Roberto', s=1)
plt.scatter(franco_mass, franco_sfr, color='yellow', marker="*", label='Franco et. al.', s=8)
plt.plot(np.unique(mass), np.poly1d(np.polyfit(mass, sfr, 1))(np.unique(mass)), color='red', label='All Fit')
plt.plot(np.unique(rob_mass), np.poly1d(np.polyfit(rob_mass, rob_sfr, 1))(np.unique(rob_mass)), color='green',
         label='Roberto fit')
plt.plot(np.unique(spec_mass), np.poly1d(np.polyfit(spec_mass, spec_sfr, 1))(np.unique(spec_mass)), color='blue',
         label='Spec Fit')
plt.legend(loc='best')
plt.ylabel('Log(SFR)')
plt.xlabel('Log(M*) [Msun]')
plt.title("M* vs SFR")
plt.savefig("MstarVsSFR.png", bbox_inches='tight', dpi=300)
plt.show()

# plt.errorbar(mass, sfr,yerr=sfr_error, xerr=mass_error, ecolor='lightgrey', fmt='.', ms=1, mfc='blue', label='All', elinewidth=1)
plt.errorbar(spec_mass, spec_sfr, yerr=spec_sfr_error, xerr=spec_mass_error, ecolor='lightgrey', fmt='.', mec='yellow',
             ms=1, label='Spectroscopic', elinewidth=1, zorder=-32)
plt.scatter(rob_mass, rob_sfr, color='black', label='Roberto', s=1)
plt.plot(np.unique(mass), np.poly1d(np.polyfit(mass, sfr, 1))(np.unique(mass)), color='red', label='All Fit')
# plt.plot(np.unique(rob_mass), np.poly1d(np.polyfit(rob_mass, rob_sfr, 1))(np.unique(rob_mass)), color='yellow1)
plt.plot(np.unique(spec_mass), np.poly1d(np.polyfit(spec_mass, spec_sfr, 1))(np.unique(spec_mass)), color='black',
         label='Spec Fit')
plt.legend(loc='best')
plt.ylabel('Log(SFR)')
plt.xlabel('Log(M*) [Msun]')
plt.title('M* vs SFR')
plt.savefig("MstarVsSFR_AllSpec.png", bbox_inches='tight', dpi=300)
plt.show()

# Only Franco
plt.scatter(franco_mass, franco_sfr, label='Franco')
plt.title("Franco et al. Stellar Mass vs SFR N = " + str(len(franco_mass)))
plt.ylabel("Log SFR")
plt.xlabel("Log M*")
plt.legend(loc='best')
plt.savefig("FrancoMstarVsSFR.png", bbox_inches='tight', dpi=300)
plt.show()

# sSFR

ssfr, ssfr_error = create_points_and_error("sSFR", full_catalog)
spec_ssfr, spec_ssfr_error = create_points_and_error("sSFR", spec_catalog)
rob_ssfr, rob_ssfr_error = create_points_and_error("sSFR", roberto_catalog)

plt.hist(ssfr, label="All", histtype='step', color='lightgrey')
plt.hist(spec_ssfr, label='Spec', histtype='step', color='yellow')
plt.hist(rob_ssfr, label='Roberto', histtype='step', color='black')
plt.ylabel("Count")
plt.xlabel("Log sSFR")
plt.legend(loc='best')
plt.title("sSFR Distribution")
plt.yscale('log')
plt.savefig("sSFR.png", bbox_inches='tight', dpi=300)
plt.show()

plt.hist(sfr, label="All", histtype='step', color='lightgrey')
plt.hist(spec_sfr, label='Spec', histtype='step', color='yellow')
plt.hist(rob_sfr, label='Roberto', histtype='step', color='black')
plt.ylabel("Count")
plt.xlabel("Log SFR")
plt.legend(loc='best')
plt.title("SFR Distribution")
plt.yscale('log')
plt.savefig("SFR.png", bbox_inches='tight', dpi=300)
plt.show()

plt.hist(mass, label="All", histtype='step', color='lightgrey')
plt.hist(spec_mass, label='Spec', histtype='step', color='yellow')
plt.hist(rob_mass, label='Roberto', histtype='step', color='black')
plt.ylabel("Count")
plt.xlabel("Log Mass [Msun]")
plt.legend(loc='best')
plt.title("Mass Distribution")
plt.yscale('log')
plt.savefig("SFR.png", bbox_inches='tight', dpi=300)
plt.show()

# Print difference between Roberto and Spec Z
spec_ids = spec_catalog['id']
diff_ids = []
for id in roberto_catalog['id']:
    if id not in spec_ids:
        diff_ids.append(id)
print("Fraction of differing galaxies (Roberto): " + str(len(diff_ids) / len(roberto_catalog)))
print("Fraction of differing galaxies (Spec): " + str(len(diff_ids) / len(spec_catalog)))

# Plot distribution of galaxies by redshift
pos_z_mask = full_catalog["z"] > 0.001

z = full_catalog[pos_z_mask]["z"]
spec_z = spec_catalog["z"]
rob_z = roberto_catalog["z"]

plt.hist(z, label="All", histtype='step', color='lightgrey')
plt.hist(spec_z, label='Spec', histtype='step', color='yellow')
plt.hist(rob_z, label='Roberto', histtype='step', color='black')
plt.ylabel("Count")
plt.yscale("log")
plt.legend(loc='best')
plt.xlabel("Redshift")
plt.title("Redshift Distribution")
plt.savefig("RedshiftDistribution.png", bbox_inches='tight', dpi=300)
plt.show()

# Plot the histograms that Leinhard has with the percentages

muse_low_z_sfr, low_z_sfr_error, low_z_sfr_z = create_points_and_error_by_z("SFR", muse_catalog, 0, 1)
muse_mid_z_sfr, mid_z_sfr_error, mid_z_sfr_z = create_points_and_error_by_z("SFR", muse_catalog, 1, 2)
muse_high_z_sfr, high_z_sfr_error, high_z_sfr_z = create_points_and_error_by_z("SFR", muse_catalog, 2, 3)
muse_vhigh_z_sfr, vhigh_z_sfr_error, vhigh_z_sfr_z = create_points_and_error_by_z("SFR", muse_catalog, 3, 4)

muse_low_z_mass, low_z_mass_error, low_z_mass_z = create_points_and_error_by_z("Mstar", muse_catalog, 0, 1)
muse_mid_z_mass, mid_z_mass_error, mid_z_mass_z = create_points_and_error_by_z("Mstar", muse_catalog, 1, 2)
muse_high_z_mass, high_z_mass_error, high_z_mass_z = create_points_and_error_by_z("Mstar", muse_catalog, 2, 3)
muse_vhigh_z_mass, vhigh_z_mass_error, vhigh_z_mass_z = create_points_and_error_by_z("Mstar", muse_catalog, 3, 4)

leinhardt_low_z_sfr, low_z_sfr_error, low_z_sfr_z = create_points_and_error_by_z("SFR", leinhardt_catalog, 0, 1)
leinhardt_mid_z_sfr, mid_z_sfr_error, mid_z_sfr_z = create_points_and_error_by_z("SFR", leinhardt_catalog, 1, 2)
leinhardt_high_z_sfr, high_z_sfr_error, high_z_sfr_z = create_points_and_error_by_z("SFR", leinhardt_catalog, 2, 3)
leinhardt_vhigh_z_sfr, vhigh_z_sfr_error, vhigh_z_sfr_z = create_points_and_error_by_z("SFR", leinhardt_catalog, 3, 4)

leinhardt_low_z_mass, low_z_mass_error, low_z_mass_z = create_points_and_error_by_z("Mstar", leinhardt_catalog, 0, 1)
leinhardt_mid_z_mass, mid_z_mass_error, mid_z_mass_z = create_points_and_error_by_z("Mstar", leinhardt_catalog, 1, 2)
leinhardt_high_z_mass, high_z_mass_error, high_z_mass_z = create_points_and_error_by_z("Mstar", leinhardt_catalog, 2, 3)
leinhardt_vhigh_z_mass, vhigh_z_mass_error, vhigh_z_mass_z = create_points_and_error_by_z("Mstar", leinhardt_catalog, 3, 4)


f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='all')
ax1.hist(low_z_mass, color='lightgrey', histtype='step', bins=10)
ax1.set_yscale('log')
ax1.set_title('0 < Z < 1')
ax1.hist(spec_low_z_mass, color='yellow', histtype='step', bins=10)
ax1.hist(rob_low_z_mass, color='black', histtype='step', bins=10)
ax1.hist(muse_low_z_mass, color='green', histtype='step', bins=10)
ax1.hist(leinhardt_low_z_mass, color='orange', histtype='step', bins=10)
ax2.set_title('0 < Z < 1')
ax2.hist(low_z_sfr, color='lightgrey', histtype='step', bins=10)
ax2.hist(spec_low_z_sfr, color='yellow', histtype='step', bins=10)
ax2.hist(rob_low_z_sfr, color='black', histtype='step', bins=10)
ax2.hist(muse_low_z_sfr, color='green', histtype='step', bins=10)
ax2.hist(leinhardt_low_z_sfr, color='orange', histtype='step', bins=10)
ax3.set_title('1 < Z < 2')
ax3.hist(mid_z_mass, color='lightgrey', histtype='step', bins=10)
ax3.hist(spec_mid_z_mass, color='yellow', histtype='step', bins=10)
ax3.hist(rob_mid_z_mass, color='black', histtype='step', bins=10)
ax3.hist(muse_mid_z_mass, color='green', histtype='step', bins=10)
ax3.hist(leinhardt_mid_z_mass, color='orange', histtype='step', bins=10)
ax4.set_title('1 < Z < 2')
ax4.hist(mid_z_sfr, color='lightgrey', histtype='step', label='All', bins=10)
ax4.hist(spec_mid_z_sfr, color='yellow', histtype='step', label='Spectroscopic', bins=10)
ax4.hist(rob_mid_z_sfr, color='black', histtype='step', label='Roberto', bins=10)
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
ax1.hist(spec_high_z_mass, color='yellow', histtype='step', bins=10)
ax1.hist(rob_high_z_mass, color='black', histtype='step', bins=10)
ax1.hist(muse_high_z_mass, color='green', histtype='step', bins=10)
ax1.hist(leinhardt_high_z_mass, color='orange', histtype='step', bins=10)
ax2.set_title('2 < Z < 3')
ax2.hist(high_z_sfr, color='lightgrey', histtype='step', bins=10)
ax2.hist(spec_high_z_sfr, color='yellow', histtype='step', bins=10)
ax2.hist(rob_high_z_sfr, color='black', histtype='step', bins=10)
ax2.hist(muse_high_z_sfr, color='green', histtype='step', bins=10)
ax2.hist(leinhardt_high_z_sfr, color='orange', histtype='step', bins=10)
ax3.set_title('3 < Z < 4')
ax3.hist(vhigh_z_mass, color='lightgrey', histtype='step', bins=10)
ax3.hist(spec_vhigh_z_mass, color='yellow', histtype='step', bins=10)
ax3.hist(rob_vhigh_z_mass, color='black', histtype='step', bins=10)
ax3.hist(muse_vhigh_z_mass, color='green', histtype='step', bins=10)
ax3.hist(leinhardt_vhigh_z_mass, color='orange', histtype='step', bins=10)
ax4.set_title('3 < Z < 4')
ax4.hist(vhigh_z_sfr, color='lightgrey', histtype='step', label='All', bins=10)
ax4.hist(spec_vhigh_z_sfr, color='yellow', histtype='step', label='Spectroscopic', bins=10)
ax4.hist(rob_vhigh_z_sfr, color='black', histtype='step', label='Roberto', bins=10)
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

# plot MStar vs sSFR

plt.errorbar(mass, ssfr, yerr=ssfr_error, xerr=mass_error, ecolor='lightgrey', fmt='.', mec='lightblue', ms=1,
             label='All', elinewidth=1, zorder=-32)
plt.plot(np.unique(mass), np.poly1d(np.polyfit(mass, ssfr, 1))(np.unique(mass)), color='r')
# plt.errorbar(mass_star, sfr, xerr=mass_star_error, yerr=sfr_error)
plt.title("Log Stellar Mass vs Log Specific Star Formation Rate")
plt.xlabel("Log Stellar Mass (Log(Msun))")
plt.ylabel("Log sSFR")
plt.legend(loc='best')
# plt.xlim(5,np.max(mass_star))
# plt.ylim(np.min(sfr),np.max(sfr))
plt.savefig("MstarVssSFR.png", bbox_inches='tight', dpi=300)
plt.show()

low_z_ssfr, low_z_ssfr_error, low_z_ssfr_z = create_points_and_error_by_z("sSFR", full_catalog, 0, 1)
mid_z_ssfr, mid_z_ssfr_error, mid_z_ssfr_z = create_points_and_error_by_z("sSFR", full_catalog, 1, 2)
high_z_ssfr, high_z_ssfr_error, high_z_ssfr_z = create_points_and_error_by_z("sSFR", full_catalog, 2, 3)
vhigh_z_ssfr, vhigh_z_ssfr_error, vhigh_z_ssfr_z = create_points_and_error_by_z("sSFR", full_catalog, 3, 4)

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='all', sharey='all')
ax1.errorbar(low_z_mass, low_z_ssfr, yerr=low_z_ssfr_error, xerr=low_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='pink', elinewidth=1)
ax1.plot(np.unique(low_z_mass), np.poly1d(np.polyfit(low_z_mass, low_z_ssfr, 1))(np.unique(low_z_mass)),
         label='All fit', color='r', zorder=10)
ax1.set_title('0 < Z < 1')
ax2.errorbar(mid_z_mass, mid_z_ssfr, yerr=mid_z_ssfr_error, xerr=mid_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='blue', elinewidth=1)
ax2.plot(np.unique(mid_z_mass), np.poly1d(np.polyfit(mid_z_mass, mid_z_ssfr, 1))(np.unique(mid_z_mass)), color='r',
         zorder=10)
ax2.set_title('1 < Z < 2')
ax3.errorbar(high_z_mass, high_z_ssfr, yerr=high_z_ssfr_error, xerr=high_z_mass_error, ecolor='lightgrey', fmt='.',
             ms=1, mec='green', elinewidth=1)
ax3.plot(np.unique(high_z_mass), np.poly1d(np.polyfit(high_z_mass, high_z_ssfr, 1))(np.unique(high_z_mass)), color='r',
         zorder=10)
ax3.set_title('2 < Z < 3')
ax4.errorbar(vhigh_z_mass, vhigh_z_ssfr, yerr=vhigh_z_ssfr_error, xerr=vhigh_z_mass_error, ecolor='lightgrey', fmt='.',
             ms=1, mec='orange', elinewidth=1)
ax4.plot(np.unique(vhigh_z_mass), np.poly1d(np.polyfit(vhigh_z_mass, vhigh_z_ssfr, 1))(np.unique(vhigh_z_mass)),
         color='r', zorder=10)
ax4.set_title('3 < Z < 4')
handles, labels = ax1.get_legend_handles_labels()
f.legend(loc='best', handles=handles, labels=labels, prop={'size': 6})
f.text(0.5, 0.01, 'Log Stellar Mass (Mstar)', ha='center')
f.text(0.01, 0.5, 'Log sSFR', va='center', rotation='vertical')
f.savefig("AllMstarVssFRbyZ.png", bbox_inches='tight', dpi=300)
f.show()

# All of them

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='all', sharey='all')
ax1.errorbar(low_z_mass, low_z_ssfr, yerr=low_z_ssfr_error, xerr=low_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='pink', elinewidth=1)
ax1.plot(np.unique(low_z_mass), np.poly1d(np.polyfit(low_z_mass, low_z_ssfr, 1))(np.unique(low_z_mass)),
         label='All fit', color='r', zorder=10)
ax1.set_title('0 < Z < 1')
ax2.errorbar(mid_z_mass, mid_z_ssfr, yerr=mid_z_ssfr_error, xerr=mid_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='blue', elinewidth=1)
ax2.plot(np.unique(mid_z_mass), np.poly1d(np.polyfit(mid_z_mass, mid_z_ssfr, 1))(np.unique(mid_z_mass)), color='r',
         zorder=10)
ax2.set_title('1 < Z < 2')
ax3.errorbar(high_z_mass, high_z_ssfr, yerr=high_z_ssfr_error, xerr=high_z_mass_error, ecolor='lightgrey', fmt='.',
             ms=1, mec='green', elinewidth=1)
ax3.plot(np.unique(high_z_mass), np.poly1d(np.polyfit(high_z_mass, high_z_ssfr, 1))(np.unique(high_z_mass)), color='r',
         zorder=10)
ax3.set_title('2 < Z < 3')
ax4.errorbar(vhigh_z_mass, vhigh_z_ssfr, yerr=vhigh_z_ssfr_error, xerr=vhigh_z_mass_error, ecolor='lightgrey', fmt='.',
             ms=1, mec='orange', elinewidth=1)
ax4.plot(np.unique(vhigh_z_mass), np.poly1d(np.polyfit(vhigh_z_mass, vhigh_z_ssfr, 1))(np.unique(vhigh_z_mass)),
         color='r', zorder=10)
ax4.set_title('3 < Z < 4')

low_z_ssfr, low_z_ssfr_error, low_z_ssfr_z = create_points_and_error_by_z("sSFR", roberto_catalog, 0, 1)
mid_z_ssfr, mid_z_ssfr_error, mid_z_ssfr_z = create_points_and_error_by_z("sSFR", roberto_catalog, 1, 2)
high_z_ssfr, high_z_ssfr_error, high_z_ssfr_z = create_points_and_error_by_z("sSFR", roberto_catalog, 2, 3)
vhigh_z_ssfr, vhigh_z_ssfr_error, vhigh_z_ssfr_z = create_points_and_error_by_z("sSFR", roberto_catalog, 3, 4)

low_z_mass, low_z_mass_error, low_z_mass_z = create_points_and_error_by_z("Mstar", roberto_catalog, 0, 1)
mid_z_mass, mid_z_mass_error, mid_z_mass_z = create_points_and_error_by_z("Mstar", roberto_catalog, 1, 2)
high_z_mass, high_z_mass_error, high_z_mass_z = create_points_and_error_by_z("Mstar", roberto_catalog, 2, 3)
vhigh_z_mass, vhigh_z_mass_error, vhigh_z_mass_z = create_points_and_error_by_z("Mstar", roberto_catalog, 3, 4)

ax1.errorbar(low_z_mass, low_z_ssfr, yerr=low_z_ssfr_error, xerr=low_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='black', label='Roberto', elinewidth=1)
ax1.plot(np.unique(low_z_mass), np.poly1d(np.polyfit(low_z_mass, low_z_ssfr, 1))(np.unique(rob_low_z_mass)),
         label='Roberto fit', color='g', zorder=10)
ax1.set_title('0 < Z < 1')
ax2.errorbar(mid_z_mass, mid_z_ssfr, yerr=mid_z_ssfr_error, xerr=mid_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='black', elinewidth=1)
ax2.plot(np.unique(mid_z_mass), np.poly1d(np.polyfit(mid_z_mass, mid_z_ssfr, 1))(np.unique(mid_z_mass)), color='g',
         zorder=10)
ax2.set_title('1 < Z < 2')
ax3.errorbar(high_z_mass, high_z_ssfr, yerr=high_z_ssfr_error, xerr=high_z_mass_error, ecolor='lightgrey', fmt='.',
             ms=1, mec='black', elinewidth=1)
ax3.plot(np.unique(high_z_mass), np.poly1d(np.polyfit(high_z_mass, high_z_ssfr, 1))(np.unique(high_z_mass)), color='g',
         zorder=10)
ax3.set_title('2 < Z < 3')
ax4.errorbar(vhigh_z_mass, vhigh_z_ssfr, yerr=vhigh_z_ssfr_error, xerr=vhigh_z_mass_error, ecolor='lightgrey', fmt='.',
             ms=1, mec='black', elinewidth=1)
ax4.plot(np.unique(vhigh_z_mass), np.poly1d(np.polyfit(vhigh_z_mass, vhigh_z_ssfr, 1))(np.unique(vhigh_z_mass)),
         color='g', zorder=10)
ax4.set_title('3 < Z < 4')

low_z_ssfr, low_z_ssfr_error, low_z_ssfr_z = create_points_and_error_by_z("sSFR", spec_catalog, 0, 1)
mid_z_ssfr, mid_z_ssfr_error, mid_z_ssfr_z = create_points_and_error_by_z("sSFR", spec_catalog, 1, 2)
high_z_ssfr, high_z_ssfr_error, high_z_ssfr_z = create_points_and_error_by_z("sSFR", spec_catalog, 2, 3)
vhigh_z_ssfr, vhigh_z_ssfr_error, vhigh_z_ssfr_z = create_points_and_error_by_z("sSFR", spec_catalog, 3, 4)

low_z_mass, low_z_mass_error, low_z_mass_z = create_points_and_error_by_z("Mstar", spec_catalog, 0, 1)
mid_z_mass, mid_z_mass_error, mid_z_mass_z = create_points_and_error_by_z("Mstar", spec_catalog, 1, 2)
high_z_mass, high_z_mass_error, high_z_mass_z = create_points_and_error_by_z("Mstar", spec_catalog, 2, 3)
vhigh_z_mass, vhigh_z_mass_error, vhigh_z_mass_z = create_points_and_error_by_z("Mstar", spec_catalog, 3, 4)

ax1.errorbar(low_z_mass, low_z_ssfr, yerr=low_z_ssfr_error, xerr=low_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='yellow', label='Spec', elinewidth=1)
ax1.plot(np.unique(low_z_mass), np.poly1d(np.polyfit(low_z_mass, low_z_ssfr, 1))(np.unique(low_z_mass)),
         label='Spec fit', color='b', zorder=10)
ax1.set_title('0 < Z < 1')
ax2.errorbar(mid_z_mass, mid_z_ssfr, yerr=mid_z_ssfr_error, xerr=mid_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='yellow', elinewidth=1)
ax2.plot(np.unique(mid_z_mass), np.poly1d(np.polyfit(mid_z_mass, mid_z_ssfr, 1))(np.unique(mid_z_mass)), color='b',
         zorder=10)
ax2.set_title('1 < Z < 2')
ax3.errorbar(high_z_mass, high_z_ssfr, yerr=high_z_ssfr_error, xerr=high_z_mass_error, ecolor='lightgrey', fmt='.',
             ms=1, mec='yellow', elinewidth=1)
ax3.plot(np.unique(high_z_mass), np.poly1d(np.polyfit(high_z_mass, high_z_ssfr, 1))(np.unique(high_z_mass)), color='b',
         zorder=10)
ax3.set_title('2 < Z < 3')
ax4.errorbar(vhigh_z_mass, vhigh_z_ssfr, yerr=vhigh_z_ssfr_error, xerr=vhigh_z_mass_error, ecolor='lightgrey', fmt='.',
             ms=1, mec='yellow', elinewidth=1)
ax4.plot(np.unique(vhigh_z_mass), np.poly1d(np.polyfit(vhigh_z_mass, vhigh_z_ssfr, 1))(np.unique(vhigh_z_mass)),
         color='b', zorder=10)
ax4.set_title('3 < Z < 4')

handles, labels = ax1.get_legend_handles_labels()
f.legend(loc='best', handles=handles, labels=labels, prop={'size': 6})
f.text(0.5, 0.01, 'Log Stellar Mass (Mstar)', ha='center')
f.text(0.01, 0.5, 'Log sSFR', va='center', rotation='vertical')
f.savefig("MstarVssFRbyZ.png", bbox_inches='tight', dpi=300)
f.show()

# Now only Roberto


low_z_ssfr, low_z_ssfr_error, low_z_ssfr_z = create_points_and_error_by_z("sSFR", roberto_catalog, 0, 1)
mid_z_ssfr, mid_z_ssfr_error, mid_z_ssfr_z = create_points_and_error_by_z("sSFR", roberto_catalog, 1, 2)
high_z_ssfr, high_z_ssfr_error, high_z_ssfr_z = create_points_and_error_by_z("sSFR", roberto_catalog, 2, 3)
vhigh_z_ssfr, vhigh_z_ssfr_error, vhigh_z_ssfr_z = create_points_and_error_by_z("sSFR", roberto_catalog, 3, 4)

low_z_mass, low_z_mass_error, low_z_mass_z = create_points_and_error_by_z("Mstar", roberto_catalog, 0, 1)
mid_z_mass, mid_z_mass_error, mid_z_mass_z = create_points_and_error_by_z("Mstar", roberto_catalog, 1, 2)
high_z_mass, high_z_mass_error, high_z_mass_z = create_points_and_error_by_z("Mstar", roberto_catalog, 2, 3)
vhigh_z_mass, vhigh_z_mass_error, vhigh_z_mass_z = create_points_and_error_by_z("Mstar", roberto_catalog, 3, 4)

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='all', sharey='all')
ax1.errorbar(low_z_mass, low_z_ssfr, yerr=low_z_ssfr_error, xerr=low_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='pink', label='Roberto', elinewidth=1)
ax1.plot(np.unique(low_z_mass), np.poly1d(np.polyfit(low_z_mass, low_z_ssfr, 1))(np.unique(rob_low_z_mass)),
         label='Roberto fit', color='r', zorder=10)
ax1.set_title('0 < Z < 1')
ax2.errorbar(mid_z_mass, mid_z_ssfr, yerr=mid_z_ssfr_error, xerr=mid_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='blue', elinewidth=1)
ax2.plot(np.unique(mid_z_mass), np.poly1d(np.polyfit(mid_z_mass, mid_z_ssfr, 1))(np.unique(mid_z_mass)), color='r',
         zorder=10)
ax2.set_title('1 < Z < 2')
ax3.errorbar(high_z_mass, high_z_ssfr, yerr=high_z_ssfr_error, xerr=high_z_mass_error, ecolor='lightgrey', fmt='.',
             ms=1, mec='green', elinewidth=1)
ax3.plot(np.unique(high_z_mass), np.poly1d(np.polyfit(high_z_mass, high_z_ssfr, 1))(np.unique(high_z_mass)), color='r',
         zorder=10)
ax3.set_title('2 < Z < 3')
ax4.errorbar(vhigh_z_mass, vhigh_z_ssfr, yerr=vhigh_z_ssfr_error, xerr=vhigh_z_mass_error, ecolor='lightgrey', fmt='.',
             ms=1, mec='orange', elinewidth=1)
ax4.plot(np.unique(vhigh_z_mass), np.poly1d(np.polyfit(vhigh_z_mass, vhigh_z_ssfr, 1))(np.unique(vhigh_z_mass)),
         color='r', zorder=10)
ax4.set_title('3 < Z < 4')
handles, labels = ax1.get_legend_handles_labels()
f.legend(loc='best', handles=handles, labels=labels, prop={'size': 6})
f.text(0.5, 0.01, 'Log Stellar Mass (Mstar)', ha='center')
f.text(0.5, 0.98, 'Roberto Stellar Mass vs sSFR', ha='center')
f.text(0.01, 0.5, 'Log sSFR', va='center', rotation='vertical')
f.savefig("AllMstarVssFRbyZ_Roberto.png", bbox_inches='tight', dpi=300)
f.show()

# Now Only Spectroscopic


low_z_ssfr, low_z_ssfr_error, low_z_ssfr_z = create_points_and_error_by_z("sSFR", spec_catalog, 0, 1)
mid_z_ssfr, mid_z_ssfr_error, mid_z_ssfr_z = create_points_and_error_by_z("sSFR", spec_catalog, 1, 2)
high_z_ssfr, high_z_ssfr_error, high_z_ssfr_z = create_points_and_error_by_z("sSFR", spec_catalog, 2, 3)
vhigh_z_ssfr, vhigh_z_ssfr_error, vhigh_z_ssfr_z = create_points_and_error_by_z("sSFR", spec_catalog, 3, 4)

low_z_mass, low_z_mass_error, low_z_mass_z = create_points_and_error_by_z("Mstar", spec_catalog, 0, 1)
mid_z_mass, mid_z_mass_error, mid_z_mass_z = create_points_and_error_by_z("Mstar", spec_catalog, 1, 2)
high_z_mass, high_z_mass_error, high_z_mass_z = create_points_and_error_by_z("Mstar", spec_catalog, 2, 3)
vhigh_z_mass, vhigh_z_mass_error, vhigh_z_mass_z = create_points_and_error_by_z("Mstar", spec_catalog, 3, 4)

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='all', sharey='all')
ax1.errorbar(low_z_mass, low_z_ssfr, yerr=low_z_ssfr_error, xerr=low_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='pink', label='Spec', elinewidth=1)
ax1.plot(np.unique(low_z_mass), np.poly1d(np.polyfit(low_z_mass, low_z_ssfr, 1))(np.unique(low_z_mass)),
         label='Spec fit', color='r', zorder=10)
ax1.set_title('0 < Z < 1')
ax2.errorbar(mid_z_mass, mid_z_ssfr, yerr=mid_z_ssfr_error, xerr=mid_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='blue', elinewidth=1)
ax2.plot(np.unique(mid_z_mass), np.poly1d(np.polyfit(mid_z_mass, mid_z_ssfr, 1))(np.unique(mid_z_mass)), color='r',
         zorder=10)
ax2.set_title('1 < Z < 2')
ax3.errorbar(high_z_mass, high_z_ssfr, yerr=high_z_ssfr_error, xerr=high_z_mass_error, ecolor='lightgrey', fmt='.',
             ms=1, mec='green', elinewidth=1)
ax3.plot(np.unique(high_z_mass), np.poly1d(np.polyfit(high_z_mass, high_z_ssfr, 1))(np.unique(high_z_mass)), color='r',
         zorder=10)
ax3.set_title('2 < Z < 3')
ax4.errorbar(vhigh_z_mass, vhigh_z_ssfr, yerr=vhigh_z_ssfr_error, xerr=vhigh_z_mass_error, ecolor='lightgrey', fmt='.',
             ms=1, mec='orange', elinewidth=1)
ax4.plot(np.unique(vhigh_z_mass), np.poly1d(np.polyfit(vhigh_z_mass, vhigh_z_ssfr, 1))(np.unique(vhigh_z_mass)),
         color='r', zorder=10)
ax4.set_title('3 < Z < 4')
handles, labels = ax1.get_legend_handles_labels()
f.legend(loc='best', handles=handles, labels=labels, prop={'size': 6})
f.text(0.5, 0.01, 'Log Stellar Mass (Mstar)', ha='center')
f.text(0.5, 0.98, 'Spec Stellar Mass vs sSFR', ha='center')
f.text(0.01, 0.5, 'Log sSFR', va='center', rotation='vertical')
f.savefig("AllMstarVssFRbyZ_Spec.png", bbox_inches='tight', dpi=300)
f.show()

# Plot Av
