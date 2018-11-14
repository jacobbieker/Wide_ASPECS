import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import vstack
from astropy.table import Table, hstack
from scipy import stats
from numpy.polynomial.polynomial import polyfit


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

spec_catalog = full_catalog[spec_z_mask]

F850LP = 2238.1


def janksy_to_AB(ab):
    # Given in microJansky, get Jansky for AB Magnitude
    return 25.0 -2.5 * np.log10((ab))

def convert_flux_to_mag(flux):
    # Flux is in ergs
    return -2.5 * np.log10(flux) + -21.10

count = 0
for row in full_catalog:
    if row['f850lp'] > 0.00001 and row['f160w'] > 0.00001:
        if janksy_to_AB(row['f850lp']) < 27.5 and janksy_to_AB(row['f160w']) < 27.5:
            count += 1
print("Count: ", count)
exit()

# Now the quality cuts
roberto_catalog = perform_cuts(roberto_catalog)
spec_catalog = perform_cuts(spec_catalog)
franco_catalog = perform_cuts(franco_catalog)
full_catalog = perform_cuts(full_catalog)
count = 0
for row in full_catalog:
    if row['f850lp'] > 0.00001 and row['f160w'] > 0.00001:
        if janksy_to_AB(row['f850lp']) < 27.5 and janksy_to_AB(row['f160w']) < 27.5:
            count += 1
print("Count: ", count)

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
             mec='pink', elinewidth=1)
ax1.plot(np.unique(low_z_mass), np.poly1d(np.polyfit(low_z_mass, low_z_sfr, 1))(np.unique(low_z_mass)), label='All fit',
         color='r', zorder=10)
ax1.set_title('0 < Z < 1')
ax2.errorbar(mid_z_mass, mid_z_sfr, yerr=mid_z_sfr_error, xerr=mid_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='blue', elinewidth=1)
ax2.plot(np.unique(mid_z_mass), np.poly1d(np.polyfit(mid_z_mass, mid_z_sfr, 1))(np.unique(mid_z_mass)), color='r',
         zorder=10)
ax2.set_title('1 < Z < 2')
ax3.errorbar(high_z_mass, high_z_sfr, yerr=high_z_sfr_error, xerr=high_z_mass_error, ecolor='lightgrey', fmt='.', ms=1,
             mec='green', elinewidth=1)
ax3.plot(np.unique(high_z_mass), np.poly1d(np.polyfit(high_z_mass, high_z_sfr, 1))(np.unique(high_z_mass)), color='r',
         zorder=10)
ax3.set_title('2 < Z < 3')
ax4.errorbar(vhigh_z_mass, vhigh_z_sfr, yerr=vhigh_z_sfr_error, xerr=vhigh_z_mass_error, ecolor='lightgrey', fmt='.',
             ms=1, mec='orange', elinewidth=1)
ax4.plot(np.unique(vhigh_z_mass), np.poly1d(np.polyfit(vhigh_z_mass, vhigh_z_sfr, 1))(np.unique(vhigh_z_mass)),
         color='r', zorder=10)
# ax4.errorbar(vhigh_z_mass, vhigh_z_sfr, yerr=vhigh_z_sfr_error, xerr=vhigh_z_mass_error)
ax4.set_title('3 < Z < 4')
handles, labels = ax1.get_legend_handles_labels()
f.legend(loc='best', handles=handles, labels=labels, prop={'size': 6})
f.text(0.5, 0.01, 'Log Stellar Mass (Mstar)', ha='center')
f.text(0.01, 0.5, 'Log Star Formation Rate', va='center', rotation='vertical')
f.savefig("AllMstarVsSFRbyZ04.png", bbox_inches='tight', dpi=300)
f.show()

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
# ax1.scatter(franco_low_z_mass, franco_low_z_sfr, color='yellow1, marker='*', s=8)
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
# ax2.scatter(franco_mid_z_mass, franco_mid_z_sfr, color='yellow1, marker='*', s=8)
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
# ax3.scatter(franco_high_z_mass, franco_high_z_sfr, color='yellow1, marker='*', s=8)
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
# ax4.scatter(franco_vhigh_z_mass, franco_vhigh_z_sfr, color='yellow1, marker='*', s=8)
# ax4.errorbar(vhigh_z_mass, vhigh_z_sfr, yerr=vhigh_z_sfr_error, xerr=vhigh_z_mass_error)
ax4.set_title('3 < Z < 4')
handles, labels = ax1.get_legend_handles_labels()
f.legend(loc='best', handles=handles, labels=labels, prop={'size': 6})
f.text(0.5, 0.01, 'Log Stellar Mass (Mstar)', ha='center')
f.text(0.01, 0.5, 'Log Star Formation Rate', va='center', rotation='vertical')
f.savefig("MstarVsSFRbyZ04.png", bbox_inches='tight', dpi=300)
f.show()

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
# plt.scatter(franco_mass, franco_sfr, color='yellow1, marker="*", label='Franco et. al.', s=8)
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

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='all')
ax1.hist(low_z_mass, color='lightgrey', histtype='step', bins=10)
ax1.set_yscale('log')
ax1.set_title('0 < Z < 1')
ax1.hist(spec_low_z_mass, color='yellow', histtype='step', bins=10)
ax1.hist(rob_low_z_mass, color='black', histtype='step', bins=10)
ax2.set_title('0 < Z < 1')
ax2.hist(low_z_sfr, color='lightgrey', histtype='step', bins=10)
ax2.hist(spec_low_z_sfr, color='yellow', histtype='step', bins=10)
ax2.hist(rob_low_z_sfr, color='black', histtype='step', bins=10)
ax3.set_title('1 < Z < 2')
ax3.hist(mid_z_mass, color='lightgrey', histtype='step', bins=10)
ax3.hist(spec_mid_z_mass, color='yellow', histtype='step', bins=10)
ax3.hist(rob_mid_z_mass, color='black', histtype='step', bins=10)
ax4.set_title('1 < Z < 2')
ax4.hist(mid_z_sfr, color='lightgrey', histtype='step', label='All', bins=10)
ax4.hist(spec_mid_z_sfr, color='yellow', histtype='step', label='Spectroscopic', bins=10)
ax4.hist(rob_mid_z_sfr, color='black', histtype='step', label='Roberto', bins=10)
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
ax2.set_title('2 < Z < 3')
ax2.hist(high_z_sfr, color='lightgrey', histtype='step', bins=10)
ax2.hist(spec_high_z_sfr, color='yellow', histtype='step', bins=10)
ax2.hist(rob_high_z_sfr, color='black', histtype='step', bins=10)
ax3.set_title('3 < Z < 4')
ax3.hist(vhigh_z_mass, color='lightgrey', histtype='step', bins=10)
ax3.hist(spec_vhigh_z_mass, color='yellow', histtype='step', bins=10)
ax3.hist(rob_vhigh_z_mass, color='black', histtype='step', bins=10)
ax4.set_title('3 < Z < 4')
ax4.hist(vhigh_z_sfr, color='lightgrey', histtype='step', label='All', bins=10)
ax4.hist(spec_vhigh_z_sfr, color='yellow', histtype='step', label='Spectroscopic', bins=10)
ax4.hist(rob_vhigh_z_sfr, color='black', histtype='step', label='Roberto', bins=10)
handles, labels = ax4.get_legend_handles_labels()
f.legend(loc='best', handles=handles, labels=labels, prop={'size': 6})
f.text(0.25, 0.01, 'Log Stellar Mass (Mstar)', ha='center')
f.text(0.75, 0.01, 'Log SFR (Mstar)', ha='center')
f.text(0.01, 0.25, 'Count', va='center', rotation='vertical')
f.text(0.01, 0.75, 'Count', va='center', rotation='vertical')
f.savefig("LeinhardHist24.png", bbox_inches='tight', dpi=300)
f.show()

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
