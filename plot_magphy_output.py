import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats
from numpy.polynomial.polynomial import polyfit

def create_points_and_error(column_base_name, full_catalog):
    centerpoints = full_catalog[str(column_base_name + "_50")]
    lower_error = full_catalog[str(column_base_name + "_16")]
    upper_error = full_catalog[str(column_base_name + "_84")]
    centerpoints = np.nan_to_num(centerpoints)
    zero_mask = centerpoints != 0.0
    centerpoints = centerpoints[zero_mask]
    lower_error = lower_error[zero_mask]
    upper_error = upper_error[zero_mask]
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
    lower_error = lower_error[zero_mask]
    upper_error = upper_error[zero_mask]
    z_values = z_values[zero_mask]
    error_bars = [lower_error, upper_error]
    return centerpoints, error_bars, z_values

hdu_list = fits.open("data/jacob_aspecs_catalog_fixed_magphys_jcb3.fits")
print(hdu_list.info())
print(hdu_list[1].header)
print(hdu_list[1].columns.info())
print(hdu_list[1].data)

full_catalog = hdu_list[1].data

# Get Redshifts
spec_redshift_cols = ["z_spec_3dh", "zm_vds", "zm_coeS", "zs_mor", "zm_ina", "zm_her"]
photo_redshift_cols = ["zm_s12", "zm_z13", "zm_m12", "z_m2", "zm_b15", "zm_coe"]

# Mask if any of these have a redshift larger than 0.001
spec_z_mask = (full_catalog["z_spec_3dh"] > 0.001) | (full_catalog["zm_vds"] > 0.001) | (full_catalog["zm_coeS"] > 0.001) \
              | (full_catalog["zs_mor"] > 0.001) | (full_catalog["zm_ina"] > 0.001) | (full_catalog["zm_her"] > 0.001)

photo_z_mask = (full_catalog["zm_s12"] > 0.001) | (full_catalog["zm_z13"] > 0.001) |(full_catalog["zm_m12"] > 0.001) \
               | (full_catalog["z_m2"] > 0.001) | (full_catalog["zm_b15"] > 0.001) | (full_catalog["zm_coe"] > 0.001)

only_photo_z_mask = photo_z_mask & ((full_catalog["z_spec_3dh"] < 0.001) & (full_catalog["zm_vds"] < 0.001) & (full_catalog["zm_coeS"] < 0.001)
                                    & (full_catalog["zs_mor"] < 0.001) & (full_catalog["zm_ina"] < 0.001) & (full_catalog["zm_her"] < 0.001))

only_spec_z_mask = spec_z_mask & ((full_catalog["zm_s12"] < 0.001) & (full_catalog["zm_z13"] < 0.001) &(full_catalog["zm_m12"] < 0.001)
                                  & (full_catalog["z_m2"] < 0.001) & (full_catalog["zm_b15"] < 0.001) & (full_catalog["zm_coe"] < 0.001))

print(len(full_catalog[spec_z_mask]))


pos_z_mask = full_catalog["z"] > 0.001

from copy import deepcopy

all_catalog = deepcopy(full_catalog)

full_catalog = full_catalog[spec_z_mask]

# Quality control cuts
#quality_cuts = (full_catalog["Q0"] < 2.0) & (full_catalog["Q2"] < 1.0) & ((full_catalog["Mstar_50"] - full_catalog["Mstar_16"]) < 0.5) &\
#               ((full_catalog["SFR_50"] - full_catalog["SFR_16"]) < 0.5)

#full_catalog = full_catalog[quality_cuts]


# Now get the SFR

sfr_cols = ["SFR", "SFR_2.5", "SFR_16", "SFR_50", "SFR_84", "SFR_97.5"]
lgM_Lh_cols = ["lgM_Lh"]
age_cols = ["age_M", "age_M_2.5", "age_M_16", "age_M_50", "age_M_84", "age_M_97.5"]
Ldust_cols = ["Ldust"]
Mdust_cols = ["Mdust"]
fmuSFH_cols = ["fmuSFH"]
Mstar_cols = ["Mstar", "Mstar_2.5", "Mstar_16", "Mstar_50", "Mstar_84", "Mstar_97.5"]
sSFR_cols = ["sSFR"]
Tdust_cols = ["Tdust"]

full_catalog = np.nan_to_num(full_catalog)



# Now stellar mass distribution
spec_mass_star, spec_mass_star_error = create_points_and_error("Mstar", full_catalog)
mass_star, mass_star_error = create_points_and_error("Mstar", all_catalog)
plt.hist(mass_star, label='Stellar Mass', bins=100, histtype='step')
plt.hist(spec_mass_star, label="Spec Stellar Mass", bins=100, histtype='step')
plt.title("Stellar Mass Distribution")
plt.xlabel("Log Stellar Mass (Log(Msun))")
plt.ylabel("Count")
plt.legend(loc='best')
plt.show()

# Star Formation vs Stellar Mass

sfr, sfr_error = create_points_and_error("SFR", all_catalog)
spec_sfr, spec_sfr_error = create_points_and_error("SFR", full_catalog)
plt.scatter(mass_star, sfr, s=2, color='lightgrey')
plt.scatter(spec_mass_star, spec_sfr, s=2, color='black')
#plt.errorbar(mass_star, sfr, xerr=mass_star_error, yerr=sfr_error, ecolor='grey', capsize=0, alpha=0.2, fmt='o')
plt.title("Log Stellar Mass vs Log Star Formation Rate")
plt.xlabel("Log Stellar Mass (Log(Msun))")
plt.ylabel("Log Star Formation Rate")
plt.legend(loc='best')
plt.show()

# Star Formation vs Stellar Mass split by Z

low_z_sfr, low_z_sfr_error, low_z_sfr_z = create_points_and_error_by_z("SFR", full_catalog, 0, 1)
mid_z_sfr, mid_z_sfr_error, mid_z_sfr_z = create_points_and_error_by_z("SFR", full_catalog, 1, 2)
high_z_sfr, high_z_sfr_error, high_z_sfr_z = create_points_and_error_by_z("SFR", full_catalog, 2, 3)
vhigh_z_sfr, vhigh_z_sfr_error, vhigh_z_sfr_z = create_points_and_error_by_z("SFR", full_catalog, 3, 4)

low_z_mass, low_z_mass_error, low_z_mass_z = create_points_and_error_by_z("Mstar", full_catalog, 0, 1)
mid_z_mass, mid_z_mass_error, mid_z_mass_z = create_points_and_error_by_z("Mstar", full_catalog, 1, 2)
high_z_mass, high_z_mass_error, high_z_mass_z = create_points_and_error_by_z("Mstar", full_catalog, 2, 3)
vhigh_z_mass, vhigh_z_mass_error, vhigh_z_mass_z = create_points_and_error_by_z("Mstar", full_catalog, 3, 4)

all_low_z_sfr, low_z_sfr_error, low_z_sfr_z = create_points_and_error_by_z("SFR", all_catalog, 0, 1)
all_mid_z_sfr, mid_z_sfr_error, mid_z_sfr_z = create_points_and_error_by_z("SFR", all_catalog, 1, 2)
all_high_z_sfr, high_z_sfr_error, high_z_sfr_z = create_points_and_error_by_z("SFR", all_catalog, 2, 3)
all_vhigh_z_sfr, vhigh_z_sfr_error, vhigh_z_sfr_z = create_points_and_error_by_z("SFR", all_catalog, 3, 4)

all_low_z_mass, low_z_mass_error, low_z_mass_z = create_points_and_error_by_z("Mstar", all_catalog, 0, 1)
all_mid_z_mass, mid_z_mass_error, mid_z_mass_z = create_points_and_error_by_z("Mstar", all_catalog, 1, 2)
all_high_z_mass, high_z_mass_error, high_z_mass_z = create_points_and_error_by_z("Mstar", all_catalog, 2, 3)
all_vhigh_z_mass, vhigh_z_mass_error, vhigh_z_mass_z = create_points_and_error_by_z("Mstar", all_catalog, 3, 4)


f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='all', sharey='all')
ax1.scatter(all_low_z_mass, all_low_z_sfr, color='pink', s=1)
ax1.scatter(low_z_mass, low_z_sfr, color='black', s=1)
ax1.plot(np.unique(low_z_mass), np.poly1d(np.polyfit(low_z_mass, low_z_sfr, 1))(np.unique(low_z_mass)), color='r')
ax1.plot(np.unique(all_low_z_mass), np.poly1d(np.polyfit(all_low_z_mass, all_low_z_sfr, 1))(np.unique(all_low_z_mass)), color='lightgrey')
ax1.set_title('0 < Z < 1')
ax2.scatter(all_mid_z_mass, all_mid_z_sfr, color='blue', s=1)
ax2.scatter(mid_z_mass, mid_z_sfr, color='black', s=1)
ax2.plot(np.unique(mid_z_mass), np.poly1d(np.polyfit(mid_z_mass, mid_z_sfr, 1))(np.unique(mid_z_mass)), color='r')
ax2.plot(np.unique(all_mid_z_mass), np.poly1d(np.polyfit(all_mid_z_mass, all_mid_z_sfr, 1))(np.unique(all_mid_z_mass)), color='lightgrey')
ax2.set_title('1 < Z < 2')
ax3.scatter(all_high_z_mass, all_high_z_sfr, color='green', s=1)
ax3.scatter(high_z_mass, high_z_sfr, color='black', s=1)
ax3.plot(np.unique(high_z_mass), np.poly1d(np.polyfit(high_z_mass, high_z_sfr, 1))(np.unique(high_z_mass)), color='r')
ax3.plot(np.unique(all_high_z_mass), np.poly1d(np.polyfit(all_high_z_mass, all_high_z_sfr, 1))(np.unique(all_high_z_mass)), color='lightgrey')
ax3.set_title('2 < Z < 3')
#ax3.set_ylim(-2, 4)
ax4.scatter(all_vhigh_z_mass, all_vhigh_z_sfr, color='orange', s=1)
ax4.scatter(vhigh_z_mass, vhigh_z_sfr, color='black', s=1)
ax4.plot(np.unique(vhigh_z_mass), np.poly1d(np.polyfit(vhigh_z_mass, vhigh_z_sfr, 1))(np.unique(vhigh_z_mass)), color='r')
ax4.plot(np.unique(all_vhigh_z_mass), np.poly1d(np.polyfit(all_vhigh_z_mass, all_vhigh_z_sfr, 1))(np.unique(all_vhigh_z_mass)), color='lightgrey')
#ax4.errorbar(vhigh_z_mass, vhigh_z_sfr, yerr=vhigh_z_sfr_error, xerr=vhigh_z_mass_error, ecolor='grey', capsize=0, alpha=0.2, fmt='o')
ax4.set_title('3 < Z < 4')
f.text(0.5, 0.01, 'Log Stellar Mass (M*)', ha='center')
f.text(0.01, 0.5, 'Log Star Formation Rate', va='center', rotation='vertical')
f.show()


# Specific Star Formation Rate

ssfr, ssfr_error = create_points_and_error("sSFR", full_catalog)

plt.hist(ssfr, label="sSFR", bins=100)
plt.ylabel("Count")
plt.xlabel("Log sSFR")
plt.title("sSFR Distribution")
plt.show()

plt.scatter(spec_mass_star, ssfr, s=2)
plt.plot(np.unique(spec_mass_star), np.poly1d(np.polyfit(spec_mass_star, ssfr, 1))(np.unique(spec_mass_star)), color='r')
#plt.errorbar(mass_star, sfr, xerr=mass_star_error, yerr=sfr_error)
plt.title("Log Stellar Mass vs Log Specific Star Formation Rate")
plt.xlabel("Log Stellar Mass (Log(Msun))")
plt.ylabel("Log sSFR")
#plt.xlim(5,np.max(mass_star))
#plt.ylim(np.min(sfr),np.max(sfr))
plt.show()

# Now do it by redshift

low_z_ssfr, low_z_ssfr_error, low_z_ssfr_z = create_points_and_error_by_z("sSFR", full_catalog, 0, 1)
mid_z_ssfr, mid_z_ssfr_error, mid_z_ssfr_z = create_points_and_error_by_z("sSFR", full_catalog, 1, 2)
high_z_ssfr, high_z_ssfr_error, high_z_ssfr_z = create_points_and_error_by_z("sSFR", full_catalog, 2, 3)
vhigh_z_ssfr, vhigh_z_ssfr_error, vhigh_z_ssfr_z = create_points_and_error_by_z("sSFR", full_catalog, 3, 4)


f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='all', sharey='all')
ax1.scatter(low_z_mass, low_z_ssfr, color='pink', s=1)
ax1.plot(np.unique(low_z_mass), np.poly1d(np.polyfit(low_z_mass, low_z_ssfr, 1))(np.unique(low_z_mass)), color='r')
ax1.set_title('0 < Z < 1')
ax2.scatter(mid_z_mass, mid_z_ssfr, color='blue', s=1)
ax2.plot(np.unique(mid_z_mass), np.poly1d(np.polyfit(mid_z_mass, mid_z_ssfr, 1))(np.unique(mid_z_mass)), color='r')
ax2.set_title('1 < Z < 2')
ax3.scatter(high_z_mass, high_z_ssfr, color='green', s=1)
ax3.plot(np.unique(high_z_mass), np.poly1d(np.polyfit(high_z_mass, high_z_ssfr, 1))(np.unique(high_z_mass)), color='r')
ax3.set_title('2 < Z < 3')
ax4.scatter(vhigh_z_mass, vhigh_z_ssfr, color='orange', s=1)
ax4.plot(np.unique(vhigh_z_mass), np.poly1d(np.polyfit(vhigh_z_mass, vhigh_z_ssfr, 1))(np.unique(vhigh_z_mass)), color='r')
#ax4.errorbar(vhigh_z_mass, vhigh_z_sfr, yerr=vhigh_z_sfr_error, xerr=vhigh_z_mass_error)
ax4.set_title('3 < Z < 4')
f.text(0.5, 0.01, 'Log Stellar Mass (M*)', ha='center')
f.text(0.01, 0.5, 'Log Specific Star Formation Rate', va='center', rotation='vertical')
f.show()

# Plot the Age distribution
age, age_error = create_points_and_error("age_M", full_catalog)

plt.hist(age, label="Log(Age)", bins=100)
plt.ylabel("Count")
plt.xlabel("Log(Age)")
plt.title("Age Distribution")
plt.show()

# Plot the Dust distribution and Temp
dust_mass, dust_mass_error = create_points_and_error("Mdust", full_catalog)
plt.hist(dust_mass, label="Log Dust Mass", bins=100)
plt.ylabel("Count")
plt.xlabel("Log Dust Mass")
plt.title("Dust Mass Distribution")
plt.show()

dust_temp, dust_temp_error = create_points_and_error("Tdust", full_catalog)
plt.hist(dust_temp, label="Log Dust Temperature", bins=100)
plt.ylabel("Count")
plt.xlabel("Log Dust Mass")
plt.title("Dust Temperature Distribution")
plt.show()

# Dust Mass vs Dust Temp
plt.scatter(dust_mass, dust_temp, s=1)
plt.plot(np.unique(dust_mass), np.poly1d(np.polyfit(dust_mass, dust_temp, 1))(np.unique(dust_mass)), color='r')
#plt.errorbar(dust_mass, dust_temp, xerr=dust_mass_error, yerr=dust_temp_error)
plt.ylabel("Log Dust Temperature")
plt.xlabel("Log Dust Mass")
plt.title("Dust Mass vs Temperature")
plt.show()

# Now dust Luminosity vs Dust temp, should show something

dust_lum, dust_lum_error = create_points_and_error("Ldust", full_catalog)

plt.scatter(dust_mass, dust_lum, s=1)
plt.plot(np.unique(dust_mass), np.poly1d(np.polyfit(dust_mass, dust_lum, 1))(np.unique(dust_mass)), color='r')
#plt.errorbar(dust_mass, dust_lum, xerr=dust_mass_error, yerr=dust_lum_error)
plt.ylabel("Log Dust Luminosity")
plt.xlabel("Log Dust Mass")
plt.title("Dust Mass vs Luminosity")
plt.show()

# Now from Leinhard, Log MStar and Log SFR number and percentage of the total with MUSE vs Non MUSE
all_mstar_12, vhigh_z_ssfr_error, vhigh_z_ssfr_z = create_points_and_error_by_z("MStar", all_catalog, 1, 2)
all_sfr_12, vhigh_z_ssfr_error, vhigh_z_ssfr_z = create_points_and_error_by_z("SFR", all_catalog, 1, 2)

full_mstar_12, vhigh_z_ssfr_error, vhigh_z_ssfr_z = create_points_and_error_by_z("MStar", full_catalog, 1, 2)
full_sfr_12, vhigh_z_ssfr_error, vhigh_z_ssfr_z = create_points_and_error_by_z("SFR", full_catalog, 1, 2)

all_mstar_23, vhigh_z_ssfr_error, vhigh_z_ssfr_z = create_points_and_error_by_z("MStar", all_catalog, 2, 3)
all_sfr_23, vhigh_z_ssfr_error, vhigh_z_ssfr_z = create_points_and_error_by_z("SFR", all_catalog, 2, 3)

full_mstar_23, vhigh_z_ssfr_error, vhigh_z_ssfr_z = create_points_and_error_by_z("MStar", full_catalog, 2, 3)
full_sfr_23, vhigh_z_ssfr_error, vhigh_z_ssfr_z = create_points_and_error_by_z("SFR", full_catalog, 2, 3)

all_mstar_34, vhigh_z_ssfr_error, vhigh_z_ssfr_z = create_points_and_error_by_z("MStar", all_catalog, 3, 4)
all_sfr_34, vhigh_z_ssfr_error, vhigh_z_ssfr_z = create_points_and_error_by_z("SFR", all_catalog, 3, 4)

full_mstar_34, vhigh_z_ssfr_error, vhigh_z_ssfr_z = create_points_and_error_by_z("MStar", full_catalog, 3, 4)
full_sfr_34, vhigh_z_ssfr_error, vhigh_z_ssfr_z = create_points_and_error_by_z("SFR", full_catalog, 3, 4)

f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex='col', sharey='all')
ax1.hist(all_mstar_12, color='black', histtype='step', bins=10)
ax1.set_yscale('log')
ax1.set_title('1 < Z < 2 LogMStar')
ax1.hist(full_mstar_12, color='blue', histtype='step', bins=10)
ax2.set_title('1 < Z < 2 LogSFR')
ax2.hist(all_sfr_12, color='black', histtype='step', bins=10)
ax2.hist(full_sfr_12, color='blue', histtype='step', bins=10)
ax3.set_title('2 < Z < 3 LogMStar')
ax3.hist(all_mstar_23, color='black', histtype='step', bins=10)
ax3.hist(full_mstar_23, color='blue', histtype='step', bins=10)
ax4.set_title('2 < Z < 3 LogSFR')
ax4.hist(all_sfr_23, color='black', histtype='step', bins=10)
ax4.hist(full_sfr_23, color='blue', histtype='step', bins=10)

ax6.set_title('3 < Z < 4 LogSFR')
ax6.hist(all_sfr_34, color='black', histtype='step', bins=10)
ax6.hist(full_sfr_34, color='blue', histtype='step', bins=10)

ax5.set_title('3 < Z < 4 LogMStar')
ax5.hist(all_mstar_34, color='black', histtype='step', bins=10)
ax5.hist(full_mstar_34, color='blue', histtype='step', bins=10)

f.show()

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='all')
ax1.hist(all_mstar_12, color='black', histtype='step', bins=10)
ax1.set_yscale('log')
ax1.set_title('1 < Z < 2 LogMStar')
ax1.set_ylabel('Count')
ax1.set_xlabel('Log M*')
ax1.hist(full_mstar_12, color='blue', histtype='step', bins=10)
ax2.set_title('1 < Z < 2 LogSFR')
ax2.set_ylabel('Count')
ax2.set_xlabel('Log SFR')
ax2.hist(all_sfr_12, color='black', histtype='step', bins=10)
ax2.hist(full_sfr_12, color='blue', histtype='step', bins=10)
ax3.set_title('2 < Z < 3 LogMStar')
ax3.set_ylabel('Count')
ax3.set_xlabel('Log M*')
ax3.hist(all_mstar_23, color='black', histtype='step', bins=10)
ax3.hist(full_mstar_23, color='blue', histtype='step', bins=10)
ax4.set_title('2 < Z < 3 LogSFR')
ax4.set_ylabel('Count')
ax4.set_xlabel('Log SFR')
ax4.hist(all_sfr_23, color='black', histtype='step', bins=10)
ax4.hist(full_sfr_23, color='blue', histtype='step', bins=10)

f.show()

f, ((ax1, ax2)) = plt.subplots(1, 2, sharex='col', sharey='all')
ax1.hist(all_mstar_12, color='black', histtype='step', bins=10)
ax1.set_yscale('log')
ax1.set_title('1 < Z < 2 LogMStar')
ax1.set_ylabel('Count')
ax1.set_xlabel('Log M*')
ax1.hist(full_mstar_12, color='blue', histtype='step', bins=10)
ax2.set_title('1 < Z < 2 LogSFR')
ax2.set_ylabel('Count')
ax2.set_xlabel('Log SFR')
ax2.hist(all_sfr_12, color='black', histtype='step', bins=10)
ax2.hist(full_sfr_12, color='blue', histtype='step', bins=10)
f.show()

f, ((ax3, ax4)) = plt.subplots(1, 2, sharex='col', sharey='all')
ax3.set_yscale('log')
ax3.set_title('2 < Z < 3 LogMStar')
ax3.set_ylabel('Count')
ax3.set_xlabel('Log M*')
ax3.hist(all_mstar_23, color='black', histtype='step', bins=10)
ax3.hist(full_mstar_23, color='blue', histtype='step', bins=10)
ax4.set_title('2 < Z < 3 LogSFR')
ax4.set_ylabel('Count')
ax4.set_xlabel('Log SFR')
ax4.hist(all_sfr_23, color='black', histtype='step', bins=10)
ax4.hist(full_sfr_23, color='blue', histtype='step', bins=10)

f.show()

# Now try the tV