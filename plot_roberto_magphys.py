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

# List of the fits files to open
roberto_magphys_output = ["/home/jacob/Research/MAGPHYS/roberto2/catalog.fits",
                          "/home/jacob/Research/MAGPHYS/roberto3/catalog.fits", "/home/jacob/Research/MAGPHYS/roberto4/catalog.fits",
                          "/home/jacob/Research/MAGPHYS/roberto5/catalog.fits", "/home/jacob/Research/MAGPHYS/roberto6/catalog.fits",
                          "/home/jacob/Research/MAGPHYS/roberto7/catalog.fits", "/home/jacob/Research/MAGPHYS/roberto8/catalog.fits",
                          "/home/jacob/Research/MAGPHYS/roberto9/catalog.fits"]

full_catalog = Table.read("/home/jacob/Research/MAGPHYS/roberto1/catalog.fits", format='fits')

for filename in roberto_magphys_output:
    hdu_list = Table.read(filename, format='fits')
    full_catalog = vstack([full_catalog, hdu_list])

print(full_catalog)
# Get Redshifts
spec_redshift_cols = ["z_spec_3dh", "zm_vds", "zm_coeS", "zs_mor", "zm_ina", "zm_her"]
photo_redshift_cols = ["zm_s12", "zm_z13", "zm_m12", "z_m2", "zm_b15", "zm_coe"]

# Mask if any of these have a redshift larger than 0.001

pos_z_mask = full_catalog["z"] > 0.001
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
mass_star, mass_star_error = create_points_and_error("Mstar", full_catalog)
plt.hist(mass_star, label='Stellar Mass', bins=100)
plt.title("Stellar Mass Distribution")
plt.xlabel("Log Stellar Mass (Log(Msun))")
plt.ylabel("Count")
plt.show()

# Star Formation vs Stellar Mass

sfr, sfr_error = create_points_and_error("SFR", full_catalog)
plt.scatter(mass_star, sfr, s=2)
plt.errorbar(mass_star, sfr, xerr=mass_star_error, yerr=sfr_error)
plt.title("Log Stellar Mass vs Log Star Formation Rate")
plt.xlabel("Log Stellar Mass (Log(Msun))")
plt.ylabel("Log Star Formation Rate")
#plt.xlim(5,np.max(mass_star))
#plt.ylim(np.min(sfr),np.max(sfr))
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


f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
ax1.scatter(low_z_mass, low_z_sfr, color='pink', s=1)
ax1.plot(np.unique(low_z_mass), np.poly1d(np.polyfit(low_z_mass, low_z_sfr, 1))(np.unique(low_z_mass)), color='r')
ax1.set_title('0 < Z < 1')
ax2.scatter(mid_z_mass, mid_z_sfr, color='blue', s=1)
ax2.plot(np.unique(mid_z_mass), np.poly1d(np.polyfit(mid_z_mass, mid_z_sfr, 1))(np.unique(mid_z_mass)), color='r')
ax2.set_title('1 < Z < 2')
ax3.scatter(high_z_mass, high_z_sfr, color='green', s=1)
ax3.plot(np.unique(high_z_mass), np.poly1d(np.polyfit(high_z_mass, high_z_sfr, 1))(np.unique(high_z_mass)), color='r')
ax3.set_title('2 < Z < 3')
ax4.scatter(vhigh_z_mass, vhigh_z_sfr, color='orange', s=1)
ax4.plot(np.unique(vhigh_z_mass), np.poly1d(np.polyfit(vhigh_z_mass, vhigh_z_sfr, 1))(np.unique(vhigh_z_mass)), color='r')
#ax4.errorbar(vhigh_z_mass, vhigh_z_sfr, yerr=vhigh_z_sfr_error, xerr=vhigh_z_mass_error)
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

plt.scatter(mass_star, ssfr, s=2)
plt.plot(np.unique(mass_star), np.poly1d(np.polyfit(mass_star, ssfr, 1))(np.unique(mass_star)), color='r')
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


f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
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

