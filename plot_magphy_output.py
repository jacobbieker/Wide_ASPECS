import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats

def create_points_and_error(column_base_name, full_catalog):
    centerpoints = full_catalog[str(column_base_name + "_50")]
    lower_error = full_catalog[str(column_base_name + "_16")]
    upper_error = full_catalog[str(column_base_name + "_84")]
    centerpoints = np.nan_to_num(centerpoints)
    lower_error = np.nan_to_num(lower_error)
    upper_error = np.nan_to_num(upper_error)
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

hdu_list = fits.open("data/jacob_aspecs_catalog_fixed_magphys_jcb.fits")
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

pos_z_mask = full_catalog["z"] > 0.001

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
x = range(len(mass_star))
plt.scatter(x, mass_star)
plt.errorbar(x, mass_star, yerr=mass_star_error)
plt.show()
plt.hist(mass_star, label='Stellar Mass', log=True, bins=100)
plt.title("Stellar Mass Distribution")
plt.xlabel("Log Stellar Mass (Log(Msun))")
plt.ylabel("Count")
plt.show()

# Star Formation vs Stellar Mass

sfr, sfr_error = create_points_and_error("SFR", full_catalog)
plt.scatter(mass_star, sfr, s=2)
#plt.errorbar(mass_star, sfr, xerr=mass_star_error, yerr=sfr_error)
plt.title("Log Stellar Mass vs Log Star Formation Rate")
plt.xlabel("Log Stellar Mass (Log(Msun))")
plt.ylabel("Log Star Formation Rate")
plt.xlim(5,np.max(mass_star))
plt.ylim(np.min(sfr),np.max(sfr))
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
ax1.set_title('0 < Z < 1')
ax2.scatter(mid_z_mass, mid_z_sfr, color='blue', s=1)
ax2.set_title('1 < Z < 2')
ax3.scatter(high_z_mass, high_z_sfr, color='green', s=1)
ax3.set_title('2 < Z < 3')
ax4.scatter(vhigh_z_mass, vhigh_z_sfr, color='orange', s=1)
#ax4.errorbar(vhigh_z_mass, vhigh_z_sfr, yerr=vhigh_z_sfr_error, xerr=vhigh_z_mass_error)
ax4.set_title('3 < Z < 4')
f.text(0.5, 0.01, 'Log Stellar Mass (Mstar)', ha='center')
f.text(0.01, 0.5, 'Log Star Formation Rate', va='center', rotation='vertical')
f.show()