import matplotlib.pyplot as plt
from catalog_builder import build_catalog
import numpy as np

# First plot the skymap of the data
full_catalog = build_catalog()

plt.scatter(full_catalog['ra'], full_catalog['dc'], label='All galaxies', s=2)
plt.ylabel("Dec")
plt.xlabel("RA")

# Get the ones with spectroscopic redshift

spec_redshift_cols = ["z_spec_3dh", "zm_vds", "zm_coeS", "zs_mor", "zm_ina", "zm_her"]

# Mask if any of these have a redshift larger than 0.001
spec_z_mask = full_catalog["z_spec_3dh"] > 0.001

plt.scatter(full_catalog[spec_z_mask]['ra'], full_catalog[spec_z_mask]['dc'], label='Spectrographic redshift', s=2)
plt.legend()
plt.show()

# Now colorcoded by redshift

low_z = full_catalog["z"] < 1
mid_z = (full_catalog["z"] > 1) & (full_catalog["z"] < 2)
high_z = (full_catalog["z"] > 2) & (full_catalog["z"] < 3)
highest_z = full_catalog["z"] > 3

plt.scatter(full_catalog[low_z]['ra'], full_catalog[low_z]['dc'], alpha=0.15, label="0 < Z < 1", color='purple', s=1)
plt.scatter(full_catalog[mid_z]['ra'], full_catalog[mid_z]['dc'], alpha=0.15,label="1 < Z < 2", color='blue', s=1)
plt.scatter(full_catalog[high_z]['ra'], full_catalog[high_z]['dc'], alpha=0.15,label="2 < Z < 3", color='green', s=1)
plt.scatter(full_catalog[highest_z]['ra'], full_catalog[highest_z]['dc'],alpha=0.15, label="3 < Z", color='red', s=1)
plt.ylabel('Dec')
plt.xlabel('RA')
plt.title("Galaxies By Redshift")
plt.legend()
plt.show()


plt.scatter(full_catalog[highest_z]['ra'], full_catalog[highest_z]['dc'],alpha=0.15, label="3 < Z", color='red', s=1)
plt.scatter(full_catalog[high_z]['ra'], full_catalog[high_z]['dc'], alpha=0.15,label="2 < Z < 3", color='green',s=1)
plt.scatter(full_catalog[mid_z]['ra'], full_catalog[mid_z]['dc'], alpha=0.15,label="1 < Z < 2", color='blue',s=1)
plt.scatter(full_catalog[low_z]['ra'], full_catalog[low_z]['dc'], alpha=0.15, label="0 < Z < 1", color='purple', s=1)
plt.ylabel('Dec')
plt.xlabel('RA')
plt.title("Galaxies By Redshift (Reversed)")
plt.legend()
plt.show()


# Now stellar mass distribution

plt.hist(full_catalog['lmass'], label='Stellar Mass')
plt.title("Stellar Mass Distribution")
plt.xscale("log")
plt.xlabel("Stellar Mass (Msun)")
plt.ylabel("Count")
plt.xlim(np.log(np.max(full_catalog['lmass'])))
plt.show()

# Star Formation vs Stellar Mass

plt.scatter(full_catalog["lmass"], full_catalog['lsfr'], s=2)
plt.title("Stellar Mass vs Star Formation Rate")
plt.xlabel("Stellar Mass (Msun)")
plt.ylabel("Star Formation Rate")
plt.xscale("log")
plt.yscale("log")
plt.xlim(np.log(np.max(full_catalog['lmass'])))
plt.show()

# Star Formation vs Stellar Mass split by Z

plt.scatter(full_catalog[low_z]["lmass"], full_catalog[low_z]['lsfr'], s=2)
plt.title("Stellar Mass vs Star Formation Rate Z < 1")
plt.xlabel("Log Stellar Mass (Msun)")
plt.ylabel("Log Star Formation Rate")
plt.xscale("log")
plt.yscale("log")
plt.xlim(np.log(np.max(full_catalog[low_z]['lmass'])))
plt.show()

plt.scatter(full_catalog[mid_z]["lmass"], full_catalog[mid_z]['lsfr'], s=2)
plt.title("Stellar Mass vs Star Formation Rate 1 < Z < 2")
plt.xlabel("Log Stellar Mass (Msun)")
plt.ylabel("Log Star Formation Rate")
plt.xscale("log")
plt.yscale("log")
plt.xlim(np.log(np.max(full_catalog[mid_z]['lmass'])))
plt.show()

plt.scatter(full_catalog[high_z]["lmass"], full_catalog[high_z]['lsfr'], s=2)
plt.title("Stellar Mass vs Star Formation Rate 2 < Z < 3")
plt.xlabel("Log Stellar Mass (Msun)")
plt.ylabel("Log Star Formation Rate")
plt.xscale("log")
plt.yscale("log")
plt.xlim(np.log(np.max(full_catalog[high_z]['lmass'])))
plt.show()

plt.scatter(full_catalog[highest_z]["lmass"], full_catalog[highest_z]['lsfr'], s=2)
plt.title("Stellar Mass vs Star Formation Rate 3 < Z")
plt.xlabel("Log Stellar Mass (Msun)")
plt.ylabel("Log Star Formation Rate")
plt.xscale("log")
plt.yscale("log")
plt.xlim(np.log(np.max(full_catalog[highest_z]['lmass'])))
plt.show()

# Distribution of IRAC fluxes

plt.hist(full_catalog['fIRAC1'])
plt.title("fIRAC 3.6mm Distribution")
plt.show()

# Ra and Dec colorcoded by redshift distance

# 1.1mm in RA and Dec

# IRAC fluxes in skymap

