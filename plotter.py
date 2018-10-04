import matplotlib.pyplot as plt
from catalog_builder import build_catalog

# First plot the skymap of the data
full_catalog = build_catalog()

plt.scatter(full_catalog['ra'], full_catalog['dc'], label='All galaxies')
plt.ylabel("Dec")
plt.xlabel("RA")

# Get the ones with spectroscopic redshift

spec_redshift_cols = ["z_spec_3dh", "zm_vds", "zm_coeS", "zs_mor", "zm_ina", "zm_her"]

# Mask if any of these have a redshift larger than 0.001
spec_z_mask = full_catalog["z_spec_3dh"] > 0.001

plt.scatter(full_catalog[spec_z_mask]['ra'], full_catalog[spec_z_mask]['dc'], label='Spectrographic redshift')
plt.legend()
plt.show()


# Now stellar mass distribution

# Star Formation vs Stellar Mass

# Star Formation vs Stellar Mass split by Z

# Distribution of IRAC fluxes

# Ra and Dec colorcoded by redshift distance

# 1.1mm in RA and Dec

# IRAC fluxes in skymap

