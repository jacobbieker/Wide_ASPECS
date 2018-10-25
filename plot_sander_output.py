import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

hdu_list = fits.open("data/cat_fitted_with_MUSE.fits")
print(hdu_list.info())
print(hdu_list[1].header)
print(hdu_list[1].columns.info())
print(hdu_list[1].data)

full_catalog = hdu_list[1].data

# Get Redshifts

pos_z_mask = full_catalog["z"] > 0.001

# Now get the SFR

sfr_cols = ["SFR", "SFR_2.5", "SFR_16", "SFR_50", "SFR_84", "SFR_97.5"]
lgM_Lh_cols = ["lgM_Lh"]
age_cols = ["age_M"]
Ldust_cols = ["Ldust"]
Mdust_cols = ["Mdust"]
fmuSFH_cols = ["fmuSFH"]
Mstar_cols = ["Mstar"]
sSFR_cols = ["sSFR"]
Tdust_cols = ["Tdust"]
