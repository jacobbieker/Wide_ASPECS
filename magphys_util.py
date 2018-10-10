from astropy.io import fits
import numpy as np


hdu_list = fits.open("data/magphys_in.fits")
print(hdu_list.info())

print(hdu_list[1].header)
print(hdu_list[1].data)

