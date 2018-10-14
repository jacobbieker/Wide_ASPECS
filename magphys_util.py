from astropy.io import fits
import numpy as np
from catalog_builder import build_catalog
from astropy.table import Table

hdu_list = fits.open("data/magphys_in.fits")
print(hdu_list.info())

print(hdu_list[1].header)
print(hdu_list[1].data)

data_field = hdu_list[1].data

print(data_field)

list_of_used_ones = []
for element in data_field:
    list_of_used_ones.append(element[0])

print(list_of_used_ones)


full_catalog = build_catalog()
for column in full_catalog.columns:
    if "id" not in column and "ra" not in column and "dc" not in column and "z" not in column and "flag" not in column:
        # Convert to Jansky
        full_catalog[column] *= 0.000001
        continue
    # Now need to convert to Jy from uJy units and write out to FITS file for MAGPHYS

print(full_catalog)
print(data_field[3])




