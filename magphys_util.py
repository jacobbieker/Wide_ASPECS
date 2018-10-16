from astropy.io import fits
import numpy as np
from catalog_builder import build_catalog
from astropy.table import Table

hdu_list = fits.open("data/magphys_in.fits")
print(hdu_list.info())
print(hdu_list[1].header)
print(hdu_list[1].columns.info())
print(hdu_list[1].data)

data_field = hdu_list[1].data

list_of_used_ones = []
for element in data_field:
    list_of_used_ones.append(element[0])

full_catalog = build_catalog()

#print(len(full_catalog.columns))
#print(len(data_field[0]))

#print(hdu_list[1].columns.names)
print(full_catalog.columns)
print("In full catalog but not in FITS ---------------------")
for name in full_catalog.columns:
    if name not in hdu_list[1].columns.names:
        print(name)

print("In FITS but not in Full Catalog ---------------------")
for name in hdu_list[1].columns.names:
    if name not in full_catalog.columns:
        print(name)


# Now build the FITS file for use everywhere

fits.writeto("full_catalog_two.fits", np.array(full_catalog))
