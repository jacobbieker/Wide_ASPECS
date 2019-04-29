import numpy as np
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube
import regions
import matplotlib.pyplot as plt
import astropy.units as u
from spectral_cube import SpectralCube
from astropy.coordinates import SkyCoord, Angle, SkyOffsetFrame, ICRS, Distance
from astropy.table import Table, hstack, join


def load_table(ascii_table, header=0, start=1):
    ascii_table_data = Table.read(ascii_table, format="ascii", header_start=header, data_start=start)
    return ascii_table_data


def make_skycoords(source, ra='ra', dec='dec', distance=None):
    """
    Makes and returns a SkyCoord array from given source
    :param source: Source with information
    :param ra: Key for RA
    :param dec: Key for Dec
    :return: SkyCoord list
    """
    try:
        if distance is None:
            skycoords = SkyCoord(source[ra] * u.deg, source[dec] * u.deg, frame='icrs')
        else:
            distances = Distance(z=source[distance])
            skycoords = SkyCoord(source[ra] * u.deg, source[dec] * u.deg, distance=distances, frame='icrs')
    except:
        if distance is None:
            skycoords = SkyCoord(source[ra], source[dec], unit=(u.hour, u.deg), frame='icrs')
        else:
            distances = Distance(z=source[distance])
            skycoords = SkyCoord(source[ra], source[dec], unit=(u.hour, u.deg), distance=distances, frame='icrs')

    return skycoords


"""

Steps:

1. Get the Table with units of the SN, rx, ry, potentially SkyCoord, width, and SN

3. Get the slice with the closest to rounded GHz value to the number

4. Take (width - 1)/2. points in either direction

5. For each point, get the flux

6. Try with adding SN of each slice, and then as SN of all the included channels to see if it gets
close

7. Check if the summed up flux matches what is in rflux as well



"""

cubes = ["A1", "A2", "B1", "B2", "C1", "C2"]

for cube_name in cubes:
    cube = SpectralCube.read("/home/jacob/Research/Wide_ASPECS/Data/gs_{}_2chn.fits".format(cube_name))

    region_list = regions.read_ds9('line_search_P3_wa_crop.reg')
    #sub_cube = cube.subcube_from_regions(region_list)
    real_catalog = load_table("line_search_P3_wa_crop.out")
    sn_cut = 9.
    real_catalog = real_catalog[real_catalog['rsnrrbin'] > sn_cut]
    real_catalog_coords = make_skycoords(real_catalog, ra='rra', dec='rdc')
    real_catalog_freq = real_catalog['rfreq'] * u.GHz
    cube = cube.with_spectral_unit(u.GHz)
    width_of_channel = (cube.spectral_axis[-1] - cube.spectral_axis[0]) / cube.shape[0]

    for index, coord in enumerate(real_catalog):
        lower_ghz_bound = real_catalog_freq[index]-(width_of_channel*((coord['width']-1)/2))
        upper_ghz_bound = real_catalog_freq[index]+(width_of_channel*((coord['width']-1)/2))
        if lower_ghz_bound > cube.spectral_axis[0] and upper_ghz_bound < cube.spectral_axis[-1]:
            if coord['width'] < 5:
                print("Cube: {}".format(cube_name))
                sub_cube = cube.subcube(zlo=real_catalog_freq[index]-(width_of_channel*((coord['width']-1)/2)),
                                        zhi=real_catalog_freq[index]+(width_of_channel*((coord['width']-1)/2)))
                # Try to get the SN of the cube now...
                rms_values = []
                for i in range(sub_cube.shape[0]):
                    rms_cube = cube.unitless.unmasked_data[i,:,:]
                    rms_noise = np.nanstd(sub_cube)
                    rms_values.append(rms_noise)

                # Now cut down to right around the source
                sub_cube = sub_cube.subcube(xlo=int(coord['rx']-1), xhi=int(coord['rx']),
                                            ylo=int(coord['ry']-1), yhi=int(coord['ry']))

                sub_vals = []
                print(sub_cube.shape)
                for i in range(sub_cube.shape[0]):
                    sub_val_cube = sub_cube[i,:,:]
                    #print(sub_val_cube)
                    sub_cube_values = sub_val_cube[~np.isnan(sub_val_cube)]
                    sub_cube_values /= rms_values[i]
                    print(sub_cube_values)
                    sub_vals.append(sub_cube_values)
                sub_vals = np.asarray(sub_vals)
                print("\nSource: {} Width: {}".format(coord['rsnrrbin'], coord['width']))
                print("Flux: {} Peak Flux: {}".format(coord['rflux'], coord['rpeak']))
                print("Sum Axis = 1: {}".format(np.sum(sub_vals, axis=1)))
                print("Sum Axis = 0: {}".format(np.sum(sub_vals, axis=0)))
                print("Sum: {}\n".format(np.sum(sub_vals)))

                #sub_cube[0,:,:].quicklook() # uses aplpy
                # using wcsaxes
                #wcs = sub_cube[0,:,:].wcs
                #fig = plt.figure()
                #ax = fig.add_axes([0.1,0.1,0.8,0.8])
                #ax.imshow(sub_cube.unitless[0,:,:]) # you may need cube[5,:,:].value depending on mpl version
                #fig.show()
