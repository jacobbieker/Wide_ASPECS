import numpy as np
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube
import regions
import matplotlib.pyplot as plt
from astropy.stats import mad_std
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

cubes = ["A1","A2", "B1", "B2", "C1", "C2"]

for cube_name in cubes:
    cube = SpectralCube.read("/home/jacob/Research/Wide_ASPECS/Data/gs_{}_2chn.fits".format(cube_name))

    #region_list = regions.read_ds9('line_search_P3_wa_crop.reg')
    #sub_cube = cube.subcube_from_regions(region_list)
    real_catalog = load_table("line_search_P3_wa_crop.out")
    sn_cut = 9.0
    real_catalog = real_catalog[real_catalog['rsnrrbin'] > sn_cut]
    real_catalog_coords = make_skycoords(real_catalog, ra='rra', dec='rdc')
    real_catalog_freq = real_catalog['rfreq'] * u.GHz
    cube = cube.with_spectral_unit(u.GHz)
    width_of_channel = (cube.spectral_axis[-1] - cube.spectral_axis[0]) / cube.shape[0]

    for index, coord in enumerate(real_catalog):
        lower_ghz_bound = real_catalog_freq[index]-(width_of_channel*((coord['width']-1)/2))
        upper_ghz_bound = real_catalog_freq[index]+(width_of_channel*((coord['width']-1)/2))
        if lower_ghz_bound > cube.spectral_axis[0] and upper_ghz_bound < cube.spectral_axis[-1]:
            print("Cube: {}".format(cube_name))
            print("Source SN: {}".format(coord['rsnrrbin']))
            print("Source Limits: Low: {}\n High: {}".format(lower_ghz_bound, upper_ghz_bound))
            subcube = cube.subcube(zlo=real_catalog_freq[index]-(width_of_channel*((coord['width']-1)/2)),
                                    zhi=real_catalog_freq[index]+(width_of_channel*((coord['width']-1)/2)))
            # Try to get the SN of the cube now...
            subcube.allow_huge_operations=True

            # Now cut down to right around the source
            rms_sub_cube = subcube.subcube(xlo=int(coord['rx']+50-15), xhi=int(coord['rx']+50+15),
                                            ylo=int(coord['ry']+50-15), yhi=int(coord['ry']+50+15))

            sub_cube = subcube.subcube(xlo=int(coord['rx']-15), xhi=int(coord['rx']+15),
                                        ylo=int(coord['ry']-15), yhi=int(coord['ry']+15))

            sub_cube = sub_cube.mean(axis=0)
            try:
                rms_sub_cube = subcube.mean(axis=0)
            except:
                rms_sub_cube = subcube.subcube(xlo=int(coord['rx']-50-15), xhi=int(coord['rx']-50+15),
                                               ylo=int(coord['ry']-50-15), yhi=int(coord['ry']-50+15))
                rms_sub_cube = subcube.mean(axis=0)

            print(sub_cube.shape)
            rms_cube = cube.unitless.unmasked_data[:,:]
            rms_noise = mad_std(sub_cube, ignore_nan=True)
            print(rms_noise)

            rms_cube = sub_cube
            print("RMS Shape {}".format(rms_cube.shape))
            rms_noise = mad_std(sub_cube)
            print("RMS Noise around Source: {}".format(rms_noise))

            rms_rms_cube = rms_sub_cube
            rms_noise = np.nanstd(rms_rms_cube.value)
            print("RMS Noise around Random Space: {}".format(rms_noise))
            try:
                max_val = np.max(sub_cube.value/rms_noise)
                print("Max SN vs Cat SN: {}".format(coord['rsnrrbin']/max_val))

                print("\nSource: {} Width: {}".format(coord['rsnrrbin'], coord['width']))
                print("Flux: {} Peak Flux: {}".format(coord['rflux'], coord['rpeak']))

                # now plot the contours
                class nf(float):
                    def __repr__(self):
                        str = '%.1f' % (self.__float__(),)
                        if str[-1] == '0':
                            return '%.0f' % self.__float__()
                        else:
                            return '%.1f' % self.__float__()

                wcs = sub_cube.wcs
                fig = plt.figure()
                ax = fig.add_axes([0.1,0.1,0.8,0.8], projection=wcs)
                im = ax.imshow(sub_cube.value/rms_noise, origin='lower')
                fig.colorbar(im)
                cs = ax.contour(sub_cube.value/rms_noise, levels=[-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12], colors='white', alpha=0.5)
                cs.levels = [nf(val) for val in cs.levels]
                if plt.rcParams["text.usetex"]:
                    fmt = r'%r'
                else:
                    fmt = '%r'

                ax.clabel(cs, cs.levels, inline=True, fmt=fmt, fontsize=10)
                fig.suptitle("S/N: {} Cube: {} Pix (X,Y): ({},{}) Freq: {} GHz\nCube Range: {} - {}".format(coord['rsnrrbin'], cube_name,
                                                                                   coord['rx'], coord['ry'],
                                                                                   coord['rfreq'],
                                                                                    np.round(real_catalog_freq[index]-(width_of_channel*((coord['width']-1)/2)), 4),
                                                                                    np.round(real_catalog_freq[index]+(width_of_channel*((coord['width']-1)/2)), 4)))
                plt.savefig("/home/jacob/Development/Wide_ASPECS/May_Output/Contours_All_SN_{}_Cube_{}.png".format(coord['rsnrrbin'], cube_name), dpi=300)
            except:
                continue

            #sub_cube[0,:,:].quicklook() # uses aplpy
            # using wcsaxes
            #wcs = sub_cube[0,:,:].wcs
            #fig = plt.figure()
            #ax = fig.add_axes([0.1,0.1,0.8,0.8])
            #ax.imshow(sub_cube.unitless[0,:,:]) # you may need cube[5,:,:].value depending on mpl version
            #fig.show()
