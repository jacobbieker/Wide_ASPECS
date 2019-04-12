import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from spectral_cube import SpectralCube
from astropy.coordinates import match_coordinates_sky, search_around_sky
from astropy.coordinates import SkyCoord, Angle, SkyOffsetFrame, ICRS, Distance
from astropy.table import Table, hstack, join
import astropy.io.fits as fits


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


def angular_distance():
    raise NotImplementedError


def generate_random_catalog(number_of_points, filename):
    """
    Generates random number of points in position of Wide ASPECS  and within the mask



    :param number_of_points:
    :param mask:
    :return:
    """
    cube = SpectralCube.read(filename)
    cube[0,:,:].quicklook() # uses aplpy
    # using wcsaxes
    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.8], projection=cube[0,:,:].wcs)
    ax.imshow(cube.unitless[0,:,:]) # you may need cube[5,:,:].value depending on mpl version
    #fig.show()

    # Have the mask now, so can create from random pixel coordinates, generate SkyCoord

    # Just do it 10 times the number of points, hoping enough are in the valid area
    rand_x = np.random.randint(0, cube.unitless[0,:,:].shape[0], 100*number_of_points)
    rand_y = np.random.randint(0, cube.unitless[0,:,:].shape[0], 100*number_of_points)

    random_catalog_coords = []
    mask = cube.unitless[0,:,:]

    for x,y in zip(rand_x,rand_y):
        if not np.isnan(mask[x,y]):
            random_catalog_coords.append(SkyCoord.from_pixel(x,y,wcs=mask.wcs))
            #print(SkyCoord.from_pixel(x,y,wcs=mask.wcs))
            if len(random_catalog_coords) == number_of_points:
                random_catalog_coords = SkyCoord(random_catalog_coords)
                break

    ax.scatter(random_catalog_coords.ra.hour, random_catalog_coords.dec.degree, c='r')
    fig.show()
    return random_catalog_coords


random_catalog = generate_random_catalog(100, "/home/jacob/Research/Wide_ASPECS/Data/gs_A1_2chn.fits")
random_catalog2 = generate_random_catalog(100, "/home/jacob/Research/Wide_ASPECS/Data/gs_A1_2chn.fits")
real_catalog = load_table("line_search_P3_wa_crop.out")
real_catalog = real_catalog[real_catalog['rsnrrbin'] > 5.8]
real_catalog = make_skycoords(real_catalog, ra='rra', dec='rdc')
print(real_catalog.shape)


def angular_correlation_function(data_catalog, random_catalog):
    """
    Calculates the arrays for the data, random, and data_random for w(theta)
    :param data_catalog:
    :param random_catalog:
    :return:
    """
    distance_bins = np.logspace(0,np.log10(3000), 51)
    # First create the data data one
    data_data = np.zeros(len(data_catalog)*(len(data_catalog)-1))

    for i, element in enumerate(data_catalog):
        for j in range(i+1, len(data_catalog)):
            separation = element.separation(data_catalog[j]).arcsecond
            data_data[i+j-1] = separation
    print(data_data)
    data_data, _ = np.histogram(data_data[np.nonzero(data_data)], bins=distance_bins)

    print("Done with Data Data")

    random_random = np.zeros(len(random_catalog)*(len(random_catalog)-1))
    for i, element in enumerate(random_catalog):
        for j in range(i+1, len(random_catalog)):
            separation = element.separation(random_catalog[j]).arcsecond
            random_random[i+j-1] = separation
    random_random, _ = np.histogram(random_random[np.nonzero(random_random)], bins=distance_bins)

    print("Done with Random Random")

    data_random = np.zeros(len(data_catalog)*len(random_catalog))
    for i, element in enumerate(data_catalog):
        for j, rand_element in enumerate(random_catalog):
            separation = element.separation(rand_element).arcsecond
            data_random[i+j] = separation
    data_random, _ = np.histogram(data_random[np.nonzero(data_random)], bins=distance_bins)

    print("Done with Data Random")

    return data_data, data_random, random_random



def xi_r(data_array, data_random_array, random_array):
    """

    (DD/RR - 2 DR/RR +1)
    :param data_array:
    :param data_random_array:
    :param random_array:
    :return:
    """
    data_array = data_array/np.sum(data_array)
    data_random_array = data_random_array/np.sum(data_random_array)
    random_array = random_array/np.sum(random_array)
    print(sum(random_array))
    print(sum(data_random_array))
    print(sum(data_array))
    return data_array / random_array - ((2 * data_random_array / random_array) + 1)


def xi_r_error(omega_theta, data_array):
    """
    Data array is not normalized
    :param omega_theta:
    :param data_array:
    :return:
    """

    return (1 + omega_theta) / np.sqrt(data_array)


dd, dr, rr = angular_correlation_function(real_catalog, random_catalog)

omega_w = xi_r(dd, dr, rr)
e_omega_w = xi_r_error(omega_w, dd)
distance_bins = np.logspace(0,np.log10(3000), len(omega_w))

plt.cla()
plt.errorbar(x=distance_bins, y=omega_w, yerr=e_omega_w, fmt='o')
plt.xscale("log")
plt.title("Data vs Random")
plt.xlabel("Angular Distance (arcseconds)")
plt.ylabel("$\omega(\\theta)$")
#plt.yscale("log")
plt.tight_layout()
#plt.savefig("Data_vs_Random.png", dpi=300)
plt.show()

dd, dr, rr = angular_correlation_function(random_catalog2, random_catalog)

omega_w = xi_r(dd, dr, rr)
e_omega_w = xi_r_error(omega_w, dd)
distance_bins = np.logspace(0,np.log10(3000), len(omega_w))

plt.cla()
plt.errorbar(x=distance_bins, y=omega_w, yerr=e_omega_w, fmt='o')
plt.xscale("log")
plt.title("Random vs Random")
plt.xlabel("Angular Distance (arcseconds)")
plt.ylabel("$\omega(\\theta)$")
#plt.yscale("log")
plt.tight_layout()
#plt.savefig("Random_vs_Random.png", dpi=300)
plt.show()
