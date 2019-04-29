import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from spectral_cube import SpectralCube
from astropy.coordinates import SkyCoord, Angle, SkyOffsetFrame, ICRS, Distance
from astropy.table import Table, hstack, join
from astropy.stats import histogram


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
    wcs = cube[0,:,:].wcs
    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.8], projection=wcs)
    #ax.imshow(cube.unitless[0,:,:]) # you may need cube[5,:,:].value depending on mpl version
    #fig.show()

    # Have the mask now, so can create from random pixel coordinates, generate SkyCoord

    # Just do it 10 times the number of points, hoping enough are in the valid area
    rand_x = np.random.randint(0, cube.unitless[0,:,:].shape[0], 100*number_of_points)
    rand_y = np.random.randint(0, cube.unitless[0,:,:].shape[1], 100*number_of_points)


    random_catalog_coords = []
    r_c_x = []
    r_c_y = []
    mask = cube.unitless[0,:,:]

    for x,y in zip(rand_x,rand_y):
        if not np.isnan(mask[x,y]):
            random_catalog_coords.append(SkyCoord.from_pixel(x,y,wcs=wcs, origin=1))
            r_c_x.append(x)
            r_c_y.append(y)
            #if len(r_c_x) == number_of_points:
            #    break
            #print(SkyCoord.from_pixel(x,y,wcs=mask.wcs))
            if len(random_catalog_coords) == number_of_points:
                random_catalog_coords = SkyCoord(random_catalog_coords)
                break

    #random_catalog_coords = SkyCoord(ras_sim, decs_sim)

    values = random_catalog_coords.to_pixel(wcs=wcs)
    r_c_x, r_c_y = np.fliplr(np.flipud(values))
    random_catalog_coords = SkyCoord.from_pixel(r_c_x, r_c_y, wcs=wcs, origin=1)
    #ax.scatter(r_c_x, r_c_y, c='r')
    #ax.imshow(mask)
    #x,y = real_catalog.to_pixel(wcs=wcs)
    #ax.scatter(x,y, c='b')

    ax.scatter(random_catalog_coords.ra.hour, random_catalog_coords.dec.degree, c='r')
    ax.scatter(real_catalog.ra.hour, real_catalog.dec.degree, c='b')
    fig.show()
    return random_catalog_coords

real_catalog = load_table("line_search_P3_wa_crop.out")
sn_cut = 9.0
real_catalog = real_catalog[real_catalog['rsnrrbin'] > sn_cut]
real_catalog = make_skycoords(real_catalog, ra='rra', dec='rdc')
random_catalog = generate_random_catalog(1000, "/home/jacob/Research/Wide_ASPECS/Data/gs_A1_2chn.fits")
random_catalog2 = generate_random_catalog(1000, "/home/jacob/Research/Wide_ASPECS/Data/gs_A1_2chn.fits")
print(real_catalog.shape)
#exit()

def angular_correlation_function(data_catalog, random_catalog):
    """
    Calculates the arrays for the data, random, and data_random for w(theta)
    :param data_catalog:
    :param random_catalog:
    :return:
    """
    distance_bins = np.logspace(0,np.log10(3000), 51) * u.arcsecond
    # First create the data data one
    seps = None


    # Get it for each one that is not the current ones
    for i, element in enumerate(data_catalog):
        #print(element)
        sep2d = element.separation(data_catalog[i+1:]).arcsecond
        #print(sep2d)
        if seps is None:
            seps = sep2d
        else:
            seps = np.concatenate((seps, sep2d))
    print(seps.shape)
    data_data, _ = histogram(seps, bins=distance_bins)
    print(data_data.shape)
    print(data_data)

    print("Done with Data Data")

    seps = None

    # Get it for each one that is not the current ones
    for i, element in enumerate(random_catalog):
        sep2d = element.separation(random_catalog[i+1:]).arcsecond
        if seps is None:
            seps = sep2d
        else:
            seps = np.concatenate((seps, sep2d))
    random_random, _ = np.histogram(seps, bins=distance_bins)

    print("Done with Random Random")

    seps = None

    # Get it for each one that is not the current ones
    for i, element in enumerate(data_catalog):
        sep2d = element.separation(random_catalog).arcsecond
        if seps is None:
            seps = sep2d
        else:
            seps = np.concatenate((seps, sep2d))
    data_random, _ = np.histogram(seps, bins=distance_bins)

    print("Done with Data Random")

    return data_data, data_random, random_random



def xi_r(data_array, data_random_array, random_array, real_catalog, random_catalog):
    """

    (DD/RR - 2 DR/RR +1)
    :param data_array:
    :param data_random_array:
    :param random_array:
    :return:
    """
    data_array_norm = (real_catalog.shape[0]*(real_catalog.shape[0]-1))/2.
    data_random_array_norm = (real_catalog.shape[0]*random_catalog.shape[0])
    random_array_norm = (random_catalog.shape[0]*(random_catalog.shape[0]-1))/2.

    data_array = data_array / data_array_norm
    data_random_array =  data_random_array / data_random_array_norm
    random_array = random_array / random_array_norm

    print(sum(random_array))
    print(sum(data_random_array))
    print(sum(data_array))

    # return 2 * (5000/real_catalog.shape[0])*(data_array/data_random_array) - 1

    # return 4 * (data_array*random_array)/(data_random_array)**2 - 1

    return data_array / random_array - (2 * (data_random_array / random_array)) + 1


def xi_r_error(omega_theta, data_array):
    """
    Data array is not normalized
    :param omega_theta:
    :param data_array:
    :return:
    """

    return (1 + omega_theta) / np.sqrt(data_array)


dd, dr, rr = angular_correlation_function(real_catalog, real_catalog)

omega_w = xi_r(dd, dr, rr, real_catalog, real_catalog)
e_omega_w = xi_r_error(omega_w, dd)
distance_bins = np.logspace(0,np.log10(3000), len(omega_w))

plt.cla()
plt.errorbar(x=distance_bins, y=omega_w, yerr=e_omega_w, fmt='o')
plt.xscale("log")
plt.title("Data vs Data")
plt.xlabel("Angular Distance (arcseconds)")
plt.ylabel("$\omega(\\theta)$")
#plt.yscale("log")
plt.tight_layout()
plt.savefig("Data_vs_Data_{}.png".format(sn_cut), dpi=300)
plt.show()

dd, dr, rr = angular_correlation_function(real_catalog, random_catalog)

omega_w = xi_r(dd, dr, rr, real_catalog, random_catalog)
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
plt.savefig("Data_vs_Random_{}.png".format(sn_cut), dpi=300)
plt.show()

dd, dr, rr = angular_correlation_function(random_catalog2, random_catalog)

omega_w = xi_r(dd, dr, rr,  random_catalog2, random_catalog)
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
plt.savefig("Random_vs_Random_{}.png".format(sn_cut), dpi=300)
plt.show()
