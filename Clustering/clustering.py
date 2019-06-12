import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from spectral_cube import SpectralCube
from astropy.coordinates import SkyCoord, Angle, SkyOffsetFrame, ICRS, Distance
from astropy.table import Table, hstack, join
from astropy.stats import histogram
from scipy.spatial.distance import cdist
from scipy.special import gamma
from astropy.cosmology import FlatLambdaCDM
from astropy import constants as const
from scipy.interpolate import interp2d, interp1d
from scipy.stats import norm
from scipy.optimize import leastsq, curve_fit

import matplotlib.mlab as mlab
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

def calc_gamma(beta):
    return 1 + beta


def H_gamma(gam):
    """
    Returns H_gamma for the r0 calculation
    :param gam: gamma
    :return:
    """
    return (gamma(0.5)*gamma(0.5*(gam-1)))/gamma(0.5*gam)


def H_z(z):
    """
    Calculates H(z)
    :param z:
    :return:
    """
    # return np.sqrt(cosmo.H0**2 * (0.3*(1+z)**3 + 0.7))
    return cosmo.H(z)


def E_z(z):
    """
    Calculates E_z for a given redshift
    :param z:
    :return:
    """
    return H_z(z) / const.c.to(u.km/u.s)

def Chi(z):
    """
    Returns the radial comoving distance, Dc, X in the equation
    :param z:
    :return:
    """

    return cosmo.comoving_distance(z)

def round_of_rating(number):
    """Round a number to the closest half integer.
    1.5
    2.5
    3.0
    4.0"""

    return round(number * 2) / 2

def redshift_distribution(table, use_matched=False):
    """
    Calculates redshift distribution in 0.1 redshift increments for all the sources

    If its noisy with empty bins, then need to make_curve and use that instead

    Get this from the matched ones


    :param table:
    :param bins:
    :return: A curve that captures the redshift distribution of the sample
    """
    # Get min and max Z

    min_z = np.min(table['Z (CO)'])
    max_z = np.max(table['Z (CO)'])

    bins = np.arange(1.5, 3.5, 0.5)
    if use_matched:
        only_matched = (table['Roberto ID'] > 0)
    else:
        only_matched = (table['Roberto ID'] > -1000000)
    values, bins = np.histogram(table[only_matched]['Z (CO)'], bins=bins)
    print(sum(values))
    #plt.hist(table[only_matched]['Z (CO)'], bins=bins)
    #plt.show()
    bin_centers = 0.5*(bins[1:]+bins[:-1]) # convert to centers
    mask = (values > 0)
    interp_values = values[mask]
    interp_bins = bin_centers[mask]
    interp_values = np.concatenate((np.asarray([0]), interp_values))
    interp_bins = np.concatenate((np.asarray([0]), interp_bins))

    def func(x, a, b, c):
        return a*x**2 + b*x + c

    popt, pcov = curve_fit(func, interp_bins, interp_values)

    xdata = np.linspace(0.1, np.max(interp_bins), 10000)
    f = interp1d(interp_bins, interp_values, kind='slinear')

    plt.hist(table[only_matched]['Z (CO)'], bins=bins)
    plt.title("SN > {}".format(np.round(np.min(table['S/N']), 2)))
    plt.ylabel("Count")
    plt.xlabel("Redshift (z)")
    plt.plot(xdata, func(xdata, *popt))
    plt.plot(xdata, f(xdata))
    plt.savefig("SN_{}_Redshift_Distribution.png".format(np.round(np.min(table['S/N']), 2)), dpi=300)
    plt.cla()
    return f, xdata


def calculate_r0(a, beta, table):
    """
    Calculates r0 given a beta and a, along with other things

    Calculates H_gamma*\frac{}{}

    :param a:
    :param beta:
    :param table: The Astropy Table from the matching that is used to get the
    Redshift Distribution for the sample
    :return:
    """

    # Convert A to radians to the 1, so need to raise to positive beta
    # A should be in units of arcseconds
    a_rad = a.radian ** (-1.*beta)

    print(a_rad)
    # Need to calc redshift distribution
    # Got that from the linear interpolation
    z_dist_func, zs = redshift_distribution(table)

    # Need to calc redshift vector E as a vector for every redshift Ez = Hz/c
    Ez = E_z(zs)
    print("Ez")
    print(Ez)
    # Need to calculate X -> Dc = DH as a vector for all redshifts too, can use Astropy distance
    # Radial comoving Distance, so I think all of the comoving distance
    # http://docs.astropy.org/en/stable/cosmology/

    X = Chi(zs)
    print("X")
    print(X)
    # So now have everything to get r0
    # Divide over the equation to have r0 on its own so
    # But write parts in normal Equation 16 order here
    # Need Integral here

    # Need to do elementwise multiplication of all of these together, then sum in integral
    top = z_dist_func(zs) * z_dist_func(zs) * Ez * X**(0.8)
    diff = zs[1] - zs[0]
    print("Top")
    print(top)
    print("Top Integrand")
    top_integrand = diff * sum(top)
    top = top_integrand
    # Need integral here
    bottom = z_dist_func(zs)**2
    print("Bottom")
    print(bottom)
    print("Bottom Integrated")
    bottom_integral = diff * sum(bottom)
    bottom = bottom_integral
    # Front
    front = H_gamma(calc_gamma(beta))
    print("Front")
    print(front)
    # Now put together with swap to other side
    print("Top, Bottom, Front")
    print(bottom)
    print(top)
    print(front)

    print("Top/Bottom*front")
    print(front*(top/bottom))

    r0_gamma = a_rad * (bottom/top) * (1/front)
    print(r0_gamma)
    r0 = r0_gamma**(-1.8)
    return r0


def load_table(ascii_table, header=0, start=1):
    ascii_table_data = Table.read(ascii_table, format="ascii", header_start=header, data_start=start)
    return ascii_table_data

sn8_table = Table.read("/home/jacob/Development/Wide_ASPECS/Final_Output/ASPECS_Line_Candidates_cleaned_all_closest_Sep_1.0_SN_5.5.ecsv")
#sn85_table = Table.read("/home/jacob/Development/Wide_ASPECS/Final_Output/ASPECS_Line_Candidates_cleaned_all_closest_Sep_1.0_SN_5.5.ecsv")
#sn9_table = Table.read("/home/jacob/Development/Wide_ASPECS/Final_Output/ASPECS_Line_Candidates_cleaned_all_closest_Sep_1.0_SN_6.15.ecsv")
#sn95_table = Table.read("/home/jacob/Development/Wide_ASPECS/Final_Output/ASPECS_Line_Candidates_cleaned_all_closest_Sep_1.0_SN_6.25.ecsv")
redshift_distribution(sn8_table)
print(calculate_r0(Angle(0.0005 * u.arcsecond), 0.8, sn8_table))
#redshift_distribution(sn85_table)
#redshift_distribution(sn9_table)
#redshift_distribution(sn95_table)
exit()
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
    #rand_x = np.random.randint(0, cube.unitless[0,:,:].shape[0]+1, 100*number_of_points)
    #rand_y = np.random.randint(0, cube.unitless[0,:,:].shape[1]+1, 100*number_of_points)
    random_catalog_coords = []
    r_c_x = []
    r_c_y = []
    r_c_y2 = []
    mask = cube.filled_data[0:30,:,:]
    mask = mask.mean(axis=0)

    xs = []
    ys = []
    xs2 = []
    ys2 = []

    while len(random_catalog_coords) < number_of_points:
        ra_random = np.random.uniform(low=53.037, high=53.213) * u.degree
        dec_random = (np.random.uniform(low=-27.8713, high=-27.737)) * u.degree
        c = SkyCoord(ra=ra_random, dec=dec_random)
        #print(c)

        x, y = c.to_pixel(wcs=wcs)
        if not np.isnan(mask[int(y),int(x)]):
            xs.append(x)
            ys.append(y)
            r_c_x.append(np.asarray([x,y]))
            random_catalog_coords.append(c)


    random_catalog_coords = SkyCoord(random_catalog_coords)
    plt.cla()
    plt.scatter(xs, ys, s=1)
    plt.title("X vs Y for RandInt")
    plt.savefig("X_V_Y_{}.png".format(len(random_catalog_coords)))
    plt.cla()

    #values = random_catalog_coords.to_pixel(wcs=wcs)
    #r_c_x, r_c_y = np.fliplr(np.flipud(values))
    #random_catalog_coords = SkyCoord.from_pixel(r_c_x, r_c_y, wcs=wcs, origin=1)
    #ax.scatter(r_c_x, r_c_y, c='r')
    #ax.imshow(mask)
    #x,y = real_catalog.to_pixel(wcs=wcs)
    #ax.scatter(x,y, c='b')

    ax.scatter(random_catalog_coords.ra.degree, random_catalog_coords.dec.degree, c='r', s=1,  transform=ax.get_transform('world'))
    ax.scatter(real_catalog.ra.degree, real_catalog.dec.degree, c='b', s=3, transform=ax.get_transform('world'))
    fig.show()
    #random_catalog_coords = np.asarray(r_c_x)
    return random_catalog_coords, np.asarray(r_c_x)

#exit()

def random_tester(random_cat, random_cat2):
    """
    Takes list of tuples of x,y coordinates and computes distances between them
    Should be a Numpy array of points
    :param random_cat:
    :param random_cat2:
    :return:
    """
    # First create the data data one
    data_data = None


    # Get it for each one that is not the current ones
    for i, element in enumerate(random_cat[:-1]):
        sep2d = cdist(element.reshape((1,2)), random_cat[i+1:], 'euclidean')
        #print(sep2d)
        if data_data is None:
            data_data = sep2d
        else:
            data_data = np.concatenate((data_data.reshape((-1,1)), sep2d.reshape((-1,1))))
    min_dist = np.min(data_data)
    print("Min Distance: {}".format(min_dist))
    min_dist = min_dist
    max_dist = np.max(data_data)

    print("Done with Data Data")

    random_random = None

    # Get it for each one that is not the current ones
    for i, element in enumerate(random_cat2[:-1]):
        sep2d = cdist(element.reshape((1,2)), random_cat[i+1:], 'euclidean')
        #print(sep2d)
        if random_random is None:
            random_random = sep2d
        else:
            random_random = np.concatenate((random_random.reshape((-1,1)), sep2d.reshape((-1,1))))

    m_dist = np.max(random_random)
    if m_dist > max_dist:
        max_dist = m_dist
    print("Done with Random Random")

    data_random = None

    # Get it for each one that is not the current ones
    for i, element in enumerate(random_cat):
        sep2d = cdist(element.reshape((1,2)), random_cat2.reshape((-1,2)), 'euclidean')
        if data_random is None:
            data_random = sep2d
        else:
            data_random = np.concatenate((data_random.reshape((-1,1)), sep2d.reshape((-1,1))))

    print("Done with Data Random")
    m_dist = np.max(data_random)
    if m_dist > max_dist:
        max_dist = m_dist

    return data_data, data_random, random_random, min_dist, max_dist



def angular_correlation_function(data_catalog, random_catalog):
    """
    Calculates the arrays for the data, random, and data_random for w(theta)
    :param data_catalog:
    :param random_catalog:
    :return:
    """
    # First create the data data one
    data_data = None


    # Get it for each one that is not the current ones
    for i, element in enumerate(data_catalog):
        #print(element)
        sep2d = element.separation(data_catalog[i+1:]).arcsecond
        #print(sep2d)
        if data_data is None:
            data_data = sep2d
        else:
            data_data = np.concatenate((data_data, sep2d))
    min_dist = np.min(data_data)
    print("Min Distance: {}".format(min_dist))
    min_dist = 13.5
    max_dist = np.max(data_data)

    print("Done with Data Data")

    random_random = None

    # Get it for each one that is not the current ones
    for i, element in enumerate(random_catalog):
        sep2d = element.separation(random_catalog[i+1:]).arcsecond
        if random_random is None:
            random_random = sep2d
        else:
            random_random = np.concatenate((random_random, sep2d))

    # Plot distribution of those to make sure random

    m_dist = np.max(random_random)
    if m_dist > max_dist:
        max_dist = m_dist
    print("Done with Random Random")

    data_random = None

    # Get it for each one that is not the current ones
    for i, element in enumerate(data_catalog):
        sep2d = element.separation(random_catalog).arcsecond
        if data_random is None:
            data_random = sep2d
        else:
            data_random = np.concatenate((data_random, sep2d))

    print("Done with Data Random")
    m_dist = np.max(data_random)
    if m_dist > max_dist:
        max_dist = m_dist

    return data_data, data_random, random_random, min_dist, max_dist


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

    # return 2 * (5000/real_catalog.shape[0])*(data_array/data_random_array) - 1

    # return 4 * (data_array*random_array)/(data_random_array)**2 - 1
    #print("Data-Data: {}".format(data_array))
    #print("RR: {}".format(random_array))
    #print("DR: {}".format(data_random_array))
    #print("DD/RR: {}".format(data_array/random_array))
    #print("DR/RR: {}".format(data_random_array/random_array))
    return (data_array / random_array) - 2 * (data_random_array / random_array) + 1


def xi_r_error(omega_theta, data_array):
    """
    Data array is not normalized

    w(theta) = A*theta^-beta with beta = 0.8

    cross correlation with the galaxy ones

    Plot for multiple SN and compare the slopes

    Same binning for all of them

    Small Poisson Error, works with small sample

    > 30 use sqrt number

    Gehrels 1986

    :param omega_theta:
    :param data_array:
    :return:
    """

    lower = [0.000, 0.173, 0.708, 1.367, 2.086, 2.840, 3.620, 4.419, 5.232, 6.057, 6.891, 7.734, 8.585, 9.441, 10.30, 11.17,
             12.04, 12.92, 13.80, 14.68, 15.57, 16.45, 17.35, 18.24, 19.14, 20.03, 20.93, 21.84, 22.74, 23.65, 24.55]

    upper = [1.841, 3.300, 4.638, 5.918, 7.163, 8.382, 9.584, 10.77, 11.95, 13.11, 14.27, 15.42, 16.56, 17.70, 18.83, 19.96,
             21.08, 22.20, 23.32, 24.44, 25.55, 26.66, 27.76, 28.87, 29.97, 31.07, 32.16, 33.26, 34.35, 35.45, 36.54]

    l_errors = []
    u_errors = []
    for index, element in enumerate(data_array):
        if element < 30:
            l_errors.append(lower[element])
            u_errors.append(upper[element])
        else:
            l_errors.append((1 + omega_theta[index]) / np.sqrt(element))
            u_errors.append((1 + omega_theta[index]) / np.sqrt(element))

    l_errors = np.asarray(l_errors)
    u_errors = np.asarray(u_errors)
    return l_errors, u_errors


def correlation_function(x, a):
    return a*(x**(-0.8))


from scipy.optimize import leastsq
fitfunc = lambda a, x: a*(x**(-0.8))
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

pinit = [1.0]

from scipy.optimize import curve_fit
negative = True
num_points = 7500
if negative:
    real_catalog = load_table("line_search_N3_wa_crop.out")
else:
    real_catalog = load_table("line_search_P3_wa_crop.out")
real_catalog = real_catalog[real_catalog['rsnrrbin'] > 8.5]
real_catalog = make_skycoords(real_catalog, ra='rra', dec='rdc')
np.random.seed(5227)
random_catalog, r_pixels = generate_random_catalog(num_points, "/media/jacob/A6548D38548D0BED/gs_A1_2chn.fits")
#np.random.seed(5227)
random_catalog2, r2_pixels = generate_random_catalog(num_points, "/media/jacob/A6548D38548D0BED/gs_A1_2chn.fits")

data_data, data_random, random_random, min_dist, max_dist = angular_correlation_function(random_catalog, random_catalog2)
from astroML.correlation import two_point
for bin_num in [5,6,7,8,9,10]:
    distance_bins = np.logspace(np.log10(min_dist-0.001),np.log10(max_dist+1), bin_num)
    distance_bins = np.concatenate((np.asarray([0.0]), distance_bins))
    dd, bins = histogram(data_data, bins=distance_bins)
    dr, _ = histogram(data_random, bins=distance_bins)
    rr, _ = histogram(random_random, bins=distance_bins)
    x_vals = []
    omega_w = xi_r(dd, dr, rr,  random_catalog, random_catalog2)
    le_omega_w, ue_omega_w = xi_r_error(omega_w, dd)
    distance_bins1 = 0.5*(distance_bins[1:]+distance_bins[:-1])
    #distance_bins1 = np.concatenate((np.asarray([1.5]), distance_bins1))
    # Best fit to the data
    pinit = [1.]
    out = leastsq(errfunc, pinit,
                  args=(distance_bins1, omega_w, (le_omega_w+ue_omega_w)/2.), full_output=True)
    a = out[0][0]
    s_sq = (errfunc(out[0][0], distance_bins1, omega_w, (le_omega_w+ue_omega_w)/2.)**2).sum()/(len(distance_bins1)-len(pinit))
    cov_matrix = out[1] * s_sq

    a_error = np.absolute(cov_matrix[0][0])**0.5
#    print("Value for A: {}".format(a))
    xmin=min_dist
    xmax=max_dist+0.1*max_dist
    x_fit = np.linspace(xmin, xmax, 10000)

    plt.cla()
    plt.errorbar(x=distance_bins1, y=omega_w, yerr=(le_omega_w, ue_omega_w), fmt='o')
    #plt.plot(x_fit, correlation_function(x_fit, params[0][0]), '-', c='r', label='Fit Normal: A = {}'.format(np.round( params[0][0],4)))
    plt.plot(x_fit, correlation_function(x_fit, a), '-', c='r', label='Fit: A = {}+-{}'.format(np.round(a,4), np.round(a_error,5)))
    #plt.scatter(x=distance_bins1, y=corr, c='r', label='Normal')
    #plt.scatter(x=distance_bins1, y=corr2, c='g', label='Other Random')
    plt.legend(loc='best')
    plt.xscale("log")
    plt.title("Random vs Random")
    plt.xlabel("Angular Distance (arcseconds)")
    plt.ylabel("$\omega(\\theta)$")
    #plt.yscale("log")
    #plt.tight_layout()
    plt.savefig("final/Random_vs_Random_{}NoParenFlip_bin{}.png".format(num_points, bin_num), dpi=300)

    plt.cla()
    plt.errorbar(x=distance_bins1, y=omega_w, yerr=(le_omega_w, ue_omega_w), fmt='o')
    #plt.plot(x_fit, correlation_function(x_fit, params[0][0]), '-', c='r', label='Fit Normal: A = {}'.format(np.round(a,4)))
    plt.plot(x_fit, correlation_function(x_fit, a), '-', c='r', label='Fit: A = {}+-{}'.format(np.round(a,4), np.round(a_error,5)))
    #plt.scatter(x=distance_bins1, y=corr, c='r', label='r')
    #plt.scatter(x=distance_bins1, y=corr2, c='g', label='Other Random')
    plt.legend(loc='best')
    plt.xscale("log")
    plt.title("Random vs Random")
    plt.xlabel("Angular Distance (arcseconds)")
    plt.ylabel("$\omega(\\theta)$")
    plt.yscale("log")
    #plt.tight_layout()
    plt.savefig("final/Log_Random_vs_Random_{}NoParenFlip_bin{}.png".format(num_points, bin_num), dpi=300)

#exit()

"""
4 panel plot

Show correlation between different SN cuts

Add taking into account the error on the A fitting

"""

dds = {}
drs = {}
rrs = {}
dist_bns = {}
dist_binners = {}
snners = [6.25, 6.1, 5.9, 5.5]

for sn_cut in snners:
    if negative:
        real_catalog = load_table("line_search_N3_wa_crop.out")
    else:
        real_catalog = load_table("line_search_P3_wa_crop.out")
    real_catalog = real_catalog[real_catalog['rsnrrbin'] > sn_cut]
    real_catalog = make_skycoords(real_catalog, ra='rra', dec='rdc')
    print(real_catalog.shape)
    dds[sn_cut] = []
    drs[sn_cut] = []
    rrs[sn_cut] = []
    data_data, data_random, random_random, min_dist, max_dist = angular_correlation_function(real_catalog, random_catalog)
    for bin_num in [5,6,7,8,9,10]:
        distance_bins = np.logspace(np.log10(min_dist-0.001),np.log10(max_dist+1), bin_num+1)
        dist_bns[bin_num] = distance_bins
        dd, _ = histogram(data_data, bins=distance_bins)
        dr, _ = histogram(data_random, bins=distance_bins)
        rr, _ = histogram(random_random, bins=distance_bins)
        dds[sn_cut].append(dd)
        drs[sn_cut].append(dr)
        rrs[sn_cut].append(rr)
        omega_w = xi_r(dd, dr, rr, real_catalog, random_catalog)
        le_omega_w, ue_omega_w = xi_r_error(omega_w, dd)
        distance_bins = 0.5*(distance_bins[1:]+distance_bins[:-1])
        #distance_bins = np.logspace(np.log10(min_dist),np.log10(max_dist), len(omega_w))
        dist_binners[bin_num] = distance_bins
        # Best fit to the data
        pinit = [1.]
        out = leastsq(errfunc, pinit,
                      args=(distance_bins, omega_w, (le_omega_w+ue_omega_w)/2.), full_output=True)
        a = out[0][0]
        s_sq = (errfunc(out[0][0], distance_bins, omega_w, (le_omega_w+ue_omega_w)/2.)**2).sum()/(len(distance_bins)-len(pinit))
        cov_matrix = out[1] * s_sq

        a_error = np.absolute(cov_matrix[0][0])**0.5
        print("Value for A: {}".format(a))

        # Now get one for only the positive points
        pos_mask = (omega_w > 0.)
        pos_ue = ue_omega_w[pos_mask]
        pos_le = le_omega_w[pos_mask]
        pos_bins = distance_bins[pos_mask]
        pos_omega = omega_w[pos_mask]

        pinit = [1.]
        out = leastsq(errfunc, pinit,
                      args=(pos_bins, pos_omega, (pos_le+pos_ue)/2.), full_output=True)
        pos_a = out[0][0]
        s_sq = (errfunc(out[0][0], pos_bins, pos_omega, (pos_le+pos_ue)/2.)**2).sum()/(len(pos_bins)-len(pinit))
        cov_matrix = out[1] * s_sq

        pos_a_error = np.absolute(cov_matrix[0][0])**0.5
        print("Value for A: {}".format(a))

        xmin=distance_bins[0]-0.001
        xmax=distance_bins[-1]+0.1*distance_bins[-1]
        x_fit = np.linspace(xmin, xmax, 10000)
        plt.cla()
        plt.errorbar(x=distance_bins, y=omega_w, yerr=(le_omega_w,ue_omega_w), fmt='o')
        plt.plot(x_fit, correlation_function(x_fit, a), '-', c='r', label='Fit: A = {}+-{}'.format(np.round(a,4), np.round(a_error,5)))
        plt.plot(x_fit, correlation_function(x_fit, pos_a), '--', c='g', label='Pos Fit: A = {}+-{}'.format(np.round(pos_a,4), np.round(pos_a_error,5)))
        plt.xscale("log")
        plt.title("Data vs Random")
        plt.legend(loc='best')
        plt.xlabel("Angular Distance (arcseconds)")
        plt.ylabel("$\omega(\\theta)$")
        plt.yscale("log")
        #plt.tight_layout()
        plt.savefig("final/Log_Data_vs_Random_{}_bin{}_sn{}.png".format(num_points, bin_num, sn_cut), dpi=300)
        #plt.show()

        plt.cla()
        plt.errorbar(x=distance_bins, y=omega_w, yerr=(le_omega_w,ue_omega_w), fmt='o')
        plt.plot(x_fit, correlation_function(x_fit, a), '-', c='r', label='Fit: A = {}+-{}'.format(np.round(a,4), np.round(a_error,5)))
        plt.plot(x_fit, correlation_function(x_fit, pos_a), '--', c='g', label='Pos Fit: A = {}+-{}'.format(np.round(pos_a,4), np.round(pos_a_error,5)))
        plt.xscale("log")
        plt.title("Data vs Random")
        plt.legend(loc='best', fontsize='5')
        plt.xlabel("Angular Distance (arcseconds)")
        plt.ylabel("$\omega(\\theta)$")
        #plt.yscale("log")
        #plt.tight_layout()
        plt.savefig("final/Data_vs_Random_{}_bin{}_sn{}.png".format(num_points, bin_num, sn_cut), dpi=300)
        #plt.show()

# Now have dictionary of lists of datas


def plot_four(dd, dr, rr, distance_bins, distance_bins1, use_log=True):
    """
    Take set of 4 S/N cut lists for DD, DR, and RR and plot a 4 panel plot, has to be same binning for it
    Assumes dd, dr, and rr are all in same order, 9.5, 9.0, 8.5, 8.0 SN
    :param dd:
    :param dr:
    :param rr:
    :return:
    """
    plt.cla()
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='all', sharey='all', figsize=(10,10))
    sn_cut = snners
    for index, data_data in enumerate(dd):
        if negative:
            real_catalog = load_table("line_search_N3_wa_crop.out")
        else:
            real_catalog = load_table("line_search_P3_wa_crop.out")
        real_catalog = real_catalog[real_catalog['rsnrrbin'] > sn_cut[index]]
        real_catalog = make_skycoords(real_catalog, ra='rra', dec='rdc')
        omega_w = xi_r(dd[index], dr[index], rr[index], real_catalog, random_catalog)
        le_omega_w, ue_omega_w = xi_r_error(omega_w, dd[index])
        pinit = [1.]
        out = leastsq(errfunc, pinit,
                      args=(distance_bins1, omega_w, (le_omega_w+ue_omega_w)/2.), full_output=True)
        a = out[0][0]
        s_sq = (errfunc(out[0][0], distance_bins1, omega_w, (le_omega_w+ue_omega_w)/2.)**2).sum()/(len(distance_bins1)-len(pinit))
        cov_matrix = out[1] * s_sq

        a_error = np.absolute(cov_matrix[0][0])**0.5
        print("Value for A: {}".format(a))
        # Now get one for only the positive points
        pos_mask = (omega_w > 0.)
        pos_ue = ue_omega_w[pos_mask]
        pos_le = le_omega_w[pos_mask]
        pos_bins = distance_bins1[pos_mask]
        pos_omega = omega_w[pos_mask]

        pinit = [1.]
        out = leastsq(errfunc, pinit,
                      args=(pos_bins, pos_omega, (pos_le+pos_ue)/2.), full_output=True)
        pos_a = out[0][0]
        s_sq = (errfunc(out[0][0], pos_bins, pos_omega, (pos_le+pos_ue)/2.)**2).sum()/(len(pos_bins)-len(pinit))
        cov_matrix = out[1] * s_sq

        pos_a_error = np.absolute(cov_matrix[0][0])**0.5
        print("Value for A: {}".format(a))
        xmin=min_dist-0.001
        xmax=max_dist+0.1*max_dist
        x_fit = np.linspace(xmin, xmax, 10000)
        if index == 0:
            ax1.errorbar(x=distance_bins1, y=omega_w, yerr=(le_omega_w,ue_omega_w), fmt='o')
            ax1.plot(x_fit, correlation_function(x_fit, a), '-', c='r', label='Fit: A = {}+-{}'.format(np.round(a,4), np.round(a_error, 5)))
            ax1.plot(x_fit, correlation_function(x_fit, pos_a), '--', c='g', label='Pos Fit: A = {}+-{}'.format(np.round(pos_a,4), np.round(pos_a_error,5)))
            ax1.set_xscale("log")
            ax1.set_title("S/N > {}".format(snners[index]))
            #plt.title("Data vs Random")
            ax1.legend(loc='best', fontsize='8')
            ax1.set_ylabel("$\omega(\\theta)$")
            #ax1.set_ylabel("$\omega(\\theta)$")
            if use_log:
                ax1.set_yscale("log")
            #plt.savefig("final/Log_Data_vs_Random_{}_bin{}_sn{}.png".format(num_points, bin_num, sn_cut), dpi=300)
        if index == 1:
            ax2.errorbar(x=distance_bins1, y=omega_w, yerr=(le_omega_w,ue_omega_w), fmt='o')
            ax2.plot(x_fit, correlation_function(x_fit, a), '-', c='r', label='Fit: A = {}+-{}'.format(np.round(a,4), np.round(a_error, 5)))
            ax2.plot(x_fit, correlation_function(x_fit, pos_a), '--', c='g', label='Pos Fit: A = {}+-{}'.format(np.round(pos_a,4), np.round(pos_a_error,5)))
            ax2.set_xscale("log")
            ax2.set_title("S/N > {}".format(snners[index]))
            #plt.title("Data vs Random")
            ax2.legend(loc='best', fontsize='8')
            #plt.xlabel("Angular Distance (arcseconds)")
            #plt.ylabel("$\omega(\\theta)$")
            if use_log:
                ax2.set_yscale("log")
            #plt.savefig("final/Log_Data_vs_Random_{}_bin{}_sn{}.png".format(num_points, bin_num, sn_cut), dpi=300)
        if index == 2:
            ax3.errorbar(x=distance_bins1, y=omega_w, yerr=(le_omega_w,ue_omega_w), fmt='o')
            ax3.plot(x_fit, correlation_function(x_fit, a), '-', c='r', label='Fit: A = {}+-{}'.format(np.round(a,4), np.round(a_error, 5)))
            ax3.plot(x_fit, correlation_function(x_fit, pos_a), '--', c='g', label='Pos Fit: A = {}+-{}'.format(np.round(pos_a,4), np.round(pos_a_error,5)))
            ax3.set_xscale("log")
            ax3.set_title("S/N > {}".format(snners[index]))
            #plt.title("Data vs Random")
            ax3.legend(loc='best', fontsize='8')
            ax3.set_xlabel("Angular Distance (arcseconds)")
            ax3.set_ylabel("$\omega(\\theta)$")
            if use_log:
                ax3.set_yscale("log")
            #plt.savefig("final/Log_Data_vs_Random_{}_bin{}_sn{}.png".format(num_points, bin_num, sn_cut), dpi=300)
        if index == 3:
            ax4.errorbar(x=distance_bins1, y=omega_w, yerr=(le_omega_w,ue_omega_w), fmt='o')
            ax4.plot(x_fit, correlation_function(x_fit, a), '-', c='r', label='Fit: A = {}+-{}'.format(np.round(a,4), np.round(a_error, 5)))
            ax4.plot(x_fit, correlation_function(x_fit, pos_a), '--', c='g', label='Pos Fit: A = {}+-{}'.format(np.round(pos_a,4), np.round(pos_a_error,5)))
            ax4.set_xscale("log")
            ax4.set_title("S/N > {}".format(snners[index]))
            #plt.title("Data vs Random")
            ax4.legend(loc='best', fontsize='8')
            #plt.xlabel("Angular Distance (arcseconds)")
            if use_log:
                ax4.set_yscale("log")
            ax4.set_xlabel("Angular Distance (arcseconds)")
            #plt.savefig("final/Log_Data_vs_Random_{}_bin{}_sn{}.png".format(num_points, bin_num, sn_cut), dpi=300)
            # Now add the plot to the others
            ax1.plot(x_fit, correlation_function(x_fit, a), '--', c='orange', label='S/N>8. Fit')
            ax2.plot(x_fit, correlation_function(x_fit, a), '--', c='orange', label='S/N>8. Fit')
            ax3.plot(x_fit, correlation_function(x_fit, a), '--', c='orange', label='S/N>8. Fit')

    # Actually plot it all now
    f.align_xlabels()
    f.align_ylabels()
    f.suptitle("Data vs Random")
    if use_log:
        plt.savefig("final/Log_4Panel_Data_Vs_Random_bin{}_N{}.png".format(len(distance_bins1), negative), dpi=300)
    else:
        plt.savefig("final/4Panel_Data_Vs_Random_bin{}_N{}.png".format(len(distance_bins1), negative), dpi=300)

num_ones = len(dds[snners[2]])

for i in range(num_ones):
    dd_plot = []
    dr_plot = []
    rr_plot = []

    for key in dds.keys():
        dd_plot.append(dds[key][i])
        dr_plot.append(drs[key][i])
        rr_plot.append(rrs[key][i])

    # Number of bins
    print(dd_plot)
    print(np.asarray(dd_plot).shape)
    bin_num = len(dd_plot[0])
    distance_bins = dist_bns[bin_num]
    distance_bins1 = dist_binners[bin_num]
    plot_four(dd_plot, dr_plot, rr_plot, distance_bins, distance_bins1)
    plot_four(dd_plot, dr_plot, rr_plot, distance_bins, distance_bins1, use_log=False)
