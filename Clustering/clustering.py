import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from spectral_cube import SpectralCube
from astropy.coordinates import SkyCoord, Angle, SkyOffsetFrame, ICRS, Distance
from astropy.table import Table, hstack, join, Column, vstack
from astropy.stats import histogram
from scipy.spatial.distance import cdist
from scipy.special import gamma
from astropy.io import ascii
from astropy.cosmology import FlatLambdaCDM
from astropy import constants as const
from scipy.interpolate import interp2d, interp1d
from scipy.stats import norm
from astropy import modeling
from scipy.integrate import newton_cotes
from scipy.optimize import leastsq, curve_fit

import matplotlib.mlab as mlab


def idl_calc_r0(a, beta, z_dist):
    """
    a is alpha
    beta is beta
    z_dist is 2D array with z in the first vector and N in the second one.
    :param a:
    :param beta:
    :param z_dist:
    :return:
    """
    clight = 2.9979246e5 * (u.km/u.s)
    HORIZON_MPC = 2.9979246e3 # * (u.Mpc / u.littleh)

    # Compute Hg
    gamma_par = beta + 1.
    Hg = (gamma(1.0/2.0)*gamma((gamma_par-1.0)/2.0))/gamma(gamma_par/2.0)

    # Computing Integral of redshift distribution
    int_zdist1 = idl_tabulate(z_dist[:,0], z_dist[:,1])

    # Compute Ez
    Hz = (100.0 * (u.littleh * u.km / u.s / u.Mpc))
    Ez = E_z(z_dist[:,0])

    # Compute comoving distance, should be cMpc/h
    zcom = HORIZON_MPC * cosmo.comoving_distance(z_dist[:,0])

    # Compute Big Integral
    expr = z_dist[:,1]*z_dist[:,1]*Ez*zcom**(1.-gamma_par)
    int_big = idl_tabulate(z_dist[:,0], expr)

    # Comput r0 using inverse Limber equation
    A_param_rad = a * ((np.pi/180.0)/3600.) ** beta
    r0_param = ((A_param_rad * int_zdist1 * int_zdist1) / (Hg * int_big)) ** (1.0/gamma_par) # in cMpc/h units

    return r0_param




def combine_catalogs(catalog_one, catalog_two):
    combined = join(catalog_one, catalog_two, keys='id')
    return combined


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

initial_catalog = Table.read(
    "/home/jacob/Development/Wide_ASPECS/independent/jacob_mapghys_in_nov2018_all_jcb4_magphys_jcb4.fits",
    format='fits')  # hdu_list[1].data
roberto_catalog = Table.read("roberto_catalog_muse_skelton_matched_manFix.fits", format='fits')

combined_catalog = combine_catalogs(initial_catalog, roberto_catalog)

gal_catalog = make_skycoords(combined_catalog, ra='ra', dec='dc')


cosmo = FlatLambdaCDM(H0=100, Om0=0.3, Tcmb0=2.725)


def load_table(ascii_table, header=0, start=1):
    ascii_table_data = Table.read(ascii_table, format="ascii", header_start=header, data_start=start)
    return ascii_table_data


def fid(neg, pos):
    return 1. - (len(neg) / len(pos))


neg_catalog = load_table("line_search_N3_wa_crop.out")
pos_catalog = load_table("line_search_P3_wa_crop.out")

line_widths = [i for i in range(3, 21, 2)]

sn_cuts = np.arange(5., 8.1, 0.05)
six_fids = []
for width in line_widths:
    neg_widths = neg_catalog[neg_catalog['width'] == width]
    pos_widths = pos_catalog[pos_catalog['width'] == width]
    print("Neg Width Lines: {}".format(len(neg_widths)))
    print("Pos Width Lines: {}".format(len(pos_widths)))
    print("Width {} MaxSN: {}".format(width, np.max(neg_widths['rsnrrbin'])))
    fid_width = []
    sn_vals = []
    six_fid = -1
    for sn in sn_cuts:
        neg_sn = neg_widths[neg_widths['rsnrrbin'] >= sn]
        pos_sn = pos_widths[pos_widths['rsnrrbin'] >= sn]
        # print("SN: {} Len: {}".format(sn, len(pos_sn)))
        if len(pos_sn) > 0:
            fid_width.append(fid(neg_sn, pos_sn))
            sn_vals.append(sn)
            if six_fid < 0 and fid_width[-1] >= 0.6:
                six_fid = sn
        elif len(neg_sn) == 0:
            fid_width.append(1)
            sn_vals.append(sn)
            if six_fid < 0 and fid_width[-1] >= 0.6:
                six_fid = sn
        six_fids.append(six_fid)


def construct_fid_mask(catalog):
    """
    Constructs the fidelity mask based off my results, not Robertos
    :param catalog:
    :return:
    """
    masks = []
    six_fids = [6.25, 6.2, 6.1, 6.1, 6.1, 6.15, 6.1, 6.15, 6.05]
    #six_fids = [6.35, 6.25, 6.15, 6.15, 6.15, 6.25, 6.15, 6.25, 6.05]
    #six_fids = [6.3, 6.2, 6.1, 6.15, 6.1, 6.20, 6.1, 6.20, 6.05]
    for index, width in enumerate(line_widths):
        print(six_fids[index])
        masks.append(catalog[((catalog['width'] == width) & (catalog['rsnrrbin'] >= six_fids[index]))])

    total = masks[0]
    t_sum = 0
    for mask in masks[1:]:
        t_sum += len(mask)
        total = vstack((total, mask))

    print("Total One: {}".format(len(total)))
    return total

print(" Num above Fid: 0.6: " + str(len(construct_fid_mask(pos_catalog))))
print( " Num above 6.25: " + str(len(pos_catalog[pos_catalog['rsnrrbin'] > 6.25])))
print( " Num above 6.1: " + str(len(pos_catalog[pos_catalog['rsnrrbin'] > 6.1])))


def idl_tabulate(x, f, p=5) :
    def idl_newton_cotes(x, f) :
        if x.shape[0] < 2 :
            return 0
        rn = (x.shape[0] - 1) * (x - x[0]) / (x[-1] - x[0])
        weights = newton_cotes(rn)[0]
        return (x[-1] - x[0]) / (x.shape[0] - 1) * np.dot(weights, f)
    ret = 0
    for idx in range(0, x.shape[0], p - 1) :
        ret += idl_newton_cotes(x[idx:idx + p], f[idx:idx + p])
    return ret


def calc_gamma(beta):
    return beta + 1


def H_gamma(gam):
    """
    Returns H_gamma for the r0 calculation
    :param gam: gamma
    :return:

    Hg = (gamma(0.5)*gamma((gamma_par-1.0)/2.0))/gamma(gamma_par/2.0)
    """
    return (gamma(0.5) * gamma(0.5 * (gam - 1))) / gamma(0.5 * gam)


def H_z(z):
    """
    Calculates H(z)
    :param z:
    :return:
    """
    # return np.sqrt(cosmo.H0**2 * (0.3*(1+z)**3 + 0.7))
    # Returns the same as the equation above from the paper
    print(" H(z): {}".format(cosmo.H(z)))
    return cosmo.H(z)


def E_z(z):
    """
    Calculates E_z for a given redshift

    clight = double(2.9979246e5) ;km/s

    HORIZON_MPC = double(2.9979246e3) ; HORIZON in Mpc/h

    int_zdist1 = int_tabulated(z_dist1[*,0], z_dist1[*,1], /DOUBLE)
    Integrates from the beginning to the end of it all


    ;Computing Ez
    Hz = 100.0*bigh(z_dist1[*,0], omega_m, omega_v, w)
    ; Hz=H0*E(z) for the array of z, with H0=100h km/s/Mpc.
    bigh is a function that compute the value of E(z) for the provided z and cosmology
    Ez = (Hz/clight) ; in units of cMpc/h

    :param z:
    :return:
    """
    print(" Ez(z): {}".format(H_z(z) / const.c.to(u.km / u.s)))
    return H_z(z) / const.c.to(u.km / u.s)


def Chi(z):
    """
    Returns the radial comoving distance, Dc, X in the equation

    ;Computing comoving distance
    a = 1.0d/(1.0d + z_dist1[*,0])
    zcom = HORIZON_MPC*dofa(a, omega_m, omega_v, w)
    ; in units of cMpc/h. dofa is a function that compute
    the comoving distance for the provided z (or alternatively a)

    :param z:
    :return:
    """
    HORIZON_MPC = 2.9979246e3 * (u.Mpc / u.littleh)
    return 1 * cosmo.comoving_distance(z)


def round_of_rating(number):
    """Round a number to the closest half integer.
    1.5
    2.5
    3.0
    4.0"""

    return round(number * 2) / 2

from scipy import arange, array, exp

def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return array(list(map(pointwise, array(xs))))

    return ufunclike

def redshift_distribution(table, use_matched=False, plot=False):
    """
    Calculates redshift distribution in 0.1 redshift increments for all the sources

    If its noisy with empty bins, then need to make_curve and use that instead

    Get this from the matched ones

    ;Computing big integral
    expr = z_dist1[*,1]*z_dist2[*,1]*Ez*zcom^(1.0-gamma_par)
    int_big = int_tabulated(z_dist1[*,0], expr, /DOUBLE)

    :param table:
    :param bins:
    :return: A curve that captures the redshift distribution of the sample
    """
    # Get min and max Z

    min_z = np.min(table['Z (CO)'])
    max_z = np.max(table['Z (CO)'])
    width_bin = 0.5
    #print("Max Z: {}".format(max_z))
    bins = np.arange(0, max_z+0.6, width_bin)
    if use_matched:
        only_matched = (table['Roberto ID'] > 0)
    else:
        only_matched = (table['Roberto ID'] > -1000000)
    values, bins = np.histogram(table[only_matched]['Z (CO)'], bins=bins)#, density=True)
    # Multiply by the bin width to get the values to sum to 1
    #values = values * width_bin
    # plt.hist(table[only_matched]['Z (CO)'], bins=bins)
    # plt.show()
    bin_centers = 0.5 * (bins[1:] + bins[:-1])  # convert to centers
    mask = (values > 0)
    interp_values = values[mask]
    #print(sum(interp_values))
    interp_bins = bin_centers[mask]

    #interp_values = np.concatenate((np.asarray([1.5]), interp_values))
    #interp_bins = np.concatenate((np.asarray([values[0]]), interp_bins))
    #interp_values = np.concatenate((interp_values, np.asarray([3.5])))
    #interp_bins = np.concatenate((interp_bins, np.asarray([values[-1]])))
    #print("Min Intern bins: {}".format(min(interp_bins)))
    #print("Max Intern bins: {}".format(max(interp_bins)))

    def func(x, a, b, c):
        return a * x ** 2 + b * x + c

    popt, pcov = curve_fit(func, interp_bins, interp_values)

    xdata = np.linspace(1.5, 3.5, 10000)
    f = interp1d(interp_bins, interp_values, kind='slinear')
    f = extrap1d(f)

    # Fit Gaussian
    fitter = modeling.fitting.LevMarLSQFitter()
    model = modeling.models.Gaussian1D()  # depending on the data you need to give some initial values
    fitted_model = fitter(model, interp_bins, interp_values)
    if plot:
        plt.plot(xdata, fitted_model(xdata))
        _, _, _ = plt.hist(bin_centers, weights=values, bins=bins)
        # values = values * ((bins[0] + bins[1]) * 0.5)
        plt.title("SN > {}".format(np.round(np.min(table['S/N']), 2)))
        plt.ylabel("Count")
        plt.xlabel("Redshift (z)")
        plt.plot(xdata, func(xdata, *popt))
        plt.plot(xdata, f(xdata))
        plt.savefig("SN_{}_Redshift_Distribution_60.png".format(np.round(np.min(table['S/N']), 2)), dpi=300)
        plt.cla()
        t = Table()
        t['Z'] = Column(xdata, description='Redshift')
        t['Num_Gal'] = Column(f(xdata), description='Num Gal at Redshift')
        ascii.write(t, "redshift_distribution_interpolated_60.txt")

        t = Table()
        t['Z'] = Column(xdata, description='Redshift')
        t['Num_Gal'] = Column(fitted_model(xdata), description='Num Gal at Redshift')
        ascii.write(t, "redshift_distribution_interpolated_gauss_60.txt")

        t = Table()
        t['Z'] = Column(bin_centers, description='Redshift')
        t['Num_Gal'] = Column(values, description='Num Gal at Redshift')
        ascii.write(t, "redshift_distribution_points_60.txt")

    return f, xdata, fitted_model


def calculate_r0(a, beta, table, use_gauss=False, plot=False):
    """
    Calculates r0 given a beta and a, along with other things

    Calculates H_gamma*\frac{}{}

    ;Computing r0 using the inverse of Limber equation
    A_param_rad = A_param*((!dpi/180.0d)/3600.0d)^(beta_param) ; A in units of radians^beta
    r0_param = ((A_param_rad * int_zdist1 * int_zdist2) / (Hg * int_big))^(1.0/gamma_par) ;in cMpc/h units

    :param a:
    :param beta:
    :param table: The Astropy Table from the matching that is used to get the
    Redshift Distribution for the sample
    :return:
    """

    # Convert A to radians to the 1, so need to raise to positive beta
    # A should be in units of arcseconds
    a_rad = a.arcsec * ((np.pi/180.0)/3600.) ** beta
    print("A Rad: {}".format(a_rad))
    #a_rad = a.radian ** (beta)
    #print("A Rad 2: {}".format(a_rad))
    # Need to calc redshift distribution
    # Got that from the linear interpolation
    z_dist_func, zs, gaussian = redshift_distribution(table, plot=plot)

    # Need to calc redshift vector E as a vector for every redshift Ez = Hz/c
    Ez = E_z(zs)
    # Need to calculate X -> Dc = DH as a vector for all redshifts too, can use Astropy distance
    # Radial comoving Distance, so I think all of the comoving distance
    # http://docs.astropy.org/en/stable/cosmology/
    #print("E(z): {}".format(Ez))

    X = Chi(zs)
    #print("X: {}".format(X))
    # So now have everything to get r0
    # Divide over the equation to have r0 on its own so
    # But write parts in normal Equation 16 order here
    # Need Integral here

    # Need to do elementwise multiplication of all of these together, then sum in integral
    """
    
        ;Computing big integral
    expr = z_dist1[*,1]*z_dist2[*,1]*Ez*zcom^(1.0-gamma_par)
    int_big = int_tabulated(z_dist1[*,0], expr, /DOUBLE)
    
    r0_param = ((A_param_rad * int_zdist1 * int_zdist2) / (Hg * int_big))^(1.0/gamma_par) ;in cMpc/h units
    """
    chi_pow = 1 - calc_gamma(beta)
    if use_gauss:
        top = gaussian(zs) * gaussian(zs) * Ez * np.power(X, chi_pow)
    else:
        top = z_dist_func(zs) * z_dist_func(zs) * Ez * np.power(X, chi_pow)
    # Now integrate over the Z_distribution with the expression
    print("Top: {}".format(top))
    top_unit = top.unit
    top_integrand = idl_tabulate(zs, top) * top.unit
    print("Top: {}".format(top_integrand))
    #diff = zs[1] - zs[0]
    #top_integrand = diff * sum(top)
    top = top_integrand
    # Need integral here
    if use_gauss:
        bottom_integral = 0.5 * sum(gaussian(zs))
    else:
        bottom_integral = 0.5 * sum(z_dist_func(zs))

    bottom_integral = idl_tabulate(zs, gaussian(zs))
    print("Bottom: {}".format(bottom_integral))
    bottom = bottom_integral ** 2
    # Front
    front = H_gamma(calc_gamma(beta))
    # Now put together with swap to other side
    r0_gamma = a_rad * (bottom / (front * top))
    #print("R0Gamm: {}".format(r0_gamma))
    r0 = r0_gamma ** (1 / calc_gamma(beta))
    # Convert to cMpc/h
    print("Before Convert: {}".format(r0))
    r0 = r0.to((u.Mpc / u.littleh), u.with_H0(cosmo.H0))
    return r0

sn8_table = Table.read(
    "/home/jacob/Development/Wide_ASPECS/Final_Output/ASPECS_Line_Candidates_cleaned_all_closest_Sep_1.0_SN_fid_0.6.ecsv")
# sn85_table = Table.read("/home/jacob/Development/Wide_ASPECS/Final_Output/ASPECS_Line_Candidates_cleaned_all_closest_Sep_1.0_SN_5.5.ecsv")
# sn9_table = Table.read("/home/jacob/Development/Wide_ASPECS/Final_Output/ASPECS_Line_Candidates_cleaned_all_closest_Sep_1.0_SN_6.15.ecsv")
sn95_table = Table.read(
    "/home/jacob/Development/Wide_ASPECS/Final_Output/ASPECS_Line_Candidates_cleaned_all_closest_Sep_1.0_SN_fid_0.6.ecsv")
print(calculate_r0(Angle(1.1703 * u.arcsecond), 0.8, sn95_table))
print(calculate_r0(Angle((1.1703+ 0.16) * u.arcsecond), 0.8, sn95_table, use_gauss=False, plot=True))
print(calculate_r0(Angle((1.1703 - 0.16) * u.arcsecond), 0.8, sn95_table, use_gauss=False))
print(calculate_r0(Angle(1.1703 * u.arcsecond), 0.8, sn95_table, use_gauss=True, plot=True))
print(calculate_r0(Angle((1.1703 + 0.16) * u.arcsecond), 0.8, sn95_table, use_gauss=True))
print(calculate_r0(Angle((1.1703 - 0.16) * u.arcsecond), 0.8, sn95_table, use_gauss=True))

# redshift_distribution(sn85_table)
# redshift_distribution(sn9_table)
# redshift_distribution(sn95_table)
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


def trim_galaxy(points, filename):
    """
    Removes galaxies that do not fall within the Wide ASPECS footprint in position of Wide ASPECS  and within the mask

    :param number_of_points:
    :param mask:
    :return:
    """
    cube = SpectralCube.read(filename)
    cube[0, :, :].quicklook()  # uses aplpy
    # using wcsaxes
    wcs = cube[0, :, :].wcs
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=wcs)
    # ax.imshow(cube.unitless[0,:,:]) # you may need cube[5,:,:].value depending on mpl version
    # fig.show()

    # Have the mask now, so can create from random pixel coordinates, generate SkyCoord

    # Just do it 10 times the number of points, hoping enough are in the valid area
    # rand_x = np.random.randint(0, cube.unitless[0,:,:].shape[0]+1, 100*number_of_points)
    # rand_y = np.random.randint(0, cube.unitless[0,:,:].shape[1]+1, 100*number_of_points)
    random_catalog_coords = []
    r_c_x = []
    r_c_y = []
    r_c_y2 = []
    mask = cube.filled_data[0:30, :, :]
    mask = mask.mean(axis=0)

    xs = []
    ys = []
    xs2 = []
    ys2 = []

    for point in points:
        if 53.037 <= point.ra.degree <= 53.213:
            if -27.8713 <= point.dec.degree <= -27.737:
                c = point
                # print(c)

                x, y = c.to_pixel(wcs=wcs)
                if x < mask.shape[0] and y < mask.shape[1]:
                    if not np.isnan(mask[int(y), int(x)]):
                        xs.append(x)
                        ys.append(y)
                        r_c_x.append(np.asarray([x, y]))
                        random_catalog_coords.append(c)

    random_catalog_coords = SkyCoord(random_catalog_coords)

    plt.cla()
    plt.scatter(xs, ys, s=1)
    plt.title("X vs Y for RandInt")
    plt.savefig("X_V_Y_{}.png".format(len(random_catalog_coords)))
    plt.cla()

    # values = random_catalog_coords.to_pixel(wcs=wcs)
    # r_c_x, r_c_y = np.fliplr(np.flipud(values))
    # random_catalog_coords = SkyCoord.from_pixel(r_c_x, r_c_y, wcs=wcs, origin=1)
    # ax.scatter(r_c_x, r_c_y, c='r')
    # ax.imshow(mask)
    # x,y = real_catalog.to_pixel(wcs=wcs)
    # ax.scatter(x,y, c='b')

    ax.scatter(random_catalog_coords.ra.degree, random_catalog_coords.dec.degree, c='r', s=1,
               transform=ax.get_transform('world'))
    ax.scatter(real_catalog.ra.degree, real_catalog.dec.degree, c='b', s=3, transform=ax.get_transform('world'))
    fig.show()
    # random_catalog_coords = np.asarray(r_c_x)
    return random_catalog_coords, np.asarray(r_c_x)


def generate_random_catalog(number_of_points, filename):
    """
    Generates random number of points in position of Wide ASPECS  and within the mask

    :param number_of_points:
    :param mask:
    :return:
    """
    cube = SpectralCube.read(filename)
    cube[0, :, :].quicklook()  # uses aplpy
    # using wcsaxes
    wcs = cube[0, :, :].wcs
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=wcs)
    # ax.imshow(cube.unitless[0,:,:]) # you may need cube[5,:,:].value depending on mpl version
    # fig.show()

    # Have the mask now, so can create from random pixel coordinates, generate SkyCoord

    # Just do it 10 times the number of points, hoping enough are in the valid area
    # rand_x = np.random.randint(0, cube.unitless[0,:,:].shape[0]+1, 100*number_of_points)
    # rand_y = np.random.randint(0, cube.unitless[0,:,:].shape[1]+1, 100*number_of_points)
    random_catalog_coords = []
    r_c_x = []
    r_c_y = []
    r_c_y2 = []
    mask = cube.filled_data[0:30, :, :]
    mask = mask.mean(axis=0)

    xs = []
    ys = []
    xs2 = []
    ys2 = []

    while len(random_catalog_coords) < number_of_points:
        ra_random = np.random.uniform(low=53.037, high=53.213) * u.degree
        dec_random = (np.random.uniform(low=-27.8713, high=-27.737)) * u.degree
        c = SkyCoord(ra=ra_random, dec=dec_random)
        # print(c)

        x, y = c.to_pixel(wcs=wcs)
        if not np.isnan(mask[int(y), int(x)]):
            xs.append(x)
            ys.append(y)
            r_c_x.append(np.asarray([x, y]))
            random_catalog_coords.append(c)

    random_catalog_coords = SkyCoord(random_catalog_coords)

    plt.cla()
    plt.scatter(xs, ys, s=1)
    plt.title("X vs Y for RandInt")
    plt.savefig("X_V_Y_{}.png".format(len(random_catalog_coords)))
    plt.cla()

    # values = random_catalog_coords.to_pixel(wcs=wcs)
    # r_c_x, r_c_y = np.fliplr(np.flipud(values))
    # random_catalog_coords = SkyCoord.from_pixel(r_c_x, r_c_y, wcs=wcs, origin=1)
    # ax.scatter(r_c_x, r_c_y, c='r')
    # ax.imshow(mask)
    # x,y = real_catalog.to_pixel(wcs=wcs)
    # ax.scatter(x,y, c='b')

    ax.scatter(random_catalog_coords.ra.degree, random_catalog_coords.dec.degree, c='r', s=1,
               transform=ax.get_transform('world'))
    ax.scatter(real_catalog.ra.degree, real_catalog.dec.degree, c='b', s=3, transform=ax.get_transform('world'))
    fig.show()
    # random_catalog_coords = np.asarray(r_c_x)
    return random_catalog_coords, np.asarray(r_c_x)


# exit()


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
        # print(element)
        sep2d = element.separation(data_catalog[i + 1:]).arcsecond
        # print(sep2d)
        if data_data is None:
            data_data = sep2d
        else:
            data_data = np.concatenate((data_data, sep2d))
    min_dist = np.min(data_data)
    print("Min Distance: {}".format(min_dist))
    min_dist = 8.39
    # min_dist = 18.75 # For Negative Ones
    max_dist = np.max(data_data)

    random_random = None

    # Get it for each one that is not the current ones
    for i, element in enumerate(random_catalog):
        sep2d = element.separation(random_catalog[i + 1:]).arcsecond
        if random_random is None:
            random_random = sep2d
        else:
            random_random = np.concatenate((random_random, sep2d))

    # Plot distribution of those to make sure random

    m_dist = np.max(random_random)
    if m_dist > max_dist:
        max_dist = m_dist

    data_random = None

    # Get it for each one that is not the current ones
    for i, element in enumerate(data_catalog):
        sep2d = element.separation(random_catalog).arcsecond
        if data_random is None:
            data_random = sep2d
        else:
            data_random = np.concatenate((data_random, sep2d))

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
    data_array_norm = (real_catalog.shape[0] * (real_catalog.shape[0] - 1)) / 2.
    data_random_array_norm = (real_catalog.shape[0] * random_catalog.shape[0])
    random_array_norm = (random_catalog.shape[0] * (random_catalog.shape[0] - 1)) / 2.

    data_array = data_array / data_array_norm
    data_random_array = data_random_array / data_random_array_norm
    random_array = random_array / random_array_norm

    # return 2 * (5000/real_catalog.shape[0])*(data_array/data_random_array) - 1

    # return 4 * (data_array*random_array)/(data_random_array)**2 - 1
    # print("Data-Data: {}".format(data_array))
    # print("RR: {}".format(random_array))
    # print("DR: {}".format(data_random_array))
    # print("DD/RR: {}".format(data_array/random_array))
    # print("DR/RR: {}".format(data_random_array/random_array))
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

    lower = [0.000, 0.173, 0.708, 1.367, 2.086, 2.840, 3.620, 4.419, 5.232, 6.057, 6.891, 7.734, 8.585, 9.441, 10.30,
             11.17,
             12.04, 12.92, 13.80, 14.68, 15.57, 16.45, 17.35, 18.24, 19.14, 20.03, 20.93, 21.84, 22.74, 23.65, 24.55]

    upper = [1.841, 3.300, 4.638, 5.918, 7.163, 8.382, 9.584, 10.77, 11.95, 13.11, 14.27, 15.42, 16.56, 17.70, 18.83,
             19.96,
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
    return a * (x ** (-0.8))


from scipy.optimize import leastsq

fitfunc = lambda a, x: a * (x ** (-0.8))
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

pinit = [1.0]


def cross_correlation_method(co_lines, random_catalog, galaxies, num_gals, num_co, num_random):
    """
    Need it in Skycoords, to do angular correlation function for each of them
    :param co_lines:
    :param random:
    :param galaxies:
    :param num_gals:
    :param num_co:
    :param num_random:
    :return:
    """
    front = num_random / num_gals

    _, co_gal, _, min_dist, max_dist = angular_correlation_function(co_lines, galaxies)
    _, co_random, _, min_dist, max_dist = angular_correlation_function(co_lines, random_catalog)

    co_gal_array_norm = (galaxies.shape[0] * co_lines.shape[0])
    co_random_array_norm = (random_catalog.shape[0] * co_lines.shape[0])
    co_gal = co_gal / co_gal_array_norm
    co_random = co_random / co_random_array_norm

    return front * (co_gal / co_random) - 1


from scipy.optimize import curve_fit

negative = True
num_points = 20000
if negative:
    real_catalog = load_table("line_search_N3_wa_crop.out")
else:
    real_catalog = load_table("line_search_P3_wa_crop.out")
real_catalog = real_catalog[real_catalog['rsnrrbin'] > 6.25]
real_catalog = make_skycoords(real_catalog, ra='rra', dec='rdc')
# Trim galaxy coordinates
#gal_catalog, _ = trim_galaxy(gal_catalog, "/media/jacob/A6548D38548D0BED/gs_A1_2chn.fits")
np.random.seed(5227)
random_catalog, r_pixels = generate_random_catalog(num_points, "/media/jacob/A6548D38548D0BED/gs_A1_2chn.fits")
# np.random.seed(5227)
random_catalog2, r2_pixels = generate_random_catalog(num_points, "/media/jacob/A6548D38548D0BED/gs_A1_2chn.fits")
# Now Write out the coordinates for the random ones
# convert to a table of ra and dec
ra = []
dec = []
t = Table()
for coord in random_catalog:
    ra.append(coord.ra.degree)
    dec.append(coord.dec.degree)

t['ra'] = Column(ra, unit='degree', description='RA')
t['dec'] = Column(dec, unit='degree', description='DEC')
ascii.write(t, "random_catalog60.txt")
# ascii.write(random_catalog2, "random_catalog2.txt")
# exit()

"""
data_data, data_random, random_random, min_dist, max_dist = angular_correlation_function(random_catalog,
                                                                                         random_catalog2)
from astroML.correlation import two_point

for bin_num in [5, 6, 7, 8, 9, 10]:
    distance_bins = np.logspace(np.log10(min_dist - 0.001), np.log10(max_dist + 1), bin_num)
    distance_bins = np.concatenate((np.asarray([0.0]), distance_bins))
    dd, bins = histogram(data_data, bins=distance_bins)
    dr, _ = histogram(data_random, bins=distance_bins)
    rr, _ = histogram(random_random, bins=distance_bins)
    x_vals = []
    omega_w = xi_r(dd, dr, rr, random_catalog, random_catalog2)
    le_omega_w, ue_omega_w = xi_r_error(omega_w, dd)
    distance_bins1 = 0.5 * (distance_bins[1:] + distance_bins[:-1])
    # distance_bins1 = np.concatenate((np.asarray([1.5]), distance_bins1))
    # Best fit to the data
    pinit = [1.]
    out = leastsq(errfunc, pinit,
                  args=(distance_bins1, omega_w, (le_omega_w + ue_omega_w) / 2.), full_output=True)
    a = out[0][0]
    s_sq = (errfunc(out[0][0], distance_bins1, omega_w, (le_omega_w + ue_omega_w) / 2.) ** 2).sum() / (
                len(distance_bins1) - len(pinit))
    cov_matrix = out[1] * s_sq

    a_error = np.absolute(cov_matrix[0][0]) ** 0.5
    #    print("Value for A: {}".format(a))
    xmin = min_dist
    xmax = max_dist + 0.1 * max_dist
    x_fit = np.linspace(xmin, xmax, 10000)

    plt.cla()
    plt.errorbar(x=distance_bins1, y=omega_w, yerr=(le_omega_w, ue_omega_w), fmt='o')
    # plt.plot(x_fit, correlation_function(x_fit, params[0][0]), '-', c='r', label='Fit Normal: A = {}'.format(np.round( params[0][0],4)))
    plt.plot(x_fit, correlation_function(x_fit, a), '-', c='r',
             label='Fit: A = {}+-{}'.format(np.round(a, 4), np.round(a_error, 5)))
    # plt.scatter(x=distance_bins1, y=corr, c='r', label='Normal')
    # plt.scatter(x=distance_bins1, y=corr2, c='g', label='Other Random')
    plt.legend(loc='best')
    plt.xscale("log")
    plt.title("Random vs Random")
    plt.xlabel("Angular Distance (arcseconds)")
    plt.ylabel("$\omega(\\theta)$")
    # plt.yscale("log")
    # plt.tight_layout()
    plt.savefig("final/Random_vs_Random_{}NoParenFlip_bin{}.png".format(num_points, 60), dpi=300)

    plt.cla()
    plt.errorbar(x=distance_bins1, y=omega_w, yerr=(le_omega_w, ue_omega_w), fmt='o')
    # plt.plot(x_fit, correlation_function(x_fit, params[0][0]), '-', c='r', label='Fit Normal: A = {}'.format(np.round(a,4)))
    plt.plot(x_fit, correlation_function(x_fit, a), '-', c='r',
             label='Fit: A = {}+-{}'.format(np.round(a, 4), np.round(a_error, 5)))
    # plt.scatter(x=distance_bins1, y=corr, c='r', label='r')
    # plt.scatter(x=distance_bins1, y=corr2, c='g', label='Other Random')
    plt.legend(loc='best')
    plt.xscale("log")
    plt.title("Random vs Random")
    plt.xlabel("Angular Distance (arcseconds)")
    plt.ylabel("$\omega(\\theta)$")
    plt.yscale("log")
    # plt.tight_layout()
    plt.savefig("final/Log_Random_vs_Random_{}NoParenFlip_bin{}.png".format(num_points, 60), dpi=300)
"""
# exit()

"""
4 panel plot

Show correlation between different SN cuts

Add taking into account the error on the A fitting

"""
plt.close()
dds = {}
drs = {}
rrs = {}
grs = {}
dist_bns = {}
dist_binners = {}
snners = [6.25, 6.1, 5.9, 5.5]
snners = [60]

for sn_cut in snners:
    if negative:
        real_catalog = load_table("line_search_N3_wa_crop.out")
    else:
        real_catalog = load_table("line_search_P3_wa_crop.out")

    fidelity_sixty = ((real_catalog['width'] == 3) & (real_catalog['rsnrrbin'] >= 6.54)) | (
            (real_catalog['width'] == 5) & (real_catalog['rsnrrbin'] >= 6.83)) | \
                     ((real_catalog['width'] == 7) & (real_catalog['rsnrrbin'] >= 6.18)) | (
                             (real_catalog['width'] == 9) & (real_catalog['rsnrrbin'] >= 6.49)) | \
                     ((real_catalog['width'] == 11) & (real_catalog['rsnrrbin'] >= 6.61)) | (
                             (real_catalog['width'] == 13) & (real_catalog['rsnrrbin'] >= 6.54)) | \
                     ((real_catalog['width'] == 15) & (real_catalog['rsnrrbin'] >= 6.89)) | (
                             (real_catalog['width'] == 17) & (real_catalog['rsnrrbin'] >= 6.83)) | \
                     ((real_catalog['width'] == 19) & (real_catalog['rsnrrbin'] >= 6.1))
    # real_catalog = real_catalog[real_catalog['rsnrrbin'] > sn_cut]
    # real_catalog = real_catalog[fidelity_sixty]
    real_catalog = construct_fid_mask(real_catalog)
    real_catalog = make_skycoords(real_catalog, ra='rra', dec='rdc')
    print(real_catalog.shape)
    dds[sn_cut] = []
    drs[sn_cut] = []
    rrs[sn_cut] = []
    grs[sn_cut] = []

    data_data, data_random, random_random, min_dist, max_dist = angular_correlation_function(real_catalog,
                                                                                         random_catalog)

    for bin_num in [5, 6, 7, 8, 9, 10]:
        distance_bins = np.logspace(np.log10(min_dist - 0.001), np.log10(max_dist + 1), bin_num + 1)
        dist_bns[bin_num] = distance_bins
        dd, _ = histogram(data_data, bins=distance_bins)
        dr, _ = histogram(data_random, bins=distance_bins)
        rr, _ = histogram(random_random, bins=distance_bins)
        dds[sn_cut].append(dd)
        drs[sn_cut].append(dr)
        rrs[sn_cut].append(rr)
        omega_w = xi_r(dd, dr,  rr, real_catalog, random_catalog)
        le_omega_w, ue_omega_w = xi_r_error(omega_w, dd)
        distance_bins = 0.5 * (distance_bins[1:] + distance_bins[:-1])
        # distance_bins = np.logspace(np.log10(min_dist),np.log10(max_dist), len(omega_w))
        dist_binners[bin_num] = distance_bins
        # Best fit to the data
        pinit = [1.]
        out = leastsq(errfunc, pinit,
                      args=(distance_bins, omega_w, (le_omega_w + ue_omega_w) / 2.), full_output=True)
        a = out[0][0]
        s_sq = (errfunc(out[0][0], distance_bins, omega_w, (le_omega_w + ue_omega_w) / 2.) ** 2).sum() / (
                len(distance_bins) - len(pinit))
        cov_matrix = out[1] * s_sq

        a_error = np.absolute(cov_matrix[0][0]) ** 0.5
        print("Value for A: {}".format(a))

        # Now get one for only the positive points
        pos_mask = (omega_w > 0.)
        pos_ue = ue_omega_w[pos_mask]
        pos_le = le_omega_w[pos_mask]
        pos_bins = distance_bins[pos_mask]
        pos_omega = omega_w[pos_mask]

        pinit = [1.]
        out = leastsq(errfunc, pinit,
                      args=(pos_bins, pos_omega, (pos_le + pos_ue) / 2.), full_output=True)
        pos_a = out[0][0]
        s_sq = (errfunc(out[0][0], pos_bins, pos_omega, (pos_le + pos_ue) / 2.) ** 2).sum() / (
                len(pos_bins) - len(pinit))
        cov_matrix = out[1] * s_sq

        pos_a_error = np.absolute(cov_matrix[0][0]) ** 0.5
        print("Value for A: {}".format(a))

        xmin = distance_bins[0] - 0.001
        xmax = distance_bins[-1] + 0.1 * distance_bins[-1]
        x_fit = np.linspace(xmin, xmax, 10000)
        """
        plt.cla()
        plt.errorbar(x=distance_bins, y=omega_w, yerr=(le_omega_w, ue_omega_w), fmt='o')
        plt.plot(x_fit, correlation_function(x_fit, a), '-', c='r',
                 label='Fit: A = {}+-{}'.format(np.round(a, 4), np.round(a_error, 5)))
        plt.plot(x_fit, correlation_function(x_fit, pos_a), '--', c='g',
                 label='Pos Fit: A = {}+-{}'.format(np.round(pos_a, 4), np.round(pos_a_error, 5)))
        plt.xscale("log")
        plt.title("Data vs Random")
        plt.legend(loc='best')
        plt.xlabel("Angular Distance (arcseconds)")
        plt.ylabel("$\omega(\\theta)$")
        plt.yscale("log")
        #plt.tight_layout()
        plt.savefig("fidelity/Log_Data_vs_Random_{}_bin{}_sn{}_N{}.png".format(num_points, bin_num, sn_cut, negative), dpi=300)
        # plt.show()
        """
        r0 = calculate_r0(Angle(a * u.arcsecond), 0.8, sn8_table)
        guass_r0 = calculate_r0(Angle(a * u.arcsecond), 0.8, sn8_table, use_gauss=True)

        plt.errorbar(x=distance_bins, y=omega_w, yerr=(le_omega_w, ue_omega_w), fmt='o')
        plt.plot(x_fit, correlation_function(x_fit, a), '-', c='r',
                 label='Fit: A = {}+-{}'.format(np.round(a, 4), np.round(a_error, 5)))
        plt.plot(x_fit, correlation_function(x_fit, pos_a), '--', c='g',
                 label='Pos Fit: A = {}+-{}'.format(np.round(pos_a, 4), np.round(pos_a_error, 5)))
        plt.xscale("log")
        plt.title("Data vs Random")
        plt.legend(loc='best', fontsize='5')
        plt.xlabel("Angular Distance (arcseconds)")
        plt.ylabel("$\omega(\\theta)$")
        plt.yscale("log")
        plt.title("r0: {} Gauss r0: {}".format(np.round(r0, 3), np.round(guass_r0, 3)))
        #plt.tight_layout()
        plt.savefig("fidelity/Data_vs_Random_{}_bin{}_sn{}_N{}.png".format(num_points, bin_num, sn_cut, negative), dpi=300)
        # plt.show()
        plt.cla()


exit()


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
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='all', sharey='all', figsize=(10, 10))
    sn_cut = snners
    for index, data_data in enumerate(dd):
        if negative:
            real_catalog = load_table("line_search_N3_wa_crop.out")
        else:
            real_catalog = load_table("line_search_P3_wa_crop.out")
        if negative:
            open_file = open("negative_r0_bin{}_sn60.txt".format(len(data_data)), "a")
        else:
            open_file = open("positive_r0_bin{}_sn60.txt".format(len(data_data)), "a")
        fidelity_sixty = ((real_catalog['width'] == 3) & (real_catalog['rsnrrbin'] >= 6.54)) | (
                    (real_catalog['width'] == 5) & (real_catalog['rsnrrbin'] >= 6.83)) | \
                         ((real_catalog['width'] == 7) & (real_catalog['rsnrrbin'] >= 6.18)) | (
                                     (real_catalog['width'] == 9) & (real_catalog['rsnrrbin'] >= 6.49)) | \
                         ((real_catalog['width'] == 11) & (real_catalog['rsnrrbin'] >= 6.61)) | (
                                     (real_catalog['width'] == 13) & (real_catalog['rsnrrbin'] >= 6.54)) | \
                         ((real_catalog['width'] == 15) & (real_catalog['rsnrrbin'] >= 6.89)) | (
                                     (real_catalog['width'] == 17) & (real_catalog['rsnrrbin'] >= 6.83)) | \
                         ((real_catalog['width'] == 19) & (real_catalog['rsnrrbin'] >= 6.1))
        # real_catalog = real_catalog[real_catalog['rsnrrbin'] > sn_cut[index]]
        # real_catalog = real_catalog[fidelity_sixty]
        real_catalog = construct_fid_mask(real_catalog)
        real_catalog = make_skycoords(real_catalog, ra='rra', dec='rdc')
        # real_catalog = make_skycoords(real_catalog, ra='rra', dec='rdc')
        omega_w = xi_r(dd[index], dr[index], rr[index], real_catalog, random_catalog)
        le_omega_w, ue_omega_w = xi_r_error(omega_w, dd[index])
        pinit = [1.]
        out = leastsq(errfunc, pinit,
                      args=(distance_bins1, omega_w, (le_omega_w + ue_omega_w) / 2.), full_output=True)
        a = out[0][0]
        s_sq = (errfunc(out[0][0], distance_bins1, omega_w, (le_omega_w + ue_omega_w) / 2.) ** 2).sum() / (
                len(distance_bins1) - len(pinit))
        cov_matrix = out[1] * s_sq

        a_error = np.absolute(cov_matrix[0][0]) ** 0.5
        # Now get one for only the positive points
        pos_mask = (omega_w > 0.)
        pos_ue = ue_omega_w[pos_mask]
        pos_le = le_omega_w[pos_mask]
        pos_bins = distance_bins1[pos_mask]
        pos_omega = omega_w[pos_mask]

        pinit = [1.]
        out = leastsq(errfunc, pinit,
                      args=(pos_bins, pos_omega, (pos_le + pos_ue) / 2.), full_output=True)
        pos_a = out[0][0]
        s_sq = (errfunc(out[0][0], pos_bins, pos_omega, (pos_le + pos_ue) / 2.) ** 2).sum() / (
                len(pos_bins) - len(pinit))
        cov_matrix = out[1] * s_sq

        pos_a_error = np.absolute(cov_matrix[0][0]) ** 0.5
        xmin = min_dist - 0.001
        xmax = max_dist + 0.1 * max_dist
        x_fit = np.linspace(xmin, xmax, 10000)
        sn8_table = Table.read(
            "/home/jacob/Development/Wide_ASPECS/Final_Output/ASPECS_Line_Candidates_cleaned_all_closest_Sep_1.0_SN_fid_0.6.ecsv")  # .format(sn_cut[index]))
        if index == 0:
            open_file.write("{} SN\n".format(sn_cut[index]))
            open_file.write("A: {}, r0: {}\n".format(a, calculate_r0(Angle(a * u.arcsecond), 0.8, sn8_table)))
            open_file.write("A+: {}, r0: {}\n".format(a + a_error,
                                                      calculate_r0(Angle((a + a_error) * u.arcsecond), 0.8, sn8_table)))
            open_file.write("A-: {}, r0: {}\n".format(a - a_error,
                                                      calculate_r0(Angle((a - a_error) * u.arcsecond), 0.8, sn8_table)))
            open_file.write(
                "pos A: {}, r0: {}\n".format(pos_a, calculate_r0(Angle(pos_a * u.arcsecond), 0.8, sn8_table)))
            open_file.write("pos A+: {}, r0: {}\n".format(pos_a + pos_a_error,
                                                          calculate_r0(Angle((pos_a + pos_a_error) * u.arcsecond), 0.8,
                                                                       sn8_table)))
            open_file.write("pos A-: {}, r0: {}\n".format(pos_a - pos_a_error,
                                                          calculate_r0(Angle((pos_a - pos_a_error) * u.arcsecond), 0.8,
                                                                       sn8_table)))

            open_file.write(
                "Gauss A: {}, r0: {}\n".format(a, calculate_r0(Angle(a * u.arcsecond), 0.8, sn8_table, use_gauss=True)))
            open_file.write("Gauss A+: {}, r0: {}\n".format(a + a_error,
                                                            calculate_r0(Angle((a + a_error) * u.arcsecond), 0.8,
                                                                         sn8_table, use_gauss=True)))
            open_file.write("Gauss A-: {}, r0: {}\n".format(a - a_error,
                                                            calculate_r0(Angle((a - a_error) * u.arcsecond), 0.8,
                                                                         sn8_table, use_gauss=True)))
            open_file.write("Gauss pos A: {}, r0: {}\n".format(pos_a,
                                                               calculate_r0(Angle(pos_a * u.arcsecond), 0.8, sn8_table,
                                                                            use_gauss=True)))
            open_file.write("Gauss pos A+: {}, r0: {}\n".format(pos_a + pos_a_error,
                                                                calculate_r0(Angle((pos_a + pos_a_error) * u.arcsecond),
                                                                             0.8, sn8_table, use_gauss=True)))
            open_file.write("Gauss pos A-: {}, r0: {}\n".format(pos_a - pos_a_error,
                                                                calculate_r0(Angle((pos_a - pos_a_error) * u.arcsecond),
                                                                             0.8, sn8_table, use_gauss=True)))
            ax1.errorbar(x=distance_bins1, y=omega_w, yerr=(le_omega_w, ue_omega_w), fmt='o')
            ax1.plot(x_fit, correlation_function(x_fit, a), '-', c='r',
                     label='Fit: A = {}+-{}'.format(np.round(a, 4), np.round(a_error, 5)))
            ax1.plot(x_fit, correlation_function(x_fit, pos_a), '--', c='g',
                     label='Pos Fit: A = {}+-{}'.format(np.round(pos_a, 4), np.round(pos_a_error, 5)))
            ax1.set_xscale("log")
            ax1.set_title("S/N > {}".format(snners[index]))
            # plt.title("Data vs Random")
            ax1.legend(loc='best', fontsize='8')
            ax1.set_ylabel("$\omega(\\theta)$")
            # ax1.set_ylabel("$\omega(\\theta)$")
            if use_log:
                ax1.set_yscale("log")
            # plt.savefig("final/Log_Data_vs_Random_{}_bin{}_sn{}.png".format(num_points, bin_num, sn_cut), dpi=300)
        if index == 1:
            open_file.write("{} SN\n".format(sn_cut[index]))
            open_file.write("A: {}, r0: {}\n".format(a, calculate_r0(Angle(a * u.arcsecond), 0.8, sn8_table)))
            open_file.write("A+: {}, r0: {}\n".format(a + a_error,
                                                      calculate_r0(Angle((a + a_error) * u.arcsecond), 0.8, sn8_table)))
            open_file.write("A-: {}, r0: {}\n".format(a - a_error,
                                                      calculate_r0(Angle((a - a_error) * u.arcsecond), 0.8, sn8_table)))
            open_file.write(
                "pos A: {}, r0: {}\n".format(pos_a, calculate_r0(Angle(pos_a * u.arcsecond), 0.8, sn8_table)))
            open_file.write("pos A+: {}, r0: {}\n".format(pos_a + pos_a_error,
                                                          calculate_r0(Angle((pos_a + pos_a_error) * u.arcsecond), 0.8,
                                                                       sn8_table)))
            open_file.write("pos A-: {}, r0: {}\n".format(pos_a - pos_a_error,
                                                          calculate_r0(Angle((pos_a - pos_a_error) * u.arcsecond), 0.8,
                                                                       sn8_table)))

            open_file.write(
                "Gauss A: {}, r0: {}\n".format(a, calculate_r0(Angle(a * u.arcsecond), 0.8, sn8_table, use_gauss=True)))
            open_file.write("Gauss A+: {}, r0: {}\n".format(a + a_error,
                                                            calculate_r0(Angle((a + a_error) * u.arcsecond), 0.8,
                                                                         sn8_table, use_gauss=True)))
            open_file.write("Gauss A-: {}, r0: {}\n".format(a - a_error,
                                                            calculate_r0(Angle((a - a_error) * u.arcsecond), 0.8,
                                                                         sn8_table, use_gauss=True)))
            open_file.write("Gauss pos A: {}, r0: {}\n".format(pos_a,
                                                               calculate_r0(Angle(pos_a * u.arcsecond), 0.8, sn8_table,
                                                                            use_gauss=True)))
            open_file.write("Gauss pos A+: {}, r0: {}\n".format(pos_a + pos_a_error,
                                                                calculate_r0(Angle((pos_a + pos_a_error) * u.arcsecond),
                                                                             0.8, sn8_table, use_gauss=True)))
            open_file.write("Gauss pos A-: {}, r0: {}\n".format(pos_a - pos_a_error,
                                                                calculate_r0(Angle((pos_a - pos_a_error) * u.arcsecond),
                                                                             0.8, sn8_table, use_gauss=True)))
            ax2.errorbar(x=distance_bins1, y=omega_w, yerr=(le_omega_w, ue_omega_w), fmt='o')
            ax2.plot(x_fit, correlation_function(x_fit, a), '-', c='r',
                     label='Fit: A = {}+-{}'.format(np.round(a, 4), np.round(a_error, 5)))
            ax2.plot(x_fit, correlation_function(x_fit, pos_a), '--', c='g',
                     label='Pos Fit: A = {}+-{}'.format(np.round(pos_a, 4), np.round(pos_a_error, 5)))
            ax2.set_xscale("log")
            ax2.set_title("S/N > {}".format(snners[index]))
            # plt.title("Data vs Random")
            ax2.legend(loc='best', fontsize='8')
            # plt.xlabel("Angular Distance (arcseconds)")
            # plt.ylabel("$\omega(\\theta)$")
            if use_log:
                ax2.set_yscale("log")
            # plt.savefig("final/Log_Data_vs_Random_{}_bin{}_sn{}.png".format(num_points, bin_num, sn_cut), dpi=300)
        if index == 2:
            open_file.write("{} SN\n".format(sn_cut[index]))
            open_file.write("A: {}, r0: {}\n".format(a, calculate_r0(Angle(a * u.arcsecond), 0.8, sn8_table)))
            open_file.write("A+: {}, r0: {}\n".format(a + a_error,
                                                      calculate_r0(Angle((a + a_error) * u.arcsecond), 0.8, sn8_table)))
            open_file.write("A-: {}, r0: {}\n".format(a - a_error,
                                                      calculate_r0(Angle((a - a_error) * u.arcsecond), 0.8, sn8_table)))
            open_file.write(
                "pos A: {}, r0: {}\n".format(pos_a, calculate_r0(Angle(pos_a * u.arcsecond), 0.8, sn8_table)))
            open_file.write("pos A+: {}, r0: {}\n".format(pos_a + pos_a_error,
                                                          calculate_r0(Angle((pos_a + pos_a_error) * u.arcsecond), 0.8,
                                                                       sn8_table)))
            open_file.write("pos A-: {}, r0: {}\n".format(pos_a - pos_a_error,
                                                          calculate_r0(Angle((pos_a - pos_a_error) * u.arcsecond), 0.8,
                                                                       sn8_table)))

            open_file.write(
                "Gauss A: {}, r0: {}\n".format(a, calculate_r0(Angle(a * u.arcsecond), 0.8, sn8_table, use_gauss=True)))
            open_file.write("Gauss A+: {}, r0: {}\n".format(a + a_error,
                                                            calculate_r0(Angle((a + a_error) * u.arcsecond), 0.8,
                                                                         sn8_table, use_gauss=True)))
            open_file.write("Gauss A-: {}, r0: {}\n".format(a - a_error,
                                                            calculate_r0(Angle((a - a_error) * u.arcsecond), 0.8,
                                                                         sn8_table, use_gauss=True)))
            open_file.write("Gauss pos A: {}, r0: {}\n".format(pos_a,
                                                               calculate_r0(Angle(pos_a * u.arcsecond), 0.8, sn8_table,
                                                                            use_gauss=True)))
            open_file.write("Gauss pos A+: {}, r0: {}\n".format(pos_a + pos_a_error,
                                                                calculate_r0(Angle((pos_a + pos_a_error) * u.arcsecond),
                                                                             0.8, sn8_table, use_gauss=True)))
            open_file.write("Gauss pos A-: {}, r0: {}\n".format(pos_a - pos_a_error,
                                                                calculate_r0(Angle((pos_a - pos_a_error) * u.arcsecond),
                                                                             0.8, sn8_table, use_gauss=True)))
            ax3.errorbar(x=distance_bins1, y=omega_w, yerr=(le_omega_w, ue_omega_w), fmt='o')
            ax3.plot(x_fit, correlation_function(x_fit, a), '-', c='r',
                     label='Fit: A = {}+-{}'.format(np.round(a, 4), np.round(a_error, 5)))
            ax3.plot(x_fit, correlation_function(x_fit, pos_a), '--', c='g',
                     label='Pos Fit: A = {}+-{}'.format(np.round(pos_a, 4), np.round(pos_a_error, 5)))
            ax3.set_xscale("log")
            ax3.set_title("S/N > {}".format(snners[index]))
            # plt.title("Data vs Random")
            ax3.legend(loc='best', fontsize='8')
            ax3.set_xlabel("Angular Distance (arcseconds)")
            ax3.set_ylabel("$\omega(\\theta)$")
            if use_log:
                ax3.set_yscale("log")
            # plt.savefig("final/Log_Data_vs_Random_{}_bin{}_sn{}.png".format(num_points, bin_num, sn_cut), dpi=300)
        if index == 3:
            open_file.write("{} SN\n".format(sn_cut[index]))
            open_file.write("A: {}, r0: {}\n".format(a, calculate_r0(Angle(a * u.arcsecond), 0.8, sn8_table)))
            open_file.write("A+: {}, r0: {}\n".format(a + a_error,
                                                      calculate_r0(Angle((a + a_error) * u.arcsecond), 0.8, sn8_table)))
            open_file.write("A-: {}, r0: {}\n".format(a - a_error,
                                                      calculate_r0(Angle((a - a_error) * u.arcsecond), 0.8, sn8_table)))
            open_file.write(
                "pos A: {}, r0: {}\n".format(pos_a, calculate_r0(Angle(pos_a * u.arcsecond), 0.8, sn8_table)))
            open_file.write("pos A+: {}, r0: {}\n".format(pos_a + pos_a_error,
                                                          calculate_r0(Angle((pos_a + pos_a_error) * u.arcsecond), 0.8,
                                                                       sn8_table)))
            open_file.write("pos A-: {}, r0: {}\n".format(pos_a - pos_a_error,
                                                          calculate_r0(Angle((pos_a - pos_a_error) * u.arcsecond), 0.8,
                                                                       sn8_table)))

            open_file.write(
                "Gauss A: {}, r0: {}\n".format(a, calculate_r0(Angle(a * u.arcsecond), 0.8, sn8_table, use_gauss=True)))
            open_file.write("Gauss A+: {}, r0: {}\n".format(a + a_error,
                                                            calculate_r0(Angle((a + a_error) * u.arcsecond), 0.8,
                                                                         sn8_table, use_gauss=True)))
            open_file.write("Gauss A-: {}, r0: {}\n".format(a - a_error,
                                                            calculate_r0(Angle((a - a_error) * u.arcsecond), 0.8,
                                                                         sn8_table, use_gauss=True)))
            open_file.write("Gauss pos A: {}, r0: {}\n".format(pos_a,
                                                               calculate_r0(Angle(pos_a * u.arcsecond), 0.8, sn8_table,
                                                                            use_gauss=True)))
            open_file.write("Gauss pos A+: {}, r0: {}\n".format(pos_a + pos_a_error,
                                                                calculate_r0(Angle((pos_a + pos_a_error) * u.arcsecond),
                                                                             0.8, sn8_table, use_gauss=True)))
            open_file.write("Gauss pos A-: {}, r0: {}\n".format(pos_a - pos_a_error,
                                                                calculate_r0(Angle((pos_a - pos_a_error) * u.arcsecond),
                                                                             0.8, sn8_table, use_gauss=True)))
            ax4.errorbar(x=distance_bins1, y=omega_w, yerr=(le_omega_w, ue_omega_w), fmt='o')
            ax4.plot(x_fit, correlation_function(x_fit, a), '-', c='r',
                     label='Fit: A = {}+-{}'.format(np.round(a, 4), np.round(a_error, 5)))
            ax4.plot(x_fit, correlation_function(x_fit, pos_a), '--', c='g',
                     label='Pos Fit: A = {}+-{}'.format(np.round(pos_a, 4), np.round(pos_a_error, 5)))
            ax4.set_xscale("log")
            ax4.set_title("S/N > {}".format(snners[index]))
            # plt.title("Data vs Random")
            ax4.legend(loc='best', fontsize='8')
            # plt.xlabel("Angular Distance (arcseconds)")
            if use_log:
                ax4.set_yscale("log")
            ax4.set_xlabel("Angular Distance (arcseconds)")
            # plt.savefig("final/Log_Data_vs_Random_{}_bin{}_sn{}.png".format(num_points, bin_num, sn_cut), dpi=300)
            # Now add the plot to the others
            ax1.plot(x_fit, correlation_function(x_fit, a), '--', c='orange', label='S/N>8. Fit')
            ax2.plot(x_fit, correlation_function(x_fit, a), '--', c='orange', label='S/N>8. Fit')
            ax3.plot(x_fit, correlation_function(x_fit, a), '--', c='orange', label='S/N>8. Fit')

    # Actually plot it all now
    f.align_xlabels()
    f.align_ylabels()
    f.suptitle("Data vs Random")
    if use_log:
        plt.savefig("final/Log_4Panel_Data_Vs_Random_bin{}_N{}_sn{}.png".format(len(distance_bins1), negative, sn_cut), dpi=300)
    else:
        plt.savefig("final/4Panel_Data_Vs_Random_bin{}_N{}_sn{}.png".format(len(distance_bins1), negative, sn_cut), dpi=300)

    open_file.close()


num_ones = len(dds[snners[0]])

for i in range(num_ones):
    dd_plot = []
    dr_plot = []
    rr_plot = []

    for key in dds.keys():
        dd_plot.append(dds[key][i])
        dr_plot.append(drs[key][i])
        rr_plot.append(rrs[key][i])

    # Number of bins
    bin_num = len(dd_plot[0])
    distance_bins = dist_bns[bin_num]
    distance_bins1 = dist_binners[bin_num]
    plot_four(dd_plot, dr_plot, rr_plot, distance_bins, distance_bins1)
    plot_four(dd_plot, dr_plot, rr_plot, distance_bins, distance_bins1, use_log=False)
