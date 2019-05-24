import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from spectral_cube import SpectralCube
from astropy.coordinates import SkyCoord, Angle, SkyOffsetFrame, ICRS, Distance
from astropy.table import Table, hstack, join
from astropy.stats import histogram
from scipy.spatial.distance import cdist

print(SkyCoord("3:32:41.5041", "-27:44:13.553", unit=(u.hourangle, u.deg), frame='icrs'))
print(SkyCoord("3:32:08.9056", "-27:47:02.788", unit=(u.hourangle, u.deg), frame='icrs'))
print(SkyCoord("3:32:51.1022", "-27:49:27.011", unit=(u.hourangle, u.deg), frame='icrs'))
print(SkyCoord("3:32:18.9546", "-27:52:16.837", unit=(u.hourangle, u.deg), frame='icrs'))

#exit()
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
    min_dist = 3.
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
    print("Data-Data: {}".format(data_array))
    print("RR: {}".format(random_array))
    print("DR: {}".format(data_random_array))
    print("DD/RR: {}".format(data_array/random_array))
    print("DR/RR: {}".format(data_random_array/random_array))
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
num_points = 7500
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
    distance_bins1 = np.logspace(np.log10(min_dist),np.log10(max_dist+1), len(omega_w))
    #distance_bins1 = np.concatenate((np.asarray([1.5]), distance_bins1))
    # Best fit to the data
    out = leastsq(errfunc, pinit,
                  args=(distance_bins1, omega_w, (le_omega_w+ue_omega_w)/2.), full_output=True)
    a = out[0][0]
#    print("Value for A: {}".format(a))
    xmin=min_dist
    xmax=max_dist+0.1*max_dist
    x_fit = np.linspace(xmin, xmax, 10000)

    plt.cla()
    plt.errorbar(x=distance_bins1, y=omega_w, yerr=(le_omega_w, ue_omega_w), fmt='o')
    #plt.plot(x_fit, correlation_function(x_fit, params[0][0]), '-', c='r', label='Fit Normal: A = {}'.format(np.round( params[0][0],4)))
    plt.plot(x_fit, correlation_function(x_fit, a), '-', c='g', label='Fit: A = {}'.format(np.round(a,4)))
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
    plt.plot(x_fit, correlation_function(x_fit, a), '-', c='g', label='Fit: A = {}'.format(np.round(a,4)))
    #plt.scatter(x=distance_bins1, y=corr, c='r', label='Normal')
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


for sn_cut in [9.5, 9.0, 8.5, 8.0]:
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
        distance_bins = np.logspace(np.log10(min_dist),np.log10(max_dist), len(omega_w))
        dist_binners[bin_num] = distance_bins
        # Best fit to the data
        params = curve_fit(correlation_function, distance_bins, omega_w)
        a = params[0][0]
        print("Value for A: {}".format(a))
        xmin=min_dist-0.001
        xmax=max_dist+0.1*max_dist
        x_fit = np.linspace(xmin, xmax, 10000)
        plt.cla()
        plt.errorbar(x=distance_bins, y=omega_w, yerr=(le_omega_w,ue_omega_w), fmt='o')
        plt.plot(x_fit, correlation_function(x_fit, a), '-', c='r', label='Fit: A = {}'.format(np.round(a,4)))
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
        plt.plot(x_fit, correlation_function(x_fit, a), '-', c='r', label='Fit: A = {}'.format(np.round(a,4)))
        plt.xscale("log")
        plt.title("Data vs Random")
        plt.legend(loc='best')
        plt.xlabel("Angular Distance (arcseconds)")
        plt.ylabel("$\omega(\\theta)$")
        #plt.yscale("log")
        #plt.tight_layout()
        plt.savefig("final/Data_vs_Random_{}_bin{}_sn{}.png".format(num_points, bin_num, sn_cut), dpi=300)
        #plt.show()

# Now have dictionary of lists of datas


def plot_four(dd, dr, rr, distance_bins, distance_bins1):
    """
    Take set of 4 S/N cut lists for DD, DR, and RR and plot a 4 panel plot, has to be same binning for it
    Assumes dd, dr, and rr are all in same order, 9.5, 9.0, 8.5, 8.0 SN
    :param dd:
    :param dr:
    :param rr:
    :return:
    """
    plt.cla()
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='all', sharey='all')
    sn_cut = [9.5,9.0,8.5,8.0]
    for index, data_data in enumerate(dd):
        real_catalog = load_table("line_search_P3_wa_crop.out")
        real_catalog = real_catalog[real_catalog['rsnrrbin'] > sn_cut[index]]
        real_catalog = make_skycoords(real_catalog, ra='rra', dec='rdc')
        omega_w = xi_r(dd[index], dr[index], rr[index], real_catalog, random_catalog)
        le_omega_w, ue_omega_w = xi_r_error(omega_w, dd[index])
        params = curve_fit(correlation_function, distance_bins, omega_w)
        a = params[0][0]
        print("Value for A: {}".format(a))
        xmin=min_dist-0.001
        xmax=max_dist+0.1*max_dist
        x_fit = np.linspace(xmin, xmax, 10000)
        if index == 0:
            ax1.errorbar(x=distance_bins1, y=omega_w, yerr=(le_omega_w,ue_omega_w), fmt='o')
            ax1.plot(x_fit, correlation_function(x_fit, a), '-', c='r', label='Fit: A = {}'.format(np.round(a,4)))
            ax1.xscale("log")
            ax1.set_title("S/N > 9.5")
            #plt.title("Data vs Random")
            ax1.set_legend(loc='best')
            ax1.set_xlabel("Angular Distance (arcseconds)")
            #ax1.set_ylabel("$\omega(\\theta)$")
            ax1.yscale("log")
            #plt.savefig("final/Log_Data_vs_Random_{}_bin{}_sn{}.png".format(num_points, bin_num, sn_cut), dpi=300)
        if index == 1:
            ax2.errorbar(x=distance_bins1, y=omega_w, yerr=(le_omega_w,ue_omega_w), fmt='o')
            ax2.plot(x_fit, correlation_function(x_fit, a), '-', c='r', label='Fit: A = {}'.format(np.round(a,4)))
            ax2.xscale("log")
            ax2.set_title("S/N > 9.0")
            #plt.title("Data vs Random")
            ax2.set_legend(loc='best')
            #plt.xlabel("Angular Distance (arcseconds)")
            #plt.ylabel("$\omega(\\theta)$")
            ax2.yscale("log")
            #plt.savefig("final/Log_Data_vs_Random_{}_bin{}_sn{}.png".format(num_points, bin_num, sn_cut), dpi=300)
        if index == 2:
            ax3.errorbar(x=distance_bins1, y=omega_w, yerr=(le_omega_w,ue_omega_w), fmt='o')
            ax3.plot(x_fit, correlation_function(x_fit, a), '-', c='r', label='Fit: A = {}'.format(np.round(a,4)))
            ax3.xscale("log")
            ax3.set_title("S/N > 8.5")
            #plt.title("Data vs Random")
            ax3.set_legend(loc='best')
            ax3.set_xlabel("Angular Distance (arcseconds)")
            ax3.set_ylabel("$\omega(\\theta)$")
            ax3.yscale("log")
            #plt.savefig("final/Log_Data_vs_Random_{}_bin{}_sn{}.png".format(num_points, bin_num, sn_cut), dpi=300)
        if index == 3:
            ax4.errorbar(x=distance_bins1, y=omega_w, yerr=(le_omega_w,ue_omega_w), fmt='o')
            ax4.plot(x_fit, correlation_function(x_fit, a), '-', c='r', label='Fit: A = {}'.format(np.round(a,4)))
            ax4.xscale("log")
            ax4.set_title("S/N > 8.0")
            #plt.title("Data vs Random")
            ax4.set_legend(loc='best')
            #plt.xlabel("Angular Distance (arcseconds)")
            ax4.set_ylabel("$\omega(\\theta)$")
            ax4.yscale("log")
            #plt.savefig("final/Log_Data_vs_Random_{}_bin{}_sn{}.png".format(num_points, bin_num, sn_cut), dpi=300)

    # Actually plot it all now
    f.align_xlabels()
    f.align_ylabels()
    f.sup_title("Data vs Random")
    plt.show()



num_ones = len(dds[9.5])

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
    bin_num = len(dd_plot[i])
    distance_bins = dist_bns[bin_num]
    distance_bins1 = dist_binners[bin_num]
    plot_four(dd_plot, dr_plot, rr_plot, distance_bins, distance_bins1)

"""
        dd, dr, rr, min_dist, max_dist = angular_correlation_function(real_catalog, real_catalog, number_of_bins=bin_num)
    
        omega_w = xi_r(dd, dr, rr, real_catalog, real_catalog)
        e_omega_w = xi_r_error(omega_w, dd)
        distance_bins = np.logspace(np.log10(min_dist),np.log10(max_dist), len(omega_w))
    
        plt.cla()
        plt.errorbar(x=distance_bins, y=omega_w, yerr=e_omega_w, fmt='o')
        plt.xscale("log")
        plt.title("Data vs Data")
        plt.xlabel("Angular Distance (arcseconds)")
        plt.ylabel("$\omega(\\theta)$")
        plt.yscale("log")
        plt.tight_layout()
        plt.savefig("Data_vs_Data_{}_bin{}.png".format(sn_cut, bin_num), dpi=300)
        plt.show()
"""
