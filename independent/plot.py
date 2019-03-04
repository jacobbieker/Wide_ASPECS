"""

This is focused on the plotting, both of cutouts, SN, and others

"""
import numpy as np
import matplotlib.pyplot as plt


def plot_signal_to_noise():
    return NotImplementedError


def plot_single_cutout():
    return NotImplementedError


def plot_multiple_cutouts():
    return NotImplementedError


def plot_mstar_vs_sfr(type='square'):
    return NotImplementedError


def plot_leinerdt():
    return NotImplementedError


def whitaker_main_sequence(z, mass_start, mass_end):
    """
    Returns the main sequence fit for Mstar between mass_start and mass_end for a given Z from
    Whitaker et al. 2014

    :param z:
    :param mass_start:
    :param mass_end:
    :return:
    """
    mass_range = np.linspace(mass_start, mass_end, 100)

    if z > 2.0: # Valid till 2.5
        sfr = -19.99 + 3.44 * mass_range + -0.13 * mass_range ** 2
    elif 1.5 < z <= 2.0:
        sfr = -24.04 + 4.17 * mass_range + -0.16 * mass_range ** 2
    elif 1.0 < z <= 1.5:
        sfr = -26.03 + 4.62 * mass_range + -0.19 * mass_range ** 2
    elif z <= 1.0: #Valid to 0.5
        sfr = -27.4 + 5.02 * mass_range + -0.22 * mass_range ** 2

    return mass_range, sfr


def schrieber_main_sequence(z, mass_start, mass_end):
    """"
    Returns the x and y values for the Schreiber C. et al. 2015 Main Sequence for galaxies

    Because M* and SFR are given in log10 space anyway, need to do 10^value for linspace

    But not for SFR output, because it equals log10(SFR) which is what is given
    """

    r = np.log10(1 + z)

    mass_range = np.linspace(mass_start, mass_end, 100)

    print(mass_start)
    print(np.min(mass_range))
    print(mass_end)
    print(np.max(mass_range))

    m_0 = 0.5  # +-0.07
    a_0 = 1.5  # +- 0.15
    a_1 = 0.3  # +- 0.08
    m_1 = 0.36  # +- 0.3
    a_2 = 2.5  # +- 0.6

    m_s = mass_range - 9.0
    print(np.min(m_s))
    print(np.max(m_s))

    sfr = []
    for m in m_s:
        sfr.append((m - m_0 + a_0 * r - a_1 * (np.max([0, (m - m_1 - a_2 * r)]) ** 2)))

    return mass_range, sfr


def create_points_and_error_by_z(column_base_name, initial_catalog, low_z, high_z):
    z_mask = (initial_catalog["z"] >= low_z) & (initial_catalog["z"] <= high_z)
    centerpoints = initial_catalog[z_mask][str(column_base_name + "_50")]
    lower_error = initial_catalog[z_mask][str(column_base_name + "_16")]
    upper_error = initial_catalog[z_mask][str(column_base_name + "_84")]
    z_values = initial_catalog[z_mask]["z"]
    centerpoints = np.nan_to_num(centerpoints)
    lower_error = np.nan_to_num(lower_error)
    upper_error = np.nan_to_num(upper_error)
    zero_mask = centerpoints != 0.0
    centerpoints = centerpoints[zero_mask]
    lower_error = centerpoints - lower_error[zero_mask]
    upper_error = upper_error[zero_mask] - centerpoints
    z_values = z_values[zero_mask]
    error_bars = [lower_error, upper_error]
    return centerpoints, error_bars, z_values


def create_points_and_error(column_base_name, initial_catalog):
    centerpoints = initial_catalog[str(column_base_name + "_50")]
    lower_error = initial_catalog[str(column_base_name + "_16")]
    upper_error = initial_catalog[str(column_base_name + "_84")]
    centerpoints = np.nan_to_num(centerpoints)
    zero_mask = centerpoints != 0.0
    centerpoints = centerpoints[zero_mask]
    lower_error = centerpoints - lower_error[zero_mask]
    upper_error = upper_error[zero_mask] - centerpoints
    error_bars = [lower_error, upper_error]
    print(error_bars[0].shape)

    return centerpoints, error_bars
