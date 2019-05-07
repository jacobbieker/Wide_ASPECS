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


def plot_mstar_vs_sfr_third(aspecs_catalog, only_aspecs_match_catalog, matched_catalog, snr_limit, max_z=0.3, labels=('All', 'ASPECS', 'ASPECS CO\n Magphys'), z_lows=(0.0, 1.1, 2.2, 3), z_highs=(0.4, 1.8, 3, 4.4), colors=('lightgrey', 'red', 'orange', 'green'), filename="", z_name="z_co", added_thing="", type='square'):
    """
    Given a set of catalogs, plot them with labels in a M* by SFR overlaid with best fits

    Does perform cuts and errors on the catalogs


    :param z_lows:
    :param z_highs:
    :param type:
    :return:
    """

    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='all', sharey='all')

    for index, z_range in enumerate(zip(z_lows, z_highs)):
        use_labels = False
        if index == 0:
            ax = ax1
            use_labels = True
        elif index == 1:
            ax = ax2
            use_labels = False
        elif index == 2:
            ax = ax3
            use_labels = False
        else:
            ax = ax4
            use_labels  = False
        # Calculate and make the Main Sequence Plots
        # Do it for the mid point of the range
        mid_range_z = (z_range[1] + z_range[0]) / 2.
        whitaker_dotted = False
        whitaker_mass, whitaker_sfr = whitaker_main_sequence(mid_range_z, 6., 12.)
        schrieber_mass, schrieber_sfr = schrieber_main_sequence(mid_range_z, 6., 12.)
        if mid_range_z < 0.5 or mid_range_z > 2.5:
            whitaker_dotted = True
        if use_labels:
            ax.plot(schrieber_mass, schrieber_sfr, color='green', label='S15', zorder=20)
            if whitaker_dotted:
                ax.plot(whitaker_mass, whitaker_sfr, color='orange', label='W14', linestyle='dashed',
                        zorder=20)
            else:
                ax.plot(whitaker_mass, whitaker_sfr, color='orange', label='W14', zorder=20)

        ax.plot(schrieber_mass, schrieber_sfr, color='green', zorder=20)
        if whitaker_dotted:
            ax.plot(whitaker_mass, whitaker_sfr, color='orange', linestyle='dashed', zorder=20)
        else:
            ax.plot(whitaker_mass, whitaker_sfr, color='orange', zorder=20)
        sfr, sfr_error, sfr_z = create_points_and_error_by_z("SFR", matched_catalog, z_range[0], z_range[1], added="_1")
        mstar, mstar_error, mstar_z = create_points_and_error_by_z("Mstar", matched_catalog, z_range[0], z_range[1], added="_1")
        if use_labels:
            ax.errorbar(mstar, sfr, yerr=sfr_error, xerr=mstar_error, ecolor=colors[0], label=labels[0],  mec='darkgrey',  fmt='.', ms=1, elinewidth=1)
        else:
            ax.errorbar(mstar, sfr, yerr=sfr_error, xerr=mstar_error, ecolor=colors[0], fmt='.', mec='darkgrey', ms=1, elinewidth=1)

        sfr, sfr_error, sfr_z = create_points_and_error_by_z("SFR", aspecs_catalog, z_range[0], z_range[1], z=z_name, added="_1")
        mstar, mstar_error, mstar_z = create_points_and_error_by_z("Mstar", aspecs_catalog, z_range[0], z_range[1], z=z_name, added="_1")
        if use_labels:
            ax.errorbar(mstar, sfr, yerr=sfr_error, xerr=mstar_error, ecolor=colors[1], label=labels[1], fmt='.', ms=5, mec='red', zorder=20, elinewidth=1)
        else:
            ax.errorbar(mstar, sfr, yerr=sfr_error, xerr=mstar_error, ecolor=colors[1], fmt='.', ms=5, mec='red', zorder=20, elinewidth=1)

        sfr, sfr_error, sfr_z = create_points_and_error_by_z("SFR", only_aspecs_match_catalog, z_range[0], z_range[1], added=added_thing, z=z_name)
        mstar, mstar_error, mstar_z = create_points_and_error_by_z("Mstar", only_aspecs_match_catalog, z_range[0], z_range[1], added=added_thing, z=z_name)
        if use_labels:
            ax.errorbar(mstar, sfr, yerr=sfr_error, xerr=mstar_error, ecolor='blue', label=labels[2], fmt='.', ms=5, mec='blue', zorder=20, elinewidth=1)
        else:
            ax.errorbar(mstar, sfr, yerr=sfr_error, xerr=mstar_error, ecolor='blue', fmt='.', ms=5, mec='blue', zorder=20, elinewidth=1)


        ax.set_title(str(np.round(z_range[0], 1)) + ' < Z < ' + str(np.round(z_range[1], 1)))

        if use_labels:
            handles, ax_labels = ax.get_legend_handles_labels()
            f.legend(loc='best', handles=handles, labels=ax_labels, prop={'size': 6})
    f.text(0.5, 0.01, 'Log(M*)', ha='center')
    f.text(0.01, 0.5, 'Log(SFR)', va='center', rotation='vertical')
    f.savefig(filename, bbox_inches='tight', dpi=300)
    f.show()

def plot_mstar_vs_sfr(aspecs_catalog, matched_catalog, snr_limit, max_z=0.3, labels=('All', 'ASPECS'), z_lows=(0.0, 1.1, 2.2, 3), z_highs=(0.4, 1.8, 3, 4.4), colors=('lightgrey', 'red', 'orange', 'green'), filename="", type='square'):
    """
    Given a set of catalogs, plot them with labels in a M* by SFR overlaid with best fits

    Does perform cuts and errors on the catalogs


    :param z_lows:
    :param z_highs:
    :param type:
    :return:
    """

    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='all', sharey='all')

    for index, z_range in enumerate(zip(z_lows, z_highs)):
        use_labels = False
        if index == 0:
            ax = ax1
            use_labels = True
        elif index == 1:
            ax = ax2
            use_labels = False
        elif index == 2:
            ax = ax3
            use_labels = False
        else:
            ax = ax4
            use_labels  = False
        # Calculate and make the Main Sequence Plots
        # Do it for the mid point of the range
        mid_range_z = (z_range[1] + z_range[0]) / 2.
        whitaker_dotted = False
        whitaker_mass, whitaker_sfr = whitaker_main_sequence(mid_range_z, 6., 12.)
        schrieber_mass, schrieber_sfr = schrieber_main_sequence(mid_range_z, 6., 12.)
        if mid_range_z < 0.5 or mid_range_z > 2.5:
            whitaker_dotted = True
        if use_labels:
            ax.plot(schrieber_mass, schrieber_sfr, color='green', label='S15', zorder=20)
            if whitaker_dotted:
                ax.plot(whitaker_mass, whitaker_sfr, color='orange', label='W14', linestyle='dashed',
                        zorder=20)
            else:
                ax.plot(whitaker_mass, whitaker_sfr, color='orange', label='W14', zorder=20)

        ax.plot(schrieber_mass, schrieber_sfr, color='green', zorder=20)
        if whitaker_dotted:
            ax.plot(whitaker_mass, whitaker_sfr, color='orange', linestyle='dashed', zorder=20)
        else:
            ax.plot(whitaker_mass, whitaker_sfr, color='orange', zorder=20)
        sfr, sfr_error, sfr_z = create_points_and_error_by_z("SFR", matched_catalog, z_range[0], z_range[1])
        mstar, mstar_error, mstar_z = create_points_and_error_by_z("Mstar", matched_catalog, z_range[0], z_range[1])
        if use_labels:
            ax.errorbar(mstar, sfr, yerr=sfr_error, xerr=mstar_error, ecolor=colors[0], label=labels[0],  mec='darkgrey',  fmt='.', ms=1, elinewidth=1)
        else:
            ax.errorbar(mstar, sfr, yerr=sfr_error, xerr=mstar_error, ecolor=colors[0], fmt='.', mec='darkgrey', ms=1, elinewidth=1)

        sfr, sfr_error, sfr_z = create_points_and_error_by_z("SFR", aspecs_catalog, z_range[0]-max_z, z_range[1]+max_z)
        mstar, mstar_error, mstar_z = create_points_and_error_by_z("Mstar", aspecs_catalog, z_range[0]-max_z, z_range[1]+max_z)
        if use_labels:
            ax.errorbar(mstar, sfr, yerr=sfr_error, xerr=mstar_error, ecolor=colors[1], label=labels[1], fmt='.', ms=5, mec='red', zorder=20, elinewidth=1)
        else:
            ax.errorbar(mstar, sfr, yerr=sfr_error, xerr=mstar_error, ecolor=colors[1], fmt='.', ms=5, mec='red', zorder=20, elinewidth=1)

        ax.set_title(str(np.round(z_range[0], 1)) + ' < Z < ' + str(np.round(z_range[1], 1)))

        if use_labels:
            handles, ax_labels = ax.get_legend_handles_labels()
            f.legend(loc='best', handles=handles, labels=ax_labels, prop={'size': 6})
    f.text(0.5, 0.01, 'Log(M*)', ha='center')
    f.text(0.01, 0.5, 'Log(SFR)', va='center', rotation='vertical')
    f.savefig(filename, bbox_inches='tight', dpi=300)
    f.show()


def plot_mstar_vs_sfr_specz(spec_z_catalog, matched_catalog, no_spec_z_catalog, snr_limit, max_z=0.3, labels=('All', 'Spec Z', 'ASPECS'), z_lows=(0.0, 1.1, 2.2, 3), z_highs=(0.4, 1.8, 3, 4.4), colors=('lightgrey', 'red', 'blue', 'orange', 'green'), filename="", type='square'):
    """
    Given a set of catalogs, plot them with labels in a M* by SFR overlaid with best fits

    Does perform cuts and errors on the catalogs


    :param z_lows:
    :param z_highs:
    :param type:
    :return:
    """

    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='all', sharey='all')

    for index, z_range in enumerate(zip(z_lows, z_highs)):
        use_labels = False
        if index == 0:
            ax = ax1
            use_labels = True
        elif index == 1:
            ax = ax2
            use_labels = False
        elif index == 2:
            ax = ax3
            use_labels = False
        else:
            ax = ax4
            use_labels  = False
        # Calculate and make the Main Sequence Plots
        # Do it for the mid point of the range
        mid_range_z = (z_range[1] + z_range[0]) / 2.
        whitaker_dotted = False
        whitaker_mass, whitaker_sfr = whitaker_main_sequence(mid_range_z, 6., 12.)
        schrieber_mass, schrieber_sfr = schrieber_main_sequence(mid_range_z, 6., 12.)
        if mid_range_z < 0.5 or mid_range_z > 2.5:
            whitaker_dotted = True
        if use_labels:
            ax.plot(schrieber_mass, schrieber_sfr, color='green', label='S15', zorder=20)
            if whitaker_dotted:
                ax.plot(whitaker_mass, whitaker_sfr, color='orange', label='W14', linestyle='dashed',
                        zorder=20)
            else:
                ax.plot(whitaker_mass, whitaker_sfr, color='orange', label='W14', zorder=20)

        ax.plot(schrieber_mass, schrieber_sfr, color='green', zorder=20)
        if whitaker_dotted:
            ax.plot(whitaker_mass, whitaker_sfr, color='orange', linestyle='dashed', zorder=20)
        else:
            ax.plot(whitaker_mass, whitaker_sfr, color='orange', zorder=20)

        sfr, sfr_error, sfr_z = create_points_and_error_by_z("SFR", matched_catalog, z_range[0], z_range[1])
        mstar, mstar_error, mstar_z = create_points_and_error_by_z("Mstar", matched_catalog, z_range[0], z_range[1])
        if use_labels:
            ax.errorbar(mstar, sfr, yerr=sfr_error, xerr=mstar_error, ecolor=colors[0], label=labels[0],  mec='darkgrey',  fmt='.', ms=1, elinewidth=1)
        else:
            ax.errorbar(mstar, sfr, yerr=sfr_error, xerr=mstar_error, ecolor=colors[0], fmt='.', mec='darkgrey', ms=1, elinewidth=1)

        sfr, sfr_error, sfr_z = create_points_and_error_by_z("SFR", spec_z_catalog, z_range[0] - max_z, z_range[1] + max_z)
        mstar, mstar_error, mstar_z = create_points_and_error_by_z("Mstar", spec_z_catalog, z_range[0] - max_z, z_range[1] + max_z)

        if use_labels:
            ax.errorbar(mstar, sfr, yerr=sfr_error, xerr=mstar_error, ecolor=colors[1], label=labels[1], fmt='.', ms=5, mec='red', zorder=20, elinewidth=1)
        else:
            ax.errorbar(mstar, sfr, yerr=sfr_error, xerr=mstar_error, ecolor=colors[1], fmt='.', ms=5, mec='red', zorder=20, elinewidth=1)

        sfr, sfr_error, sfr_z = create_points_and_error_by_z("SFR", no_spec_z_catalog, z_range[0] - max_z, z_range[1] + max_z)
        mstar, mstar_error, mstar_z = create_points_and_error_by_z("Mstar", no_spec_z_catalog, z_range[0] - max_z, z_range[1] + max_z)

        if use_labels:
            ax.errorbar(mstar, sfr, yerr=sfr_error, xerr=mstar_error, ecolor=colors[2], label=labels[2], fmt='.', ms=5, mec='blue', zorder=20, elinewidth=1)
        else:
            ax.errorbar(mstar, sfr, yerr=sfr_error, xerr=mstar_error, ecolor=colors[2], fmt='.', ms=5, mec='blue', zorder=20, elinewidth=1)

        ax.set_title(str(np.round(z_range[0], 1)) + ' < Z < ' + str(np.round(z_range[1], 1)))

        if use_labels:
            handles, ax_labels = ax.get_legend_handles_labels()
            f.legend(loc='best', handles=handles, labels=ax_labels, prop={'size': 6})
    f.text(0.5, 0.01, 'Log(M*)', ha='center')
    f.text(0.01, 0.5, 'Log(SFR)', va='center', rotation='vertical')
    f.savefig(filename, bbox_inches='tight', dpi=300)
    f.show()

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

    m_0 = 0.5  # +-0.07
    a_0 = 1.5  # +- 0.15
    a_1 = 0.3  # +- 0.08
    m_1 = 0.36  # +- 0.3
    a_2 = 2.5  # +- 0.6

    m_s = mass_range - 9.0

    sfr = []
    for m in m_s:
        sfr.append((m - m_0 + a_0 * r - a_1 * (np.max([0, (m - m_1 - a_2 * r)]) ** 2)))

    return mass_range, sfr


def create_points_and_error_by_z(column_base_name, initial_catalog, low_z, high_z, added="_1", z="z_1"):
    z_mask = (initial_catalog[z] >= low_z) & (initial_catalog[z] <= high_z)
    centerpoints = initial_catalog[z_mask][str(column_base_name + "_50"+added)]
    lower_error = initial_catalog[z_mask][str(column_base_name + "_16"+added)]
    upper_error = initial_catalog[z_mask][str(column_base_name + "_84"+added)]
    z_values = initial_catalog[z_mask][z]
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
    centerpoints = initial_catalog[str(column_base_name + "_50_1")]
    lower_error = initial_catalog[str(column_base_name + "_16_1")]
    upper_error = initial_catalog[str(column_base_name + "_84_1")]
    centerpoints = np.nan_to_num(centerpoints)
    zero_mask = centerpoints != 0.0
    centerpoints = centerpoints[zero_mask]
    lower_error = centerpoints - lower_error[zero_mask]
    upper_error = upper_error[zero_mask] - centerpoints
    error_bars = [lower_error, upper_error]

    return centerpoints, error_bars
