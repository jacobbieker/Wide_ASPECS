from match_sources import match_lines_to_catalog, make_skycoords, match_lines_to_catalog_pilot
from plot import plot_mstar_vs_sfr, plot_mstar_vs_sfr_third
from data import load_table, combine_catalogs, perform_cuts, save_catalog
import astropy.io.ascii as ascii
from astropy.coordinates import match_coordinates_sky, search_around_sky
import astropy.units as u

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, hstack, join

initial_catalog = Table.read(
    "/home/jacob/Development/Wide_ASPECS/independent/jacob_mapghys_in_nov2018_all_jcb4_magphys_jcb4.fits",
    format='fits')  # hdu_list[1].data
roberto_catalog = Table.read("roberto_catalog_muse_skelton_matched_manFix.fits", format='fits')
combined_catalog = combine_catalogs(initial_catalog, roberto_catalog)

aspecs_lines = load_table("ASPECS_pilot.txt")
# matched_locs = load_table("ASPECS_Pilot_Matches.txt")
matched_coords = make_skycoords(combined_catalog, ra='ra', dec='dc')
aspecs_coords = make_skycoords(aspecs_lines, ra='rra', dec='rdc')
idxc, idxcatalog, d2d, d3d = search_around_sky(matched_coords, aspecs_coords, 1.0 * u.arcsecond)
# special_ones = ["/home/jacob/Development/Wide_ASPECS/independent/matches/sn59_sep15.fits", "/home/jacob/Development/Wide_ASPECS/independent/matches/sn6_sep15.fits"]

combined_catalog = combine_catalogs(initial_catalog, roberto_catalog)
for option_name in ["all_closest"]:
    for separation in [1.6]:
        for snrlimit in [6.0, ]:
            aspecs_table, aspecs_catalog = match_lines_to_catalog_pilot(aspecs_lines, combined_catalog,
                                                                        method=option_name, max_sep=separation,
                                                                        snr_limit=snrlimit, max_redshift=0.3)
            save_catalog(aspecs_catalog,
                         "/home/jacob/Development/Wide_ASPECS/Pilot/aspecs_pilot_zco_3mm_catalog_SN{}_method_{}_Sep_{}".format(
                             snrlimit, option_name, separation))
            ascii.write(aspecs_table,
                        "/home/jacob/Development/Wide_ASPECS/Pilot/ASPECS_Pilot_Line_3mm_Candidates_{}_Sep_{}_SN_{}.txt".format(
                            option_name, separation, snrlimit), format='fixed_width', bookend=False, delimiter=None,
                        overwrite=True, formats={'RA (J2000)': '', 'DEC (J2000)': '%2.6f', 'Roberto RA': '%2.6f',
                                                 'Roberto DEC': '%2.6f', 'Observed CO (GHz)': '%3.4f',
                                                 'Restframe CO (GHz)': '%3.4f', 'Z (CO)': '%2.3f',
                                                 'Z (Matched)': '%2.3f',
                                                 'Delta Z': '%2.3f', 'Delta V (Km/s)': '%4.3f', 'Km/s': '%4.3f',
                                                 'Separation (Arcsecond)': '%2.4f', 'S/N': '%2.3f',
                                                 'Flux Density at Peak (Jy/beam)': '%2.4f',
                                                 'Integrated Flux (Jy km/s)': '%2.4f', 'Cosmic Volume (Mpc^3)': '%8.0f',
                                                 'Log(M*)': '%2.4f', 'Error Log(M*)': '%2.4f', 'Log(SFR)': '%2.4f',
                                                 'Error Log(SFR)': '%2.4f'})
            ascii.write(aspecs_table,
                        "/home/jacob/Development/Wide_ASPECS/Pilot/ASPECS_Pilot_Line_3mm_Candidates_{}_Sep_{}_SN_{}.ecsv".format(
                            option_name, separation, snrlimit), format='ecsv', overwrite=True)
            # Check if matched in the aspecs catalog
            #
            # Do the new plotting
            plt.errorbar(np.log(aspecs_lines['mstar']*9), np.log(aspecs_lines['sfr']),
                         yerr=(np.log(aspecs_lines['sl_error']), np.log(aspecs_lines['su_error'])),
                         xerr=(np.log(aspecs_lines['ml_error']*9), np.log(aspecs_lines['mu_error']*9)), fmt='.',
                         label='Pilot', c='g')
            plt.errorbar(aspecs_catalog['Mstar_50_1'], aspecs_catalog['SFR_50_1'], yerr=(
                (aspecs_catalog['SFR_50_1'] - aspecs_catalog['SFR_16_1']),
                (aspecs_catalog['SFR_84_1'] - aspecs_catalog['SFR_50_1'])
            ),
                         xerr=((aspecs_catalog['Mstar_50_1'] - aspecs_catalog['Mstar_16_1']),
                               (aspecs_catalog['Mstar_84_1'] - aspecs_catalog['Mstar_50_1'])
                               ), fmt='.',
                         label='MAGPHYS (Original)')
            plt.errorbar(aspecs_catalog['Mstar_50_2'], aspecs_catalog['SFR_50_2'], yerr=(
                (aspecs_catalog['SFR_50_2'] - aspecs_catalog['SFR_16_2']),
            (aspecs_catalog['SFR_84_2'] - aspecs_catalog['SFR_50_2'])
            ),
                         xerr=((aspecs_catalog['Mstar_50_2'] - aspecs_catalog['Mstar_16_2']),
                               (aspecs_catalog['Mstar_84_2'] - aspecs_catalog['Mstar_50_2'])
                               ), fmt='.',
                         label='MAGPHYS (MUSE)')
            plt.xlabel("Log(Mstar)")
            plt.ylabel("Log(SFR)")
            plt.title("Mstar vs SFR Blind ASPECS Pilot")
            # plt.ylim(0.0,3.2)
            # plt.xlim(8.0,12.)
            for i, element in enumerate(aspecs_lines):
                plt.annotate(element['name'], (aspecs_lines[i]['mstar'], aspecs_lines[i]['sfr']))
            plt.legend(loc='best')
            plt.savefig("Comparison_Blind_Pilot_Sources.png", dpi=300)
            plt.show()
            # combined_catalog = perform_cuts(combined_catalog)
            # plot_mstar_vs_sfr(aspecs_catalog, combined_catalog, snr_limit=6.,z_lows=(0.0, 1.0, 2.0, 4), z_highs=(1.0, 2., 4, 9.), filename="/home/jacob/Development/Wide_ASPECS/Pilot/Pilot_Mstar_vs_SFR_{}_sep_{}_sn_{}.png".format(option_name, separation, snrlimit))
