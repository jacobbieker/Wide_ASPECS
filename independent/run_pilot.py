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
roberto_catalog = Table.read("/home/jacob/Development/Wide_ASPECS/data/jacob_aspecs_catalog_fixed_magphys_jcb3.fits", format='fits')
combined_catalog = combine_catalogs(initial_catalog, roberto_catalog)

aspecs_lines = load_table("ASPECS_pilot.txt")
matched_locs = load_table("ASPECS_Pilot_Matches.txt")
matched_coords = make_skycoords(matched_locs, ra='rra', dec='rdc')
aspecs_coords = make_skycoords(aspecs_lines, ra='rra', dec='rdc')
#idxc, idxcatalog, d2d, d3d = search_around_sky(matched_coords, aspecs_coords, 1.0 * u.arcsecond)
# special_ones = ["/home/jacob/Development/Wide_ASPECS/independent/matches/sn59_sep15.fits", "/home/jacob/Development/Wide_ASPECS/independent/matches/sn6_sep15.fits"]

for option_name in ["all_closest"]:
    for separation in [2.0]:
        for snrlimit in [0.0, ]:
            aspecs_table, aspecs_catalog = match_lines_to_catalog_pilot(aspecs_lines, combined_catalog,
                                                                        method=option_name, max_sep=separation,
                                                                        max_redshift=0.3)
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
            new_cat = make_skycoords(aspecs_catalog, ra='ra', dec='dc')
            # Get the new ones
            idx, d2d, d3d = match_coordinates_sky(matched_coords, new_cat)
            aspecs_catalog = aspecs_catalog[idx]
            #
            # Do the new plotting
            plt.errorbar(aspecs_lines['mstar'], aspecs_lines['sfr'],
                         yerr=(aspecs_lines['sl_error'], aspecs_lines['su_error']),
                         xerr=(aspecs_lines['ml_error'], aspecs_lines['mu_error']), fmt='.',
                         label='Pilot', c='g')
            '''
                        plt.errorbar(np.log10(aspecs_lines['mstar'])+9, np.log10(aspecs_lines['sfr']),
                         yerr=(np.log10(aspecs_lines['sl_error']), np.log10(aspecs_lines['su_error'])),
                         xerr=(np.log10(aspecs_lines['ml_error']), np.log10(aspecs_lines['mu_error'])), fmt='.',
                         label='Pilot', c='g')
            plt.errorbar(aspecs_catalog['Mstar_50_1'], aspecs_catalog['SFR_50_1'], yerr=(
                (aspecs_catalog['SFR_50_1'] - aspecs_catalog['SFR_16_1']),
                (aspecs_catalog['SFR_84_1'] - aspecs_catalog['SFR_50_1'])
            ),
                         xerr=((aspecs_catalog['Mstar_50_1'] - aspecs_catalog['Mstar_16_1']),
                               (aspecs_catalog['Mstar_84_1'] - aspecs_catalog['Mstar_50_1'])
                               ), fmt='.',
                         label='MAGPHYS')
            '''
            plt.errorbar(10**aspecs_catalog['Mstar_50_2']/(10**9), 10**aspecs_catalog['SFR_50_2'], yerr=(
                (10**aspecs_catalog['SFR_50_2'] - 10**aspecs_catalog['SFR_16_2']),
                (10**aspecs_catalog['SFR_84_2'] - 10**aspecs_catalog['SFR_50_2'])
            ),
                         xerr=((10**aspecs_catalog['Mstar_50_2']/(10**9) - 10**aspecs_catalog['Mstar_16_2']/(10**9)),
                               (10**aspecs_catalog['Mstar_84_2']/(10**9) - 10**aspecs_catalog['Mstar_50_2']/(10**9))
                               ), fmt='.',
                         label='MAGPHYS')
            #plt.errorbar(aspecs_catalog['Mstar_50_2'], aspecs_catalog['SFR_50_2'], yerr=(
            #    (aspecs_catalog['SFR_50_2'] - aspecs_catalog['SFR_16_2']),
            #(aspecs_catalog['SFR_84_2'] - aspecs_catalog['SFR_50_2'])
            #),
            #             xerr=((aspecs_catalog['Mstar_50_2'] - aspecs_catalog['Mstar_16_2']),
            #                   (aspecs_catalog['Mstar_84_2'] - aspecs_catalog['Mstar_50_2'])
            #                   ), fmt='.',
            #             label='MAGPHYS (Original)')
            plt.xlabel("Mstar (x10^9 Msun)")
            plt.ylabel("SFR")
            plt.title("Mstar vs SFR Blind ASPECS Pilot")
            # plt.ylim(0.0,3.2)
            # plt.xlim(8.0,12.)
            for i, element in enumerate(aspecs_lines):
                plt.annotate(element['name'], (aspecs_lines[i]['mstar'], aspecs_lines[i]['sfr']))
            names = ["31","32","33","35"]
            for i, row in enumerate(aspecs_catalog):
                plt.annotate(names[i], (10**row['Mstar_50_1']/(10**9), 10**row['SFR_50_1']))
                print("Mstar: {} SFR: {} Z: {}".format(10**row['Mstar_50_1']/(10**9), 10**row['SFR_50_1'], row['z_1']))
                print("Mstar2: {} SF2: {} Z2: {}".format(10**row['Mstar_50_2']/(10**9), 10**row['SFR_50_2'], row['z_2']))
            #plt.annotate("34", (33.113,16.8655))
            #plt.annotate("32", (407.3802, 135.5189))
            #plt.annotate("31", (2.79898, 32.885))
            #plt.annotate("33", (0.00238, 0.010046))

            for index, line in enumerate(idx):
                #print("\nLine: {} RA: {} DEC: {}".format(line['name'], aspecs_coords[index].ra.to_string(unit=u.hour, sep=':'), aspecs_coords[index].dec.to_string(unit=u.degree, sep=':')))
                print("Matched Line: RA: {} DEC: {}\n".format(new_cat[idx[index]].ra.to_string(unit=u.hour, sep=':'), new_cat[idx[index]].dec.to_string(unit=u.degree, sep=':')))
            plt.legend(loc='best')
            #plt.show()
            plt.savefig("Comparison_Blind_Pilot_Sources_2.png", dpi=300)
            # combined_catalog = perform_cuts(combined_catalog)
            # plot_mstar_vs_sfr(aspecs_catalog, combined_catalog, snr_limit=6.,z_lows=(0.0, 1.0, 2.0, 4), z_highs=(1.0, 2., 4, 9.), filename="/home/jacob/Development/Wide_ASPECS/Pilot/Pilot_Mstar_vs_SFR_{}_sep_{}_sn_{}.png".format(option_name, separation, snrlimit))
