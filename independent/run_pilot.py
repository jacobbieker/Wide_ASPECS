from match_sources import match_lines_to_catalog
from plot import plot_mstar_vs_sfr, plot_mstar_vs_sfr_third
from data import load_table, combine_catalogs, perform_cuts, save_catalog
import astropy.io.ascii as ascii


import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, hstack, join

initial_catalog = Table.read("/home/jacob/Development/Wide_ASPECS/independent/jacob_mapghys_in_nov2018_all_jcb4_magphys_jcb4.fits", format='fits')  # hdu_list[1].data
roberto_catalog = Table.read("roberto_catalog_muse_skelton_matched_manFix.fits", format='fits')

aspecs_lines = load_table("ASPECS_pilot.txt")
print(aspecs_lines)
#special_ones = ["/home/jacob/Development/Wide_ASPECS/independent/matches/sn59_sep15.fits", "/home/jacob/Development/Wide_ASPECS/independent/matches/sn6_sep15.fits"]


combined_catalog = combine_catalogs(initial_catalog, roberto_catalog)
for option_name in ["all_closest"]:
    for separation in [0.25, 0.5, 1.0, 1.5]:
        for snrlimit in [6.0,]:
            aspecs_table, aspecs_catalog = match_lines_to_catalog(aspecs_lines, combined_catalog, method=option_name, max_sep=separation, snr_limit=snrlimit, max_redshift=0.3)
            save_catalog(aspecs_catalog, "/home/jacob/Development/Wide_ASPECS/Pilot/aspecs_pilot_zco_catalog_SN{}_method_{}_Sep_{}".format(snrlimit, option_name, separation))
            ascii.write(aspecs_table, "/home/jacob/Development/Wide_ASPECS/Pilot/ASPECS_Pilot_Line_Candidates_{}_Sep_{}_SN_{}.txt".format(option_name, separation, snrlimit), format='fixed_width', bookend=False, delimiter=None, formats={'Roberto RA': '%2.6f', 'Roberto DEC': '%2.6f','Observed CO (GHz)': '%3.4f', 'Restframe CO (GHz)': '%3.4f', 'Z (CO)': '%2.3f', 'Z (Matched)': '%2.3f',
                                                                                                                                                                                            'Delta Z': '%2.3f', 'Delta V (Km/s)': '%4.3f', 'Km/s': '%4.3f', 'Separation (Arcsecond)': '%2.4f', 'S/N': '%2.3f', 'Flux Density at Peak (Jy/beam)': '%2.4f',
                                                                                                                                                                                            'Integrated Flux (Jy km/s)': '%2.4f', 'Cosmic Volume (Mpc^3)': '%8.0f', 'Log(M*)': '%2.4f', 'Error Log(M*)': '%2.4f', 'Log(SFR)': '%2.4f', 'Error Log(SFR)': '%2.4f'})
            ascii.write(aspecs_table, "/home/jacob/Development/Wide_ASPECS/Pilot/ASPECS_Pilot_Line_Candidates_{}_Sep_{}_SN_{}.ecsv".format(option_name, separation, snrlimit), format='ecsv')
            #combined_catalog = perform_cuts(combined_catalog)
            #plot_mstar_vs_sfr(aspecs_catalog, combined_catalog, snr_limit=6., filename="Mstar_vs_SFR_{}_sep_{}_sn_{}.png".format(option_name, separation, snrlimit))