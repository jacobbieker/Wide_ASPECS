from match_sources import match_lines_to_catalog, make_skycoords, match_lines_to_catalog_pilot
from plot import plot_mstar_vs_sfr, plot_mstar_vs_sfr_third
from data import load_table, combine_catalogs, perform_cuts, save_catalog
import astropy.io.ascii as ascii
from astropy.coordinates import match_coordinates_sky, search_around_sky
import astropy.units as u


import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, hstack, join

initial_catalog = Table.read("/home/jacob/Development/Wide_ASPECS/independent/jacob_mapghys_in_nov2018_all_jcb4_magphys_jcb4.fits", format='fits')  # hdu_list[1].data
roberto_catalog = Table.read("/home/jacob/Development/Wide_ASPECS/data/jacob_aspecs_catalog_fixed_magphys_jcb2.fits", format='fits')
combined_catalog = combine_catalogs(initial_catalog, roberto_catalog)

aspecs_lines = load_table("ASPECS_Pilot_C_Matches.txt")
#matched_locs = load_table("ASPECS_Pilot_Matches.txt")
matched_coords = make_skycoords(combined_catalog, ra='ra', dec='dc')
aspecs_coords = make_skycoords(aspecs_lines, ra='rra', dec='rdc')
idxc, idxcatalog, d2d, d3d  = search_around_sky(matched_coords, aspecs_coords, 1.0 * u.arcsecond)
print(d2d)
properties = load_table("ASPECS_Pilot_C_Properties.txt")

#print(aspecs_lines[idxcatalog])
#print(combined_catalog[idxc]['id', 'Mstar_50_1', 'Mstar_50_2', 'SFR_50_1', 'SFR_50_2', 'z_1', 'z_2'])
matches = {}
matches_2 = {}

for element in properties:
    matches[element['name']] = (0,0,-1,0,0,0,0)
    matches_2[element['name']] = (0,0,-1,0,0,0,0)

for index, row in enumerate(aspecs_lines[idxcatalog]):
    print("\nID: {} \nZ: {} \nMstar_1: {} SFR_1: {} \n Mstar_2: {} SFR_2: {}".format(row['name'], np.round(combined_catalog[idxc[index]]['z_2'], 3), (combined_catalog[idxc[index]]['Mstar_50_1']),
                                                                                     combined_catalog[idxc[index]]['SFR_50_1'], (combined_catalog[idxc[index]]['Mstar_50_2']),
                                                                                     combined_catalog[idxc[index]]['SFR_50_2']))
    for prop in properties:
        if prop['name'] == row['name']:
            if np.abs(prop['sfr'] - combined_catalog[idxc[index]]['SFR_50_1']) < np.abs(prop['sfr'] - matches[row['name']][1]):
                if np.abs(prop['mstar'] - combined_catalog[idxc[index]]['Mstar_50_1']) < np.abs(prop['mstar'] - matches[row['name']][0]):
                    matches[row['name']] = (combined_catalog[idxc[index]]['Mstar_50_1'], combined_catalog[idxc[index]]['SFR_50_1'], idxc[index],
                                            (combined_catalog[idxc[index]]['SFR_50_1']-combined_catalog[idxc[index]]['SFR_16_1']),
                                            (combined_catalog[idxc[index]]['SFR_84_1']-combined_catalog[idxc[index]]['SFR_50_1']),
                                            (combined_catalog[idxc[index]]['Mstar_50_1']-combined_catalog[idxc[index]]['Mstar_16_1']),
                                            (combined_catalog[idxc[index]]['Mstar_84_1']-combined_catalog[idxc[index]]['Mstar_50_1']))
            if np.abs(prop['sfr'] - combined_catalog[idxc[index]]['SFR_50_2']) < np.abs(prop['sfr'] - matches_2[row['name']][1]):
                if np.abs(prop['mstar'] - combined_catalog[idxc[index]]['Mstar_50_2']) < np.abs(prop['mstar'] - matches_2[row['name']][0]):
                    matches_2[row['name']] = (combined_catalog[idxc[index]]['Mstar_50_2'], combined_catalog[idxc[index]]['SFR_50_2'], idxc[index],
                                              (combined_catalog[idxc[index]]['SFR_50_2']-combined_catalog[idxc[index]]['SFR_16_2']),
                                              (combined_catalog[idxc[index]]['SFR_84_2']-combined_catalog[idxc[index]]['SFR_50_2']),
                                              (combined_catalog[idxc[index]]['Mstar_50_2']-combined_catalog[idxc[index]]['Mstar_16_2']),
                                              (combined_catalog[idxc[index]]['Mstar_84_2']-combined_catalog[idxc[index]]['Mstar_50_2']))

# Get the closest matches to those values and the closest ones

# now have them together
one_sfr = []
one_mstar = []
one_mu_error = []
one_ml_error = []
one_su_error = []
one_sl_error = []
two_sfr = []
two_mstar = []
t_mu_error = []
t_ml_error = []
t_su_error = []
t_sl_error = []
labels = []

for key in matches.keys():
    labels.append(key)
    one_sfr.append(matches[key][1])
    one_mstar.append(matches[key][0])
    one_sl_error.append(matches[key][3])
    one_su_error.append(matches[key][4])
    one_ml_error.append(matches[key][5])
    one_mu_error.append(matches[key][6])
    two_sfr.append(matches_2[key][1])
    two_mstar.append(matches_2[key][0])
    t_sl_error.append(matches_2[key][3])
    t_su_error.append(matches_2[key][4])
    t_ml_error.append(matches_2[key][5])
    t_mu_error.append(matches_2[key][6])


plt.errorbar(properties['mstar'], properties['sfr'], yerr=(properties['sl_error'],properties['su_error']),
             xerr=(properties['ml_error'],properties['mu_error']), fmt='.', label='Pilot', c='g')
plt.errorbar(two_mstar, two_sfr, yerr=(t_sl_error,t_su_error),
             xerr=(t_ml_error,t_mu_error), fmt='.', label='MAGPHYS (MUSE)', c='b')
plt.errorbar(one_mstar, one_sfr, yerr=(one_sl_error,one_su_error),
             xerr=(one_ml_error,one_mu_error), fmt='.', label='MAGPHYS (Original)', c='orange')
plt.xlabel("Log(Mstar)")
plt.ylabel("Log(SFR)")
plt.title("Mstar vs SFR Continuum ASPECS Pilot")
#plt.ylim(0.0,3.2)
#plt.xlim(8.0,12.)
for i, element in enumerate(properties):
    plt.annotate(element['name'], (properties[i]['mstar'], properties[i]['sfr']))
for i, element in enumerate(labels):
    plt.annotate(element, (one_mstar[i], one_sfr[i]))
    plt.annotate(element, (two_mstar[i], two_sfr[i]))
plt.legend(loc='best')
#plt.savefig("Comparison_Continuum_Pilot_Sources.png", dpi=300)
plt.show()
exit()
#special_ones = ["/home/jacob/Development/Wide_ASPECS/independent/matches/sn59_sep15.fits", "/home/jacob/Development/Wide_ASPECS/independent/matches/sn6_sep15.fits"]

#combined_catalog = combine_catalogs(initial_catalog, roberto_catalog)
for option_name in ["all_closest"]:
    for separation in [1.6]:
        for snrlimit in [6.0,]:
            aspecs_table, aspecs_catalog = match_lines_to_catalog_pilot(aspecs_lines, combined_catalog, method=option_name, max_sep=separation, snr_limit=snrlimit, max_redshift=0.3)
            save_catalog(aspecs_catalog, "/home/jacob/Development/Wide_ASPECS/Pilot/aspecs_pilot_zco_3mm_catalog_SN{}_method_{}_Sep_{}".format(snrlimit, option_name, separation))
            ascii.write(aspecs_table, "/home/jacob/Development/Wide_ASPECS/Pilot/ASPECS_Pilot_Line_3mm_Candidates_{}_Sep_{}_SN_{}.txt".format(option_name, separation, snrlimit), format='fixed_width', bookend=False, delimiter=None, overwrite=True, formats={'RA (J2000)': '', 'DEC (J2000)': '%2.6f', 'Roberto RA': '%2.6f', 'Roberto DEC': '%2.6f','Observed CO (GHz)': '%3.4f', 'Restframe CO (GHz)': '%3.4f', 'Z (CO)': '%2.3f', 'Z (Matched)': '%2.3f',
                                                                                                                                                                                                                                                                'Delta Z': '%2.3f', 'Delta V (Km/s)': '%4.3f', 'Km/s': '%4.3f', 'Separation (Arcsecond)': '%2.4f', 'S/N': '%2.3f', 'Flux Density at Peak (Jy/beam)': '%2.4f',
                                                                                                                                                                                                                                                                'Integrated Flux (Jy km/s)': '%2.4f', 'Cosmic Volume (Mpc^3)': '%8.0f', 'Log(M*)': '%2.4f', 'Error Log(M*)': '%2.4f', 'Log(SFR)': '%2.4f', 'Error Log(SFR)': '%2.4f'})
            ascii.write(aspecs_table, "/home/jacob/Development/Wide_ASPECS/Pilot/ASPECS_Pilot_Line_3mm_Candidates_{}_Sep_{}_SN_{}.ecsv".format(option_name, separation, snrlimit), format='ecsv', overwrite=True)
            # Check if matched in the aspecs catalog

            #combined_catalog = perform_cuts(combined_catalog)
            #plot_mstar_vs_sfr(aspecs_catalog, combined_catalog, snr_limit=6.,z_lows=(0.0, 1.0, 2.0, 4), z_highs=(1.0, 2., 4, 9.), filename="/home/jacob/Development/Wide_ASPECS/Pilot/Pilot_Mstar_vs_SFR_{}_sep_{}_sn_{}.png".format(option_name, separation, snrlimit))