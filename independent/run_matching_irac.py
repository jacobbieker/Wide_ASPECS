from match_sources import match_lines_to_catalog
from plot import plot_mstar_vs_sfr, plot_mstar_vs_sfr_third, plot_mstar_vs_sfr_specz
from data import load_table, combine_catalogs, perform_cuts, save_catalog
import astropy.io.ascii as ascii


import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, hstack, join

initial_catalog = Table.read("/home/jacob/Development/Wide_ASPECS/independent/jacob_mapghys_in_nov2018_all_jcb4_magphys_jcb4.fits", format='fits')  # hdu_list[1].data
roberto_catalog = Table.read("roberto_catalog_muse_skelton_matched_manFix.fits", format='fits')

aspecs_lines = load_table("/home/jacob/Development/Wide_ASPECS/independent/line_search_P3_wa_crop.out")

special_ones = ["/home/jacob/Development/Wide_ASPECS/independent/matches/sn59_sep15.fits", "/home/jacob/Development/Wide_ASPECS/independent/matches/sn6_sep15.fits"]

combined_catalog = combine_catalogs(initial_catalog, roberto_catalog)
cat_mask = (combined_catalog['fIRAC1_1'] >= (np.median(combined_catalog['fIRAC1_1']) + 1*np.std(combined_catalog['fIRAC1_1'])))

# Need to convert AB magnitude to Jansky, select those with an AB magntiude larger than
# either 17 ( low end of spectrum) , or 21.8 median
# Detection Limit is 22.8 magAB
# Decarli uses those selected galaxies
# mAB= 8.9−2.5 log(F) from : https://arxiv.org/pdf/1609.08772.pdf where F is in Jansky
# Same group assumes mag_3.6um < 26, above is undetected, potential < 25.3 as well
# mag 22.8-22.8 for m_3.6 for 3.6 um selected catalog as training sample

# Another is > 2 uJy from: https://arxiv.org/pdf/1212.1234.pdf which is about < 41 mag according to the other one
# But they say the limit for IRAC is 1 uJy which corresponds to 23.9 magAB

#print(combined_catalog['fIRAC1_1'])

def jansky_to_ab(jansky):
    return -2.5 * np.log10(jansky) + 8.9

jansky_to_ab_limit = jansky_to_ab(3631)

print(jansky_to_ab_limit)
irac_above_0 = (combined_catalog['fIRAC1_1'] > 0.0000001)
print(np.mean(combined_catalog['fIRAC1_1']))
print("Jansky Max ")
print(jansky_to_ab(np.max(combined_catalog['fIRAC1_1'])))
print("Jansky Min ")
print(jansky_to_ab(np.min(combined_catalog[irac_above_0]['fIRAC1_1'])))

converted_ones = jansky_to_ab(combined_catalog[irac_above_0]['fIRAC1_1'])

# Converted Ones
#plt.hist(jansky_to_ab(combined_catalog['fIRAC1_1']), bins=20)
#plt.show()

low_irac_mask = (converted_ones <= 19.3)
high_irac_mask = (converted_ones <= 24.5)

print(np.max(combined_catalog[irac_above_0][low_irac_mask]['fIRAC1_1']))
print(np.min(combined_catalog[irac_above_0][low_irac_mask]['fIRAC1_1']))
print(np.mean(jansky_to_ab(combined_catalog[irac_above_0][low_irac_mask]['fIRAC1_1'])))#-np.std(jansky_to_ab(combined_catalog[irac_above_0]['fIRAC1_1'])))

print("Low Mask: {}".format(len(combined_catalog[irac_above_0][low_irac_mask])))
print("High Mask: {}".format(len(combined_catalog[irac_above_0][high_irac_mask])))
print("No Mask: {}".format(len(combined_catalog)))
print("Mean Mask: {}".format(len(combined_catalog[cat_mask])))
combined_catalog = combined_catalog[irac_above_0][low_irac_mask]
exit()
for option_name in ["all_closest"]:
    for separation in [2.]:
        for snrlimit in [6., 5.85, 5.75, 5.65, 5.55, 5.45]:
            aspecs_table, aspecs_catalog, spec_catalog, no_spec_catalog = match_lines_to_catalog(aspecs_lines, combined_catalog, method=option_name, max_sep=separation, snr_limit=snrlimit)
            save_catalog(aspecs_catalog, "/home/jacob/Development/Wide_ASPECS/IRAC_Output/22.8/22.8aspecs_zco_IRAC1STD_catalog_SN{}_method_{}_Sep_{}".format(snrlimit, option_name, separation))
            ascii.write(aspecs_table, "/home/jacob/Development/Wide_ASPECS/IRAC_Output/22.8/IRAC22.8_ASPECS_Line_Candidates_Only_Matched_{}_Sep_{}_SN_{}.txt".format(option_name, separation, snrlimit), format='fixed_width', bookend=False, delimiter=None, formats={'RA (J2000)': '%2.6f', 'DEC (J2000)': '%2.6f', 'Roberto RA': '%2.6f', 'Roberto DEC': '%2.6f','Observed CO (GHz)': '%3.4f', 'Restframe CO (GHz)': '%3.4f', 'Z (CO)': '%2.3f', 'Z (Matched)': '%2.3f',
                                                                                                                                                                                                         'Delta Z': '%2.3f', 'Delta V (Km/s)': '%4.3f', 'Km/s': '%4.3f', 'Separation (Arcsecond)': '%2.4f', 'S/N': '%2.3f', 'Flux Density at Peak (Jy/beam)': '%2.4f',
                                                                                                                                                                                                         'Integrated Flux (Jy km/s)': '%2.4f', 'Cosmic Volume (Mpc^3)': '%8.0f', 'Log(M*)': '%2.4f', 'Error Log(M*)': '%2.4f', 'Log(SFR)': '%2.4f', 'Error Log(SFR)': '%2.4f'})
            ascii.write(aspecs_table, "/home/jacob/Development/Wide_ASPECS/IRAC_Output/22.8/IRAC22.8_ASPECS_Line_Candidates_{}_Sep_{}_SN_{}.ecsv".format(option_name, separation, snrlimit), format='ecsv')
            plot_mstar_vs_sfr_specz(spec_catalog, combined_catalog, no_spec_catalog, snr_limit=snrlimit, filename="/home/jacob/Development/Wide_ASPECS/IRAC_Output/22.8/IRAC22.8_Mstar_vs_SFR_{}_sep_{}_sn_{}.png".format(option_name, separation, snrlimit))
            combined_catalog = perform_cuts(combined_catalog)
            plot_mstar_vs_sfr(spec_catalog, combined_catalog, snr_limit=snrlimit, filename="/home/jacob/Development/Wide_ASPECS/IRAC_Output/22.8/IRAC22.8_0nlySpecZ_Mstar_vs_SFR_{}_sep_{}_sn_{}.png".format(option_name, separation, snrlimit))
            plot_mstar_vs_sfr_specz(spec_catalog, combined_catalog, no_spec_catalog, snr_limit=snrlimit, filename="/home/jacob/Development/Wide_ASPECS/IRAC_Output/22.8/IRAC22.8_Cleaned_Mstar_vs_SFR_{}_sep_{}_sn_{}.png".format(option_name, separation, snrlimit))

            combined_catalog = combine_catalogs(initial_catalog, roberto_catalog)
            # Only Select IRAC Sources above... Median?
            combined_catalog = combined_catalog[irac_above_0][low_irac_mask]

exit()
combined_catalog = combine_catalogs(initial_catalog, roberto_catalog)
# Only Select IRAC Sources above... Median?
combined_catalog = combined_catalog[high_irac_mask]
for option_name in ["all_closest"]:
    for separation in [2.]:
        for snrlimit in [10., 9.5, 9.0, 8.5, 8.]:
            aspecs_table, aspecs_catalog, spec_catalog, no_spec_catalog = match_lines_to_catalog(aspecs_lines, combined_catalog, method=option_name, max_sep=separation, snr_limit=snrlimit)
            save_catalog(aspecs_catalog, "/home/jacob/Development/Wide_ASPECS/IRAC_Output/22.8/22.8aspecs_zco_IRAC1STD_catalog_SN{}_method_{}_Sep_{}".format(snrlimit, option_name, separation))
            ascii.write(aspecs_table, "/home/jacob/Development/Wide_ASPECS/IRAC_Output/22.8/IRAC22.8_ASPECS_Line_Candidates_Only_Matched_{}_Sep_{}_SN_{}.txt".format(option_name, separation, snrlimit), format='fixed_width', bookend=False, delimiter=None, formats={'RA (J2000)': '%2.6f', 'DEC (J2000)': '%2.6f', 'Roberto RA': '%2.6f', 'Roberto DEC': '%2.6f','Observed CO (GHz)': '%3.4f', 'Restframe CO (GHz)': '%3.4f', 'Z (CO)': '%2.3f', 'Z (Matched)': '%2.3f',
                                                                                                                                                                                                                                                                  'Delta Z': '%2.3f', 'Delta V (Km/s)': '%4.3f', 'Km/s': '%4.3f', 'Separation (Arcsecond)': '%2.4f', 'S/N': '%2.3f', 'Flux Density at Peak (Jy/beam)': '%2.4f',
                                                                                                                                                                                                                                                                  'Integrated Flux (Jy km/s)': '%2.4f', 'Cosmic Volume (Mpc^3)': '%8.0f', 'Log(M*)': '%2.4f', 'Error Log(M*)': '%2.4f', 'Log(SFR)': '%2.4f', 'Error Log(SFR)': '%2.4f'})
            ascii.write(aspecs_table, "/home/jacob/Development/Wide_ASPECS/IRAC_Output/22.8/IRAC22.8_ASPECS_Line_Candidates_{}_Sep_{}_SN_{}.ecsv".format(option_name, separation, snrlimit), format='ecsv')
            plot_mstar_vs_sfr_specz(spec_catalog, combined_catalog, no_spec_catalog, snr_limit=snrlimit, filename="/home/jacob/Development/Wide_ASPECS/IRAC_Output/22.8/IRAC22.8_Mstar_vs_SFR_{}_sep_{}_sn_{}.png".format(option_name, separation, snrlimit))
            combined_catalog = perform_cuts(combined_catalog)
            plot_mstar_vs_sfr(spec_catalog, combined_catalog, snr_limit=snrlimit, filename="/home/jacob/Development/Wide_ASPECS/IRAC_Output/22.8/IRAC22.8_0nlySpecZ_Mstar_vs_SFR_{}_sep_{}_sn_{}.png".format(option_name, separation, snrlimit))
            plot_mstar_vs_sfr_specz(spec_catalog, combined_catalog, no_spec_catalog, snr_limit=snrlimit, filename="/home/jacob/Development/Wide_ASPECS/IRAC_Output/22.8/IRAC22.8_Cleaned_Mstar_vs_SFR_{}_sep_{}_sn_{}.png".format(option_name, separation, snrlimit))

            combined_catalog = combine_catalogs(initial_catalog, roberto_catalog)
            # Only Select IRAC Sources above... Median?
            combined_catalog = combined_catalog[high_irac_mask]