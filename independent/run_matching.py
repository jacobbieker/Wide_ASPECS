from match_sources import match_lines_to_catalog
from plot import plot_mstar_vs_sfr, plot_mstar_vs_sfr_third, plot_mstar_vs_sfr_specz
from data import load_table, combine_catalogs, perform_cuts, save_catalog
import astropy.io.ascii as ascii


import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, hstack, join, Column

t = Table.read(
    "/home/jacob/Development/Wide_ASPECS/Final_Output/No_Cut_ASPECS_Line_Candidates_all_closest_Sep_1.0_SN_fid_60.ecsv")

t.sort("S/N")
t.reverse()
print(ascii.write(t, format='latex'))
#ascii.write(t, "Ordered_Fid60.ecsv", format="ecsv")
#exit()

t = t['RA (J2000)', 'DEC (J2000)', 'Observed CO (GHz)',
      'Transition', 'Z (CO)',
      'Delta Z', 'Spec Z', 'S/N', 'Cosmic Volume (Mpc^3)', 'Log(M*)', 'Roberto ID',
      'Separation (Arcsecond)', 'Log(SFR)', 'Z (Matched)', 'Error Log(SFR)', 'Error Log(M*)' ]

counterparts = (t["Roberto ID"] > 0.0)
t = t[counterparts]

t = t['RA (J2000)', 'DEC (J2000)', 'Observed CO (GHz)',
      'Z (CO)', 'Z (Matched)',
      'Delta Z', 'Spec Z', 'S/N', 'Separation (Arcsecond)', 'Log(SFR)', 'Log(M*)', 'Error Log(SFR)', 'Error Log(M*)' ]

ids = ["ID.{}".format(i+1) for i in range(len(t))]
t.sort("S/N")
t.reverse()
aa = Column(ids, name='ID')
t.add_column(aa, index=0)
#counterparts = (t["Roberto ID"] > 0.0)
print(ascii.write(t, format='latex'))
exit()
plt.hist(t['Z (CO)'], bins=10)
plt.title("Redshift Distribution of CO Line Candidates")
plt.ylabel("Number of Sources")
plt.xlabel("Redshift")
plt.savefig("redshift_catalog.png", dpi=300)

exit()

initial_catalog = Table.read("/home/jacob/Development/Wide_ASPECS/independent/jacob_mapghys_in_nov2018_all_jcb4_magphys_jcb4.fits", format='fits')  # hdu_list[1].data
roberto_catalog = Table.read("roberto_catalog_muse_skelton_matched_manFix.fits", format='fits')

negative = False
if negative:
    aspecs_lines = load_table("line_search_N3_wa_crop.out")
else:
    aspecs_lines = load_table("line_search_P3_wa_crop.out")

special_ones = ["/home/jacob/Development/Wide_ASPECS/independent/matches/sn59_sep15.fits", "/home/jacob/Development/Wide_ASPECS/independent/matches/sn6_sep15.fits"]


combined_catalog = combine_catalogs(initial_catalog, roberto_catalog)
for option_name in ["all_closest"]:
    for separation in [1.0]:
        for snrlimit in ["fid_40"]:
            aspecs_table, aspecs_catalog, spec_catalog, no_spec_catalog = match_lines_to_catalog(aspecs_lines, combined_catalog, method=option_name, max_sep=separation, snr_limit=snrlimit)
            save_catalog(aspecs_catalog, "/home/jacob/Development/Wide_ASPECS/Final_Output/No_Cut_aspecs_zco_catalog_SN{}_method_{}_Sep_{}".format(snrlimit, option_name, separation))
            ascii.write(aspecs_table, "/home/jacob/Development/Wide_ASPECS/Final_Output/No_Cut_ASPECS_Line_Candidates_Only_Matched_{}_Sep_{}_SN_{}.txt".format(option_name, separation, snrlimit), format='fixed_width', bookend=False, delimiter=None, formats={'RA (J2000)': '%2.6f', 'DEC (J2000)': '%2.6f', 'Roberto RA': '%2.6f', 'Roberto DEC': '%2.6f','Observed CO (GHz)': '%3.4f', 'Restframe CO (GHz)': '%3.4f', 'Z (CO)': '%2.3f', 'Z (Matched)': '%2.3f',
                                                                                                                                 'Delta Z': '%2.3f', 'Delta V (Km/s)': '%4.3f', 'Km/s': '%4.3f', 'Separation (Arcsecond)': '%2.4f', 'S/N': '%2.3f', 'Flux Density at Peak (Jy/beam)': '%2.4f',
                                                                                                                                 'Integrated Flux (Jy km/s)': '%2.4f', 'Cosmic Volume (Mpc^3)': '%8.0f', 'Log(M*)': '%2.4f', 'Error Log(M*)': '%2.4f', 'Log(SFR)': '%2.4f', 'Error Log(SFR)': '%2.4f'})
            ascii.write(aspecs_table, "/home/jacob/Development/Wide_ASPECS/Final_Output/No_Cut_ASPECS_Line_Candidates_{}_Sep_{}_SN_{}.ecsv".format(option_name, separation, snrlimit), format='ecsv')
            combined_catalog = perform_cuts(combined_catalog)
            #exit()
            plot_mstar_vs_sfr(spec_catalog, combined_catalog, snr_limit=snrlimit, filename="/home/jacob/Development/Wide_ASPECS/Final_Output/No_Cut_OnlySpecZ_Mstar_vs_SFR_{}_sep_{}_sn_{}.png".format(option_name, separation, snrlimit))
            plot_mstar_vs_sfr_specz(spec_catalog, combined_catalog, no_spec_catalog, snr_limit=snrlimit, filename="/home/jacob/Development/Wide_ASPECS/Final_Output/No_Cut_Mstar_vs_SFR_{}_sep_{}_sn_{}.png".format(option_name, separation, snrlimit))
            aspecs_table, aspecs_catalog, spec_catalog, no_spec_catalog = match_lines_to_catalog(aspecs_lines, combined_catalog, method=option_name, max_sep=separation, snr_limit=snrlimit)
            save_catalog(aspecs_catalog, "/home/jacob/Development/Wide_ASPECS/Final_Output/No_Cut_aspecs_zco_cleaned_catalog_SN{}_method_{}_Sep_{}".format(snrlimit, option_name, separation))
            ascii.write(aspecs_table, "/home/jacob/Development/Wide_ASPECS/Final_Output/No_Cut_ASPECS_Line_Candidates_cleaned_{}_Sep_{}_SN_{}.txt".format(option_name, separation, snrlimit), format='fixed_width', bookend=False, delimiter=None, formats={'RA (J2000)': '%2.6f', 'DEC (J2000)': '%2.6f', 'Roberto RA': '%2.6f', 'Roberto DEC': '%2.6f','Observed CO (GHz)': '%3.4f', 'Restframe CO (GHz)': '%3.4f', 'Z (CO)': '%2.3f', 'Z (Matched)': '%2.3f',
                                                                                                                                                                                            'Delta Z': '%2.3f', 'Delta V (Km/s)': '%4.3f', 'Km/s': '%4.3f', 'Separation (Arcsecond)': '%2.4f', 'S/N': '%2.3f', 'Flux Density at Peak (Jy/beam)': '%2.4f',
                                                                                                                                                                                            'Integrated Flux (Jy km/s)': '%2.4f', 'Cosmic Volume (Mpc^3)': '%8.0f', 'Log(M*)': '%2.4f', 'Error Log(M*)': '%2.4f', 'Log(SFR)': '%2.4f', 'Error Log(SFR)': '%2.4f'})
            ascii.write(aspecs_table, "/home/jacob/Development/Wide_ASPECS/Final_Output/No_Cut_ASPECS_Line_Candidates_cleaned_{}_Sep_{}_SN_{}.ecsv".format(option_name, separation, snrlimit), format='ecsv')
            plot_mstar_vs_sfr_specz(spec_catalog, combined_catalog, no_spec_catalog, snr_limit=snrlimit, filename="/home/jacob/Development/Wide_ASPECS/Final_Output/No_Cut_Mstar_vs_SFR_cleaned_{}_sep_{}_sn_{}.png".format(option_name, separation, snrlimit))
            temp = aspecs_table['RA (J2000)', 'DEC (J2000)', 'Observed CO (GHz)',
                                'Transition', 'Z (Matched)', 'Z (CO)',
                                'Delta Z', 'Delta V (Km/s)', 'S/N']
            combined_catalog = combine_catalogs(initial_catalog, roberto_catalog)

