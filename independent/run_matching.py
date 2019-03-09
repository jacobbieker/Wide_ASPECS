from match_sources import match_lines_to_catalog
from plot import plot_mstar_vs_sfr
from data import load_table, combine_catalogs, perform_cuts
import astropy.io.ascii as ascii


import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, hstack, join

initial_catalog = Table.read("magphys_catalog.fits", format='fits')  # hdu_list[1].data
roberto_catalog = Table.read("roberto_catalog_muse_skelton_matched_manFix.fits", format='fits')

combined_catalog = combine_catalogs(initial_catalog, roberto_catalog)
aspecs_lines = load_table("line_search_P3_wa_crop.out")

aspecs_table, aspecs_catalog = match_lines_to_catalog(aspecs_lines, combined_catalog, method='all', max_sep=2., snr_limit=6., max_redshift=0.3)

print(len(aspecs_table))
print(aspecs_catalog)
print("Len Matches: {}".format(len(aspecs_catalog)))
ascii.write(aspecs_table, "ASPECS_Line_Candidates_test.txt", format='fixed_width', bookend=False, delimiter=None, formats={'RA (J2000)': '%2.6f', 'DEC (J2000)': '%2.6f', 'Roberto RA': '%2.6f', 'Roberto DEC': '%2.6f','Observed CO (GHz)': '%3.4f', 'Restframe CO (GHz)': '%3.4f', 'Z (CO)': '%2.3f', 'Z (Matched)': '%2.3f',
                                                                                                                             'Delta Z': '%2.3f', 'Delta V (Km/s)': '%4.3f', 'Km/s': '%4.3f', 'Separation (Arcsecond)': '%2.4f', 'S/N': '%2.3f', 'Flux Density at Peak (Jy/beam)': '%2.4f',
                                                                                                                             'Integrated Flux (Jy km/s)': '%2.4f', 'Cosmic Volume (Mpc^3)': '%8.0f', 'Log(M*)': '%2.4f', 'Error Log(M*)': '%2.4f', 'Log(SFR)': '%2.4f', 'Error Log(SFR)': '%2.4f'})
#exit()

combined_catalog = perform_cuts(combined_catalog)

plot_mstar_vs_sfr(aspecs_catalog, combined_catalog, snr_limit=6.)

