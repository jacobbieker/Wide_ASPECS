import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from spectral_cube import SpectralCube
from astropy.coordinates import SkyCoord, Angle, SkyOffsetFrame, ICRS, Distance
from astropy.table import Table, hstack, join, Column
from astropy.stats import histogram
from scipy.spatial.distance import cdist
from scipy.special import gamma
from astropy.io import ascii
from astropy.cosmology import FlatLambdaCDM
from astropy import constants as const
from scipy.interpolate import interp2d, interp1d
from scipy.stats import norm
from astropy import modeling
from scipy.optimize import leastsq, curve_fit

def load_table(ascii_table, header=0, start=1):
    ascii_table_data = Table.read(ascii_table, format="ascii", header_start=header, data_start=start)
    return ascii_table_data

def fid(neg, pos):
    return 1 - (len(neg)/len(pos))

neg_catalog = load_table("line_search_N3_wa_crop.out")
pos_catalog = load_table("line_search_P3_wa_crop.out")

line_widths = [i for i in range(3, 21, 2)]
print(line_widths)

sn_cuts = np.arange(5., 8.1, 0.01)
print(sn_cuts)
for width in line_widths:
    neg_widths = neg_catalog[neg_catalog['width'] == width]
    pos_widths = pos_catalog[pos_catalog['width'] == width]
    print("Neg Width Lines: {}".format(len(neg_widths)))
    print("Pos Width Lines: {}".format(len(pos_widths)))
    print("Width {} MaxSN: {}".format(width, np.max(neg_widths['rsnrrbin'])))
    fid_width = []
    sn_vals = []
    six_fid = -1
    for sn in sn_cuts:
        neg_sn = neg_widths[neg_widths['rsnrrbin'] >= sn]
        pos_sn = pos_widths[pos_widths['rsnrrbin'] >= sn]
        #print("SN: {} Len: {}".format(sn, len(pos_sn)))
        if len(pos_sn) > 0:
            fid_width.append(fid(neg_sn, pos_sn))
            sn_vals.append(sn)
            #print(fid_width[-1])
            if six_fid < 0 and fid_width[-1] >= 0.6:
                six_fid = sn
        elif len(neg_sn) == 0:
            fid_width.append(1)
            sn_vals.append(sn)
            if six_fid < 0 and fid_width[-1] >= 0.6:
                six_fid = sn


    # Now plot it
    plt.plot(sn_vals, fid_width, label="{}".format(width))
    #plt.axvline(x=six_fid, c='black', ls="--")
    print("Width {} 0.6 Fid: {}".format(width, six_fid))

plt.title("Fidelity vs SN")
plt.axhline(y=0.6, c='r', ls='--')
plt.ylabel("Fidelity")
plt.xlabel("S/N")
plt.legend(loc="best")
plt.show()