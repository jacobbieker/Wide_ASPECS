import numpy as np
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube

cubes = ["A1", "A2"]

number_sn_bins = 65

all_sn_summed = np.asarray([0. for i in range(number_sn_bins-1)])
all_fids = []
for i in range(len(cubes)):
    total_rms = 0.0
    cube = SpectralCube.read("/home/jacob/Research/Wide_ASPECS/Data/gs_{}_2chn.fits".format(cubes[i]))
    sn_summed = np.asarray([0. for i in range(number_sn_bins-1)])
    sn_bins = None
    print(cube)

    num_bins = 200
    # Now go through each slice and do it for all of them

    sn_max = 0.0
    sn_min = 0.0

    avg_fidelity_above_65 = 0.0
    avg_fidelity_above_59 = 0.0
    avg_fidelity_above_53 = 0.0
    count_65 = 0
    count_59 = 0
    count_53 = 0
    for slice in range(480):
        sub_cube = cube.unitless.unmasked_data[slice,:,:]
        rms_noise = np.nanstd(sub_cube)
        total_rms += rms_noise
        #print(rms_noise)

    print("Mean RMS per channel: {}".format(total_rms/480.))