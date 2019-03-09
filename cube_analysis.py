import numpy as np
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube

cubes = ["A1", "A2"]

number_sn_bins = 100
sn_summed = np.asarray([0. for i in range(number_sn_bins-1)])

for i in range(len(cubes)):
    cube = SpectralCube.read("/home/jacob/Research/Wide_ASPECS/Data/gs_{}_2chn.fits".format(cubes[i]))

    print(cube)

    num_bins = 200

    # Now do it all!!!
    '''
    sub_cube = cube.unitless.unmasked_data[:,:,:]
    rms_noise = np.nanstd(sub_cube)

    #print(np.unique(sub_cube[~np.isnan(sub_cube)], return_counts=True))

    #print(sub_cube)

    sub_cube = sub_cube[~np.isnan(sub_cube)].flatten()

    print(sub_cube)

    _, bins, _ = plt.hist(sub_cube, bins=num_bins)
    plt.xlabel("Jy/beam")
    plt.yscale("log")
    plt.ylabel("Count (Log)")
    plt.savefig("Feb_Output/Fluxes_{}_2chn_All.png".format(cubes[i]))
    plt.cla()

    # Now flip the axis by splitting into < 0 and > 0 and subtracting one from other

    neg_sub_cube = sub_cube[sub_cube < 0]
    pos_sub_cube = sub_cube[sub_cube > 0]

    neg_sub_cube *= -1

    # Now digitize it and subtract sum of each bin together

    neg_bins = np.digitize(neg_sub_cube, bins)
    pos_bins = np.digitize(pos_sub_cube, bins)

    _, bins, _ = plt.hist(pos_bins, bins=num_bins)
    plt.hist(neg_bins, bins=bins)
    plt.title("Pos vs Neg Flux")
    plt.xlabel("S/N")
    plt.ylabel("Count (Log)")
    plt.yscale("log")
    plt.savefig("Feb_Output/Pos_Vs_Neg_{}_2chn_All.png".format(cubes[i]))
    plt.cla()

    _, bins, _ = plt.hist(pos_sub_cube, bins=num_bins)
    plt.hist(neg_sub_cube, bins=bins)
    plt.title("Pos vs Neg Flux")
    plt.xlabel("Flux (Jy/beam)")
    plt.ylabel("Count (Log)")
    plt.yscale("log")
    plt.savefig("Feb_Output/Pos_Vs_Neg_Values_{}_2chn_All.png".format(cubes[i]))
    plt.cla()

    # Now calculate S/N by making a few bins and dividing it out

    num_bins = 10

    sub_cube /= rms_noise
    print(np.min(sub_cube))
    print(np.max(sub_cube))
    print(len(sub_cube[sub_cube < -5.5]))
    print(len(sub_cube[sub_cube > 5.5]))

    # Divide the false by the true to get the fidelity, do the on its own


    # Bins
    bins = np.linspace(-6.5, 6.5, number_sn_bins)
    #bins = [-6.5,-5.5,-4.5,-3.5-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5]

    values, bins, _ = plt.hist(sub_cube, bins=bins)
    plt.title("S/N")
    plt.xlabel("S/N")
    plt.ylabel("Count (Log)")
    plt.yscale("log")
    plt.savefig("Feb_Output/SN_{}_2chn_All.png".format(cubes[i]))
    plt.cla()

    # Now do the Fidelity
    # First by working backwards
    fid = []
    for j in range(int(len(values)/2.)):
        fidel = values[j] / values[int(len(values))-j-1]
        fid.append((bins[j], fidel))

    # So goes in reverse order, False / Pos from -6.5 to 0
    fid = list(reversed(fid))
    print("Fidelity Fractions")
    print(fid)

    # Divide positive by false
    neg_bins = values[0:int(number_sn_bins/2.)-1]
    neg_bins = list(reversed(neg_bins))
    pos_bins = values[int(number_sn_bins/2.):]

    print(len(neg_bins))
    print(len(pos_bins))
    print(neg_bins)
    print(pos_bins)
    print(bins)

    fidelity = np.asarray(neg_bins) / pos_bins

    fidelity = list(fidelity)

    fidelity = fidelity[0:int(len(fidelity)/2.)] + [1] + fidelity[int(len(fidelity)/2.):]

    # Get from

    plt.hist(fidelity, bins=number_sn_bins)
    plt.title("Fidelity (False/Positive)")
    plt.xlabel("Fraction of False/True")
    plt.ylabel("Count (Log)")
    plt.yscale('log')
    plt.savefig("Feb_Output/Fidelity_{}_2chn_All.png".format(cubes[i]))
    print(fidelity)
    # TODO Fidelity
    '''
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

        #print(np.unique(sub_cube[~np.isnan(sub_cube)], return_counts=True))

        #print(sub_cube)

        sub_cube = sub_cube[~np.isnan(sub_cube)].flatten()

        _, bins, _ = plt.hist(sub_cube, bins=num_bins)
        plt.xlabel("Jy/beam")
        plt.yscale("log")
        plt.ylabel("Count (Log)")
        plt.savefig("Feb_Output/Fluxes/Fluxes_{}_2chn_slice{}.png".format(cubes[i],slice))
        plt.cla()

        # Now flip the axis by splitting into < 0 and > 0 and subtracting one from other

        neg_sub_cube = sub_cube[sub_cube < 0]
        pos_sub_cube = sub_cube[sub_cube > 0]

        neg_sub_cube *= -1

        sub_cube /= rms_noise
        if np.min(sub_cube) < sn_min:
            sn_min = np.min(sub_cube)
        if np.max(sub_cube) > sn_max:
            sn_max = np.max(sub_cube)
        num_bins = 10
        bins = np.linspace(-6.5, 6.5, number_sn_bins)

        values, bins, _ = plt.hist(sub_cube, bins=bins)
        sn_summed = sn_summed + np.asarray(values)
        plt.title("S/N")
        plt.xlabel("S/N")
        plt.ylabel("Count (Log)")
        plt.yscale("log")
        plt.savefig("Feb_Output/Noise/SN_{}_2chn_slice{}.png".format(cubes[i],slice))
        plt.cla()

        fid = []
        for j in range(int(len(values)/2.)):
            if np.isclose(values[int(len(values))-j-1], 0.0):
                fidel = np.nan
            else:
                fidel = values[j] / values[int(len(values))-j-1]
            fid.append(fidel)

        print(fid[0])
        if fid[2] is not np.nan:
            avg_fidelity_above_53 += fid[2]
            count_53 += 1
        if fid[1] is not np.nan:
            avg_fidelity_above_59 += fid[1]
            count_59 += 1
        if fid[0] is not np.nan:
            avg_fidelity_above_65 += fid[0]
            count_65 += 1

    print(sn_max)
    print(sn_min)
    print(count_65)
    #print(avg_fidelity_above_65 / count_65)
    print(count_59)
    #print(avg_fidelity_above_59 / count_59)
    print(count_53)
    #print(avg_fidelity_above_53 / count_53)
    plt.hist(sn_summed)
    plt.title("S/N Summed Cube")
    plt.yscale("log")
    plt.show()
    plt.cla()
    del sub_cube

