import numpy as np
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube

cube = SpectralCube.read("/home/jacob/Research/Wide_ASPECS/Data/gs_A1_2chn.fits")

print(cube)

num_bins = 200

# Now do it all!!!

sub_cube = cube.unitless.unmasked_data[:,:,:]

#print(np.unique(sub_cube[~np.isnan(sub_cube)], return_counts=True))

#print(sub_cube)

kept_sub_cube = sub_cube[~np.isnan(sub_cube)].flatten()

print(kept_sub_cube)

_, bins, _ = plt.hist(kept_sub_cube, bins=num_bins)
plt.xlabel("Jy/beam")
plt.yscale("log")
plt.ylabel("Count (Log)")
plt.savefig("Feb_Output/Fluxes_A1_2chn_All.png")
plt.cla()

# Now flip the axis by splitting into < 0 and > 0 and subtracting one from other

neg_sub_cube = kept_sub_cube[kept_sub_cube < 0]
pos_sub_cube = kept_sub_cube[kept_sub_cube > 0]

neg_sub_cube *= -1


# Now digitize it and subtract sum of each bin together

neg_bins = np.digitize(neg_sub_cube, bins)
pos_bins = np.digitize(pos_sub_cube, bins)

plt.hist(pos_bins, bins=num_bins)
plt.hist(neg_bins, bins=num_bins)
plt.title("Pos vs Neg Flux")
plt.xlabel("S/N")
plt.ylabel("Count (Log)")
plt.yscale("log")
plt.savefig("Feb_Output/Pos_Vs_Neg_A1_2chn_All.png")

exit()
print(max(pos_bins))
print(max(neg_bins))

bin_sums = [0 for i in range(num_bins+1)]

for index, element in enumerate(neg_sub_cube):
    bin_sums[neg_bins[index]] -= element

for index, element in enumerate(pos_sub_cube):
    try:
        bin_sums[pos_bins[index]] += element
    except IndexError:
        print(index)

plt.hist(bin_sums, bins=bins)
plt.title("S/N")
plt.xlabel("S/N")
plt.ylabel("Count (Log)")
plt.yscale("log")
plt.savefig("Feb_Output/SN_A1_2chn_All.png")

# Now go through each slice and do it for all of them
for slice in range(480):
    sub_cube = cube.unitless.unmasked_data[slice,:,:]

    #print(np.unique(sub_cube[~np.isnan(sub_cube)], return_counts=True))

    #print(sub_cube)

    kept_sub_cube = sub_cube[~np.isnan(sub_cube)].flatten()

    print(kept_sub_cube)

    _, bins, _ = plt.hist(kept_sub_cube, bins=num_bins)
    plt.xlabel("Jy/beam")
    plt.yscale("log")
    plt.ylabel("Count (Log)")
    plt.savefig("Feb_Output/Fluxes/Fluxes_A2_2chn_slice{}.png".format(slice))
    plt.cla()

    # Now flip the axis by splitting into < 0 and > 0 and subtracting one from other

    neg_sub_cube = kept_sub_cube[kept_sub_cube < 0]
    pos_sub_cube = kept_sub_cube[kept_sub_cube > 0]

    neg_sub_cube *= -1


    # Now digitize it and subtract sum of each bin together

    neg_bins = np.digitize(neg_sub_cube, bins)
    pos_bins = np.digitize(pos_sub_cube, bins)

    print(max(pos_bins))
    print(max(neg_bins))

    bin_sums = [0 for i in range(num_bins+1)]

    for index, element in enumerate(neg_sub_cube):
        try:
            bin_sums[neg_bins[index]] += element
        except IndexError:
            print(index)
    for index, element in enumerate(pos_sub_cube):
        try:
            bin_sums[pos_bins[index]] += element
        except IndexError:
            print(index)

    plt.hist(bin_sums, bins=bins)
    plt.title("S/N")
    plt.xlabel("S/N")
    plt.ylabel("Count (Log)")
    plt.yscale("log")
    plt.savefig("Feb_Output/Noise/SN_A2_2chn_slice{}.png".format(slice))
    plt.cla()


