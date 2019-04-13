import numpy as np
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube
from scipy.optimize import curve_fit
import pickle

cubes = ["A1", "A2", "B1", "B2", "C1", "C2"]

number_sn_bins = 260

def gaussian(x, mu, sig, height):
    return (np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))))

all_sn_summed = np.asarray([0. for i in range(number_sn_bins-1)])
all_fids = []
for i in range(len(cubes)):
    total_rms = 0.0
    cube = SpectralCube.read("/media/jacob/WDRed8Tb1/gs_{}_2chn.fits".format(cubes[i]))
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
        print(rms_noise)

        sub_cube = sub_cube[~np.isnan(sub_cube)].flatten()
        print((len(sub_cube)*0.09)/3600)
        # Fideltiy is 1 - Neg(SN)/Pos(SN)
        # Walter 2016

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
        bins = np.linspace(-12, 12, number_sn_bins)

        values, bins, _ = plt.hist(sub_cube, bins=bins)
        plt.cla()
        sn_summed = sn_summed + np.asarray(values)
        all_sn_summed = all_sn_summed + np.asarray(values)
        sn_bins = bins

        neg_summed = np.asarray(list(reversed(values[0:int(len(values)/2)])))
        pos_summed = np.asarray(values[int(len(values)/2)+1:])
        fid = 1. - neg_summed / pos_summed
        all_fids.append(fid)

    print(sn_max)
    print(sn_min)
    print(count_65)
    #print(avg_fidelity_above_65 / count_65)
    print(count_59)
    #print(avg_fidelity_above_59 / count_59)
    print(count_53)
    #print(avg_fidelity_above_53 / count_53)
    sn_summed = sn_summed
    print(sn_summed)
    n = len(sn_summed)                          #the number of data
    mean = np.mean(sn_summed)                   #note this correction
    sigma = 1.        #note this correction

    def gaus(x,a,x0,sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))

    #popt,pcov = curve_fit(gaus,np.linspace(-6.5, 6.5, 129),sn_summed,p0=[1,mean,sigma])
    plt.hist(np.linspace(-12, 12, number_sn_bins-1), weights=sn_summed, bins=sn_bins, histtype='step')
    #plt.plot(np.linspace(-6.5, 6.5, 129),gaus(np.linspace(-6.5, 6.5, 129),*popt),c='r',label='fit')
    plt.title("S/N Summed Cube {}".format(cubes[i]))
    plt.yscale("log")
    plt.xlabel("S/N")
    plt.ylabel("Count")
    plt.savefig("SN_Summed_Cube_{}_stepped.png".format(cubes[i]), dpi=300)
    plt.show()
    plt.cla()
    fid = []
    neg_summed = np.asarray(list(reversed(sn_summed[0:int(len(sn_summed)/2)])))
    pos_summed = np.asarray(sn_summed[int(len(sn_summed)/2)+1:])
    fid = 1. - neg_summed / pos_summed
    print(fid)
    print(sn_summed)
    print(fid)
    plt.hist(np.linspace(0., 12., len(fid)), weights=fid, bins=sn_bins[int(len(sn_bins)/2.):], histtype='step')
    plt.title("Fidelity Summed Cube {}".format(cubes[i]))
    plt.xlabel("SN")
    plt.ylabel("1- Neg(S/N)/Pos(S/N)")
    plt.savefig("Fidelity_Summed_Cube_Hist_{}_stepped.png".format(cubes[i]), dpi=300)
    plt.show()
    plt.cla()
    with open("Fidelity_Summed_Cube_{}.p", "wb") as pickle_fle:
        pickle.dump(fid, pickle_fle)
    with open("All_Fidelity_Summed_Cube_{}.p", "wb") as pickle_fle:
        pickle.dump(all_fids, pickle_fle)
    with open("SN_Summed_Cube_{}.p", "wb") as pickle_fle:
        pickle.dump(sn_summed, pickle_fle)
    with open("All_SN_Summed_Cube_{}.p", "wb") as pickle_fle:
        pickle.dump(all_sn_summed, pickle_fle)

    del cube

n = len(all_sn_summed)                          #the number of data
mean = np.mean(all_sn_summed)                   #note this correction
sigma = 1.        #note this correction

def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

#popt,pcov = curve_fit(gaus,np.linspace(-6.5, 6.5, 129),all_sn_summed,p0=[1,mean,sigma])
plt.hist(np.linspace(-12, 12, number_sn_bins-1), weights=all_sn_summed, bins=sn_bins, histtype='step')
#plt.plot(np.linspace(-6.5, 6.5, 129), gaus(np.linspace(-6.5, 6.5, 129),*popt), c='r')
plt.title("S/N Summed Cube Both Cubes")
plt.yscale("log")
plt.xlabel("S/N")
plt.ylabel("Count")
plt.savefig("SN_Summed_Cube_Both_Cubes_stepped.png", dpi=300)
plt.show()
plt.cla()
fid = []
neg_summed = np.asarray(list(reversed(all_sn_summed[0:int(len(all_sn_summed)/2)])))
pos_summed = np.asarray(all_sn_summed[int(len(all_sn_summed)/2)+1:])
fid = 1. - neg_summed / pos_summed
print(fid)
print(all_sn_summed)
print(fid)
plt.hist(np.linspace(0., 12., len(fid)), weights=fid, bins=sn_bins[int(len(sn_bins)/2.):], histtype='step')
plt.title("Fidelity Summed Cube Both")
plt.xlabel("SN")
plt.ylabel("1- Neg(S/N)/Pos(S/N)")
plt.savefig("Fidelity_Summed_Cube_Hist_All_stepped.png", dpi=300)
plt.show()
plt.cla()
plt.plot(np.linspace(0., 12., len(fid)), fid)#, bins=sn_bins[int(len(sn_bins)/2.):])
plt.title("Fidelity Summed Cube Both")
plt.xlabel("SN")
plt.ylabel("1- Neg(S/N)/Pos(S/N)")
plt.savefig("Fidelity_Summed_Cube_All_stepped.png", dpi=300)
plt.show()
plt.cla()

# Now go through the full fidelity
all_fids = np.asarray(all_fids)
print(all_fids.shape)
mean_fid = np.nanmean(all_fids, axis=0, dtype=np.float128)
median_fid = np.nanmedian(all_fids, axis=0)

print("Median")
print(median_fid)
print("Mean")
print(mean_fid)

plt.plot(np.linspace(0., 12., len(mean_fid)), mean_fid)#, bins=sn_bins[int(len(sn_bins)/2.):])
plt.title("Fidelity Mean Summed Cube")
plt.xlabel("SN")
plt.ylabel("1- Neg(S/N)/Pos(S/N)")
plt.savefig("Fidelity_Summed_Cube_Mean_All_stepped.png", dpi=300)
plt.show()
plt.cla()

plt.plot(np.linspace(0., 12., len(median_fid)), median_fid)#, bins=sn_bins[int(len(sn_bins)/2.):])
plt.title("Fidelity Median Summed Cube")
plt.xlabel("SN")
plt.ylabel("1- Neg(S/N)/Pos(S/N)")
plt.savefig("Fidelity_Summed_Cube_Median_All_stepped.png", dpi=300)
plt.show()
plt.cla()