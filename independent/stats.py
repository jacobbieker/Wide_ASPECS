"""

This is focused on providing the statistics used in the other analysis in one place

"""
import numpy as np
from spectral_cube import SpectralCube
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from astropy import constants as const
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

wdth = 7.813302339355E+06 *u.Hz
print(wdth)

transitions = {"1-0": [0.0030, 0.3694, 115.271, 0.2801, 89],
               "2-1": [1.0059, 1.7387, 230.538, 1.4277, 1920],
               "3-2": [2.0088, 3.1080, 345.796, 2.6129, 3363],
               "4-3": [3.0115, 4.4771, 461.041, 3.8030, 4149], }

transitions1 = {"1-0": [0.0030, 0.3694, 115.271, 0.2801, 89],
                   "2-1": [1.0059, 1.7387, 230.538, 1.4277, 1920],
                   "3-2": [2.0088, 3.1080, 345.796, 2.6129, 3363],
                   "4-3": [3.0115, 4.4771, 461.041, 3.8030, 4149],
                   "5-4": [4.0142, 5.8460, 576.268, 4.9933, 4571],
                   "6-5": [5.0166, 7.2146, 691.473, 6.1843, 4809],
                   "7-6": [6.0188, 8.5829, 806.652, 7.3750, 4935],}
transitions1 = {"2-1": [0.0, 0.0873, 230.538, 0.0656, 1.4],
               "3-2": [0.2713, 0.6309, 345.796, 0.4858, 314],
               "4-3": [0.6950, 1.1744, 461.041, 0.9543, 1028],
               "5-4": [1.1186, 1.7178, 576.268, 1.4297, 1759],
               "6-5": [1.5422, 2.2612, 691.473, 1.9078, 2376],
               "7-6": [1.9656, 2.8044, 806.652, 2.3859, 2864],}

temp = {
    "C1mm": [0.8094, 1.3212, 492.161, 1.0828, 1233],
    "C1_2-1mm": [1.9755, 2.8171, 809.342, 2.3973, 2875],
    "C2": [5.9873, 7.9635, 1900.548, 6.9408, 4431],
    "C1": [3.2823, 4.8468, 492.161, 4.1242, 4287],
    "C1_2-1": [6.0422, 8.6148, 809.342, 7.4031, 4936],
}
def comoving_volume(start_z, end_z, sq_arcminutes):
    start_vol = cosmo.comoving_volume(start_z)
    end_vol = cosmo.comoving_volume(end_z)
    fraction_sky = sq_arcminutes / 148510800.
    comoving_volume = (end_vol - start_vol) * fraction_sky
    return comoving_volume

def fidelity(spectral_cube):
    return NotImplementedError

def signal_to_noise(spectral_cube):
    return NotImplementedError

def fluxes(spectral_cube):
    return NotImplementedError

def get_kms(channels, observed_ghz):
    channel_width = 7.813302339355E+06 *u.Hz
    #channel_width = 3.9025 * u.MHz
    channel_width = channel_width.to(u.GHz)
    observed_ghz = observed_ghz * u.GHz
    width = channels * channel_width

    #299792.458 = speed of light in m/s
    fwhm = (width * const.c.to('km/s')) / observed_ghz

    #fwhm = fwhm * (u.km / u.s)

    return fwhm

print(get_kms(1, 96))
print(get_kms(3, 96))
print(get_kms(19, 96))

def get_co_z(observed_ghz, transition):
    observed_ghz = observed_ghz * u.GHz
    emitted_ghz = transitions[transition][2] * u.GHz

    z_co = np.round((emitted_ghz)/ observed_ghz - 1, 3)
    return z_co

def convert_deltaZ_to_kms(delta_z, z):
    # C*delta_Z / (1+z) = Delta V only valid when delta Z is small
    delta_v = (delta_z * const.c.to('km/s')) / (1+z)
    return delta_v


def has_spec_z(source):
    spec_z_mask = (source["z_spec_3dh"] > 0.0001) | (source["zm_vds"] > 0.0001) | (
            source["zm_coeS"] > 0.0001) \
                  | (source["zs_mor"] > 0.0001) | (source["zm_ina"] > 0.0001) | (source["zm_her"] > 0.0001)
                  #| (source['muse_wide_z'] > 0.0001)
    return spec_z_mask

