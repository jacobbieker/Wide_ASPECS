"""

This is focused on providing the statistics used in the other analysis in one place

"""
import numpy as np
from spectral_cube import SpectralCube
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from astropy import constants as const
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

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

def get_co_z(observed_ghz, transition):
    observed_ghz = observed_ghz * u.GHz
    emitted_ghz = transitions[transition][2] * u.GHz

    z_co = np.round((emitted_ghz)/ observed_ghz - 1, 3)
    print(z_co)
    return z_co

def convert_deltaZ_to_kms(delta_z):
    delta_v = delta_z * const.c.to('km/s')
    return delta_v


def has_spec_z(source):
    spec_z_mask = (source["z_spec_3dh"] > 0.0001) | (source["zm_vds"] > 0.0001) | (
            source["zm_coeS"] > 0.0001) \
                  | (source["zs_mor"] > 0.0001) | (source["zm_ina"] > 0.0001) | (source["zm_her"] > 0.0001) \
                  | (source['muse_wide_z'] > 0.0001)
    return spec_z_mask

