"""

This is focused on providing the statistics used in the other analysis in one place

"""
import numpy as np
from spectral_cube import SpectralCube

from astropy.cosmology import FlatLambdaCDM
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

def get_kms():
    return NotImplementedError

