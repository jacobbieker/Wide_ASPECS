"""

This is focused on providing the statistics used in the other analysis in one place

"""
import numpy as np
from spectral_cube import SpectralCube

def fidelity(spectral_cube):
    return NotImplementedError

def signal_to_noise(spectral_cube):
    return NotImplementedError

def fluxes(spectral_cube):
    return NotImplementedError

def comoving_volume():
    return NotImplementedError

def get_kms():
    return NotImplementedError

