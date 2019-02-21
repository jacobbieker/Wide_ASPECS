import numpy as np
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube

cube = SpectralCube.read("/home/jacob/Research/gs_A1_2chn.fits")

print(cube)

sub_cube = cube[0:]

print(cube.header)