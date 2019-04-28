import numpy as np
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube

"""

Steps:

1. Get the Table with units of the SN, rx, ry, potentially SkyCoord, width, and SN

2. cube3 = cube.with_spectral_unit(u.GHz, velocity_convention='radio',
...                                 rest_value=200 * u.GHz)

3. Get the slice with the closest to rounded GHz value to the number

4. Take (width - 1)/2. points in either direction

5. For each point, get the flux

6. Try with adding SN of each slice, and then as SN of all the included channels to see if it gets
close



"""