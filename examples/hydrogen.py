"""
This example calculates two H atoms, reads the charge density from abinit,
converts it to real space by FFT and plots it using matplotlib.
"""

from numpy import empty
from scipy import fft

from qsnake import Atom, Atoms
from qsnake.calculators import Abinit

a1 = Atom("H", (0, 0, -1))
a2 = Atom("H", (0, 0, 1))
atoms = Atoms([a1, a2])
abinit = Abinit(atoms)
result = abinit.calculate()
d = result["density"]

# convert the density to real space
d = abs(fft(d))

# plot it:
max_value = d.max()
from pylab import imshow, grid, show, subplot, colorbar, title
subplot(2, 2, 1)
imshow(d[0, :, :], vmax=max_value, interpolation="nearest")
colorbar()
title("x=0")
grid(True)
subplot(2, 2, 2)
colorbar()
imshow(d[1, :, :], vmax=max_value, interpolation="nearest")
title("x=1")
grid(True)
subplot(2, 2, 3)
colorbar()
imshow(d[2, :, :], vmax=max_value, interpolation="nearest")
title("x=2")
grid(True)
subplot(2, 2, 4)
imshow(d[3, :, :], vmax=max_value, interpolation="nearest")
title("x=3")
grid(True)
colorbar()
grid(True)
show()
