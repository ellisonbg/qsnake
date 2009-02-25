"""
This example calculates three H atoms, reads the charge density from abinit,
converts it to real space by FFT and plots it using matplotlib.
"""

from math import sqrt

from numpy import empty
from scipy.fftpack import fftn

from qsnake import Atom, Atoms
from qsnake.calculators import Abinit

atoms = Atoms([
    Atom("H", (0, 0, -1)),
    Atom("H", (0, 0, 1)),
    Atom("H", (0, 0.5, 0)),
    ])
abinit = Abinit(atoms)
result = abinit.calculate(verbose=True)
d = result["density"]

# convert the density to real space
d = abs(fftn(d))

# plot it:
max_value = d.max()
from pylab import imshow, grid, show, subplot, colorbar, title, clf
clf()
x_max = 30
i_max = int(sqrt(x_max))+1
j_max = int(sqrt(x_max))+1
for i in range(x_max):
    subplot(i_max, j_max, i+1)
    imshow(d[i, :, :], vmax=max_value, interpolation="nearest")
    #colorbar()
    title("x=%d" % i)
    grid(True)
show()
