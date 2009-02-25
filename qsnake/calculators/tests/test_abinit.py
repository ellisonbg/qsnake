from qsnake import Atom, Atoms
from qsnake.calculators import Abinit

from qsnake.calculators.abinit import abinit_present
if not abinit_present():
    disabled = True

def feq(a, b, eps = 1e-10):
    """
    Compares two floating point numbers to a specified precision.
    """
    return abs(a-b) < eps

def test_basic1():
    a1 = Atom("H", (0, 0, 0))
    a2 = Atom("H", (0, 0, 1))
    a3 = Atom("H", (1, 0, 1))
    atoms = Atoms([a1, a2, a3])
    abinit = Abinit(atoms)
    r = abinit.calculate()
    assert feq(r["header"]["etotal"], -1.26918021328, 1e-5)
    assert feq(r["header"]["fermie"], -0.091979753480, 1e-5)
