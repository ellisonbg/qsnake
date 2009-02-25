from qsnake import Atom, Atoms
from qsnake.calculators import Abinit

from numpy import array

from qsnake.calculators.abinit import abinit_present
if not abinit_present():
    disabled = True

def feq(a, b, eps = 1e-10):
    """
    Compares two floating point numbers to a specified precision.
    """
    return abs(a-b) < eps

def feq_array(a, b, eps = 1e-10):
    """
    Compares two floating point numbers to a specified precision.
    """
    return (abs(a-b) < eps).all()

def test_basic1():
    a1 = Atom("H", (0, 0, 0))
    a2 = Atom("H", (0, 0, 1))
    a3 = Atom("H", (1, 0, 1))
    atoms = Atoms([a1, a2, a3])
    abinit = Abinit(atoms)
    r = abinit.calculate()
    assert feq(r["header"]["etotal"], -1.26918021328, 1e-5)
    assert feq(r["header"]["fermie"], -0.091979753480, 1e-5)
    assert feq_array(r["wf"]["occ"], array([2., 1., 0.]), 1e-10)
    assert not feq_array(r["wf"]["occ"], array([1., 1., 1.]), 1e-10)
    assert feq_array(r["wf"]["eigen"], array([-0.8806383, -0.09197975,
        -0.08961579]), 1e-5)
    assert not feq_array(r["wf"]["eigen"], array([-0.8906383, -0.09197975,
        -0.08961579]), 1e-5)
