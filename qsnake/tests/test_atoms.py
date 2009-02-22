from qsnake import Atom, Atoms

def test_basic1():
    a1 = Atom("H", (0, 0, 0))
    a2 = Atom("H", (0, 0, 1))
    a3 = Atom("O", (1, 0, 1))
    atoms = Atoms([a1, a2, a3])
    assert len(atoms) == 3
    assert atoms[0] == a1
    assert atoms[1] == a2
    assert atoms[2] == a3
    assert atoms.get_atomic_numbers() == [1, 1, 8]
