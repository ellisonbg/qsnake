from qsnake import Atom

def test_basic1():
    a = Atom(1)
    assert a.number == 1
    assert a.symbol == "H"
    assert a.name == "Hydrogen"
    a = Atom("H")
    assert a.number == 1
    assert a.symbol == "H"
    assert a.name == "Hydrogen"

def test_basic2():
    a = Atom(5)
    assert a.number == 5
    assert a.symbol == "B"
    assert a.name == "Boron"
    a = Atom("B")
    assert a.number == 5
    assert a.symbol == "B"
    assert a.name == "Boron"

def test_basic3():
    a = Atom(82)
    assert a.number == 82
    assert a.symbol == "Pb"
    assert a.name == "Lead"
    a = Atom("Pb")
    assert a.number == 82
    assert a.symbol == "Pb"
    assert a.name == "Lead"
