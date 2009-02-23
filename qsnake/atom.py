from numpy import array
from data import symbol2number, number2symbol, number2name

class Atom(object):

    def __init__(self, symbol_or_Z=None, position=(0, 0, 0)):
        try:
            self._Z = int(symbol_or_Z)
        except ValueError:
            self._Z = symbol2number(str(symbol_or_Z))
        self._position = array(position)

    @property
    def number(self):
        """
        Returns the atomic number.

        Example:
        >>> Atom("B").number
        5
        """
        return self._Z

    @property
    def symbol(self):
        """
        Returns the symbol of the atom.

        Example:
        >>> Atom("B").symbol
        'B'
        """
        return number2symbol(self._Z)

    @property
    def name(self):
        """
        Returns the atom's name.

        Example:
        >>> Atom("B").name
        'Boron'
        """
        return number2name(self._Z)

    def __str__(self):
        return "<%s: %s>" % (number2symbol(self._Z), self._position)

    @property
    def position(self):
        return self._position
