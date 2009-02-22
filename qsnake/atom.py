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
        return self._Z

    @property
    def symbol(self):
        return number2symbol(self._Z)

    @property
    def name(self):
        return number2name(self._Z)

    def __str__(self):
        return "<%d: %s>" % (number2symbol(self._Z), self._position)
