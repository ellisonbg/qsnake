class Atoms(object):

    def __init__(self, atom_list):
        self._atoms = atom_list

    def __getitem__(self, key):
        return self._atoms[key]

    def __len__(self):
        return len(self._atoms)

    def get_atomic_numbers(self):
        """
        Returns the list of atomic numbers.

        Example:
        >>> from qsnake import Atom
        >>> a1 = Atom("H", (0, 0, 0))
        >>> a2 = Atom("H", (0, 0, 1))
        >>> a3 = Atom("O", (1, 0, 1))
        >>> atoms = Atoms([a1, a2, a3])
        >>> atoms.get_atomic_numbers()
        [1, 1, 8]

        """
        return [a.number for a in self]

    def __str__(self):
        s = ""
        for a in self._atoms:
            s += str(a)
        return s
