class Atoms(object):

    def __init__(self, atom_list):
        self._atoms = atom_list

    def __getitem__(self, key):
        return self._atoms[key]

    def __len__(self):
        return len(self._atoms)

    def get_atomic_numbers(self):
        return [a.number for a in self]

    def __str__(self):
        s = ""
        for a in self._atoms:
            s += str(a)
