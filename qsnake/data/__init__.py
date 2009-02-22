_numbers = None
def symbol2number(symbol):
    """
    Converts an element symbol to an atomic number.

    Example:
    >>> symbol2number("H")
    1
    >>> symbol2number("He")
    2
    >>> symbol2number("Br")
    35
    """
    global _numbers
    if _numbers is None:
        _numbers = {}
        # initialize _symbols:
        number2symbol(1)
        for n, name in enumerate(_symbols):
            _numbers[name] = n+1
    return _numbers[symbol]

_symbols = None
def number2symbol(Z):
    """
    Converts an atomic number to an element name.

    Example:
    >>> number2symbol(1)
    'H'
    >>> number2symbol(2)
    'He'
    >>> number2symbol(35)
    'Br'
    """
    global _symbols
    if _symbols is None:
        _symbols = ["H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F", "Ne",
                "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K", "Ca",
                "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
                "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
                "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
                "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
                "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
                "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
                "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
                "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
                "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
                "Rg", "Uub", "Uut", "Uuq", "Uup", "Uuh", "Uus", "UUo"]
    return _symbols[Z-1]

_names = None
def number2name(Z):
    """
    Converts an atomic number to an element name.

    Example:
    >>> number2name(1)
    'Hydrogen'
    >>> number2name(2)
    'Helium'
    >>> number2name(35)
    'Bromine'
    """
    global _names
    if _names is None:
        _names = ["Hydrogen", "Helium", "Lithium", "Beryllium", "Boron",
                "Carbon", "Nitrogen", "Oxygen", "Fluorine", "Neon", "Sodium",
                "Magnesium", "Aluminium", "Silicon", "Phosphorus", "Sulfur",
                "Chlorine", "Argon", "Potassium", "Calcium", "Scandium",
                "Titanium", "Vanadium", "Chromium", "Manganese", "Iron",
                "Cobalt", "Nickel", "Copper", "Zinc", "Gallium", "Germanium",
                "Arsenic", "Selenium", "Bromine", "Krypton", "Rubidium",
                "Strontium", "Yttrium", "Zirconium", "Niobium", "Molybdenum",
                "Technetium", "Ruthenium", "Rhodium", "Palladium", "Silver",
                "Cadmium", "Indium", "Tin", "Antimony", "Tellurium", "Iodine",
                "Xenon", "Caesium", "Barium", "Lanthanum", "Cerium",
                "Praseodymium", "Neodymium", "Promethium", "Samarium",
                "Europium", "Gadolinium", "Terbium", "Dysprosium", "Holmium",
                "Erbium", "Thulium", "Ytterbium", "Lutetium", "Hafnium",
                "Tantalum", "Tungsten", "Rhenium", "Osmium", "Iridium",
                "Platinum", "Gold", "Mercury", "Thallium", "Lead", "Bismuth",
                "Polonium", "Astatine", "Radon", "Francium", "Radium",
                "Actinium", "Thorium", "Protactinium", "Uranium", "Neptunium",
                "Plutonium", "Americium", "Curium", "Berkelium", "Californium",
                "Einsteinium", "Fermium", "Mendelevium", "Nobelium",
                "Lawrencium", "Unnilquadium", "Unnilpentium", "Unnilhexium"]
    return _names[Z-1]
