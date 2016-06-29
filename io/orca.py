import numpy as np

import ase.units
from ase.data import chemical_symbols
from ase.atoms import Atoms
from ase.atom import Atom
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io.orca_reader import OrcaReader
from ase.calculators.orca import orca

def read_orca_out(filename, index=-1, quantity='atoms'):
    """"Interface to orcaReader and returns various quantities"""
    energy = 0.0

    data = OrcaReader(filename)[index]
    if isinstance(data, list):
        msg = 'Cannot parse multiple images from orca out files at this'
        msg += ' time.  Please select a single image.'
        raise RuntimeError(msg)

    atomic_numbers = data['Atomic_numbers']
    formula = str()
    for number in atomic_numbers:
        formula += chemical_symbols[number]

    positions = np.array(data['Positions'])
    method = data['Method']
    version = data['Version']
    charge = data['Charge']
    multiplicity = data['Multiplicity']
    energy = data['Energy'] * ase.units.Hartree
    if data['Gradient'] is None:
        forces = None
    else:
        forces = [-1*np.array(d) for d in data['Gradient']]
        convert = ase.units.Hartree / ase.units.Bohr
        forces = np.array(forces) * convert

    atoms = Atoms(formula, positions=positions)

    calc = SinglePointCalculator(atoms, energy=energy, forces=forces)
    atoms.set_calculator(calc)

    if (quantity == 'energy'):
        return energy
    elif (quantity == 'forces'):
        return forces
    elif (quantity == 'dipole'):
        return np.array(data['Dipole'])
    elif (quantity == 'atoms'):
        return atoms
    elif (quantity == 'version'):
        return version
    elif (quantity == 'multiplicity'):
        return multiplicity
    elif (quantity == 'charge'):
        return charge


def read_orca(filename):
    """Reads an orca input file"""
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    atoms = Atoms()
    for n, line in enumerate(lines):
        if ('*xyz' in line):
            i = 0
            while (lines[n + i + 1] != '\n'):
                info = lines[n + i + 1].split()
                symbol = info[0]
                position = [float(info[1]), float(info[2]), float(info[3])]
                atoms += Atom(symbol, position=position)
                i += 1
    return atoms


def write_orca(filename, atoms):
    """Writes a basic orca input file"""
# Since orca prints the geometry directly into the input file, we'll just
# the write_input method from the orca calculator, and just use the
# default settings
    calc = orca()
    calc.initialize(atoms)
    calc.write_input(filename, atoms)
