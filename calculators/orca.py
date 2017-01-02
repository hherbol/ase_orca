import os

from ase.calculators.calculator import FileIOCalculator, Parameters, ReadError


def which(program):
    '''
    A function to return the full path of a system executable.

    **Parameters**

        program: *str*
            The name of the system executable to find.

    **Returns**

        path: *str or None*
            The path to the system executable. If none exists, then None.

    **Refs**

        * http://stackoverflow.com/a/377028
    '''

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


class orca(FileIOCalculator):
    """
    orca calculator
    """
    name = which("orca")
    if name is None:
        raise Exception("Unable to find orca installation.")

    implemented_properties = ['energy', 'forces']
    command = name + ' PREFIX.inp > PREFIX.log'

    default_parameters = {'charge': 0,
                          'multiplicity': 1,
                          'method': 'HF-3c',
                          'basis': None,
                          'RunTyp': None,
                          'COSMO': None}

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='orca', atoms=None, scratch=None,
                 basisfile=None, extra=None, **kwargs):

        """Constructs a orca-calculator object.

        extra: any extra text to be included in the input card

        """

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        if restart is not None:
            try:
                self.read(restart)
            except ReadError:
                if ignore_bad_restart_file:
                    self.reset()
                else:
                    raise

        self.scratch = scratch
        self.basisfile = basisfile

        # store extra parameters
        self.extra = extra

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()
        return changed_parameters

    # Not sure what this function is for yet
    def check_state(self, atoms):
        system_changes = FileIOCalculator.check_state(self, atoms)

        ignore = ['cell', 'pbc']
        for change in system_changes:
            if change in ignore:
                system_changes.remove(change)

        return system_changes

    def write_input(self, atoms, properties=None, system_changes=None):
        """Writes the input file"""
        FileIOCalculator.write_input(self, atoms, properties, system_changes)

        self.parameters.write(self.label + '.ase')

        filename = self.label + '.inp'
        inputfile = open(filename, 'w')

        if self.extra is None:
            self.extra = ""
        route = '!'
        route += " %s" % self.parameters["method"]
        if self.parameters["basis"] is not None:
            route += " %s" % self.parameters["basis"]
        if self.parameters["COSMO"] is not None:
            route += " COSMO"
            if "cosmo" not in self.extra.lower():
                self.extra += "\n%cosmo\n SMD true\n solvent \""
                self.extra += self.parameters["COSMO"] + "\"\n end"

        if self.parameters["RunTyp"] is not None:
            if self.parameters["RunTyp"].upper() == "SP":
                self.parameters["RunTyp"] = "Energy"
            elif self.parameters["RunTyp"].upper() == "OPT":
                self.parameters["RunTyp"] = "Opt"
            elif self.parameters["RunTyp"].upper() == "ENGRAD":
                self.parameters["RunTyp"] = "Gradient"
            s = "\n%method\n RunTyp " + self.parameters["RunTyp"] + "\n end"
            if s not in self.extra:
                self.extra += s

        inputfile.write(route)
        inputfile.write(self.extra)

        charge = self.parameters["charge"]
        multiplicity = self.parameters["multiplicity"]
        inputfile.write("\n*xyz %d %d\n" % (charge, multiplicity))

        # Input coordinates
        symbols = atoms.get_chemical_symbols()
        coordinates = atoms.get_positions()
        for i in range(len(atoms)):
            inputfile.write('%-10s' % symbols[i])
            for j in range(3):
                inputfile.write('%20.10f' % coordinates[i, j])
            inputfile.write('\n')
        inputfile.write("*")

        inputfile.close()

    def read(self, label):
        """Used to read the results of a previous calculation if restarting"""
        FileIOCalculator.read(self, label)

        from ase.io.orca import read_orca_out
        filename = self.label

        if not os.path.isfile(filename + ".log"):
            raise ReadError

        self.atoms = read_orca_out(filename, quantity='atoms')
        self.parameters = Parameters.read(self.label + '.ase')
        self.read_results()

    def read_results(self):
        """Reads the output file using orcaReader"""
        from ase.io.orca import read_orca_out

        filename = self.label + ".log"

        if not os.path.isfile(filename):
            raise ReadError

        quantities = ['energy', 'forces']
        with open(filename, 'r') as fileobj:
            for quant in quantities:
                self.results[quant] = read_orca_out(fileobj,
                                                    quantity=quant)

    def clean(self):
        """Cleans up from a previous run"""
        extensions = ['.chk', '.log']

        for ext in extensions:
            f = self.label + ext
            try:
                if (self.directory is not None):
                    os.remove(os.path.join(self.directory, f))
                else:
                    os.remove(f)
            except OSError:
                pass

    def get_version(self):
        return self.read_output(self.label + '.log', 'version')
