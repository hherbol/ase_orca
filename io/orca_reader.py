from __future__ import print_function

class OrcaReader:

    def auto_type(self, data):
        """ tries to determine type"""
        try:
            return float(data)
        except ValueError:
            pass

        try:
            ds = data.split(",")
            array = []

            for d in ds:
                array.append(float(d))

            return array
        except ValueError:
            pass

        return data

    def __init__(self, filename):
        """filename is optional; if not set, use parse to set the content"""
        if isinstance(filename, str):
            fileobj = open(filename+".log", 'r')
        else:
            fileobj = filename
            fileobj.seek(0)  # Re-wind fileobj
        content = fileobj.read()
        content = content.replace("\r\n", "\n")
        self.parse(content)

    def parse(self, content):
        from ase.data import atomic_numbers
        self.data = []
        seq_count = 0

        try:
            route = [line[5:] for line in content.split('\n') if line.startswith('|  1>')][0]
        except IndexError:
            raise IOError('Could not find route line in %s: job most likely crashed.' % input_path)

        version = content.split("Program Version")[1].split("\n")[0]
        method = route.split("!")[1].split()[0]
        charge = int(content.split("Total Charge")[1].split('\n')[0].split()[-1])
        multiplicity = int(content.split("Multiplicity")[1].split('\n')[0].split()[-1])

        section, s_position = content, "CARTESIAN COORDINATES (ANGSTROEM)"
        s_energy = "FINAL SINGLE POINT ENERGY"
        s_gradient = "CARTESIAN GRADIENT"

        while s_position in section:
            section = section[section.find(s_position)+len(s_position):]
            atom_block = section[:section.find('\n\n')].split('\n')[2:]
            atoms, positions = [], []
            for line in atom_block:
                a = line.split()
                atoms.append(atomic_numbers[a[0]])
                positions.append([float(a[1]),float(a[2]),float(a[3])])
            energy = float(section[section.find(s_energy):].split("\n")[0].split()[-1])

            
            if s_gradient in section:
                grad_block = section[section.find(s_gradient):].split("\n\n")[1].split("\n")
                gradient = []
                for line in grad_block:
                    a = line.split()
                    gradient.append([float(b) for b in a[3:]])
            else:
                gradient = None

            new_dict = {}
            self.data.append(new_dict)

            new_dict["Method"] = method
            new_dict["Sequence number"] = seq_count
            new_dict["Charge"] = charge
            new_dict["Multiplicity"] = multiplicity
            new_dict["Atomic_numbers"] = atoms
            new_dict["Positions"] = positions
            new_dict["Energy"] = energy
            new_dict["Gradient"] = gradient
            new_dict["Version"] = version

    def __iter__(self):
        """returns an iterator that iterates over all keywords"""
        return self.data.__iter__()

    def __len__(self):
        return len(self.data)

    def __getitem__(self, pos):
        return self.data[pos]
