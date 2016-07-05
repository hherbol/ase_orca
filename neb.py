# -*- coding: utf-8 -*-
import threading
from math import sqrt

import numpy as np

import ase.parallel as mpi
from ase.build import minimize_rotation_and_translation, orthogonal_procrustes
from ase.calculators.calculator import Calculator
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read
from ase.optimize import BFGS
from ase.geometry import find_mic
import ase.units as units

class NEB:
    def __init__(self, images, k=0.1, climb=False, parallel=False,
                 remove_rotation_and_translation=False, procrustes=False, world=None, clancelot=False):
        """Nudged elastic band.

        images: list of Atoms objects
            Images defining path from initial to final state.
        k: float or list of floats
            Spring constant(s) in eV/Ang.  One number or one for each spring.
        climb: bool
            Use a climbing image (default is no climbing image).
        parallel: bool
            Distribute images over processors.
        remove_rotation_and_translation: bool
            TRUE actives NEB-TR for removing translation and
            rotation during NEB. By default applied non-periodic
            systems
        """
        self.images = images
        self.climb = climb
        self.parallel = parallel
        self.natoms = len(images[0])
        self.nimages = len(images)
        self.emax = np.nan
        self.use_clancelot = clancelot
        self.procrustes = procrustes
        
        self.remove_rotation_and_translation = remove_rotation_and_translation
        
        if isinstance(k, (float, int)):
            k = [k] * (self.nimages - 1)
        self.k = list(k)

        if world is None:
            world = mpi.world
        self.world = world

        if parallel:
            assert world.size == 1 or world.size % (self.nimages - 2) == 0
        if self.procrustes:
            self.remove_rotation_and_translation = True

    def rotate(self, image, R):
        pos = image.get_positions()
        for i,atom in enumerate(pos):
            pos[i] = np.dot(atom,R)
        return pos

    def interpolate(self, method='linear', mic=False):
        if self.remove_rotation_and_translation:
            if self.procrustes:
                rotation_matrix, _ = orthogonal_procrustes(self.images[-1], self.images[0])
                self.images[-1].set_positions(self.images[-1].get_positions()*rotation_matrix)
            else:
                minimize_rotation_and_translation(self.images[0], self.images[-1])

        
        interpolate(self.images, mic)
                 
        if method == 'idpp':
            self.idpp_interpolate(traj=None, log=None, mic=mic)

    def idpp_interpolate(self, traj='idpp.traj', log='idpp.log', fmax=0.1,
                         optimizer=BFGS, mic=False):
        d1 = self.images[0].get_all_distances(mic=mic)
        d2 = self.images[-1].get_all_distances(mic=mic)
        d = (d2 - d1) / (self.nimages - 1)
        old = []
        for i, image in enumerate(self.images):
            old.append(image.calc)
            image.calc = IDPP(d1 + i * d, mic=mic)
        opt = optimizer(self, trajectory=traj, logfile=log)
        opt.run(fmax=fmax)
        for image, calc in zip(self.images, old):
            image.calc = calc

    def get_positions(self):
        positions = np.empty(((self.nimages - 2) * self.natoms, 3))
        n1 = 0
        for image in self.images[1:-1]:
            n2 = n1 + self.natoms
            positions[n1:n2] = image.get_positions()
            n1 = n2
        return positions

    def set_positions(self, positions):
        n1 = 0
        for image in self.images[1:-1]:
            n2 = n1 + self.natoms
            image.set_positions(positions[n1:n2])
            n1 = n2

            # Parallel NEB with Jacapo needs this:
            try:
                image.get_calculator().set_atoms(image)
            except AttributeError:
                pass

    def get_forces(self):
        """Evaluate and return the forces."""
        images = self.images
        forces = np.empty((self.nimages, self.natoms, 3))
        
        energies = np.empty(self.nimages)

        if self.remove_rotation_and_translation:
            # Remove translation and rotation between
            # images before computing forces:
            for i in range(1, self.nimages):
                if self.procrustes:
                    rotation_matrix, _ = orthogonal_procrustes(self.images[i], self.images[i-1])
                    self.images[i].set_positions(self.rotate(self.images[i], rotation_matrix))
                else:
                    minimize_rotation_and_translation(images[i - 1], images[i])

        if not self.parallel:
            # Do all images - one at a time:
            for i in range(self.nimages):
                energies[i] = images[i].get_potential_energy()
                forces[i] = images[i].get_forces()
            forces = forces[1:-1]
            if not self.use_clancelot:
                energies_0 = energies[0]
                energies = energies[1:-1]

        elif self.world.size == 1:
            def run(image, energies, forces):
                energies[:] = image.get_potential_energy()
                forces[:] = image.get_forces()

            threads = [threading.Thread(target=run,
                                        args=(images[i],
                                              energies[i:i+1],
                                              forces[i:i+1]))
                       for i in range(self.nimages)]
            for i,thread in enumerate(threads):
                thread.start()
            for thread in threads:
                thread.join()
            forces = forces[1:-1]
            if not self.use_clancelot:
                energies_0 = energies[0]
                energies = energies[1:-1]
        else:
            # Parallelize over images:
            i = self.world.rank * (self.nimages) // self.world.size + 1
            try:
                energies[i] = images[i].get_potential_energy()
                forces[i] = images[i].get_forces()
            except:
                # Make sure other images also fail:
                error = self.world.sum(1.0)
                raise
            else:
                error = self.world.sum(0.0)
                if error:
                    raise RuntimeError('Parallel NEB failed!')

            for i in range(self.nimages):
                root = (i) * self.world.size // (self.nimages)
                self.world.broadcast(energies[i], root)
                self.world.broadcast(forces[i], root)
            forces = forces[1:-1]
            if not self.use_clancelot:
                energies_0 = energies[0]
                energies = energies[1:-1]

        if self.use_clancelot:
            V = energies.copy()
            # Output max energy in units of kT at 300 K
            self.emax = (max(V)-V[0]) / (units.kB * 300.0)
            for i in range(1, self.nimages-1):
                a = images[i-1].get_positions()
                b = images[i].get_positions()
                c = images[i+1].get_positions()
                real_force = forces[i - 1]
                k = self.k[i]

                # Find tangent
                tplus = c - b
                tminus = b - a
                dVmin = min( abs(V[i+1] - V[i]), abs(V[i-1]-V[i]) )
                dVmax = max( abs(V[i+1] - V[i]), abs(V[i-1]-V[i]) )
                if V[i+1] > V[i] and V[i] > V[i-1]:
                    tangent = tplus.copy()
                elif V[i+1] < V[i] and V[i] < V[i-1]:
                    tangent = tminus.copy()
                elif V[i+1] > V[i-1]:
                    tangent = tplus*dVmax + tminus*dVmin
                else:
                    tangent = tplus*dVmin + tminus*dVmax

                # Normalize tangent
                tangent_norm = np.sqrt(np.vdot(tangent,tangent))
                if tangent_norm != 0: tangent /= tangent_norm

                F_spring_parallel = k * (np.linalg.norm(tplus) - np.linalg.norm(tminus)) * tangent 

                F_real_parallel = (np.vdot(real_force,tangent)*tangent).reshape((-1,3))
                F_real_perpendicular = real_force - F_real_parallel

                # Set NEB forces
                forces[i-1] = F_spring_parallel + F_real_perpendicular
        else:
            imax = 1 + np.argsort(energies)[-1]
            #self.emax = energies[imax - 1]
            self.emax = (max(energies)-energies_0) / (units.kB * 300.0)

            tangent1 = find_mic(images[1].get_positions() -
                                images[0].get_positions(),
                                images[0].get_cell(), images[0].pbc)[0]
            for i in range(1, self.nimages - 1):
                tangent2 = find_mic(images[i + 1].get_positions() -
                                    images[i].get_positions(),
                                    images[i].get_cell(),
                                    images[i].pbc)[0]

                if i < imax: # tplus
                    tangent = tangent2
                elif i > imax: # tminus
                    tangent = tangent1
                else:
                    tangent = tangent1 + tangent2

                tt = np.vdot(tangent, tangent)

                f = forces[i - 1]
                real_force = forces[i-1].copy()

                ft = np.vdot(f, tangent)
                if i == imax and self.climb:
                    f -= 2 * ft / tt * tangent
                else:
                    F_real_parallel = ft / tt * tangent
                    f -= F_real_parallel
                    F_spring_parallel = np.vdot(tangent2 * self.k[i - 1] -
                                 tangent1 * self.k[i], tangent) / tt * tangent
                    f += F_spring_parallel

                tangent1 = tangent2

        return forces.reshape((-1, 3))

    def get_potential_energy(self):
        return self.emax

    def __len__(self):
        return (self.nimages - 2) * self.natoms


class IDPP(Calculator):
    """Image dependent pair potential.

    See:

        Improved initial guess for minimum energy path calculations.

        Søren Smidstrup, Andreas Pedersen, Kurt Stokbro and Hannes Jónsson

        Chem. Phys. 140, 214106 (2014)
    """

    implemented_properties = ['energy', 'forces']

    def __init__(self, target, mic):
        Calculator.__init__(self)
        self.target = target
        self.mic = mic

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        P = atoms.get_positions()
        d = []
        D = []
        for p in P:
            Di = P - p
            if self.mic:
                Di, di = find_mic(Di, atoms.get_cell(), atoms.get_pbc())
            else:
                di = np.sqrt((Di**2).sum(1))
            d.append(di)
            D.append(Di)
        d = np.array(d)
        D = np.array(D)

        dd = d - self.target
        d.ravel()[::len(d) + 1] = 1  # avoid dividing by zero
        d4 = d**4
        e = 0.5 * (dd**2 / d4).sum()
        f = -2 * ((dd * (1 - 2 * dd / d) / d**5)[..., np.newaxis] * D).sum(0)
        self.results = {'energy': e, 'forces': f}


class SingleCalculatorNEB(NEB):
    def __init__(self, images, k=0.1, climb=False):
        if isinstance(images, str):
            # this is a filename
            images = read(images)

        NEB.__init__(self, images, k, climb, False)
        self.calculators = [None] * self.nimages
        self.energies_ok = False
        self.first = True

    def interpolate(self, initial=0, final=-1, mic=False):
        """Interpolate linearly between initial and final images."""
        if final < 0:
            final = self.nimages + final
        n = final - initial
        pos1 = self.images[initial].get_positions()
        pos2 = self.images[final].get_positions()
        dist = (pos2 - pos1)
        if mic:
            cell = self.images[initial].get_cell()
            assert((cell == self.images[final].get_cell()).all())
            pbc = self.images[initial].get_pbc()
            assert((pbc == self.images[final].get_pbc()).all())
            dist, D_len = find_mic(dist, cell, pbc)
        dist /= n
        for i in range(1, n):
            self.images[initial + i].set_positions(pos1 + i * dist)

    def refine(self, steps=1, begin=0, end=-1, mic=False):
        """Refine the NEB trajectory."""
        if end < 0:
            end = self.nimages + end
        j = begin
        n = end - begin
        for i in range(n):
            for k in range(steps):
                self.images.insert(j + 1, self.images[j].copy())
                self.calculators.insert(j + 1, None)
            self.k[j:j + 1] = [self.k[j] * (steps + 1)] * (steps + 1)
            self.nimages = len(self.images)
            self.interpolate(j, j + steps + 1, mic=mic)
            j += steps + 1

    def set_positions(self, positions):
        # new positions -> new forces
        if self.energies_ok:
            # restore calculators
            self.set_calculators(self.calculators[1:-1])
        NEB.set_positions(self, positions)

    def get_calculators(self):
        """Return the original calculators."""
        calculators = []
        for i, image in enumerate(self.images):
            if self.calculators[i] is None:
                calculators.append(image.get_calculator())
            else:
                calculators.append(self.calculators[i])
        return calculators

    def set_calculators(self, calculators):
        """Set new calculators to the images."""
        self.energies_ok = False
        self.first = True

        if not isinstance(calculators, list):
            calculators = [calculators] * self.nimages

        n = len(calculators)
        if n == self.nimages:
            for i in range(self.nimages):
                self.images[i].set_calculator(calculators[i])
        elif n == self.nimages - 2:
            for i in range(1, self.nimages - 1):
                self.images[i].set_calculator(calculators[i - 1])
        else:
            raise RuntimeError(
                'len(calculators)=%d does not fit to len(images)=%d'
                % (n, self.nimages))

    def get_energies_and_forces(self):
        """Evaluate energies and forces and hide the calculators"""
        if self.energies_ok:
            return

        self.emax = -1.e32

        def calculate_and_hide(i):
            image = self.images[i]
            calc = image.get_calculator()
            if self.calculators[i] is None:
                self.calculators[i] = calc
            if calc is not None:
                if not isinstance(calc, SinglePointCalculator):
                    self.images[i].set_calculator(
                        SinglePointCalculator(
                            image,
                            energy=image.get_potential_energy(),
                            forces=image.get_forces()))
                self.emax = min(self.emax, image.get_potential_energy())

        if self.first:
            calculate_and_hide(0)

        # Do all images - one at a time:
        for i in range(1, self.nimages - 1):
            calculate_and_hide(i)

        if self.first:
            calculate_and_hide(-1)
            self.first = False

        self.energies_ok = True

    def get_forces(self):
        self.get_energies_and_forces()
        return NEB.get_forces(self)

    def n(self):
        return self.nimages

    def write(self, filename):
        from ase.io.trajectory import Trajectory
        traj = Trajectory(filename, 'w', self)
        traj.write()
        traj.close()

    def __add__(self, other):
        for image in other:
            self.images.append(image)
        return self


def fit0(E, F, R, cell=None, pbc=None):
    """Constructs curve parameters from the NEB images."""
    E = np.array(E) - E[0]
    n = len(E)
    Efit = np.empty((n - 1) * 20 + 1)
    Sfit = np.empty((n - 1) * 20 + 1)

    s = [0]
    dR = np.zeros_like(R)
    for i in range(n):
        if i < n - 1:
            dR[i] = R[i + 1] - R[i]
            if cell is not None and pbc is not None:
                dR[i], _ = find_mic(dR[i], cell, pbc)
            s.append(s[i] + sqrt((dR[i]**2).sum()))
        else:
            dR[i] = R[i] - R[i - 1]
            if cell is not None and pbc is not None:
                dR[i], _ = find_mic(dR[i], cell, pbc)

    lines = []
    dEds0 = None
    for i in range(n):
        d = dR[i]
        if i == 0:
            ds = 0.5 * s[1]
        elif i == n - 1:
            ds = 0.5 * (s[-1] - s[-2])
        else:
            ds = 0.25 * (s[i + 1] - s[i - 1])

        d = d / sqrt((d**2).sum())
        dEds = -(F[i] * d).sum()
        x = np.linspace(s[i] - ds, s[i] + ds, 3)
        y = E[i] + dEds * (x - s[i])
        lines.append((x, y))

        if i > 0:
            s0 = s[i - 1]
            s1 = s[i]
            x = np.linspace(s0, s1, 20, endpoint=False)
            c = np.linalg.solve(np.array([(1, s0, s0**2, s0**3),
                                          (1, s1, s1**2, s1**3),
                                          (0, 1, 2 * s0, 3 * s0**2),
                                          (0, 1, 2 * s1, 3 * s1**2)]),
                                np.array([E[i - 1], E[i], dEds0, dEds]))
            y = c[0] + x * (c[1] + x * (c[2] + x * c[3]))
            Sfit[(i - 1) * 20:i * 20] = x
            Efit[(i - 1) * 20:i * 20] = y

        dEds0 = dEds

    Sfit[-1] = s[-1]
    Efit[-1] = E[-1]
    return s, E, Sfit, Efit, lines


class NEBtools:
    """Class to make many of the common tools for NEB analysis available to
    the user. Useful for scripting the output of many jobs. Initialize with
    list of images which make up a single band."""

    def __init__(self, images):
        self._images = images

    def get_barrier(self, fit=True, raw=False):
        """Returns the barrier estimate from the NEB, along with the
        Delta E of the elementary reaction. If fit=True, the barrier is
        estimated based on the interpolated fit to the images; if
        fit=False, the barrier is taken as the maximum-energy image
        without interpolation. Set raw=True to get the raw energy of the
        transition state instead of the forward barrier."""
        s, E, Sfit, Efit, lines = self.get_fit()
        dE = E[-1] - E[0]
        if fit:
            barrier = max(Efit)
        else:
            barrier = max(E)
        if raw:
            barrier += self._images[0].get_potential_energy()
        return barrier, dE

    def plot_band(self, ax=None):
        """Plots the NEB band on matplotlib axes object 'ax'. If ax=None
        returns a new figure object."""
        if not ax:
            from matplotlib import pyplot
            fig = pyplot.figure()
            ax = fig.add_subplot(111)
        else:
            fig = None
        s, E, Sfit, Efit, lines = self.get_fit()
        ax.plot(s, E, 'o')
        for x, y in lines:
            ax.plot(x, y, '-g')
        ax.plot(Sfit, Efit, 'k-')
        ax.set_xlabel('path [$\AA$]')
        ax.set_ylabel('energy [eV]')
        Ef = max(Efit) - E[0]
        Er = max(Efit) - E[-1]
        dE = E[-1] - E[0]
        ax.set_title('$E_\mathrm{f} \\approx$ %.3f eV; '
                     '$E_\mathrm{r} \\approx$ %.3f eV; '
                     '$\\Delta E$ = %.3f eV'
                     % (Ef, Er, dE))
        return fig

    def get_fmax(self, **kwargs):
        """Returns fmax, as used by optimizers with NEB."""
        neb = NEB(self._images, **kwargs)
        forces = neb.get_forces()
        return np.sqrt((forces**2).sum(axis=1).max())

    def get_fit(self):
        """Returns the parameters for fitting images to band."""
        images = self._images
        if not hasattr(images, 'repeat'):
            from ase.gui.images import Images
            images = Images(images)
        N = images.repeat.prod()
        natoms = images.natoms // N
        R = images.P[:, :natoms]
        E = images.E
        F = images.F[:, :natoms]
        s, E, Sfit, Efit, lines = fit0(E, F, R, images.A[0], images.pbc)
        return s, E, Sfit, Efit, lines


def interpolate(images, mic=False):
    """Given a list of images, linearly interpolate the positions of the
    interior images."""
    pos1 = images[0].get_positions()
    pos2 = images[-1].get_positions()
    d = pos2 - pos1
    if mic:
        d = find_mic(d, images[0].get_cell(), images[0].pbc)[0]
    d /= (len(images) - 1.0)
    for i in range(1, len(images) - 1):
        images[i].set_positions(pos1 + i * d)
        # Parallel NEB with Jacapo needs this:
        try:
            images[i].get_calculator().set_atoms(images[i])
        except AttributeError:
            pass
