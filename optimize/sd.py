# My custom implementation of steepest descent
# -*- coding: utf-8 -*-
import numpy as np

from ase.optimize.optimize import Optimizer
import ase.units as units


class SD(Optimizer):
    def __init__(self, atoms, restart=None, logfile='-', trajectory=None,
                 maxstep=None, master=None):
        """SD optimizer.

        Parameters:

        atoms: Atoms object
            The Atoms object to relax.

        restart: string
            Pickle file used to store hessian matrix. If set, file with
            such a name will be searched and hessian matrix stored will
            be used, if the file exists.

        trajectory: string
            Pickle file used to store trajectory of atomic movement.

        logfile: file object or str
            If *logfile* is a string, a file with that name will be opened.
            Use '-' for stdout.

        maxstep: float
            Used to set the maximum distance an atom can move per
            iteration (default value is 0.04 Å).

        master: boolean
            Defaults to None, which causes only rank 0 to save files.  If
            set to true,  this rank will save files.
        """
        Optimizer.__init__(self, atoms, restart, logfile, trajectory, master)

        if maxstep is not None:
            if maxstep > 1.0:
                raise ValueError('You are using a much too large value for ' +
                                 'the maximum step size: %.1f Å' % maxstep)
            self.maxstep = maxstep
        self.alpha = self.maxstep

    def initialize(self):
        self.r0 = None
        self.f0 = None
        self.maxstep = 0.1
        self.alpha = 0.1

    def read(self):
        self.H, self.r0, self.f0, self.maxstep = self.load()

    def step(self, f_eV):
        atoms = self.atoms
        r = atoms.get_positions()
        # To maintain similarity in comparison with clancelot, convert force
        # to units of hartree
        f = f_eV.copy() / units.Hartree
        f = f.reshape(-1).reshape((-1, 3))

        max_step_length = np.sqrt(((f)**2).sum(axis=1).max())
        # Scale if max step is larger than alpha
        if max_step_length > 1.0:
            dr = f * self.alpha / max_step_length
        else:
            dr = f * self.alpha

        self.r0 = (r + dr).copy()

        self.f0 = f.reshape(-1).copy()

        atoms.set_positions(self.r0)
        self.dump((self.r0, self.f0, self.maxstep, self.alpha))

    def replay_trajectory(self, traj):
        """Initialize hessian from old trajectory."""
        if isinstance(traj, str):
            from ase.io.trajectory import Trajectory
            traj = Trajectory(traj, 'r')

        atoms = traj[0]
        r0 = atoms.get_positions().ravel()
        f0 = atoms.get_forces().ravel()
        for atoms in traj:
            r = atoms.get_positions().ravel()
            f = atoms.get_forces().ravel()
            r0 = (r + f).copy()
            f0 = f.copy()

        self.r0 = r0
        self.f0 = f0
