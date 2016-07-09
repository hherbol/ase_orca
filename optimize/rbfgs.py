# -*- coding: utf-8 -*-
import numpy as np
from numpy.linalg import eigh

import ase.units
from ase.optimize.optimize import Optimizer
from ase.build import orthogonal_procrustes

from scipy.linalg import block_diag

class RBFGS(Optimizer):
    def __init__(self, atoms, restart=None, logfile='-', trajectory=None,
                 maxstep=None, master=None):
        """RBFGS optimizer.

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

    def initialize(self):
        self.H = None
        self.r0 = None
        self.f0 = None
        self.maxstep = 0.04
        self.alpha = 0.1
        self.beta = 0.1
        self.fnorm0 = None
        self.prev = (None, None, None, None, None)

    def read(self):
        self.H, self.r0, self.f0, self.maxstep = self.load()

    def step_ase(self, f):
        atoms = self.atoms
        r = atoms.get_positions()
        f = f.reshape(-1)
        self.update_hessian_bfgs(r.flat, f, self.r0, self.f0)
        omega, V = eigh(self.H)
	# I think this is a "Preconditioning" step... 
        dr = np.dot(V, np.dot(f, V) / np.fabs(omega)).reshape((-1, 3))
        steplengths = (dr**2).sum(1)**0.5
        dr = self.determine_step(dr, steplengths)
        atoms.set_positions(r + dr)
        self.r0 = r.flat.copy()
        self.f0 = f.copy()
        self.dump((self.H, self.r0, self.f0, self.maxstep))

    def step_clancelot(self, f):
        atoms = self.atoms
        r = atoms.get_positions()
        f = f.reshape(-1) / ase.units.Hartree
        self.update_inv_hessian_bfgs(r.flat, f, self.r0, self.f0)
        dr = np.dot(self.H,f).reshape((-1,3))
        steplengths = (dr**2).sum(1)**0.5
        dr = self.determine_step(dr, steplengths)
        atoms.set_positions(r + dr)
        self.r0 = r.flat.copy()
        self.f0 = f.copy()
        self.dump((self.H, self.r0, self.f0, self.maxstep))

    def step_clancelot_backtrack(self, f):
        atoms = self.atoms
        r = atoms.get_positions()
        f = f.reshape(-1) # / (f**2).sum(axis=1).max()**0.5
        backtrack = self.check_fnorm(np.linalg.norm(f))
        if backtrack:
            for current, previous in zip([r,f,self.H,self.r0,self.f0], self.prev):
                try:
                    current = previous.copy()
                except:
                    current = None
        else:
            H = self.H.copy() if self.H is not None else None
            r0 = self.r0.copy() if self.r0 is not None else None
            f0 = self.f0.copy() if self.f0 is not None else None
            self.prev = [r.copy(), f.copy(), H, r0, f0]
            self.update_inv_hessian_bfgs(r.flat, f, self.r0, self.f0)
        dr = np.dot(self.H,f).reshape((-1,3))
        dr = self.determine_step(dr, (dr**2).sum(1)**0.5, method=2)
        print self.alpha, (dr**2).sum(1).max()**0.5
        atoms.set_positions(r + dr)
        self.r0 = r.flat.copy()
        self.f0 = f.copy()
        self.dump((self.H, self.r0, self.f0, self.maxstep))

    def step_clancelot_SD_backtrack(self, f):
        atoms = self.atoms
        I = np.eye(3 * len(atoms))
        r = atoms.get_positions()
        f = f.reshape(-1) # / (f**2).sum(axis=1).max()**0.5
        backtrack = self.check_fnorm(np.linalg.norm(f))
        if backtrack:
            r = self.prev[0].copy()
            f = self.prev[1].copy()
        else:
            self.prev = [r.copy(), f.copy()]
        dr = np.dot(I,f).reshape((-1,3))
        dr = self.determine_step(dr, (dr**2).sum(1)**0.5, method=2)
        atoms.set_positions(r + dr)
        self.dump((r.copy(), f.copy(), self.alpha))

    step = step_clancelot_SD_backtrack

    def determine_step(self, dr, steplengths, method=1):
        """Determine step to take according to maxstep
        
        Normalize all steps as the largest step. This way
        we still move along the eigendirection.
        """
        if method == 2:
            maxsteplength = np.max(steplengths)
            dr *= self.alpha / maxsteplength
        else:
            maxsteplength = np.max(steplengths)
            if maxsteplength >= self.maxstep:
                dr *= self.maxstep / maxsteplength
        
        return dr

    def check_fnorm(self,fnorm):
        fnorm_flag = False
        if self.fnorm0 is not None and fnorm > self.fnorm0:
            self.alpha *= self.beta
            fnorm_flag = True
        self.fnorm0 = fnorm
        if 1e-7 > min([self.alpha, self.maxstep]):
            raise Exception("Lowered step size too far!")
        return fnorm_flag

    def update_hessian_bfgs(self, r, f, r0, f0):
        # https://en.wikipedia.org/wiki/Quasi-Newton_method
        if self.H is None:
            scalar = np.linalg.norm(f) / self.maxstep
            self.H = np.eye(3 * len(self.atoms)) * 70.0 # * scalar# * 70.0
            return
        dr = r - r0

        if np.abs(dr).max() < 1e-7:
            # Same configuration again (maybe a restart):
            return

        df = f - f0
        a = np.dot(dr, df)
        dg = np.dot(self.H, dr)
        b = np.dot(dr, dg)
        self.H -= np.outer(df, df) / a + np.outer(dg, dg) / b

    def update_inv_hessian_dfp(self, r, f, r0, f0):
        # A rapidly convergent descent method for minimization - Fletcher and Powell
        # http://bioinfo.ict.ac.cn/~dbu/AlgorithmCourses/Lectures/Fletcher-Powell.pdf
        if self.H is None:
            scalar = np.linalg.norm(f) / self.maxstep
            self.H = np.eye(3 * len(self.atoms))# * 70.0 # * scalar# * 70.0
            return
        dr = r - r0

        if np.abs(dr).max() < 1e-7:
            # Same configuration again (maybe a restart):
            return

        df = f - f0
        dg = np.dot(self.H, df)
        a = np.dot(df, dg)
        b = np.dot(dr, df)
        self.H -= np.outer(df, df) / a + np.outer(dr, dr) / b

    def update_inv_hessian_bfgs(self, r, f, r0, f0):
        if self.H is None:
            scalar = np.linalg.norm(f) / self.maxstep
            self.H = np.eye(3 * len(self.atoms)) * scalar# * 70.0
            return

        I = np.eye(3 * len(self.atoms))
        s = r - r0
        y = -(f - f0)

        if np.abs(s).max() < 1e-7:
            # Same configuration again (maybe a restart):
            return

        A1 = I - np.outer(s,y) / np.dot(y,s)
        A2 = I - np.outer(y,s) / np.dot(y,s)
        self.H = np.dot(A1, np.dot(self.H, A2)) + np.outer(s,s)/np.dot(y,s)

    def replay_trajectory(self, traj):
        """Initialize hessian from old trajectory."""
        if isinstance(traj, str):
            from ase.io.trajectory import Trajectory
            traj = Trajectory(traj, 'r')
        self.H = None
        atoms = traj[0]
        r0 = atoms.get_positions().ravel()
        f0 = atoms.get_forces().ravel()
        for atoms in traj:
            r = atoms.get_positions().ravel()
            f = atoms.get_forces().ravel()
            self.update(r, f, r0, f0)
            r0 = r
            f0 = f

        self.r0 = r0
        self.f0 = f0

    # Rotate coordinates
    def get_rotation(self, dr=None):
        atoms = self.atoms
        r_start = atoms.get_positions_full().reshape((-1, atoms.natoms, 3))
        r = r_start[1:-1].copy().reshape((-1,3))
        r_start = [r_start[0]]

        full_rotation = np.array([np.empty((3,3)) for i in range((atoms.nimages-2) * atoms.natoms)])
        if dr is not None:
            r += dr

        images = np.append([r_full[0]], r).reshape((-1, atoms.natoms, 3))
        for i in range(1,len(images)):
            rotation_matrix, _ = orthogonal_procrustes(images[i], images[i-1])
            for j,atom in enumerate(images[i]):
                full_rotation[(i-1)*3+j] = rotation_matrix
        
        return block_diag(*full_rotation)
