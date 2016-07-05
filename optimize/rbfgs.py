# -*- coding: utf-8 -*-
import numpy as np
from numpy.linalg import eigh

from ase.optimize.optimize import Optimizer
from ase.build import orthogonal_procrustes

from scipy.linalg import block_diag

class RBFGS(Optimizer):
    def __init__(self, atoms, restart=None, logfile='-', trajectory=None,
                 maxstep=None, master=None, procrustes=True):
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
        self.maxstep = 0.1
        self.beta = 0.5
        self.prev = (float("inf"),None,None,None)

    def read(self):
        self.H, self.r0, self.f0, self.maxstep = self.load()

    def step(self, f):
        # Get variables
        atoms = self.atoms
        procrustes = atoms.procrustes
        r_full = atoms.get_positions_full().reshape((-1, atoms.natoms, 3))
        r = r_full[1:-1].copy().reshape((-1,3))
        r_full = r_full

        fmax = np.sqrt((f**2).sum(axis=1).max())
        f = f.reshape(-1)

        # Backtrack if needed
        if fmax > self.prev[0]:
            r,f,self.H = self.prev[1].copy(), self.prev[2].copy(), self.prev[3].copy() if self.prev[3] is not None else None
            #r,f,self.H = self.prev[1].copy(), self.prev[2].copy(), None
            self.maxstep *= self.beta
            #print "Backtrack to %lg with sum(H) = %lg" % (self.maxstep, self.H.sum() if self.H is not None else -1)
        else:
            self.prev = (fmax, r.copy(), f.copy(), self.H if self.H is not None else None)
            #print "Stored %lg with sum(H) = %lg" % (self.maxstep, self.H.sum() if self.H is not None else -1)

        # Update inv hessian
        self.update(r.flat, f, self.r0, self.f0)
        
        # Take step
        omega, V = eigh(self.H)
        dr = np.dot(V, np.dot(f, V) / np.fabs(omega)).reshape((-1, 3))
        steplengths = (dr**2).sum(1)**0.5
        dr = self.determine_step(dr, steplengths)


        #atoms.set_positions(r + dr)

        # Rotate coordinates
        full_rotation = np.array([np.empty((3,3)) for i in range((atoms.nimages-2) * atoms.natoms)])
        
        #images = atoms.get_positions_full().reshape((-1, atoms.natoms, 3))[:-1]
        images = np.append([r_full[0]], r+dr).reshape((-1, atoms.natoms, 3))

        for i in range(1,len(images)):
            rotation_matrix, _ = orthogonal_procrustes(images[i], images[i-1])
            for j,atom in enumerate(images[i]):
                full_rotation[(i-1)*3+j] = rotation_matrix
        R = block_diag(*full_rotation)

        r = np.dot(r.flatten(),R)
        dr = np.dot(dr.flatten(),R)
        new_pos = (r+dr).reshape((-1,3))
        atoms.set_positions(new_pos)

        # Rotate the new and old gradient, as well as old positions
        
        new_forces = np.dot(f.flatten(),R)
        if self.r0 is not None:
            old_forces = np.dot(self.f0,R)
            old_pos = np.dot(self.r0.flatten(),R)

        # Rotate the hessian
        self.H = np.dot(np.dot(R,self.H), R.T)
    
        # Store previous results
        self.r0 = r.flatten().copy()
        self.f0 = new_forces.copy()

        self.dump((self.H, self.r0, self.f0, self.maxstep))

    def step_1(self, f):

        # Get variables
        atoms = self.atoms
        procrustes = atoms.procrustes
        r = atoms.get_positions()
        f = f.reshape(-1)

        # If our first time, initialize H to identity
        if self.H is None:
            self.update(r.flat, f, self.r0, self.f0)

        # Take step towards minimum
        omega, V = eigh(self.H)
        dr = np.dot(V, np.dot(f, V) / np.fabs(omega)).reshape((-1, 3))
        steplengths = (dr**2).sum(1)**0.5
        dr = self.determine_step(dr, steplengths)
        atoms.set_positions(r + dr)

        # Rotate coordinates
        full_rotation = np.array([np.empty((3,3)) for i in range((atoms.nimages-2) * atoms.natoms)])
        images = atoms.get_positions_full()
        images = images.reshape((-1, atoms.natoms, 3))[:-1]
        for i in range(1,len(images)):
            rotation_matrix, _ = orthogonal_procrustes(images[i], images[i-1])
            for j,atom in enumerate(images[i]):
                images[i][j] = np.dot(atom,rotation_matrix)
                full_rotation[(i-1)*3+j] = rotation_matrix
        new_pos = images[1:].flatten().reshape((-1,3))
        atoms.set_positions(new_pos)

        # Rotate the new and old gradient, as well as old positions
        R = block_diag(*full_rotation)
        new_forces = np.dot(f.flatten(),R)
        if self.r0 is not None:
            old_forces = np.dot(self.f0,R)
            old_pos = np.dot(self.r0.flatten(),R)

        # Rotate the hessian
        self.H = np.dot(np.dot(R,self.H), R.T)

        if self.r0 is not None:
            self.update(new_pos.flatten(), new_forces, old_pos, old_forces)
        self.r0 = new_pos.flatten().copy()
        self.f0 = new_forces.copy()

        self.dump((self.H, self.r0, self.f0, self.maxstep))

    def step_2(self, f):
        # Get variables
        atoms = self.atoms
        procrustes = atoms.procrustes
        r = atoms.get_positions()
        f = f.reshape(-1)

        self.update(r.flat, f, self.r0, self.f0)
        
        omega, V = eigh(self.H)
        dr = np.dot(V, np.dot(f, V) / np.fabs(omega)).reshape((-1, 3))
        steplengths = (dr**2).sum(1)**0.5
        dr = self.determine_step(dr, steplengths)
        atoms.set_positions(r + dr)

        # ROTATE
        full_rotation = np.array([np.empty((3,3)) for i in range((atoms.nimages-2) * atoms.natoms)])
        images = atoms.get_positions_full().reshape((-1, atoms.natoms, 3))[:-1]
        for i in range(1,len(images)):
            rotation_matrix, _ = orthogonal_procrustes(images[i], images[i-1])
            for j,atom in enumerate(images[i]):
                images[i][j] = np.dot(atom,rotation_matrix)
                full_rotation[(i-1)*3+j] = rotation_matrix
        new_pos = images[1:].flatten().reshape((-1,3))
        atoms.set_positions(new_pos)

        # Rotate the new and old gradient, as well as old positions
        R = block_diag(*full_rotation)
        new_forces = np.dot(f.flatten(),R)
        if self.r0 is not None:
            old_forces = np.dot(self.f0,R)
            old_pos = np.dot(self.r0.flatten(),R)

        # Rotate the hessian
        self.H = np.dot(np.dot(R,self.H), R.T)
        
        self.r0 = new_pos.flatten().copy()
        self.f0 = new_forces.copy()
        self.dump((self.H, self.r0, self.f0, self.maxstep))


    def determine_step(self, dr, steplengths):
        """Determine step to take according to maxstep
        
        Normalize all steps as the largest step. This way
        we still move along the eigendirection.
        """
        maxsteplength = np.max(steplengths)
        #if maxsteplength >= self.maxstep:
        dr *= self.maxstep / maxsteplength
        
        return dr

    def update(self, r, f, r0, f0):
        if self.H is None:
            self.H = np.eye(3 * len(self.atoms)) * 70.0
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

class oldRBFGS(RBFGS):
    def determine_step(self, dr, steplengths):
        """Old RBFGS behaviour for scaling step lengths

        This keeps the behaviour of truncating individual steps. Some might
        depend of this as some absurd kind of stimulated annealing to find the
        global minimum.
        """
        dr /= np.maximum(steplengths / self.maxstep, 1.0).reshape(-1, 1)
        return dr
