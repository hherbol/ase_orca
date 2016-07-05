import numpy as np


def rotation_matrix_from_points(m0, m1):
    """Returns a rigid transformation/rotation matrix that minimizes the
    RMSD between two set of points.
    
    m0 and m1 should be (3, npoints) numpy arrays with
    coordinates as columns::

        (x1  x2   x3   ... xN
         y1  y2   y3   ... yN
         z1  z2   z3   ... zN)

    The centeroids should be set to origin prior to
    computing the rotation matrix.

    The rotation matrix is computed using quaternion
    algebra as detailed in::
        
        Melander et al. J. Chem. Theory Comput., 2015, 11,1055
    """

    v0 = np.copy(m0)
    v1 = np.copy(m1)

    # compute the rotation quaternion

    R11, R22, R33 = np.sum(v0 * v1, axis=1)
    R12, R23, R31 = np.sum(v0 * np.roll(v1, -1, axis=0), axis=1)
    R13, R21, R32 = np.sum(v0 * np.roll(v1, -2, axis=0), axis=1)

    f = [[R11 + R22 + R33, R23 - R32, R31 - R13, R12 - R21],
         [R23 - R32, R11 - R22 - R33, R12 + R21, R13 + R31],
         [R31 - R13, R12 + R21, -R11 + R22 - R33, R23 + R32],
         [R12 - R21, R13 + R31, R23 + R32, -R11 - R22 + R33]]

    F = np.array(f)

    w, V = np.linalg.eigh(F)
    # eigenvector corresponding to the most
    # positive eigenvalue
    q = V[:, np.argmax(w)]

    # Rotation matrix from the quaternion q

    R = quaternion_to_matrix(q)

    return R

    
def quaternion_to_matrix(q):
    """Returns a rotation matrix.
    
    Computed from a unit quaternion Input as (4,) numpy array.
    """

    q0, q1, q2, q3 = q
    R_q = [[q0**2 + q1**2 - q2**2 - q3**2,
            2 * (q1 * q2 - q0 * q3),
            2 * (q1 * q3 + q0 * q2)],
           [2 * (q1 * q2 + q0 * q3),
            q0**2 - q1**2 + q2**2 - q3**2,
            2 * (q2 * q3 - q0 * q1)],
           [2 * (q1 * q3 - q0 * q2),
            2 * (q2 * q3 + q0 * q1),
            q0**2 - q1**2 - q2**2 + q3**2]]
    return np.array(R_q)


def minimize_rotation_and_translation(target, atoms):
    """Minimize RMSD between atoms and target.
    
    Rotate and translate atoms to best match target.  For more details, see::
        
        Melander et al. J. Chem. Theory Comput., 2015, 11,1055
    """

    p = atoms.get_positions()
    p0 = target.get_positions()

    # centeroids to origin
    c = np.mean(p, axis=0)
    p -= c
    c0 = np.mean(p0, axis=0)
    p0 -= c0

    # Compute rotation matrix
    R = rotation_matrix_from_points(p.T, p0.T)

    atoms.set_positions(np.dot(p, R.T) + c0)

def orthogonal_procrustes(A, ref_matrix, reflection=False):
    # Adaptation of scipy.linalg.orthogonal_procrustes -> https://github.com/scipy/scipy/blob/v0.16.0/scipy/linalg/_procrustes.py#L14
    # Info here: http://compgroups.net/comp.soft-sys.matlab/procrustes-analysis-without-reflection/896635
    # goal is to find unitary matrix R with det(R) > 0 such that ||A*R - ref_matrix||^2 is minimized
    from scipy.linalg.decomp_svd import svd # Singular Value Decomposition, factors matrices
    from scipy.linalg import det
    import numpy as np

    A = np.asarray_chkfinite(A)
    ref_matrix = np.asarray_chkfinite(ref_matrix)

    if A.ndim != 2:
        raise ValueError('expected ndim to be 2, but observed %s' % A.ndim)
    if A.shape != ref_matrix.shape:
        raise ValueError('the shapes of A and ref_matrix differ (%s vs %s)' % (A.shape, ref_matrix.shape))


    u, w, vt = svd(ref_matrix.T.dot(A).T)

    # Goal: minimize ||A*R - ref||^2, switch to trace
    # trace((A*R-ref).T*(A*R-ref)), now we distribute
    # trace(R'*A'*A*R) + trace(ref.T*ref) - trace((A*R).T*ref) - trace(ref.T*(A*R)), trace doesn't care about order, so re-order
    # trace(R*R.T*A.T*A) + trace(ref.T*ref) - trace(R.T*A.T*ref) - trace(ref.T*A*R), simplify
    # trace(A.T*A) + trace(ref.T*ref) - 2*trace(ref.T*A*R)
    # Thus, to minimize we want to maximize trace(ref.T * A * R)

    # u*w*v.T = (ref.T*A).T
    # ref.T * A = w * u.T * v
    # trace(ref.T * A * R) = trace (w * u.T * v * R)
    # differences minimized when trace(ref.T * A * R) is maximized, thus when trace(u.T * v * R) is maximized
    # This occurs when u.T * v * R = I (as u, v and R are all unitary matrices so max is 1)
    # R is a rotation matrix so R.T = R^-1
    # u.T * v * I = R^-1 = R.T
    # R = u * v.T
    # Thus, R = u.dot(vt)

    R = u.dot(vt) # Get the rotation matrix, including reflections
    if not reflection and det(R) < 0: # If we don't want reflection
        # To remove reflection, we change the sign of the rightmost column of u (or v) and the scalar associated
        # with that column
        u[:,-1] *= -1
        w[-1] *= -1
        R = u.dot(vt)

    scale = w.sum() # Get the scaled difference

    return R,scale