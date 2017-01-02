import print_helper

ALPHA = 0.1
FMAX = 0.03
FRMS = 0.001
MAX_STEP = 0.2  # 0.2 is better, but to keep with ASE we use the same
RIGID_ROTATION = True


def method_ASE(RIGID_ROTATION, FPTR="../xyz/CNH_HCN.xyz"):
    from ase.io import read, write
    from ase.neb import NEB
    from ase.calculators.orca import orca
    from ase.optimize import SD
    from os import mkdir, chdir
    from os.path import isdir

    # Read in frames
    if not isdir("ASE"):
        mkdir("ASE")
    chdir("ASE")

    start, images = [], []
    FPTR = "../" + FPTR
    if not FPTR.endswith(".xyz"):
        FPTR += ".xyz"

    for a in read(FPTR, ':'):
        start += [a]
        images += [a.copy()]

    # Setup neb
    neb = NEB(images, parallel=False,
              remove_rotation_and_translation=RIGID_ROTATION,
              method="improvedtangent")
    # Set the calculator
    for image in images:
        image.set_calculator(orca(label="CNH_HCN", RunTyp="Gradient"))

    # Run neb
    dyn = SD(neb, trajectory='CNH_HCN.traj')
    dyn.run(fmax=FMAX)

    # Save?
    traj = read("CNH_HCN.traj", ':-1')
    write("../CNH_HCN_opt_ASE.xyz", traj)

    chdir("..")
    print("\nDONE WITH ASE SIMULATION...\n")


def method_CLANCELOT(opt_method="LBFGS"):
    from files import read_xyz, write_xyz
    import neb
    from units import convert, convert_energy

    FPTR = "./../xyz/CNH_HCN.xyz"
    frames = read_xyz(FPTR)
    route = '! HF-3c'

    if RIGID_ROTATION:
        is_on = "ON"
    else:
        is_on = "OFF"

    print("\nRUNNING CLANCELOT SIMULATION WITH RIGID_ROTATION %s...\n" % is_on)

    run_name = 'CNH_HCN_c_' + opt_method
    new_opt_params = {'step_size': ALPHA,
                      'step_size_adjustment': 0.5,
                      'max_step': MAX_STEP,
                      'linesearch': 'backtrack',
                      'accelerate': True,
                      'reset_step_size': 5,
                      'g_rms': convert("eV/Ang", "Ha/Ang", 0.03),
                      'g_max': convert("eV/Ang", "Ha/Ang", FMAX)}

    opt = neb.NEB(run_name,
                  frames,
                  route,
                  k=convert_energy("eV", "Ha", 0.1),
                  opt=opt_method,
                  new_opt_params=new_opt_params)
    output = opt.optimize()
    frames = output[-1]

    write_xyz(frames, "CNH_HCN_opt_%s" % opt_method)

    print("\nDONE WITH CLANCELOT SIMULATION...\n")


# We run ASE, then Clancelot, then we take the last frame in clancelot and
# validate it within ASE. If it converges in 1 iteration, and Clancelot
# converged faster than ASE originally, then Clancelot is validated as being
# faster.
s = print_helper.color_set("VERIFYING CLANCELOT AGAINST ASE", "BLUE")
s = print_helper.color_set(s, "BOLD")

method_ASE(RIGID_ROTATION)
method_CLANCELOT(opt_method="SD")
print(s)
method_ASE(False, "CNH_HCN_opt_SD")
