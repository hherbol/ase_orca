ALPHA = 0.1
FMAX = 0.1
RIGID_ROTATION = True


def method_ASE():
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
    FPTR = "../../xyz/CNH_HCN.xyz"

    for a in read(FPTR, ':'):
        start += [a]
        images += [a.copy()]

    # Setup neb
    neb = NEB(images, parallel=False,
              remove_rotation_and_translation=RIGID_ROTATION)
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

    run_name = 'CNH_HCN_clancelot'
    new_opt_params = {'step_size': ALPHA,
                      'step_size_adjustment': 0.5,
                      'max_step': 0.2,
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
    opt.optimize()
    write_xyz(frames, "CNH_HCN_opt_clancelot")

    print("\nDONE WITH CLANCELOT SIMULATION...\n")


method_ASE()
method_CLANCELOT(opt_method="SD")
