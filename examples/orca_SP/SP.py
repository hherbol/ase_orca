from ase.io import read, write
from ase.calculators.orca import orca
from ase.optimize import BFGS


# Optimize using ASE optimization method and single point calculations
def method_1():
    print("Running Method 1 - Optimization via ASE optimizer.")
    # Read in second to last frame
    image = read("../xyz/CNH_HCN.xyz", '-2')

    # Set the calculator
    image.set_calculator(orca(label="CNH", RunTyp="Gradient"))

    # Run optimizer
    dyn = BFGS(image, trajectory='CNH.traj')
    dyn.run(fmax=0.1)

    # Final output
    print "Energy via Method 1 = ", image.get_potential_energy(), "eV"
    # Resave trajectory to xyz
    traj = read("CNH.traj", ':')
    write("CNH.xyz", traj)


# Use OPT in orca
def method_2():
    print("Running Method 2 - Optimization via Orca.")
    # Read in second to last frame
    image = read("../xyz/CNH_HCN.xyz", '-2')

    # Set the calculator
    image.set_calculator(orca(label="CNH", RunTyp="Opt"))

    # Run the simulation to get the potential energy
    # Note - it does not recalculate anything as the system itself
    #        hasn't changed between calls to get_potential_energy()
    image.get_potential_energy()
    print "Energy via Method 2 = ", image.get_potential_energy(), "eV"


method_1()
method_2()
