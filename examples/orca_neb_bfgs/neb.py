from ase.io import read
from ase.neb import NEB
from ase.calculators.orca import orca
from ase.optimize import BFGS

#from ase.io.extxyz import write_extxyz as write_xyz

# Read in frames
start, images = [], []
FPTR = "../xyz/CNH_HCN.xyz"
#FPTR = "xyz/CNH_HCN_opt.xyz"
for a in read(FPTR,':'):
    start += [a]
    images += [a.copy()]

# Setup neb
neb = NEB(images, parallel=False, remove_rotation_and_translation=True)
# Set the calculator
for image in images:
    image.set_calculator(orca(label="CNH_HCN", RunTyp="Gradient"))

# Run neb
dyn = BFGS(neb, trajectory='CNH_HCN.traj')
dyn.run(fmax=0.1)

# Save?
for image in images:
    print(image.get_distance(1, 2), image.get_potential_energy())
