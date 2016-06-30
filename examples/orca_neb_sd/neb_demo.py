from sys import argv

FMAX = 0.2

def method_ASE():
	from ase.io import read, write
	from ase.neb import NEB
	from ase.calculators.orca import orca
	from ase.optimize import SD
	from os import mkdir, chdir
	from os.path import isdir

	# Read in frames
	if not isdir("ASE"): mkdir("ASE")
	chdir("ASE")

	start, images = [], []
	FPTR = "../../xyz/CNH_HCN.xyz"
	#FPTR = "xyz/CNH_HCN_opt.xyz"
	for a in read(FPTR,':'):
	    start += [a]
	    images += [a.copy()]

	# Setup neb
	neb = NEB(images, parallel=False, remove_rotation_and_translation=False)
	# Set the calculator
	for image in images:
	    image.set_calculator(orca(label="CNH_HCN", RunTyp="Gradient"))

	# Run neb
	dyn = SD(neb, trajectory='CNH_HCN.traj')
	dyn.run(fmax=FMAX)

	# Save?
	traj = read("CNH_HCN.traj",':-1')
	write("CNH_HCN_opt_ASE.xyz",traj)

	chdir("..")
	print("\nDONE WITH ASE SIMULATION...\n")

def method_CLANCELOT():
	from files import read_xyz, write_xyz
	from neb import neb
	from units import convert, convert_energy

	FPTR = "./../xyz/CNH_HCN.xyz"
	frames = read_xyz(FPTR)
	route = '! HF-3c'

	run_name = 'CNH_HCN_clancelot'
	neb(run_name, frames, route, k=convert_energy("eV","Ha",0.1), opt="SD", maxiter=1000, gtol=convert("eV/Ang","Ha/Ang",0.03), fmax=convert("eV/Ang","Ha/Ang",FMAX), DFT='orca', alpha=0.1, fit_rigid=False, reset=100)
	write_xyz(frames,"CNH_HCN_opt_clancelot")

	print("\nDONE WITH CLANCELOT SIMULATION...\n")

if len(argv) > 1 and argv[1].lower() == "ase":
	print("Running ASE NEB of CNH Isomerization.\n")
	method_ASE()
elif len(argv) > 1 and argv[1].lower() == "clancelot":
	print("Running Clancelot NEB of CNH Isomerization.\n")
	method_CLANCELOT()
elif len(argv) > 1 and argv[1].lower() == "all":
	print("Running both ASE and Clancelot NEB of CNH Isomerization.\n")
	method_ASE()
	method_CLANCELOT()
else:
	print("Running default: ASE NEB of CNH Isomerization.")
	print("\tNote: Other options are (1) ase, (2) clancelot and (3) all.\n")
	method_ASE()