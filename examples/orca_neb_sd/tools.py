from ase.io import read, write
from os import mkdir, chdir, system
from os.path import isdir

def traj_to_xyzs():
	if not isdir("ASE_traj"): mkdir("ASE_traj")
	chdir("ASE_traj")
	system("rm *.xyz")

	images = read("../ASE/CNH_HCN.traj",':')
	original = read("../../xyz/CNH_HCN.xyz",':')
	N = len(original)

	xyz_counter, i, frame = 0, 0, []
	for image in images:
		if i < N:
			frame.append(image)
			i += 1
		else:
			write("CNH_HCN_%d.xyz" % xyz_counter, frame)
			frame = []
			xyz_counter += 1
			i = 0
	write("CNH_HCN_%d.xyz" % xyz_counter, frame)
	chdir("..")

traj_to_xyzs()