import os, sys
import argparse
import numpy as np
from ase import io, units
from ttci import make_LMP_script, compute_LMP_L2  # assumes ttci.py includes these

def compute_Lxyz(lmp_data_file, potfile, mass, T, dx=2.0, numcell=50):
	"""Compute log partition terms in x, y, z directions using L2 quadrature."""
	L = []
	for j in range(3):
		lmp_script = make_LMP_script(lmp_data_file, potfile, mass)
		Ldim, _ = compute_LMP_L2(lmp_script, T, dx, numcell, dim=[0, j])
		L.append(Ldim)
	return np.asarray(L)

def run_phase_diagram(structures, potfile, mass, presrange, temprange, ncell_dict, input_dir, temp_dir, output_file, dx=2.0, numcell=50, sp=1e-6):
	kB = units.kB

	Structs = list(structures.keys())
	P = np.zeros((len(Structs), len(presrange), len(temprange)))
	F = np.zeros_like(P)
	G = np.zeros_like(P)
	V = np.zeros_like(P)
	N_atoms = np.zeros_like(P)

	for iden, p in enumerate(presrange):
		for itemp, T in enumerate(temprange):
			beta = 1 / (kB * T)

			for ist, struct in enumerate(Structs):
				basename = f"{struct}_p{p:03d}"
				structfile = os.path.join(input_dir, f"{basename}.xyz")
				CF = io.read(structfile)

				V1 = CF.get_volume()
				N = CF.get_global_number_of_atoms()
				N_atoms[ist, iden, itemp] = N

				temp_lmp_file = os.path.join(temp_dir, f"{basename}_{N}.data")
				os.makedirs(temp_dir, exist_ok=True)
				io.write(temp_lmp_file, CF, format='lammps-data')

				L, U0 = compute_Lxyz(temp_lmp_file, potfile, mass, T, dx, numcell)
				lnQ = -beta * U0 + N * np.sum(np.log(L))
				F[ist, iden, itemp] = lnQ

				# Pressure computation
				CF.set_cell(CF.get_cell() * (1 + sp), scale_atoms=True)
				V2 = CF.get_volume()
				temp_lmp_file_2 = os.path.join(temp_dir, f"CF2_{basename}.data")
				io.write(temp_lmp_file_2, CF, format='lammps-data')

				L2, U2 = compute_Lxyz(temp_lmp_file_2, potfile, mass, T, dx, numcell)
				lnQ2 = -beta * U2 + N * np.sum(np.log(L2))
				Pcalc = kB * T * (lnQ2 - lnQ) / (V2 - V1)
				P[ist, iden, itemp] = Pcalc
				V[ist, iden, itemp] = V1

				# Gibbs Energy
				G[ist, iden, itemp] = (-lnQ + V1 * Pcalc) / N

	np.savez(output_file, F=F, G=G, P=P, V=V, presrange=presrange, temprange=temprange, structs=Structs, ncell=ncell_dict, mass=mass, sp=sp)

def parse_args():
	parser = argparse.ArgumentParser(description="Compute phase diagram using TTCI")

	parser.add_argument('--potfile', required=True, help='Path to LAMMPS potential file')
	parser.add_argument('--mass', type=float, required=True, help='Atomic mass of the element')
	parser.add_argument('--input_dir', required=True, help='Directory with structure .xyz files')
	parser.add_argument('--output_file', required=True, help='Path to save .npz results')
	parser.add_argument('--temp_dir', default='./tempfile/', help='Directory to store temporary LAMMPS files')

	parser.add_argument('--pres', type=int, nargs=3, metavar=('start', 'stop', 'step'), default=[20, 60, 1],
    help='Pressure range in units of 0.1 GPa, e.g. 20 60 1 => 2.0 to 6.0 GPa')

	parser.add_argument('--temp', type=int, nargs=3, metavar=('start', 'stop', 'step'), default=[100, 600, 20],
    help='Temperature range in Kelvin')

	parser.add_argument('--structures', nargs='+', required=True,
											help='List of structure names (must match basename in input_dir)')
	parser.add_argument('--ncells', nargs='+', type=int, required=True,
											help='Matching list of supercell sizes for each structure')

	return parser.parse_args()

if __name__ == "__main__":
	args = parse_args()

	assert len(args.structures) == len(args.ncells), "Each structure must have a matching ncell size."

	struct_dict = dict(zip(args.structures, args.ncells))
	presrange = range(*args.pres)
	temprange = range(*args.temp)

	run_phase_diagram(
		structures=struct_dict,
		potfile=args.potfile,
		mass=args.mass,
		presrange=presrange,
		temprange=temprange,
		ncell_dict=struct_dict,
		input_dir=args.input_dir,
		temp_dir=args.temp_dir,
		output_file=args.output_file
	)