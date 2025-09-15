# ttci.py

import os
import numpy as np
from ase import io, units
from lammps import lammps, LAMMPS_INT, LMP_STYLE_GLOBAL
import quadrature_rules as quadclass  # assumes methods like Gauss8 exist

# ------------------------
# Quadrature Helpers
# ------------------------
def get_w_and_grid(method, dx, num_interval):
	qd = getattr(quadclass, method)()
	W, disp_grid = qd.get_weight_and_grid(-dx, dx, num_interval)
	return W, disp_grid


# ------------------------
# LAMMPS SCRIPT UTILITIES
# ------------------------
def make_LMP_script(structfile, potentialpath, mass):
	return f"""
	units           metal
	boundary        p p p
	atom_modify     map hash
	read_data       {structfile} nocoeff
	mass            1 {mass:f}
	include         {potentialpath}
	"""

def compute_LMP_potential_energy(lmp, coord):
	lmp.scatter_atoms("x", 1, 3, coord.ctypes)
	lmp.command('run 0')
	return float(lmp.extract_compute("thermo_pe", LMP_STYLE_GLOBAL, LAMMPS_INT))

def LMP_U_new(lmp, disp, Q0, U0):
	return compute_LMP_potential_energy(lmp, disp + Q0) - U0


# ------------------------
# TT-CUR CORE DECOMPOSITION
# ------------------------
def compute_Qval_3cores(L, alpha, d):
	Lval = alpha * L[0]
	for i in range(1, d - 1):
		Lval = alpha * Lval @ L[i]
	Lval = alpha * Lval @ L[d - 1]
	return np.squeeze(Lval)


def compute_tt_Qval(fun, grid, l2rindices, r2lindices, Ns, w, alpha=None):
	Gs = compute_TT_CUR_cores(single_to_multi(fun), grid, l2rindices, r2lindices, Ns)
	L = [np.tensordot(w, x, (0, 1)) for x in Gs]
	d = len(grid)

	if alpha is None:
		alpha = 1.0
		Qval = compute_Qval_3cores(L, alpha, d)
		if abs(Qval) > 1e3:
			asign = -1
		elif abs(Qval) < 1e-3:
			asign = 1
		while abs(Qval) <= 1.0:
			alpha += asign * 0.01
			Qval = compute_Qval_3cores(L, alpha, d)
	else:
		Qval = compute_Qval_3cores(L, alpha, d)

	return abs(Qval), alpha


# ------------------------
# TTCI ENERGY & PRESSURE
# ------------------------
def ttci_lmp_energy(structfile, curT, potfile, mass, dx, numcell, method, r, inds_step):
	args = ["-screen", "none", "-log", "none"]
	lmpScript = make_LMP_script(structfile, potfile, mass)
	lmp = lammps(cmdargs=args)
	lmp.commands_string(lmpScript)

	Q0 = np.array(lmp.gather_atoms("x", 1, 3))
	N = int(Q0.shape[0] / 3)
	U0 = compute_LMP_potential_energy(lmp, Q0)
	beta = 1 / (units.kB * curT)

	w, disp_grid = get_w_and_grid(method, dx, numcell)
	i0 = np.argmin(np.abs(disp_grid))

	inds = [i0 + x for x in inds_step]
	d = 3 * N
	grid = [disp_grid] * d
	Ns = [len(disp_grid)] * d
	ranks = [1] + [len(inds)] * (d - 1) + [1]

	l2rindices = nested_indices(Ns, ranks, inds)
	r2lindices = nested_indices(Ns, ranks, inds, -1)

	normfactor = 100.0

	def fun1(x):
		Udisp = LMP_U_new(lmp, np.squeeze(x), Q0, U0)
		return np.exp(-beta * Udisp + normfactor)

	Qval1, alpha1 = compute_tt_Qval(fun1, grid, l2rindices, r2lindices, Ns, w)

	beta2 = beta + 1e-3

	def fun2(x):
		Udisp = LMP_U_new(lmp, np.squeeze(x), Q0, U0)
		return np.exp(-beta2 * Udisp + normfactor)

	Qval2, alpha2 = compute_tt_Qval(fun2, grid, l2rindices, r2lindices, Ns, w)

	lnQ1 = np.log(Qval1) - d * np.log(alpha1)
	lnQ2 = np.log(Qval2) - d * np.log(alpha2)

	E = 1.5 * N * units.kB * curT + U0 - (lnQ2 - lnQ1) / (beta2 - beta)
	return E

def ttci_lmp_pressure(structfile, basename, potfile, mass, curT, dx, sp, numcell, method, r, inds_step):
	w, disp_grid = get_w_and_grid(method, dx, numcell)
	i0 = np.argmin(np.abs(disp_grid))

	CF = io.read(structfile, format='lammps-data')
	Q1 = CF.get_positions()
	V1 = CF.get_volume()
	beta = 1 / (units.kB * curT)

	inds = [i0 + x for x in inds_step]
	d = 3 * Q1.shape[0]

	grid = [disp_grid] * d
	Ns = [len(disp_grid)] * d
	ranks = [1] + [len(inds)] * (d - 1) + [1]
	l2rindices = nested_indices(Ns, ranks, inds)
	r2lindices = nested_indices(Ns, ranks, inds, -1)

	lmpScript1 = make_LMP_script(structfile, potfile, mass)
	lmp1 = lammps(cmdargs=["-screen", "none", "-log", "none"])
	lmp1.commands_string(lmpScript1)
	U1 = compute_LMP_potential_energy(lmp1, Q1)

	def fun1(x):
		return np.exp(-beta * LMP_U_new(lmp1, np.reshape(x, (-1, 3)), Q1, U1))

	Qval1, alpha1 = compute_tt_Qval(fun1, grid, l2rindices, r2lindices, Ns, w)

	CF.set_cell(CF.get_cell() * (1 - sp), scale_atoms=True)
	Q2 = CF.get_positions()
	V2 = CF.get_volume()

	structfile2 = f'./tempfile/CF2_{basename}.data'
	os.makedirs(os.path.dirname(structfile2), exist_ok=True)
	io.write(structfile2, CF, format='lammps-data')

	lmpScript2 = make_LMP_script(structfile2, potfile, mass)
	lmp2 = lammps(cmdargs=["-screen", "none", "-log", "none"])
	lmp2.commands_string(lmpScript2)
	U2 = compute_LMP_potential_energy(lmp2, Q2)

	def fun2(x):
		return np.exp(-beta * LMP_U_new(lmp2, np.reshape(x, (-1, 3)), Q2, U2))

	Qval2, alpha2 = compute_tt_Qval(fun2, grid, l2rindices, r2lindices, Ns, w)

	lnQval1 = np.log(Qval1) - d * np.log(alpha1)
	lnQval2 = np.log(Qval2) - d * np.log(alpha2)

	P = -(U2 - U1) / (V2 - V1) + 1 / beta * (lnQval2 - lnQval1) / (V2 - V1)
	P *= units.eV / (0.1 * units.GPa)
	return P


# ------------------------
# PHASE DIAGRAM UTILITIES
# ------------------------

def compute_LMP_L2(lmp_script, T, dx, numcell, dim=[0, 0], method='Gauss8'):
	args = ["-screen", "none", "-log", "none"]
	lmp = lammps(cmdargs=args)
	lmp.commands_string(lmp_script)

	atomdimIdx = dim[0] * 3 + dim[1]

	Q0 = np.array(lmp.gather_atoms("x", 1, 3))
	U0 = compute_LMP_potential_energy(lmp, Q0)
	beta = 1 / (units.kB * T)

	W, disp_grid = get_w_and_grid(method, dx, numcell)

	temp_disp = np.zeros_like(Q0)
	Li = np.zeros_like(disp_grid)

	for l, delta in enumerate(disp_grid):
		temp_disp[atomdimIdx] = delta
		Li[l] = np.exp(-beta * LMP_U_new(lmp, temp_disp, Q0, U0))

	Lval = np.dot(W, Li)
	return [Lval, U0]


def compute_Lxyz(lmp_data_file, potfile, mass, T, dx=2.0, numcell=50):
	L = []
	for j in range(3):
		lmp_script = make_LMP_script(lmp_data_file, potfile, mass)
		Ldim, _ = compute_LMP_L2(lmp_script, T, dx, numcell, dim=[0, j])
		L.append(Ldim)
	return np.asarray(L)


def compute_lnQ(U0, Lvec, N, beta):
	return -beta * U0 + N * np.sum(np.log(Lvec))


def compute_pressure_L2(lnQ1, lnQ2, V1, V2, T):
	kB = units.kB
	P = kB * T * (lnQ2 - lnQ1) / (V2 - V1)
	return P * units.eV / (0.1 * units.GPa)


def compute_gibbs_energy(lnQ, V, P, N):
	return (-lnQ + V * P) / N

# ----------------
# Math Utilities
# ----------------
def nested_indices(Ns, ranks, inds, direction=1):
	if direction == -1:
		return [np.flipud(x) for x in nested_indices(Ns[::-1], ranks[::-1], inds, direction=1)[::-1]]
	else:
		indices = [np.empty(shape=(0, 0), dtype='int64')]
		for n, r, rprev in zip(Ns[:-1], ranks[1:-1], ranks[0:-2]):
			indices.append(l2rnestindices(inds, indices[-1], (rprev, n)))
		indices.append(np.empty(shape=(len(Ns), 0), dtype='int64'))
		return indices

def grideval(grid, indices):
	assert len(grid) == len(indices)
	return [g[i] for g, i in zip(grid, indices)]

def meshindices(indices):
	s0 = (1,) * len(indices)
	return [c.reshape(s0[:i] + (-1,) + s0[i + 1:]) for i, inds in enumerate(indices) if inds.size > 0 for c in np.atleast_2d(inds)]

def l2rnestindices(unnested_indices, prev_inds, shape):
	tmp = np.unravel_index(unnested_indices, shape)
	prev_inds = prev_inds.reshape(prev_inds.shape[0], shape[0])
	nested_indices = np.vstack([prev_inds[:, tmp[0]], tmp[1][None, :]])
	return nested_indices

def single_to_multi(single_point_function):
	def inner(x):
		y = np.array(np.broadcast_arrays(*x))
		return np.apply_along_axis(single_point_function, 0, y)
	return inner

def compute_TT_CUR_cores(fun, grid, l2rindices, r2lindices, Ns):
	d = len(grid)
	maxvol_mats = [fun(grideval(grid, meshindices([l2rindices[i+1], r2lindices[i+1]]))) for i in range(1)]
	Vol = np.abs(np.linalg.det(maxvol_mats[0]))
	print(f'CUR Volume = {Vol:0.3e}')
	fibers = [fun(grideval(grid, meshindices([l2rindices[i], np.arange(Ns[i]), r2lindices[i+1]]))) for i in [0, 1, d-1]]
	cores = [np.linalg.lstsq(maxvol_mats[0].T, f.T.reshape(f.shape[2], -1), rcond=None)[0].reshape(f.shape[2], f.shape[1], f.shape[0]).T for f in fibers[:2]] + [fibers[2]]
	return cores