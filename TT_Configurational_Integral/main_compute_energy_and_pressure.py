# main_compute.py

import sys, os
import argparse
import numpy as np
from ttci import ttci_lmp_energy, ttci_lmp_pressure

def main():
    parser = argparse.ArgumentParser(description="Compute energy and pressure using TTCI method")

    # Required inputs
    parser.add_argument('--structfile', required=True, help='LAMMPS structure file (e.g., .data)')
    parser.add_argument('--potfile', required=True, help='LAMMPS potential file (e.g., .include)')
    parser.add_argument('--temp', type=float, required=True, help='Temperature in Kelvin')

    # Physical and system parameters
    parser.add_argument('--mass', type=float, default=63.546, help='Atomic mass (default: Cu = 63.546)')
    parser.add_argument('--atom_type', type=int, default=1, help='LAMMPS atom type number (default: 1)')

    # Volume perturbation for pressure calculation
    parser.add_argument('--dV', type=float, default=0.005, help='Volume shrinkage ratio for pressure computation')

    # Quadrature and TT parameters
    parser.add_argument('--dx', type=float, default=2.0, help='Displacement grid spacing')
    parser.add_argument('--numcell', type=int, default=51, help='Number of quadrature grid points')
    parser.add_argument('--method', type=str, default='Gauss3', help='Quadrature method (Gauss3, Gauss8, etc.)')

    # Tensor Train decomposition configuration
    parser.add_argument('--inds', type=int, nargs='+', default=[-1, 1],
				help='Displacement index steps for TT decomposition. '
						'The number of indices determines the TT rank.')

    args = parser.parse_args()

    structfile = args.structfile
    potfile = args.potfile
    T = args.temp
    mass = args.mass
    dx = args.dx
    numcell = args.numcell
    method = args.method
    dV = args.dV
    inds_step = args.inds
    r = len(inds_step)  # TT rank is determined by the length of inds_step

    basename = os.path.splitext(os.path.basename(structfile))[0]

    print(f"\nComputing for structure: {structfile} at temperature {T} K")
    print(f"TTCI rank: {r} (from inds_step = {inds_step})")

    print("Computing Energy...")
    energy = ttci_lmp_energy(structfile, T, potfile, mass, dx, numcell, method, r, inds_step)
    print(f"Energy: {energy:.6f} eV")

    print("Computing Pressure...")
    pressure = ttci_lmp_pressure(structfile, basename, potfile, mass, T, dx, dV, numcell, method, r, inds_step)
    print(f"Pressure: {pressure:.2f} bar")

if __name__ == "__main__":
    main()