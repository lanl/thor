# Tensor Train Configurational Integral

This repository contains code to compute energy, pressure, and phase diagrams of atomic systems using the **Tensor-Train Configurational Integral (TTCI)** method.

## Dependencies

You will need:

- LAMMPS (with Python bindings)
- ASE (Atomic Simulation Environment)
- Potential files (e.g. .include, .eam)
- Structure files in LAMMPS format (e.g. .data, .xyz → converted to .data)

---

## Reproduce Plots

### Cu with SMATB Energy and Pressure

    python Fig3_Cu.p

### Ar (HIPNN model) Energy and Pressure

    python Fig4_Ar_HIPNN.py

### Sn Phase Diagram (TTCI vs MD)

    python Fig5_phase_diagram.py

---

## Notes

- TTCI rank is controlled by the `--inds` argument (e.g., `-1 1` → rank 2)
- Quadrature rule controlled via `--method` (e.g., `Gauss3`, `Gauss8`)
- Ensure LAMMPS and ASE are correctly configured to find potentials and structure files
