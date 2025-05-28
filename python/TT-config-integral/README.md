# Tensor Train Configurational Integral

This repository contains code to compute energy, pressure, and phase diagrams of atomic systems using the **Tensor-Train Configurational Integral (TTCI)** method.

It integrates:
- Tensor-train quadrature and CUR decomposition
- ASE and LAMMPS for structure and energy evaluation
- Scripts to reproduce figures and export data

---

## Directory Structure

code_for_submission/
├── data_and_figures/                   # Data files (.txt, .npz) and figures (.pdf)
│   ├── Cu_FCC_smatb_energy_pressure.txt
│   ├── Sn-phase_transition_data.txt
│   ├── phase_transition_plot.pdf
│   └── Fig_*.pdf (saved figures)
├── Fig3_Cu.py                          # Script: Cu energy & pressure comparison
├── Fig4_Ar_HIPNN.py                    # Script: Ar using HIPNN model
├── Fig5_phase_diagram.py              # Script: Sn MEAM phase diagram
├── main_compute_energy_and_pressure.py # CLI: compute TTCI energy and pressure
├── main_compute_phase_diagram.py       # CLI: compute phase diagrams via TTCI
├── quadrature_rules.py                # Gaussian quadrature logic (Gauss3, Gauss8)
├── ttci.py                             # Core TTCI routines (energy, pressure, TT CUR)

---

## Dependencies

You also need:
- LAMMPS (with Python bindings)
- ASE (Atomic Simulation Environment)
- Potential files (e.g. .include, .eam)
- Structure files in LAMMPS format (e.g. .data, .xyz → converted to .data)

---

## 📊 Reproduce Plots

### Cu with SMATB Energy and Pressure

    python Fig3_Cu.py

### Ar (HIPNN model) Energy and Pressure

    python Fig4_Ar_HIPNN.py

### Sn Phase Diagram (TTCI vs MD)

    python Fig5_phase_diagram.py

---

## Notes

- TTCI rank is controlled by the `--inds` argument (e.g., `-1 1` → rank 2)
- Quadrature rule controlled via `--method` (e.g., `Gauss3`, `Gauss8`)
- Ensure LAMMPS and ASE are correctly configured to find potentials and structure files
