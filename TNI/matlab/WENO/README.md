## Information

    Developer   : Mustafa Engin Danis (@LANL, T-1)
    Date        : 10/26/2023
    Last Update : 04/28/2025
    Purpose     : Full-tensor Finite Difference Weighted Essentially Non-oscillatory (WENO) Solver for Hyperbolic PDEs
    Publication : Danis, M. Engin et al., Tensor-Train WENO Scheme for Compressible Flows, Journal of Computational Physics 529.C (2025), https://doi.org/10.1016/j.jcp.2025.113891

    Please cite the above paper if you use any part of this WENO library

## Contents

`src/ft`   : Full-tensor solver
`src/tt`   : Tensor-train solver
`src/misc` : Other helper functions 

## Usage:
  - In startup.m, do: `setenv("TTWENOLIB","/path/of/WENO")`
  - In the main MATLAB script to run simulations, 
    - for TT solver, do: `run(getenv("TTWENOLIB")+"/tt_setup.m")`
    - for FT solver, do: `run(getenv("TTWENOLIB")+"/ft_setup.m")`
  -  running in a terminal might require exporting TTWENOLIB as an environment variable: `export TTWENOLIB=/path/of/WENO`
