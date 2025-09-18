# Tensor Network Integrators

This folder contains applications in Matlab and Fortran, for a variety of scientific applications.

## Capabilities

The following capabilities are provided in [TNI](TNI):
- [```HDI```](examples/HDI): multidimensional numerical integration using a low-rank TT ([Alexandrov et al. 2023](https://doi.org/10.3390/math11030534));
- [```CI```](../TT_Configurational_Integral): Solving configurational integrals for crystalline solids ([Truong et al. 2025](https://arxiv.org/abs/2505.21826));
- [```Linear-STSC```](examples/Linear-STSC): TT space-time spectral collocation method for convection-diffusion-reaction
           (CDR) equation ([Adak et al. 2024](https://arxiv.org/abs/2402.18073));
           with variable coefficients ([Adak et al. 2025](https://www.mdpi.com/2227-7390/13/14/2277));
- [```Non-linear-STSC```](examples/Non-linear-STSC): space-time spectral collocation method for the nonlinear convection diffusion equation ([Adak et al. 2025](https://arxiv.org/abs/2406.02505));
- [```Maxwell-Mimetic```](examples/Maxwell-Mimetic): mimetic finite difference method for 3D wave propagation ([Manzini et al. 2023](https://doi.org/10.1016/j.matcom.2023.03.026));
- [```NTE```](examples/NTE): time-independent Boltzmann neutron transport equation ([Truong et al. 2024](https://www.sciencedirect.com/science/article/pii/S002199912400192X));
- [```RK-1```](examples/RK-1): explicit and implicit Runge-Kutta integrators (Chinomona et al. 2025, LA-UR-25-25351);
- [```SWE```](examples/SWE): high-order finite volume methods for shallow water equations with TT ([Danis et al. 2025](https://doi.org/10.1175/MWR-D-24-0165.1));
- [```TENO```](examples/TENO): TT TENO scheme for compressible flows ([Danis et al. 2025](https://doi.org/10.2514/6.2025-0304));
- [```WENO```](examples/WENO): TT WENO scheme for compressible flows ([Danis et al. 2025](https://doi.org/10.1016/j.jcp.2025.113891));
- [```OpInf```](examples/OpInf): TT operator inference ([Danis et al._2025](https://doi.org/10.48550/arXiv.2509.08071));
- [```IGA```](examples/IGA): isogeometric framework for complex geometries ([Tran et al._2025](https://doi.org/10.48550/arXiv.2509.13224));
- [```Multimat```](examples/Multimat): Multimaterial hydrodynamics (Truong et al., in prep.).

## Installation

When cloning the git repo, make sure to pull it with all of its submodules:
```sh
  git clone --recursive https://github.com/lanl/thor.git
```
Alternatively, you can pull the submodules by running the following git command:
```sh
  git submodule update --init --recursive
```

### Running MATLAB inside Jupyter Notebooks

The [Jupyter Notebook](https://jupyter.org/) is a server-client application, very popular in Python community, that allows editing and running interactive notebooks with notes interspersed with the code, via a web browser. Jupyter can be installed on a wide variety of platforms, and run not only with Python but also with other programming languages, including MATLAB.  

All of the examples in the [```examples```](examples) subdirectory are documented in Jupyter notebooks. 
To have Jupyter on your local system with an interactive MATLAB code, please follows these steps:

1. Have a working MATLAB with `matlab` command in your path.
2. Use [this link](https://jupyter.org/install) for the up-to-date instructions on how to install Jupyter Notebook on your system.
2. Install [Jupyter-MATLAB proxy](https://github.com/mathworks/jupyter-matlab-proxy).

### [FORTRAN](fortran), serial version

* Clone the repo and change to [`fortran`](fortran) directory.
* In the subdirectory `config`, find a configuration that most closely matches your system.
* Copy it to [`fortran`](fortran) under the name `make.inc`. This is a system-dependent part of the Makefile.
* Modify it to match your compiler and the location of BLAS/LAPACK libraries.
* Build the code using `make`. Test using `make test`.

```shell
git clone git@github.com:lanl/thor.git
cd thor/TNI/fortran
cp config/gfortran.inc make.inc
# <edit make.inc>
make -j
make test # optional
```

You can also follow the (# LANL HPC Installation Notes) for some of the existing LANL systems (current as of September 2025).
The code in this folder is derived from the following third-party libraries and packages:
 - TT-Toolbox : https://github.com/oseledets/TT-Toolbox
 - tt-fort    : https://github.com/oseledets/tt-fort
 - ttcross    : https://github.com/savostyanov/ttcross

### LANL HPC Installation Notes

#### Darwin
```shell
salloc -p general
cd fortran
. load_env.sh intel
make -j
make test # optional
```

#### Chicoma
```shell
salloc --qos=debug --reservation=debug --partition=debug
cd fortran
. load_env.sh chicoma
make -j
make test # optional
```

