# Tensors for High-dimensional Object Representation(THOR) / [Wiki](https://TODO)

This repository contains a novel, multi-GPU implementation of tensors in tensor train format.  It is based on the modern Fortran MPI/GPU communication library (Thunder), and a CUDA-enabled algebra of distributed arrays (DRay).

The THOR Project (Tensors for High-dimensional Object Representation) aims to advance the state-of-the-art in tensor calculations, manipulation, and research. 
We strive to provide a high-performance tensor library for various scientific applications, containing ready-to-use utilities and applicaions in Fortran, Matlab, and Python. 

Check our [wiki page](https://TODO) for introduction, getting started, etc.


## Installation

### MATLAB

Please download [TT-toolbox](https://github.com/oseledets/TT-Toolbox) package at the same directory level as THOR. Refer to individual application in for additional tools.

### FORTRAN (serial version in TNI/fortran)

* Clone the repo and change to `TNI/fortran` directory.
* In the subdirectory `config`, find a configuration that most closely matches your system.
* Copy it to `TNI/fortran` under the name `make.inc`. This is a system-dependent part of the Makefile.
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

## Capabilities

List capabilities here

## Examples

Examples go here

## How to Cite THOR

Publications:
 * Boureima et al. (in prep)
 * Truong et al. (in prep)

## Authors
- [Boian S. Alexandrov](mailto:boian@lanl.gov): Theoretical Division, Los Alamos National Laboratory
- [Ismael Boureima](mailto:iboureima@lanl.gov): Theoretical Division, Los Alamos National Laboratory
- [Rujeko Chinomona](mailto:crujeko@lanl.gov): Theoretical Division, Los Alamos National Laboratory
- [William Dai](mailto:dai@lanl.gov): Theoretical Division, Los Alamos National Laboratory
- [Engin Danis](mailto:danis@lanl.gov): Theoretical Division, Los Alamos National Laboratory
- [Maksim Ekin Eren](mailto:maksim@lanl.gov): Information Systems and Modeling Group, Los Alamos National Laboratory ([Website](https://www.maksimeren.com/))
- [Oleg Korobkin](mailto:korobkin@lanl.gov): Theoretical Division, Los Alamos National Laboratory
- [Rahul Somasundaram](mailto:rahul@lanl.gov): Theoretical Division, Los Alamos National Laboratory
- [Kim Rasmussen](mailto:kor@lanl.gov): Theoretical Division, Los Alamos National Laboratory
- [Duc Truong](mailto:dptruong@lanl.gov): Theoretical Division, Los Alamos National Laboratory

## Maintainers
- [Ismael Boureima](mailto:iboureima@lanl.gov): Theoretical Division, Los Alamos National Laboratory
- [Oleg Korobkin](mailto:korobkin@lanl.gov): Theoretical Division, Los Alamos National Laboratory
- [Rahul Somasundaram](mailto:rahul@lanl.gov): Theoretical Division, Los Alamos National Laboratory
- [Duc Truong](mailto:dptruong@lanl.gov): Theoretical Division, Los Alamos National Laboratory

## Copyright Notice
>© 2025. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit
others to do so.

**LANL Copyright Assertion #O4849 **

## License
This program is open source under the BSD-3 License.
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and
the following disclaimer.
 
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
and the following disclaimer in the documentation and/or other materials provided with the
distribution.
 
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


## Developer Test Suite
Developer test suites are located under [```tests/```](tests/) directory.

## LANL HPC Installation Notes

### Darwin
```shell
salloc -p general
cd fortran
. load_env.sh intel
make -j
make test # optional
```

### Chicoma
```shell
salloc --qos=debug --reservation=debug --partition=debug
cd fortran
. load_env.sh chicoma
make -j
make test # optional
```
