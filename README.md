# Tensors for High-dimensional Object Representation (THOR)

This repository contains a novel, multi-GPU implementation of tensors in tensor train format.  It is based on the modern Fortran MPI/GPU communication library (Thunder), and a CUDA-enabled algebra of distributed arrays (DRay).

The THOR Project (Tensors for High-dimensional Object Representation) aims to advance the state-of-the-art in tensor calculations, manipulation, and research. 
We strive to provide a high-performance tensor library for various scientific applications, containing ready-to-use utilities and applicaions in Fortran, Matlab, and Python. 

For the collection of published examples, navigate to [```TNI```](TNI) subfolder.

## How to Cite THOR

If you used THOR for scientific publications, please include the following citation:
```
@techreport{alexandrov2024thor,
  title        = {Tensors Optimized for High-level Research (THOR): an efficient and easy-to-use library for tensor networks},
  author       = {Alexandrov, Boian and Boureima, Ismael Djibrilla and Korobkin, Oleg and Danis, Mustafa Engin},
  institution  = {Los Alamos National Laboratory},
  number       = {LA-UR-24-24375},
  year         = {2024},
  month        = {may},
  day          = {6},
  note         = {Approved for public release; distribution is unlimited},
  type         = {Report}
}
```
Upcoming publications:
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
Developer test suites are located under [```TNI/fortran/tests/```](TNI/fortran/tests/) directory.

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
