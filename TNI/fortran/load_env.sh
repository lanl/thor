#!/bin/bash

if [ $# -ne 1 ]; then
  echo  "Usage: source load_env.sh (gnu|mac|intel|cuda)"
  return
fi

option=$1
case $option in
gnu)
   echo "Loading GNU environment:"
   module purge
   module load openmpi 
   #module load gcc/8.2 openmpi/3.1.1-gcc_8.2 
   ln -sf config/gfortran.inc make.inc
   export LD_LIBRARY_PATH=$HOME/local/share/openblas/lib:$LD_LIBRARY_PATH
   module list
   ;;
mac)
   echo "Loading GNU environment:"
   ln -sf config/mac_gnu.inc make.inc
   ;;
intel)
   echo "Loading Intel environment (on Darwin):"
   module purge
   module load intel-mpi/2021.9.0  intel-mkl/2021.3.0
   ln -sf config/darwin_ifort.inc make.inc
   module list
   ;;
ch-intel)
   echo "Loading Intel environment (on Chicoma):"
   module swap PrgEnv-cray PrgEnv-intel/8.5.0
   module load intel-mkl/2023.2.0
   ln -sf config/chicoma_intel.inc make.inc
   module list
   ;;
intel2024)
   echo "Loading intel-2024 environment (on Darwin):"
   module purge
   module load openmpi/5.0.2-intel_2024.0.0
   ln -sf config/darwin_intel2024.inc make.inc
   module list
   ;;
gnu132)
   echo "Loading GNU+OpenMPI environment (on Darwin):"
   module purge
   module load openmpi/5.0.2-gcc_13.2.0
   ln -sf config/darwin_gcc13_2.inc make.inc
   module list
   ;;
aocc)
   echo "Loading Clang env (on Darwin)"
   module purge
   module load openmpi/4.1.0-aocc_3.0.0
   ln -sf config/darwin_aocc3.inc make.inc
   module list
   ;;
cuda)
   echo "Loading CUDA environment:"
   module purge
   module load nvhpc/23.3 openmpi
   ln -sf config/darwin_nvhpc.inc make.inc
   module list
   ;;
*)
   echo "Unsupported environment: $1" >> /dev/stderr
   echo "Use one of (gnu|intel|cuda)" >> /dev/stderr
   return
esac

ls -l make.inc
echo "Environment setup successfully"

