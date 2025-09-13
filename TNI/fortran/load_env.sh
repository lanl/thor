#!/bin/bash

if [ $# -ne 1 ]; then
  echo  "Usage: source load_env.sh (gnu|mac|intel|ch-intel|gnu132)"
  return
fi

option=$1
case $option in
gnu) # Oleg's laptop
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
   module load openmpi/5.0.2-intel_2024.0.0 intel-mkl/2024.0.0
   ln -sf config/darwin_intel.inc make.inc
   module list
   ;;
ch-intel)
   echo "Loading Intel environment (on Chicoma):"
   module swap PrgEnv-cray PrgEnv-intel/8.5.0
   module load intel-mkl/2023.2.0
   ln -sf config/chicoma_intel.inc make.inc
   module list
   ;;
gnu132)
   echo "Loading GNU+OpenMPI environment (on Darwin):"
   module purge
   module load openmpi/5.0.2-gcc_13.2.0
   ln -sf config/darwin_gcc13_2.inc make.inc
   module list
   ;;
*)
   echo "Unsupported environment: $1" >> /dev/stderr
   echo "Use one of (gnu|intel|cuda)" >> /dev/stderr
   return
esac

ls -l make.inc
echo "Environment setup successfully"

