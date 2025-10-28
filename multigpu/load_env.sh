#!/bin/bash

if [ $# -ne 1 ]; then
  echo  "Usage: . load_env.sh (gnu|intel|cuda)"
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
   ;;
cuda|cuda12)
   echo "Loading CUDA environment:"
   module purge
   module load nvhpc/25.1 cuda/12.6.3
   ln -sf config/darwin_cuda12.inc make.inc
   ;;
cuda11)
   echo "Loading CUDA environment:"
   module purge
   module load nvhpc/24.7 cuda/11.8.0
   ln -sf config/darwin_cuda11.inc make.inc
   ;;
*)
   echo "Unsupported environment: $1" >> /dev/stderr
   echo "Use one of (gnu|intel|cuda)" >> /dev/stderr
   return
esac

module list
ls -l make.inc
alias mrun="mpirun -x LD_LIBRARY_PATH"
echo "Alias  mrun='mpirun -x LD_LIBRARY_PATH' was created. Use it in place of mpirun"
echo "Environment setup successfully"
