#!/bin/bash
source /usr/local/packages/Modules/init/bash
module load intel/14.0
source /usr/local/packages/intel/14.0/bin/ifortvars.sh intel64
export TAU_MAKEFILE=/usr/local/packages/tau-2.22.2/x86_64/lib/Makefile.tau-coarray-icpc-mpi
export TAU_COMM_MATRIX=1
export PATH=/usr/local/packages/tau-2.22.2/x86_64/bin:$PATH
module load java
export PATH=/usr/local/packages/intel/14.0/mpirt/bin/intel64/:$PATH
#export FOR_COARRAY_NUM_IMAGES=1
echo "Using FOR_COARRAY_NUM_IMAGES = " $FOR_COARRAY_NUM_IMAGES
mpiexec --mca btl_tcp_if_include eth2 -np 16 ./burgers_caf
#./burgers_caf
