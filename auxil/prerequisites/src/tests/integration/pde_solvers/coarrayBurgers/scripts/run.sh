#!/bin/sh
echo "Building the code:"
export TAU_MAKEFILE=/opt/paratools/tau/x86_64/lib/Makefile.tau-mpi-pdt
export TAU_OPTIONS="-optVerbose -optCompInst"
make clean
make -f Makefile.inst

# Specify TAU parameters here:
export TAU_CALLPATH=1
export TAU_CALLPATH_DEPTH=100
#export TAU_SAMPLING=1

echo "Running the code:"
mpiexec -np 4 ./burgers

echo "Running the pprof command:"
pprof
echo "Running the TAU paraprof analyzer command:"
paraprof &
