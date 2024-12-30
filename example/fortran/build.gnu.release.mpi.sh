#!/usr/bin/env sh
rm main.exe
gfortran -cpp -ffree-line-length-none -O3 -Wl,-rpath,../../../lib -I../../../inc main.F90 ../../../lib/libparamonte* -o main.exe
mpiexec -n 3 ./main.exe