#!/usr/bin/env sh
rm main.exe
ifort -fpp -standard-semantics -O3 -Wl,-rpath,../../../lib -I../../../inc main.F90 ../../../lib/libparamonte_fortran_*_intel* -o main.exe
./main.exe