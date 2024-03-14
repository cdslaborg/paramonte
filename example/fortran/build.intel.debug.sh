#!/usr/bin/env sh
rm main.exe
ifort -fpp -standard-semantics -O0 -g3 -CB -debug full -traceback -check all -fpe0 -Wl,-rpath,../../../lib -I../../../inc main.F90 ../../../lib/libparamonte_fortran_*_intel* -o main.exe
./main.exe