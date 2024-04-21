#!/usr/bin/env sh
rm main.exe
gfortran -cpp -ffree-line-length-none -O0 -g -fcheck=all -fbacktrace -Wl,-rpath,../../../lib -I../../../inc main.F90 ../../../lib/libparamonte* -o main.exe
./main.exe
