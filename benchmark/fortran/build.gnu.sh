#!/bin/bash
gfortran -cpp -O3 -ffree-line-length-none -Wl,-rpath,../../../../lib -I../../../../inc -o main.exe main.F90 ../../../../lib/libparamonte_fortran_*_gnu*.so && ./main.exe