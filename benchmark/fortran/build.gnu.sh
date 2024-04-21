#!/bin/bash
gfortran -cpp -O3 -ffree-line-length-none -Wl,-rpath,../../../../lib -I../../../../inc -o main.exe main.F90 ../../../../lib/libparamonte*.so && ./main.exe