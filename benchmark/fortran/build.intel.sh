#!/bin/bash
ifort -fpp -standard-semantics -O3 -Wl,-rpath,../../../../lib -I../../../../inc -o main.exe main.F90 ../../../../lib/libparamonte_fortran_*_intel*.so && ./main.exe || exit 1