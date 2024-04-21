#!/bin/bash
ifort -fpp -standard-semantics -O3 -Wl,-rpath,../../../../lib -I../../../../inc -o main.exe main.F90 ../../../../lib/libparamonte*.so && ./main.exe || exit 1