#!/USR/BIN/ENV sh
rm main.exe
ifort -fpp -standard-semantics -O0 -g3 -CB -debug full -traceback -check all -fpe0 -Wl,-rpath,../../../lib -I../../../inc main.F90 ../../../lib/libparamonte* -o main.exe
mpiexec -n 3 ./main.exe