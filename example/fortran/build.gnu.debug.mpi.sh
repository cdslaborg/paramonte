#!/USR/BIN/ENV sh
rm main.exe
gfortran -cpp -ffree-line-length-none -O0 -g -fcheck=all -fbacktrace -Wl,-rpath,../../../lib -I../../../inc main.F90 ../../../lib/libparamonte* -o main.exe
mpiexec -n 3 ./main.exe
