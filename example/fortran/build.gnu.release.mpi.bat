del main.exe
set PATH=..\..\..\lib;%PATH%
gfortran -cpp -ffree-line-length-none -O3 -Wl,-rpath,../../../lib -I../../../inc main.F90 ../../../lib/libparamonte_fortran_*_gnu* -o main.exe
mpiexec -n 3 main.exe
