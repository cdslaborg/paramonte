del main.exe
set PATH=..\..\..\lib;%PATH%
gfortran -cpp -ffree-line-length-none -O0 -g -fcheck=all -fbacktrace -Wl,-rpath,../../../lib -I../../../inc main.F90 ../../../lib/libparamonte_fortran_*_gnu* -o main.exe
main.exe
