del main.exe
set PATH=..\..\..\lib;%PATH%
gfortran -cpp -ffree-line-length-none -O3 -Wl,-rpath,../../../lib -I../../../inc main.F90 ../../../lib/libparamonte* -o main.exe
main.exe
