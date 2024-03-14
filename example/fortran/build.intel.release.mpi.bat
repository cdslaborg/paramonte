del main.exe
set PATH=..\..\..\lib;%PATH%
ifort /fpp /standard-semantics /O3 /I:..\..\..\include main.F90 ..\..\..\lib\libparamonte_fortran_*_intel*.lib /exe:main.exe
mpiexec -localonly -n 3 main.exe
