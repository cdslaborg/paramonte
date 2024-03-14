del main.exe
set PATH=..\..\..\lib;%PATH%
ifort /fpp /standard-semantics /Od /debug:full /CB /Qinit:snan,arrays /warn:all /gen-interfaces /traceback /check:all /fpe-all:0 /Qtrapuv /I:..\..\..\include main.F90 ..\..\..\lib\libparamonte_fortran_*_intel*.lib /exe:main.exe
mpiexec -localonly -n 3 main.exe
