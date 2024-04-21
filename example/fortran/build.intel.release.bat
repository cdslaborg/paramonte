del main.exe
set PATH=..\..\..\lib;%PATH%
ifort /fpp /standard-semantics /O3 /I:..\..\..\include main.F90 ..\..\..\lib\libparamonte*.lib /exe:main.exe
main.exe
