::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::                                                                                                                            ::::
::::    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           ::::
::::                                                                                                                            ::::
::::    Copyright (C) 2012-present, The Computational Data Science Lab                                                          ::::
::::                                                                                                                            ::::
::::    This file is part of the ParaMonte library.                                                                             ::::
::::                                                                                                                            ::::
::::    LICENSE                                                                                                                 ::::
::::                                                                                                                            ::::
::::       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          ::::
::::                                                                                                                            ::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

@echo off
set ERRORLEVEL=0

@echo program getIntelVersion > getIntelVersion.F90
@echo use iso_fortran_env, only: compiler_version >> getIntelVersion.F90
@echo integer :: fileunit, startindex, endindex >> getIntelVersion.F90
@echo character(:), allocatable :: string >> getIntelVersion.F90
@echo string = trim(adjustl(compiler_version())) >> getIntelVersion.F90
@echo startindex = index(string,"Version") + 8 >> getIntelVersion.F90
@echo endindex = index(string,"Build") - 2 >> getIntelVersion.F90
@echo string = trim(adjustl(string(startindex:endindex))) >> getIntelVersion.F90
@echo open(newunit=fileunit,file="getIntelVersion.tmp") >> getIntelVersion.F90
@echo write(fileunit,"(A)") trim(adjustl(string)) >> getIntelVersion.F90
@echo end program getIntelVersion >> getIntelVersion.F90

ifort /nologo /standard-semantics /stand:f18 getIntelVersion.F90 -o getIntelVersion.exe | find "ERROR" >nul2>nul
if not ERRORLEVEL 1 (
    echo. 
    echo. -- ParaMonte - Fatal Error: Failed to compile the program to fetch the Intel Compiler suite version.
    echo. -- ParaMonte - Fatal Error: Please report this issue at https://github.com/cdslaborg/paramonte/issues
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)

timeout 1 >nul 2>nul

copy getIntelVersion.exe getIntelVersionCopy.exe >NUL
getIntelVersionCopy.exe | find "ERROR" >nul2>nul
if not ERRORLEVEL 1 (
    echo. -- ParaMonte - Fatal Error: Failed to run the program to fetch the Intel Compiler suite version.
    echo. -- ParaMonte - Fatal Error: Please report this issue at https://github.com/cdslaborg/paramonte/issues
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)

set /p COMPILER_VERSION=<getIntelVersion.tmp
REM echo. - COMPILER_VERSION=%COMPILER_VERSION%
del getIntelVersion.F90 getIntelVersion.obj getIntelVersion.exe getIntelVersionCopy.exe getIntelVersion.tmp