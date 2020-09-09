::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::
::::   MIT License
::::
::::   ParaMonte: plain powerful parallel Monte Carlo library.
::::
::::   Copyright (C) 2012-present, The Computational Data Science Lab
::::
::::   This file is part of the ParaMonte library.
::::
::::   Permission is hereby granted, free of charge, to any person obtaining a 
::::   copy of this software and associated documentation files (the "Software"), 
::::   to deal in the Software without restriction, including without limitation 
::::   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
::::   and/or sell copies of the Software, and to permit persons to whom the 
::::   Software is furnished to do so, subject to the following conditions:
::::
::::   The above copyright notice and this permission notice shall be 
::::   included in all copies or substantial portions of the Software.
::::
::::   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
::::   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
::::   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
::::   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
::::   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
::::   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
::::   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
::::
::::   ACKNOWLEDGMENT
::::
::::   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
::::   As per the ParaMonte library license agreement terms, if you use any parts of 
::::   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
::::   work (education/research/industry/development/...) by citing the ParaMonte 
::::   library as described on this page:
::::
::::       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

@echo off
set ERRORLEVEL=0

@echo program getIntelVersion > getIntelVersion.f90
@echo use iso_fortran_env, only: compiler_version >> getIntelVersion.f90
@echo integer :: fileunit, startindex, endindex >> getIntelVersion.f90
@echo character(:), allocatable :: string >> getIntelVersion.f90
@echo string = trim(adjustl(compiler_version())) >> getIntelVersion.f90
@echo startindex = index(string,"Version") + 8 >> getIntelVersion.f90
@echo endindex = index(string,"Build") - 2 >> getIntelVersion.f90
@echo string = trim(adjustl(string(startindex:endindex))) >> getIntelVersion.f90
@echo open(newunit=fileunit,file="getIntelVersion.tmp") >> getIntelVersion.f90
@echo write(fileunit,"(A)") trim(adjustl(string)) >> getIntelVersion.f90
@echo end program getIntelVersion >> getIntelVersion.f90

ifort /nologo /standard-semantics /stand:f08 getIntelVersion.f90 -o getIntelVersion.exe | find "ERROR" >nul2>nul
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
del getIntelVersion.f90 getIntelVersion.obj getIntelVersion.exe getIntelVersionCopy.exe getIntelVersion.tmp