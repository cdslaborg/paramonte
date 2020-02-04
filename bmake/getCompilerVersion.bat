::**********************************************************************************************************************************
::**********************************************************************************************************************************
::
::  ParaMonte: plain powerful parallel Monte Carlo library.
::
::  Copyright (C) 2012-present, The Computational Data Science Lab
::
::  This file is part of ParaMonte library. 
::
::  ParaMonte is free software: you can redistribute it and/or modify
::  it under the terms of the GNU Lesser General Public License as published by
::  the Free Software Foundation, version 3 of the License.
::
::  ParaMonte is distributed in the hope that it will be useful,
::  but WITHOUT ANY WARRANTY; without even the implied warranty of
::  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
::  GNU Lesser General Public License for more details.
::
::  You should have received a copy of the GNU Lesser General Public License
::  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
::
::**********************************************************************************************************************************
::**********************************************************************************************************************************

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