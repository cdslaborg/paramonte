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

:: This is a simple standalone build Batch script for building C/Fortran applications that use ParaMonte library on Windows OS.
::
:: Usage:
::
::      build.bat COMPILER_SUITE PARALLELISM PARAMONTE_LIBRARY_TYPE BUILD_TYPE
::
:: where the passed arguments can be one of the following options:
::
::      COMPILER_SUITE: intel
::      PARALLELISM: serial / mpi
::      PARAMONTE_LIBRARY_TYPE: static / dynamic
::      BUILD_TYPE: release / testing / debug
::
:: All of the above choices except the first depend on the way the ParaMonte library has been built in your case.  
::
:: Example Usage:
::
::      build.bat intel serial static release

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: build ParaMonteExample objects and executable
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

@echo off
set ERRORLEVEL=0
cd %~dp0

REM if defined "%2" echo %2
REM if defined "%3" echo %3
REM if defined "%4" echo %4

REM set COMPILER_SUITE="%1"
REM set PARALLELISM="%2"
REM set LTYPE="%3"
REM set BTYPE="%4"

setlocal EnableDelayedExpansion

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: get ParaMonte library filename
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set ParaMonte_LIB_NAME=
for /f tokens^=* %%a in ('where .:libparamonte_*.dll') do (
    set ParaMonte_LIB_NAME=%%~na
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: determine library type
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if defined ParaMonte_LIB_NAME (
    set LTYPE=dynamic
) else (
    for /f tokens^=* %%a in ('where .:libparamonte_*.lib') do (
        set ParaMonte_LIB_NAME=%%~na
    )
    if defined ParaMonte_LIB_NAME (
        set LTYPE=static
    ) else (
        echo.
        echo. -- ParaMonte - Fatal Error: ParaMonte library type could not be recognized.
        echo. -- ParaMonte - build failed. exiting...
        echo.
        cd %~dp0
        set ERRORLEVEL=1
        exit /B 1
    )
)

echo. -- ParaMonte - ParaMonte library name: !ParaMonte_LIB_NAME!
echo. -- ParaMonte - ParaMonte library type: !LTYPE!

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: determine library build
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set BTYPE=
echo !ParaMonte_LIB_NAME!|find "release" >nul
if errorlevel 1 (
    echo !ParaMonte_LIB_NAME!|find "testing" >nul
    if errorlevel 1 (
        echo !ParaMonte_LIB_NAME!|find "debug" >nul
        if errorlevel 1 (
            echo.
            echo. -- ParaMonte - Fatal Error: ParaMonte library build could not be recognized.
            echo. -- ParaMonte - build failed. exiting...
            echo.
            cd %~dp0
            set ERRORLEVEL=1
            exit /B 1
        ) else (
            set BTYPE=debug
        )
    ) else (
        set BTYPE=testing
    )
) else (
    set BTYPE=release
)

echo. -- ParaMonte - ParaMonte library build: !BTYPE!

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: determine library parallelism
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set PTYPE=
echo !ParaMonte_LIB_NAME!|find "mpi" >nul
if errorlevel 1 (
    echo !ParaMonte_LIB_NAME!|find "cafshared" >nul
    if errorlevel 1 (
        echo !ParaMonte_LIB_NAME!|find "cafsingle" >nul
        if errorlevel 1 (
            set PTYPE=serial
        ) else (
            set PTYPE=cafsingle
        )
    ) else (
        set PTYPE=cafshared
    )
) else (
    set PTYPE=mpi
)

echo. -- ParaMonte - ParaMonte library parallelism: !PTYPE!

if !PTYPE!==mpi (
    set MPICC_PATH=
    for /f tokens^=* %%a in ('where mpicc') do (
        set MPICC_PATH=%%~fa
        goto LABEL_MPI_FOUND
    )
    :LABEL_MPI_FOUND
    if defined MPICC_PATH (
        echo. -- ParaMonte - Intel C MPI library wrapper at: "!MPICC_PATH!"
    ) else (
        echo.
        echo. -- ParaMonte - Fatal Error: The build script cannot find MPI library wrappers.
        echo. -- ParaMonte - Fatal Error: To build MPI-parallel ParaMonte examples, 
        echo. -- ParaMonte - Fatal Error: you need Intel's MPI library installed on your system.
        echo. -- ParaMonte - Fatal Error: you can download the latest Intel MPI library from their website.
        echo. -- ParaMonte - Fatal Error: or visit www.cdslab.org/pm for more instructions on how to build.
        echo.
        cd %~dp0
        set ERRORLEVEL=1
        exit /B 1
    )
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: determine library target language
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set TARGET_LANG=
echo !ParaMonte_LIB_NAME!|find "_fortran_" >nul
if errorlevel 1 (
    echo !ParaMonte_LIB_NAME!|find "_c_" >nul
    if errorlevel 1 (
        echo.
        echo. -- ParaMonte - Fatal Error: ParaMonte library target lanugage could not be recognized.
        echo. -- ParaMonte - build failed. exiting...
        echo.
        cd %~dp0
        set ERRORLEVEL=1
        exit /B 1
    ) else (
        set TARGET_LANG=C
    )
) else (
    set TARGET_LANG=Fortran
)

echo. -- ParaMonte - ParaMonte library target language: !TARGET_LANG!

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: determine compiler/linker specs
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:: set default number of processes to be used in parallel mode

set NUM_PROCESS=3
set EXE_NAME=runExample.exe
set LINKER_FLAGS=/link /out:!EXE_NAME!

set COMPILER_NAME=
if !TARGET_LANG!==Fortran (

    set SEXT=f90
    set COMPILER_NAME=ifort
    set COMPILER_SUITE=intel

    if !PTYPE!==mpi (
        set COMPILER_NAME=mpiifort -fc=ifort
    )

    if !PTYPE!==cafshared (
        set FOR_COARRAY_NUM_IMAGES=!NUM_PROCESS!
    )

    set COMPILER_FLAGS=
    if !BTYPE!==debug   set COMPILER_FLAGS=/debug:full /CB /Od /Qinit:snan,arrays /warn:all /gen-interfaces /traceback /check:all /check:bounds /fpe-all:0 /Qdiag-error-limit:10 /Qdiag-disable:5268 /Qdiag-disable:7025 /Qtrapuv
    if !BTYPE!==release set COMPILER_FLAGS=/O3 /Qip /Qipo /Qunroll /Qunroll-aggressive /Qinline-dllimport
    if !BTYPE!==testing set COMPILER_FLAGS=/Od
    if !LTYPE!==dynamic set COMPILER_FLAGS=!COMPILER_FLAGS! /fpp /DIS_COMPATIBLE_COMPILER

)

if !TARGET_LANG!==C (

    set SEXT=c
    set COMPILER_NAME=icl
    set COMPILER_SUITE=intel

    if !PTYPE!==mpi set COMPILER_NAME=mpiicc.bat -cc=!COMPILER_NAME!.exe

    set COMPILER_FLAGS=
    if !BTYPE!==debug   set COMPILER_FLAGS=/debug:full /Od /Wall /traceback /Qcheck-pointers:rw /Qcheck-pointers-undimensioned /Qdiag-error-limit:10 /Qtrapuv
    if !BTYPE!==release set COMPILER_FLAGS=/O3 /Qip /Qipo /Qunroll /Qunroll-aggressive /Ob2 /Qparallel /Qinline-dllimport
    if !BTYPE!==testing set COMPILER_FLAGS=/Od
    if !LTYPE!==dynamic set COMPILER_FLAGS=!COMPILER_FLAGS! /DDLL_ENABLED

    if "%1"=="msvc" (

        set COMPILER_NAME=cl
        set COMPILER_SUITE=microsoft visual C++
        if !PTYPE!==mpi set COMPILER_NAME=mpicl -cc=!COMPILER_NAME!

        set COMPILER_FLAGS=
        if !BTYPE!==debug   set COMPILER_FLAGS=/Od /Z7
        if !BTYPE!==release set COMPILER_FLAGS=/O2
        if !BTYPE!==testing set COMPILER_FLAGS=/Od
        if !LTYPE!==dynamic set COMPILER_FLAGS=!COMPILER_FLAGS! /DDLL_ENABLED

    )

)

echo. -- ParaMonte - ParaMonte example compiler suite: !COMPILER_SUITE!
echo. -- ParaMonte - ParaMonte example compiler name/wrapper: !COMPILER_NAME!
echo. -- ParaMonte - ParaMonte example compiler flags: !COMPILER_FLAGS!
echo. -- ParaMonte - ParaMonte example linker flags: !LINKER_FLAGS!

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: build ParaMonte example
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if not defined ParaMonteExample_EXE_ENABLED set ParaMonteExample_EXE_ENABLED=true

if !ParaMonteExample_EXE_ENABLED! NEQ true goto LABEL_ParaMonteExample_RUN_ENABLED

if !PTYPE!==mpi (
    echo. -- ParaMonte - ParaMonte example compiling: "!COMPILER_NAME! !COMPILER_FLAGS! logfunc.!SEXT! main.!SEXT! !ParaMonte_LIB_NAME!.lib !LINKER_FLAGS!"
    call !COMPILER_NAME! !COMPILER_FLAGS! logfunc.!SEXT! main.!SEXT! !ParaMonte_LIB_NAME!.lib !LINKER_FLAGS!
) else (
    !COMPILER_NAME! !COMPILER_FLAGS! logfunc.!SEXT! main.!SEXT! !ParaMonte_LIB_NAME!.lib !LINKER_FLAGS!
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: enlarge stack size for ParaMonte example
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:LABEL_ParaMonteExample_RUN_ENABLED

if not defined ParaMonteExample_RUN_ENABLED set ParaMonteExample_RUN_ENABLED=true

if !ParaMonteExample_RUN_ENABLED! NEQ true goto LABEL_EOF

if !LTYPE!==static (
    editbin /STACK:99999999 !EXE_NAME!
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: Let Windows digest the info
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

TIMEOUT 1 >nul 2>nul

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: run ParaMonte example
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

REM https://software.intel.com/en-us/mpi-developer-reference-windows-compilation-commands

if !PTYPE!==mpi (
    echo. 
    echo. -- ParaMonte - running MPI-parallelized ParaMonte example on !NUM_PROCESS! processes...
    echo. 
    mpiexec -localonly -n !NUM_PROCESS! !EXE_NAME!
) else (
    echo. 
    echo. -- ParaMonte - running serial ParaMonte example on 1 process...
    echo. 
    !EXE_NAME!
)

goto LABEL_EOF

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: return
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:LABEL_EOF

cd %~dp0

exit /B 0
