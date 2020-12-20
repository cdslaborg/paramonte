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
::::       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::
:: This is a simple standalone build Batch script for building C/Fortran applications that use ParaMonte library on Windows OS.
::
:: Usage - building C/C++ applications via Intel Visual C/C++ compiler:
::
::      build.bat
::
:: Usage - building C/C++ applications via Microsoft Visual C/C++ compiler:
::
::      build.bat msvc
::
:: Usage - building Fortran applications via Intel Visual Fortran compiler:
::
::      build.bat

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: build ParaMonteExample objects and executable
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

@echo off
set ERRORLEVEL=0
cd %~dp0

echo.

setlocal EnableDelayedExpansion

set BUILD_NAME=ParaMonte

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
        echo. -- !BUILD_NAME! - FATAL: ParaMonte library type could not be recognized.
        echo. -- !BUILD_NAME! - FATAL: Build failed. exiting...
        echo.
        cd %~dp0
        set ERRORLEVEL=1
        exit /B 1
    )
)

echo. -- !BUILD_NAME! - ParaMonte library name: !ParaMonte_LIB_NAME!
echo. -- !BUILD_NAME! - ParaMonte library type: !LTYPE!

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: determine library build
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set BTYPE=
echo !ParaMonte_LIB_NAME! | find "release" >nul
if errorlevel 1 (
    echo !ParaMonte_LIB_NAME! | find "testing" >nul
    if errorlevel 1 (
        echo !ParaMonte_LIB_NAME! | find "debug" >nul
        if errorlevel 1 (
            echo.
            for %%f in (find.exe) do @echo. -- !BUILD_NAME! - FATAL: Exhausted all builds in search of the ParaMonte library using command %%~dpfnx$PATH:f
            echo. -- !BUILD_NAME! - FATAL: Please make sure the above path command to the Windows find.exe application.
            echo. -- !BUILD_NAME! - FATAL: ParaMonte library build could not be recognized.
            echo. -- !BUILD_NAME! - FATAL: Build failed. exiting...
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

echo. -- !BUILD_NAME! - ParaMonte library build: !BTYPE!

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

echo. -- !BUILD_NAME! - ParaMonte library parallelism: !PTYPE!

if !PTYPE!==mpi (
    set MPICC_PATH=
    for /f tokens^=* %%a in ('where mpicc') do (
        set MPICC_PATH=%%~fa
        goto LABEL_MPI_FOUND
    )
    :LABEL_MPI_FOUND
    if defined MPICC_PATH (
        echo. -- !BUILD_NAME! - Intel C MPI library wrapper at: "!MPICC_PATH!"
    ) else (
        echo.
        echo. -- !BUILD_NAME! - FATAL: The build-script cannot find MPI library wrappers.
        echo. -- !BUILD_NAME! - FATAL: To build MPI-parallel ParaMonte examples, 
        echo. -- !BUILD_NAME! - FATAL: you need Intel's MPI library installed on your system.
        echo. -- !BUILD_NAME! - FATAL: you can download the latest Intel MPI library from their website.
        echo. -- !BUILD_NAME! - FATAL: or visit www.cdslab.org/pm for more instructions on how to build.
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
        echo !ParaMonte_LIB_NAME!|find "_cpp_" >nul
        if errorlevel 1 (
            echo.
            echo. -- !BUILD_NAME! - FATAL: ParaMonte library target lanugage could not be recognized.
            echo. -- !BUILD_NAME! - FATAL: Build failed. exiting...
            echo.
            cd %~dp0
            set ERRORLEVEL=1
            exit /B 1
        ) else (
            set TARGET_LANG_IS_CCPP=true
            set TARGET_LANG=C++
            set SEXT=cpp
        )
    ) else (
        set TARGET_LANG_IS_CCPP=true
        set TARGET_LANG=C
        set SEXT=c
    )
) else (
    set TARGET_LANG_IS_CCPP=false
    set TARGET_LANG=Fortran
    set SEXT=f90
)

echo. -- !BUILD_NAME! - ParaMonte library target language: !TARGET_LANG!

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: determine compiler/linker specs
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:: set default number of processes to be used in parallel mode

if not defined nproc set nproc=3
set SRC_FILES=logfunc.!SEXT! main.!SEXT!
set EXE_NAME=main.exe

set COMPILER_NAME=
if !TARGET_LANG!==Fortran (

    set COMPILER_NAME=ifort
    set COMPILER_SUITE=intel
    set LINKER_FLAGS=/link /out:!EXE_NAME!
    set SRC_FILES=paramonte.!SEXT! !SRC_FILES!

    if !PTYPE!==mpi (
        if !LTYPE!==static set COMPILER_NAME=mpiifort -fc=ifort
    )

    set COMPILER_FLAGS=/fpp
    REM set COMPILER_FLAGS=/fpp /DIS_COMPATIBLE_COMPILER
    if !BTYPE!==debug   set COMPILER_FLAGS=!COMPILER_FLAGS! /debug:full /CB /Od /Qinit:snan,arrays /warn:all /gen-interfaces /traceback /check:all /check:bounds /fpe-all:0 /Qdiag-error-limit:10 /Qdiag-disable:5268 /Qdiag-disable:7025 /Qtrapuv
    if !BTYPE!==release set COMPILER_FLAGS=!COMPILER_FLAGS! /O3 /Qip /Qipo /Qunroll /Qunroll-aggressive /Qinline-dllimport
    if !BTYPE!==testing set COMPILER_FLAGS=!COMPILER_FLAGS! /Od

    set CAF_ENABLED=false
    if !PTYPE!==cafsingle (
        set CAFTYPE=single
        set CAF_ENABLED=true
    )
    if !PTYPE!==cafshared (
        set CAFTYPE=shared
        set CAF_ENABLED=true
    )
    if !PTYPE!==cafdistributed (
        set CAFTYPE=distributed
        set CAF_ENABLED=true
    )
    if !CAF_ENABLED!==true (
        set FOR_COARRAY_NUM_IMAGES=!nproc!
        set "COARRAY_FLAGS=/Qcoarray:!CAFTYPE! /Qcoarray-num-images:!nproc!"
        set "LINKER_FLAGS=!COARRAY_FLAGS! !LINKER_FLAGS!"
        set "COMPILER_FLAGS=!COMPILER_FLAGS! !COARRAY_FLAGS!"
    )

)

if !TARGET_LANG_IS_CCPP!==true (

    set COMPILER_NAME=icl
    set COMPILER_SUITE=intel
    set LINKER_FLAGS=/link /out:!EXE_NAME!

    if !PTYPE!==mpi set COMPILER_NAME=mpiicc.bat -cc=!COMPILER_NAME!.exe

    set COMPILER_FLAGS=
    if !BTYPE!==debug   set COMPILER_FLAGS=/debug:full /Od /Wall /traceback /Qcheck-pointers:rw /Qcheck-pointers-undimensioned /Qdiag-error-limit:10 /Qtrapuv
    if !BTYPE!==release set COMPILER_FLAGS=/O3 /Qip /Qipo /Qunroll /Qunroll-aggressive /Ob2 /Qparallel /Qinline-dllimport
    if !BTYPE!==testing set COMPILER_FLAGS=/Od
    if !LTYPE!==dynamic set COMPILER_FLAGS=!COMPILER_FLAGS! /DDLL_ENABLED

    if "%1"=="msvc" (

        set COMPILER_NAME=cl
        set COMPILER_SUITE=microsoft visual C++
        REM set LINKER_FLAGS=/link /Fe:!EXE_NAME!
        if !PTYPE!==mpi set COMPILER_NAME=mpicl -cc=!COMPILER_NAME!

        set COMPILER_FLAGS=
        if !BTYPE!==debug   set COMPILER_FLAGS=/Od /Z7
        if !BTYPE!==release set COMPILER_FLAGS=/O2
        if !BTYPE!==testing set COMPILER_FLAGS=/Od
        if !LTYPE!==dynamic set COMPILER_FLAGS=!COMPILER_FLAGS! /DDLL_ENABLED

    )

)

echo. -- !BUILD_NAME! - ParaMonte example compiler suite: !COMPILER_SUITE!
echo. -- !BUILD_NAME! - ParaMonte example compiler name/wrapper: !COMPILER_NAME!
echo. -- !BUILD_NAME! - ParaMonte example compiler flags: !COMPILER_FLAGS!
echo. -- !BUILD_NAME! - ParaMonte example linker flags: !LINKER_FLAGS!

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: build ParaMonte example
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if not defined ParaMonteExample_EXE_ENABLED set ParaMonteExample_EXE_ENABLED=true

if !ParaMonteExample_EXE_ENABLED!==false goto LABEL_ParaMonteExample_RUN_ENABLED

echo. -- !BUILD_NAME! - example compile command:    "!COMPILER_NAME! !COMPILER_FLAGS! !SRC_FILES! !ParaMonte_LIB_NAME!.lib !LINKER_FLAGS!"
if !PTYPE!==mpi (
                                                call !COMPILER_NAME! !COMPILER_FLAGS! !SRC_FILES! !ParaMonte_LIB_NAME!.lib !LINKER_FLAGS! || (
                                                     echo. 
                                                     echo. -- !BUILD_NAME! - FATAL: Failed to compile and link the MPI-parallel application. 
                                                     echo. -- !BUILD_NAME! - FATAL: If you are using the Microsoft Visual Studio C/C++ Compiler to 
                                                     echo. -- !BUILD_NAME! - FATAL: compile a C/C++ application try:
                                                     echo. -- !BUILD_NAME! - FATAL: 
                                                     echo. -- !BUILD_NAME! - FATAL:     build.bat msvc
                                                     echo. -- !BUILD_NAME! - FATAL: 
                                                     echo. -- !BUILD_NAME! - FATAL: Exiting the ParaMonte example build script...
                                                     echo. 
                                                     cd %~dp0
                                                     set ERRORLEVEL=1
                                                     exit /B 1
                                                )
) else (
                                                     !COMPILER_NAME! !COMPILER_FLAGS! !SRC_FILES! !ParaMonte_LIB_NAME!.lib !LINKER_FLAGS! || (
                                                     echo. 
                                                     echo. -- !BUILD_NAME! - FATAL: Failed to compile and link the application. 
                                                     echo. -- !BUILD_NAME! - FATAL: If you are using the Microsoft Visual Studio C/C++ Compiler to 
                                                     echo. -- !BUILD_NAME! - FATAL: compile a C/C++ application try:
                                                     echo. -- !BUILD_NAME! - FATAL: 
                                                     echo. -- !BUILD_NAME! - FATAL:     build.bat msvc
                                                     echo. -- !BUILD_NAME! - FATAL: 
                                                     echo. -- !BUILD_NAME! - FATAL: Exiting the ParaMonte example build script...
                                                     echo. 
                                                     cd %~dp0
                                                     set ERRORLEVEL=1
                                                     exit /B 1
                                                )
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: enlarge stack size for ParaMonte example
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:LABEL_ParaMonteExample_RUN_ENABLED

if not defined ParaMonteExample_RUN_ENABLED set ParaMonteExample_RUN_ENABLED=true

if !ParaMonteExample_RUN_ENABLED!==false goto LABEL_EOF

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
    echo. -- !BUILD_NAME! - running MPI-parallelized ParaMonte example on !nproc! processes...
    echo. 
    mpiexec -localonly -n !nproc! !EXE_NAME! || (
        echo. 
        echo. -- !BUILD_NAME! - FATAL: Failed to run the MPI-parallel application. exiting...
        echo. 
        cd %~dp0
        set ERRORLEVEL=1
        exit /B 1
    )
) else (
    echo. 
    echo. -- !BUILD_NAME! - running serial ParaMonte example on 1 process...
    echo. 
    !EXE_NAME! || (
        echo. 
        echo. -- !BUILD_NAME! - FATAL: Failed to run the serial application. exiting...
        echo. 
        cd %~dp0
        set ERRORLEVEL=1
        exit /B 1
    )
)

goto LABEL_EOF

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: return
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:LABEL_EOF

if !TARGET_LANG!==Fortran set TARGET_LANG_WEB=fortran
if !TARGET_LANG!==C++ set TARGET_LANG_WEB=cpp
if !TARGET_LANG!==C set TARGET_LANG_WEB=c

echo. -- !BUILD_NAME! - NOTE: For information on building ParaMonte simulations in !TARGET_LANG!, visit:
echo. -- !BUILD_NAME! - NOTE: 
echo. -- !BUILD_NAME! - NOTE:     https://www.cdslab.org/paramonte/notes/build/!TARGET_LANG_WEB!
echo. -- !BUILD_NAME! - NOTE: 
echo. -- !BUILD_NAME! - NOTE: For information on running ParaMonte simulations in !TARGET_LANG!, visit:
echo. -- !BUILD_NAME! - NOTE: 
echo. -- !BUILD_NAME! - NOTE:     https://www.cdslab.org/paramonte/notes/run/default
echo.

cd %~dp0

exit /B 0
