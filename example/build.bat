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
        echo. -- ParaMonte - FATAL: ParaMonte library type could not be recognized.
        echo. -- ParaMonte - FATAL: Build failed. exiting...
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
            echo. -- ParaMonte - FATAL: ParaMonte library build could not be recognized.
            echo. -- ParaMonte - FATAL: Build failed. exiting...
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
        echo. -- ParaMonte - FATAL: The build-script cannot find MPI library wrappers.
        echo. -- ParaMonte - FATAL: To build MPI-parallel ParaMonte examples, 
        echo. -- ParaMonte - FATAL: you need Intel's MPI library installed on your system.
        echo. -- ParaMonte - FATAL: you can download the latest Intel MPI library from their website.
        echo. -- ParaMonte - FATAL: or visit www.cdslab.org/pm for more instructions on how to build.
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
            echo. -- ParaMonte - FATAL: ParaMonte library target lanugage could not be recognized.
            echo. -- ParaMonte - FATAL: Build failed. exiting...
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

echo. -- ParaMonte - ParaMonte library target language: !TARGET_LANG!

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
        set COMPILER_NAME=mpiifort -fc=ifort
    )

    set COMPILER_FLAGS=/fpp /DIS_COMPATIBLE_COMPILER
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

echo. -- ParaMonte - ParaMonte example compiler suite: !COMPILER_SUITE!
echo. -- ParaMonte - ParaMonte example compiler name/wrapper: !COMPILER_NAME!
echo. -- ParaMonte - ParaMonte example compiler flags: !COMPILER_FLAGS!
echo. -- ParaMonte - ParaMonte example linker flags: !LINKER_FLAGS!

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: build ParaMonte example
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if not defined ParaMonteExample_EXE_ENABLED set ParaMonteExample_EXE_ENABLED=true

if !ParaMonteExample_EXE_ENABLED!==false goto LABEL_ParaMonteExample_RUN_ENABLED

if !PTYPE!==mpi (
    echo. -- ParaMonte - ParaMonte example compiling: "!COMPILER_NAME! !COMPILER_FLAGS! !SRC_FILES! !ParaMonte_LIB_NAME!.lib !LINKER_FLAGS!"
                                                  call !COMPILER_NAME! !COMPILER_FLAGS! !SRC_FILES! !ParaMonte_LIB_NAME!.lib !LINKER_FLAGS! || (
                                                        echo. 
                                                        echo. -- !BUILD_SCRIPT_NAME! - FATAL: Failed to compile and link the MPI-parallel application. exiting...
                                                        echo. -- !BUILD_SCRIPT_NAME! - FATAL: I you are using the Microsoft Visual Studio C/C++ Compiler to 
                                                        echo. -- !BUILD_SCRIPT_NAME! - FATAL: compile a C/C++ application try:
                                                        echo. -- !BUILD_SCRIPT_NAME! - FATAL: 
                                                        echo. -- !BUILD_SCRIPT_NAME! - FATAL:     build.bat msvc
                                                        echo. -- !BUILD_SCRIPT_NAME! - FATAL: 
                                                        echo. -- !BUILD_SCRIPT_NAME! - FATAL: Exiting the ParaMonte example build script...
                                                        echo. 
                                                        cd %~dp0
                                                        set ERRORLEVEL=1
                                                        exit /B 1
                                                  )
) else (
                                                       !COMPILER_NAME! !COMPILER_FLAGS! !SRC_FILES! !ParaMonte_LIB_NAME!.lib !LINKER_FLAGS! || (
                                                        echo. 
                                                        echo. -- !BUILD_SCRIPT_NAME! - FATAL: Failed to compile and link the application. exiting...
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
    echo. -- ParaMonte - running MPI-parallelized ParaMonte example on !nproc! processes...
    echo. 
    mpiexec -localonly -n !nproc! !EXE_NAME! || (
        echo. 
        echo. -- !BUILD_SCRIPT_NAME! - FATAL: failed to run the MPI-parallel application. exiting...
        echo. 
        cd %~dp0
        set ERRORLEVEL=1
        exit /B 1
    )
) else (
    echo. 
    echo. -- ParaMonte - running serial ParaMonte example on 1 process...
    echo. 
    !EXE_NAME! || (
        echo. 
        echo. -- !BUILD_SCRIPT_NAME! - FATAL: failed to run the serial application. exiting...
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

cd %~dp0

exit /B 0
