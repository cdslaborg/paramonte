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

:: NOTE: Do not change the contents of this file unless you know what the consequences are.
:: This is the batch file that builds objects, dynamic libraries, as well as the test and example binaries of the ParaMonte library on Windows systems.
:: Upon invocation of this file from Intel Parallel Studio's Windows command line (Intel's cmd.exe, which is provided via Microsoft Visual Studio,
:: and is preloaded with Intel suite environmental flags), this file will first call the configuration file configParaMonte.bat to read the user's
:: requested configuration for building the ParaMonte library.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set build type: release, debug, testing :: set library type: static, dynamic
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:: setlocal EnableDelayedExpansion

:: silence cmd output
@echo off

:: setlocal EnableDelayedExpansion

:: define variables only locally
::setlocal

set ParaMonte_BLD_ROOT_DIR=%cd%

:: change directory to the folder containing this batch file 
cd %~dp0

:: fetch ParaMonte library kernel version

set "ParaMonteVersion="

cd .\bmake\
for /f "tokens=*" %%i in ('head.bat 1 "..\.VERSION"') do set "ParaMonteVersion=%%i"
cd %~dp0

set "FPP_PARAMONTE_VERSION_FLAG="
if defined ParaMonteVersion (
    set FPP_PARAMONTE_VERSION_FLAG=/define:PARAMONTE_VERSION='!ParaMonteVersion!'
)

:: fetch ParaMonte library kernel release date

set dayName=!date:~0,3!
set year=!date:~10,4!
set day=!date:~7,2!

set m=100
for %%m in (Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec) do (
    set /a m+=1
    set month[!m:~-2!]=%%m
)
set monthNow=%date:~3,3%
set monthNow=%monthNow: =%
set monthName=!month[%monthNow%]!
set ParaMonteRelease=!dayName!.!monthName!.!day!.!year!

set "FPP_PARAMONTE_RELEASE_FLAG="
set FPP_PARAMONTE_RELEASE_FLAG=/define:PARAMONTE_RELEASE='!ParaMonteRelease!'

set SERIAL_ENABLED=true

REM echo. 
REM type .\auxil\.ParaMonteBanner
REM echo. 

echo. 
echo. ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
echo. ::::                                                                                                                            ::::
echo.                                       ParaMonte library version !ParaMonteVersion! build on Windows
echo.                                                         !ParaMonteRelease!
echo. ::::                                                                                                                            ::::
echo. ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
echo.

set BUILD_SCRIPT_NAME=ParaMonteBuild

:: set up flags via buildFlags.bat
echo. -- !BUILD_SCRIPT_NAME! - configuring ParaMonte build...
call configParaMonte.bat || (
    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: ParaMonte library build-flag setup failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)

if !ERRORLEVEL!==1 (
    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: ParaMonte library build-flag setup failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
cd %~dp0

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: check for dynamic caf parallelism
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

setlocal EnableDelayedExpansion
set TEMP_DYNAMIC_OR_CFI=
if !LTYPE!==dynamic set TEMP_DYNAMIC_OR_CFI=true
if !CFI_ENABLED!==true set TEMP_DYNAMIC_OR_CFI=true
if !TEMP_DYNAMIC_OR_CFI!==true (
    if !CAFTYPE! NEQ none (
        echo. 
        echo. -- !BUILD_SCRIPT_NAME! - ParaMonte CAF parallelism with dynamic library build or C-interface is not supported.
        echo. -- !BUILD_SCRIPT_NAME! - skipping...
        echo. 
        cd %~dp0
        set ERRORLEVEL=0
        exit /B 0
    )
)
endlocal

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set language interface preprocessor's flag
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set "FPP_CFI_FLAG="
if !CFI_ENABLED!==true (
    set FPP_CFI_FLAG=/define:CFI_ENABLED
    REM XXX INTERFACE_LANGUAGE probably needs special care here
    if not defined INTERFACE_LANGUAGE set INTERFACE_LANGUAGE=c
) else (
    if not defined INTERFACE_LANGUAGE set INTERFACE_LANGUAGE=fortran
)

set "FPP_LANG_FLAG="
if !INTERFACE_LANGUAGE!==c set FPP_LANG_FLAG=/define:C_ENABLED
if !INTERFACE_LANGUAGE!==c++ set FPP_LANG_FLAG=/define:CPP_ENABLED
if !INTERFACE_LANGUAGE!==fortran set FPP_LANG_FLAG=/define:FORTRAN_ENABLED
if !INTERFACE_LANGUAGE!==matlab set FPP_LANG_FLAG=/define:MATLAB_ENABLED
if !INTERFACE_LANGUAGE!==python set FPP_LANG_FLAG=/define:PYTHON_ENABLED

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: report build spec and setup flags
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:: echo. operating system
echo. -- !BUILD_SCRIPT_NAME! - operating system / platform: !OS! / !PLATFORM!

:: set compiler suite
echo. -- !BUILD_SCRIPT_NAME! - COMPILER_SUITE=!COMPILER_SUITE!
if !COMPILER_SUITE! NEQ intel (
    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: the specified compiler suite !COMPILER_SUITE! is not supported.
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: only COMPILER_SUITE=intel corresponding to Intel Parallel Studio for Windows is currently supported.
    echo. -- !BUILD_SCRIPT_NAME! - gracefully exiting.
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)

:: set up compiler version

if not defined COMPILER_VERSION (

    echo. -- !BUILD_SCRIPT_NAME! - Detecting intel compiler version...
    cd .\bmake\
    call getCompilerVersion.bat || (
        echo. 
        echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: failed to fetch the Fortran compiler version.
        echo. -- !BUILD_SCRIPT_NAME! - gracefully exiting.
        echo. 
        cd %~dp0
        set ERRORLEVEL=1
        exit /B 1
    )
    cd %~dp0

)

echo. -- !BUILD_SCRIPT_NAME! - COMPILER_VERSION: !COMPILER_VERSION!
echo. -- !BUILD_SCRIPT_NAME! - TEST_RUN_ENABLED: !ParaMonteTest_RUN_ENABLED!
echo. -- !BUILD_SCRIPT_NAME! - CFI_ENABLED: !CFI_ENABLED!
echo. -- !BUILD_SCRIPT_NAME! - build type: !BTYPE!
echo. -- !BUILD_SCRIPT_NAME! - link type: !LTYPE!

:: set shared library Fortran linker flags

REM set FC_LIB_FLAGS=/libs:dll /threads %= these flags are actually included by default in recent ifort implementations =%
set FC_LIB_FLAGS=/threads /libs:static
set FL_LIB_FLAGS=/threads /libs:static
set FPP_DLL_FLAGS=

set MULTITHREADING=
echo.!FL_LIB_FLAGS! | find /I "threads">Nul && ( set "MULTITHREADING=mt" )

if !ParaMonte_LIB_ENABLED!==true (
    if !LTYPE!==dynamic (
        set FPP_DLL_FLAGS=/define:DLL_ENABLED
        REM set FC_LIB_FLAGS=!FC_LIB_FLAGS! /libs:dll
        set FL_LIB_FLAGS=!FL_LIB_FLAGS! /dll
    ) else (
        REM set FC_LIB_FLAGS=!FC_LIB_FLAGS! /libs:static
        REM if !BTYPE!==debug (
        REM     set FC_LIB_FLAGS=!FC_LIB_FLAGS! %= /dbglibs not added automatically by Intel compiler =%
        REM )
        if !BTYPE!==release (
            set FL_LIB_FLAGS=!FL_LIB_FLAGS! /Qipo-c %= generate single optimized object file if needed =%
        )
    )
)
echo. -- !BUILD_SCRIPT_NAME! - ParaMonte library linker link flags: !FL_LIB_FLAGS!
echo. -- !BUILD_SCRIPT_NAME! - ParaMonte library compiler link flags: !FC_LIB_FLAGS!
echo.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set required root directories
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

echo.
echo. -- !BUILD_SCRIPT_NAME! - setting up required root directories...

:: define ParaMonte_ROOT_DIR: contains the last backward slash

REM if not defined ParaMonte_ROOT_DIR set ParaMonte_ROOT_DIR=%~dp0
set ParaMonte_ROOT_DIR=%~dp0
echo. -- !BUILD_SCRIPT_NAME! - project root directory: !ParaMonte_ROOT_DIR!

:: set the ParaMonte source files directory

set ParaMonte_SRC_DIR=!ParaMonte_ROOT_DIR!src\ParaMonte
if exist !ParaMonte_SRC_DIR! (
    echo. -- !BUILD_SCRIPT_NAME! - source files directory: !ParaMonte_SRC_DIR!
) else (
    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: source files directory does not exist: !ParaMonte_SRC_DIR!
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
echo.

:: set the ParaMonte interface source files directory and loop over them to ensure their existence

set       ParaMonteInterface_SRC_DIR=!ParaMonte_ROOT_DIR!src\interface
set      ParaMonteInterfaceC_SRC_DIR=!ParaMonteInterface_SRC_DIR!\C
set    ParaMonteInterfaceCPP_SRC_DIR=!ParaMonteInterface_SRC_DIR!\C++
set ParaMonteInterfaceMATLAB_SRC_DIR=!ParaMonteInterface_SRC_DIR!\MATLAB
set ParaMonteInterfacePython_SRC_DIR=!ParaMonteInterface_SRC_DIR!\Python

echo.
echo. -- !BUILD_SCRIPT_NAME! - interface source files directories: !ParaMonteInterface_SRC_DIR!
for %%A in (
    !ParaMonteInterface_SRC_DIR!
    !ParaMonteInterfaceC_SRC_DIR!
    !ParaMonteInterfaceCPP_SRC_DIR!
    !ParaMonteInterfaceMATLAB_SRC_DIR!
    !ParaMonteInterfacePython_SRC_DIR!
    ) do (  if exist %%A (
                echo. -- !BUILD_SCRIPT_NAME! - %%A exists.
            ) else (
                echo. 
                echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: interface source files directory does not exist: %%A
                echo. 
                cd %~dp0
                set ERRORLEVEL=1
                exit /B 1
            )
)
echo.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set preprocessor build flags
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set FPP_BUILD_FLAGS=/define:OS_IS_WINDOWS
if !BTYPE!==debug set FPP_BUILD_FLAGS=!FPP_BUILD_FLAGS! /define:DBG_ENABLED

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set default C/CPP/Fortran compilers/linkers
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set FPP_FCL_FLAGS=
if !COMPILER_SUITE!==intel (
    set CCL=icl
    set FCL=ifort
    set FPP_FCL_FLAGS=/define:IFORT_ENABLED
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set up preprocessor flags
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
REM echo. FPP_FLAGS_EXTRA = !FPP_FLAGS_EXTRA!
REM /define:IS_ENABLED
set FPP_FLAGS=/fpp /define:PARAMONTE_VERSION=^"'!ParaMonteVersion!'^" !FPP_CFI_FLAG! !FPP_LANG_FLAG! !FPP_BUILD_FLAGS! !FPP_FCL_FLAGS! !FPP_DLL_FLAGS! !USER_PREPROCESSOR_MACROS! !FPP_FLAGS_EXTRA!
REM set FPP_FLAGS=/fpp !FPP_CFI_FLAG! !FPP_LANG_FLAG! !FPP_BUILD_FLAGS! !FPP_FCL_FLAGS! !FPP_DLL_FLAGS! !USER_PREPROCESSOR_MACROS! !FPP_FLAGS_EXTRA!
:: to save the intermediate files use this on the command line: FPP /Qsave_temps <original file> <intermediate file>

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set up coarray flags
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

echo.
echo. -- !BUILD_SCRIPT_NAME! - setting up Coarray Fortran (CAF) parallelization model. Options: single, shared, distributed
echo. -- !BUILD_SCRIPT_NAME! - requested CAF: !CAFTYPE!

set CAF_ENABLED=false
if !CAFTYPE!==single set CAF_ENABLED=true
if !CAFTYPE!==shared set CAF_ENABLED=true
if !CAFTYPE!==distributed set CAF_ENABLED=true

if !CAF_ENABLED!==true (
    set SERIAL_ENABLED=false
    echo. -- !BUILD_SCRIPT_NAME! - enabling Coarray Fortran syntax via preprocesor flag /define:CAF_ENABLED
    set FPP_FLAGS=!FPP_FLAGS! /define:CAF_ENABLED
    set CAF_FLAGS=/Qcoarray=!CAFTYPE!
    if not defined FOR_COARRAY_NUM_IMAGES set FOR_COARRAY_NUM_IMAGES=3
    echo. -- !BUILD_SCRIPT_NAME! - number of Coarray images: !FOR_COARRAY_NUM_IMAGES!
) else (
    echo. -- !BUILD_SCRIPT_NAME! - ignoring Coarray Fortran parallelization.
    set CAF_FLAGS=
    set CAFTYPE=
)

echo. -- !BUILD_SCRIPT_NAME! - Coarray Fortran flags: !CAF_FLAGS!
echo.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set non-coarray parallelization flags and definitions to be passed to the preprocessors
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set MPI_FLAGS=
if !MPI_ENABLED!==true (
    set SERIAL_ENABLED=false
    if not defined CAFTYPE (
        set FPP_FLAGS=!FPP_FLAGS! /define:MPI_ENABLED
        REM set MPI_FLAGS=-fast
        set FCL=mpiifort.bat -fc=ifort
        set CCL=mpicc -cc=icl.exe
    ) else (
        echo.
        echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: Coarray Fortran cannot be mixed with MPI.
        echo. -- !BUILD_SCRIPT_NAME! - CAFTYPE: !CAFTYPE!
        echo. -- !BUILD_SCRIPT_NAME! - MPI_ENABLED: !MPI_ENABLED!
        echo. -- !BUILD_SCRIPT_NAME! - set MPI_ENABLED and CAFTYPE to appropriate values in the ParaMonte config file and rebuild.
        echo.
        cd %~dp0
        set ERRORLEVEL=1
        exit /B 1
    )
)

set OMP_FLAGS=
if !OMP_ENABLED!==true (
    set OMP_FLAGS=/Qopenmp
    set SERIAL_ENABLED=false
)

set FCL_PARALLELIZATION_FLAGS=!CAF_FLAGS! !MPI_FLAGS! !OMP_FLAGS!
if !SERIAL_ENABLED!==true (
    REM if !BTYPE!==testing set FCL_PARALLELIZATION_FLAGS=!FCL_PARALLELIZATION_FLAGS! /Qparallel
    REM if !BTYPE!==release set FCL_PARALLELIZATION_FLAGS=!FCL_PARALLELIZATION_FLAGS! /Qparallel
)
echo. -- !BUILD_SCRIPT_NAME! - all compiler/linker parallelization flags: !FCL_PARALLELIZATION_FLAGS!
echo.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set set default Fortran compiler flags in different build modes.
:: Complete list of intel compiler options:
:: https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-alphabetical-list-of-compiler-options
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if !COMPILER_SUITE!==intel (

    ::  /QxHost
    set FCL_FLAGS_DEFAULT=/nologo /standard-semantics /F0x1000000000
    
    if !BTYPE!==debug set FCL_BUILD_FLAGS=!INTEL_FORTRAN_DEBUG_FLAGS! /stand:f08
    
    if !BTYPE!==release set FCL_BUILD_FLAGS=!INTEL_FORTRAN_RELEASE_FLAGS!
    
    :: set Fortran linker flags for release mode
    if !BTYPE!==release set FL_FLAGS=/Qopt-report:2
    if !BTYPE!==testing set FL_FLAGS=
    if !BTYPE!==debug   set FL_FLAGS=
    REM /Qipo-c:
    REM      Tells the compiler to optimize across multiple files and generate a single object file ipo_out.obj without linking
    REM      info at: https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-ipo-c-qipo-c
    REM
    
    if !BTYPE!==testing set FCL_BUILD_FLAGS=!INTEL_FORTRAN_TESTING_FLAGS!

) else (

    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: No compiler other than Intel Parallel Studio is suppoerted on Windows. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1

)

set FCL_FLAGS=!FCL_FLAGS_DEFAULT! !FCL_PARALLELIZATION_FLAGS! !FCL_BUILD_FLAGS!

if !HEAP_ARRAY_ENABLED!==true (
    set FCL_FLAGS=!FCL_FLAGS! /heap-arrays
)

echo.
echo. -- !BUILD_SCRIPT_NAME! - Fortran preprocessor flags: !FPP_FLAGS!
echo. -- !BUILD_SCRIPT_NAME! - Fortran linker library flags: !FL_LIB_FLAGS!
echo. -- !BUILD_SCRIPT_NAME! - Fortran compiler library flags: !FC_LIB_FLAGS!
echo. -- !BUILD_SCRIPT_NAME! - Fortran compiler/linker all flags: !FCL_FLAGS!
echo. -- !BUILD_SCRIPT_NAME! - Fortran compiler/linker default flags: !FCL_FLAGS_DEFAULT!
echo. -- !BUILD_SCRIPT_NAME! - Fortran compiler/linker flags in !BTYPE! build mode: !FCL_BUILD_FLAGS!
echo.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: check MATLAB's existence
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:: it is imperative to nullify these MATLAB variables for all language builds at all times

set "MATLAB_ROOT_DIR="
set "MATLAB_EXE_PATH="
set "MATLAB_BIN_DIR="
set "MATLAB_LIB_DIR="
set "MATLAB_INC_DIR=."
set "MATLAB_LIBMX_FILE="
set "MATLAB_LIBMEX_FILE="
set "MATLAB_LIBMAT_FILE="
set "MATLAB_VERSION_FILE="
REM set "MATLAB_INC_DIR_FLAG="

if not !INTERFACE_LANGUAGE!==matlab goto LABEL_BUILD_ParaMonte

echo. 
echo. -- !BUILD_SCRIPT_NAME! - searching for a MATLAB installation on your system...

set "INSTALL_LOC_LIST=C:\Program Files\MATLAB\/C:\Program Files (x86)\MATLAB\"
set MATLAB_VERSION_LIST=R2025b/R2025a/R2024b/R2024a/R2023b/R2023a/R2022b/R2022a/R2021b/R2021a/R2020b/R2020a/R2019b/R2019a/R2018b/R2018a/R2017b/R2017a

for %%D in ("!INSTALL_LOC_LIST:/=" "!") do (
    for %%V in ("!MATLAB_VERSION_LIST:/=" "!") do (
        set "MATLAB_ROOT_DIR_TEMP=%%~D%%~V"
        set "MATLAB_BIN_DIR_TEMP=!MATLAB_ROOT_DIR_TEMP!\bin"
        set "MATLAB_EXE_PATH_TEMP=!MATLAB_BIN_DIR_TEMP!\matlab.exe"
        if exist !MATLAB_EXE_PATH_TEMP! (
            set "MATLAB_ROOT_DIR=!MATLAB_ROOT_DIR_TEMP!"
            set "MATLAB_EXE_PATH=!MATLAB_EXE_PATH_TEMP!"
            set "MATLAB_BIN_DIR=!MATLAB_BIN_DIR_TEMP!"
            set "MATLAB_INC_DIR=!MATLAB_ROOT_DIR!\extern\include"
            set "MATLAB_LIB_DIR=!MATLAB_ROOT_DIR!\extern\lib\win64\microsoft"
            set "MATLAB_LIBMX_FILE=!MATLAB_LIB_DIR!\libmx.lib"
            set "MATLAB_LIBMEX_FILE=!MATLAB_LIB_DIR!\libmex.lib"
            set "MATLAB_LIBMAT_FILE=!MATLAB_LIB_DIR!\libmat.lib"
            set "MATLAB_VERSION_FILE=!MATLAB_ROOT_DIR!\extern\version\fortran_mexapi_version.F"
            set FPP_FLAGS=!FPP_FLAGS! /define:MEXPRINT_ENABLED /define:MATLAB_MEX_FILE
            REM set "MATLAB_INC_DIR_FLAG=/I:!MATLAB_INC_DIR!"
            echo. -- !BUILD_SCRIPT_NAME! - MATLAB %%~V installation detected at: !MATLAB_EXE_PATH!
            echo. 
            goto :LABEL_continue
        )
    )
)

echo. -- !BUILD_SCRIPT_NAME! - WARNING: Exhausted all possible search paths for a MATLAB installation, but failed to find MATLAB.
echo. -- !BUILD_SCRIPT_NAME! - WARNING: The ParaMonte MATLAB kernel will not be functional without building the required DLL libraries.
echo. -- !BUILD_SCRIPT_NAME! - WARNING: Please add MATLAB to your environmental variable PATH and rerun the install script.
echo. -- !BUILD_SCRIPT_NAME! - WARNING: For example, on your current Windows command-line, try:
echo. -- !BUILD_SCRIPT_NAME! - WARNING: 
echo. -- !BUILD_SCRIPT_NAME! - WARNING:     set "PATH=PATH_TO_MATLAB_BIN_DIR;!PATH!
echo. -- !BUILD_SCRIPT_NAME! - WARNING: 
echo. -- !BUILD_SCRIPT_NAME! - WARNING: where PATH_TO_MATLAB_BIN_DIR must be replaced with path to the bin folder of the current 
echo. -- !BUILD_SCRIPT_NAME! - WARNING: installation of MATLAB on your system. Typical MATLAB bin installation path on a 64-bit Windows 
echo. -- !BUILD_SCRIPT_NAME! - WARNING: Operating Systems is a string like the following:
echo. -- !BUILD_SCRIPT_NAME! - WARNING: 
echo. -- !BUILD_SCRIPT_NAME! - WARNING:     C:\Program Files\MATLAB\2020a\bin\
echo. -- !BUILD_SCRIPT_NAME! - WARNING: 
echo. -- !BUILD_SCRIPT_NAME! - WARNING: where 2020a in the path points to the MATLAB 2020a version installation on the system. You can also 
echo. -- !BUILD_SCRIPT_NAME! - WARNING: find the installation location of MATLAB by typing the following command in your MATLAB session:
echo. -- !BUILD_SCRIPT_NAME! - WARNING: 
echo. -- !BUILD_SCRIPT_NAME! - WARNING:     matlabroot
echo. -- !BUILD_SCRIPT_NAME! - WARNING: 
echo. -- !BUILD_SCRIPT_NAME! - WARNING: skipping the ParaMonte MATLAB build...

REM cd %~dp0
REM set ERRORLEVEL=1
REM exit /B 1

:LABEL_continue

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: generate ParaMonte library build directories, object files, and dynamic libraries
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:LABEL_BUILD_ParaMonte

call !ParaMonte_SRC_DIR!\buildParaMonteSource.bat || (
    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: build failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
if !ERRORLEVEL!==1 (
    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: build failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
cd %~dp0

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: generate ParaMonte library TEST build directories and object files
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:LABEL_ParaMonteTest

set ParaMonteTest_ENABLED=false
if !ParaMonteTest_OBJ_ENABLED!==true set ParaMonteTest_ENABLED=true
if !ParaMonteTest_EXE_ENABLED!==true set ParaMonteTest_ENABLED=true
if !ParaMonteTest_RUN_ENABLED!==true set ParaMonteTest_ENABLED=true

:: set path to ParaMonte test source files

set ParaMonteTest_SRC_DIR=!ParaMonte_ROOT_DIR!src\test
if exist !ParaMonteTest_SRC_DIR! (
    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - ParaMonte library test source files directory: !ParaMonteTest_SRC_DIR!
) else (
    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: ParaMonte library test source files directory does not exist: !ParaMonteTest_SRC_DIR!
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
echo.

if !ParaMonteTest_ENABLED! NEQ true goto LABEL_ParaMonteInterface

call !ParaMonteTest_SRC_DIR!\buildParaMonteTest.bat || (
    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: the ParaMonte library test failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
if !ERRORLEVEL!==1 (
    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: the ParaMonte library test failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
cd %~dp0

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: generate ParaMonte library interface build directories and object files
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:LABEL_ParaMonteInterface

REM if !ParaMonteExample_ENABLED! NEQ true goto LABEL_EOF
set ParaMonteInterfaceC_SRC_DIR=!ParaMonte_ROOT_DIR!src\interface\C
set ParaMonteInterfaceCPP_SRC_DIR=!ParaMonte_ROOT_DIR!src\interface\C++
set ParaMonteInterfaceFortran_SRC_DIR=!ParaMonte_ROOT_DIR!src\interface\Fortran

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: generate ParaMonte library MATLAB build directories and object files
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if !LTYPE!==static goto LABEL_ParaMontePython
if !CFI_ENABLED! NEQ true goto LABEL_ParaMontePython
if !INTERFACE_LANGUAGE! NEQ matlab goto LABEL_ParaMontePython

:: setup MATLAB library source files directory

set ParaMonteInterfaceMATLAB_SRC_DIR=!ParaMonte_ROOT_DIR!src\interface\MATLAB
if exist !ParaMonteInterfaceMATLAB_SRC_DIR! (
    echo. -- !BUILD_SCRIPT_NAME! - MATLAB source files directory: !ParaMonteInterfaceMATLAB_SRC_DIR!
) else (
    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: MATLAB source files directory does not exist: !ParaMonteInterfaceMATLAB_SRC_DIR!
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
echo.

call !ParaMonteInterfaceMATLAB_SRC_DIR!\buildParaMonteMATLAB.bat || (
    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: the ParaMonte library MATLAB build failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
if !ERRORLEVEL!==1 (
    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: the ParaMonte library MATLAB build failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
cd %~dp0

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: generate ParaMonte library Python build directories and object files
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:LABEL_ParaMontePython

if !LTYPE!==static goto LABEL_ParaMonteExamples
if !CFI_ENABLED! NEQ true goto LABEL_ParaMonteExamples
if !INTERFACE_LANGUAGE! NEQ python goto LABEL_ParaMonteExamples

:: setup Python library source files directory

set ParaMonteInterfacePython_SRC_DIR=!ParaMonte_ROOT_DIR!src\interface\Python
if exist !ParaMonteInterfacePython_SRC_DIR! (
    echo. -- !BUILD_SCRIPT_NAME! - Python source files directory: !ParaMonteInterfacePython_SRC_DIR!
) else (
    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: Python source files directory does not exist: !ParaMonteInterfacePython_SRC_DIR!
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
echo.

call !ParaMonteInterfacePython_SRC_DIR!\buildParaMontePython.bat || (
    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: the ParaMonte library Python build failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
if !ERRORLEVEL!==1 (
    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: the ParaMonte library Python build failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
cd %~dp0

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: generate ParaMonte library Example build directories, object files, and executables
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:LABEL_ParaMonteExamples

REM set ParaMonteExample_ENABLED=false
REM if !ParaMonteExample_EXE_ENABLED!==true set ParaMonteExample_ENABLED=true
REM if !ParaMonteExample_RUN_ENABLED!==true set ParaMonteExample_ENABLED=true

:: set path to ParaMonte example source files

set ParaMonteExample_SRC_DIR=!ParaMonte_ROOT_DIR!example
if exist !ParaMonteExample_SRC_DIR! (
    echo. -- !BUILD_SCRIPT_NAME! - ParaMonte library example source files directory: !ParaMonteExample_SRC_DIR!
) else (
    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: ParaMonte library example source files directory does not exist: !ParaMonteExample_SRC_DIR!
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
echo.

:: build ParaMonte examples

echo.
call !ParaMonteExample_SRC_DIR!\buildParaMonteExample.bat || (
    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: the ParaMonte library examples build failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
if !ERRORLEVEL!==1 (
    echo. 
    echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: the ParaMonte library examples build failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
if !ERRORLEVEL!==0 (
    echo.
    echo.
    echo. -- !BUILD_SCRIPT_NAME! - build successful. 
    echo.
)
echo.

cd %~dp0

:LABEL_EOF

:: go back to the root folder containing this batch file 

cd %~dp0

:: undefine all configuration environmental flags

if !ParaMonte_FLAG_CLEANUP_ENABLED!==true (
    call unconfigParaMonte.bat || (
        echo. 
        echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: the ParaMonte library cleanup failed. exiting...
        echo. 
        cd %~dp0
        set ERRORLEVEL=1
        exit /B 1
    )
)

cd !ParaMonte_BLD_ROOT_DIR!

:: endlocal

exit /B 0



