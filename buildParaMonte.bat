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

:: define variables only locally
::setlocal

set ParaMonte_BLD_ROOT_DIR=%cd%

:: change directory to the folder containing this batch file 
cd %~dp0

:: fetch ParaMonte library version

cd .\bmake\
for /f "tokens=*" %%i in ('head.bat 1 "..\.VERSION"') do set "ParaMonteVersion=%%i"
cd %~dp0

REM echo. 
REM type .\auxil\ParaMonteBanner.txt
REM echo. 

echo. 
echo. ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
echo. ::::                                                                                                                            ::::
echo.                                       ParaMonte library version !ParaMonteVersion! build on Windows
echo. ::::                                                                                                                            ::::
echo. ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
echo.

:: set up flags via buildFlags.bat
echo. -- ParaMonte - configuring ParaMonte build...
call configParaMonte.bat

if !ERRORLEVEL!==1 (
    echo. 
    echo. -- ParaMonte - Fatal Error: ParaMonte library build-flag setup failed. exiting...
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
        echo. -- ParaMonte - ParaMonte CAF parallelism with dynamic library build or C-interface is not supported.
        echo. -- ParaMonte - skipping...
        echo. 
        cd %~dp0
        set ERRORLEVEL=0
        exit /B 0
    )
)
endlocal

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set interoperability mode preprocessor's flag
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set FPP_CFI_FLAG=
set INTERFACE_LANGUAGE=fortran
if !CFI_ENABLED!==true (
    set FPP_CFI_FLAG=/define:CFI_ENABLED
    set INTERFACE_LANGUAGE=c
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: report build spec and setup flags
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:: echo. operating system
echo. -- ParaMonte - operating system / platform: !OS! / !PLATFORM!

:: set compiler suite
echo. -- ParaMonte - COMPILER_SUITE=!COMPILER_SUITE!
if !COMPILER_SUITE! NEQ intel (
    echo. 
    echo. -- ParaMonte - Fatal Error: the specified compiler suite !COMPILER_SUITE! is not supported.
    echo. -- ParaMonte - Fatal Error: only COMPILER_SUITE=intel corresponding to Intel Parallel Studio for Windows is currently supported.
    echo. -- ParaMonte - gracefully exiting.
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)

:: set up compiler version

if not defined COMPILER_VERSION (

    echo. -- ParaMonte - Detecting intel compiler version...
    cd .\bmake\
    call getCompilerVersion.bat
    cd %~dp0

)

echo. -- ParaMonte - COMPILER_VERSION: !COMPILER_VERSION!
echo. -- ParaMonte - TEST_RUN_ENABLED: !ParaMonteTest_RUN_ENABLED!
echo. -- ParaMonte - CFI_ENABLED: !CFI_ENABLED!
echo. -- ParaMonte - build type: !BTYPE!
echo. -- ParaMonte - link type: !LTYPE!

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
echo. -- ParaMonte - ParaMonte library linker link flags: !FL_LIB_FLAGS!
echo. -- ParaMonte - ParaMonte library compiler link flags: !FC_LIB_FLAGS!
echo.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set required root directories
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

echo.
echo. -- ParaMonte - setting up required root directories...

:: define ParaMonte_ROOT_DIR: contains the last backward slash

if not defined ParaMonte_ROOT_DIR set ParaMonte_ROOT_DIR=%~dp0
echo. -- ParaMonte - project root directory: !ParaMonte_ROOT_DIR!

:: set the ParaMonte source files directory

set ParaMonte_SRC_DIR=!ParaMonte_ROOT_DIR!src\ParaMonte
if exist !ParaMonte_SRC_DIR! (
    echo. -- ParaMonte - source files directory: !ParaMonte_SRC_DIR!
) else (
    echo. 
    echo. -- ParaMonte - Fatal Error: source files directory does not exist: !ParaMonte_SRC_DIR!
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
echo.

:: set the ParaMonte interface source files directory and loop over them to ensure their existence

set       ParaMonteInterface_SRC_DIR=!ParaMonte_ROOT_DIR!src\interface
set      ParaMonteInterfaceC_SRC_DIR=!ParaMonteInterface_SRC_DIR!\C
set ParaMonteInterfacePython_SRC_DIR=!ParaMonteInterface_SRC_DIR!\Python

echo.
echo. -- ParaMonte - interface source files directories: !ParaMonteInterface_SRC_DIR!
for %%A in (
    !ParaMonteInterface_SRC_DIR!
    !ParaMonteInterfaceC_SRC_DIR!
    !ParaMonteInterfacePython_SRC_DIR!
    ) do (  if exist %%A (
                echo. -- ParaMonte - %%A exists.
            ) else (
                echo. 
                echo. -- ParaMonte - Fatal Error: interface source files directory does not exist: %%A
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

set FPP_BUILD_FLAGS=
if !BTYPE!==debug set FPP_BUILD_FLAGS=/define:DBG_ENABLED

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

set FPP_FLAGS=/fpp !FPP_CFI_FLAG! !FPP_BUILD_FLAGS! !FPP_FCL_FLAGS! !FPP_DLL_FLAGS! !USER_PREPROCESSOR_MACROS!
:: to save the intermediate files use this on the command line: FPP /Qsave_temps <original file> <intermediate file>

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set up coarray flags
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

echo.
echo. -- ParaMonte - setting up Coarray Fortran (CAF) parallelization model. Options: single, shared, distributed
echo. -- ParaMonte - requested CAF: !CAFTYPE!

set CAF_ENABLED=false
if !CAFTYPE!==single set CAF_ENABLED=true
if !CAFTYPE!==shared set CAF_ENABLED=true
if !CAFTYPE!==distributed set CAF_ENABLED=true

if !CAF_ENABLED!==true (
    echo. -- ParaMonte - enabling Coarray Fortran syntax via preprocesor flag /define:CAF_ENABLED
    set FPP_FLAGS=!FPP_FLAGS! /define:CAF_ENABLED
    set CAF_FLAGS=/Qcoarray=!CAFTYPE!
    if not defined FOR_COARRAY_NUM_IMAGES set FOR_COARRAY_NUM_IMAGES=3
    echo. -- ParaMonte - number of Coarray images: !FOR_COARRAY_NUM_IMAGES!
) else (
    echo. -- ParaMonte - ignoring Coarray Fortran parallelization.
    set CAF_FLAGS=
    set CAFTYPE=
)

echo. -- ParaMonte - Coarray Fortran flags: !CAF_FLAGS!
echo.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set non-coarray parallelization flags and definitions to be passed to the preprocessors
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set MPI_FLAGS=
if !MPI_ENABLED!==true (
    if not defined CAFTYPE (
        set FPP_FLAGS=!FPP_FLAGS! /define:MPI_ENABLED
        REM set MPI_FLAGS=-fast
        set FCL=mpiifort.bat -fc=ifort
        set CCL=mpicc -cc=icl.exe
    ) else (
        echo.
        echo. -- ParaMonte - Fatal Error: Coarray Fortran cannot be mixed with MPI.
        echo. -- ParaMonte - CAFTYPE: !CAFTYPE!
        echo. -- ParaMonte - MPI_ENABLED: !MPI_ENABLED!
        echo. -- ParaMonte - set MPI_ENABLED and CAFTYPE to appropriate values in the ParaMonte config file and rebuild.
        echo.
        cd %~dp0
        set ERRORLEVEL=1
        exit /B 1
    )
)

set OMP_FLAGS=
if !OMP_ENABLED!==true set OMP_FLAGS=/Qopenmp
set FCL_PARALLELIZATION_FLAGS=!CAF_FLAGS! !MPI_FLAGS! !OMP_FLAGS!
echo. -- ParaMonte - all compiler/linker parallelization flags: !FCL_PARALLELIZATION_FLAGS!
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
    echo. -- ParaMonte - Fatal Error: No compiler other than Intel Parallel Studio is suppoerted on Windows. exiting...
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
echo. -- ParaMonte - Fortran preprocessor flags: !FPP_FLAGS!
echo. -- ParaMonte - Fortran linker library flags: !FL_LIB_FLAGS!
echo. -- ParaMonte - Fortran compiler library flags: !FC_LIB_FLAGS!
echo. -- ParaMonte - Fortran compiler/linker all flags: !FCL_FLAGS!
echo. -- ParaMonte - Fortran compiler/linker default flags: !FCL_FLAGS_DEFAULT!
echo. -- ParaMonte - Fortran compiler/linker flags in !BTYPE! build mode: !FCL_BUILD_FLAGS!
echo.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: generate ParaMonte library build directories, object files, and dynamic libraries
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

call !ParaMonte_SRC_DIR!\buildParaMonteSource.bat
if !ERRORLEVEL!==1 (
    echo. 
    echo. -- ParaMonte - Fatal Error: build failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
cd %~dp0

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: generate ParaMonte library Python build directories and object files
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if !LTYPE!==static goto LABEL_ParaMonteTest
if !CFI_ENABLED! NEQ true goto LABEL_ParaMonteTest

:: setup Python library source files directory

set ParaMontePython_SRC_DIR=!ParaMonte_ROOT_DIR!src\interface\Python
if exist !ParaMontePython_SRC_DIR! (
    echo. -- ParaMonte - Python source files directory: !ParaMontePython_SRC_DIR!
) else (
    echo. 
    echo. -- ParaMonte - Fatal Error: Python source files directory does not exist: !ParaMontePython_SRC_DIR!
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
echo.

call !ParaMontePython_SRC_DIR!\buildParaMontePython.bat

if !ERRORLEVEL!==1 (
    echo. 
    echo. -- ParaMonte - Fatal Error: build failed. exiting...
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
    echo. -- ParaMonte - ParaMonte library test source files directory: !ParaMonteTest_SRC_DIR!
) else (
    echo. 
    echo. -- ParaMonte - Fatal Error: ParaMonte library test source files directory does not exist: !ParaMonteTest_SRC_DIR!
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
echo.

if !ParaMonteTest_ENABLED! NEQ true goto LABEL_ParaMonteExamples

call !ParaMonteTest_SRC_DIR!\buildParaMonteTest.bat

if !ERRORLEVEL!==1 (
    echo. 
    echo. -- ParaMonte - Fatal Error: ParaMonte library test build failed. exiting...
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

set ParaMonteExample_ENABLED=false
if !ParaMonteExample_EXE_ENABLED!==true set ParaMonteExample_ENABLED=true
if !ParaMonteExample_RUN_ENABLED!==true set ParaMonteExample_ENABLED=true

:: set path to ParaMonte example source files

set ParaMonteExample_SRC_DIR=!ParaMonte_ROOT_DIR!example
if exist !ParaMonteExample_SRC_DIR! (
    echo. -- ParaMonte - ParaMonte library example source files directory: !ParaMonteExample_SRC_DIR!
) else (
    echo. 
    echo. -- ParaMonte - Fatal Error: ParaMonte library example source files directory does not exist: !ParaMonteExample_SRC_DIR!
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
echo.

REM if !ParaMonteExample_ENABLED! NEQ true goto LABEL_EOF

set ParaMonteC_SRC_DIR=!ParaMonte_ROOT_DIR!src\interface\C

:: build ParaMonte examples

echo.
call !ParaMonteExample_SRC_DIR!\buildParaMonteExample.bat
echo.

if !ERRORLEVEL!==1 (
    echo. 
    echo. -- ParaMonte - Fatal Error: ParaMonte library example build failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
if !ERRORLEVEL!==0 (
    echo.
    echo.
    echo. -- ParaMonte - build successful. 
    echo.
)

cd %~dp0

:LABEL_EOF

:: go back to the root folder containing this batch file 

cd %~dp0

:: undefine all configuration environmental flags

if !ParaMonte_FLAG_CLEANUP_ENABLED!==true (
    call unconfigParaMonte.bat 
)

cd !ParaMonte_BLD_ROOT_DIR!

:: endlocal

exit /B 0



