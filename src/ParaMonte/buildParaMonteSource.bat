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

:: NOTE: This is not a standalone build-script. It must only be called by buildParaMonte.bat script in the root directory of the project.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: generate ParaMonte library build directories
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:: setlocal EnableDelayedExpansion

cd %~dp0

set PARALLELIZATION_DIR=
if !OMP_ENABLED!==true set PARALLELIZATION_DIR=!PARALLELIZATION_DIR!omp
if !MPI_ENABLED!==true set PARALLELIZATION_DIR=!PARALLELIZATION_DIR!mpi
if !CAF_ENABLED!==true set PARALLELIZATION_DIR=!PARALLELIZATION_DIR!caf!CAFTYPE!
if not defined PARALLELIZATION_DIR set PARALLELIZATION_DIR=serial

echo.
echo. -- ParaMonte - generating build directories...
echo.

set MEMORY_ALLOCATION=stack
if !HEAP_ARRAY_ENABLED!==true set MEMORY_ALLOCATION=heap

:: define ParaMonte_BLD_DIR (win64, win32), compiler (intel), build type (debug, release, testing), lib type (static, dynamic), parallelization mode
set ParaMonte_BLD_DIR=!ParaMonte_BLD_ROOT_DIR!\build\win!PLATFORM!\!COMPILER_SUITE!\!COMPILER_VERSION!\!BTYPE!\!LTYPE!\!MEMORY_ALLOCATION!\!PARALLELIZATION_DIR!
if !CFI_ENABLED!==true (
    set ParaMonte_BLD_DIR=!ParaMonte_BLD_DIR!\C
)else (
    set ParaMonte_BLD_DIR=!ParaMonte_BLD_DIR!\Fortran
)
if exist !ParaMonte_BLD_DIR! (
    echo. -- ParaMonte - !ParaMonte_BLD_DIR! already exists. skipping...
) else (
    echo. -- ParaMonte - generating !BTYPE! library directory: !ParaMonte_BLD_DIR!
    mkdir !ParaMonte_BLD_DIR!
)
echo. -- ParaMonte - all generated build files will be stored at: %~dp0!ParaMonte_BLD_DIR!

:: generate object files directory
set ParaMonte_OBJ_DIR=!ParaMonte_BLD_DIR!\obj
if exist !ParaMonte_OBJ_DIR! (
    echo. -- ParaMonte - !ParaMonte_OBJ_DIR! already exists. skipping...
) else (
    echo. -- ParaMonte - generating object files directory: !ParaMonte_OBJ_DIR!
    mkdir !ParaMonte_OBJ_DIR!
)

:: generate modules files directory
set ParaMonte_MOD_DIR=!ParaMonte_BLD_DIR!\mod
if exist !ParaMonte_MOD_DIR! (
    echo. -- ParaMonte - !ParaMonte_MOD_DIR! already exists. skipping...
) else (
    echo. -- ParaMonte - generating module files directory: !ParaMonte_MOD_DIR!
    mkdir !ParaMonte_MOD_DIR!
)

:: generate library files directory
set ParaMonte_LIB_DIR=!ParaMonte_BLD_DIR!\lib
if exist !ParaMonte_LIB_DIR! (
    echo. -- ParaMonte - !ParaMonte_LIB_DIR! already exists. skipping...
) else (
    echo. -- ParaMonte - generating output library files directory: !ParaMonte_LIB_DIR!
    mkdir !ParaMonte_LIB_DIR!
)
echo.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: generate ParaMonte library object files
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

cd !ParaMonte_OBJ_DIR!

:: Read the name of each file from the ordered list of filenames in filelist.txt to compile

if !ParaMonte_OBJ_ENABLED! NEQ true goto LABEL_ParaMonte_LIB_ENABLED

echo.
echo. -- ParaMonte - generating object files at: !ParaMonte_OBJ_DIR!

:: First verify the source filelist exists

set ParaMonte_FILE_LIST=!ParaMonte_SRC_DIR!\filelist.txt
if not exist !ParaMonte_FILE_LIST! (
    echo.
    echo. -- ParaMonte - Fatal Error: The filelist.txt containing the ParaMonte source filenames does not exist. Path: !ParaMonte_FILE_LIST!
    echo. -- ParaMonte - build failed. exiting...
    echo.
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)

:: generate ParaMonte library object files

if !MPI_ENABLED!==true (
    for /F "eol=! tokens=*" %%A in (!ParaMonte_FILE_LIST!) do (
        echo. -- ParaMonte - generating object file for %%A
        call !FCL! !FCL_FLAGS! !FPP_FLAGS! !FC_LIB_FLAGS! ^
        /module:!ParaMonte_MOD_DIR!     %=path to output ParaMonte example module files=% ^
        /I:!ParaMonte_MOD_DIR!          %=path to input ParaMonte module files=%  ^
        /c !ParaMonte_SRC_DIR!\%%A      %=path to input ParaMonte example source files=%
        if !ERRORLEVEL! NEQ 0 (
            echo.
            echo. -- ParaMonte - Fatal Error: compilation of the object file for %%A failed.
            echo. -- ParaMonte - build failed. exiting...
            echo.
            cd %~dp0
            set ERRORLEVEL=1
            exit /B 1
        )
    )
    @echo off
    echo.
) else (
    for /F "eol=! tokens=*" %%A in (!ParaMonte_FILE_LIST!) do (
        echo. -- ParaMonte - generating object file for %%A
        !FCL! !FCL_FLAGS! !FPP_FLAGS! !FCL_IPOC_FLAG! !FC_LIB_FLAGS! ^
        /module:!ParaMonte_MOD_DIR!     %=path to output ParaMonte example module files=% ^
        /I:!ParaMonte_MOD_DIR!          %=path to input ParaMonte module files=%  ^
        /c !ParaMonte_SRC_DIR!\%%A      %=path to input ParaMonte example source files=% ^
        || (
            echo.
            echo. -- ParaMonte - Fatal Error: compilation of the object file for %%A failed.
            echo. -- ParaMonte - build failed. exiting...
            echo.
            cd %~dp0
            set ERRORLEVEL=1
            exit /B 1
        )
    )
    @echo off
    echo.
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: generate dynamic library files if requested
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:LABEL_ParaMonte_LIB_ENABLED

set PMLIB_NAME=libparamonte_!LTYPE!_!MEMORY_ALLOCATION!_!BTYPE!_!COMPILER_SUITE!_!INTERFACE_LANGUAGE!
if !CAF_ENABLED!==true set PMLIB_NAME=!PMLIB_NAME!_caf!CAFTYPE!
if !MPI_ENABLED!==true set PMLIB_NAME=!PMLIB_NAME!_mpi
if !OMP_ENABLED!==true set PMLIB_NAME=!PMLIB_NAME!_omp
set PMLIB_NAME=!PMLIB_NAME!_windows_!PLATFORM!
if defined MULTITHREADING set PMLIB_NAME=!PMLIB_NAME!_!MULTITHREADING!


:: if !ParaMonte_LIB_ENABLED! NEQ true goto LABEL_ParaMonteTest_OBJ_ENABLED
if !ParaMonte_LIB_ENABLED! NEQ true goto LABEL_EOF

REM if !LTYPE! NEQ dynamic goto LABEL_EOF

cd !ParaMonte_LIB_DIR!
echo.
echo. -- ParaMonte - building ParaMonte !LTYPE! library files at: !ParaMonte_LIB_DIR!
echo. -- ParaMonte - !LTYPE! library filename: !PMLIB_NAME!
echo.

if !LTYPE!==dynamic (

    !FCL! !FCL_FLAGS! !FL_FLAGS! !FL_LIB_FLAGS! ^
    /module:!ParaMonte_MOD_DIR!                                 %=path to output ParaMonte example module files=% ^
    /I:!ParaMonte_MOD_DIR! !ParaMonte_OBJ_DIR!\*.obj            %=path to input ParaMonte module files=%  ^
    /exe:!PMLIB_NAME!.dll

) else (

    REM https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-creating-static-libraries

    if !BTYPE!==release (

        !FCL! !FCL_FLAGS! !FL_FLAGS! !FL_LIB_FLAGS! ^
        /module:!ParaMonte_MOD_DIR!                                 %=path to output ParaMonte example module files=% ^
        /I:!ParaMonte_MOD_DIR! !ParaMonte_OBJ_DIR!\*.obj            %=path to input ParaMonte module files=%

        xilib /out:!PMLIB_NAME!.lib ipo_out.obj
        REM del !ParaMonte_OBJ_DIR!\ipo_out.obj

    ) else (

        xilib /out:!PMLIB_NAME!.lib !ParaMonte_OBJ_DIR!\*.obj

    )

)

if !ERRORLEVEL!==0 goto LABEL_EOF

echo. 
echo. -- ParaMonte - Fatal Error: ParaMonte dynamic library (DLL) generation failed.
echo. -- ParaMonte - build failed. exiting...
echo. 
cd %~dp0
set ERRORLEVEL=1
exit /B 1
echo.

:LABEL_EOF

cd %~dp0

exit /B 0
