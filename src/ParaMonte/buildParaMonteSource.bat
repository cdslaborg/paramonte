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

if !INTERFACE_LANGUAGE!==c set ParaMonte_BLD_DIR=!ParaMonte_BLD_DIR!\C
if !INTERFACE_LANGUAGE!==c++ set ParaMonte_BLD_DIR=!ParaMonte_BLD_DIR!\C++
if !INTERFACE_LANGUAGE!==matlab set ParaMonte_BLD_DIR=!ParaMonte_BLD_DIR!\MATLAB
if !INTERFACE_LANGUAGE!==python set ParaMonte_BLD_DIR=!ParaMonte_BLD_DIR!\Python
if !INTERFACE_LANGUAGE!==fortran set ParaMonte_BLD_DIR=!ParaMonte_BLD_DIR!\Fortran

if exist !ParaMonte_BLD_DIR! (
    echo. -- ParaMonte - !ParaMonte_BLD_DIR! already exists. skipping...
) else (
    echo. -- ParaMonte - generating !BTYPE! library directory: !ParaMonte_BLD_DIR!
    mkdir !ParaMonte_BLD_DIR!
)
echo. -- ParaMonte - all generated build files will be stored at:
echo. -- ParaMonte -     %~dp0!ParaMonte_BLD_DIR!

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

if not !ParaMonte_OBJ_ENABLED!==true goto LABEL_ParaMonte_LIB_ENABLED

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
        /module:"!ParaMonte_MOD_DIR!"   %=path to output ParaMonte example module files=% ^
        /I:"!ParaMonte_MOD_DIR!"        %=path to input ParaMonte module files=%  ^
        /I:"!MATLAB_INC_DIR!"           %=path to the MATLAB include files=% ^
        /c "!ParaMonte_SRC_DIR!\%%A"    %=path to input ParaMonte example source files=% ^
        || (
            REM if not !ERRORLEVEL!==0 (
            echo.
            echo. -- ParaMonte - Fatal Error: compilation of the object file for %%A failed.
            echo. -- ParaMonte - build failed. exiting...
            echo.
            cd %~dp0
            set ERRORLEVEL=1
            exit /B 1
        )
    )
    if !INTERFACE_LANGUAGE!==matlab (
        echo. -- ParaMonte - generating object file for !MATLAB_VERSION_FILE!
        call !FCL! !FCL_FLAGS! !FPP_FLAGS! !FC_LIB_FLAGS! ^
        /module:"!ParaMonte_MOD_DIR!"   %=path to output ParaMonte example module files=% ^
        /I:"!ParaMonte_MOD_DIR!"        %=path to input ParaMonte module files=%  ^
        /I:"!MATLAB_INC_DIR!"           %=path to the MATLAB include files=% ^
        /c "!MATLAB_VERSION_FILE!"      %=path to input MATLAB source file=% ^
        || (
            REM if not !ERRORLEVEL!==0 (
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
        /module:"!ParaMonte_MOD_DIR!"   %=path to output ParaMonte example module files=% ^
        /I:"!ParaMonte_MOD_DIR!"        %=path to input ParaMonte module files=%  ^
        /I:"!MATLAB_INC_DIR!"           %=path to the MATLAB include files=% ^
        /c "!ParaMonte_SRC_DIR!\%%A"    %=path to input ParaMonte example source files=% ^
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
    if !INTERFACE_LANGUAGE!==matlab (
        echo. -- ParaMonte - generating object file for !MATLAB_VERSION_FILE!
        !FCL! !FCL_FLAGS! !FPP_FLAGS! !FC_LIB_FLAGS! ^
        /module:"!ParaMonte_MOD_DIR!"   %=path to output ParaMonte example module files=% ^
        /I:"!ParaMonte_MOD_DIR!"        %=path to input ParaMonte module files=%  ^
        /I:"!MATLAB_INC_DIR!"           %=path to the MATLAB include files=% ^
        /c "!MATLAB_VERSION_FILE!"      %=path to input MATLAB source file=% ^
        || (
            REM if not !ERRORLEVEL!==0 (
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

REM set PMLIB_NAME=libparamonte_!LTYPE!_!MEMORY_ALLOCATION!_!BTYPE!_!COMPILER_SUITE!_!INTERFACE_LANGUAGE!
REM if !CAF_ENABLED!==true set PMLIB_NAME=!PMLIB_NAME!_caf!CAFTYPE!
REM if !MPI_ENABLED!==true set PMLIB_NAME=!PMLIB_NAME!_mpi
REM if !OMP_ENABLED!==true set PMLIB_NAME=!PMLIB_NAME!_omp
REM set PMLIB_NAME=!PMLIB_NAME!_windows_!PLATFORM!
REM if defined MULTITHREADING set PMLIB_NAME=!PMLIB_NAME!_!MULTITHREADING!

if !INTERFACE_LANGUAGE!==c++ (
    set LANG_ABBR=cpp
) else (
    set LANG_ABBR=!INTERFACE_LANGUAGE!
)

set PMLIB_NAME=libparamonte_!LANG_ABBR!_windows_!PLATFORM!_!COMPILER_SUITE!_!BTYPE!_!LTYPE!_!MEMORY_ALLOCATION!
if !CAF_ENABLED!==true set PMLIB_NAME=!PMLIB_NAME!_caf!CAFTYPE!
if !MPI_ENABLED!==true set PMLIB_NAME=!PMLIB_NAME!_mpi
if !OMP_ENABLED!==true set PMLIB_NAME=!PMLIB_NAME!_omp

:: if not !ParaMonte_LIB_ENABLED!==true goto LABEL_ParaMonteTest_OBJ_ENABLED
if not !ParaMonte_LIB_ENABLED!==true goto LABEL_EOF

REM if not !LTYPE!==dynamic goto LABEL_EOF

cd !ParaMonte_LIB_DIR! && del *.exp *.lib *.dll *.pdb *.obj *.mexw64 *.optrpt
echo.
echo. -- ParaMonte - building ParaMonte !LTYPE! library files at: !ParaMonte_LIB_DIR!
echo. -- ParaMonte - !LTYPE! library filename: !PMLIB_NAME!
echo.

if !LTYPE!==dynamic (

    !FCL! !FCL_FLAGS! !FL_FLAGS! !FL_LIB_FLAGS! ^
    /module:!ParaMonte_MOD_DIR!                         %=path to output ParaMonte example module files=% ^
    /I:!ParaMonte_MOD_DIR!                              %=path to input ParaMonte module files=% ^
    "!ParaMonte_OBJ_DIR!"\*.obj                         %=path to the ParaMonte object files=% ^
    "!MATLAB_LIBMX_FILE!"                               %=path to the MATLAB libmx=% ^
    "!MATLAB_LIBMEX_FILE!"                              %=path to the MATLAB libmex=% ^
    "!MATLAB_LIBMAT_FILE!"                              %=path to the MATLAB libmmat=% ^
    /exe:!PMLIB_NAME!.dll

) else (

    REM https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-creating-static-libraries

    if !BTYPE!==release (

        !FCL! !FCL_FLAGS! !FL_FLAGS! !FL_LIB_FLAGS! ^
        /module:"!ParaMonte_MOD_DIR!"                               %=path to output ParaMonte example module files=% ^
        /I:"!ParaMonte_MOD_DIR!"                                    %=path to input ParaMonte module files=% ^
        "!ParaMonte_OBJ_DIR!"\*.obj                                 %=path to input ParaMonte object files=%

        xilib /out:!PMLIB_NAME!.lib ipo_out.obj
        REM del !ParaMonte_OBJ_DIR!\ipo_out.obj

    ) else (

        xilib /out:!PMLIB_NAME!.lib "!ParaMonte_OBJ_DIR!"\*.obj

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
