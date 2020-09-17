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

:: SETLOCAL EnableDelayedExpansion

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: check MATLAB's existence
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set ERRORLEVEL=0

REM echo. 
REM echo. -- ParaMonteMATLAB - searching for a MATLAB installation on your system...
REM 
REM set "MATLAB_ROOT_DIR="
REM set "MATLAB_EXE_PATH="
REM set "MATLAB_BIN_DIR="
REM 
REM set "INSTALL_LOC_LIST=C:\Program Files\MATLAB\/C:\Program Files (x86)\MATLAB\"
REM set MATLAB_VERSION_LIST=R2025b/R2025a/R2024b/R2024a/R2023b/R2023a/R2022b/R2022a/R2021b/R2021a/R2020b/R2020a/R2019b/R2019a/R2018b/R2018a/R2017b/R2017a
REM 
REM for %%D in ("!INSTALL_LOC_LIST:/=" "!") do (
REM     for %%V in ("!MATLAB_VERSION_LIST:/=" "!") do (
REM         set "MATLAB_ROOT_DIR_TEMP=%%~D%%~V"
REM         set "MATLAB_BIN_DIR_TEMP=!MATLAB_ROOT_DIR_TEMP!\bin"
REM         set "MATLAB_EXE_PATH_TEMP=!MATLAB_BIN_DIR_TEMP!\matlab.exe"
REM         if exist !MATLAB_EXE_PATH_TEMP! (
REM             set "MATLAB_ROOT_DIR=!MATLAB_ROOT_DIR_TEMP!"
REM             set "MATLAB_EXE_PATH=!MATLAB_EXE_PATH_TEMP!"
REM             set "MATLAB_BIN_DIR=!MATLAB_BIN_DIR_TEMP!"
REM             echo. -- ParaMonteMATLAB - MATLAB %%~V installation detected at: !MATLAB_EXE_PATH!
REM             echo. 
REM             goto :LABEL_continue
REM         )
REM     )
REM )
REM 
REM echo. -- ParaMonteMATLAB - WARNING: Exhausted all possible search paths for a MATLAB installation, but failed to find MATLAB.
REM echo. -- ParaMonteMATLAB - WARNING: The ParaMonte MATLAB kernel will not be functional without building the required DLL libraries.
REM echo. -- ParaMonteMATLAB - WARNING: Please add MATLAB to your environmental variable PATH and rerun the install script.
REM echo. -- ParaMonteMATLAB - WARNING: For example, on your current Windows command-line, try:
REM echo. -- ParaMonteMATLAB - WARNING: 
REM echo. -- ParaMonteMATLAB - WARNING:     set "PATH=PATH_TO_MATLAB_BIN_DIR;!PATH!
REM echo. -- ParaMonteMATLAB - WARNING: 
REM echo. -- ParaMonteMATLAB - WARNING: where PATH_TO_MATLAB_BIN_DIR must be replaced with path to the bin folder of the current 
REM echo. -- ParaMonteMATLAB - WARNING: installation of MATLAB on your system. Typical MATLAB bin installation path on a 64-bit Windows 
REM echo. -- ParaMonteMATLAB - WARNING: Operating Systems is a string like the following:
REM echo. -- ParaMonteMATLAB - WARNING: 
REM echo. -- ParaMonteMATLAB - WARNING:     C:\Program Files\MATLAB\2020a\bin\
REM echo. -- ParaMonteMATLAB - WARNING: 
REM echo. -- ParaMonteMATLAB - WARNING: where 2020a in the path points to the MATLAB 2020a version installation on the system. You can also 
REM echo. -- ParaMonteMATLAB - WARNING: find the installation location of MATLAB by typing the following command in your MATLAB session:
REM echo. -- ParaMonteMATLAB - WARNING: 
REM echo. -- ParaMonteMATLAB - WARNING:     matlabroot
REM echo. -- ParaMonteMATLAB - WARNING: 
REM echo. -- ParaMonteMATLAB - WARNING: skipping the ParaMonte MATLAB build...
REM 
REM REM cd %~dp0
REM REM set ERRORLEVEL=1
REM REM exit /B 1
REM 
REM :LABEL_continue

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: Build ParaMonte MATLAB
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set ParaMonteMATLAB_BLD_LIB_DIR=!ParaMonte_BLD_DIR!\lib
if defined MATLAB_BIN_DIR (

    call set PMLIB_MATLAB_NAME=!PMLIB_NAME:_matlab_=_!
    set "MATLAB_BUILD_FLAGS=/DDLL_ENABLED "

    REM /subsystem:windows 
    REM mex -setup:"C:\Program Files\MATLAB\R2019a\bin\win64\mexopts\intel_c_19_vs2017.xml" C
    REM if !BTYPE!==debug   set "MATLAB_BUILD_FLAGS=!MATLAB_BUILD_FLAGS!/Od /Z7"
    REM if !BTYPE!==testing set "MATLAB_BUILD_FLAGS=!MATLAB_BUILD_FLAGS!/O2"
    REM if !BTYPE!==release set "MATLAB_BUILD_FLAGS=!MATLAB_BUILD_FLAGS!/Od"

    REM if !BTYPE!==debug   set "MATLAB_BUILD_FLAGS=!MATLAB_BUILD_FLAGS!!INTEL_CPP_DEBUG_FLAGS!"
    REM if !BTYPE!==testing set "MATLAB_BUILD_FLAGS=!MATLAB_BUILD_FLAGS!!INTEL_CPP_TESTING_FLAGS!"
    REM if !BTYPE!==release set "MATLAB_BUILD_FLAGS=!MATLAB_BUILD_FLAGS!!INTEL_CPP_RELEASE_FLAGS!"

    set "MEX_FLAGS=-v -nojvm"
    if !BTYPE!==debug set "MEX_FLAGS=!MEX_FLAGS! -g"
    if !BTYPE!==release set "MEX_FLAGS=!MEX_FLAGS! -O"
    echo. -- ParaMonteMATLAB - generating the ParaMonte MATLAB dynamic library: !ParaMonteMATLAB_BLD_LIB_DIR!\!PMLIB_MATLAB_NAME!
    echo. -- ParaMonteMATLAB - compiler options: !MATLAB_BUILD_FLAGS!
    echo. -- ParaMonteMATLAB - compiler command: "!MATLAB_BIN_DIR!\mex.bat" !MEX_FLAGS! "!ParaMonte_SRC_DIR!\paramonte.m.c" !PMLIB_NAME!.lib -output !PMLIB_MATLAB_NAME!
    cd !ParaMonteMATLAB_BLD_LIB_DIR!
    if !FRESH_RUN!==true (
        REM CC=icl COMPFLAGS="!MATLAB_BUILD_FLAGS!"
        call "!MATLAB_BIN_DIR!\mex.bat" !MEX_FLAGS! "!ParaMonte_SRC_DIR!\paramonte.m.c" !PMLIB_NAME!.lib -output !PMLIB_MATLAB_NAME! && (
        REM if !ERRORLEVEL!==0 (
            echo.
            echo. -- ParaMonteMATLAB - the ParaMonte MATLAB dynamic library build appears to have succeeded.
            echo.
        ) || (
        REM ) else (
            echo. 
            echo. -- ParaMonteMATLAB - Fatal Error: The ParaMonte MATLAB library build failed.
            echo. -- ParaMonteMATLAB - Please make sure you have the following components installed
            echo. -- ParaMonteMATLAB - on your system before rerunning the installation script:
            echo. -- ParaMonteMATLAB - 
            echo. -- ParaMonteMATLAB -     -- MATLAB, including MATLAB compilers.
            echo. -- ParaMonteMATLAB -     -- Intel Parallel Studio icl/ifort compilers 2018 or newer.
            echo. -- ParaMonteMATLAB - 
            echo. -- ParaMonteMATLAB - Once you are sure of the existence of these components in your 
            echo. -- ParaMonteMATLAB - Windows command line environment, run the following command:
            echo. -- ParaMonteMATLAB - 
            echo. -- ParaMonteMATLAB -     "!MATLAB_BIN_DIR!\mex.bat" -setup C
            echo. -- ParaMonteMATLAB - 
            echo. -- ParaMonteMATLAB - Among the options displayed, you should see the command to setup
            echo. -- ParaMonteMATLAB - the Intel Parallel Studio icl compiler on your system.
            echo. -- ParaMonteMATLAB - This command should look similar to the following,
            echo. -- ParaMonteMATLAB - 
            echo. -- ParaMonteMATLAB -     "!MATLAB_BIN_DIR_TEMP!\mex.bat" -setup:"C:\Program Files\MATLAB\R2019a\bin\win64\mexopts\intel_c_19_vs2017.xml" C
            echo. -- ParaMonteMATLAB - 
            echo. -- ParaMonteMATLAB - with minor differences depending on your specific installations of 
            echo. -- ParaMonteMATLAB - 
            echo. -- ParaMonteMATLAB -     -- the Intel Parallel Studio version
            echo. -- ParaMonteMATLAB -     -- the Microsoft Visual Studio version
            echo. -- ParaMonteMATLAB -     -- the MATLAB version
            echo. -- ParaMonteMATLAB - 
            echo. -- ParaMonteMATLAB - Copy and paste this command in your terminal, run it, and then rerun the ParaMonte MATLAB installation script.
            echo. -- ParaMonteMATLAB - 
            echo. -- ParaMonteMATLAB - Please report this error at: 
            echo. -- ParaMonteMATLAB - 
            echo. -- ParaMonteMATLAB -     https://github.com/cdslaborg/paramonte/issues
            echo. -- ParaMonteMATLAB - 
            echo. -- ParaMonteMATLAB - gracefully exiting The ParaMonte build script.
            echo. 
            cd %~dp0
            set ERRORLEVEL=1
            exit /B 1
        )
    ) else (
        echo. -- ParaMonteMATLAB - dryrun mode enabled. skippting the MATLAB dynamic library build...
    )
    cd %~dp0
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set ParaMonte test source directory
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set ParaMonteMATLABTest_SRC_DIR=!ParaMonteInterfaceMATLAB_SRC_DIR!\test

:: generate MATLAB library directory

set MATLAB_TEST_FILENAME=testParaMonte_!BTYPE!
if !CAF_ENABLED!==true set MATLAB_TEST_FILENAME=!MATLAB_TEST_FILENAME!_!CAFTYPE!
if !MPI_ENABLED!==true set MATLAB_TEST_FILENAME=!MATLAB_TEST_FILENAME!_mpi
if !OMP_ENABLED!==true set MATLAB_TEST_FILENAME=!MATLAB_TEST_FILENAME!_omp
set MATLAB_TEST_FILENAME=!MATLAB_TEST_FILENAME!.m

set ParaMonteMATLABTest_BLD_DIR=!ParaMonte_BLD_DIR!\test\MATLAB
if exist !ParaMonteMATLABTest_BLD_DIR! (
    echo. -- ParaMonteMATLAB - !ParaMonteMATLABTest_BLD_DIR! already exists. skipping...
) else (
    echo. -- ParaMonteMATLAB - generating MATLAB files directory: !ParaMonteMATLABTest_BLD_DIR!
    mkdir !ParaMonteMATLABTest_BLD_DIR!
)
echo.

:: The ParaMonte library auxil files

echo. -- ParaMonteMATLAB - copying the ParaMonte library auxiliary files
echo. -- ParaMonteMATLAB - from: !ParaMonteInterface_SRC_DIR!\auxil                 %= no need for final slash here =%
echo. -- ParaMonteMATLAB -   to: !ParaMonteMATLABTest_BLD_DIR!\paramonte\auxil\     %= final slash tells this is folder =%
xcopy /s /Y "!ParaMonteInterface_SRC_DIR!\auxil" "!ParaMonteMATLABTest_BLD_DIR!\paramonte\auxil\" || goto LABEL_copyErrorOccured
echo.

:: copy necessary ParaMonte MATLAB library files in MATLAB's source directory

echo. -- ParaMonteMATLAB - copying the paramonte library source files to the MATLAB build directory
echo. -- ParaMonteMATLAB - from: !ParaMonteInterfaceMATLAB_SRC_DIR!\paramonte   %= no need for final slash here =%
echo. -- ParaMonteMATLAB -   to: !ParaMonteMATLABTest_BLD_DIR!\paramonte\       %= final slash tells this is folder =%
xcopy /s /Y "!ParaMonteInterfaceMATLAB_SRC_DIR!\paramonte" "!ParaMonteMATLABTest_BLD_DIR!\paramonte\" || goto LABEL_copyErrorOccured
echo.

:: copy necessary ParaMonte MATLAB DLL files in MATLAB's build directory

echo. -- ParaMonteMATLAB - copying the paramonte DLL files to the MATLAB build directory
echo. -- ParaMonteMATLAB - from: !ParaMonte_LIB_DIR!\!PMLIB_NAME!.* %= no need for final slash here =%
echo. -- ParaMonteMATLAB -   to: !ParaMonteMATLABTest_BLD_DIR!\paramonte\   %= final slash tells this is folder =%
xcopy /s /Y /i "!ParaMonte_LIB_DIR!\*.*" "!ParaMonteMATLABTest_BLD_DIR!\paramonte\lib" || goto LABEL_copyErrorOccured
echo.

:: copy necessary ParaMonte MATLAB library files in MATLAB's directory

echo. -- ParaMonteMATLAB - copying the paramonte library test files to the MATLAB build directory
echo. -- ParaMonteMATLAB - from: !ParaMonteMATLABTest_SRC_DIR!\!MATLAB_TEST_FILENAME!  %= no need for final slash here =%
echo. -- ParaMonteMATLAB -   to: !ParaMonteMATLABTest_BLD_DIR!\                        %= final slash tells this is folder =%
xcopy /s /Y "!ParaMonteMATLABTest_SRC_DIR!\!MATLAB_TEST_FILENAME!" "!ParaMonteMATLABTest_BLD_DIR!\" || goto LABEL_copyErrorOccured
echo.

:: copy necessary input files in MATLAB's directory

echo. -- ParaMonteMATLAB - copying the test input files to the paramonte MATLAB build directory
echo. -- ParaMonteMATLAB - from: !ParaMonteTest_SRC_DIR!\input   %= no need for final slash here =%
echo. -- ParaMonteMATLAB -   to: !ParaMonteMATLABTest_BLD_DIR!\input\  %= final slash tells this is folder =%
xcopy /s /Y "!ParaMonteTest_SRC_DIR!\input" "!ParaMonteMATLABTest_BLD_DIR!\input\" || goto LABEL_copyErrorOccured
echo.

cd %~dp0
exit /B 0

:LABEL_copyErrorOccured

echo. 
echo. -- ParaMonteExample!LANG_NAME! - Fatal Error: failed to copy contents. exiting...
echo. 
cd %~dp0
set ERRORLEVEL=1
exit /B 1
