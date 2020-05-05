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

:: SETLOCAL EnableDelayedExpansion

:: set ParaMonte test source directory

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

:: copy necessary ParaMonte MATLAB library files in MATLAB's directory

echo. -- ParaMonteMATLAB - copying the paramonte library source files to the MATLAB build directory
echo. -- ParaMonteMATLAB - from: !ParaMonteInterfaceMATLAB_SRC_DIR!\paramonte   %= no need for final slash here =%
echo. -- ParaMonteMATLAB -   to: !ParaMonteMATLABTest_BLD_DIR!\paramonte\       %= final slash tells this is folder =%
xcopy /s /Y "!ParaMonteInterfaceMATLAB_SRC_DIR!\paramonte" "!ParaMonteMATLABTest_BLD_DIR!\paramonte\"
echo.

:: copy necessary ParaMonte MATLAB DLL files in MATLAB's directory

echo. -- ParaMonteMATLAB - copying the paramonte DLL files to the MATLAB build directory
echo. -- ParaMonteMATLAB - from: !ParaMonte_LIB_DIR!\!PMLIB_NAME!.* %= no need for final slash here =%
echo. -- ParaMonteMATLAB -   to: !ParaMonteMATLABTest_BLD_DIR!\paramonte\   %= final slash tells this is folder =%
xcopy /s /Y "!ParaMonte_LIB_DIR!\!PMLIB_NAME!.*" "!ParaMonteMATLABTest_BLD_DIR!\paramonte\"
echo.

:: copy necessary ParaMonte MATLAB library files in MATLAB's directory

echo. -- ParaMonteMATLAB - copying the paramonte library test files to the MATLAB build directory
echo. -- ParaMonteMATLAB - from: !ParaMonteMATLABTest_SRC_DIR!\!MATLAB_TEST_FILENAME!  %= no need for final slash here =%
echo. -- ParaMonteMATLAB -   to: !ParaMonteMATLABTest_BLD_DIR!\                        %= final slash tells this is folder =%
xcopy /s /Y "!ParaMonteMATLABTest_SRC_DIR!\!MATLAB_TEST_FILENAME!" "!ParaMonteMATLABTest_BLD_DIR!\"
echo.

:: copy necessary input files in MATLAB's directory

echo. -- ParaMonteMATLAB - copying the test input files to the paramonte MATLAB build directory
echo. -- ParaMonteMATLAB - from: !ParaMonteTest_SRC_DIR!\input   %= no need for final slash here =%
echo. -- ParaMonteMATLAB -   to: !ParaMonteMATLABTest_BLD_DIR!\input\  %= final slash tells this is folder =%
xcopy /s /Y "!ParaMonteTest_SRC_DIR!\input" "!ParaMonteMATLABTest_BLD_DIR!\input\"
echo.
