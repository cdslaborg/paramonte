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

:: NOTE: This is not a standalone build script. It must only be called by buildParaMonte.bat script in the root directory of the project.

:: SETLOCAL EnableDelayedExpansion

:: set ParaMonte test source directory

set ParaMontePythonTest_SRC_DIR=!ParaMontePython_SRC_DIR!\test

:: generate Python library directory

set PYTHON_TEST_FILENAME=testParaMonte_!BTYPE!
if !CAF_ENABLED!==true set PYTHON_TEST_FILENAME=!PYTHON_TEST_FILENAME!_!CAFTYPE!
if !MPI_ENABLED!==true set PYTHON_TEST_FILENAME=!PYTHON_TEST_FILENAME!_mpi
if !OMP_ENABLED!==true set PYTHON_TEST_FILENAME=!PYTHON_TEST_FILENAME!_omp
set PYTHON_TEST_FILENAME=!PYTHON_TEST_FILENAME!.py

set ParaMontePythonTest_BLD_DIR=!ParaMonte_BLD_DIR!\test\Python
if exist !ParaMontePythonTest_BLD_DIR! (
    echo. -- ParaMontePython - !ParaMontePythonTest_BLD_DIR! already exists. skipping...
) else (
    echo. -- ParaMontePython - generating Python files directory: !ParaMontePythonTest_BLD_DIR!
    mkdir !ParaMontePythonTest_BLD_DIR!
)
echo.

:: copy necessary ParaMonte Python library files in Python's directory

echo. -- ParaMontePython - copying paramonte library files to the Python directory
echo. -- ParaMontePython - from: !ParaMontePython_SRC_DIR!\paramonte       %= no need for final slash here =%
echo. -- ParaMontePython -   to: !ParaMontePythonTest_BLD_DIR!\paramonte\  %= final slash tells this is folder =%
xcopy /s /Y "!ParaMontePython_SRC_DIR!\paramonte" "!ParaMontePythonTest_BLD_DIR!\paramonte\"
echo.

:: copy necessary ParaMonte Python DLL files in Python's directory

echo. -- ParaMontePython - copying paramonte DLL files to the Python directory
echo. -- ParaMontePython - from: !ParaMonte_LIB_DIR!\!PMLIB_NAME!.* %= no need for final slash here =%
echo. -- ParaMontePython -   to: !ParaMontePythonTest_BLD_DIR!\paramonte\   %= final slash tells this is folder =%
xcopy /s /Y "!ParaMonte_LIB_DIR!\!PMLIB_NAME!.*" "!ParaMontePythonTest_BLD_DIR!\paramonte\"
echo.

:: copy necessary ParaMonte Python library files in Python's directory

echo. -- ParaMontePython - copying paramonte library test files to the Python directory
echo. -- ParaMontePython - from: !ParaMontePythonTest_SRC_DIR!\!PYTHON_TEST_FILENAME!  %= no need for final slash here =%
echo. -- ParaMontePython -   to: !ParaMontePythonTest_BLD_DIR!\                        %= final slash tells this is folder =%
xcopy /s /Y "!ParaMontePythonTest_SRC_DIR!\!PYTHON_TEST_FILENAME!" "!ParaMontePythonTest_BLD_DIR!\"
echo.

:: copy necessary input files in Python's directory

echo. -- ParaMontePython - copying input files to the Python directory
echo. -- ParaMontePython - from: !ParaMonteTest_SRC_DIR!\input   %= no need for final slash here =%
echo. -- ParaMontePython -   to: !ParaMontePythonTest_BLD_DIR!\input\  %= final slash tells this is folder =%
xcopy /s /Y "!ParaMonteTest_SRC_DIR!\input" "!ParaMontePythonTest_BLD_DIR!\input\"
echo.
