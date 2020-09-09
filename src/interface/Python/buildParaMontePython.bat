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

:: set ParaMonte test source directory

set ParaMontePythonTest_SRC_DIR=!ParaMonteInterfacePython_SRC_DIR!\test

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

:: The ParaMonte library auxil files

echo. -- ParaMontePython - copying the ParaMonte library auxiliary files
echo. -- ParaMontePython - from: !ParaMonteInterface_SRC_DIR!\auxil                 %= no need for final slash here =%
echo. -- ParaMontePython -   to: !ParaMontePythonTest_BLD_DIR!\paramonte\auxil\     %= final slash tells this is folder =%
xcopy /s /Y "!ParaMonteInterface_SRC_DIR!\auxil" "!ParaMontePythonTest_BLD_DIR!\paramonte\auxil\" || goto LABEL_copyErrorOccured
echo.

:: copy necessary ParaMonte Python library files in Python's directory

echo. -- ParaMontePython - copying ParaMonte library source files to the Python build directory
echo. -- ParaMontePython - from: !ParaMonteInterfacePython_SRC_DIR!\paramonte   %= no need for final slash here =%
echo. -- ParaMontePython -   to: !ParaMontePythonTest_BLD_DIR!\paramonte\       %= final slash tells this is folder =%
xcopy /s /Y "!ParaMonteInterfacePython_SRC_DIR!\paramonte" "!ParaMontePythonTest_BLD_DIR!\paramonte\" || goto LABEL_copyErrorOccured
echo.

:: copy necessary ParaMonte Python DLL files in Python's directory

echo. -- ParaMontePython - copying ParaMonte DLL files to the Python build directory
echo. -- ParaMontePython - from: !ParaMonte_LIB_DIR!\!PMLIB_NAME!.* %= no need for final slash here =%
echo. -- ParaMontePython -   to: !ParaMontePythonTest_BLD_DIR!\paramonte\   %= final slash tells this is folder =%
xcopy /s /Y "!ParaMonte_LIB_DIR!\!PMLIB_NAME!.*" "!ParaMontePythonTest_BLD_DIR!\paramonte\" || goto LABEL_copyErrorOccured
echo.

:: copy necessary ParaMonte Python library files in Python's directory

echo. -- ParaMontePython - copying the ParaMonte library test files to the Python build directory
echo. -- ParaMontePython - from: !ParaMontePythonTest_SRC_DIR!\!PYTHON_TEST_FILENAME!  %= no need for final slash here =%
echo. -- ParaMontePython -   to: !ParaMontePythonTest_BLD_DIR!\                        %= final slash tells this is folder =%
xcopy /s /Y "!ParaMontePythonTest_SRC_DIR!\!PYTHON_TEST_FILENAME!" "!ParaMontePythonTest_BLD_DIR!\" || goto LABEL_copyErrorOccured
echo.

:: copy necessary input files in Python's directory

echo. -- ParaMontePython - copying the test input files to the ParaMonte Python build directory
echo. -- ParaMontePython - from: !ParaMonteTest_SRC_DIR!\input          %= no need for final slash here =%
echo. -- ParaMontePython -   to: !ParaMontePythonTest_BLD_DIR!\input\   %= final slash tells this is folder =%
xcopy /s /Y "!ParaMonteTest_SRC_DIR!\input" "!ParaMontePythonTest_BLD_DIR!\input\" || goto LABEL_copyErrorOccured
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
