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

:: This batch file UNconfigures the flags required for building ParaMonte library, tests, and examples on Windows-Operating Systems.
:: It is used only upon exit from the build.bat script to remove the locally-set Windows CMD environmental configuration flags.

REM cd %~dp0

set ERRORLEVEL=1

set               ParaMonte_FLAG_CLEANUP_ENABLED=
set                           ParaMonte_ROOT_DIR=
set                        ParaMonte_OBJ_ENABLED=
set                        ParaMonte_LIB_ENABLED=
set                    ParaMonteTest_OBJ_ENABLED=
set                    ParaMonteTest_EXE_ENABLED=
set                    ParaMonteTest_RUN_ENABLED=
set                 ParaMonteExample_EXE_ENABLED=
set                 ParaMonteExample_RUN_ENABLED=
set                           INTERFACE_LANGUAGE=
set                                        BTYPE=
set                                        LTYPE=
set                           HEAP_ARRAY_ENABLED=
set                                  CFI_ENABLED=
set                                      CAFTYPE=
set                       FOR_COARRAY_NUM_IMAGES=
set                                  MPI_ENABLED=
set                                  OMP_ENABLED=
set                               COMPILER_SUITE=
set                             COMPILER_VERSION=
set                    INTEL_FORTRAN_DEBUG_FLAGS=
set                  INTEL_FORTRAN_RELEASE_FLAGS=
set                  INTEL_FORTRAN_TESTING_FLAGS=
set            INTEL_FORTRAN_PREPROCESSOR_MACROS=
set                        INTEL_CPP_DEBUG_FLAGS=
set                      INTEL_CPP_RELEASE_FLAGS=
set                      INTEL_CPP_TESTING_FLAGS=
set                                  Python_PATH=
set                                    FILE_LIST=

set ERRORLEVEL=0

exit /B 0