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

:: This batch file UNconfigures the flags required for building ParaMonte library, tests, and examples on Windows-Operating Systems.
:: It is used only upon exit from the build.bat script to remove the locally-set Windows CMD environmental configuration flags.

REM cd %~dp0

set ERRORLEVEL=1

set               ParaMonte_FLAG_CLEANUP_ENABLED=
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