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

setlocal EnableDelayedExpansion

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: build ParaMonte library example objects and executable
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

echo.
echo. ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
echo. ::::                                                                                                                            ::::
echo.                                                   ParaMonte Library Examples Build
echo. ::::                                                                                                                            ::::
echo. ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
echo.

:: set and make example directories

set ParaMonteExample_BLD_DIR=!ParaMonte_BLD_DIR!\example

:: select examples to build

set EXAM_LIST=mvn

:: select languages for which to build

set LANG_LIST=
if !CFI_ENABLED!==true (
    set LANG_LIST=!LANG_LIST!C;
    if !LTYPE!==dynamic (
        set LANG_LIST=!LANG_LIST!Python
    )
) else (
    set LANG_LIST=!LANG_LIST!Fortran
)

:: make bin directory

set ParaMonte_BIN_DIR=!ParaMonte_BLD_ROOT_DIR!\bin
echo. -- ParaMonte - ParaMonte binaries directory: !ParaMonte_BIN_DIR!
if not exist !ParaMonte_BIN_DIR! (
    mkdir "!ParaMonte_BIN_DIR!\"
)

set ParaMontePython_BIN_ROOT_DIR=!ParaMonte_BIN_DIR!\Python
REM set ParaMontePython_BIN_PMCORE_DIR=!ParaMontePython_BIN_ROOT_DIR!\pmcore
echo. -- ParaMonte - ParaMonte Python binaries root directory: !ParaMontePython_BIN_ROOT_DIR!
if not exist !ParaMontePython_BIN_ROOT_DIR! (
    mkdir "!ParaMontePython_BIN_ROOT_DIR!\"
)

:: build examples

for %%l in (!LANG_LIST!) do ( 

    set ParaMonteExample_LNG_DIR=!ParaMonteExample_BLD_DIR!\%%l

    echo. 
    echo. -- ParaMonteExample%%l - generating ParaMonte examples in %%l language...
    echo. -- ParaMonteExample%%l - ParaMonte %%l examples directory: !ParaMonteExample_LNG_DIR!

    set IS_Python_LANG=false
    if %%l==Python set IS_Python_LANG=true

    set IS_FortranC_LANG=false
    if %%l==C (
        set LANG_FILE_EXT=c
        set IS_FortranC_LANG=true
    )
    if %%l==Fortran (
        set LANG_FILE_EXT=f90
        set IS_FortranC_LANG=true
    )

    for %%e in (!EXAM_LIST!) do ( 

        REM set the binaries directory

        REM set ParaMonteExample_BIN_DIR_CURRENT=!ParaMonte_BIN_DIR!\!PMLIB_NAME!_!ParaMonteVersion!
        set ParaMonteExample_BIN_DIR_CURRENT=!ParaMonte_BIN_DIR!\!PMLIB_NAME!
        if !IS_Python_LANG!==true set ParaMonteExample_BIN_DIR_CURRENT=!ParaMontePython_BIN_ROOT_DIR!\paramonte
        echo. -- ParaMonte - The ParaMonte library binaries directory: !ParaMonteExample_BIN_DIR_CURRENT!

        set ParaMonteExample_BLD_DIR_CURRENT=!ParaMonteExample_LNG_DIR!\%%e
        if not exist !ParaMonteExample_BLD_DIR_CURRENT! (
            echo. -- ParaMonteExample%%l - The ParaMonte library %%e example directory: !ParaMonteExample_BLD_DIR_CURRENT!
            mkdir "!ParaMonteExample_BLD_DIR_CURRENT!\"
        )

        REM ParaMonte library files

        set ParaMonteExample_LIB_DIR_CURRENT=!ParaMonteExample_BLD_DIR_CURRENT!
        if !IS_Python_LANG!==true set ParaMonteExample_LIB_DIR_CURRENT=!ParaMonteExample_LIB_DIR_CURRENT!\paramonte
        echo. -- ParaMonteExample%%l - copying the ParaMonte library files...
        echo. -- ParaMonteExample%%l - from: !ParaMonte_LIB_DIR!                %= no need for final slash here =%
        echo. -- ParaMonteExample%%l -   to: !ParaMonteExample_LIB_DIR_CURRENT! %= final slash tells this is folder =%
        xcopy /s /Y "!ParaMonte_LIB_DIR!\!PMLIB_NAME!.*" "!ParaMonteExample_LIB_DIR_CURRENT!\"
        xcopy /s /Y "!ParaMonte_LIB_DIR!\!PMLIB_NAME!.*" "!ParaMonteExample_BIN_DIR_CURRENT!\"

        REM ParaMonte library example build files

        echo. -- ParaMonteExample%%l - copying the ParaMonte library %%e example build files in %%l language...

        if !IS_FortranC_LANG!==true (
            echo. -- ParaMonteExample%%l - from: !ParaMonteExample_SRC_DIR!         %= no need for final slash here =%
            echo. -- ParaMonteExample%%l -   to: !ParaMonteExample_BLD_DIR_CURRENT! %= final slash tells this is folder =%
            xcopy /s /Y "!ParaMonteExample_SRC_DIR!\build.bat" "!ParaMonteExample_BLD_DIR_CURRENT!\"
            xcopy /s /Y "!ParaMonteExample_SRC_DIR!\build.bat" "!ParaMonteExample_BIN_DIR_CURRENT!\"
        )

        if !IS_Python_LANG!==true (

            echo. -- ParaMonteExample%%l - from: !ParaMonteInterfacePython_SRC_DIR!\paramonte   %= no need for final slash here =%
            echo. -- ParaMonteExample%%l -   to: !ParaMonteExample_BLD_DIR_CURRENT!\paramonte\  %= final slash tells this is folder =%
            xcopy /s /Y "!ParaMonteInterfacePython_SRC_DIR!\paramonte" "!ParaMonteExample_BLD_DIR_CURRENT!\paramonte\"
            xcopy /s /Y "!ParaMonteInterfacePython_SRC_DIR!\paramonte" "!ParaMonteExample_BIN_DIR_CURRENT!\"

            REM PyPI build - The ParaMonte library Python setup files

            echo. -- ParaMonteExample%%l - copying the ParaMonte library Python setup files...
            echo. -- ParaMonteExample%%l - from: !ParaMonteInterfacePython_SRC_DIR!\setup   %= no need for final slash here =%
            echo. -- ParaMonteExample%%l -   to: !ParaMontePython_BIN_ROOT_DIR!\            %= final slash tells this is folder =%
            xcopy /s /Y "!ParaMonteInterfacePython_SRC_DIR!\setup" "!ParaMontePython_BIN_ROOT_DIR!\"

            REM PyPI build - The ParaMonte library license files

            echo. -- ParaMonteExample%%l - copying the ParaMonte library license file...
            echo. -- ParaMonteExample%%l - from: !ParaMonte_ROOT_DIR!\LICENSE %= no need for final slash here =%
            echo. -- ParaMonteExample%%l -   to: !ParaMontePython_BIN_ROOT_DIR!\LICENSE %= final slash tells this is folder =%
            copy "!ParaMonte_ROOT_DIR!\LICENSE" "!ParaMontePython_BIN_ROOT_DIR!\LICENSE"

            REM PyPI build - The ParaMonte library CHANGES.md files

            echo. -- ParaMonteExample%%l - copying the ParaMonte library CHANGES.md file...
            echo. -- ParaMonteExample%%l - from: !ParaMonte_ROOT_DIR!\CHANGES.md %= no need for final slash here =%
            echo. -- ParaMonteExample%%l -   to: !ParaMonteExample_BLD_DIR_CURRENT!\CHANGES.md
            REM copy "!ParaMonte_ROOT_DIR!\CHANGES.md" "!ParaMontePython_BIN_ROOT_DIR!\CHANGES.md"

        )

        REM ParaMonte library example input files

        set ParaMonteExample_INP_DIR_CURRENT=!ParaMonteExample_SRC_DIR!\%%e\input
        echo. -- ParaMonteExample%%l - copying the ParaMonte library %%e example input files in %%l language...
        echo. -- ParaMonteExample%%l - from: !ParaMonteExample_INP_DIR_CURRENT!     %= no need for final slash here =%
        echo. -- ParaMonteExample%%l -   to: !ParaMonteExample_BLD_DIR_CURRENT!\    %= final slash tells this is folder =%
        xcopy /s /Y "!ParaMonteExample_INP_DIR_CURRENT!" "!ParaMonteExample_BLD_DIR_CURRENT!\"
        if !IS_Python_LANG!==true (
            xcopy /s /Y "!ParaMonteExample_INP_DIR_CURRENT!" "!ParaMontePython_BIN_ROOT_DIR!\"
        )
        if !IS_FortranC_LANG!==true (
            xcopy /s /Y "!ParaMonteExample_INP_DIR_CURRENT!" "!ParaMonteExample_BIN_DIR_CURRENT!\"
        )

        REM ParaMonte library example source files

        set ParaMonteExample_SRC_DIR_CURRENT=!ParaMonteExample_SRC_DIR!\%%e\%%l
        echo. -- ParaMonteExample%%l - copying the ParaMonte library %%e example source files in %%l language...

        if !IS_FortranC_LANG!==true (

            echo. -- ParaMonteExample%%l - from: !ParaMonteExample_SRC_DIR_CURRENT! %= no need for final slash here =%
            echo. -- ParaMonteExample%%l -   to: !ParaMonteExample_BLD_DIR_CURRENT! %= final slash tells this is folder =%
            echo. -- ParaMonteExample%%l -   to: !ParaMonteExample_BIN_DIR_CURRENT! %= final slash tells this is folder =%
            xcopy /s /Y "!ParaMonteExample_SRC_DIR_CURRENT!" "!ParaMonteExample_BLD_DIR_CURRENT!\"
            xcopy /s /Y "!ParaMonteExample_SRC_DIR_CURRENT!" "!ParaMonteExample_BIN_DIR_CURRENT!\"

            echo. -- ParaMonteExample%%l - from: !ParaMonteExample_SRC_DIR!\main.!LANG_FILE_EXT! %= no need for final slash here =%
            echo. -- ParaMonteExample%%l -   to: !ParaMonteExample_BLD_DIR_CURRENT! %= final slash tells this is folder =%
            echo. -- ParaMonteExample%%l -   to: !ParaMonteExample_BIN_DIR_CURRENT! %= final slash tells this is folder =%
            xcopy /s /Y "!ParaMonteExample_SRC_DIR!\main.!LANG_FILE_EXT!" "!ParaMonteExample_BLD_DIR_CURRENT!\"
            xcopy /s /Y "!ParaMonteExample_SRC_DIR!\main.!LANG_FILE_EXT!" "!ParaMonteExample_BIN_DIR_CURRENT!\"

        )

        if !IS_Python_LANG!==true (

            echo. -- ParaMonteExample%%l - from: !ParaMonteExample_SRC_DIR_CURRENT!\README.md %= no need for final slash here =%
            echo. -- ParaMonteExample%%l -   to: !ParaMonteExample_BLD_DIR_CURRENT!\ %= final slash tells this is folder =%
            xcopy /s /Y "!ParaMonteExample_SRC_DIR_CURRENT!\README.md" "!ParaMonteExample_BLD_DIR_CURRENT!\"

            set PythonScriptFileName=main.py
            if !MPI_ENABLED!==true set PythonScriptFileName=main_mpi.py
            xcopy /s /Y "!ParaMonteExample_SRC_DIR_CURRENT!\!PythonScriptFileName!" "!ParaMonteExample_BLD_DIR_CURRENT!\"
            xcopy /s /Y "!ParaMonteExample_SRC_DIR_CURRENT!\!PythonScriptFileName!" "!ParaMontePython_BIN_ROOT_DIR!\"

            REM create pmcore, this must be the last step, only after the full creation of paramonte

            REM echo. -- ParaMonteExample%%l - generating ParaMonte pmcore Python library from the clone of paramonte...
            REM echo. -- ParaMonteExample%%l - from: !ParaMonteExample_INP_DIR_CURRENT!     %= no need for final slash here =%
            REM echo. -- ParaMonteExample%%l -   to: !ParaMonteExample_BLD_DIR_CURRENT!\    %= final slash tells this is folder =%
            REM xcopy /s /Y "!ParaMonteExample_BIN_DIR_CURRENT!" "!ParaMontePython_BIN_PMCORE_DIR!\"

        )

        REM ParaMonte library example header/module files

        if %%l==Fortran (
            echo. -- ParaMonteExample%%l - copying the ParaMonte library Fortran module file paradram_mod.mod...
            echo. -- ParaMonteExample%%l - from: !ParaMonte_MOD_DIR! %= no need for final slash here =%
            echo. -- ParaMonteExample%%l -   to: !ParaMonteExample_BLD_DIR_CURRENT! %= final slash tells this is folder =%
            echo. -- ParaMonteExample%%l -   to: !ParaMonteExample_BIN_DIR_CURRENT! %= final slash tells this is folder =%
            xcopy /s /Y "!ParaMonte_MOD_DIR!\paradram_mod.mod" "!ParaMonteExample_BLD_DIR_CURRENT!\"
            xcopy /s /Y "!ParaMonte_MOD_DIR!\paradram_mod.mod" "!ParaMonteExample_BIN_DIR_CURRENT!\"
            echo. -- ParaMonteExample%%l - copying the ParaMonte library Fortran module file paramonte.f90...
            echo. -- ParaMonteExample%%l - from: !ParaMonteInterfaceFortran_SRC_DIR! %= no need for final slash here =%
            echo. -- ParaMonteExample%%l -   to: !ParaMonteExample_BLD_DIR_CURRENT!\ %= final slash tells this is folder =%
            echo. -- ParaMonteExample%%l -   to: !ParaMonteExample_BIN_DIR_CURRENT!\ %= final slash tells this is folder =%
            xcopy /s /Y "!ParaMonteInterfaceFortran_SRC_DIR!" "!ParaMonteExample_BLD_DIR_CURRENT!\"
            xcopy /s /Y "!ParaMonteInterfaceFortran_SRC_DIR!" "!ParaMonteExample_BIN_DIR_CURRENT!\"
        )
        if %%l==C (
            echo. -- ParaMonteExample%%l - copying the ParaMonte library C header file paramonte.h...
            echo. -- ParaMonteExample%%l - from: !ParaMonteInterfaceC_SRC_DIR!          %= no need for final slash here =%
            echo. -- ParaMonteExample%%l -   to: !ParaMonteExample_BLD_DIR_CURRENT!\    %= final slash tells this is folder =%
            echo. -- ParaMonteExample%%l -   to: !ParaMonteExample_BIN_DIR_CURRENT!\    %= final slash tells this is folder =%
            xcopy /s /Y "!ParaMonteInterfaceC_SRC_DIR!" "!ParaMonteExample_BLD_DIR_CURRENT!\"
            xcopy /s /Y "!ParaMonteInterfaceC_SRC_DIR!" "!ParaMonteExample_BIN_DIR_CURRENT!\"
        )

        REM ParaMonte library example build and run if requested

        if !ParaMonteExample_RUN_ENABLED!==true (

            cd !ParaMonteExample_BLD_DIR_CURRENT!
            if !IS_FortranC_LANG!==true (
                call build.bat
            )
            cd %~dp0
        )

    )

    echo. 

)

if %ERRORLEVEL%==1 (
    echo. 
    echo. -- ParaMonteExample - Fatal Error: build failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
if %ERRORLEVEL%==0 (
    echo.
    echo.
    echo. -- ParaMonteExample - The ParaMonte library example build successful. 
    echo.
)

cd %~dp0

exit /B 0