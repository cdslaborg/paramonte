::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::
::   ParaMonte: plain powerful parallel Monte Carlo library.
::
::   Copyright (C) 2012-present, The Computational Data Science Lab
::
::   This file is part of the ParaMonte library.
::
::   ParaMonte is free software: you can redistribute it and/or modify it
::   under the terms of the GNU Lesser General Public License as published
::   by the Free Software Foundation, version 3 of the License.
::
::   ParaMonte is distributed in the hope that it will be useful,
::   but WITHOUT ANY WARRANTY; without even the implied warranty of
::   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
::   GNU Lesser General Public License for more details.
::
::   You should have received a copy of the GNU Lesser General Public License
::   along with the ParaMonte library. If not, see,
::
::       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
::
::   ACKNOWLEDGMENT
::
::   As per the ParaMonte library license agreement terms,
::   if you use any parts of this library for any purposes,
::   we ask you to acknowledge the use of the ParaMonte library
::   in your work (education/research/industry/development/...)
::   by citing the ParaMonte library as described on this page:
::
::       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: make bin directory
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set ParaMonte_BIN_DIR=!ParaMonte_BLD_ROOT_DIR!\bin
echo. -- ParaMonte - The ParaMonte binaries directory: !ParaMonte_BIN_DIR!
if not exist !ParaMonte_BIN_DIR! (
    mkdir "!ParaMonte_BIN_DIR!\"
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: setup examples' interface language
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set LANG_NAME=
set LANG_IS_C=false
set LANG_IS_MATLAB=false
set LANG_IS_Python=false
set LANG_IS_Fortran=false
set LANG_IS_DYNAMIC=false
set LANG_IS_COMPILED=false

if !INTERFACE_LANGUAGE!==matlab (
    set LANG_IS_DYNAMIC=true
    set LANG_IS_MATLAB=true
    set LANG_FILE_EXT=m
    set LANG_NAME=MATLAB
    set ParaMonteMATLAB_BIN_ROOT_DIR=!ParaMonte_BIN_DIR!\MATLAB
    echo. -- ParaMonte - The ParaMonte MATLAB binaries root directory: !ParaMonteMATLAB_BIN_ROOT_DIR!
    if not exist !ParaMonteMATLAB_BIN_ROOT_DIR! (
        mkdir "!ParaMonteMATLAB_BIN_ROOT_DIR!\"
    )
)

if !INTERFACE_LANGUAGE!==python (
    set LANG_IS_DYNAMIC=true
    set LANG_IS_Python=true
    set LANG_FILE_EXT=py
    set LANG_NAME=Python
    set ParaMontePython_BIN_ROOT_DIR=!ParaMonte_BIN_DIR!\Python
    echo. -- ParaMonte - The ParaMonte Python binaries root directory: !ParaMontePython_BIN_ROOT_DIR!
    if not exist !ParaMontePython_BIN_ROOT_DIR! (
        mkdir "!ParaMontePython_BIN_ROOT_DIR!\"
    )
)

if !INTERFACE_LANGUAGE!==fortran (
    set LANG_IS_COMPILED=true
    set LANG_IS_Fortran=true
    set LANG_FILE_EXT=f90
    set LANG_NAME=Fortran
)

if !INTERFACE_LANGUAGE!==c (
    set LANG_IS_COMPILED=true
    set LANG_FILE_EXT=c
    set LANG_IS_C=true
    set LANG_NAME=C
)

if not defined LANG_NAME (
    echo. 
    echo. -- ParaMonteExample - Fatal Error: unrecognized or no language specified. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: build examples
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:: select examples to build

set EXAM_LIST=mvn

:: set and make example directories

set ParaMonteExample_BLD_DIR=!ParaMonte_BLD_DIR!\example
set ParaMonteInterface_SRC_DIR_CURRENT=!ParaMonteInterface_SRC_DIR!\!LANG_NAME!

echo. 
echo. -- ParaMonteExample!LANG_NAME! - generating the ParaMonte library examples in !LANG_NAME! language...
echo. -- ParaMonteExample!LANG_NAME! - The ParaMonte !LANG_NAME! examples directory: !ParaMonteExample_BLD_DIR!

for %%e in (!EXAM_LIST!) do ( 

    set EXAM_NAME=%%e

    set ParaMonteExample_BLD_DIR_CURRENT=!ParaMonteExample_BLD_DIR!\!EXAM_NAME!
    echo. -- ParaMonteExample!LANG_NAME! - The ParaMonte library !EXAM_NAME! example directory: !ParaMonteExample_BLD_DIR_CURRENT!
    if exist !ParaMonteExample_BLD_DIR_CURRENT! (
        echo. -- ParaMonteExample!LANG_NAME! - previous example build detected. deleting the old contents...
        rmdir /S /Q !ParaMonteExample_BLD_DIR_CURRENT! || goto LABEL_rmdirErrorOccured
        REM rd /S /Q !ParaMonteExample_BLD_DIR_CURRENT!
        REM /S  Removes all directories and files in the specified directory in addition to the directory itself. Used to remove a directory tree.
        REM /Q  Quiet mode, do not ask if ok to remove a directory tree with /S
        echo. -- ParaMonteExample!LANG_NAME! - regenerating the ParaMonte library !EXAM_NAME! example directory: !ParaMonteExample_BLD_DIR_CURRENT!
    )
    mkdir "!ParaMonteExample_BLD_DIR_CURRENT!\"

    REM The ParaMonte library kernel files

    set ParaMonteExample_LIB_DIR_CURRENT=!ParaMonteExample_BLD_DIR_CURRENT!
    if !LANG_IS_Python!==true set ParaMonteExample_LIB_DIR_CURRENT=!ParaMonteExample_LIB_DIR_CURRENT!\paramonte
    if !LANG_IS_MATLAB!==true set ParaMonteExample_LIB_DIR_CURRENT=!ParaMonteExample_LIB_DIR_CURRENT!\paramonte\lib

    echo. -- ParaMonteExample!LANG_NAME! - copying the ParaMonte library files...
    echo. -- ParaMonteExample!LANG_NAME! - from: !ParaMonte_LIB_DIR!                %= no need for final slash here =%
    echo. -- ParaMonteExample!LANG_NAME! -   to: !ParaMonteExample_LIB_DIR_CURRENT! %= final slash tells this is folder =%
    xcopy /s /Y "!ParaMonte_LIB_DIR!\libparamonte*.*" "!ParaMonteExample_LIB_DIR_CURRENT!\" || goto LABEL_copyErrorOccured
    xcopy /s /Y "!ParaMonte_LIB_DIR!\libparamonte*.*" "!ParaMonteExample_LIB_DIR_CURRENT!\" || goto LABEL_copyErrorOccured
    xcopy /s /Y "!ParaMonte_LIB_DIR!\libparamonte*.*" "!ParaMonteExample_LIB_DIR_CURRENT!\" || goto LABEL_copyErrorOccured
    xcopy /s /Y "!ParaMonte_LIB_DIR!\libparamonte*.*" "!ParaMonteExample_LIB_DIR_CURRENT!\" || goto LABEL_copyErrorOccured

    REM The ParaMonte library example required files

    echo. -- ParaMonteExample!LANG_NAME! - copying the ParaMonte library !EXAM_NAME! example required files in !LANG_NAME! language...

    if !LANG_IS_COMPILED!==true (

        echo. -- ParaMonteExample!LANG_NAME! - from: !ParaMonteExample_SRC_DIR!         %= no need for final slash here =%
        echo. -- ParaMonteExample!LANG_NAME! -   to: !ParaMonteExample_BLD_DIR_CURRENT! %= final slash tells this is folder =%
        xcopy /s /Y "!ParaMonteExample_SRC_DIR!\build.bat" "!ParaMonteExample_BLD_DIR_CURRENT!\" || goto LABEL_copyErrorOccured

        REM The ParaMonte library example header/module files

        echo. -- ParaMonteExample!LANG_NAME! - copying the ParaMonte library interface files...
        echo. -- ParaMonteExample!LANG_NAME! - from: !ParaMonteInterface_SRC_DIR_CURRENT!   %= no need for final slash here =%
        echo. -- ParaMonteExample!LANG_NAME! -   to: !ParaMonteExample_BLD_DIR_CURRENT!\    %= final slash tells this is folder =%
        xcopy /s /Y "!ParaMonteInterface_SRC_DIR_CURRENT!" "!ParaMonteExample_BLD_DIR_CURRENT!\" || goto LABEL_copyErrorOccured

        if !LANG_IS_Fortran!==true (

            echo. -- ParaMonteExample!LANG_NAME! - copying the ParaMonte library Fortran module file compiled as paradram_mod.mod...
            echo. -- ParaMonteExample!LANG_NAME! - from: !ParaMonte_MOD_DIR! %= no need for final slash here =%
            echo. -- ParaMonteExample!LANG_NAME! -   to: !ParaMonteExample_BLD_DIR_CURRENT! %= final slash tells this is folder =%
            xcopy /s /Y "!ParaMonte_MOD_DIR!\paradram_mod.mod" "!ParaMonteExample_BLD_DIR_CURRENT!\" || goto LABEL_copyErrorOccured

        )

        REM if !LANG_IS_C!==true (
        REM )

    )

    if !LANG_IS_DYNAMIC!==true (

        REM The ParaMonte library README.md file

        echo. -- ParaMonteExample!LANG_NAME! - copying the ParaMonte library README.md file...
        echo. -- ParaMonteExample!LANG_NAME! - from: !ParaMonteInterface_SRC_DIR_CURRENT!\README.md
        echo. -- ParaMonteExample!LANG_NAME! -   to: !ParaMonteExample_BLD_DIR_CURRENT!\paramonte\README.md
        copy /y "!ParaMonteInterface_SRC_DIR_CURRENT!\README.md" "!ParaMonteExample_BLD_DIR_CURRENT!\README.md" || goto LABEL_copyErrorOccured

        REM The ParaMonte library CHANGES.md file

        echo. -- ParaMonteExample!LANG_NAME! - copying the ParaMonte library CHANGES.md file...
        echo. -- ParaMonteExample!LANG_NAME! - from: !ParaMonteInterface_SRC_DIR_CURRENT!\CHANGES.md
        echo. -- ParaMonteExample!LANG_NAME! -   to: !ParaMonteExample_BLD_DIR_CURRENT!\paramonte\CHANGES.md
        copy "!ParaMonteInterface_SRC_DIR_CURRENT!\CHANGES.md" "!ParaMonteExample_BLD_DIR_CURRENT!\CHANGES.md" || goto LABEL_copyErrorOccured

        REM The ParaMonte library license file

        echo. -- ParaMonteExample!LANG_NAME! - copying the ParaMonte library license file...
        echo. -- ParaMonteExample!LANG_NAME! - from: !ParaMonte_ROOT_DIR!\LICENSE
        echo. -- ParaMonteExample!LANG_NAME! -   to: !ParaMonteExample_BLD_DIR_CURRENT!\LICENSE
        copy "!ParaMonte_ROOT_DIR!\LICENSE" "!ParaMonteExample_BLD_DIR_CURRENT!\LICENSE" || goto LABEL_copyErrorOccured

        REM The ParaMonte library interface files

        echo. -- ParaMonteExample!LANG_NAME! - from: !ParaMonteInterface_SRC_DIR_CURRENT!\paramonte %= no need for final slash here =%
        echo. -- ParaMonteExample!LANG_NAME! -   to: !ParaMonteExample_BLD_DIR_CURRENT!\paramonte\  %= final slash tells this is folder =%
        xcopy /s /Y "!ParaMonteInterface_SRC_DIR_CURRENT!\paramonte" "!ParaMonteExample_BLD_DIR_CURRENT!\paramonte\" || goto LABEL_copyErrorOccured

        REM The ParaMonte library banner file

        echo. -- ParaMonteExample!LANG_NAME! - copying the ParaMonte library auxiliary files
        echo. -- ParaMonteExample!LANG_NAME! - from: !ParaMonteInterface_SRC_DIR!\auxil                     %= no need for final slash here =%
        echo. -- ParaMonteExample!LANG_NAME! -   to: !ParaMonteExample_BLD_DIR_CURRENT!\paramonte\auxil\    %= final slash tells this is folder =%
        xcopy /s /Y "!ParaMonteInterface_SRC_DIR!\auxil" "!ParaMonteExample_BLD_DIR_CURRENT!\paramonte\auxil\" || goto LABEL_copyErrorOccured
        echo.

        REM The ParaMonte library kernel version file (must appear only after the above)

        echo. -- ParaMonteExample!LANG_NAME! - copying the ParaMonte library kernel version file...
        echo. -- ParaMonteExample!LANG_NAME! - from: !ParaMonte_ROOT_DIR!\.VERSION %= no need for final slash here =%
        echo. -- ParaMonteExample!LANG_NAME! -   to: !ParaMonteExample_BLD_DIR_CURRENT!\paramonte\auxil\.VERSION_KERNEL
        copy /y "!ParaMonte_ROOT_DIR!\.VERSION" "!ParaMonteExample_BLD_DIR_CURRENT!\paramonte\auxil\.VERSION_KERNEL" || goto LABEL_copyErrorOccured

        REM The ParaMonte library interface version file (must appear only after the above)

        echo. -- ParaMonteExample!LANG_NAME! - copying the ParaMonte library interface version file...
        echo. -- ParaMonteExample!LANG_NAME! - from: !ParaMonteInterface_SRC_DIR_CURRENT!\.VERSION %= no need for final slash here =%
        echo. -- ParaMonteExample!LANG_NAME! -   to: !ParaMonteExample_BLD_DIR_CURRENT!\paramonte\auxil\.VERSION_INTERFACE
        copy /y "!ParaMonteInterface_SRC_DIR_CURRENT!\.VERSION" "!ParaMonteExample_BLD_DIR_CURRENT!\paramonte\auxil\.VERSION_INTERFACE" || goto LABEL_copyErrorOccured

        if !LANG_IS_Python!==true (

            REM The ParaMonte Python library setup files

            echo. -- ParaMonteExample!LANG_NAME! - copying the ParaMonte library Python setup files...
            echo. -- ParaMonteExample!LANG_NAME! - from: !ParaMonteInterface_SRC_DIR_CURRENT!\setup   %= no need for final slash here =%
            echo. -- ParaMonteExample!LANG_NAME! -   to: !ParaMontePython_BIN_ROOT_DIR!\            %= final slash tells this is folder =%
            xcopy /s /Y "!ParaMonteInterface_SRC_DIR_CURRENT!\setup" "!ParaMonteExample_BLD_DIR_CURRENT!\" || goto LABEL_copyErrorOccured

        )

    )

    REM The ParaMonte library example input files

    set ParaMonteExample_INP_DIR_CURRENT=!ParaMonteExample_SRC_DIR!\!EXAM_NAME!\input
    echo. -- ParaMonteExample!LANG_NAME! - copying the ParaMonte library !EXAM_NAME! example input files in !LANG_NAME! language...
    echo. -- ParaMonteExample!LANG_NAME! - from: !ParaMonteExample_INP_DIR_CURRENT!     %= no need for final slash here =%
    echo. -- ParaMonteExample!LANG_NAME! -   to: !ParaMonteExample_BLD_DIR_CURRENT!\    %= final slash tells this is folder =%
    xcopy /s /Y "!ParaMonteExample_INP_DIR_CURRENT!" "!ParaMonteExample_BLD_DIR_CURRENT!\" || goto LABEL_copyErrorOccured

    REM The ParaMonte library example source files

    set ParaMonteExample_SRC_DIR_CURRENT=!ParaMonteExample_SRC_DIR!\!EXAM_NAME!\!LANG_NAME!
    echo. -- ParaMonteExample!LANG_NAME! - copying the ParaMonte library !EXAM_NAME! example source files in !LANG_NAME! language...

    set mainFileName=main.!LANG_FILE_EXT!
    if !LANG_IS_DYNAMIC!==true (
        if !MPI_ENABLED!==true set mainFileName=main_mpi.!LANG_FILE_EXT!
    )
    copy "!ParaMonteExample_SRC_DIR!\!mainFileName!" "!ParaMonteExample_BLD_DIR_CURRENT!\" || goto LABEL_copyErrorOccured

    echo. -- ParaMonteExample!LANG_NAME! - from: !ParaMonteExample_SRC_DIR_CURRENT! %= no need for final slash here =%
    echo. -- ParaMonteExample!LANG_NAME! -   to: !ParaMonteExample_BLD_DIR_CURRENT! %= final slash tells this is folder =%
    xcopy /s /Y "!ParaMonteExample_SRC_DIR_CURRENT!" "!ParaMonteExample_BLD_DIR_CURRENT!\" || goto LABEL_copyErrorOccured

)

echo. 

if %ERRORLEVEL%==1 (
    echo. 
    echo. -- ParaMonteExample!LANG_NAME! - Fatal Error: build failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
if %ERRORLEVEL%==0 (
    echo.
    echo.
    echo. -- ParaMonteExample!LANG_NAME! - The ParaMonte library example build successful. 
    echo.
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: copy the first example to the bin directory
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set ParaMonteExample_BLD_DIR_CURRENT=!ParaMonteExample_BLD_DIR!\mvn

set ParaMonteExample_BIN_DIR_CURRENT=!ParaMonte_BIN_DIR!\!PMLIB_NAME!
if !LANG_IS_DYNAMIC!==true set ParaMonteExample_BIN_DIR_CURRENT=!ParaMonte_BIN_DIR!\!LANG_NAME!

echo. -- ParaMonteExample!LANG_NAME! - The ParaMonte !LANG_NAME! library binary directory: !ParaMonteExample_BIN_DIR_CURRENT!

if not exist !ParaMonteExample_BIN_DIR_CURRENT! (
    REM echo. -- ParaMonteExample!LANG_NAME! - previous binary directory detected. deleting the old contents...
    REM rmdir /S /Q !ParaMonteExample_BLD_DIR_CURRENT! 
    REM rd /S /Q !ParaMonteExample_BLD_DIR_CURRENT!
    REM REM /S  Removes all directories and files in the specified directory in addition to the directory itself. Used to remove a directory tree.
    REM REM /Q  Quiet mode, do not ask if ok to remove a directory tree with /S
    REM echo. -- ParaMonteExample!LANG_NAME! - regenerating the ParaMonte library !EXAM_NAME! example directory: !ParaMonteExample_BLD_DIR_CURRENT!
    mkdir "!ParaMonteExample_BIN_DIR_CURRENT!\"
)

echo. -- ParaMonteExample!LANG_NAME! - copying the ParaMonte library files to the bin folder...
echo. -- ParaMonteExample!LANG_NAME! - from: !ParaMonteExample_BLD_DIR_CURRENT! %= no need for final slash here =%
echo. -- ParaMonteExample!LANG_NAME! -   to: !ParaMonteExample_BIN_DIR_CURRENT! %= final slash tells this is folder =%
REM /s: Specifies to include subdirectories. Excludes empty subdirectories
REM /e: Copies all subdirectories, even if they are empty
REM /i: specifies the destination is a folder (Otherwise it prompts you)
xcopy /s /Y /e /v /i "!ParaMonteExample_BLD_DIR_CURRENT!" "!ParaMonteExample_BIN_DIR_CURRENT!" || goto LABEL_copyErrorOccured

if %ERRORLEVEL%==1 (
    echo. 
    echo. -- ParaMonteExample!LANG_NAME! - Fatal Error: build failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
if %ERRORLEVEL%==0 (
    echo.
    echo.
    echo. -- ParaMonteExample!LANG_NAME! - The ParaMonte library example build successful. 
    echo.
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: copy the first example to the bin directory
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

for %%e in (!EXAM_LIST!) do ( 

    set EXAM_NAME=%%e
    echo. -- ParaMonteExample!LANG_NAME! - Building/running the ParaMonte library's !EXAM_NAME! example.

    REM The ParaMonte library example build and run if requested

    if !ParaMonteExample_RUN_ENABLED!==true (

        set ParaMonteExample_BLD_DIR_CURRENT=!ParaMonteExample_BLD_DIR!\!EXAM_NAME!
        cd !ParaMonteExample_BLD_DIR_CURRENT!
        if !LANG_IS_COMPILED!==true (
            call build.bat || (
                echo. 
                echo. -- ParaMonteExample!LANG_NAME! - Fatal Error: The ParaMonte library example build/run failed. exiting...
                echo. 
                cd %~dp0
                set ERRORLEVEL=1
                exit /B 1
            )
        )
        cd %~dp0
    )

)

if %ERRORLEVEL%==1 (
    echo. 
    echo. -- ParaMonteExample!LANG_NAME! - Fatal Error: The ParaMonte library example build/run failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
if %ERRORLEVEL%==0 (
    echo.
    echo.
    echo. -- ParaMonteExample!LANG_NAME! - The ParaMonte library example build/run appears to have succeeded.
    echo.
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: quit
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

cd %~dp0

exit /B 0

:LABEL_copyErrorOccured

echo. 
echo. -- ParaMonteExample!LANG_NAME! - Fatal Error: failed to copy contents. exiting...
echo. 
cd %~dp0
set ERRORLEVEL=1
exit /B 1

:LABEL_rmdirErrorOccured

echo. 
echo. -- ParaMonteExample!LANG_NAME! - Fatal Error: failed to delete old contents. exiting...
echo. 
cd %~dp0
set ERRORLEVEL=1
exit /B 1
