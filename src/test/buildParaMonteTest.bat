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

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: build ParaMonte library test objects and executable
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

cd %~dp0

echo.
echo. ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
echo. ::::                                                                                                                            ::::
echo.                                                     ParaMonte Library Test Build
echo. ::::                                                                                                                            ::::
echo. ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
echo.

:: copy necessary lib files in the executable's directory

if !LTYPE!==dynamic (

    echo. 
    echo. -- ParaMonteTest - testing ParaMonte library in dynamic build mode is currently unsupported...
    echo. -- ParaMonteTest - skipping ...
    echo. 
    cd %~dp0
    set ERRORLEVEL=0
    exit /B 0

    REM echo. -- ParaMonteTest - copying library files to the bin directory
    REM echo. -- ParaMonteTest - from: !ParaMonte_LIB_DIR!\   %= no need for final slash here =%
    REM echo. -- ParaMonteTest -   to: !ParaMonteTest_BIN_DIR!\  %= final slash tells this is folder =%
    REM xcopy /s /Y !ParaMonte_LIB_DIR!\!PMLIB_NAME!.* !ParaMonteTest_BIN_DIR!\
    REM echo.

)

:: set and make test directories

set ParaMonteTest_BLD_DIR=!ParaMonte_BLD_DIR!\test
set ParaMonteTest_BIN_DIR=!ParaMonteTest_BLD_DIR!\bin
set ParaMonteTest_MOD_DIR=!ParaMonteTest_BLD_DIR!\mod
set ParaMonteTest_OBJ_DIR=!ParaMonteTest_BLD_DIR!\obj

:: loop over test directories and generate them
echo.
for %%A in (
    !ParaMonteTest_BLD_DIR!
    !ParaMonteTest_BIN_DIR!
    !ParaMonteTest_MOD_DIR!
    !ParaMonteTest_OBJ_DIR!
    ) do (  if exist %%A (
                echo. -- ParaMonteTest - %%A already exists. skipping...
            ) else (
                echo. -- ParaMonteTest - generating test directory: %%A
                mkdir %%A
            )
)
echo.

if !ParaMonteTest_OBJ_ENABLED! NEQ true (
    echo.
    echo. -- ParaMonteTest - Warning: skipping ParaMonte library test object files build...
    echo.
    goto LABEL_ParaMonteTest_EXE_ENABLED
)

:: Read the name of each file from the ordered list of filenames in filelist.txt to compile

cd !ParaMonteTest_OBJ_DIR!
echo.
echo. -- ParaMonteTest - building ParaMonte tests...

:: First verify the source filelist exists

set ParaMonteTest_FILE_LIST=!ParaMonteTest_SRC_DIR!\filelist.txt
if not exist !ParaMonteTest_FILE_LIST! (
    echo.
    echo. -- ParaMonteTest - Fatal Error: The filelist.txt containing the ParaMonte Test source filenames does not exist. Path: !ParaMonteTest_FILE_LIST!
    echo. -- ParaMonteTest - build failed. exiting...
    echo.
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)

:: generate object files

if !MPI_ENABLED!==true (
    echo. -- ParaMonteTest - building test as MPI application...
    echo.
    for /F "eol=! tokens=*" %%A in (!ParaMonteTest_FILE_LIST!) do (
        echo. -- ParaMonteTest - generating object file for %%A
        call !FCL! !FCL_FLAGS! !FPP_FLAGS! ^
        /module:!ParaMonteTest_MOD_DIR!     %=path to output ParaMonte Test module files=% ^
        /I:!ParaMonteTest_MOD_DIR!          %=path to required ParaMonte Test module files=%  ^
        /I:!ParaMonte_MOD_DIR!              %=path to input ParaMonte module files=%  ^
        /c !ParaMonteTest_SRC_DIR!\%%A      %=path to input ParaMonte Test source file=%
        if !ERRORLEVEL! NEQ 0 (
            echo. 
            echo. -- ParaMonteTest - Fatal Error: compilation of the object file for %%A failed.
            echo. -- ParaMonteTest - build failed. exiting...
            echo. 
            set ERRORLEVEL=1
            exit /B
        )
    )
) else (
    for /F "eol=! tokens=*" %%A in (!ParaMonteTest_FILE_LIST!) do (

        echo. -- ParaMonteTest - generating object file for %%A

        !FCL! !FCL_FLAGS! !FPP_FLAGS! ^
        /module:!ParaMonteTest_MOD_DIR!     %=path to output ParaMonte Test module files=% ^
        /I:!ParaMonteTest_MOD_DIR!          %=path to required ParaMonte Test module files=%  ^
        /I:!ParaMonte_MOD_DIR!              %=path to input ParaMonte module files=%  ^
        /c !ParaMonteTest_SRC_DIR!\%%A      %=path to input ParaMonte Test source file=%  ^
        || (
            echo. 
            echo. -- ParaMonteTest - Fatal Error: compilation of the object file for %%A failed.
            echo. -- ParaMonteTest - build failed. exiting...
            echo. 
            set ERRORLEVEL=1
            exit /B
        )
    )
    echo.
)
@echo off
echo.


:LABEL_ParaMonteTest_EXE_ENABLED

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: generate test executable
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set ParaMonteTest_EXE_NAME=testParaMonte.exe

if !ParaMonteTest_EXE_ENABLED! NEQ true (
    echo.
    echo. -- ParaMonteTest - Warning: skipping ParaMonte library test exectuable build...
    echo.
    goto LABEL_ParaMonteTest_RUN_ENABLED
)

echo.

echo. -- ParaMonteTest - generating !LTYPE!-link ParaMonte library test executable at: !ParaMonteTest_BIN_DIR!

REM if !ParaMonte_LIB_ENABLED!==true (
    set ParaMonteTest_EXE_REQUIRED_OBJ_FILES=!ParaMonteTest_OBJ_DIR!\*.obj !ParaMonte_LIB_DIR!\!PMLIB_NAME!.lib
REM ) else (
REM     set ParaMonteTest_EXE_REQUIRED_OBJ_FILES=!ParaMonteTest_OBJ_DIR!\*.obj !ParaMonte_OBJ_DIR!\*.obj
REM )

:: delete the old executable first

echo. -- ParaMonteTest - deleting old executable (if any) at: !ParaMonteTest_BIN_DIR!\!ParaMonteTest_EXE_NAME!

cd !ParaMonteTest_BIN_DIR!
del !ParaMonteTest_EXE_NAME!
if !ERRORLEVEL!==1 (
    echo. 
    echo. -- ParaMonteTest - Fatal Error: deletion of the old executable at !ParaMonteTest_BIN_DIR!\!ParaMonteTest_EXE_NAME! failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)

:: build the executable

if !MPI_ENABLED!==true (

    call !FCL! !FCL_FLAGS! !FL_FLAGS! ^
    /module:!ParaMonteTest_MOD_DIR! ^
    /I:!ParaMonteTest_MOD_DIR! ^
    /I:!ParaMonte_MOD_DIR! ^
    !ParaMonteTest_EXE_REQUIRED_OBJ_FILES! ^
    /exe:!ParaMonteTest_BIN_DIR!\!ParaMonteTest_EXE_NAME!

) else (

    !FCL! !FCL_FLAGS! !FL_FLAGS! ^
    /module:!ParaMonteTest_MOD_DIR! ^
    /I:!ParaMonteTest_MOD_DIR! ^
    /I:!ParaMonte_MOD_DIR! ^
    !ParaMonteTest_EXE_REQUIRED_OBJ_FILES! ^
    /exe:!ParaMonteTest_BIN_DIR!\!ParaMonteTest_EXE_NAME!

)

if !ERRORLEVEL!==1 ( 
    echo. 
    echo. -- ParaMonteTest - Fatal Error: linking of the test object files may have likely failed.
    echo. -- ParaMonteTest - build may have likely failed. continuing...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
echo.


::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: run ParaMonte library test executable
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:LABEL_ParaMonteTest_RUN_ENABLED

:: run ParaMonte test
:: if !ParaMonteTest_RUN_ENABLED! NEQ true goto LABEL_EXAMPLE_ENABLED
if !ParaMonteTest_RUN_ENABLED! NEQ true (
    echo.
    echo. -- ParaMonteTest - Warning: skipping ParaMonte library test run...
    echo.
    goto :eof
)

echo. 
echo. ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
echo. ::::                                                                                                                            ::::
echo.                                                    Running ParaMonte Library Tests
echo. ::::                                                                                                                            ::::
echo. ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
echo.

:: copy necessary input files in the executable's directory

echo. -- ParaMonteTest - copying input files to the bin directory
echo. -- ParaMonteTest - from: !ParaMonteTest_SRC_DIR!\input   %= no need for final slash here =%
echo. -- ParaMonteTest -   to: !ParaMonteTest_BIN_DIR!\input\  %= final slash tells this is folder =%
xcopy /s /Y "!ParaMonteTest_SRC_DIR!\input" "!ParaMonteTest_BIN_DIR!\input\"
echo.

cd !ParaMonteTest_BIN_DIR!

if !MPI_ENABLED!==true (
    mpiexec -localonly -n !FOR_COARRAY_NUM_IMAGES! !ParaMonteTest_EXE_NAME!
) else (
    !ParaMonteTest_EXE_NAME!
)

if !ERRORLEVEL!==0 ( 
    echo.
    echo.
    echo. -- ParaMonteTest - ParaMonte library test successful. 
    echo.
) else ( 
    echo.
    echo.
    echo. -- ParaMonteTest - ParaMonte library test failed. exiting...
    echo.
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)

exit /B 0

