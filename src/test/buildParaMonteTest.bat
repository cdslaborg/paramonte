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
        /c !ParaMonteTest_SRC_DIR!\%%A      %=path to input ParaMonte Test source file=% ^
        || (
            REM if !ERRORLEVEL! NEQ 0 (
            echo.
            echo. -- ParaMonteTest - Fatal Error: compilation of the object file for %%A failed.
            echo. -- ParaMonteTest - build failed. exiting...
            echo.
            cd %~dp0
            set ERRORLEVEL=1
            exit /B 1
        )
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
        /I:!ParaMonteTest_MOD_DIR!          %=path to required ParaMonte Test module files=% ^
        /I:!ParaMonte_MOD_DIR!              %=path to input ParaMonte module files=% ^
        /c !ParaMonteTest_SRC_DIR!\%%A      %=path to input ParaMonte Test source file=% ^
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
del !ParaMonteTest_EXE_NAME! || (
    echo. 
    echo. -- ParaMonteTest - Fatal Error: deletion of the old executable at !ParaMonteTest_BIN_DIR!\!ParaMonteTest_EXE_NAME! failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
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
    /exe:!ParaMonteTest_BIN_DIR!\!ParaMonteTest_EXE_NAME! ^
    || (
        echo. 
        echo. -- ParaMonteTest - Fatal Error: linking of the test object files may have likely failed.
        echo. -- ParaMonteTest - build may have likely failed. continuing...
        echo. 
        cd %~dp0
        set ERRORLEVEL=1
        exit /B 1
    )

) else (

    !FCL! !FCL_FLAGS! !FL_FLAGS! ^
    /module:!ParaMonteTest_MOD_DIR! ^
    /I:!ParaMonteTest_MOD_DIR! ^
    /I:!ParaMonte_MOD_DIR! ^
    !ParaMonteTest_EXE_REQUIRED_OBJ_FILES! ^
    /exe:!ParaMonteTest_BIN_DIR!\!ParaMonteTest_EXE_NAME! ^
    || (
        echo. 
        echo. -- ParaMonteTest - Fatal Error: linking of the test object files may have likely failed.
        echo. -- ParaMonteTest - build may have likely failed. continuing...
        echo. 
        cd %~dp0
        set ERRORLEVEL=1
        exit /B 1
    )

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
xcopy /s /Y "!ParaMonteTest_SRC_DIR!\input" "!ParaMonteTest_BIN_DIR!\input\" || goto LABEL_copyErrorOccured
echo.

cd !ParaMonteTest_BIN_DIR!

if !MPI_ENABLED!==true (
    mpiexec -localonly -n !FOR_COARRAY_NUM_IMAGES! !ParaMonteTest_EXE_NAME! || (
        echo.
        echo.
        echo. -- ParaMonteTest - the ParaMonte library MPI-parallel test run failed. exiting...
        echo.
        cd %~dp0
        set ERRORLEVEL=1
        exit /B 1
    )
) else (
    !ParaMonteTest_EXE_NAME! || (
        echo.
        echo.
        echo. -- ParaMonteTest - the ParaMonte library test run failed. exiting...
        echo.
        cd %~dp0
        set ERRORLEVEL=1
        exit /B 1
    )
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

:LABEL_copyErrorOccured

echo. 
echo. -- ParaMonteExample!LANG_NAME! - Fatal Error: failed to copy contents. exiting...
echo. 
cd %~dp0
set ERRORLEVEL=1
exit /B 1
