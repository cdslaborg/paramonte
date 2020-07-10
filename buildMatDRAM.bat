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

:: NOTE: This is not a standalone build-script. It must only be called by install.bat script in the root directory of the project.

set ParaMonte_ROOT_DIR=%~dp0
set ParaMonteExample_SRC_DIR=!ParaMonte_ROOT_DIR!example
set ParaMonteInterface_SRC_DIR=!ParaMonte_ROOT_DIR!src\interface
echo. -- !BUILD_SCRIPT_NAME! - project root directory: !ParaMonte_ROOT_DIR!

setlocal EnableDelayedExpansion

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: build MatDRAM library example objects and executable
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

echo.
echo. ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
echo. ::::                                                                                                                            ::::
echo.                                                   MatDRAM MATLAB Library Build
echo. ::::                                                                                                                            ::::
echo. ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
echo.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: make bin directory
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set ParaMonte_BIN_DIR=!ParaMonte_ROOT_DIR!\bin
echo. -- MatDRAM - The MatDRAM binaries directory: !ParaMonte_BIN_DIR!
if not exist !ParaMonte_BIN_DIR! (
    mkdir "!ParaMonte_BIN_DIR!\"
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: setup examples' interface language
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set LANG_FILE_EXT=m
set LANG_NAME=MATLAB

if not defined LANG_NAME (
    echo. 
    echo. -- MatDRAMExample - Fatal Error: unrecognized or no language specified. exiting...
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

set MatDRAM_BLD_DIR=!ParaMonte_ROOT_DIR!build\MatDRAM\example

set ParaMonteInterface_SRC_DIR_CURRENT=!ParaMonteInterface_SRC_DIR!\!LANG_NAME!

echo. 
echo. -- MatDRAM - generating the MatDRAM library examples in !LANG_NAME! language...
echo. -- MatDRAM - The MatDRAM !LANG_NAME! examples directory: !MatDRAM_BLD_DIR!

for %%e in (!EXAM_LIST!) do ( 

    set EXAM_NAME=%%e

    set MatDRAM_BLD_DIR_CURRENT=!MatDRAM_BLD_DIR!\!EXAM_NAME!
    echo. -- MatDRAM - The MatDRAM library !EXAM_NAME! example directory: !MatDRAM_BLD_DIR_CURRENT!
    if exist !MatDRAM_BLD_DIR_CURRENT! (
        echo. -- MatDRAM - previous example build detected. deleting the old contents...
        rmdir /S /Q !MatDRAM_BLD_DIR_CURRENT! || goto LABEL_rmdirErrorOccured
        echo. -- MatDRAM - regenerating the MatDRAM library !EXAM_NAME! example directory: !MatDRAM_BLD_DIR_CURRENT!
    )
    mkdir "!MatDRAM_BLD_DIR_CURRENT!\"

    REM The MatDRAM library example required files

    echo. -- MatDRAM - copying the MatDRAM library !EXAM_NAME! example required files in !LANG_NAME! language...

    REM The MatDRAM library README.md file

    REM echo. -- MatDRAM - copying the MatDRAM library README.md file...
    REM echo. -- MatDRAM - from: !ParaMonteInterface_SRC_DIR_CURRENT!\README.md
    REM echo. -- MatDRAM -   to: !MatDRAM_BLD_DIR_CURRENT!\paramonte\README.md
    REM copy /y "!ParaMonteInterface_SRC_DIR_CURRENT!\README.md" "!MatDRAM_BLD_DIR_CURRENT!\README.md" || goto LABEL_copyErrorOccured

    REM The MatDRAM library CHANGES.md file

    REM echo. -- MatDRAM - copying the MatDRAM library CHANGES.md file...
    REM echo. -- MatDRAM - from: !ParaMonteInterface_SRC_DIR_CURRENT!\CHANGES.md
    REM echo. -- MatDRAM -   to: !MatDRAM_BLD_DIR_CURRENT!\paramonte\CHANGES.md
    REM copy "!ParaMonteInterface_SRC_DIR_CURRENT!\CHANGES.md" "!MatDRAM_BLD_DIR_CURRENT!\CHANGES.md" || goto LABEL_copyErrorOccured

    REM The MatDRAM library license file

    echo. -- MatDRAM - copying the MatDRAM library license file...
    echo. -- MatDRAM - from: !ParaMonte_ROOT_DIR!\LICENSE
    echo. -- MatDRAM -   to: !MatDRAM_BLD_DIR_CURRENT!\LICENSE
    copy "!ParaMonte_ROOT_DIR!\LICENSE" "!MatDRAM_BLD_DIR_CURRENT!\LICENSE" || goto LABEL_copyErrorOccured

    REM The MatDRAM library interface files

    echo. -- MatDRAM - from: !ParaMonteInterface_SRC_DIR_CURRENT!\MatDRAM %= no need for final slash here =%
    echo. -- MatDRAM -   to: !MatDRAM_BLD_DIR_CURRENT!\paramonte\  %= final slash tells this is folder =%
    xcopy /s /Y "!ParaMonteInterface_SRC_DIR_CURRENT!\paramonte" "!MatDRAM_BLD_DIR_CURRENT!\paramonte\" || goto LABEL_copyErrorOccured

    REM The MatDRAM library banner file

    echo. -- MatDRAM - copying the MatDRAM library auxiliary files
    echo. -- MatDRAM - from: !ParaMonteInterface_SRC_DIR!\auxil                     %= no need for final slash here =%
    echo. -- MatDRAM -   to: !MatDRAM_BLD_DIR_CURRENT!\paramonte\auxil\    %= final slash tells this is folder =%
    xcopy /s /Y "!ParaMonteInterface_SRC_DIR!\auxil" "!MatDRAM_BLD_DIR_CURRENT!\paramonte\auxil\" || goto LABEL_copyErrorOccured
    echo.

    REM The MatDRAM library kernel version file (must appear only after the above)

    echo. -- MatDRAM - copying the MatDRAM library kernel version file...
    echo. -- MatDRAM - from: !ParaMonte_ROOT_DIR!\.VERSION %= no need for final slash here =%
    echo. -- MatDRAM -   to: !MatDRAM_BLD_DIR_CURRENT!\paramonte\auxil\.VERSION_KERNEL
    copy /y "!ParaMonte_ROOT_DIR!\.VERSION" "!MatDRAM_BLD_DIR_CURRENT!\paramonte\auxil\.VERSION_KERNEL" || goto LABEL_copyErrorOccured

    REM The MatDRAM library interface version file (must appear only after the above)

    echo. -- MatDRAM - copying the MatDRAM library interface version file...
    echo. -- MatDRAM - from: !ParaMonteInterface_SRC_DIR_CURRENT!\.VERSION %= no need for final slash here =%
    echo. -- MatDRAM -   to: !MatDRAM_BLD_DIR_CURRENT!\paramonte\auxil\.VERSION_INTERFACE
    copy /y "!ParaMonteInterface_SRC_DIR_CURRENT!\.VERSION" "!MatDRAM_BLD_DIR_CURRENT!\paramonte\auxil\.VERSION_INTERFACE" || goto LABEL_copyErrorOccured

    REM The MatDRAM library example input files

    REM set ParaMonteExample_INP_DIR_CURRENT=!ParaMonteExample_SRC_DIR!\!EXAM_NAME!\input
    REM echo. -- MatDRAM - copying the MatDRAM library !EXAM_NAME! example input files in !LANG_NAME! language...
    REM echo. -- MatDRAM - from: !ParaMonteExample_INP_DIR_CURRENT!     %= no need for final slash here =%
    REM echo. -- MatDRAM -   to: !MatDRAM_BLD_DIR_CURRENT!\    %= final slash tells this is folder =%
    REM xcopy /s /Y "!ParaMonteExample_INP_DIR_CURRENT!" "!MatDRAM_BLD_DIR_CURRENT!\" || goto LABEL_copyErrorOccured

    REM The MatDRAM library example source files

    set ParaMonteExample_SRC_DIR_CURRENT=!ParaMonteExample_SRC_DIR!\!EXAM_NAME!\!LANG_NAME!
    echo. -- MatDRAM - copying the MatDRAM library !EXAM_NAME! example source files in !LANG_NAME! language...

    set mainFileName=main.!LANG_FILE_EXT!
    copy "!ParaMonteExample_SRC_DIR!\!mainFileName!" "!MatDRAM_BLD_DIR_CURRENT!\" || goto LABEL_copyErrorOccured

    echo. -- MatDRAM - from: !ParaMonteExample_SRC_DIR_CURRENT! %= no need for final slash here =%
    echo. -- MatDRAM -   to: !MatDRAM_BLD_DIR_CURRENT! %= final slash tells this is folder =%
    xcopy /s /Y "!ParaMonteExample_SRC_DIR_CURRENT!" "!MatDRAM_BLD_DIR_CURRENT!\" || goto LABEL_copyErrorOccured

)

echo. 

if %ERRORLEVEL%==1 (
    echo. 
    echo. -- MatDRAM - Fatal Error: build failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
if %ERRORLEVEL%==0 (
    echo.
    echo.
    echo. -- MatDRAM - The MatDRAM library example build successful. 
    echo.
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: copy the first example to the bin directory
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set MatDRAM_BLD_DIR_CURRENT=!MatDRAM_BLD_DIR!\mvn

set MatDRAMExample_BIN_DIR_CURRENT=!ParaMonte_BIN_DIR!\MatDRAM

echo. -- MatDRAM - The MatDRAM !LANG_NAME! library binary directory: !MatDRAMExample_BIN_DIR_CURRENT!

if not exist !MatDRAMExample_BIN_DIR_CURRENT! (
    mkdir "!MatDRAMExample_BIN_DIR_CURRENT!\"
)

echo. -- MatDRAM - copying the MatDRAM library files to the bin folder...
echo. -- MatDRAM - from: !MatDRAM_BLD_DIR_CURRENT! %= no need for final slash here =%
echo. -- MatDRAM -   to: !MatDRAMExample_BIN_DIR_CURRENT! %= final slash tells this is folder =%
xcopy /s /Y /e /v /i "!MatDRAM_BLD_DIR_CURRENT!" "!MatDRAMExample_BIN_DIR_CURRENT!" || goto LABEL_copyErrorOccured

if %ERRORLEVEL%==1 (
    echo. 
    echo. -- MatDRAM - Fatal Error: build failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
if %ERRORLEVEL%==0 (
    echo.
    echo.
    echo. -- MatDRAM - The MatDRAM library example build successful. 
    echo.
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: copy the first example to the bin directory
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if !ParaMonteExample_RUN_ENABLED!==true (

    for %%e in (!EXAM_LIST!) do ( 

        set EXAM_NAME=%%e
        echo. -- MatDRAM - Running the MatDRAM library's !EXAM_NAME! example.

        REM The MatDRAM library example build and run if requested

        set MatDRAM_BLD_DIR_CURRENT=!MatDRAM_BLD_DIR!\!EXAM_NAME!
        cd !MatDRAM_BLD_DIR_CURRENT!

        matlab -batch "main" && (
            echo.
            echo.
            echo. -- MatDRAM - The MatDRAM library example build/run appears to have succeeded.
            echo.
        ) || (
            echo. 
            echo. -- MatDRAM - Fatal Error: The MatDRAM library example build/run failed. exiting...
            echo. 
            cd %~dp0
            set ERRORLEVEL=1
            exit /B 1
        )

        cd %~dp0

    )

)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: quit
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

cd %~dp0

exit /B 0

:LABEL_copyErrorOccured

echo. 
echo. -- MatDRAM - Fatal Error: failed to copy contents. exiting...
echo. 
cd %~dp0
set ERRORLEVEL=1
exit /B 1

:LABEL_rmdirErrorOccured

echo. 
echo. -- MatDRAM - Fatal Error: failed to delete old contents. exiting...
echo. 
cd %~dp0
set ERRORLEVEL=1
exit /B 1
