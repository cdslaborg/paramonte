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
::
::
::  USAGE:
::
::      install.bat > install.bat.out 2>&1
::
::  This batch file configures the flags required for building ParaMonte library, tests, and examples on Windows Operating Systems.
::
::  Prerequisites:
::
::      See this page for illustrative instructions:
::
::          https://www.cdslab.org/recipes/programming/intel-parallel-studio-installation-windows/intel-parallel-studio-installation-windows
::
::      In sum, you need the following software/compilers:
::
::          --  A recent Microsoft Visual Studio (>2017). The community edition of Visual Studio can be downloaded and installed free of charge:
::              --  https://visualstudio.microsoft.com/vs/community/
::              --  Ensure C++ development tools are chosen to be installed at the at the time of installation as
::                  it is required for the intergation of Visual Studio with Intel Studio.
::
::          --  Install Intel Parallel Studio >2018 (after installing Microsoft Visual Studio >2017).
::              --  Intel Parallel Studio is free of charge for students, educators, and open-srouce developers.
::              --  Once Intel Studio is installed, open Intel's special Windows-command-prompt which
::                  automatically defines all of the prerequisite environmental variables.
::              --  Run this script on the command prompt: install.bat -language -build -memory

@echo off
set ERRORLEVEL=0
cd %~dp0

setlocal EnableDelayedExpansion

set "INSTALL_SCRIPT_NAME=ParaMonte install.bat"

:: parse arguments

REM type .\bmake\install_usage.txt

set FPP_FLAGS_EXTRA=
set LANG_LIST=
set BTYPE_LIST=
set LTYPE_LIST=
set MEMORY_LIST=
set PARALLELISM_LIST=
set FOR_COARRAY_NUM_IMAGES=
set ParaMonte_INSTALL_CLEANUP_ENABLED=true
set DRY_RUN=false
set MatDRAM_ENABLED=false

echo.
type .\auxil\.ParaMonteBanner
echo.

echo.
echo.-- !INSTALL_SCRIPT_NAME! - parsing input arguments...
echo.

:LABEL_parseArgLoop

set FLAG_SUPPORTED=true
set VALUE_SUPPORTED=true

if not "%1"=="" (

    echo.-- !INSTALL_SCRIPT_NAME! - processing: %1

    set FLAG=%1
    set VALUE=%2
    call :getLowerCase FLAG
    call :getLowerCase VALUE

    set FLAG_SUPPORTED=false
    set VALUE_SUPPORTED=false

    REM -lang

    if "!FLAG!"=="--lang" (
        set FLAG_SUPPORTED=true
        for %%a in ("!VALUE:/=" "!") do (
            set DELIM=
            if defined LANG_LIST set DELIM=/
            set LANG_LIST=!LANG_LIST!!DELIM!%%~a
            set VALUE_SUPPORTED=false
            for %%V in ( "c" "c++" "fortran" "matlab" "python" ) do ( if /I "%%~a"=="%%~V" set "VALUE_SUPPORTED=true" )
            if not !VALUE_SUPPORTED!==true goto LABEL_REPORT_ERR
        )
        shift
    )

    REM --build

    if "!FLAG!"=="--build" (
        set FLAG_SUPPORTED=true
        for %%a in ("!VALUE:/=" "!") do (
            set DELIM=
            if defined BTYPE_LIST set DELIM=/
            set BTYPE_LIST=!BTYPE_LIST!!DELIM!%%~a
            set VALUE_SUPPORTED=false
            for %%V in ( "release" "testing" "debug" ) do ( if /I "%%~a"=="%%~V" set "VALUE_SUPPORTED=true" )
            if not !VALUE_SUPPORTED!==true goto LABEL_REPORT_ERR
        )
        shift
    )

    REM --lib

    if "!FLAG!"=="--lib" (
        set FLAG_SUPPORTED=true
        for %%a in ("!VALUE:/=" "!") do (
            set DELIM=
            if defined LTYPE_LIST set DELIM=/
            set LTYPE_LIST=!LTYPE_LIST!!DELIM!%%~a
            set VALUE_SUPPORTED=false
            for %%V in ( "dynamic" "static" ) do ( if /I "%%~a"=="%%~V" set "VALUE_SUPPORTED=true" )
            if not !VALUE_SUPPORTED!==true goto LABEL_REPORT_ERR
        )
        shift
    )

    REM --par

    if "!FLAG!"=="--par" (
        set FLAG_SUPPORTED=true
        for %%a in ("!VALUE:/=" "!") do (
            set DELIM=
            if defined PARALLELISM_LIST set DELIM=/
            set PARALLELISM_LIST=!PARALLELISM_LIST!!DELIM!%%~a
            set VALUE_SUPPORTED=false
            for %%V in ( "none" "mpi" "cafsingle" "cafshared" ) do ( if /I "%%~a"=="%%~V" set "VALUE_SUPPORTED=true" )
            if not !VALUE_SUPPORTED!==true goto LABEL_REPORT_ERR
        )
        shift
    )

    REM --mem

    if "!FLAG!"=="--mem" (
        set FLAG_SUPPORTED=true
        for %%a in ("!VALUE:/=" "!") do (
            set DELIM=
            if defined MEMORY_LIST set DELIM=/
            set MEMORY_LIST=!MEMORY_LIST!!DELIM!%%~a
            set VALUE_SUPPORTED=false
            for %%V in ( "stack" "heap" ) do ( if /I "%%~a"=="%%~V" set "VALUE_SUPPORTED=true" )
            if not !VALUE_SUPPORTED!==true goto LABEL_REPORT_ERR
        )
        shift
    )

    REM --test_enabled

    if "!FLAG!"=="--test_enabled" (
        set FLAG_SUPPORTED=true
        set "TEST_ENABLED=!VALUE!"
        set VALUE_SUPPORTED=false
        if !TEST_ENABLED!==true set "VALUE_SUPPORTED=true"
        if !TEST_ENABLED!==false set "VALUE_SUPPORTED=true"
        if not !VALUE_SUPPORTED!==true goto LABEL_REPORT_ERR
        shift
    )

    REM --exam_enabled

    if "!FLAG!"=="--exam_enabled" (
        set FLAG_SUPPORTED=true
        set "EXAM_ENABLED=!VALUE!"
        set VALUE_SUPPORTED=false
        if !EXAM_ENABLED!==true set "VALUE_SUPPORTED=true"
        if !EXAM_ENABLED!==false set "VALUE_SUPPORTED=true"
        if not !VALUE_SUPPORTED!==true goto LABEL_REPORT_ERR
        shift
    )

    REM --matdram

    if "!FLAG!"=="--matdram" (
        set FLAG_SUPPORTED=true
        set MatDRAM_ENABLED=true
    )

    REM --nproc

    if "!FLAG!"=="--nproc" (
        set FLAG_SUPPORTED=true
        call :getUpperCase VALUE
        set "FOR_COARRAY_NUM_IMAGES=!VALUE!"
        shift
    )

    REM --clean

    if "!FLAG!"=="--clean" (
        set FLAG_SUPPORTED=true
        set "ParaMonte_INSTALL_CLEANUP_ENABLED=!VALUE!"
        set VALUE_SUPPORTED=false
        if !ParaMonte_INSTALL_CLEANUP_ENABLED!==true set "VALUE_SUPPORTED=true"
        if !ParaMonte_INSTALL_CLEANUP_ENABLED!==false set "VALUE_SUPPORTED=true"
        if not !VALUE_SUPPORTED!==true goto LABEL_REPORT_ERR
        shift
    )

    REM --clf preprocessor/Compiler/linker flags

    if "!FLAG!"=="--fpp" (
        set FLAG_SUPPORTED=true
        set FPP_FLAGS_EXTRA=%2
        shift
    )

    REM --dryrun

    if "!FLAG!"=="--dryrun" (
        set FLAG_SUPPORTED=true
        set DRY_RUN=true
    )

    REM --help

    if "!FLAG!"=="--help" (
        set FLAG_SUPPORTED=true
        type .\install.bat.usage.txt
        exit /b 0
    )

    if !FLAG_SUPPORTED! NEQ true goto LABEL_REPORT_ERR

    shift
    goto :LABEL_parseArgLoop

)

:LABEL_REPORT_ERR

REM check flag/value support

if "!FLAG_SUPPORTED!"=="true" (
    if "!VALUE_SUPPORTED!" NEQ "true" (
        echo.
        echo.-- !INSTALL_SCRIPT_NAME! - FATAL: The requested input value "!VALUE!" specified
        echo.-- !INSTALL_SCRIPT_NAME! - FATAL: with the input flag "!FLAG!" is not supported.
        goto LABEL_ERR
    )
) else (
    echo.
    echo.-- !INSTALL_SCRIPT_NAME! - FATAL: The requested input flag "!FLAG!" is not supported.
    goto LABEL_ERR
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: build MatDRAM if explicitly requested. WARNING: If true, all other builds will be disabled.
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if !MatDRAM_ENABLED!==true (
    call buildMatDRAM.bat || (
        echo.
        echo.-- !INSTALL_SCRIPT_NAME! - Fatal Error: The MatDRAM library build failed for the following configuration.
        echo.-- !INSTALL_SCRIPT_NAME! - If you cannot identify the cause of the failure, please report this error at:
        echo.-- !INSTALL_SCRIPT_NAME! -
        echo.-- !INSTALL_SCRIPT_NAME! -     https://github.com/cdslaborg/paramonte/issues
        echo.-- !INSTALL_SCRIPT_NAME! -
        echo.-- !INSTALL_SCRIPT_NAME! - gracefully exiting...
        echo.
        cd %~dp0
        set ERRORLEVEL=1
        exit /B 1
    )
    goto LABEL_EOF
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: echo warnings
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if defined PARALLELISM_LIST (
    for %%P in ("!PARALLELISM_LIST:/=" "!") do (
        set PARALLELISM=%%~P
        set INITIALS=!PARALLELISM:~0,3!
        if !INITIALS!==caf (
            if defined LANG_LIST (
                for %%G in ("!LANG_LIST:/=" "!") do (
                    if %%G NEQ "fortran" (
                        echo.
                        echo.-- !INSTALL_SCRIPT_NAME! - WARNING: The Coarray parallelism flag "--par %%~P" cannot be
                        echo.-- !INSTALL_SCRIPT_NAME! - WARNING: specified along with the %%~G language "--lang %%~G".
                        echo.-- !INSTALL_SCRIPT_NAME! - WARNING: This configuration will be ignored at build time.
                        REM goto LABEL_ERR
                    )
                )
            )
            if defined LTYPE_LIST (
                for %%L in ("!LTYPE_LIST:/=" "!") do (
                    echo.
                    echo.-- !INSTALL_SCRIPT_NAME! - WARNING: The Coarray parallelism flag "--par %%~P" cannot be
                    echo.-- !INSTALL_SCRIPT_NAME! - WARNING: specified along with the dynamic library build flag "--lib %%~L".
                    echo.-- !INSTALL_SCRIPT_NAME! - WARNING: This configuration will be ignored at build time.
                    REM goto LABEL_ERR
                )
            )
        )
    )
)

if defined LANG_LIST (
    if defined LTYPE_LIST (
        for %%G in ("!LANG_LIST:/=" "!") do (
            set LANG_IS_DYNAMIC=false
            if %%~G==matlab set LANG_IS_DYNAMIC=true
            if %%~G==python set LANG_IS_DYNAMIC=true
            if !LANG_IS_DYNAMIC!==true (
                for %%L in ("!LTYPE_LIST:/=" "!") do (
                    if %%~L==static (
                        echo.
                        echo.-- !INSTALL_SCRIPT_NAME! - WARNING: The dynamic library option "--lib %%~L" cannot be
                        echo.-- !INSTALL_SCRIPT_NAME! - WARNING: specified along with the %%~G language "--lang %%~G".
                        echo.-- !INSTALL_SCRIPT_NAME! - WARNING: This configuration will be ignored at build time.
                        REM goto LABEL_ERR
                    )
                )
            )
        )
    )
)

if defined LTYPE_LIST (
    if defined MEMORY_LIST (
        for %%L in ("!LTYPE_LIST:/=" "!") do (
            for %%M in ("!MEMORY_LIST:/=" "!") do (
                if %%~M==stack (
                    if %%~L==dynamic (
                        echo.
                        echo.-- !INSTALL_SCRIPT_NAME! - WARNING: The stack memory allocation option "--mem %%~M" cannot be
                        echo.-- !INSTALL_SCRIPT_NAME! - WARNING: specified along with the dynamic library build "--lib %%~L".
                        echo.-- !INSTALL_SCRIPT_NAME! - WARNING: This configuration will be ignored at build time.
                        REM goto LABEL_ERR
                    )
                )
            )
        )
    )
)

echo.

:: set build type

if not defined LANG_LIST        set LANG_LIST=c/c++/fortran/matlab/python
if not defined BTYPE_LIST       set BTYPE_LIST=release/testing/debug
if not defined LTYPE_LIST       set LTYPE_LIST=static/dynamic
if not defined MEMORY_LIST      set MEMORY_LIST=stack/heap
if not defined PARALLELISM_LIST set PARALLELISM_LIST=none/mpi/cafsingle/cafshared

REM remove redundancies

set TEMP=
set C_IS_MISSING=true
set CPP_IS_MISSING=true
set Fortran_IS_MISSING=true
set MATLAB_IS_MISSING=true
set Python_IS_MISSING=true
for %%G in ("!LANG_LIST:/=" "!") do (
    if %%~G==fortran (
        if !Fortran_IS_MISSING!==true (
            if not defined TEMP (
                set TEMP=%%~G
            ) else (
                set TEMP=!TEMP!/%%~G
            )
            set Fortran_IS_MISSING=false
        )
    )
    if %%~G==c (
        if !C_IS_MISSING!==true (
            if not defined TEMP (
                set TEMP=%%~G
            ) else (
                set TEMP=!TEMP!/%%~G
            )
            set C_IS_MISSING=false
        )
    )
    if %%~G==c++ (
        if !CPP_IS_MISSING!==true (
            if not defined TEMP (
                set TEMP=%%~G
            ) else (
                set TEMP=!TEMP!/%%~G
            )
            set CPP_IS_MISSING=false
        )
    )
    if %%~G==matlab (
        if !MATLAB_IS_MISSING!==true (
            if not defined TEMP (
                set TEMP=%%~G
            ) else (
                set TEMP=!TEMP!/%%~G
            )
            set MATLAB_IS_MISSING=false
        )
    )
    if %%~G==python (
        if !Python_IS_MISSING!==true (
            if not defined TEMP (
                set TEMP=%%~G
            ) else (
                set TEMP=!TEMP!/%%~G
            )
            set Python_IS_MISSING=false
        )
    )
)
set LANG_LIST=!TEMP!

REM build

if !DRY_RUN!==true (
    set FRESH_RUN=false
) else (
    set FRESH_RUN=true
)

echo. LANG_LIST=!LANG_LIST!
echo. BTYPE_LIST=!BTYPE_LIST!
echo. LTYPE_LIST=!LTYPE_LIST!
echo. MEMORY_LIST=!MEMORY_LIST!
echo. PARALLELISM_LIST=!PARALLELISM_LIST!

for %%G in ("!LANG_LIST:/=" "!") do (

    for %%B in ("!BTYPE_LIST:/=" "!") do (

        for %%L in ("!LTYPE_LIST:/=" "!") do (

            for %%M in ("!MEMORY_LIST:/=" "!") do (

                for %%P in ("!PARALLELISM_LIST:/=" "!") do (

                    set BENABLED=true

                    set ParaMonte_OBJ_ENABLED=!FRESH_RUN!
                    set ParaMonte_LIB_ENABLED=!FRESH_RUN!
                    set ParaMonteExample_EXE_ENABLED=!EXAM_ENABLED!
                    set ParaMonteExample_RUN_ENABLED=!EXAM_ENABLED!

                    set INTERFACE_LANGUAGE=%%~G
                    set BTYPE=%%~B
                    set LTYPE=%%~L
                    set HEAP_ARRAY_ENABLED=false
                    if %%~M==heap set HEAP_ARRAY_ENABLED=true

                    if %%~G==fortran (
                        set CFI_ENABLED=false
                        set ParaMonteTest_OBJ_ENABLED=!FRESH_RUN!
                        set ParaMonteTest_EXE_ENABLED=!FRESH_RUN!
                        set ParaMonteTest_RUN_ENABLED=!TEST_ENABLED!
                    ) else (
                        set CFI_ENABLED=true
                        set ParaMonteTest_OBJ_ENABLED=false
                        set ParaMonteTest_EXE_ENABLED=false
                        set ParaMonteTest_RUN_ENABLED=false
                    )

                    set CAFTYPE=
                    set CAF_ENABLED=false
                    set MPI_ENABLED=false
                    set PARALLELISM=%%~P
                    if !PARALLELISM!==mpi set MPI_ENABLED=true
                    set INITIALS=!PARALLELISM:~0,3!
                    if !INITIALS!==caf (
                        set CAF_ENABLED=true
                        set CAFTYPE=!PARALLELISM:~3!
                    )

                    if !LTYPE!==dynamic (
                        if !HEAP_ARRAY_ENABLED!==false set BENABLED=false
                        if !CAF_ENABLED!==true set BENABLED=false
                    )
                    if !CAF_ENABLED!==true (
                        if !MPI_ENABLED!==true (
                            set BENABLED=false
                        )
                        if !HEAP_ARRAY_ENABLED!==false (
                            set BENABLED=false
                        )
                        if !CFI_ENABLED!==true (
                            set BENABLED=false
                        )
                    )

                    if %%~G==matlab (
                        if !LTYPE!==static set BENABLED=false
                        if !LTYPE! NEQ dynamic set BENABLED=false
                        if !CAF_ENABLED!==true set BENABLED=false
                        if !HEAP_ARRAY_ENABLED!==false set BENABLED=false
                    )

                    if %%~G==python (
                        if !LTYPE!==static set BENABLED=false
                        if !LTYPE! NEQ dynamic set BENABLED=false
                        if !CAF_ENABLED!==true set BENABLED=false
                        if !HEAP_ARRAY_ENABLED!==false set BENABLED=false
                    )

                    if !BENABLED!==true (

                        echo.
                        echo.************************************************************************************************************************************
                        echo.**** ParaMonte - current library build: --lang %%~G --build %%~B --lib %%~L --mem %%~M --par %%~P
                        echo.************************************************************************************************************************************
                        echo.

                        call buildParaMonte.bat || (
                            echo.
                            echo.-- !INSTALL_SCRIPT_NAME! - Fatal Error: The ParaMonte library build failed for the following configuration:
                            echo.-- !INSTALL_SCRIPT_NAME! -
                            echo.-- !INSTALL_SCRIPT_NAME! -               language: %%~G
                            echo.-- !INSTALL_SCRIPT_NAME! -             build type: %%~B
                            echo.-- !INSTALL_SCRIPT_NAME! -           library type: %%~L
                            echo.-- !INSTALL_SCRIPT_NAME! -      memory allocation: %%~M
                            echo.-- !INSTALL_SCRIPT_NAME! -            parallelism: %%~P
                            echo.-- !INSTALL_SCRIPT_NAME! -
                            echo.-- !INSTALL_SCRIPT_NAME! - If you cannot identify the cause of the failure, please report this error at:
                            echo.-- !INSTALL_SCRIPT_NAME! -
                            echo.-- !INSTALL_SCRIPT_NAME! -     https://github.com/cdslaborg/paramonte/issues
                            echo.-- !INSTALL_SCRIPT_NAME! -
                            echo.-- !INSTALL_SCRIPT_NAME! - gracefully exiting...
                            echo.
                            cd %~dp0
                            set ERRORLEVEL=1
                            exit /B 1
                        )

                    ) else (

                        echo.
                        echo.-- !INSTALL_SCRIPT_NAME! - inconsistent configuration flags detected. skipping...
                        echo.

                    )

                    if !ERRORLEVEL! NEQ 0 (
                        echo.
                        echo.-- !INSTALL_SCRIPT_NAME! - Fatal Error: The ParaMonte library build failed for the following configuration:
                        echo.-- !INSTALL_SCRIPT_NAME! -
                        echo.-- !INSTALL_SCRIPT_NAME! -               language: %%~G
                        echo.-- !INSTALL_SCRIPT_NAME! -             build type: %%~B
                        echo.-- !INSTALL_SCRIPT_NAME! -           library type: %%~L
                        echo.-- !INSTALL_SCRIPT_NAME! -      memory allocation: %%~M
                        echo.-- !INSTALL_SCRIPT_NAME! -            parallelism: %%~P
                        echo.-- !INSTALL_SCRIPT_NAME! -
                        echo.-- !INSTALL_SCRIPT_NAME! - If you cannot identify the cause of the failure, please report this error at:
                        echo.-- !INSTALL_SCRIPT_NAME! -
                        echo.-- !INSTALL_SCRIPT_NAME! -     https://github.com/cdslaborg/paramonte/issues
                        echo.-- !INSTALL_SCRIPT_NAME! -
                        echo.-- !INSTALL_SCRIPT_NAME! - gracefully exiting...
                        echo.
                        cd %~dp0
                        set ERRORLEVEL=1
                        exit /B 1
                    )

                )

            )

        )

    )

    REM ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    REM :: if MATLAB, generate MatDRAM
    REM ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

REM    if %%~G==matlab (
REM
REM        echo.
REM        echo.-- !INSTALL_SCRIPT_NAME! - Generating MATLAB MatDRAM library...
REM        echo.
REM
REM        set MatDRAM_ORIGIN_PATH=.\bin\MATLAB
REM        set MatDRAM_DESTINATION_PATH=.\bin\MatDRAM
REM        echo.-- !INSTALL_SCRIPT_NAME! - copying the MatDRAM library files...
REM        echo.-- !INSTALL_SCRIPT_NAME! - from: !MatDRAM_ORIGIN_PATH!         %= no need for final slash here =%
REM        echo.-- !INSTALL_SCRIPT_NAME! -   to: !MatDRAM_DESTINATION_PATH!\   %= final slash tells this is folder =%
REM        xcopy /s /Y "!MatDRAM_ORIGIN_PATH!" "!MatDRAM_DESTINATION_PATH!\" || goto LABEL_copyErrorOccured
REM
REM        REM add the MatDRAM indicator file
REM
REM        REM xcopy /s /Y "!MatDRAM_DESTINATION_PATH!\paramonte\kernel\.MatDRAM" "!MatDRAM_DESTINATION_PATH!\auxil\" || goto LABEL_copyErrorOccured
REM
REM        REM delete the binary files
REM
REM        rd /s /q "!MatDRAM_DESTINATION_PATH!\paramonte\lib" >nul 2>&1 || goto LABEL_delErrorOccured
REM
REM        REM delete the mpi example file
REM
REM        del /s /q "!MatDRAM_DESTINATION_PATH!\main_mpi.m" >nul 2>&1 || goto LABEL_delErrorOccured
REM
REM    )

)

goto LABEL_EOF

:: subroutines

:getLowerCase
:: Subroutine to convert a variable VALUE to all lower case.
:: The argument for this subroutine is the variable NAME.
FOR %%i IN ("/=/" "+=+" "A=a" "B=b" "C=c" "D=d" "E=e" "F=f" "G=g" "H=h" "I=i" "J=j" "K=k" "L=l" "M=m" "N=n" "O=o" "P=p" "Q=q" "R=r" "S=s" "T=t" "U=u" "V=v" "W=w" "X=x" "Y=y" "Z=z") DO CALL SET "%1=%%%1:%%~i%%"
GOTO:EOF

:getUpperCase
:: Subroutine to convert a variable VALUE to all UPPER CASE.
:: The argument for this subroutine is the variable NAME.
FOR %%i IN ("/=/" "+=+" "a=A" "b=B" "c=C" "d=D" "e=E" "f=F" "g=G" "h=H" "i=I" "j=J" "k=K" "l=L" "m=M" "n=N" "o=O" "p=P" "q=Q" "r=R" "s=S" "t=T" "u=U" "v=V" "w=W" "x=X" "y=Y" "z=Z") DO CALL SET "%1=%%%1:%%~i%%"
GOTO:EOF

:getTitleCase
:: Subroutine to convert a variable VALUE to Title Case.
:: The argument for this subroutine is the variable NAME.
FOR %%i IN ("/=/" "+=+" " a= A" " b= B" " c= C" " d= D" " e= E" " f= F" " g= G" " h= H" " i= I" " j= J" " k= K" " l= L" " m= M" " n= N" " o= O" " p= P" " q= Q" " r= R" " s= S" " t= T" " u= U" " v= V" " w= W" " x= X" " y= Y" " z= Z") DO CALL SET "%1=%%%1:%%~i%%"
GOTO:EOF

:LABEL_ERR

echo.
echo.-- !INSTALL_SCRIPT_NAME! - To see the list of possible flags and associated values, try:
echo.-- !INSTALL_SCRIPT_NAME! -
echo.-- !INSTALL_SCRIPT_NAME! -     install.bat --help
echo.-- !INSTALL_SCRIPT_NAME! -
echo.-- !INSTALL_SCRIPT_NAME! - gracefully exiting the !INSTALL_SCRIPT_NAME! script.
echo.

exit /B 1

:LABEL_copyErrorOccured

echo.
echo. -- !INSTALL_SCRIPT_NAME! - Fatal Error: failed to copy contents. exiting...
echo.
cd %~dp0
set ERRORLEVEL=1
exit /B 1

:LABEL_delErrorOccured

echo.
echo. -- !INSTALL_SCRIPT_NAME! - Fatal Error: failed to delete contents. exiting...
echo.
cd %~dp0
set ERRORLEVEL=1
exit /B 1

:LABEL_EOF

:: undefine all configuration environmental flags

if !ParaMonte_INSTALL_CLEANUP_ENABLED!==true (
    echo.
    echo.-- !INSTALL_SCRIPT_NAME! - cleaning up the environment...
    call unconfigParaMonte.bat || (
        echo.
        echo. -- !BUILD_SCRIPT_NAME! - Fatal Error: the ParaMonte library cleanup failed. exiting...
        echo.
        cd %~dp0
        set ERRORLEVEL=1
        exit /B 1
    )
)

echo.
echo.-- !INSTALL_SCRIPT_NAME! - mission accomplished.
echo.

exit /B 0
