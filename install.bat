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

:: parse arguments

REM type .\bmake\install_usage.txt

set LANG_LIST=
set BTYPE_LIST=
set LTYPE_LIST=
set CAFTYPE_LIST=
set MPI_ENABLED_LIST=
set HEAP_ARRAY_ENABLED_LIST=
set FOR_COARRAY_NUM_IMAGES=
set ParaMonte_FLAG_CLEANUP_ENABLED=false
set DRY_RUN=false

echo.
type .\auxil\ParaMonteBanner.txt
echo.

echo.
echo.-- ParaMonte - parsing input arguments...
echo.

:LABEL_parseArgLoop

set FLAG_SUPPORTED=true
set VALUE_SUPPORTED=true

if not "%1"=="" (

    echo.-- ParaMonte - processing: %1

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
            for %%V in ( "c" "fortran" "python" ) do ( if /I "%%~a"=="%%~V" set "VALUE_SUPPORTED=true" )
            if !VALUE_SUPPORTED! NEQ true goto LABEL_REPORT_ERR
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
            if !VALUE_SUPPORTED! NEQ true goto LABEL_REPORT_ERR
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
            if !VALUE_SUPPORTED! NEQ true goto LABEL_REPORT_ERR
        )
        shift
    )

    REM --caf

    if "!FLAG!"=="--caf" (
        set FLAG_SUPPORTED=true
        for %%a in ("!VALUE:/=" "!") do (
            set DELIM=
            if defined CAFTYPE_LIST set DELIM=/
            set CAFTYPE_LIST=!CAFTYPE_LIST!!DELIM!%%~a
            set VALUE_SUPPORTED=false
            for %%V in ( "none" "single" "shared" ) do ( if /I "%%~a"=="%%~V" set "VALUE_SUPPORTED=true" )
            if !VALUE_SUPPORTED! NEQ true goto LABEL_REPORT_ERR
        )
        shift
    )

    REM --mpi_enabled

    if "!FLAG!"=="--mpi_enabled" (
        set FLAG_SUPPORTED=true
        for %%a in ("!VALUE:/=" "!") do (
            set DELIM=
            if defined MPI_ENABLED_LIST set DELIM=/
            set MPI_ENABLED_LIST=!MPI_ENABLED_LIST!!DELIM!%%~a
            set VALUE_SUPPORTED=false
            for %%V in ( "true" "false" ) do ( if /I "%%~a"=="%%~V" set "VALUE_SUPPORTED=true" )
            if !VALUE_SUPPORTED! NEQ true goto LABEL_REPORT_ERR
        )
        shift
    )

    REM --heap_enabled

    if "!FLAG!"=="--heap_enabled" (
        set FLAG_SUPPORTED=true
        for %%a in ("!VALUE:/=" "!") do (
            set DELIM=
            if defined HEAP_ARRAY_ENABLED_LIST set DELIM=/
            set HEAP_ARRAY_ENABLED_LIST=!HEAP_ARRAY_ENABLED_LIST!!DELIM!%%~a
            set VALUE_SUPPORTED=false
            for %%V in ( "true" "false" ) do ( if /I "%%~a"=="%%~V" set "VALUE_SUPPORTED=true" )
            if !VALUE_SUPPORTED! NEQ true goto LABEL_REPORT_ERR
        )
        shift
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
        set "ParaMonte_FLAG_CLEANUP_ENABLED=true"
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

    shift
    goto :LABEL_parseArgLoop

)

:LABEL_REPORT_ERR

REM check flag/value support

if "!FLAG_SUPPORTED!"=="true" (
    if "!VALUE_SUPPORTED!" NEQ "true" (
        echo.
        echo.-- ParaMonte - FATAL: The requested input value "!VALUE!" specified 
        echo.-- ParaMonte - FATAL: with the input flag "!FLAG!" is not supported.
        goto LABEL_ERR
    )
) else (
    echo.
    echo.-- ParaMonte - FATAL: The requested input flag "!FLAG!" is not supported.
    goto LABEL_ERR
)

REM echo warnings

if defined CAFTYPE_LIST (
    for %%c in ("!CAFTYPE_LIST:/=" "!") do (
        if %%c NEQ "none" (
            if defined MPI_ENABLED_LIST (
                for %%a in ("!MPI_ENABLED_LIST:/=" "!") do (
                    if %%~a==true (
                        echo.
                        echo.-- ParaMonte - WARNING: The Coarray flag "--caf !CAFTYPE_LIST!" cannot 
                        echo.-- ParaMonte - WARNING: be specified along with the MPI flag "--mpi_enabled %%~a".
                        echo.-- ParaMonte - WARNING: This configuration will be ignored at build time.
                        REM goto LABEL_ERR
                    )
                )
            )
            if defined LANG_LIST (
                for %%l in ("!LANG_LIST:/=" "!") do (
                    if %%l NEQ "fortran" (
                        if defined CAFTYPE_LIST (
                            echo.
                            echo.-- ParaMonte - WARNING: The Coarray flag "--caf !CAFTYPE_LIST!" cannot be 
                            echo.-- ParaMonte - WARNING: specified along with the %%~l language "--lang %%~l".
                            echo.-- ParaMonte - WARNING: This configuration will be ignored at build time.
                            REM goto LABEL_ERR
                        )
                    )
                )
            )
            if defined LANG_LIST (
                for %%l in ("!LANG_LIST:/=" "!") do (
                    if %%l NEQ "fortran" (
                        if defined CAFTYPE_LIST (
                            echo.
                            echo.-- ParaMonte - WARNING: The Coarray flag "--caf !CAFTYPE_LIST!" cannot be 
                            echo.-- ParaMonte - WARNING: specified along with the %%~l language "--lang %%~l".
                            echo.-- ParaMonte - WARNING: This configuration will be ignored at build time.
                            REM goto LABEL_ERR
                        )
                    )
                )
            )
        )
    )
)

if defined LANG_LIST (
    if defined LTYPE_LIST (
        for %%l in ("!LANG_LIST:/=" "!") do (
            if %%~l==python (
                for %%a in ("!LTYPE_LIST:/=" "!") do (
                    if %%~a==static (
                        echo.
                        echo.-- ParaMonte - WARNING: The dynamic library option "--lib %%~a" cannot be 
                        echo.-- ParaMonte - WARNING: specified along with the %%~l language "--lang %%~l".
                        echo.-- ParaMonte - WARNING: This configuration will be ignored at build time.
                        REM goto LABEL_ERR
                    )
                )
            )
        )
    )
)

if defined LTYPE_LIST (
    if defined HEAP_ARRAY_ENABLED_LIST (
        for %%l in ("!LTYPE_LIST:/=" "!") do (
            for %%h in ("!HEAP_ARRAY_ENABLED_LIST:/=" "!") do (
                if %%~h==false (
                    if %%~l==dynamic (
                        echo.
                        echo.-- ParaMonte - WARNING: The stack memory allocation option "--heap_enabled %%~h" cannot be 
                        echo.-- ParaMonte - WARNING: specified along with the dynamic library build "--lib %%~l".
                        echo.-- ParaMonte - WARNING: This configuration will be ignored at build time.
                        REM goto LABEL_ERR
                    )
                )
            )
        )
    )
)

echo.

:: set build type

if not defined LANG_LIST                        set LANG_LIST=c/fortran/python
if not defined BTYPE_LIST                       set BTYPE_LIST=release/testing/debug
if not defined LTYPE_LIST                       set LTYPE_LIST=static/dynamic
if not defined CAFTYPE_LIST                     set CAFTYPE_LIST=none/single/shared
if not defined MPI_ENABLED_LIST                 set MPI_ENABLED_LIST=true/false
if not defined HEAP_ARRAY_ENABLED_LIST          set HEAP_ARRAY_ENABLED_LIST=true/false
if not defined ParaMonte_FLAG_CLEANUP_ENABLED   set ParaMonte_FLAG_CLEANUP_ENABLED=false

REM remove redundancies

set TEMP=
set C_IS_MISSING=true
set Fortran_IS_MISSING=true
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
    if %%~G==python (
        if !C_IS_MISSING!==true (
            if not defined TEMP (
                set TEMP=%%~G
            ) else (
                set TEMP=!TEMP!/%%~G
            )
            set C_IS_MISSING=false
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

REM set PAR_TYPE_LIST=!MPI_ENABLED_LIST!/!CAFTYPE_LIST!

echo. LANG_LIST=!LANG_LIST!
echo. BTYPE_LIST=!BTYPE_LIST!
echo. LTYPE_LIST=!LTYPE_LIST!
echo. HEAP_ARRAY_ENABLED_LIST=!HEAP_ARRAY_ENABLED_LIST!
echo. CAFTYPE_LIST=!CAFTYPE_LIST!
echo. MPI_ENABLED_LIST=!MPI_ENABLED_LIST!

for %%G in ("!LANG_LIST:/=" "!") do (

    for %%B in ("!BTYPE_LIST:/=" "!") do (

        for %%L in ("!LTYPE_LIST:/=" "!") do (

            for %%H in ("!HEAP_ARRAY_ENABLED_LIST:/=" "!") do (

                for %%C in ("!CAFTYPE_LIST:/=" "!") do (

                    for %%M in ("!MPI_ENABLED_LIST:/=" "!") do (

                        set BENABLED=true

                        set ParaMonte_OBJ_ENABLED=!FRESH_RUN!
                        set ParaMonte_LIB_ENABLED=!FRESH_RUN!
                        set ParaMonteExample_EXE_ENABLED=true
                        set ParaMonteExample_RUN_ENABLED=true

                        set BTYPE=%%~B
                        set LTYPE=%%~L
                        set CAFTYPE=%%~C
                        set MPI_ENABLED=%%~M
                        set HEAP_ARRAY_ENABLED=%%~H

                        if %%~G==fortran (
                            set CFI_ENABLED=false
                            set ParaMonteTest_OBJ_ENABLED=!FRESH_RUN!
                            set ParaMonteTest_EXE_ENABLED=!FRESH_RUN!
                            set ParaMonteTest_RUN_ENABLED=!FRESH_RUN!
                        ) else (
                            set CFI_ENABLED=true
                            set ParaMonteTest_OBJ_ENABLED=false
                            set ParaMonteTest_EXE_ENABLED=false
                            set ParaMonteTest_RUN_ENABLED=false
                        )

                        set CAF_ENABLED=true
                        if !CAFTYPE!==none set CAF_ENABLED=false

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

                        if !BENABLED!==true (

                            echo.
                            echo.************************************************************************************************************************************
                            echo.**** ParaMonte - current library build: --lang %%~G --build %%~B --lib %%~L --heap_enabled %%~H --mpi_enabled %%~M --caf %%~C
                            echo.************************************************************************************************************************************
                            echo.

                            call buildParaMonte.bat

                        )

                        if !ERRORLEVEL! NEQ 0 (
                            echo.
                            echo.-- ParaMonte - Fatal Error: ParaMonte library build failed for the following configuration:
                            echo.-- ParaMonte - 
                            echo.-- ParaMonte -               language: %%~G
                            echo.-- ParaMonte -             build type: %%~B
                            echo.-- ParaMonte -           library type: %%~L
                            echo.-- ParaMonte -     heap array enabled: %%~H
                            echo.-- ParaMonte -           Coarray type: %%~C
                            echo.-- ParaMonte -            MPI enabled: %%~M
                            echo.-- ParaMonte - 
                            echo.-- ParaMonte - Please report this error at: 
                            echo.-- ParaMonte - 
                            echo.-- ParaMonte -     https://github.com/cdslaborg/paramonte/issues
                            echo.-- ParaMonte - 
                            echo.-- ParaMonte - gracefully exiting...
                            echo.
                            cd %~dp0
                            set ERRORLEVEL=1
                            exit /B 1
                        )

                    )

                )

            )

        )

    )

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
echo.-- ParaMonte - To see the list of possible flags and associated values, try:
echo.-- ParaMonte - 
echo.-- ParaMonte -     install.bat --help
echo.-- ParaMonte - 
echo.-- ParaMonte - gracefully exiting ParaMonte build. 
echo.

exit /B 1

:LABEL_EOF

echo.
echo.-- ParaMonte - mission accomplished. 
echo.

exit /B 0
