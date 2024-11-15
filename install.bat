::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::                                                                                                                            ::::
::::    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           ::::
::::                                                                                                                            ::::
::::    Copyright (C) 2012-present, The Computational Data Science Lab                                                          ::::
::::                                                                                                                            ::::
::::    This file is part of the ParaMonte library.                                                                             ::::
::::                                                                                                                            ::::
::::    LICENSE                                                                                                                 ::::
::::                                                                                                                            ::::
::::       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          ::::
::::                                                                                                                            ::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

::
:: See the file install.bat.usage in the same folder for usage guidelines of this Batch script.
::

@echo off
set ERRORLEVEL=0
setlocal EnableDelayedExpansion
set BUILD_SCRIPT_NAME=install.bat
set "script_name=install.bat"
:: change directory to the folder containing this batch file
cd "%~dp0"

REM WARNING: paramonte_dir ends with a forward slash.

set "paramonte_dir=%~dp0"
set "paramonte_src_dir=!paramonte_dir!src"
set "paramonte_auxil_dir=!paramonte_dir!auxil"
set "paramonte_example_dir=!paramonte_dir!example"
set "paramonte_external_dir=!paramonte_dir!external"
set "paramonte_benchmark_dir=!paramonte_dir!benchmark"
set "paramonte_src_fortran_dir=!paramonte_src_dir!\fortran"
set "paramonte_src_fortran_main_dir=!paramonte_src_fortran_dir!\main"
set "paramonte_src_fortran_test_dir=!paramonte_src_fortran_dir!\test"
set "paramonte_req_dir=!paramonte_external_dir!\prerequisites"

set "paramonte_web_github=https://github.com/cdslaborg/paramonte"
set "paramonte_web_github_issues=!paramonte_web_github!/issues/new/choose"
set "build_name=ParaMontes"

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: Set up platform.
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set os=!OS!
set arch=!PLATFORM!
if !arch!==x64 (
    set arch=amd64
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: Set up color coding.
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

REM \warning
REM ESC contains the escape ASCII character.
set "ESC="
set "ColorReset=!ESC![0m"
set "ColorBold=!ESC![1m"
set "Red=!ESC![31m"
set "Green=!ESC![32m"
set "Yellow=!ESC![33m"
set "Blue=!ESC![34m"
set "Magenta=!ESC![35m"
set "Cyan=!ESC![36m"
set "White=!ESC![37m"
set "BoldRed=!ESC![1;31m"
set "BoldGreen=!ESC![1;32m"
set "BoldYellow=!ESC![1;33m"
set "BoldBlue=!ESC![1;34m"
set "BoldMagenta=!ESC![1;35m"
set "BoldCyan=!ESC![1;36m"
set "BoldWhite=!ESC![1;37m"

set "pmcolor=!BoldCyan!"
set "pmattn= !pmcolor!-- !build_name! !script_name! -!ColorReset!"
set "pmnote=!pmattn! !BoldYellow!NOTE:!ColorReset!"
set "pmwarn=!pmattn! !BoldMagenta!WARNING:!ColorReset!"
set "pmfatal=!pmattn! !BoldRed!FATAL ERROR:!ColorReset!"
set "warning=!BoldMagenta!WARNING!ColorReset!"

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: Fetch and set the ParaMonte library Fortran version.
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set "ParaMonteVersion="

cd "!paramonte_auxil_dir!"
for /f "tokens=*" %%i in ('head.bat 1 "..\VERSION.md"') do set "ParaMonteVersion=%%i"
cd %~dp0

REM uncomment the following conditional block to set the ParaMonte version in the source files via the preprocessor macros.
REM This is, however, not recommended. Generating the include source file is the preferred default method of
REM the ParaMonte version to the binaries. Starting ParaMonte release 1.4.2, this is the default behavior.
REM set "FPP_PARAMONTE_VERSION_FLAG="
REM if defined ParaMonteVersion (
REM     set FPP_PARAMONTE_VERSION_FLAG=/define:PARAMONTE_VERSION='!ParaMonteVersion!'
REM )

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: fetch ParaMonte library Fortran release date
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set dayName=!date:~0,3!
set year=!date:~10,4!
set day=!date:~7,2!

set m=100
for %%m in (Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec) do (
    set /a m+=1
    set month[!m:~-2!]=%%m
)
set monthNow=%date:~3,3%
set monthNow=%monthNow: =%
set monthName=!month[%monthNow%]!
set ParaMonteRelease=!dayName!.!monthName!.!day!.!year!

set SERIAL_ENABLED=true

REM echo.
REM type "!paramonte_auxil_dir!\.paramonte.banner"
REM echo.

echo.
echo. ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
echo. ::::                                                                                  ::::
echo.                ParaMonte library version !ParaMonteVersion! build on Windows
echo.                                     !ParaMonteRelease!
echo. ::::                                                                                  ::::
echo. ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
echo.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: parse arguments
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

REM type "!paramonte_auxil_dir!\install_usage.txt"

echo.
echo.!pmnote! parsing the command-line arguments...
echo.

set bdir=
set FOR_COARRAY_NUM_IMAGES=3
set "ddir=!paramonte_dir!bin"
set "flag_ddir=-Dddir="!ddir!""
set "cmakeBuildGenerator="
set "makename="

set list_build=
set list_checking=
set list_fc=
set list_lang=
set list_lib=
set list_mem=
set list_par=

set flag_bench=
set flag_benchpp=
set flag_blas=
set flag_codecov=
set flag_cfi=
set flag_deps=
set flag_exam=
set flag_exampp=
set flag_fpp=
set flag_fresh=
set flag_G=
set flag_j=
set flag_lapack=
set flag_matlabroot=
set flag_me=
set flag_mod=
set flag_nproc=
set flag_perfprof=
set flag_pdt=
set flag_purity=
set flag_test=

set flag_ski=
set flag_iki=
set flag_lki=
set flag_cki=
set flag_rki=

REM
REM The variable `ntry` is to bypass the need for duplicate build with CMake for development and testing times.
REM The duplicate build with CMake is required to ensure the generation of FPP source files in the output package.
REM This need for duplicate builds is an issue within the current CMake build scripts of the ParaMonte library that
REM must be resolved in the future with a better solution.
REM

set "flag_dev=-Ddev_enabled=0"
set ntry=2

REM
REM MATLAB MEX variables (must be removed once CMake FindMatlab.cmake module bug for Windows is resolved.)
REM

set "matlabroot="

REM
REM Echo the ParaMonte banner.
REM

echo.
type "!paramonte_auxil_dir!\.paramonte.banner"
echo.

:LABEL_parseArgLoop

if not "%1"=="" (

    echo.!pmnote! processing: %1

    set "FLAG=%~1"
    set "VALUE=%~2"
    REM call :getLowerCase FLAG
    REM call :getLowerCase VALUE

    set FLAG_SUPPORTED=
    set VALUE_SUPPORTED=true

    REM :::::::::
    REM list args
    REM :::::::::

    REM --build

    if "!FLAG!"=="--build" (
        set FLAG_SUPPORTED=true
        set "list_build=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --checking

    if "!FLAG!"=="--checking" (
        set FLAG_SUPPORTED=true
        set "list_checking=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --fc

    if "!FLAG!"=="--fc" (
        set FLAG_SUPPORTED=true
        set "list_fc=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --lang

    if "!FLAG!"=="--lang" (
        set FLAG_SUPPORTED=true
        set "list_lang=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --lib

    if "!FLAG!"=="--lib" (
        set FLAG_SUPPORTED=true
        set "list_lib=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --mem

    if "!FLAG!"=="--mem" (
        set FLAG_SUPPORTED=true
        set "list_mem=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --par

    if "!FLAG!"=="--par" (
        set FLAG_SUPPORTED=true
        set "list_par=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM :::::::::
    REM flag args
    REM :::::::::

    REM --bench

    if "!FLAG!"=="--bench" (
        set FLAG_SUPPORTED=true
        set "flag_bench=-Dbench=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --benchpp

    if "!FLAG!"=="--benchpp" (
        set FLAG_SUPPORTED=true
        set "flag_benchpp=-Dbenchpp=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --blas

    if "!FLAG!"=="--blas" (
        set FLAG_SUPPORTED=true
        set "flag_blas=-Dblas=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --codecov

    if "!FLAG!"=="--codecov" (
        set FLAG_SUPPORTED=true
        set "flag_codecov=-Dcodecov=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --deps

    if "!FLAG!"=="--deps" (
        set FLAG_SUPPORTED=true
        set "flag_deps=-Ddeps=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --exam

    if "!FLAG!"=="--exam" (
        set FLAG_SUPPORTED=true
        set "flag_exam=-Dexam=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --exampp

    if "!FLAG!"=="--exampp" (
        set FLAG_SUPPORTED=true
        set "flag_exampp=-Dexampp=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --cfi

    if "!FLAG!"=="--cfi" (
        set FLAG_SUPPORTED=true
        set "flag_cfi=-Dcfi=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --fpp

    if "!FLAG!"=="--fpp" (
        set FLAG_SUPPORTED=true
        set "flag_fpp=-Dfpp=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --fresh

    if "!FLAG!"=="--fresh" (
        set FLAG_SUPPORTED=true
        set "flag_fresh=-Dfresh=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --lapack

    if "!FLAG!"=="--lapack" (
        set FLAG_SUPPORTED=true
        set "flag_lapack=-Dlapack=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --matlabroot

    if "!FLAG!"=="--matlabroot" (
        set FLAG_SUPPORTED=true
        set "matlabroot=!VALUE!"
        set "flag_matlabroot=-Dmatlabroot="!VALUE!""
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --me

    if "!FLAG!"=="--me" (
        set FLAG_SUPPORTED=true
        set "flag_me=-Dme=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --mod

    if "!FLAG!"=="--mod" (
        set FLAG_SUPPORTED=true
        set "flag_mod=-Dmod=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --nproc

    if "!FLAG!"=="--nproc" (
        set FLAG_SUPPORTED=true
        set "flag_nproc=-Dnproc=!VALUE!"
        set FOR_COARRAY_NUM_IMAGES=!VALUE!
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --perfprof

    if "!FLAG!"=="--perfprof" (
        set FLAG_SUPPORTED=true
        set "flag_perfprof=-Dperfprof=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --pdt

    if "!FLAG!"=="--pdt" (
        set FLAG_SUPPORTED=true
        set "flag_pdt=-Dpdt=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --purity

    if "!FLAG!"=="--purity" (
        set FLAG_SUPPORTED=true
        set "flag_purity=-Dpurity=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --test

    if "!FLAG!"=="--test" (
        set FLAG_SUPPORTED=true
        set "flag_test=-Dtest=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM ::::::::::::::::::::
    REM flag args: type kind
    REM ::::::::::::::::::::

    REM --ski

    if "!FLAG!"=="--ski" (
        set FLAG_SUPPORTED=true
        set "flag_ski=-Dski=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --iki

    if "!FLAG!"=="--iki" (
        set FLAG_SUPPORTED=true
        set "flag_iki=-Diki=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --lki

    if "!FLAG!"=="--lki" (
        set FLAG_SUPPORTED=true
        set "flag_lki=-Dlki=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --cki

    if "!FLAG!"=="--cki" (
        set FLAG_SUPPORTED=true
        set "flag_cki=-Dcki=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --rki

    if "!FLAG!"=="--rki" (
        set FLAG_SUPPORTED=true
        set "flag_rki=-Drki=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM ::::::::::
    REM other args
    REM ::::::::::

    REM --ddir

    if "!FLAG!"=="--ddir" (
        set "ddir=!VALUE!"
        set FLAG_SUPPORTED=true
        set "flag_ddir=-Dddir="!VALUE!""
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --bdir

    if "!FLAG!"=="--bdir" (
        set FLAG_SUPPORTED=true
        set "bdir=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM --help

    if "!FLAG!"=="--help" (
        set FLAG_SUPPORTED=true
        type "!paramonte_dir!install.bat.md"
        type "!paramonte_dir!install.config.md"
        exit /b 0
    )

    REM --dev

    if "!FLAG!"=="--dev" (
        set FLAG_SUPPORTED=true
        set "flag_dev=-Ddev_enabled=1"
        set ntry=1
    )

    REM -G

    if "!FLAG!"=="-G" (
        set FLAG_SUPPORTED=true
        set "flag_G=-G "!VALUE!""
        set "cmakeBuildGenerator=!VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM -j

    if "!FLAG!"=="-j" (
        set FLAG_SUPPORTED=true
        set "flag_j=-j !VALUE!"
        if "!VALUE!"=="" set "VALUE_SUPPORTED=false"
        if /i "!VALUE:~0,2!"=="--" set "VALUE_SUPPORTED=false"
        shift
    )

    REM
    REM Check for errors.
    REM

    if  defined FLAG_SUPPORTED (
        if !VALUE_SUPPORTED! NEQ true goto LABEL_REPORT_ERR
        if !FLAG_SUPPORTED! NEQ true goto LABEL_REPORT_ERR
    )

    shift
    goto :LABEL_parseArgLoop

)

:LABEL_REPORT_ERR

REM Check flag/value support

if defined FLAG_SUPPORTED (
    if "!FLAG_SUPPORTED!"=="true" (
        if "!VALUE_SUPPORTED!" NEQ "true" (
            echo.
            echo.!pmfatal! The requested input value "!VALUE!" specified
            echo.!pmfatal! with the input flag "!FLAG!" is unsupported.
            goto LABEL_ERR
        )
    ) else (
        echo.
        echo.!pmfatal! The requested input flag "!FLAG!" is unsupported.
        goto LABEL_ERR
    )
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: Set the default values for the input command line arguments.
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if not defined list_build set list_build=release
if not defined list_checking set list_checking=nocheck
if not defined list_lang set list_lang=fortran
if not defined list_lib set list_lib=shared
if not defined list_mem set list_mem=heap
if not defined list_par set list_par=serial
if not defined flag_j set "flag_j=-j"

REM Set the optional values.

if defined flag_exam (
    if not defined flag_exampp (
        set "flag_exampp="
    )
)

if defined flag_bench (
    if not defined flag_benchpp (
        set "flag_benchpp="
    )
)

REM Set the default Fortran compiler and the `list_fc` flag.

if not defined list_fc (
    for %%C in (ifort ifx gfortran) do (
        echo.!pmnote! Checking for the presence of %%~C Fortran compiler...
        call :mktemp tempfile
        where %%~C > "!tempfile!"
        set /p list_fc=<"!tempfile!"
        if exist "!list_fc!" (
            echo.!pmnote! The %%~C Fortran compiler detected in the environment.
            echo.!pmnote! fc="!list_fc!"
            goto :loopExit
        ) else (
            set list_fc=
        )
    )
)

:loopExit
if not defined list_fc (
    echo.!pmwarn! No compatible Fortran compiler detected in the environment.
    echo.!pmwarn! Proceeding without a guarantee of build success...
    set "list_fc=default"
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: Set CMake default flags.
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

REM delete the prerequisite folder if requested.

set "substring=prereq"
if not "!flag_fresh!" == "" (
    if not "!flag_fresh:!substring!=!"=="!flag_fresh!" (
        if exist "!paramonte_req_dir!" (
            echo.!pmnote! Removing the old prerequisites of the ParaMonte library build at paramonte_req_dir="!paramonte_req_dir!"
            rmdir /S /Q "!paramonte_req_dir!"
        )
    )
)

if not defined flag_G (

    REM Set the default CMake makefile generator.

    set "replacement="
    set "cmakeBuildGenerator="

    REM Above all, search for the Ninja makefile generator: ninja
    REM The ninja executable is installed either as part of Microsoft Visual Studio or Quickstart Fortran software.

    if not defined cmakeBuildGenerator (
        echo.!pmnote! Searching for the Ninja build generator in the command-line environment...
        set "NINJA_FOUND="
        for %%X in (ninja.exe) do (set NINJA_FOUND=%%~$PATH:X)
        if defined NINJA_FOUND (
            echo.!pmnote! !BoldYellow!Setting CMake makefile generator to Ninja...!ColorReset!
            set "cmakeBuildGenerator=Ninja"
        ) else (
            echo.!pmnote! Failed to detect the Ninja build generator in the command-line environment. skipping...
        )
    )

    REM Firstly, search for the CMake makefile generator: make

    if not defined cmakeBuildGenerator (
        echo.!pmnote! Searching for the GNU Make application in the command-line environment...
        set "make_version="
        set "substring=GNU Make"
        for /f "Tokens=* Delims=" %%i in ('make --version') do set make_version=!make_version!%%i
        for /f "delims=" %%S in (^""!substring!=!replacement!"^") do (set "make_version_modified=!make_version:%%~S!")
        if not "!make_version_modified!" == "!make_version!" (
            echo.!pmnote! !BoldYellow!Setting CMake makefile generator to GNU MinGW Make application...!ColorReset!
            set "cmakeBuildGenerator=MinGW Makefiles"
        ) else (
            echo.!pmnote! Failed to detect the GNU Make application in the command-line environment. skipping...
        )
    )

    REM Secondly, search for the CMake makefile generator: mingw32-make

    if not defined cmakeBuildGenerator (
        echo.!pmnote! Searching for the GNU Make application in the command-line environment...
        set "make_version="
        set "substring=GNU Make"
        for /f "Tokens=* Delims=" %%i in ('mingw32-make --version') do set make_version=!make_version!%%i
        for /f "delims=" %%S in (^""!substring!=!replacement!"^") do (set "make_version_modified=!make_version:%%~S!")
        if not "!make_version_modified!" == "!make_version!" (
            echo.!pmnote! !BoldYellow!Setting CMake makefile generator to GNU MinGW Make application...!ColorReset!
            set "cmakeBuildGenerator=MinGW Makefiles"
        ) else (
            echo.!pmnote! Failed to detect the GNU MinGW Make application in the command-line environment. skipping...
        )
    )

    REM Thirdly, search for the CMake makefile generator: NMake

    if not defined cmakeBuildGenerator (
        echo.!pmnote! Searching for the Microsoft NMake application in the command-line environment...
        set "make_version="
        set "substring=Microsoft"
        for /f "Tokens=* Delims=" %%i in ('nmake') do set make_version=!make_version!%%i
        for /f "delims=" %%S in (^""!substring!=!replacement!"^") do (set "make_version_modified=!make_version:%%~S!")
        if not "!make_version_modified!" == "!make_version!" (
            echo.!pmnote! !BoldYellow!Setting CMake makefile generator to Microsoft NMake application...!ColorReset!
            set "cmakeBuildGenerator=NMake Makefiles"
        ) else (
            echo.!pmnote! Failed to detect the Microsoft NMake application in the command-line environment. skipping...
        )
    )

    REM Revert to the Windows CMD default if no CMake makefile generator is identified.

    if not defined cmakeBuildGenerator (
        echo.
        echo.!pmwarn! Failed to infer the CMake makefile generator.
        echo.!pmwarn! Procedding with Microsoft NMake as the default makefile generator...
        echo.!pmwarn! The CMake configuration and build may fail.
        set "cmakeBuildGenerator=NMake Makefiles"
    )

    set "flag_G=-G "!cmakeBuildGenerator!""


)

if  defined cmakeBuildGenerator (
    echo.!pmnote! cmakeBuildGenerator=!cmakeBuildGenerator!
    if "!cmakeBuildGenerator!" == "NMake Makefiles" set makename=nmake
    if "!cmakeBuildGenerator!" == "MinGW Makefiles" set makename=mingw
    if "!cmakeBuildGenerator!" == "Ninja" set makename=ninja
    echo.!pmnote! makename=!makename!
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: Build the library for all requested configurations.
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

REM Here we trust the user to assign semi-colon separated items as flag values.

for %%C in ("!list_fc:;=" "!") do (

    REM
    REM Set up the CMake fc flag.
    REM

    set "fc=%%~C"
    set "flag_fc="
    if exist "!fc!" (
        set "fcpath=!fc!"
        set "flag_fc=-Dfc="!fcpath!""
    ) else (
        if not "!fc!"=="default" (
            call :mktemp tempfile
            where "!fc!" > "!tempfile!"
            set /p fcpath=<"!tempfile!"
            if exist "!fcpath!" (
                echo.!pmnote! The !fc! Fortran compiler detected in the environment.
                echo.!pmnote! fc="!fcpath!"
            ) else (
                set "fcpath=!fc!"
            )
            if exist "!fcpath!" (
                set "flag_fc=-Dfc="!fcpath!""
                echo.!pmnote! Fortran compiler path fcpath="!fcpath!"
            ) else (
                echo.!pmwarn! Failed to detect the full path for the specified compiler fc=!fc!
                echo.!pmwarn! Proceeding with the build without guarantee of success...
            )
        )
    )

    REM Get the compiler ID.

    set csid=csid
    set "replacement="
    set "substring=intel"
    set "fcpath_lower=!fcpath!"
    call :getLowerCase fcpath_lower
    set "fcpath_modified=!fcpath_lower:%substring%=!"
    for /f "delims=" %%S in (^""!substring!=!replacement!"^") do (set "fcpath_modified=!fcpath_modified:%%~S!")
    if not "!fcpath_modified!"=="!fcpath_lower!" (
        echo.!pmnote! The Fortran compiler vendor is Intel.
        set csid=intel
    ) else (
        set "substring=gfortran"
        for /f "delims=" %%S in (^""!substring!=!replacement!"^") do (set "fcpath_modified=!fcpath_modified:%%~S!")
        if not "!fcpath_modified!"=="!fcpath_lower!" (
            echo.!pmnote! The Fortran compiler vendor is GNU.
            set csid=gnu
        ) else (
            echo.!pmwarn! "Failed to detect the specified Fortran compiler ID (vendor) from its path "!fcpath!""
            echo.!pmwarn! "Processing with CMake build configuration without guarantee of success..."
        )
    )

    REM Get the compiler version.

    set csvs=csvs
    if not !csid!==csid (
        REM Get unique file name.
        cd "!tmp!"
        set "tempsrc=!tmp!\getCompilerVersion.F90"
        call :mktemp tempout "!tmp!" "getCompilerVersion"
        set "tempexe=!tempout!.exe"
        echo F | xcopy /Y "!paramonte_auxil_dir!\getCompilerVersion.F90" "!tempsrc!"
        "!fcpath!" "!tempsrc!" -o "!tempexe!"
        "!tempexe!" "!tempout!.tmp" || (
            echo.!pmwarn! "Failed to infer the specified Fortran compiler version from executable "!fcpath!""
            echo.!pmwarn! "Processing with CMake build configuration without guarantee of success..."
        )
        cd "!paramonte_auxil_dir!"
        for /f "tokens=*" %%i in ('head.bat 1 "!tempout!.tmp"') do set "csvs=%%~i"
        cd %~dp0
    )

    echo.!pmnote! compiler suite !csid!
    echo.!pmnote! compiler version !csvs!

    for %%G in ("!list_lang:;=" "!") do (

        set "flag_lang=-Dlang=%%~G"

        for %%B in ("!list_build:;=" "!") do (

            set "flag_build=-Dbuild=%%~B"

            for %%L in ("!list_lib:;=" "!") do (

                set "flag_lib=-Dlib=%%~L"

                for %%M in ("!list_mem:;=" "!") do (

                    set "flag_mem=-Dmem=%%~M"

                    for %%P in ("!list_par:;=" "!") do (

                        set "flag_par=-Dpar=%%~P"

                        for %%H in ("!list_checking:;=" "!") do (

                            set "flag_checking=-Dchecking=%%~H"

                            REM
                            REM First, determine the parallelism and MPI library name to be used in build directory.
                            REM

                            if %%~P==mpi (
                                set parname=mpi
                                for /f %%i in ('mpiexec --version') do set mpiexec_version=%%i
                                echo !mpiexec_version! | find "Intel" >nul
                                if errorlevel 0 (
                                    echo.!pmnote! Intel MPI library detected...
                                    set parname=impi
                                ) else (
                                    echo.!pmwarn! Failed to infer the MPI library vendor and version.
                                    echo.!pmwarn! The CMake configuration and library build may fail. skipping...
                                )
                            ) else (
                                if %%~P==omp (
                                    set parname=openmp
                                ) else (
                                    if %%~P==none (
                                        set parname=serial
                                    ) else (
                                        set parname=%%~P
                                    )
                                )
                            )

                            REM
                            REM Set the ParaMonte CMake build directory.
                            REM

                            if  not defined bdir (
                                set "paramonte_bld_dir=!paramonte_dir!bld\!os!\!arch!\!csid!\!csvs!\%%~B\%%~L\%%~M\!parname!\%%~H\%%~G"
                                if "!flag_perfprof!" == "-Dperfprof=all" set "paramonte_bld_dir=!paramonte_bld_dir!\perfprof"
                                if "!flag_codecov!" == "-Dcodecov=all" set "paramonte_bld_dir=!paramonte_bld_dir!\codecov"
                                if defined makename set "paramonte_bld_dir=!paramonte_bld_dir!\!makename!"
                                echo.!pmnote! The ParaMonte library build directory paramonte_bld_dir="!paramonte_bld_dir!"
                            ) else (
                                echo.!pmnote! User-specified library build directory detected bdir="!bdir!"
                                set "paramonte_bld_dir=!bdir!"
                            )

                            REM Make the build directory if needed.

                            if  not exist "!paramonte_bld_dir!" (
                                echo.!pmnote! Generating the ParaMonte build directory...
                                mkdir "!paramonte_bld_dir!"
                            )

                            REM The following loop temporarily bypasses an existing bug where the first fresh installation
                            REM does not copy the FPP source files to the deployment and installation directories.

                            set "flag_fresh_current=!flag_fresh!"

                            for /l %%x in (1, 1, !ntry!) do (

                                REM
                                REM Configure and build the library via CMake.
                                REM

                                echo.!pmnote! All generated build files will be stored at "!paramonte_bld_dir!"
                                echo.!pmnote! Changing directory to "!paramonte_bld_dir!"
                                echo.
                                echo.****************************************************************************************************
                                echo.
                                echo.!pmnote! Invoking CMake as:
                                echo.

                                echo.cd "!paramonte_bld_dir!"
                                echo.cmake !paramonte_dir! !flag_G! -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON !flag_build! !flag_checking! !flag_lib! !flag_mem! !flag_par! !flag_fc!
                                echo.!flag_ddir! !flag_bench! !flag_benchpp! !flag_blas! !flag_codecov! !flag_cfi! !flag_deps! !flag_exam! !flag_exampp! !flag_fpp! !flag_fresh_current! !flag_lapack! !flag_matlabroot!
                                echo.!flag_lang! !flag_me! !flag_mod! !flag_nproc! !flag_perfprof! !flag_pdt! !flag_purity! !flag_test! !flag_ski! !flag_iki! !flag_lki! !flag_cki! !flag_rki! !flag_dev!

                                cd "!paramonte_bld_dir!"
                                cmake !paramonte_dir! !flag_G! -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON !flag_build! !flag_checking! !flag_lib! !flag_mem! !flag_par! !flag_fc! ^
                                !flag_ddir! !flag_bench! !flag_benchpp! !flag_blas! !flag_codecov! !flag_cfi! !flag_deps! !flag_exam! !flag_exampp! !flag_fpp! !flag_fresh_current! !flag_lapack! !flag_matlabroot! ^
                                !flag_lang! !flag_me! !flag_mod! !flag_nproc! !flag_perfprof! !flag_pdt! !flag_purity! !flag_test! !flag_ski! !flag_iki! !flag_lki! !flag_cki! !flag_rki! !flag_dev! ^
                                && (
                                    echo.
                                    echo.!pmnote! !BoldGreen!ParaMonte configuration with CMake appears to have succeeded.!ColorReset!
                                ) || (
                                    echo.
                                    echo.!pmfatal! !BoldRed!ParaMonte configuration with CMake appears to have failed.!ColorReset!
                                    echo.!pmfatal! !BoldRed!This error could happen for a variety of silly reasons such Dropbox interference with the build files.!ColorReset!
                                    echo.!pmfatal! !BoldRed!Make sure to disable intensive applications such as Google Drive and Dropbox which lock file ownerships.!ColorReset!
                                    echo.!pmfatal! !BoldRed!Then retry the build. Keep in mind that this is only one probable cause of the build failure.!ColorReset!
                                    goto LABEL_ERR
                                )

                                REM
                                REM Reset the fresh flag to ensure the build is not erased during the second CMake configure cycle.
                                REM

                                set "flag_fresh_current=-Dfresh=none"

                                echo.
                                echo.****************************************************************************************************
                                echo.

                                set PATH=!PATH!;!paramonte_bld_dir!\lib
                                cd "!paramonte_bld_dir!" && cmake --build "!paramonte_bld_dir!" !flag_j! && (
                                    echo.
                                    echo.!pmnote! !BoldGreen!ParaMonte build appears to have succeeded.!ColorReset!
                                    echo.
                                ) || (
                                    echo.
                                    echo.!pmnote! !BoldRed!ParaMonte build appears to have failed.!ColorReset!
                                    echo.
                                    goto LABEL_ERR
                                )

                                cd "!paramonte_bld_dir!" && cmake --build "!paramonte_bld_dir!" --target install !flag_j! && (
                                    echo.
                                    echo.!pmnote! !BoldGreen!ParaMonte installation appears to have succeeded.!ColorReset!
                                    echo.
                                ) || (
                                    echo.
                                    echo.!pmnote! !BoldRed!ParaMonte installation appears to have failed.!ColorReset!
                                    echo.
                                    goto LABEL_ERR
                                )

                                REM
                                REM  Search for MATLAB installations and build MEX files before deploying the package.
                                REM

                                if %%~G==matlab (

                                    set "MATLAB_ROOT_DIR="
                                    set "MATLAB_EXE_PATH="
                                    set "MATLAB_BIN_DIR="
                                    set "MATLAB_LIB_DIR="
                                    set "MATLAB_INC_DIR=."
                                    set "MATLAB_LIBMX_FILE="
                                    set "MATLAB_LIBMEX_FILE="
                                    set "MATLAB_LIBMAT_FILE="
                                    set "MATLAB_VERSION_FILE="
                                    REM set "MATLAB_INC_DIR_FLAG="

                                    echo.
                                    echo.!pmnote! Searching for a MATLAB installations on your system...

                                    set "INSTALL_LOC_LIST=C:\Program Files\MATLAB\/C:\Program Files (x86)\MATLAB\"
                                    set MATLAB_VERSION_LIST=R2035b/R2035a/R2034b/R2034a/R2033b/R2033a/R2032b/R2032a/R2031b/R2031a/R2030b/R2030a/R2029b/R2029a/R2028b/R2028a/R2027b/R2027a/R2026b/R2026a
                                    set MATLAB_VERSION_LIST=!MATLAB_VERSION_LIST!/R2025b/R2025a/R2024b/R2024a/R2023b/R2023a/R2022b/R2022a/R2021b/R2021a/R2020b/R2020a/R2019b/R2019a/R2018b/R2018a/R2017b/R2017a

                                    REM
                                    REM Amir Shahmoradi Oct 25, 2024:
                                    REM The following block is currently was added despite its functionality being already implemented within CMake.
                                    REM The reason for its existence is to resolve the vicious bug that exists in CMake intrinsic module FindMatlab.cmake yielding the following runtime error:
                                    REM
                                    REM     Error using pm.sampling.Sampler/run MATLAB:mex:ErrInvalidMEXFile : Invalid MEX-file 'pm_sampling.mexw64': Gateway function is missing
                                    REM
                                    REM See also,
                                    REM
                                    REM     https://gitlab.kitware.com/cmake/cmake/-/issues/25068#note_1580985
                                    REM
                                    REM for a relevant discussion of this bug faced by others and the status of a resolution to fix it.
                                    REM Note that this CMake bug is different from another vicious MATLAB-MEX-version related bug that causes the MEX files to fail at runtime
                                    REM while the same MEX compilation and run for ParaMonte 1 succeeds with MATLAB R2022b and older.
                                    REM See
                                    REM
                                    REM     https://www.mathworks.com/matlabcentral/answers/2157360-matlab-mex-errinvalidmexfile-invalid-mex-file-the-specified-procedure-could-not-be-found?s_tid=prof_contriblnk
                                    REM
                                    REM for more relevant discussion of this bug and possible causes.
                                    REM
                                    REM As of today, both CMake and MATLAB MEX compatibility bugs remain unresolved.
                                    REM The following block can be commented out by setting the value of
                                    REM `MATLAB_FOUND` to `none` in the following `set` command.
                                    REM
                                    REM \todo
                                    REM \pvhigh
                                    REM Once the CMake bug in FindMatlab.cmake intrinsic modules is resolved, the whole shenanigan above and below for MEX compilation must be removed.
                                    REM

                                    set MATLAB_FOUND=false
                                    for %%D in ("!INSTALL_LOC_LIST:/=" "!") do (
                                        for %%V in ("!MATLAB_VERSION_LIST:/=" "!") do (

                                            if !MATLAB_FOUND!==false (

                                                if  defined matlabroot (
                                                    set "MATLAB_ROOT_DIR_TEMP=!matlabroot!"
                                                    echo.!pmnote! !BoldYellow!Searching for user-specified MATLAB installation at: !MATLAB_ROOT_DIR_TEMP! !ColorReset!
                                                ) else (
                                                    set "MATLAB_ROOT_DIR_TEMP=%%~D%%~V"
                                                )
                                                set "MATLAB_BIN_DIR_TEMP=!MATLAB_ROOT_DIR_TEMP!\bin"
                                                set "MATLAB_EXE_PATH_TEMP=!MATLAB_BIN_DIR_TEMP!\matlab.exe"

                                                if  exist !MATLAB_EXE_PATH_TEMP! (

                                                    set MATLAB_FOUND=true
                                                    set "MATLAB_ROOT_DIR=!MATLAB_ROOT_DIR_TEMP!"
                                                    set "MATLAB_EXE_PATH=!MATLAB_EXE_PATH_TEMP!"
                                                    set "MATLAB_BIN_DIR=!MATLAB_BIN_DIR_TEMP!"
                                                    set "MATLAB_INC_DIR=!MATLAB_ROOT_DIR!\extern\include"
                                                    set "MATLAB_LIB_DIR=!MATLAB_ROOT_DIR!\extern\lib\win64\microsoft"
                                                    set "MATLAB_LIBMX_FILE=!MATLAB_LIB_DIR!\libmx.lib"
                                                    set "MATLAB_LIBMEX_FILE=!MATLAB_LIB_DIR!\libmex.lib"
                                                    set "MATLAB_LIBMAT_FILE=!MATLAB_LIB_DIR!\libmat.lib"
                                                    set "MATLAB_VERSION_FILE=!MATLAB_ROOT_DIR!\extern\version\fortran_mexapi_version.F"
                                                    echo.!pmnote! !BoldYellow!MATLAB installation detected at: !MATLAB_EXE_PATH! !ColorReset!
                                                    REM set "MATLAB_INC_DIR_FLAG=/I:!MATLAB_INC_DIR!"
                                                    REM set FPP_FLAGS=/define:MATLAB_MEX_FILE

                                                    REM
                                                    REM  Build MATLAB MEX files.
                                                    REM

                                                    call set PMLIB_MATLAB_NAME=!PMLIB_NAME:_matlab_=_!

                                                    REM /subsystem:windows
                                                    REM mex -setup:"C:\Program Files\MATLAB\R2019a\bin\win64\mexopts\intel_c_19_vs2017.xml" C
                                                    REM if !BTYPE!==debug   set "MATLAB_BUILD_FLAGS=!MATLAB_BUILD_FLAGS!/Od /Z7"
                                                    REM if !BTYPE!==testing set "MATLAB_BUILD_FLAGS=!MATLAB_BUILD_FLAGS!/O2"
                                                    REM if !BTYPE!==release set "MATLAB_BUILD_FLAGS=!MATLAB_BUILD_FLAGS!/Od"

                                                    REM if !BTYPE!==debug   set "MATLAB_BUILD_FLAGS=!MATLAB_BUILD_FLAGS!!INTEL_CPP_DEBUG_FLAGS!"
                                                    REM if !BTYPE!==testing set "MATLAB_BUILD_FLAGS=!MATLAB_BUILD_FLAGS!!INTEL_CPP_TESTING_FLAGS!"
                                                    REM if !BTYPE!==release set "MATLAB_BUILD_FLAGS=!MATLAB_BUILD_FLAGS!!INTEL_CPP_RELEASE_FLAGS!"

                                                    set "MEX_FLAGS=-v -nojvm"

                                                    REM
                                                    REM If openmp is enabled, define the macro OMP_ENABLED=1.
                                                    REM \todo
                                                    REM \pvhigh
                                                    REM This is a weakness point as the input value for `--par` flag may not be completely lower case.
                                                    REM

                                                    REM ;pm_parallelism
                                                    set "list_mex=pm_sampling"

                                                    set "ismatlabomp=false"
                                                    if omp==%%~P set "ismatlabomp=true"
                                                    if OMP==%%~P set "ismatlabomp=true"
                                                    if openmp==%%~P set "ismatlabomp=true"
                                                    if OPENMP==%%~P set "ismatlabomp=true"
                                                    if !ismatlabomp!==true set "MEX_FLAGS=!MEX_FLAGS! -DOMP_ENABLED"

                                                    echo.!pmnote!!BoldYellow!Generating the ParaMonte MATLAB MEX files...!ColorReset!

                                                    for %%X in ("!list_mex:;=" "!") do (
                                                        echo.!pmnote!!BoldYellow!Compiler command: "!MATLAB_BIN_DIR!\mex.bat" !MEX_FLAGS! "!paramonte_src_dir!\matlab\xrc\pm_sampling.c" libparamonte.lib -output pm_sampling!ColorReset!
                                                        cd !paramonte_bld_dir!\lib
                                                        REM we cannot use the version variable when MATLAB directory is user-specified.
                                                        REM if not exist "%%~V" (mkdir "%%~V")
                                                        REM cd %%~V
                                                        call "!MATLAB_BIN_DIR!\mex.bat" !MEX_FLAGS! "!paramonte_src_dir!\matlab\xrc\pm_sampling.c" libparamonte.lib -output pm_sampling && (
                                                            echo.!pmnote! !BoldGreen!The ParaMonte MATLAB shared library build appears to have succeeded.!ColorReset!
                                                        ) || (
                                                            echo.
                                                            echo.!pmwarn! !BoldMagenta!The ParaMonte MATLAB library build failed.!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!Please make sure you have the following components installed!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!on your system before rerunning the installation script:!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!    -- MATLAB, including MATLAB MEX compilers.!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!    -- Intel OneAPI icx/icl and ifx/ifort compilers 2023 or newer.!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!Once you are sure of the existence of these components in your !ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!Windows command line environment, run the following command:!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!    "!MATLAB_BIN_DIR!\mex.bat" -setup C!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!Among the options displayed, you should see the command to setup!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!the Intel OneAPI icl/icx or Microsoft cl compiler for C on your system.!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!This command should look similar to the following,!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!    "!MATLAB_BIN_DIR_TEMP!\mex.bat" -setup:"C:\Program Files\MATLAB\R2024a\bin\win64\mexopts\intel_c_24_vs2022.xml" C!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!with minor differences in the xml file name depending on your specific installations of !ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!    -- the Intel OneAPI or Microsoft compiler version!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!    -- the Microsoft Visual Studio version!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!    -- the MATLAB version!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!Copy and paste this command into your terminal, run it, and then rerun the ParaMonte MATLAB installation script.!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!Please report this or any other issues at:!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!!ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!    https://github.com/cdslaborg/paramonte/issues !ColorReset!
                                                            echo.!pmwarn! !BoldMagenta!!ColorReset!
                                                            echo.
                                                            REM set ERRORLEVEL=1
                                                            REM exit /B 1
                                                        )
                                                    )
                                                    cd %~dp0

                                                )

                                                set "MATLAB_ROOT_DIR_TEMP="
                                                set "MATLAB_BIN_DIR_TEMP="
                                                set "MATLAB_EXE_PATH_TEMP="

                                            )
                                        )
                                    )

                                    if  MATLAB_FOUND==false (
                                        echo.!pmwarn! !BoldMagenta!Exhausted all possible search paths for a MATLAB installation, but failed to find MATLAB.!ColorReset!
                                        echo.!pmwarn! !BoldMagenta!The ParaMonte MATLAB kernel will not be functional without building the required DLL libraries.!ColorReset!
                                        echo.!pmwarn! !BoldMagenta!Please add MATLAB to your environmental variable PATH and rerun the install script.!ColorReset!
                                        echo.!pmwarn! !BoldMagenta!For example, on your current Windows command-line, try:!ColorReset!
                                        echo.!pmwarn! !BoldMagenta!!ColorReset!
                                        echo.!pmwarn! !BoldMagenta!    set "PATH=PATH_TO_MATLAB_BIN_DIR;!PATH!
                                        echo.!pmwarn! !BoldMagenta!!ColorReset!
                                        echo.!pmwarn! !BoldMagenta!where PATH_TO_MATLAB_BIN_DIR must be replaced with path to the bin folder of the current!ColorReset!
                                        echo.!pmwarn! !BoldMagenta!installation of MATLAB on your system. Typical MATLAB bin installation path on a 64-bit Windows!ColorReset!
                                        echo.!pmwarn! !BoldMagenta!Operating Systems is a string like the following:!ColorReset!
                                        echo.!pmwarn! !BoldMagenta!!ColorReset!
                                        echo.!pmwarn! !BoldMagenta!    C:\Program Files\MATLAB\2020a\bin\!ColorReset!
                                        echo.!pmwarn! !BoldMagenta!!ColorReset!
                                        echo.!pmwarn! !BoldMagenta!where 2020a in the path points to the MATLAB 2020a version installation on the system. You can also!ColorReset!
                                        echo.!pmwarn! !BoldMagenta!find the installation location of MATLAB by typing the following command in your MATLAB session:!ColorReset!
                                        echo.!pmwarn! !BoldMagenta!!ColorReset!
                                        echo.!pmwarn! !BoldMagenta!    matlabroot!ColorReset!
                                        echo.!pmwarn! !BoldMagenta!
                                        echo.!pmwarn! !BoldMagenta!skipping the ParaMonte MATLAB build...!ColorReset!
                                    )

                                )

                                REM
                                REM End of MATLAB MEX build.
                                REM

                                cd "!paramonte_bld_dir!" && cmake --build "!paramonte_bld_dir!" --target deploy !flag_j! && (
                                    echo.
                                    echo.!pmnote! !BoldGreen!ParaMonte deploy appears to have succeeded.!ColorReset!
                                    echo.
                                ) || (
                                    echo.
                                    echo.!pmnote! !BoldRed!ParaMonte deploy appears to have failed.!ColorReset!
                                    echo.
                                    goto LABEL_ERR
                                )

                            )

                            cd "!paramonte_bld_dir!" && cmake --build "!paramonte_bld_dir!" --target test && (
                                echo.
                                echo.!pmnote! !BoldGreen!ParaMonte test appears to have succeeded.!ColorReset!
                                echo.
                            ) || (
                                echo.
                                echo.!pmnote! !BoldRed!ParaMonte test appears to have failed.!ColorReset!
                                echo.
                                goto LABEL_ERR
                            )

                            cd "!paramonte_bld_dir!" && cmake --build "!paramonte_bld_dir!" --target example && (
                                echo.
                                echo.!pmnote! !BoldGreen!ParaMonte example appears to have succeeded.!ColorReset!
                                echo.
                            ) || (
                                echo.
                                echo.!pmnote! !BoldRed!ParaMonte example appears to have failed.!ColorReset!
                                echo.
                                goto LABEL_ERR
                            )

                            cd "!paramonte_bld_dir!" && cmake --build "!paramonte_bld_dir!" --target benchmark && (
                                echo.
                                echo.!pmnote! !BoldGreen!ParaMonte benchmark appears to have succeeded.!ColorReset!
                                echo.
                            ) || (
                                echo.
                                echo.!pmnote! !BoldRed!ParaMonte benchmark appears to have failed.!ColorReset!
                                echo.
                                goto LABEL_ERR
                            )

                            REM
                            REM Mission Accomplished.
                            REM

                            echo.
                            echo.!pmnote! !BoldGreen!All build files for the current build configurations are stored at!ColorReset! "!paramonte_bld_dir!"
                            echo.

                        )
                    )
                )
            )
        )
    )

)

echo.
echo.!pmnote! !BoldGreen!All build files for all requested build configurations are stored at!ColorReset! "!paramonte_bld_dir!"
echo.!pmnote! !BoldGreen!The installed binary files for all requested build configurations are ready to use at!ColorReset! "!ddir!"
echo.

REM
REM zip the binary folder. The application tar.exe
REM

set "zipperFound="
set zipperName=tar.exe
for %%X in (!zipperName!) do (set zipperFound=%%~$PATH:X)
if  "!zipperFound!"=="" (
    echo.
    echo.!pmwarn! !BoldMagenta!Skipping the binary archive generation because the !zipperName! application could not be found.!ColorReset!
    echo.
) else (
    echo.
    echo.!pmnote! !BoldGreen!Generating the binary archive zip file using !zipperName! at:!ColorReset! "!ddir!"
    echo.
    call :NORMALIZEPATH "!ddir!"
    if  exist "!ddir!" (
        cd "!ddir!"
        echo.
        echo. -- ParaMonte - compressing all subdirectories in the directory: "!ddir!"
        echo.
        for /f "tokens=* usebackq" %%G in (`dir /b /a:d "!ddir!"`) do (
            if exist "%%~G.zip" (
                echo.!pmwarn! !BoldMagenta!: compressed subdirectory already exists:!ColorReset! "!ddir!\%%~G.zip"
                echo.!pmwarn! !BoldMagenta!: overwriting the existing archive file...!ColorReset!
            )
            echo. -- ParaMonte - compressing subdirectory: "!ddir!%%~G"
            tar.exe -a -cf "%%~G.zip" "%%~G" || (
                echo.
                echo.!pmfatal! !BoldRed!: compression failed for subdirectory:!ColorReset! "!ddir!\%%~G"
                echo.!pmfatal! !BoldRed!: gracefully exiting.!ColorReset!
                echo.
                cd "!paramonte_dir!"
                set ERRORLEVEL=1
                exit /B 1
            )
        )
    ) else (
        echo.
        echo.!pmfatal! !BoldRed!: The final binary deployment destination directory does not exist: "!ddir!"
        echo.
        cd "!paramonte_dir!"
        set ERRORLEVEL=1
        exit /B 1
    )
)

goto LABEL_EOF

:: subroutines

:mktemp
REM Get unique random file name: mktemp tempfile mktemp_dir prefix suffix
REM where tempfile is the output and mktemp_dir and prefix and suffix are three optional input arguments.
if "%~2" == "" (
    set mktemp_dir=!tmp!
) else (
    set "mktemp_dir=%~2"
)
:loopUniq
set "%1=!mktemp_dir!\%~3!RANDOM!%~4"
if exist "%~3" goto :loopUniq
GOTO:EOF

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
echo.!pmnote! To see the list of possible flags and associated values, try:
echo.!pmnote!
echo.!pmnote!     install.bat --help
echo.!pmnote!
echo.!pmnote! gracefully exiting the !script_name! script.
echo.

exit /B 1

:LABEL_copyErrorOccured

echo.
echo.!pmfatal! Failed to copy contents. exiting...
echo.
cd %~dp0
set ERRORLEVEL=1
exit /B 1

:LABEL_delErrorOccured

echo.
echo.!pmfatal! Failed to delete contents. exiting...
echo.
cd %~dp0
set ERRORLEVEL=1
exit /B 1

:NORMALIZEPATH
cd "!paramonte_dir!"
set DESTINATION_DIR=%~dpfn1
exit /B

:LABEL_EOF

echo.
echo.!pmnote! !BoldGreen!mission accomplished.!ColorReset!
echo.

exit /B 0
