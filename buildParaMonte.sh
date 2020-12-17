#!/bin/bash
####################################################################################################################################
####################################################################################################################################
####
####   MIT License
####
####   ParaMonte: plain powerful parallel Monte Carlo library.
####
####   Copyright (C) 2012-present, The Computational Data Science Lab
####
####   This file is part of the ParaMonte library.
####
####   Permission is hereby granted, free of charge, to any person obtaining a
####   copy of this software and associated documentation files (the "Software"),
####   to deal in the Software without restriction, including without limitation
####   the rights to use, copy, modify, merge, publish, distribute, sublicense,
####   and/or sell copies of the Software, and to permit persons to whom the
####   Software is furnished to do so, subject to the following conditions:
####
####   The above copyright notice and this permission notice shall be
####   included in all copies or substantial portions of the Software.
####
####   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
####   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
####   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
####   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
####   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
####   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
####   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
####
####   ACKNOWLEDGMENT
####
####   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
####   As per the ParaMonte library license agreement terms, if you use any parts of
####   this library for any purposes, kindly acknowledge the use of ParaMonte in your
####   work (education/research/industry/development/...) by citing the ParaMonte
####   library as described on this page:
####
####       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
####
####################################################################################################################################
####################################################################################################################################
#
# NOTE: Do not change the contents of this file unless you know what the consequences are.
# This is the Bash script file that builds objects, dynamic libraries,
# as well as the test and example binaries of the ParaMonte library on non-Windows systems.
# Upon invocation of this file from a Bash command-line interface,
# this file will first call the configuration file configParaMonte.bat to read the user's
# requested configuration for building the ParaMonte library.

####################################################################################################################################
#### set minimum requirements
####################################################################################################################################

openCoarraysVersion="2.9.0"
gnuVersionOpenCoarrays="10.1.0"
mpichVersionOpenCoarrays="3.2"
gnuVersionParaMonteCompatible="8.4.0"
cmakeVersionParaMonteCompatible="3.14.0"
intelVersionParaMonteCompatible="18.0.0"

####################################################################################################################################
#### set up color coding
####################################################################################################################################

asciiEscVal=27; Esc="$(printf "\\$(printf "%o" "${asciiEscVal}")")";
ColorReset="${Esc}[m"
ColorBold="${Esc}[1m"
Red="${Esc}[31m"
Green="${Esc}[32m"
Yellow="${Esc}[33m"
Blue="${Esc}[34m"
Magenta="${Esc}[35m"
Cyan="${Esc}[36m"
White="${Esc}[37m"
BoldRed="${Esc}[1;31m"
BoldGreen="${Esc}[1;32m"
BoldYellow="${Esc}[1;33m"
BoldBlue="${Esc}[1;34m"
BoldMagenta="${Esc}[1;35m"
BoldCyan="${Esc}[1;36m"
BoldWhite="${Esc}[1;37m"

####################################################################################################################################
#### set up the main root paths
####################################################################################################################################

BUILD_NAME="ParaMonte"; export BUILD_NAME

pmcolor="${BoldCyan}"
pmattn=" ${pmcolor}-- ${BUILD_NAME} -${ColorReset}"
pmnote="${pmattn} ${BoldYellow}NOTE:${ColorReset}"
pmwarn="${pmattn} ${BoldMagenta}WARNING:${ColorReset}"
pmfatal="${pmattn} ${BoldRed}FATAL:${ColorReset}"
warning="${BoldMagenta}WARNING${ColorReset}"

workingDir="$(pwd)"
sourceFileDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
ParaMonte_ROOT_DIR="${sourceFileDir}"; export ParaMonte_ROOT_DIR
ParaMonteKernel_SRC_DIR="${ParaMonte_ROOT_DIR}/src/kernel"; export ParaMonteKernel_SRC_DIR
ParaMonteKernelVersion_SRC_INC_PATH="${ParaMonteKernel_SRC_DIR}/ParaMonte_mod@version@kernel.inc.f90"; export ParaMonteKernel_SRC_DIR
ParaMonteTest_SRC_DIR="${ParaMonte_ROOT_DIR}/src/kernel/tests"; export ParaMonteTest_SRC_DIR
#export ParaMonte_ROOT_DIR="${ParaMonte_ROOT_DIR:-${PWD%/}}"

ParaMonte_ROOT_BUILD_DIR="${ParaMonte_ROOT_DIR}/build"; export ParaMonte_ROOT_BUILD_DIR
if ! [ -d "${ParaMonte_ROOT_BUILD_DIR}" ]; then
    mkdir -p "${ParaMonte_ROOT_BUILD_DIR}"
fi

if [[ ! -f "$(pwd)/build${BUILD_NAME}.sh" ]]; then
    echo >&2
    echo >&2 "${pmfatal} Build failed."
    echo >&2 "${pmfatal} Please run this script inside the top-level ParaMonte library root directory."
    echo >&2 "${pmfatal} This is the directory which contains this file in the GitHub repository of ParaMonte."
    echo >&2
    exit 1
fi

####################################################################################################################################
#### utils
####################################################################################################################################

# Compare two version strings [$1: version string 1 (v1), $2: version string 2 (v2)]
# Return values:
#   0: v1 == v2
#   1: v1 > v2
#   2: v1 < v2
function compareVersions() {

    # Trivial v1 == v2 test based on string comparison
    [[ "$1" == "$2" ]] && return 0

    # Local variables
    local regex="^(.*)-r([0-9]*)$" va1=() vr1=0 va2=() vr2=0 len i IFS="."

    # Split version strings into arrays, extract trailing revisions
    if [[ "$1" =~ ${regex} ]]; then
        va1=(${BASH_REMATCH[1]})
        [[ -n "${BASH_REMATCH[2]}" ]] && vr1=${BASH_REMATCH[2]}
    else
        va1=($1)
    fi
    if [[ "$2" =~ ${regex} ]]; then
        va2=(${BASH_REMATCH[1]})
        [[ -n "${BASH_REMATCH[2]}" ]] && vr2=${BASH_REMATCH[2]}
    else
        va2=($2)
    fi

    # Bring va1 and va2 to same length by filling empty fields with zeros
    (( ${#va1[@]} > ${#va2[@]} )) && len=${#va1[@]} || len=${#va2[@]}
    for ((i=0; i < len; ++i)); do
        [[ -z "${va1[i]}" ]] && va1[i]="0"
        [[ -z "${va2[i]}" ]] && va2[i]="0"
    done

    # Append revisions, increment length
    va1+=($vr1)
    va2+=($vr2)
    len=$((len+1))

    # *** DEBUG ***
    #echo "TEST: '${va1[@]} (?) ${va2[@]}'"

    # Compare version elements, check if v1 > v2 or v1 < v2
    for ((i=0; i < len; ++i)); do
        if (( 10#${va1[i]} > 10#${va2[i]} )); then
            return 1
        elif (( 10#${va1[i]} < 10#${va2[i]} )); then
            return 2
        fi
    done

    # All elements are equal, thus v1 == v2
    return 0
}

verify() {
    if [ $1 -eq 0 ]; then
        echo >&2 "-- ${pmattn} ${BoldGreen}The ParaMonte $2 appears to have succeeded.${ColorReset}"
    else
        echo >&2
        echo >&2 "    -- ParaMonte - FATAL: ParaMonte $2 appears to have failed."
        echo >&2 "    -- ParaMonte - FATAL: If the source of the error cannot be identified,"
        echo >&2 "    -- ParaMonte - FATAL: consider a fresh installation of ParaMonte's required compilers by calling"
        echo >&2 "    -- ParaMonte - FATAL: "
        echo >&2 "    -- ParaMonte - FATAL:     ./install.sh --fresh"
        echo >&2 "    -- ParaMonte - FATAL: "
        echo >&2 "    -- ParaMonte - FATAL: If the error happens during the installation of ParaMonte prerequisites,"
        echo >&2 "    -- ParaMonte - FATAL: it is possible that the current existing GNU compiler collection installed"
        echo >&2 "    -- ParaMonte - FATAL: on your system cannot compile the downloaded version of GNU that is required"
        echo >&2 "    -- ParaMonte - FATAL: for ParaMonte build. In such case, make sure you have a GNU compiler collection"
        echo >&2 "    -- ParaMonte - FATAL: version ${gnuVersionParaMonteCompatible} or newer installed on your system, with an updated PATH environmental"
        echo >&2 "    -- ParaMonte - FATAL: variable, then reinstall ParaMonte."
        echo >&2 "    -- ParaMonte - FATAL: "
        echo >&2 "    -- ParaMonte - FATAL: If the error is solely due to the failures of some ParaMonte tests, then"
        echo >&2 "    -- ParaMonte - FATAL: you may want to skip the testing of the library by specifying \"-t none\" or"
        echo >&2 "    -- ParaMonte - FATAL: \"--test none\" when calling the ParaMonte installation script."
        echo >&2 "    -- ParaMonte - FATAL: "
        echo >&2 "    -- ParaMonte - FATAL:     ./install.sh --test none"
        echo >&2 "    -- ParaMonte - FATAL: "
        echo >&2 "    -- ParaMonte - FATAL: or,"
        echo >&2 "    -- ParaMonte - FATAL: "
        echo >&2 "    -- ParaMonte - FATAL:     ./install.sh --t none"
        echo >&2 "    -- ParaMonte - FATAL: "
        echo >&2 "    -- ParaMonte - FATAL: To get more help on the usage of the optional install flags, try:"
        echo >&2 "    -- ParaMonte - FATAL: "
        echo >&2 "    -- ParaMonte - FATAL:     ./install.sh --help"
        echo >&2 "    -- ParaMonte - FATAL: "
        echo >&2 "    -- ParaMonte - FATAL: If all ParaMonte installation attempts fail, please report this issue at"
        echo >&2 "    -- ParaMonte - FATAL: "
        echo >&2 "    -- ParaMonte - FATAL:     https://github.com/shahmoradi/paramonte/issues"
        echo >&2 "    -- ParaMonte - FATAL: "
        echo >&2 "    -- ParaMonte - FATAL: or by contacting the ParaMonte authors directly (e.g., shahmoradi@utexas.edu)."
        echo >&2
        echo >&2
        echo >&2 "    ${pmattn} gracefully exiting."
        echo >&2
        exit 1
    fi
}

####################################################################################################################################
#### determine the platform
####################################################################################################################################

UNAME_PLATFORM="$(uname -s)"
case "${UNAME_PLATFORM}" in
    Linux*)     PLATFORM=linux;;
    Darwin*)    PLATFORM=darwin;;
    CYGWIN*)    PLATFORM=cygwin;;
    MINGW*)     PLATFORM=mingw;;
    *)          PLATFORM="unknown:${UNAME_PLATFORM}"
esac
if [[ "$PLATFORM" =~ .*"unknown".* ]]; then
    echo >&2
    echo >&2 "-- ${BUILD_NAME} - FATAL: Build failed. unrecognized platform - ${PLATFORM}"
    echo >&2 "-- ${BUILD_NAME} - FATAL: The supported platforms include: Linux, Darwin, CYGWIN, MINGW"
    echo >&2 "-- ${BUILD_NAME} - FATAL: The ParaMonte build has been only tested on Linux and Darwin platforms."
    echo >&2
    echo >&2 "-- ${BUILD_NAME} - gracefully exiting."
    echo >&2
    exit 1
else
    export PLATFORM
fi

UNAME_PLATFORM_FULL="$(uname -a)"
if [[ "$UNAME_PLATFORM_FULL" =~ .*"Microsoft".* && "$UNAME_PLATFORM_FULL" =~ .*"Linux".* ]]; then
    isWSL=true
else
    isWSL=false
fi

isMacOS=false
isLinux=false
isMingw=false
isCygwin=false
isWindows=false
if [[ "${UNAME_PLATFORM}" =~ .*"Darwin".* ]]; then
    isMacOS=true
    OSNAME="macOS"
elif [[ "${UNAME_PLATFORM}" =~ .*"Linux".* ]]; then
    isLinux=true
    OSNAME="Linux"
elif [[ "${UNAME_PLATFORM}" =~ .*"CYGWIN".* ]]; then
    OSNAME="Cygwin"
    isWindows=true
    isCygwin=true
elif [[ "${UNAME_PLATFORM}" =~ .*"MINGW".* ]]; then
    OSNAME="MinGW"
    isWindows=true
    isMingw=true
fi
export isMacOS
export isLinux
export isWindows
export isCygwin
export isMingw
export OSNAME
export isWSL

####################################################################################################################################
#### determine the architecture
####################################################################################################################################

ARCHITECTURE=$(uname -p)
if [[ "$ARCHITECTURE" =~ .*"64".* ]]; then
    ARCHITECTURE="x64"
else
    ARCHITECTURE=$(uname -m)
    if [[ "$ARCHITECTURE" =~ .*"64".* ]]; then ARCHITECTURE="x64"; fi
fi
export ARCHITECTURE

####################################################################################################################################
#### set the ParaMonte version (to be used by cmake)
####################################################################################################################################

read -r ParaMonteVersion < .VERSION

unset fppParaMonteVersion
unset FPP_PARAMONTE_VERSION_FLAG
if ! [ -z ${ParaMonteVersion+x} ]; then

    export ParaMonteVersion

    # NOTE: uncomment the following line, if you wish to insert the ParaMonte version to the source code via Cmake via the compiler preprocessor.
    # fppParaMonteVersion="${ParaMonteVersion}"; export fppParaMonteVersion

    # write the version source include file.

    if ! [ -f "${ParaMonteKernelVersion_SRC_INC_PATH}" ] || [ "$(grep -c "${ParaMonteVersion}" "${ParaMonteKernelVersion_SRC_INC_PATH}")" = "0" ]; then

        echo "! WARNING - DO NOT CHANGE THE CONTENTS OF THIS FILE MANUALLY." > "${ParaMonteKernelVersion_SRC_INC_PATH}"
        echo "! WARNING - This file is auto-generated by the ParaMonte build scripts." >> "${ParaMonteKernelVersion_SRC_INC_PATH}"
        echo "self%version = \"${ParaMonteVersion}\"" >> "${ParaMonteKernelVersion_SRC_INC_PATH}"

    fi

fi

#echo "$(cat ./auxil/.ParaMonteBanner)"

####################################################################################################################################
#### report build setup
####################################################################################################################################

echo >&2
echo >&2 "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
echo >&2 "::::                                                                                                                            ::::"
echo >&2 "                                      ParaMonte library version ${ParaMonteVersion} build on ${OSNAME}                              "
echo >&2 "::::                                                                                                                            ::::"
echo >&2 "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
echo >&2

echo >&2
echo >&2 "-- ${BUILD_NAME} -             working directory: $(pwd)"
echo >&2 "-- ${BUILD_NAME} -         source file directory: ${sourceFileDir}"
echo >&2 "-- ${BUILD_NAME} -     current system's platform: ${OSNAME}"
echo >&2 "-- ${BUILD_NAME} - current system's architecture: ${ARCHITECTURE}"
echo >&2

echo >&2
echo >&2 "-- ${BUILD_NAME} -            ParaMonte_ROOT_DIR: ${ParaMonte_ROOT_DIR}"
echo >&2 "-- ${BUILD_NAME} -         ParaMonteTest_SRC_DIR: ${ParaMonteTest_SRC_DIR}"
echo >&2 "-- ${BUILD_NAME} -       ParaMonteKernel_SRC_DIR: ${ParaMonteKernel_SRC_DIR}"
echo >&2 "-- ${BUILD_NAME} -      ParaMonte_ROOT_BUILD_DIR: ${ParaMonte_ROOT_BUILD_DIR}"
echo >&2

####################################################################################################################################
#### usage
####################################################################################################################################

usage()
{
cat << EndOfMessage

    usage:

        buildParaMonte.sh
        -L <language: C/C++/Fortran/MATLAB/Python>
        -s <compiler suite: intel/gnu>
        -b <build mode: release/testing/debug>
        -l <library type: static/dynamic>
        -c <coarray: none/single/shared/distributed>
        -m <mpi enabled: true/false>
        -i <C-Fortran interface enabled: true/false>
        -e <heap allocation enabled: true/false>
        -t <ParaMonte test run enabled: none/all/pm/nopm>
        -x <ParaMonte example run enabled: true/false>
        -f <path to Fortran compiler>
        -M <path to mpiexec>
        -F <purge the existing prerequisite library installations and perform a fresh installation>
        -y <assume yes as answer to all installation permission inquiries>
        -B <perform GCC bootstrap installation>
        -n <default number of processors for parallel application>
        -a <clean bash variables upon exit from the script>
        -h <help on the script usage>

    example:

        buildParaMonte.sh -b release -l dynamic -c none -m true -i true -d true -n 3

    flag definitions:

        -L | --lang             : the ParaMonte library interface programming language: C, C++, Fortran, MATLAB, Python
        -s | --compiler_suite   : the ParaMonte library build compiler suite: intel, gnu
        -b | --build            : the ParaMonte library build type: release, testing, debug
        -l | --lib              : the ParaMonte library type: static, dynamic
        -c | --caf              : the ParaMonte library Coarray Fortran parallelism type: none, single, shared, distributed
        -m | --mpi_enabled      : the ParaMonte library MPI parallelism enabled?: true, false
        -i | --cfi_enabled      : the ParaMonte library C-Fortran interface enabled? must be true if the library is to be called from non-Fortran languages: true, false
        -e | --heap_enabled     : the ParaMonte library heap array allocation enabled?: true, false
        -t | --test             : the ParaMonte library test run enabled?: none, all, pm, nopm
        -x | --exam_enabled     : the ParaMonte library examples run enabled?: true, false
        -S | --shared_enabled   : release the shared library dependencies in the binary release of the ParaMonte library: true, false
        -f | --fortran          : path to Fortran compiler. If provided, the ParaMonte library will be built via the specified compiler.
        -M | --mpiexec          : path to mpiexec routine. If provided, it will be used to find the MPI library.
        -F | --fresh            : enables a fresh installation of all of the prerequisites of ParaMonte library. Applicable only to GNU compiler suite.
        -y | --yes-to-all       : if a fresh installation of all of the prerequisites is needed, automatically answer yes to all permission requests.
        -B | --bootstrap        : enables robust bootstrap build when building the required GCC version with an old GCC version. Applicable only to GNU compiler suite.
        -n | --nproc            : the default number of processes (coarray images) on which the ParaMonte examples/tests (if any) will be run: positive integer
        -a | --clean            : clean the environmental variables upon exit, if flag is provided.
        -h | --help             : help with the script usage


EndOfMessage
}

####################################################################################################################################
# Configure ParaMonte Build
####################################################################################################################################

cmakeInstallEnabled=false
cafInstallEnabled=false
gnuInstallEnabled=false
mpiInstallEnabled=false

unset PMCS
unset BTYPE
unset LTYPE
unset TTYPE
unset CAFTYPE
unset MPI_ENABLED
unset CFI_ENABLED
unset GCC_BOOTSTRAP
unset MPIEXEC_PATH_USER
unset HEAP_ARRAY_ENABLED
unset INTERFACE_LANGUAGE
unset FOR_COARRAY_NUM_IMAGES
unset Fortran_COMPILER_PATH_USER
ParaMonteExample_RUN_ENABLED=true
localInstallationEnabled=false
freshInstallEnabled=false
YES_TO_ALL_DISABLED=true
CODECOV_ENABLED=false
DRYRUN_ENABLED=false
shared_enabled=true
CLEAN=false

chmod a+x ./configParaMonte.sh
echo >&2
echo >&2 "-- ${BUILD_NAME} - configuring ParaMonte Build..."
echo >&2 "-- ${BUILD_NAME} - default configuration: "
source ./configParaMonte.sh
#output=$("./configParaMonte.sh")
#echo "${output}"

while [ "$1" != "" ]; do
    case $1 in
        -L | --lang )           shift
                                INTERFACE_LANGUAGE=$1
                                ;;
        -s | --compiler_suite ) shift
                                PMCS=$1; export PMCS;
                                ;;
        -b | --build )          shift
                                BTYPE=$1; export BTYPE
                                ;;
        -l | --lib )            shift
                                LTYPE=$1; export LTYPE
                                ;;
        -c | --caf )            shift
                                CAFTYPE=$1; export CAFTYPE
                                ;;
        -m | --mpi_enabled )    shift
                                MPI_ENABLED=$1; export MPI_ENABLED
                                ;;
        -i | --cfi_enabled )    shift
                                CFI_ENABLED=$1; export CFI_ENABLED
                                ;;
        -e | --heap_enabled )   shift
                                HEAP_ARRAY_ENABLED=$1; export HEAP_ARRAY_ENABLED
                                ;;
        -t | --test )           shift
                                TTYPE=$1; export TTYPE
                                ;;
        -x | --exam_enabled )   shift
                                ParaMonteExample_RUN_ENABLED=$1; export ParaMonteExample_RUN_ENABLED
                                ;;
        -S | --shared_enabled ) shift
                                shared_enabled=$1; export shared_enabled
                                ;;
        -f | --fortran )        shift
                                Fortran_COMPILER_PATH_USER=$1; export Fortran_COMPILER_PATH_USER
                                ;;
        -M | --mpiexec )        shift
                                MPIEXEC_PATH_USER=$1; export MPIEXEC_PATH_USER
                                ;;
        -n | --nproc )          shift
                                FOR_COARRAY_NUM_IMAGES=$1; export FOR_COARRAY_NUM_IMAGES
                                ;;
        -F | --fresh )          freshInstallEnabled=true; export freshInstallEnabled
                                ;;
        -o | --local )          localInstallationEnabled=true; export localInstallationEnabled
                                ;;
        -d | --dryrun )         DRYRUN_ENABLED=true; export DRYRUN_ENABLED
                                ;;
        -c | --codecov )        CODECOV_ENABLED=true; export CODECOV_ENABLED
                                ;;
        -y | --yes-to-all )     YES_TO_ALL_DISABLED=false; export YES_TO_ALL_DISABLED
                                ;;
        -B | --bootstrap )      GCC_BOOTSTRAP="--bootstrap"; export GCC_BOOTSTRAP
                                ;;
        -a | --clean )          shift
                                CLEAN=true
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     echo >&2 "-- ${BUILD_NAME} - FATAL: the input flag is not recognized: $1"
                                usage
                                exit 1
    esac
    shift
done

echo >&2
echo >&2 "-- ${BUILD_NAME} - current requested configuration: "
. ./configParaMonte.sh

####################################################################################################################################
# check flag consistencies
####################################################################################################################################

if [ -z ${INTERFACE_LANGUAGE+x} ]; then
    echo >&2
    echo >&2 "-- ${BUILD_NAME} - FATAL: The environmental INTERFACE_LANGUAGE must be specified as input."
    echo >&2
    usage
    echo >&2
    echo >&2 "-- ${BUILD_NAME} - gracefully exiting."
    echo >&2
    exit 1
fi

CAF_ENABLED=false
if [ "${CAFTYPE}" = "single" ] || [ "${CAFTYPE}" = "shared" ] || [ "${CAFTYPE}" = "distributed" ]; then
    CAF_ENABLED=true
fi
export CAF_ENABLED

if [ "${isMacOS}" = "true" ]; then
    #if [[ "${PMCS}" =~ .*"intel".* ]]; then
    if [[ ${PMCS} == [iI][nN][tT][eE][lL] ]]; then
        if [[ "${MPI_ENABLED}" =~ .*"true".* ]]; then
            echo >&2
            echo >&2 "-- ${BUILD_NAME} - FATAL: incompatible input flags specified by the user:"
            echo >&2 "-- ${BUILD_NAME} - FATAL:     -s | --compiler_suite : ${PMCS}"
            echo >&2 "-- ${BUILD_NAME} - FATAL:     -m | --mpi_enabled : ${MPI_ENABLED}"
            echo >&2 "-- ${BUILD_NAME} - FATAL: \"--compiler_suite ${PMCS}\" cannot be used along with \"--mpi_enabled ${MPI_ENABLED}\" on macOS."
            echo >&2 "-- ${BUILD_NAME} - FATAL: For parallel ParaMonte builds on macOS, use \"--compiler_suite gnu\" instead."
            echo >&2
            echo >&2 "-- ${BUILD_NAME} - gracefully exiting."
            echo >&2
            exit 1
        fi
        if [[ "${CAF_ENABLED}" =~ .*"true".* ]]; then
            echo >&2
            echo >&2 "-- ${BUILD_NAME} - FATAL: incompatible input flags specified by the user:"
            echo >&2 "-- ${BUILD_NAME} - FATAL:     -s | --compiler_suite : ${PMCS}"
            echo >&2 "-- ${BUILD_NAME} - FATAL:     -c | --caf_enabled : ${CAF_ENABLED}"
            echo >&2 "-- ${BUILD_NAME} - FATAL: \"--compiler_suite ${PMCS}\" cannot be used along with \"--caf_enabled ${CAF_ENABLED}\" on macOS."
            echo >&2 "-- ${BUILD_NAME} - FATAL: For parallel ParaMonte builds on macOS, use \"--compiler_suite gnu\" instead."
            echo >&2
            echo >&2 "-- ${BUILD_NAME} - gracefully exiting."
            echo >&2
            exit 1
        fi
    else
        if [ "${MPI_ENABLED}" = "true" ] || [ "${CAF_ENABLED}" = "true" ]; then PMCS="gnu"; fi
    fi
fi

if [ "${CFI_ENABLED}" = "true" ]; then
    if [ "${CAFTYPE}" = "shared" ] || [ "${CAFTYPE}" = "distributed" ]; then
        echo >&2
        echo >&2 "-- ${BUILD_NAME} - FATAL: incompatible input flags specified by the user:"
        echo >&2 "-- ${BUILD_NAME} - FATAL:     -i | --cfi_enabled : ${CFI_ENABLED}"
        echo >&2 "-- ${BUILD_NAME} - FATAL:     -c | --caf : ${CAFTYPE}"
        echo >&2 "-- ${BUILD_NAME} - FATAL: coarray parallelism is not available in C language."
        echo >&2
        echo >&2 "-- ${BUILD_NAME} - gracefully exiting."
        echo >&2
        exit 1
    fi
fi

if [ "${MPI_ENABLED}" = "true" ]; then
    if [ "${CAFTYPE}" = "shared" ] || [ "${CAFTYPE}" = "distributed" ]; then
        echo >&2
        echo >&2 "-- ${BUILD_NAME} - FATAL: incompatible input flags specified by the user:"
        echo >&2 "-- ${BUILD_NAME} - FATAL:     -m | --mpi_enabled : ${MPI_ENABLED}"
        echo >&2 "-- ${BUILD_NAME} - FATAL:     -c | --caf : ${CAFTYPE}"
        echo >&2 "-- ${BUILD_NAME} - FATAL: coarray parallelism cannot be mixed with MPI in the current version of ParaMonte."
        echo >&2
        echo >&2 "-- ${BUILD_NAME} - gracefully exiting."
        echo >&2
        exit 1
    fi
fi

if [ "${LTYPE}" = "dynamic" ]; then
    if [ "${CAFTYPE}" = "shared" ] || [ "${CAFTYPE}" = "distributed" ]; then
        echo >&2
        echo >&2 "-- ${BUILD_NAME} - FATAL: incompatible input flags specified by the user:"
        echo >&2 "-- ${BUILD_NAME} - FATAL:     -l | --lib : ${LTYPE}"
        echo >&2 "-- ${BUILD_NAME} - FATAL:     -c | --caf : ${CAFTYPE}"
        echo >&2 "-- ${BUILD_NAME} - FATAL: ParaMonte dynamic library build with coarray parallelism currently unsupported."
        echo >&2
        echo >&2 "-- ${BUILD_NAME} - gracefully exiting."
        echo >&2
        exit 1
    fi
fi

####################################################################################################################################
#### set local dependencies paths
####################################################################################################################################

ParaMonte_REQ_DIR="${ParaMonte_ROOT_DIR}/build/prerequisites"; export ParaMonte_REQ_DIR
ParaMonte_REQ_INSTALL_DIR="${ParaMonte_REQ_DIR}/prerequisites/installations"; export ParaMonte_REQ_INSTALL_DIR

#############################################
#### set up the local CMAKE installation path
#############################################

if [ -d "${ParaMonte_REQ_INSTALL_DIR}" ]; then
    CMAKE_LOCAL_INSTALLATION_PATH="$(find "${ParaMonte_REQ_INSTALL_DIR}/" -path **/bin/cmake)"
    if [ -f "${CMAKE_LOCAL_INSTALLATION_PATH}" ]; then
        CMAKE_LOCAL_INSTALLATION_BIN_DIR="$(dirname "${CMAKE_LOCAL_INSTALLATION_PATH}")"
        export CMAKE_LOCAL_INSTALLATION_PATH
        export CMAKE_LOCAL_INSTALLATION_BIN_DIR
        CMAKE_LOCAL_INSTALLATION_VERSION="$(cmake --version)"
        CMAKE_LOCAL_INSTALLATION_VERSION_ARRAY=($CMAKE_LOCAL_INSTALLATION_VERSION)
        CMAKE_LOCAL_INSTALLATION_VERSION="${CMAKE_LOCAL_INSTALLATION_VERSION_ARRAY[2]}"
        unset CMAKE_LOCAL_INSTALLATION_VERSION_ARRAY
    else
        unset CMAKE_LOCAL_INSTALLATION_PATH
    fi
fi

############################################
#### set up the local GNU installation paths
############################################

if [ -d "${ParaMonte_REQ_INSTALL_DIR}/gnu" ]; then

    GNU_LOCAL_INSTALLATION_GFORTRAN_PATH="$(find "${ParaMonte_REQ_INSTALL_DIR}/gnu" -name gfortran)"

    if [ -f "${GNU_LOCAL_INSTALLATION_GFORTRAN_PATH}" ]; then

        GNU_LOCAL_INSTALLATION_BIN_DIR="$(dirname "${GNU_LOCAL_INSTALLATION_GFORTRAN_PATH}")"
        export GNU_LOCAL_INSTALLATION_GFORTRAN_PATH
        export GNU_LOCAL_INSTALLATION_BIN_DIR

        GNU_LOCAL_INSTALLATION_LIB_DIR="${GNU_LOCAL_INSTALLATION_BIN_DIR}/../lib"
        if [ -d "${GNU_LOCAL_INSTALLATION_LIB_DIR}" ]; then
            export GNU_LOCAL_INSTALLATION_LIB_DIR
        else
            unset GNU_LOCAL_INSTALLATION_LIB_DIR
        fi

        GNU_LOCAL_INSTALLATION_LIB64_DIR="${GNU_LOCAL_INSTALLATION_BIN_DIR}/../lib64"
        if [ -d "${GNU_LOCAL_INSTALLATION_LIB64_DIR}" ]; then
            export GNU_LOCAL_INSTALLATION_LIB64_DIR
        else
            unset GNU_LOCAL_INSTALLATION_LIB64_DIR
        fi

    else

        unset GNU_LOCAL_INSTALLATION_GFORTRAN_PATH

    fi

fi

############################################
#### set up the local MPI installation paths
############################################

if [ -d "${ParaMonte_REQ_INSTALL_DIR}" ]; then

    MPI_LOCAL_INSTALLATION_MPIEXEC_PATH="$(find "${ParaMonte_REQ_INSTALL_DIR}"/ -name mpiexec)"

    if [ -f "${MPI_LOCAL_INSTALLATION_MPIEXEC_PATH}" ]; then

        MPI_LOCAL_INSTALLATION_BIN_DIR="$(dirname "${MPI_LOCAL_INSTALLATION_MPIEXEC_PATH}")"
        export MPI_LOCAL_INSTALLATION_MPIEXEC_PATH
        export MPI_LOCAL_INSTALLATION_BIN_DIR

        MPI_LOCAL_INSTALLATION_LIB_DIR="${MPI_LOCAL_INSTALLATION_BIN_DIR}/../lib"
        if [ -d "${MPI_LOCAL_INSTALLATION_LIB_DIR}" ]; then
            export MPI_LOCAL_INSTALLATION_LIB_DIR
        else
            unset MPI_LOCAL_INSTALLATION_LIB_DIR
        fi

        MPI_LOCAL_INSTALLATION_LIB64_DIR="${MPI_LOCAL_INSTALLATION_BIN_DIR}/../lib64"
        if [ -d "${MPI_LOCAL_INSTALLATION_LIB64_DIR}" ]; then
            export MPI_LOCAL_INSTALLATION_LIB64_DIR
        else
            unset MPI_LOCAL_INSTALLATION_LIB64_DIR
        fi

    else

        unset MPI_LOCAL_INSTALLATION_MPIEXEC_PATH

    fi

fi

############################################
#### set up the local CAF installation paths
############################################

if [ -d "${ParaMonte_REQ_INSTALL_DIR}/opencoarrays" ]; then

    CAF_LOCAL_INSTALLATION_WRAPPER_PATH="$(find "${ParaMonte_REQ_INSTALL_DIR}"/opencoarrays/ -path **/bin/caf)"

    if [ -f "${CAF_LOCAL_INSTALLATION_WRAPPER_PATH}" ]; then

        CAF_LOCAL_INSTALLATION_BIN_DIR="$(dirname "${CAF_LOCAL_INSTALLATION_WRAPPER_PATH}")"
        export CAF_LOCAL_INSTALLATION_BIN_DIR

        CAF_LOCAL_INSTALLATION_LIB_DIR="${CAF_LOCAL_INSTALLATION_BIN_DIR}/../lib"
        if [ -d "${CAF_LOCAL_INSTALLATION_LIB_DIR}" ]; then
            export CAF_LOCAL_INSTALLATION_LIB_DIR
        else
            unset CAF_LOCAL_INSTALLATION_LIB_DIR
        fi

        CAF_LOCAL_INSTALLATION_LIB64_DIR="${CAF_LOCAL_INSTALLATION_BIN_DIR}/../lib64"
        if [ -d "${CAF_LOCAL_INSTALLATION_LIB64_DIR}" ]; then
            export CAF_LOCAL_INSTALLATION_LIB64_DIR
        else
            unset CAF_LOCAL_INSTALLATION_LIB64_DIR
        fi

        CAF_LOCAL_INSTALLATION_SETUP_FILE="${CAF_LOCAL_INSTALLATION_BIN_DIR}/../setup.sh";
        if [ -f "${CAF_LOCAL_INSTALLATION_SETUP_FILE}" ]; then
            export CAF_LOCAL_INSTALLATION_SETUP_FILE
        else
            unset CAF_LOCAL_INSTALLATION_SETUP_FILE
        fi

    else

        unset CAF_LOCAL_INSTALLATION_WRAPPER_PATH

    fi

fi

####################################################################################################################################
#### check cmake version
####################################################################################################################################

cmakeVersion="$(cmake --version)"
cmakeVersionArray=($cmakeVersion)
cmakeVersion="${cmakeVersionArray[2]}"

if [ "${cmakeVersion}" = "" ]; then
    cmakeInstallEnabled=true
else
    cmakeInstallEnabled=false
    compareVersions "${cmakeVersion}" "${cmakeVersionParaMonteCompatible}"
    if [ "$?" = "2" ]; then
        cmakeInstallEnabled=true
        echo >&2 "-- ${BUILD_NAME} - failed to detect a ParaMonte-compatible installation of cmake!"
    else
        echo >&2 "-- ${BUILD_NAME} - the current cmake installation is ParaMonte compatible!"
    fi
fi
export cmakeInstallEnabled

echo >&2
echo >&2 "-- ${BUILD_NAME} -              cmake path: $(which cmake)"
echo >&2 "-- ${BUILD_NAME} -   cmake version current: ${cmakeVersion}"
echo >&2 "-- ${BUILD_NAME} -  cmake version required: ${cmakeVersionParaMonteCompatible}"
echo >&2 "-- ${BUILD_NAME} -     cmakeInstallEnabled: ${cmakeInstallEnabled}"

####################################################################################################################################
#### set compiler suite
####################################################################################################################################

if [ -z ${PMCS+x} ]; then

    SUITE_LIST="intel gnu"

else

    if [[ ${PMCS} == [iI][nN][tT][eE][lL] ]]; then
        PMCS=intel
        SUITE_LIST=${PMCS}
    else
        if [[ ${PMCS} == [gG][nN][uU] ]]; then
            PMCS=gnu
            SUITE_LIST=${PMCS}
        else
            echo >&2
            echo >&2 "-- ${BUILD_NAME} - FATAL: the requested compiler suite ${PMCS} is unrecognized."
            echo >&2 "-- ${BUILD_NAME} - FATAL: please choose either intel or gnu, or drop the option."
            echo >&2 "-- ${BUILD_NAME} - FATAL: The installer will automatically find the proper compiler suite."
            echo >&2
            usage
            echo >&2
            echo >&2 "-- ${BUILD_NAME} - gracefully exiting."
            echo >&2
            exit 1
        fi
    fi

fi

####################################################################################################################################
#### detect the available compiler suites, C/Fortran compilers and CAF/MPI libraries
####################################################################################################################################

gnuCCompilerName=gcc
gnuFortranCompilerName=gfortran

intelCCompilerName=icc
intelFortranCompilerName=ifort

LANG_LIST="C Fortran"

for SUITE in $SUITE_LIST
do

    echo >&2
    echo >&2 "-- ${BUILD_NAME}Compiler - checking for the ${SUITE} compilers and libraries presence..."
    echo >&2

    for LANG in $LANG_LIST
    do

        suiteLangCompilerName="${SUITE}${LANG}CompilerName"
        suiteLangCompilerPath="${SUITE}${LANG}CompilerPath"
        suiteLangCompilerVersion="${SUITE}${LANG}CompilerVersion"

        if eval "command -v ${!suiteLangCompilerName} >/dev/null 2>&1"; then

            eval "unset ${suiteLangCompilerPath}"
            eval ${suiteLangCompilerPath}='$(command -v ${!suiteLangCompilerName})'

            echo >&2 "-- ${BUILD_NAME}Compiler - the ${SUITE} ${!suiteLangCompilerName} detected at ${suiteLangCompilerPath} = ${!suiteLangCompilerPath}"

            #### get compiler version

            if [ "${LANG}" = "C" ]; then

                eval ${suiteLangCompilerVersion}="$(${!suiteLangCompilerName} -dumpversion)"
                echo >&2 "-- ${BUILD_NAME}Compiler - the ${SUITE} ${LANG} compiler version is ${suiteLangCompilerVersion} = ${!suiteLangCompilerVersion}"

            elif [ "${LANG}" = "Fortran" ]; then

                tempDir=$(mktemp -d "${TMPDIR:-/tmp}/cversion.XXXXXXXXX")
                echo >&2 "-- ${BUILD_NAME}Compiler - changing directory to: ${tempDir}"
                cd "${tempDir}" && cp "${ParaMonte_ROOT_DIR}/auxil/getCompilerVersion.f90" "./getCompilerVersion.f90"

                if ${!suiteLangCompilerPath} getCompilerVersion.f90 -o getCompilerVersion.exe; then

                    chmod +x getCompilerVersion.exe && ./getCompilerVersion.exe && {
                        eval ${suiteLangCompilerVersion}='$(head -n 1 getCompilerVersion.tmp)'
                        echo >&2 "-- ${BUILD_NAME}Compiler - the ${SUITE} ${LANG} compiler version is ${suiteLangCompilerVersion} = ${!suiteLangCompilerVersion}"
                        isParaMonteCompatibleCompiler=$(head -n 1 isParaMonteCompatibleCompiler.tmp)
                        isGfortran10=$(head -n 1 isGfortran10.tmp); export isGfortran10
                        if [ "$isParaMonteCompatibleCompiler" = "true" ]; then
                            echo >&2 "-- ${BUILD_NAME}Compiler - the ${SUITE} ${LANG} compiler is ParaMonte compatible!"
                            eval "export ${suiteLangCompilerPath}"
                        else
                            echo >&2 "-- ${BUILD_NAME}Compiler - the ${SUITE} ${LANG} compiler is not ParaMonte compatible. skipping..."
                            unset ${suiteLangCompilerPath}
                        fi
                        rm *.tmp *.exe
                    } || {
                        echo >&2 "-- ${BUILD_NAME}Compiler - failed to detect the ${SUITE} ${LANG} compiler version. skipping..."
                        unset ${suiteLangCompilerPath}
                    }

                else

                    echo >&2 "-- ${BUILD_NAME}Compiler - failed to detect the ${SUITE} ${LANG} compiler version. skipping..."
                    unset ${suiteLangCompilerPath}

                fi

                echo >&2 "-- ${BUILD_NAME}Compiler - changing directory to: ${ParaMonte_ROOT_DIR}"
                cd "${ParaMonte_ROOT_DIR}"
                rm -rf "${tempDir}"

            fi

        else

            echo >&2 "-- ${BUILD_NAME}Compiler - the ${SUITE} ${!suiteLangCompilerName} not found. skipping..."
            unset ${suiteLangCompilerPath}

        fi

    done

done

####################################################################################################################################
#### if no compiler has been identified and version extraction has failed,
#### try one last time by searching for the alternative compiler names.
####################################################################################################################################

if [[ "${SUITE_LIST}" =~ .*"gnu".* ]] && [ -z ${gnuFortranCompilerPath+x} ]; then

    echo >&2
    echo >&2 "-- ${BUILD_NAME}Compiler - searching all system paths for a the gnu Fortran compiler..."
    echo >&2

    #### return a comma-separated list of available compiler paths

    clist=$(( IFS=:; for p in $PATH; do unset lsout; lsout=$(ls -dm "$p"/${gnuFortranCompilerName}*); if ! [[ -z "${lsout// }" ]]; then echo "${lsout}, "; fi; done ) 2>/dev/null)

    for cname in $(echo $clist | sed "s/,/ /g"); do

        echo >&2 "-- ${BUILD_NAME}Compiler - checking ${cname}"

        gnuFortranCompilerVersion="$("${cname}" -dumpversion)" && \
        gnuFortranCompilerPath="$(command -v "${cname}")" && \
        gnuFortranCompilerName="$(basename "${cname}")" && \
        {

            echo >&2 "-- ${BUILD_NAME}Compiler - the gnu Fortran compiler version: ${gnuFortranCompilerVersion}"
            echo >&2 "-- ${BUILD_NAME}Compiler - the gnu Fortran compiler path: ${gnuFortranCompilerPath}"
            echo >&2 "-- ${BUILD_NAME}Compiler - the gnu Fortran compiler name: ${gnuFortranCompilerName}"

            compareVersions "${gnuFortranCompilerVersion}" "${gnuVersionParaMonteCompatible}"
            if [ "$?" = "2" ]; then

                echo >&2 "-- ${BUILD_NAME}Compiler - the gnu Fortran compiler is not ParaMonte compatible. skipping..."
                unset gnuFortranCompilerVersion
                unset gnuFortranCompilerPath

            else

                echo >&2 "-- ${BUILD_NAME}Compiler - the gnu Fortran compiler is ParaMonte compatible!"

                export gnuFortranCompilerVersion
                export gnuFortranCompilerPath
                export gnuFortranCompilerName

                #### check for a companion C processor
                #gnuFortranCompilerDir="$(dirname "${gnuFortranCompilerPath}")"

                gnuCCompanionDetected=false
                tempFilePath="${gnuFortranCompilerPath//gfortran/gcc}"
                if [ -f "${tempFilePath}" ]; then
                    tempFileVersion="$("${tempFilePath}" -dumpversion)" && \
                    {
                        gnuCCompanionDetected=true
                        gnuCCompilerVersion="${tempFileVersion}"
                        gnuCCompilerPath="${tempFilePath}"
                        gnuCCompilerName="$(basename "${tempFilePath}")"
                    }
                fi
                if [ "${gnuCCompanionDetected}" = "true" ]; then
                    echo >&2 "-- ${BUILD_NAME}Compiler - a companion gnu C processor detected!"
                    echo >&2 "-- ${BUILD_NAME}Compiler - the gnu C compiler version: ${gnuCCompilerVersion}"
                    echo >&2 "-- ${BUILD_NAME}Compiler - the gnu C compiler path: ${gnuCCompilerPath}"
                    echo >&2 "-- ${BUILD_NAME}Compiler - the gnu C compiler name: ${gnuCCompilerName}"
                    export gnuCCompilerVersion
                    export gnuCCompilerPath
                    export gnuCCompilerName
                else
                    echo >&2 "-- ${BUILD_NAME}Compiler - no companion gnu C processor detected. skipping..."
                fi

                break

            fi

        } || {

            echo >&2 "-- ${BUILD_NAME}Compiler - the gnu Fortran compiler version extraction failed. skipping..."
            unset ${suiteLangCompilerVersion}
            unset ${suiteLangCompilerPath}
            continue

        }

    done
fi

if ! [ -z ${intelFortranCompilerPath+x} ]; then
    intelFortranCompilerBinDir="$(dirname "${intelFortranCompilerPath}")"
    if [[ ":$PATH:" != *":${intelFortranCompilerBinDir}:"* ]]; then
        PATH="${intelFortranCompilerBinDir}:${PATH}"
        export PATH
    fi
fi

if ! [ -z ${gnuFortranCompilerPath+x} ]; then
    gnuFortranCompilerBinDir="$(dirname "${gnuFortranCompilerPath}")"
    if [[ ":$PATH:" != *":${gnuFortranCompilerBinDir}:"* ]]; then
        PATH="${gnuFortranCompilerBinDir}:${PATH}"
        export PATH
    fi
    if ! [ -z ${gnuFortranCompilerName+x} ]; then
        gnuFortranCompilerName="$(basename "${gnuFortranCompilerPath}")"
    fi
    if ! [ "${gnuFortranCompilerName}" = "gfortran" ]; then
        echo >&2 "gnuFortranCompilerName = ${gnuFortranCompilerName}"
        alias gfortran="${gnuFortranCompilerPath}" >> $HOME/.bashrc
        #shopt -s expand_aliases
        source $HOME/.bashrc
    fi
fi

####################################################################################################################################
#### identify the MPI wrappers
####################################################################################################################################

if [ "${MPI_ENABLED}" = "true" ] || [ "${CAF_ENABLED}" = "true" ]; then

    gnuCMpiWrapperList="mpicc"
    gnuFortranMpiWrapperList="mpifort:mpif90"

    intelCMpiWrapperList="mpiicc"
    intelFortranMpiWrapperList="mpiifort"

    for SUITE in $SUITE_LIST
    do

        echo >&2
        echo >&2 "-- ${BUILD_NAME}MPI - checking for ${SUITE} MPI wrappers and libraries presence..."
        echo >&2

        for LANG in $LANG_LIST
        do

            suiteLangMpiWrapperName="${SUITE}${LANG}MpiWrapperName"
            suiteLangMpiWrapperList="${SUITE}${LANG}MpiWrapperList"
            suiteLangMpiWrapperPath="${SUITE}${LANG}MpiWrapperPath"

            searchFailed=false
            for mpiWrapperName in ${!suiteLangMpiWrapperList//:/ }
            do

                eval ${suiteLangMpiWrapperName}="${mpiWrapperName}"

                if eval "command -v ${!suiteLangMpiWrapperName} >/dev/null 2>&1"; then

                    eval "unset ${suiteLangMpiWrapperPath}"
                    eval ${suiteLangMpiWrapperPath}='$(command -v ${!suiteLangMpiWrapperName})'

                    # ensure the MPI library is not intel

                    if [[ "${SUITE}" = "gnu" && "${!suiteLangMpiWrapperPath}" =~ .*"intel".* ]]; then

                        searchFailed=true

                    else

                        echo >&2 "-- ${BUILD_NAME}MPI - ${SUITE} ${!suiteLangMpiWrapperName} MPI wrapper detected at: ${suiteLangMpiWrapperPath} = ${!suiteLangMpiWrapperPath}"

                        if [ "${LANG}" = "Fortran" ]; then

                            # check if the compiler wrapper can compile a simple Fortran MPI test code.

                            tempDir=$(mktemp -d "${TMPDIR:-/tmp}/cversion.XXXXXXXXX")
                            echo >&2 "-- ${BUILD_NAME}MPI - changing directory to: ${tempDir}"
                            cd "${tempDir}" && cp "${ParaMonte_ROOT_DIR}/auxil/testMPI.f90" "./testMPI.f90" && \
                            {
                                ${!suiteLangMpiWrapperName} testMPI.f90 -o main.exe && mpiexec -n 1 ./main.exe
                            } && {
                            #} &> /dev/null && {
                                echo >&2 "-- ${BUILD_NAME}MPI - checking whether ${SUITE} ${!suiteLangMpiWrapperName} MPI wrapper compiles and runs a test program...yes"
                            }|| {
                                echo >&2 "-- ${BUILD_NAME}MPI - failed to compile a simple MPI test program with ${SUITE} ${!suiteLangMpiWrapperName}. skipping..."
                                unset ${suiteLangMpiWrapperPath}
                                if [ "${MPI_ENABLED}" = "true" ] ; then
                                    mpiInstallEnabled=true
                                fi
                            }
                            rm -rf "${tempDir}"
                            cd "${ParaMonte_ROOT_DIR}"

                        fi

                        break

                    fi

                else

                    searchFailed=true

                fi

                if [ "${searchFailed}" = "true" ]; then
                    echo >&2 "-- ${BUILD_NAME}MPI - failed to detect the ${SUITE} ${LANG} MPI wrapper. skipping..."
                    unset ${suiteLangMpiWrapperPath}
                fi

            done

        done

    done

fi

####################################################################################################################################
#### detect CAF wrapper
####################################################################################################################################

echo >&2
if [ "${CAF_ENABLED}" = "true" ]; then

    #unset CAF_LOCAL_INSTALLATION_WRAPPER_PATH

    if [ -z "${CAF_LOCAL_INSTALLATION_WRAPPER_PATH+x}" ] && ( [ -z ${PMCS+x} ] || [[ "${PMCS}" =~ .*"intel".* ]] ); then

        if [ -z "${intelFortranMpiWrapperPath+x}" ]; then
            echo >&2
            echo >&2 "-- ${BUILD_NAME}CAF - ${warning}: Failed to identify the Intel MPI library."
            echo >&2 "-- ${BUILD_NAME}CAF - ${warning}: The Intel MPI library is required to compile"
            echo >&2 "-- ${BUILD_NAME}CAF - ${warning}: Coarray Fortran applications via Intel compilers."
            echo >&2 "-- ${BUILD_NAME}CAF - ${warning}: The ParaMonte build will continue at the risk of failing..."
            echo >&2
        elif [ -z "${intelFortranCompilerPath+x}" ]; then
            echo >&2
            echo >&2 "-- ${BUILD_NAME}CAF - ${warning}: Failed to identify the Intel Fortran compiler path."
            echo >&2 "-- ${BUILD_NAME}CAF - ${warning}: The Intel MPI library and Fortran compiler are required"
            echo >&2 "-- ${BUILD_NAME}CAF - ${warning}: to compile Coarray Fortran applications via Intel compilers."
            echo >&2 "-- ${BUILD_NAME}CAF - ${warning}: The ParaMonte build will continue at the risk of failing..."
            echo >&2
        else
            CAF_LOCAL_INSTALLATION_WRAPPER_PATH="${intelFortranCompilerPath}"
            echo >&2 "-- ${BUILD_NAME}CAF - The inferred Coarray Fortran compiler wrapper path: ${CAF_LOCAL_INSTALLATION_WRAPPER_PATH}"
        fi

    fi

    if [ -z "${CAF_LOCAL_INSTALLATION_WRAPPER_PATH+x}" ] && ( [ -z ${PMCS+x} ] || [[ "${PMCS}" =~ .*"gnu".* ]] ); then

        # assume OpenCoarrays

        PMCS=gnu

        if command -v caf >/dev/null 2>&1; then
            CAF_LOCAL_INSTALLATION_WRAPPER_PATH="$(command -v caf)"
        fi

        if [ -f "${CAF_LOCAL_INSTALLATION_WRAPPER_PATH}" ]; then
            echo >&2 "-- ${BUILD_NAME}CAF - OpenCoarrays Fortran compiler wrapper detected at: ${CAF_LOCAL_INSTALLATION_WRAPPER_PATH}"
            cafVersion="$(caf -dumpversion)"
            cafVersionRequired="${gnuVersionParaMonteCompatible}"
            echo >&2 "-- ${BUILD_NAME}CAF - caf version: ${cafVersion}"
            echo >&2 "-- ${BUILD_NAME}CAF - caf version required: ${cafVersionRequired}"
            compareVersions "$cafVersion" "$cafVersionRequired"
            if [ "$?" = "2" ]; then
                echo >&2
                echo >&2 "-- ${BUILD_NAME}CAF - ${warning}: The existing OpenCoarrays caf compiler wrapper version is not ParaMonte compatible."
                echo >&2 "-- ${BUILD_NAME}CAF - ${warning}: A fresh installation of the OpenCoarrays library might be needed."
                echo >&2
                cafInstallEnabled=true
                #mpiInstallEnabled=true
                #gnuInstallEnabled=true
            else
                echo >&2 "-- ${BUILD_NAME}CAF - OpenCoarrays Fortran compiler wrapper is ParaMonte compatible!"
            fi
        else
            echo >&2
            echo >&2 "-- ${BUILD_NAME}CAF - ${warning}: The OpenCoarrays caf compiler wrapper was not found on your system or at least,"
            echo >&2 "-- ${BUILD_NAME}CAF - ${warning}: the path containing the OpenCoarrays caf compiler wrapper is not in the PATH environmental variable."
            echo >&2 "-- ${BUILD_NAME}CAF - ${warning}: A fresh installation of the OpenCoarrays library might be needed."
            echo >&2
            cafInstallEnabled=true
            #mpiInstallEnabled=true
            #gnuInstallEnabled=true
        fi

    fi

fi

echo >&2

####################################################################################################################################
#### set the ParaMonte compiler suite
####################################################################################################################################

if [ -z ${Fortran_COMPILER_PATH_USER+x} ]; then

    #### attempt to set up build with the Intel compiler suite

    if [ -z ${PMCS+x} ]; then

        # if no preference, then priority is with intel

        if ! ${intelFortranCompilerPath+false}; then
            if [ "${MPI_ENABLED}" = "true" ] || [ "${CAF_ENABLED}" = "true" ]; then
                if ! ${intelFortranMpiWrapperPath+false}; then
                    PMCS=intel
                    COMPILER_VERSION=${intelFortranCompilerVersion}
                    Fortran_COMPILER_PATH="${intelFortranCompilerPath}"
                    MPIEXEC_PATH="$(dirname "${intelFortranMpiWrapperPath}")/mpiexec"
                    #MPIEXEC_PATH=$(command -v mpiexec)
                else
                    PMCS=gnu
                fi
            else
                PMCS=intel
                COMPILER_VERSION=${intelFortranCompilerVersion}
                Fortran_COMPILER_PATH="${intelFortranCompilerPath}"
            fi
        else
            PMCS=gnu
        fi

    else

        errorOccurred=false
        if [ "${PMCS}" = "intel" ]; then
            if ! ${intelFortranCompilerPath+false}; then
                if [ "${MPI_ENABLED}" = "true" ] || [ "${CAF_ENABLED}" = "true" ]; then
                    if ! ${intelFortranMpiWrapperPath+false}; then
                        COMPILER_VERSION=${intelFortranCompilerVersion}
                        Fortran_COMPILER_PATH="${intelFortranCompilerPath}"
                        MPIEXEC_PATH="$(dirname "${intelFortranMpiWrapperPath}")/mpiexec"
                    else
                        if [ -z ${MPIEXEC_PATH_USER+x} ]; then
                            echo >&2
                            echo >&2 "-- ${BUILD_NAME} - ${warning}: Failed to identify the Intel MPI library."
                            echo >&2 "-- ${BUILD_NAME} - ${warning}: The Intel MPI library is required to compile"
                            echo >&2 "-- ${BUILD_NAME} - ${warning}: the parallel ParaMonte library via Intel compilers."
                            echo >&2 "-- ${BUILD_NAME} - ${warning}: The library build will continue at the risk of failing..."
                            echo >&2
                        fi
                    fi
                else
                    COMPILER_VERSION=${intelFortranCompilerVersion}
                    Fortran_COMPILER_PATH="${intelFortranCompilerPath}"
                fi
            else
                errorOccurred=true
            fi
        fi

        if [ "${errorOccurred}" = "true" ]; then
            echo >&2
            echo >&2 "-- ${BUILD_NAME} - FATAL: the Fortran compiler and wrapper components of the"
            echo >&2 "-- ${BUILD_NAME} - FATAL: requested compiler suite (${PMCS}) could not be detected."
            echo >&2 "-- ${BUILD_NAME} - FATAL: Please make sure the ${PMCS} compiler suite is installed on your system"
            echo >&2 "-- ${BUILD_NAME} - FATAL: and is available in the environmental PATH variable of your shell."
            echo >&2 "-- ${BUILD_NAME} - FATAL: Otherwise, drop the \"-s intel\" or \"--compiler_suite intel\""
            echo >&2 "-- ${BUILD_NAME} - FATAL: from the input list of flags to the build script to the let"
            echo >&2 "-- ${BUILD_NAME} - FATAL: the script build the library with other available compilers."
            echo >&2
            echo >&2 "-- ${BUILD_NAME} - gracefully exiting."
            echo >&2
            exit 1
        fi

    fi

    #### If needed, set up build with the GNU compiler suite

    if [ "${PMCS}" = "gnu" ]; then
        if ! ${gnuFortranCompilerPath+false}; then
            COMPILER_VERSION=${gnuFortranCompilerVersion}
            Fortran_COMPILER_PATH="${gnuFortranCompilerPath}"
            gnuInstallEnabled=false
        else
            gnuInstallEnabled=true
            echo >&2
            echo >&2 "-- ${BUILD_NAME} - ${warning}: The GNU Fortran compiler could not be found on your system."
            echo >&2 "-- ${BUILD_NAME} - ${warning}: If you do not have GNU compiler suite installed on your system,"
            echo >&2 "-- ${BUILD_NAME} - ${warning}: ParaMonte may be able to install the GNU compiler suite for you."
            echo >&2
        fi

        if [ "${MPI_ENABLED}" = "true" ]; then
            if [ -f "${gnuFortranMpiWrapperPath}" ]; then
                MPIEXEC_PATH="$(dirname "${gnuFortranMpiWrapperPath}")/mpiexec"
                if [ -f "${MPIEXEC_PATH}" ]; then
                    mpiInstallEnabled=false
                else
                    if command -v mpiexec >/dev/null 2>&1; then
                        MPIEXEC_PATH="$(command -v mpiexec)"
                        if [[ "${MPIEXEC_PATH}" =~ .*"intel".* ]]; then
                            mpiInstallEnabled=true
                        else
                            mpiInstallEnabled=false
                        fi
                    else
                        mpiInstallEnabled=true
                    fi
                fi
            fi
        else
            mpiInstallEnabled=false
        fi

        if [ "${mpiInstallEnabled}" = "true" ]; then
            unset MPIEXEC_PATH
            # enable usage of the local installation of the GNU compiler for MPI-GNUcompiler compatibility
            gnuInstallEnabled=true
            echo >&2
            echo >&2 "-- ${BUILD_NAME} - ${warning}: The mpiexec executable could not be found on your system."
            echo >&2 "-- ${BUILD_NAME} - ${warning}: If you do not have an MPI library installed on your system,"
            echo >&2 "-- ${BUILD_NAME} - ${warning}: ParaMonte may be able to install the MPICH MPI library for you."
            echo >&2
        fi

        if [ "${CAF_ENABLED}" = "true" ]; then
            if [ -f "${CAF_LOCAL_INSTALLATION_WRAPPER_PATH}" ]; then
                COMPILER_VERSION="unknown"
                Fortran_COMPILER_PATH="${CAF_LOCAL_INSTALLATION_WRAPPER_PATH}"
                cafInstallEnabled=false
            else
                cafInstallEnabled=true
                mpiInstallEnabled=true
                gnuInstallEnabled=true
            fi
        fi

    fi

else # user has specified the Fortran compiler path

    cafInstallEnabled=false
    #gnuInstallEnabled=false

    PMCS="unknown"
    COMPILER_VERSION="unknown"
    Fortran_COMPILER_NAME=${Fortran_COMPILER_PATH_USER##*/}
    Fortran_COMPILER_PATH="${Fortran_COMPILER_PATH_USER}"
    echo >&2
    echo >&2 "-- ${BUILD_NAME}Compiler - user-requested compiler path: ${Fortran_COMPILER_PATH}"
    echo >&2 "-- ${BUILD_NAME}Compiler - user-requested compiler name: ${Fortran_COMPILER_NAME}"
    echo >&2
    if [[ "${Fortran_COMPILER_NAME}" =~ .*"gfortran".* || "${Fortran_COMPILER_NAME}" =~ .*"caf".* ]]; then
        PMCS=gnu
    fi
    if [ "${Fortran_COMPILER_NAME}" = "ifort" ]; then
        PMCS=intel
    fi

    if [ "${MPI_ENABLED}" = "true" ] && [ -z ${MPIEXEC_PATH_USER+x} ]; then
        if [ -f "${MPIEXEC_PATH}" ]; then
            mpiInstallEnabled=false
        else
            if command -v mpiexec >/dev/null 2>&1; then
                MPIEXEC_PATH="$(command -v mpiexec)"
                mpiInstallEnabled=false
            else
                mpiInstallEnabled=true
            fi
        fi
        #mpiInstallEnabled=true
        #gnuInstallEnabled=true
    fi

fi


if ! [ -z ${MPIEXEC_PATH_USER+x} ]; then
    mpiInstallEnabled=false
    if [ -f "${MPIEXEC_PATH_USER}" ]; then
        if ! [ -z ${MPIEXEC_PATH+x} ]; then
            echo >&2
            echo >&2 "-- ${BUILD_NAME} - ${warning}: The specified mpiexec path via the input flag -M | --mpiexec"
            echo >&2 "-- ${BUILD_NAME} - ${warning}: will overwrite the inferred mpiexec path by the build script."
            echo >&2 "-- ${BUILD_NAME} - ${warning}: user-specified mpiexec path: ${MPIEXEC_PATH_USER}"
            echo >&2 "-- ${BUILD_NAME} - ${warning}:       inferred mpiexec path: ${MPIEXEC_PATH}"
            echo >&2
            echo >&2 "-- ${BUILD_NAME} - gracefully exiting."
            echo >&2
        fi
        MPIEXEC_PATH="${MPIEXEC_PATH_USER}"
    else
        usage
        echo >&2
        echo >&2 "-- ${BUILD_NAME} - FATAL: The specified mpiexec path via the input flag -M | --mpiexec"
        echo >&2 "-- ${BUILD_NAME} - FATAL: does not point to a valid file name on the system."
        echo >&2 "-- ${BUILD_NAME} - FATAL: Please provide a valid path to the mpiexec executable"
        echo >&2 "-- ${BUILD_NAME} - FATAL: or drop it from the list of input flags."
        echo >&2
        echo >&2 "-- ${BUILD_NAME} - gracefully exiting."
        echo >&2
        exit 1
    fi
fi

#if [ "${MPI_ENABLED}" = "true" ]; then
#    if ! ${intelFortranMpiWrapperPath+false} && ! ${intelFortranCompilerPath+false}; then
#        PMCS=intel
#        COMPILER_VERSION=${intelFortranCompilerVersion}
#    else
#        #if ! ${gnuFortranMpiWrapperPath+false} && ! ${gnuFortranCompilerPath+false}; then
#            PMCS=gnu
#            COMPILER_VERSION=${gnuFortranCompilerVersion}
#        #fi
#    fi
#else
#    if ! ${intelFortranCompilerPath+false}; then
#        PMCS=intel
#        COMPILER_VERSION=${intelFortranCompilerVersion}
#    else
#        if ! ${gnuFortranCompilerPath+false}; then
#            PMCS=gnu
#            COMPILER_VERSION=${gnuFortranCompilerVersion}
#        fi
#    fi
#fi

#if [ "${MPI_ENABLED}" = "true" ]; then
#    echo >&2
#    echo >&2 "-- ${BUILD_NAME}MPI - intel Fortran MPI compiler wrapper: ${intelFortranMpiWrapperPath}"
#    echo >&2 "-- ${BUILD_NAME}MPI - intel Fortran compiler: ${intelFortranCompilerPath}"
#    echo >&2
#    echo >&2 "-- ${BUILD_NAME}MPI - gnu Fortran MPI compiler wrapper: ${gnuFortranMpiWrapperPath}"
#    echo >&2 "-- ${BUILD_NAME}MPI - gnu Fortran compiler: ${gnuFortranCompilerPath}"
#    echo >&2
#    if [ -z ${intelFortranMpiWrapperPath+x} ] || [ -z ${intelFortranCompilerPath+x} ]; then
#        #if [ -z ${gnuFortranMpiWrapperPath+x} ] || [ -z ${gnuFortranCompilerPath+x} ]; then
#            mpiInstallEnabled=true
#            gnuInstallEnabled=true
#        #fi
#    fi
#fi

####################################################################################################################################
#### set up ParaMonte library prerequisites
####################################################################################################################################

if [ "${freshInstallEnabled}" = "true" ] || [ "${localInstallationEnabled}" = "true" ]; then
    #cmakeInstallEnabled=true
    cafInstallEnabled=true
    mpiInstallEnabled=true
    gnuInstallEnabled=true
fi

echo >&2 "-- ${BUILD_NAME} -  localInstallationEnabled: ${localInstallationEnabled}"
echo >&2 "-- ${BUILD_NAME} -       freshInstallEnabled: ${freshInstallEnabled}"
echo >&2 "-- ${BUILD_NAME} -       cmakeInstallEnabled: ${cmakeInstallEnabled}"
echo >&2 "-- ${BUILD_NAME} -         cafInstallEnabled: ${cafInstallEnabled}"
echo >&2 "-- ${BUILD_NAME} -         mpiInstallEnabled: ${mpiInstallEnabled}"
echo >&2 "-- ${BUILD_NAME} -         gnuInstallEnabled: ${gnuInstallEnabled}"

# chmod 777 -R "${ParaMonte_ROOT_DIR}/auxil/prerequisites"

echo >&2
if [ "${cafInstallEnabled}" = "true" ] || [ "${mpiInstallEnabled}" = "true" ] || [ "${gnuInstallEnabled}" = "true" ] || [ "${cmakeInstallEnabled}" = "true" ]; then

    #ParaMonte_REQ_DIR="${ParaMonte_ROOT_DIR}/build/prerequisites"
    #CAF_LOCAL_INSTALLATION_SETUP_FILE="${ParaMonte_REQ_DIR}/prerequisites/installations/opencoarrays/${openCoarraysVersion}/setup.sh"

    if [ "${freshInstallEnabled}" = "true" ]; then
        rm -rf "${ParaMonte_REQ_DIR}"
    fi

    answer=y
    if [ ! -d "${ParaMonte_REQ_DIR}" ] || [ ! "$(ls -A ${ParaMonte_REQ_DIR})" ]; then

        echo >&2
        echo >&2 "-- ${BUILD_NAME} - ${warning}: ParaMonte build with the requested configuration requires the installations"
        echo >&2 "-- ${BUILD_NAME} - ${warning}: of either OpenCoarrays, MPICH MPI library (on Linux) or Open-MPI MPI library (on macOS),"
        echo >&2 "-- ${BUILD_NAME} - ${warning}: GNU compilers, or CMAKE on your system."
        echo >&2 "-- ${BUILD_NAME} - ${warning}: ParaMonte can install all  prerequisites on your system from the web, if needed."
        echo >&2 "-- ${BUILD_NAME} - ${warning}: The prerequisite build objects may occupy up to 5Gb of your system's memory."
        echo >&2

        ############################################################################################################################
        #### request permission to download and unpack prereqs
        ############################################################################################################################

        if [ "${YES_TO_ALL_DISABLED}" = "true" ]; then
            answerNotGiven=true
            while [ "${answerNotGiven}" = "true" ]; do
                read -p "-- ${BUILD_NAME} - Do you wish to continue with the installation of the prerequisites (y/n)? " answer
                if [[ $answer == [yY] || $answer == [yY][eE][sS] ]]; then
                    answer=y
                    answerNotGiven=false
                fi
                if [[ $answer == [nN] || $answer == [nN][oO] ]]; then
                    answer=n
                    answerNotGiven=false
                fi
                if [ "${answerNotGiven}" = "true" ]; then
                    echo >&2 "-- ${BUILD_NAME} - please enter either y or n"
                fi
            done
        else
            echo >&2 "-- ${BUILD_NAME} - Do you wish to continue with the installation of the prerequisites (y/n)? y"
            answer=y
        fi
        echo >&2

        ############################################################################################################################
        #### download and unpack prereqs
        ############################################################################################################################

        if ! [ "${isMacOS}" = "true" ]; then
            if [ "${answer}" = "y" ]; then

                echo >&2 "-- ${BUILD_NAME} - generating directory: ${ParaMonte_REQ_DIR}/"
                if ! [ -d "${ParaMonte_REQ_DIR}" ]; then
                    mkdir -p "${ParaMonte_REQ_DIR}"
                fi

                tarFileName="OpenCoarrays-${openCoarraysVersion}.tar.gz"
                tarFileWeb="https://github.com/sourceryinstitute/OpenCoarrays/releases/download/${openCoarraysVersion}/${tarFileName}"
                if command -v wget >/dev/null 2>&1; then
                    wget -P "${ParaMonte_REQ_DIR}/.." "${tarFileWeb}"
                    verify $? "download of the prerequisites"
                else
                    echo >&2
                    echo >&2 "-- ${BUILD_NAME} - ${fatal}: The GNU Wget software is required to download the ParaMonte prerequisites."
                    echo >&2 "-- ${BUILD_NAME} - ${fatal}: Please install this wget application and rerun the install/build script."
                    echo >&2
                    echo >&2 "Gracefully exiting..."
                    exit 1
                fi

                #tarFileName="prerequisites.tar.gz"
                #cp -rv "${ParaMonte_ROOT_DIR}/auxil/${tarFileName}" "${ParaMonte_REQ_DIR}/../"
                #verify $? "installation setup of prerequisites"

                (cd "${ParaMonte_REQ_DIR}/../" && tar xvzf "${tarFileName}" -C prerequisites --strip-components 1)
                verify $? "unpacking of prerequisites"
                # chmod +x -R "${ParaMonte_REQ_DIR}"

            else

                echo >&2 "-- ${BUILD_NAME} - ${warning}: The ParaMonte installation will proceed with no guarantee of success."
                echo >&2
                #exit 1

            fi
        fi

    fi

    ################################################################################################################################
    #### install prereqs
    ################################################################################################################################

    if [ "${answer}" = "y" ]; then

        ############################################################################################################################
        #### check brew installation on macOS
        ############################################################################################################################

        if [ "${isMacOS}" = "true" ]; then
            if command -v brew >/dev/null 2>&1; then
                brewCompilerPath="$(command -v brew)"
                echo >&2 "-- ${BUILD_NAME} - Homebrew detected at: ${brewCompilerPath}"
            else
                echo >&2 "-- ${BUILD_NAME} - Homebrew missing. Installing Homebrew..."
                (xcode-select --install && /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)" </dev/null )
                (command -v brew >/dev/null 2>&1 )
                verify $? "installation of Homebrew"
            fi
        fi

        ############################################################################################################################
        #### check cmake installation
        ############################################################################################################################

        if [ "${cmakeInstallEnabled}" = "true" ]; then
            if [ -z "${CMAKE_LOCAL_INSTALLATION_PATH+x}" ] && [ -f "${CMAKE_LOCAL_INSTALLATION_PATH}" ]; then
                #\todo the cmake version extraction needs enhancement.
                echo >&2 "-- ${BUILD_NAME} - cmake ${cmakeVersionParaMonteCompatible} detected."
            else
                echo >&2 "-- ${BUILD_NAME} - cmake ${cmakeVersionParaMonteCompatible} missing."
                echo >&2 "-- ${BUILD_NAME} - installing the prerequisites...this can take a while."
                if [ "${isMacOS}" = "true" ]; then
                    (brew install cmake && brew link cmake )
                    (command -v cmake >/dev/null 2>&1 )
                    verify $? "installation of cmake"
                else
                    chmod +x "${ParaMonte_REQ_DIR}/install.sh"
                    (cd ${ParaMonte_REQ_DIR} && yes | ./install.sh --yes-to-all --package cmake --install-version ${cmakeVersionParaMonteCompatible} )
                    verify $? "installation of cmake"
                    if [ -d "${ParaMonte_REQ_INSTALL_DIR}" ]; then
                        CMAKE_LOCAL_INSTALLATION_PATH="$(find "${ParaMonte_REQ_INSTALL_DIR}/" -path **/bin/cmake)"
                        if [ -f "${CMAKE_LOCAL_INSTALLATION_PATH}" ]; then
                            CMAKE_LOCAL_INSTALLATION_BIN_DIR="$(dirname "${CMAKE_LOCAL_INSTALLATION_PATH}")"
                            export CMAKE_LOCAL_INSTALLATION_BIN_DIR
                            export CMAKE_LOCAL_INSTALLATION_PATH
                            if [[ ":$PATH:" != *":${CMAKE_LOCAL_INSTALLATION_BIN_DIR}:"* ]]; then
                                PATH="${CMAKE_LOCAL_INSTALLATION_BIN_DIR}:${PATH}"
                                export PATH
                            fi
                            cmakeVersion="$(${CMAKE_LOCAL_INSTALLATION_PATH} --version)"
                            cmakeVersionArray=($cmakeVersion)
                            cmakeVersion="${cmakeVersionArray[2]}"
                            echo >&2 "-- ${BUILD_NAME} -       cmake binary path: ${CMAKE_LOCAL_INSTALLATION_PATH}"
                            echo >&2 "-- ${BUILD_NAME} - cmake installed version: ${cmakeVersion}"
                        else
                            unset CMAKE_LOCAL_INSTALLATION_PATH
                            unset CMAKE_LOCAL_INSTALLATION_BIN_DIR
                        fi
                    fi
                fi
            fi
        fi

        ############################################################################################################################
        #### check mpi installation
        ############################################################################################################################

        localMpiInstallationDetected=false
        if [ "${isMacOS}" = "true" ]; then
            CURRENT_PKG="the Open-MPI library"
            echo >&2 "-- ${BUILD_NAME} - ${CURRENT_PKG} missing."
            echo >&2 "-- ${BUILD_NAME} - installing the prerequisites...this can take a while."
            (brew install open-mpi && brew link open-mpi )
            (command -v mpiexec >/dev/null 2>&1 )
            verify $? "installation of ${CURRENT_PKG}"
            MPIEXEC_PATH="$(command -v mpiexec)"
        else
            CURRENT_PKG="the MPICH library"
            if [ "${mpiInstallEnabled}" = "true" ]; then
                MPIEXEC_PATH="${MPI_LOCAL_INSTALLATION_BIN_DIR}/mpiexec"
                if [[ -f "${MPIEXEC_PATH}" ]]; then
                    echo >&2 "-- ${BUILD_NAME} - Local installation of ${CURRENT_PKG} detected: ${MPI_LOCAL_INSTALLATION_BIN_DIR}"
                    localMpiInstallationDetected=true
                    if [[ ":$PATH:" != *":${MPI_LOCAL_INSTALLATION_BIN_DIR}:"* ]]; then
                        PATH="${MPI_LOCAL_INSTALLATION_BIN_DIR}:${PATH}"
                    fi
                    if [[ ":$LD_LIBRARY_PATH:" != *":${MPI_LOCAL_INSTALLATION_LIB_DIR}:"* ]]; then
                        LD_LIBRARY_PATH="${MPI_LOCAL_INSTALLATION_LIB_DIR}:${LD_LIBRARY_PATH}"
                        export LD_LIBRARY_PATH
                    fi
                    PATH="${MPI_LOCAL_INSTALLATION_LIB_DIR}:${PATH}"
                    export PATH
                else
                    ##########################################################################
                    echo >&2 "-- ${BUILD_NAME} - ${CURRENT_PKG} missing."
                    echo >&2 "-- ${BUILD_NAME} - installing the prerequisites...this can take a while."
                    chmod +x "${ParaMonte_REQ_DIR}/install.sh"
                    (cd ${ParaMonte_REQ_DIR} && yes | ./install.sh --yes-to-all --package mpich --install-version ${mpichVersionOpenCoarrays} ) ||
                    {
                        if [ -z ${GCC_BOOTSTRAP+x} ]; then
                            if [ "${YES_TO_ALL_DISABLED}" = "true" ]; then
                                answerNotGiven=true
                                while [ "${answerNotGiven}" = "true" ]; do
                                    echo >&2
                                    read -p "-- ${BUILD_NAME} - ${CURRENT_PKG} installation failed. Shall I retry with bootstrap (y/n)? " answer
                                    echo >&2
                                    if [[ $answer == [yY] || $answer == [yY][eE][sS] ]]; then
                                        answer=y
                                        answerNotGiven=false
                                    fi
                                    if [[ $answer == [nN] || $answer == [nN][oO] ]]; then
                                        answer=n
                                        answerNotGiven=false
                                    fi
                                    if [ "${answerNotGiven}" = "true" ]; then
                                        echo >&2 "-- ${BUILD_NAME} - please enter either y or n"
                                    fi
                                done
                            else
                                echo >&2
                                echo >&2 "-- ${BUILD_NAME} - ${CURRENT_PKG} installation failed. Shall I retry with bootstrap (y/n)? y"
                                echo >&2
                                answer=y
                            fi
                            echo >&2
                            if [ "${answer}" = "y" ]; then
                                GCC_BOOTSTRAP="--bootstrap"
                                chmod +x "${ParaMonte_REQ_DIR}/install.sh"
                                (cd ${ParaMonte_REQ_DIR} && yes | ./install.sh --yes-to-all ${GCC_BOOTSTRAP} )
                                verify $? "installation of ${CURRENT_PKG}"
                            else
                                verify 1 "installation of ${CURRENT_PKG}"
                            fi
                        else
                            verify 1 "installation of ${CURRENT_PKG}"
                        fi
                    }
                    ##########################################################################
                fi
            fi
        fi

        ############################################################################################################################
        #### check gnu
        ############################################################################################################################

        CURRENT_PKG="the GNU compiler collection"
        if [ "${gnuInstallEnabled}" = "true" ]; then # || [ "${localMpiInstallationDetected}" = "false" ]); then
            Fortran_COMPILER_PATH="${GNU_LOCAL_INSTALLATION_BIN_DIR}/gfortran"
            if [[ -f "${Fortran_COMPILER_PATH}" ]]; then
                echo >&2 "-- ${BUILD_NAME} - Local installation of ${CURRENT_PKG} detected: ${GNU_LOCAL_INSTALLATION_BIN_DIR}"
                if [[ ":$PATH:" != *":${GNU_LOCAL_INSTALLATION_LIB_DIR}:"* ]]; then
                    PATH="${GNU_LOCAL_INSTALLATION_LIB_DIR}:${PATH}"
                    export PATH
                fi
                if [[ ":$LD_LIBRARY_PATH:" != *":${GNU_LOCAL_INSTALLATION_LIB_DIR}:"* ]]; then
                    LD_LIBRARY_PATH="${GNU_LOCAL_INSTALLATION_LIB_DIR}:${LD_LIBRARY_PATH}"
                    export LD_LIBRARY_PATH
                fi
                #gnuInstallEnabled=false
            else
                ##########################################################################
                echo >&2 "-- ${BUILD_NAME} - ${CURRENT_PKG} missing."
                echo >&2 "-- ${BUILD_NAME} - installing the prerequisites...this can take a while."
                if [ "${isMacOS}" = "true" ]; then
                    (brew install gcc && brew link gcc )
                    (command -v gfortran >/dev/null 2>&1 )
                    verify $? "installation of ${CURRENT_PKG}"
                    Fortran_COMPILER_PATH="$(command -v gfortran)"
                else
                    chmod +x "${ParaMonte_REQ_DIR}/install.sh"
                    (cd ${ParaMonte_REQ_DIR} && yes | ./install.sh --yes-to-all ${GCC_BOOTSTRAP} ) ||
                    {
                        if [ -z ${GCC_BOOTSTRAP+x} ]; then
                            if [ "${YES_TO_ALL_DISABLED}" = "true" ]; then
                                answerNotGiven=true
                                while [ "${answerNotGiven}" = "true" ]; do
                                    echo >&2
                                    read -p "-- ${BUILD_NAME} - ${CURRENT_PKG} installation failed. Shall I retry with bootstrap (y/n)? " answer
                                    echo >&2
                                    if [[ $answer == [yY] || $answer == [yY][eE][sS] ]]; then
                                        answer=y
                                        answerNotGiven=false
                                    fi
                                    if [[ $answer == [nN] || $answer == [nN][oO] ]]; then
                                        answer=n
                                        answerNotGiven=false
                                    fi
                                    if [ "${answerNotGiven}" = "true" ]; then
                                        echo >&2 "-- ${BUILD_NAME} - please enter either y or n"
                                    fi
                                done
                            else
                                echo >&2
                                echo >&2 "-- ${BUILD_NAME} - ${CURRENT_PKG} installation failed. Shall I retry with bootstrap (y/n)? y"
                                echo >&2
                                answer=y
                            fi
                            echo >&2
                            if [ "${answer}" = "y" ]; then
                                GCC_BOOTSTRAP="--bootstrap"
                                chmod +x "${ParaMonte_REQ_DIR}/install.sh"
                                (cd ${ParaMonte_REQ_DIR} && yes | ./install.sh --yes-to-all ${GCC_BOOTSTRAP} )
                                verify $? "installation of ${CURRENT_PKG}"
                            else
                                verify 1 "installation of ${CURRENT_PKG}"
                            fi
                        else
                            verify 1 "installation of ${CURRENT_PKG}"
                        fi
                    }
                fi
                ##########################################################################
            fi
        fi

        ############################################################################################################################
        #### check caf
        ############################################################################################################################

        CURRENT_PKG="the OpenCoarrays compiler wrapper"
        if [ "${CAF_ENABLED}" = "true" ] && [ "${cafInstallEnabled}" = "true" ]; then
            if [[ -f "${CAF_LOCAL_INSTALLATION_WRAPPER_PATH}" ]]; then
                echo >&2 "-- ${BUILD_NAME} - Local installation of ${CURRENT_PKG} detected: ${CAF_LOCAL_INSTALLATION_WRAPPER_PATH}"
            else
                ##########################################################################
                echo >&2 "-- ${BUILD_NAME} - ${CURRENT_PKG} missing."
                echo >&2 "-- ${BUILD_NAME} - installing the prerequisites...this can take a while."
                if [ "${isMacOS}" = "true" ]; then
                    (brew install opencoarrays && brew link opencoarrays )
                    (command -v caf >/dev/null 2>&1 )
                    verify $? "installation of ${CURRENT_PKG}"
                    CAF_LOCAL_INSTALLATION_WRAPPER_PATH="$(command -v caf)"
                else
                    chmod +x "${ParaMonte_REQ_DIR}/install.sh"
                    (cd ${ParaMonte_REQ_DIR} && yes | ./install.sh ${GCC_BOOTSTRAP} --yes-to-all) ||
                    {
                        if [ -z ${GCC_BOOTSTRAP+x} ]; then
                            if [ "${YES_TO_ALL_DISABLED}" = "true" ]; then
                                answerNotGiven=true
                                while [ "${answerNotGiven}" = "true" ]; do
                                    echo >&2
                                    read -p "-- ${BUILD_NAME} - ${CURRENT_PKG} installation failed. Shall I retry with bootstrap (y/n)? " answer
                                    echo >&2
                                    if [[ $answer == [yY] || $answer == [yY][eE][sS] ]]; then
                                        answer=y
                                        answerNotGiven=false
                                    fi
                                    if [[ $answer == [nN] || $answer == [nN][oO] ]]; then
                                        answer=n
                                        answerNotGiven=false
                                    fi
                                    if [ "${answerNotGiven}" = "true" ]; then
                                        echo >&2 "-- ${BUILD_NAME} - please enter either y or n"
                                    fi
                                done
                            else
                                echo >&2
                                echo >&2 "-- ${BUILD_NAME} - ${CURRENT_PKG} installation failed. Shall I retry with bootstrap (y/n)? y"
                                echo >&2
                                answer=y
                            fi
                            echo >&2
                            if [ "${answer}" = "y" ]; then
                                GCC_BOOTSTRAP="--bootstrap"
                                chmod +x "${ParaMonte_REQ_DIR}/install.sh"
                                (cd ${ParaMonte_REQ_DIR} && yes | ./install.sh --yes-to-all ${GCC_BOOTSTRAP} )
                                verify $? "installation of ${CURRENT_PKG}"
                            else
                                verify 1 "installation of ${CURRENT_PKG}"
                            fi
                        else
                            verify 1 "installation of ${CURRENT_PKG}"
                        fi
                    }
                fi
                ##########################################################################
            fi
            Fortran_COMPILER_PATH="${CAF_LOCAL_INSTALLATION_WRAPPER_PATH}"
        fi

        if [ -f "${CAF_LOCAL_INSTALLATION_SETUP_FILE}" ]; then
            source "${CAF_LOCAL_INSTALLATION_SETUP_FILE}"
            # source "${SETUP_FILE_PATH}"
            # echo "" >> ${SETUP_FILE_PATH}
            # echo "source ${CAF_LOCAL_INSTALLATION_SETUP_FILE}" >> ${SETUP_FILE_PATH}
            # echo "" >> ${SETUP_FILE_PATH}
        fi

    fi

fi

# check one last time if Fortran compiler path exists, if not unset

if [[ ! -f "${Fortran_COMPILER_PATH}" ]]; then
    unset Fortran_COMPILER_PATH
fi

####################################################################################################################################
#### set up PATH & LD_LIBRARY_PATH
####################################################################################################################################

if [ -d "${CMAKE_LOCAL_INSTALLATION_BIN_DIR}" ]; then
    if [ -z ${PATH+x} ]; then
        PATH="${CMAKE_LOCAL_INSTALLATION_BIN_DIR}"
    else
        if [[ ":$PATH:" != *":${CMAKE_LOCAL_INSTALLATION_BIN_DIR}:"* ]]; then
            PATH="${CMAKE_LOCAL_INSTALLATION_BIN_DIR}:${PATH}"
        fi
    fi
fi

####################################################################################################################################
##### set up setup.sh file
####################################################################################################################################

# SETUP_FILE_PATH="${ParaMonte_ROOT_DIR}/build/setup.sh"
# export SETUP_FILE_PATH
#
# {
#     echo "if ! [[ \"${PATH}\" =~ .*\"${GNU_LOCAL_INSTALLATION_BIN_DIR}\".* ]]; then"
#     echo "    export PATH=\"${GNU_LOCAL_INSTALLATION_BIN_DIR}:${PATH}\""
#     echo "fi"
#     echo ""
#     echo "if ! [[ \"${PATH}\" =~ .*\"${MPI_LOCAL_INSTALLATION_BIN_DIR}\".* ]]; then"
#     echo "    export PATH=\"${MPI_LOCAL_INSTALLATION_BIN_DIR}:${PATH}\""
#     echo "fi"
#     echo ""
#     echo "if ! [[ \"${PATH}\" =~ .*\"${CAF_LOCAL_INSTALLATION_BIN_DIR}\".* ]]; then"
#     echo "    export PATH=\"${CAF_LOCAL_INSTALLATION_BIN_DIR}:${PATH}\""
#     echo "fi"
#     echo ""
#     echo "if ! [[ \"${LD_LIBRARY_PATH}\" =~ .*\"${GNU_LOCAL_INSTALLATION_LIB_DIR}\".* ]]; then"
#     echo "    export PATH=\"${GNU_LOCAL_INSTALLATION_LIB_DIR}:${LD_LIBRARY_PATH}\""
#     echo "fi"
#     echo ""
#     echo "if ! [[ \"${LD_LIBRARY_PATH}\" =~ .*\"${MPI_LOCAL_INSTALLATION_LIB_DIR}\".* ]]; then"
#     echo "    export PATH=\"${MPI_LOCAL_INSTALLATION_LIB_DIR}:${LD_LIBRARY_PATH}\""
#     echo "fi"
#     echo ""
#     echo "if ! [[ \"${LD_LIBRARY_PATH}\" =~ .*\"${CAF_LOCAL_INSTALLATION_LIB_DIR}\".* ]]; then"
#     echo "    export PATH=\"${CAF_LOCAL_INSTALLATION_LIB_DIR}:${LD_LIBRARY_PATH}\""
#     echo "fi"
# } > ${SETUP_FILE_PATH}

if [ "${PMCS}" = "gnu" ] && ( [[ "${Fortran_COMPILER_PATH}" =~ .*"build/prerequisites/prerequisites/installations".* ]] || [[ "${MPIEXEC_PATH}" =~ .*"build/prerequisites/prerequisites/installations".* ]] ) && ! [ "${isMacOS}" = "true" ]; then

    SETUP_FILE_PATH="${ParaMonte_ROOT_DIR}/build/setup.sh"
    export SETUP_FILE_PATH

    # add the default OpenCoarrays setup file
    {
    echo "# ParaMonte runtime environment setup script."
    echo "# Source this Bash script in your Bash environment like,"
    echo "#     source ./setup.sh"
    echo "# before compiling your source files and linking with ParaMonte library."
    echo ""
    } > ${SETUP_FILE_PATH}
    chmod +x ${SETUP_FILE_PATH}

    if [ -f "${CAF_LOCAL_INSTALLATION_SETUP_FILE}" ]; then
        {
            echo ""
            echo "source ${CAF_LOCAL_INSTALLATION_SETUP_FILE}"
            echo ""
        } >> ${SETUP_FILE_PATH}
        if [ "${LTYPE}" = "dynamic" ]; then
            {
            echo "if [ -z \${LD_LIBRARY_PATH+x} ]; then"
            echo "    LD_LIBRARY_PATH=."
            echo "else"
            echo "    sourceFileDir=\"\$( cd \"\$( dirname \"\${BASH_SOURCE[0]}\" )\" >/dev/null 2>&1 && pwd )\""
            echo "    if [[ \":\$LD_LIBRARY_PATH:\" != *\":${sourceFileDir}:\"* ]]; then"
            echo "        LD_LIBRARY_PATH=\"${sourceFileDir}:\${LD_LIBRARY_PATH}\""
            echo "    fi"
            echo "fi"
            echo "export LD_LIBRARY_PATH"
            echo ""
            } >> ${SETUP_FILE_PATH}
        fi
    fi

    #####################################################
    #### add the GNU installation paths to the setup file
    #####################################################

    if [ -d "${GNU_LOCAL_INSTALLATION_BIN_DIR}" ]; then
        if [ -z ${PATH+x} ]; then
            PATH="${GNU_LOCAL_INSTALLATION_BIN_DIR}"
        else
            if [[ ":$PATH:" != *":${GNU_LOCAL_INSTALLATION_BIN_DIR}:"* ]]; then
                PATH="${GNU_LOCAL_INSTALLATION_BIN_DIR}:${PATH}"
            fi
        fi
        {
            echo "if [ -z \${PATH+x} ]; then"
            echo "    export PATH=\"${GNU_LOCAL_INSTALLATION_BIN_DIR}\""
            echo "else"
            echo "    if [[ \":\$PATH:\" != *\":${GNU_LOCAL_INSTALLATION_BIN_DIR}:\"* ]]; then"
            echo "        export PATH=\"${GNU_LOCAL_INSTALLATION_BIN_DIR}:\${PATH}\""
            echo "    fi"
            echo "fi"
        } >> ${SETUP_FILE_PATH}
    fi

    if [ -d "${GNU_LOCAL_INSTALLATION_LIB_DIR}" ]; then
        if [ -z ${LD_LIBRARY_PATH+x} ]; then
            LD_LIBRARY_PATH="${GNU_LOCAL_INSTALLATION_LIB_DIR}"
        else
            if [[ ":$LD_LIBRARY_PATH:" != *":${GNU_LOCAL_INSTALLATION_LIB_DIR}:"* ]]; then
                LD_LIBRARY_PATH="${GNU_LOCAL_INSTALLATION_LIB_DIR}:${LD_LIBRARY_PATH}"
            fi
        fi
        {
            echo "if [ -z \${LD_LIBRARY_PATH+x} ]; then"
            echo "    export LD_LIBRARY_PATH=\"${GNU_LOCAL_INSTALLATION_LIB_DIR}\""
            echo "else"
            echo "    if [[ \":\$LD_LIBRARY_PATH:\" != *\":${GNU_LOCAL_INSTALLATION_LIB_DIR}:\"* ]]; then"
            echo "        export LD_LIBRARY_PATH=\"${GNU_LOCAL_INSTALLATION_LIB_DIR}:\${LD_LIBRARY_PATH}\""
            echo "    fi"
            echo "fi"
        } >> ${SETUP_FILE_PATH}
    fi

    if [ -d "${GNU_LOCAL_INSTALLATION_LIB64_DIR}" ]; then
        if [ -z ${LD_LIBRARY_PATH+x} ]; then
            LD_LIBRARY_PATH="${GNU_LOCAL_INSTALLATION_LIB64_DIR}"
        else
            if [[ ":$LD_LIBRARY_PATH:" != *":${GNU_LOCAL_INSTALLATION_LIB64_DIR}:"* ]]; then
                LD_LIBRARY_PATH="${GNU_LOCAL_INSTALLATION_LIB64_DIR}:${LD_LIBRARY_PATH}"
            fi
        fi
        {
            echo "if [ -z \${LD_LIBRARY_PATH+x} ]; then"
            echo "    export LD_LIBRARY_PATH=\"${GNU_LOCAL_INSTALLATION_LIB64_DIR}\""
            echo "else"
            echo "    if [[ \":\$LD_LIBRARY_PATH:\" != *\":${GNU_LOCAL_INSTALLATION_LIB64_DIR}:\"* ]]; then"
            echo "        export LD_LIBRARY_PATH=\"${GNU_LOCAL_INSTALLATION_LIB64_DIR}:\${LD_LIBRARY_PATH}\""
            echo "    fi"
            echo "fi"
        } >> ${SETUP_FILE_PATH}
    fi

    ##########################
    #### add path to local MPI
    ##########################

    if [ -d "${MPI_LOCAL_INSTALLATION_BIN_DIR}" ]; then
        if [ -z ${PATH+x} ]; then
            PATH="${MPI_LOCAL_INSTALLATION_BIN_DIR}"
        else
            if [[ ":$PATH:" != *":${MPI_LOCAL_INSTALLATION_BIN_DIR}:"* ]]; then
                PATH="${MPI_LOCAL_INSTALLATION_BIN_DIR}:${PATH}"
            fi
        fi
        {
            echo "if [ -z \${PATH+x} ]; then"
            echo "    export PATH=\"${MPI_LOCAL_INSTALLATION_BIN_DIR}\""
            echo "else"
            echo "    if [[ \":\$PATH:\" != *\":${MPI_LOCAL_INSTALLATION_BIN_DIR}:\"* ]]; then"
            echo "        export PATH=\"${MPI_LOCAL_INSTALLATION_BIN_DIR}:\${PATH}\""
            echo "    fi"
            echo "fi"
        } >> ${SETUP_FILE_PATH}
    fi

    if [ -d "${MPI_LOCAL_INSTALLATION_LIB_DIR}" ]; then
        if [ -z ${LD_LIBRARY_PATH+x} ]; then
            LD_LIBRARY_PATH="${MPI_LOCAL_INSTALLATION_LIB_DIR}"
        else
            if [[ ":$LD_LIBRARY_PATH:" != *":${MPI_LOCAL_INSTALLATION_LIB_DIR}:"* ]]; then
                LD_LIBRARY_PATH="${MPI_LOCAL_INSTALLATION_LIB_DIR}:${LD_LIBRARY_PATH}"
            fi
        fi
        {
            echo "if [ -z \${LD_LIBRARY_PATH+x} ]; then"
            echo "    export LD_LIBRARY_PATH=\"${MPI_LOCAL_INSTALLATION_LIB_DIR}\""
            echo "else"
            echo "    if [[ \":\$LD_LIBRARY_PATH:\" != *\":${MPI_LOCAL_INSTALLATION_LIB_DIR}:\"* ]]; then"
            echo "        export LD_LIBRARY_PATH=\"${MPI_LOCAL_INSTALLATION_LIB_DIR}:\${LD_LIBRARY_PATH}\""
            echo "    fi"
            echo "fi"
        } >> ${SETUP_FILE_PATH}
    fi

    if [ -d "${MPI_LOCAL_INSTALLATION_LIB64_DIR}" ]; then
        if [ -z ${LD_LIBRARY_PATH+x} ]; then
            LD_LIBRARY_PATH="${MPI_LOCAL_INSTALLATION_LIB64_DIR}"
        else
            if [[ ":$LD_LIBRARY_PATH:" != *":${MPI_LOCAL_INSTALLATION_LIB64_DIR}:"* ]]; then
                LD_LIBRARY_PATH="${MPI_LOCAL_INSTALLATION_LIB64_DIR}:${LD_LIBRARY_PATH}"
            fi
        fi
        {
            echo "if [ -z \${LD_LIBRARY_PATH+x} ]; then"
            echo "    export LD_LIBRARY_PATH=\"${MPI_LOCAL_INSTALLATION_LIB64_DIR}\""
            echo "else"
            echo "    if [[ \":\$LD_LIBRARY_PATH:\" != *\":${MPI_LOCAL_INSTALLATION_LIB64_DIR}:\"* ]]; then"
            echo "        export LD_LIBRARY_PATH=\"${MPI_LOCAL_INSTALLATION_LIB64_DIR}:\${LD_LIBRARY_PATH}\""
            echo "    fi"
            echo "fi"
        } >> ${SETUP_FILE_PATH}
    fi

    ############################
    #### add path to caf wrapper
    ############################

    if [ -d "${CAF_LOCAL_INSTALLATION_BIN_DIR}" ]; then
        if [ -z ${PATH+x} ]; then
            PATH="${CAF_LOCAL_INSTALLATION_BIN_DIR}"
        else
            if [[ ":$PATH:" != *":${CAF_LOCAL_INSTALLATION_BIN_DIR}:"* ]]; then
                PATH="${CAF_LOCAL_INSTALLATION_BIN_DIR}:${PATH}"
            fi
        fi
        {
            echo "if [ -z \${PATH+x} ]; then"
            echo "    export PATH=\"${CAF_LOCAL_INSTALLATION_BIN_DIR}\""
            echo "else"
            echo "    if [[ \":\$PATH:\" != *\":${CAF_LOCAL_INSTALLATION_BIN_DIR}:\"* ]]; then"
            echo "        export PATH=\"${CAF_LOCAL_INSTALLATION_BIN_DIR}:\${PATH}\""
            echo "    fi"
            echo "fi"
        } >> ${SETUP_FILE_PATH}
    fi

    if [ -d "${CAF_LOCAL_INSTALLATION_LIB_DIR}" ]; then
        if [ -z ${LD_LIBRARY_PATH+x} ]; then
            LD_LIBRARY_PATH="${CAF_LOCAL_INSTALLATION_LIB_DIR}"
        else
            if [[ ":$LD_LIBRARY_PATH:" != *":${CAF_LOCAL_INSTALLATION_LIB_DIR}:"* ]]; then
                LD_LIBRARY_PATH="${CAF_LOCAL_INSTALLATION_LIB_DIR}:${LD_LIBRARY_PATH}"
            fi
        fi
        {
            echo "if [ -z \${LD_LIBRARY_PATH+x} ]; then"
            echo "    export LD_LIBRARY_PATH=\"${CAF_LOCAL_INSTALLATION_LIB_DIR}\""
            echo "else"
            echo "    if [[ \":\$LD_LIBRARY_PATH:\" != *\":${CAF_LOCAL_INSTALLATION_LIB_DIR}:\"* ]]; then"
            echo "        export LD_LIBRARY_PATH=\"${CAF_LOCAL_INSTALLATION_LIB_DIR}:\${LD_LIBRARY_PATH}\""
            echo "    fi"
            echo "fi"
        } >> ${SETUP_FILE_PATH}
    fi

    export PATH
    export LD_LIBRARY_PATH

fi

####################################################################################################################################
#### set Fortran compiler version one last time
####################################################################################################################################

echo >&2
echo >&2 "-- ${BUILD_NAME}Compiler - Fortran compiler path: ${Fortran_COMPILER_PATH}"
echo >&2 "-- ${BUILD_NAME}Compiler - MPI mpiexec path: ${MPIEXEC_PATH}"
echo >&2

if [ "${PMCS}" = "gnu" ] || [ "${COMPILER_VERSION}" = "unknown" ]; then

    #cd ./auxil/

    LANG=Fortran
    isUnknownVersion=false

    #tempDir=$(mktemp -d --tmpdir="${ParaMonte_ROOT_BUILD_DIR}/") || mkdir -p "${ParaMonte_ROOT_BUILD_DIR}"
    tempDir=$(mktemp -d "${TMPDIR:-/tmp}/cversion.XXXXXXXXX")
    echo >&2 "-- ${BUILD_NAME}Compiler - changing directory to: ${tempDir}"
    cd "${tempDir}"
    cp "${ParaMonte_ROOT_DIR}/auxil/getCompilerVersion.f90" "./getCompilerVersion.f90"
    if ${Fortran_COMPILER_PATH} getCompilerVersion.f90 -o getCompilerVersion.exe; then

        chmod +x getCompilerVersion.exe
        ./getCompilerVersion.exe -outdir "" && {
            COMPILER_VERSION=$(head -n 1 getCompilerVersion.tmp)
            echo >&2 "-- ${BUILD_NAME}Compiler - ${PMCS} ${LANG} compiler version: ${COMPILER_VERSION}"
            isParaMonteCompatibleCompiler=$(head -n 1 isParaMonteCompatibleCompiler.tmp)
            isGfortran10=$(head -n 1 isGfortran10.tmp); export isGfortran10
            if [ "$isParaMonteCompatibleCompiler" = "true" ]; then
                echo >&2 "-- ${BUILD_NAME}Compiler - ${PMCS} ${LANG} compiler is ParaMonte compatible!"
            else
                echo >&2 "-- ${BUILD_NAME}Compiler - ${SUITE} ${LANG} compiler is not ParaMonte compatible..."
                echo >&2 "-- ${BUILD_NAME}Compiler - ParaMonte installation failed."
                echo >&2
                echo >&2 "-- ${BUILD_NAME} - gracefully exiting."
                echo >&2
                exit 1
            fi
            rm *.tmp *.exe
        } || {
            isUnknownVersion=true
        }

    else

        isUnknownVersion=true

    fi

    echo >&2 "-- ${BUILD_NAME}Compiler - changing directory to: ${ParaMonte_ROOT_DIR}"
    cd "${ParaMonte_ROOT_DIR}"
    rm -rf "${tempDir}"

    if [ "${isUnknownVersion}" = "true" ]; then
        echo >&2 "-- ${BUILD_NAME}Compiler - failed to detect the ${PMCS} ${LANG} compiler version. skipping..."
        COMPILER_VERSION="unknown"
    fi

fi

####################################################################################################################################
#### determine the MPI library name to used in the ParaMonte library name and build directory
####################################################################################################################################

unset MPILIB_NAME
if [ -f "${MPIEXEC_PATH}" ]; then

    mpiVersionInfo="$(${MPIEXEC_PATH} --version)"

    if [[ "${mpiVersionInfo}" =~ .*"Intel".* ]]; then
        MPILIB_NAME="impi"
    elif [[ "${mpiVersionInfo}" =~ .*[oO][pP][eE][nN][rR][tT][eE].* ]] || [[ "${mpiVersionInfo}" =~ .*[oO][pP][eE][nN]-?[mM][pP][iI].* ]]; then
        MPILIB_NAME="openmpi"
    elif [[ "${mpiVersionInfo}" =~ .*[mM][pP][iI][cC][hH].* ]]; then
        MPILIB_NAME="mpich"
    else # look for mpichversion
        mpichversion_PATH="$(dirname ${MPIEXEC_PATH})"/mpichversion
        if [ -f "${mpichversion_PATH}" ] && [[ "$(mpichversion_PATH)" =~ .*[mM][pP][iI][cC][hH].* ]]; then
            MPILIB_NAME="mpich"
        fi
    fi
fi

if [ -z ${MPILIB_NAME+x} ] && [ "${MPI_ENABLED}" = "true" ]; then
    MPILIB_NAME="mpi"
    echo >&2
    echo >&2 "${pmwarn} The make of the MPI library could be identified."
    echo >&2 "${pmwarn} The MPI library's behavior does not match the Intel, MPICH, or OpenMPI libraries."
    echo >&2 "${pmwarn} The ParaMonte library name will be suffixed with the generic \"${MPILIB_NAME}\" label."
    echo >&2
fi

####################################################################################################################################
# set ParaMonte build dir
####################################################################################################################################

export PMCS
export COMPILER_VERSION
echo >&2 "-- ${BUILD_NAME} - selected compiler suite: ${PMCS}"
echo >&2 "-- ${BUILD_NAME} - selected compiler version: ${COMPILER_VERSION}"

unset PARALLELIZATION_DIR
if [ "${OMP_ENABLED}" = "true" ]; then PARALLELIZATION_DIR=${PARALLELIZATION_DIR}omp; fi
if [ "${MPI_ENABLED}" = "true" ]; then PARALLELIZATION_DIR=${PARALLELIZATION_DIR}${MPILIB_NAME}; fi
if [ "${CAF_ENABLED}" = "true" ]; then PARALLELIZATION_DIR=${PARALLELIZATION_DIR}caf${CAFTYPE}; fi
if [ -z ${PARALLELIZATION_DIR+x} ]; then PARALLELIZATION_DIR=serial; fi
export PARALLELIZATION_DIR

unset MEMORY_ALLOCATION
if [ "${HEAP_ARRAY_ENABLED}" = "true" ]; then MEMORY_ALLOCATION="heap"; fi
if [ "${HEAP_ARRAY_ENABLED}" = "false" ]; then MEMORY_ALLOCATION="stack"; fi
if [ -z ${MEMORY_ALLOCATION+x} ]; then
    echo >&2 "-- ${BUILD_NAME} - FATAL: HEAP_ARRAY_ENABLED must be set to either true or false."
    echo >&2 "-- ${BUILD_NAME} - FATAL: you have provided HEAP_ARRAY_ENABLED=${HEAP_ARRAY_ENABLED}"
    echo >&2
    echo >&2 "-- ${BUILD_NAME} - gracefully exiting."
    echo >&2
    exit 1
fi
export MEMORY_ALLOCATION

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set ParaMonte library build directories
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

echo >&2
echo >&2 "-- ${BUILD_NAME} - setting up build directories..."
echo >&2

ParaMonte_BLD_DIR="${ParaMonte_ROOT_DIR}/build/${PLATFORM}${ARCHITECTURE}/${PMCS}/${COMPILER_VERSION}/${BTYPE}/${LTYPE}/${MEMORY_ALLOCATION}/${PARALLELIZATION_DIR}"
if [ ${INTERFACE_LANGUAGE} = "c" ]; then ParaMonte_BLD_DIR="${ParaMonte_BLD_DIR}/C"; fi
if [ ${INTERFACE_LANGUAGE} = "c++" ]; then ParaMonte_BLD_DIR="${ParaMonte_BLD_DIR}/C++"; fi
if [ ${INTERFACE_LANGUAGE} = "fortran" ]; then ParaMonte_BLD_DIR="${ParaMonte_BLD_DIR}/Fortran"; fi
if [ ${INTERFACE_LANGUAGE} = "matlab" ]; then ParaMonte_BLD_DIR="${ParaMonte_BLD_DIR}/MATLAB"; fi
if [ ${INTERFACE_LANGUAGE} = "python" ]; then ParaMonte_BLD_DIR="${ParaMonte_BLD_DIR}/Python"; fi
if [ -z ${CFI_ENABLED+x} ]; then
    echo >&2 "-- ${BUILD_NAME} - FATAL: CFI_ENABLED must be set to either true or false."
    echo >&2 "-- ${BUILD_NAME} - FATAL: you have provided CFI_ENABLED=${CFI_ENABLED}"
    echo >&2
    echo >&2 "-- ${BUILD_NAME} - gracefully exiting."
    echo >&2
    exit 1
fi

if [ "${CODECOV_ENABLED}" = "true" ]; then
    ParaMonte_BLD_DIR="${ParaMonte_BLD_DIR}/codecov"
fi

export ParaMonte_BLD_DIR

echo >&2 "-- ${BUILD_NAME} - build directory: ${ParaMonte_BLD_DIR}"
if [ -d "${ParaMonte_BLD_DIR}" ]; then
    echo >&2 "-- ${BUILD_NAME} - ParaMonte build directory already exists. skipping..."
else
    echo >&2 "-- ${BUILD_NAME} - generating ParaMonte build directory..."
    mkdir -p "${ParaMonte_BLD_DIR}"
fi
echo >&2 "-- ${BUILD_NAME} - all generated build files will be stored at: ${ParaMonte_BLD_DIR}"

# set object/module/lib files directories

ParaMonte_OBJ_DIR="${ParaMonte_BLD_DIR}/obj"; export ParaMonte_OBJ_DIR
ParaMonte_MOD_DIR="${ParaMonte_BLD_DIR}/mod"; export ParaMonte_MOD_DIR
ParaMonte_LIB_DIR="${ParaMonte_BLD_DIR}/lib"; export ParaMonte_LIB_DIR

echo >&2 "-- ${BUILD_NAME} - ParaMonte object files directory: ${ParaMonte_OBJ_DIR}"
echo >&2 "-- ${BUILD_NAME} - ParaMonte module files directory: ${ParaMonte_MOD_DIR}"
echo >&2 "-- ${BUILD_NAME} - ParaMonte library files directory: ${ParaMonte_LIB_DIR}"

# make bin directory

ParaMonte_BIN_DIR="${ParaMonte_ROOT_DIR}/bin"
echo >&2 "-- ${BUILD_NAME} - ParaMonte binaries directory: ${ParaMonte_BIN_DIR}"
if [[ -d "${ParaMonte_BIN_DIR}" ]]; then
    echo >&2 "-- ${BUILD_NAME} - ParaMonte binaries directory already exists. skipping..."
else
    echo >&2 "-- ${BUILD_NAME} - generating ParaMonte binaries directory..."
    mkdir "${ParaMonte_BIN_DIR}/"
fi
export ParaMonte_BIN_DIR

####################################################################################################################################
# set ParaMonte library source directories
####################################################################################################################################

         ParaMonteExample_SRC_DIR="${ParaMonte_ROOT_DIR}/example"
       ParaMonteInterface_SRC_DIR="${ParaMonte_ROOT_DIR}/src/interface"
      ParaMonteInterfaceC_SRC_DIR="${ParaMonteInterface_SRC_DIR}/C"
 ParaMonteInterfaceMATLAB_SRC_DIR="${ParaMonteInterface_SRC_DIR}/MATLAB"
 ParaMonteInterfacePython_SRC_DIR="${ParaMonteInterface_SRC_DIR}/Python"
ParaMonteInterfaceFortran_SRC_DIR="${ParaMonteInterface_SRC_DIR}/Fortran"
      ParaMonteMATLABTest_SRC_DIR="${ParaMonteInterfaceMATLAB_SRC_DIR}/test"
      ParaMontePythonTest_SRC_DIR="${ParaMonteInterfacePython_SRC_DIR}/test"

export ParaMonteExample_SRC_DIR
export ParaMonteMATLABTest_SRC_DIR
export ParaMontePythonTest_SRC_DIR
export ParaMonteInterface_SRC_DIR
export ParaMonteInterfaceC_SRC_DIR
export ParaMonteInterfaceMATLAB_SRC_DIR
export ParaMonteInterfacePython_SRC_DIR
export ParaMonteInterfaceFortran_SRC_DIR

####################################################################################################################################
# check for MATLAB's existence
####################################################################################################################################

unset MATLAB_ROOT_DIR_OPTION # it is imperative to nullify this MATLAB variable for all language builds at all times

if [ "${INTERFACE_LANGUAGE}" = "matlab" ] && [ "${LTYPE}" = "dynamic" ] && [ "${CFI_ENABLED}" = "true" ]; then

    echo >&2
    echo >&2 "-- ${BUILD_NAME}MATLAB - searching for a MATLAB installation on your system..."

    unset MATLAB_ROOT_DIR
    unset MATLAB_EXE_PATH
    unset MATLAB_BIN_DIR

    if command -v matlab >/dev/null 2>&1; then
        MATLAB_EXE_PATH=$(command -v matlab)
        MATLAB_BIN_DIR=$(dirname $MATLAB_EXE_PATH)
        MATLAB_ROOT_DIR=$(dirname $MATLAB_BIN_DIR)
        # echo >&2 "-- ${BUILD_NAME} - MATLAB detected at: ${MATLAB_EXE_PATH}"
    else
        echo >&2 "-- ${BUILD_NAME} - ${warning}: MATLAB could not be found in among the search paths."
        echo >&2 "-- ${BUILD_NAME} - ${warning}: searching for MATLAB in the default installation directories..."
    fi

    if [ -z ${MATLAB_EXE_PATH+x} ]; then
        if [ "${isMacOS}" = "true" ]; then
            INSTALL_LOC_LIST="/Applications/MATLAB_"
        else
            INSTALL_LOC_LIST="/usr/local/MATLAB"
        fi
        MATLAB_VERSION_LIST="R2025b R2025a R2024b R2024a R2023b R2023a R2022b R2022a R2021b R2021a R2020b R2020a R2019b R2019a R2018b R2018a R2017b R2017a"

        for INSTALL_LOC in $INSTALL_LOC_LIST; do
            for MATLAB_VERSION in $MATLAB_VERSION_LIST; do
                if [ "${isMacOS}" = "true" ]; then
                    MATLAB_ROOT_DIR="${INSTALL_LOC_LIST}${MATLAB_VERSION}.app"
                else
                    MATLAB_ROOT_DIR="${INSTALL_LOC_LIST}/${MATLAB_VERSION}"
                fi
                MATLAB_BIN_DIR="${MATLAB_ROOT_DIR}/bin"
                MATLAB_EXE_PATH="${MATLAB_BIN_DIR}/matlab"
                if [[ -f "${MATLAB_EXE_PATH}" ]]; then
                    #MATLAB_ROOT_DIR_OPTION="-DMATLAB_ROOT_DIR=${MATLAB_ROOT_DIR}"
                    MATLAB_EXE_PATH="${MATLAB_EXE_PATH}"
                    MATLAB_BIN_DIR="${MATLAB_BIN_DIR}"
                    echo >&2 "-- ${BUILD_NAME}MATLAB - MATLAB ${MATLAB_VERSION} installation detected at: ${MATLAB_EXE_PATH}"
                    echo >&2
                    break 2
                fi
            done
            # if [ -z ${MATLAB_EXE_PATH+x} ]; then break; fi
        done
    fi

    if [ -z ${MATLAB_EXE_PATH+x} ]; then

        EXAMPLE_MATLAB_VERSION="R2020a"
        if [ "${isMacOS}" = "true" ]; then
            EXAMPLE_MATLAB_ROOT_DIR="/Applications/MATLAB_${EXAMPLE_MATLAB_VERSION}.app"
        else
            EXAMPLE_MATLAB_ROOT_DIR="/usr/local/MATLAB/${EXAMPLE_MATLAB_VERSION}"
        fi
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: Exhausted all possible search paths for a MATLAB installation, but failed to find MATLAB."
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: The ParaMonte MATLAB kernel will not be functional without building the required DLL libraries."
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: Please add MATLAB to your environmental variable PATH and rerun the install script."
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: For example, in your current terminal, try:"
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: "
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}:     export PATH=\"PATH_TO_MATLAB_BIN_DIR:$PATH\""
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: "
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: where PATH_TO_MATLAB_BIN_DIR must be replaced with path to the bin folder of the current"
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: installation of MATLAB on your system. Typical MATLAB bin installation path on a 64-bit ${OSNAME}"
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: Operating Systems is a string like the following:"
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: "
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}:     \"${EXAMPLE_MATLAB_ROOT_DIR}/bin\""
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: "
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: where ${EXAMPLE_MATLAB_VERSION} in the path points to the MATLAB ${EXAMPLE_MATLAB_VERSION} version installation on the system. You can also "
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: find the installation location of MATLAB by typing the following command in your MATLAB session:"
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: "
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}:     matlabroot"
        echo >&2
        if [ "${YES_TO_ALL_DISABLED}" = "true" ]; then
            answerNotGiven=true
            while [ "${answerNotGiven}" = "true" ]; do
                read -p "-- ${BUILD_NAME}MATLAB - Do you wish to continue with the rest of the installation process (y/n)? " answer
                if [[ $answer == [yY] || $answer == [yY][eE][sS] ]]; then
                    answer=y
                    answerNotGiven=false
                fi
                if [[ $answer == [nN] || $answer == [nN][oO] ]]; then
                    answer=n
                    answerNotGiven=false
                fi
                if [ "${answerNotGiven}" = "true" ]; then
                    echo >&2 "-- ${BUILD_NAME}MATLAB - please enter either y or n"
                fi
            done
        else
            echo >&2 "-- ${BUILD_NAME}MATLAB - Do you wish to continue with the rest of the installation process (y/n)? y"
            answer=y
        fi
        if [ "${answer}" = "y" ]; then
            echo >&2
            echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: skipping the ParaMonte MATLAB dynamic library build..."
            echo >&2
        else
            echo >&2
            echo >&2 "-- ${BUILD_NAME}MATLAB - exiting the ParaMonte MATLAB library build."
            echo >&2
            exit 1
        fi
    fi

fi

####################################################################################################################################
# build ParaMonte
####################################################################################################################################

# cat <<EOF >${ParaMonte_BLD_DIR}/build.sh
# cmake ${ParaMonte_ROOT_DIR} \
# -DPMCS=${PMCS} \
# -DMPI_ENABLED=${MPI_ENABLED} \
# -DCAFTYPE=${CAFTYPE} \
# -DBTYPE=${BTYPE} \
# -DLTYPE=${LTYPE} \
# -DHEAP_ARRAY_ENABLED=${HEAP_ARRAY_ENABLED} \
# -DCFI_ENABLED=${CFI_ENABLED} \
# -DOMP_ENABLED=${OMP_ENABLED}
# EOF
#
# chmod +x ${ParaMonte_BLD_DIR}/build.sh
# exec ${ParaMonte_BLD_DIR}/build.sh

#unset FC_OPTION
#unset MPIEXEC_OPTION

if [ -z ${Fortran_COMPILER_PATH+x} ]; then
    FC_OPTION=""
else
    FC_OPTION="-DCMAKE_Fortran_COMPILER=${Fortran_COMPILER_PATH}"
    export Fortran_COMPILER_PATH
fi
if [ -z ${MPIEXEC_PATH+x} ]; then
    MPIEXEC_OPTION=""
else
    MPIEXEC_OPTION="-DMPIEXEC_EXECUTABLE=${MPIEXEC_PATH}"
fi

echo >&2
echo >&2 "-- ${BUILD_NAME} - CMAKE Fortran compiler option: ${FC_OPTION}"
echo >&2 "-- ${BUILD_NAME} - CMAKE mpiexec option: ${MPIEXEC_OPTION}"
echo >&2

if [ -f "${CAF_LOCAL_INSTALLATION_SETUP_FILE}" ] && ([ "${gnuInstallEnabled}" = "true" ] || [ "${mpiInstallEnabled}" = "true" ] || [ "${cafInstallEnabled}" = "true" ]); then
    CAF_LOCAL_INSTALLATION_SETUP_FILE_CMD="source ${CAF_LOCAL_INSTALLATION_SETUP_FILE}"
else
    CAF_LOCAL_INSTALLATION_SETUP_FILE_CMD=""
fi
if [ "${isMacOS}" = "true" ]; then CAF_LOCAL_INSTALLATION_SETUP_FILE_CMD=""; fi

####################################################################################################################################
#### generate the current paramonte library name
####################################################################################################################################

if [ "${INTERFACE_LANGUAGE}" = "c++" ] || [ "${INTERFACE_LANGUAGE}" = "C++" ]; then
    LANG_ABBR="cpp"
else
    LANG_ABBR="${INTERFACE_LANGUAGE}"
fi

PMLIB_NAME="paramonte_${LANG_ABBR}_${PLATFORM}_${ARCHITECTURE}_${PMCS}_${BTYPE}_${LTYPE}_${MEMORY_ALLOCATION}"

if [ "${CAF_ENABLED}" = "true" ]; then
    PMLIB_NAME="${PMLIB_NAME}_caf${CAFTYPE}"
fi

if [ "${MPI_ENABLED}" = "true" ]; then
    PMLIB_NAME="${PMLIB_NAME}_${MPILIB_NAME}"
fi

if [ "${OMP_ENABLED}" = "true" ]; then
    PMLIB_NAME="${PMLIB_NAME}_omp"
fi

echo >&2
echo >&2 "-- ${BUILD_NAME} - NOTE: The current ParaMonte library build will tagged with: ${PMLIB_NAME}"
echo >&2

####################################################################################################################################
#### call cmake
####################################################################################################################################

if [ "${CODECOV_ENABLED}" = "true" ]; then
    CODECOV_ENABLED_FLAG="-DCODECOV_ENABLED=${CODECOV_ENABLED}"
    ParaMonteTest_RUN_ENABLED=true
    if [ -z ${TTYPE+x} ]; then
        BASIC_TEST_ENABLED_FLAG="-DBASIC_TEST_ENABLED=true"
        SAMPLER_TEST_ENABLED_FLAG="-DSAMPLER_TEST_ENABLED=true"
        export BASIC_TEST_ENABLED_FLAG
        export SAMPLER_TEST_ENABLED_FLAG
    elif [[ ${TTYPE} == [nN][oO][nN][eE] ]]; then
        echo >&2
        echo >&2 "-- ${BUILD_NAME} - FATAL: Testing must be activated to generate code coverage."
        echo >&2 "-- ${BUILD_NAME} - FATAL: Drop the \"-t none\" from the script options and rerun the build/install script."
        echo >&2
        usage
        exit 1
    fi
else
    CODECOV_ENABLED_FLAG=""
fi

if [ "${DRYRUN_ENABLED}" != "true" ]; then

    # determine whether the current system is a Windows Subsystem for Linux (WSL).
    # The ParaMonte tests that contain internal procedure calls from outside the parent procedure
    # fail with static builds on WSL for code coverage purposes.
    # This GFortran bug does not exist when the program is compiled as a shared library.
    # However, the same bug happens also with shared library files at the time of code coverage generation.
    # The only resolution left is to fence the internal function calls in the code when the code is being
    # built on WSL. This is done by enabling the preprocessor flag OS_IS_WSL.

    if [ "${isWSL}" = "true" ]; then
        OS_IS_WSL_FLAG="-DOS_IS_WSL=${isWSL}"
    else
        OS_IS_WSL_FLAG=""
    fi

    (cd ${ParaMonte_BLD_DIR} && \
    ${CAF_LOCAL_INSTALLATION_SETUP_FILE_CMD} && \
    cmake \
    --verbose=1 \
    "${FC_OPTION}" \
    "${MPIEXEC_OPTION}" \
    "${MATLAB_ROOT_DIR_OPTION}" \
    -DPMLIB_NAME=${PMLIB_NAME} \
    -DINTERFACE_LANGUAGE=${INTERFACE_LANGUAGE} \
    -DPMCS=${PMCS} \
    -DMPI_ENABLED=${MPI_ENABLED} \
    -DCAFTYPE=${CAFTYPE} \
    -DBTYPE=${BTYPE} \
    -DLTYPE=${LTYPE} \
    -DHEAP_ARRAY_ENABLED=${HEAP_ARRAY_ENABLED} \
    -DCFI_ENABLED=${CFI_ENABLED} \
    -DOMP_ENABLED=${OMP_ENABLED} \
    ${SAMPLER_TEST_ENABLED_FLAG} \
    ${BASIC_TEST_ENABLED_FLAG} \
    ${CODECOV_ENABLED_FLAG} \
    ${OS_IS_WSL_FLAG} \
    ${ParaMonte_ROOT_DIR} \
    )
    verify $? "configuration with cmake"

    (cd ${ParaMonte_BLD_DIR} && make)
    verify $? "build with make"

    (cd ${ParaMonte_BLD_DIR} && make install)
    verify $? "installation"

fi

####################################################################################################################################
#### Test the ParaMonte library if requested
####################################################################################################################################

LD_LIBRARY_PATH=${ParaMonte_BLD_DIR}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH

# set test object/module/lib files directories

ParaMonteTest_BLD_DIR="${ParaMonte_BLD_DIR}/test"
ParaMonteTest_OBJ_DIR="${ParaMonteTest_BLD_DIR}/obj"
ParaMonteTest_BIN_DIR="${ParaMonteTest_BLD_DIR}/bin"
ParaMonteTest_SRC_INPUT_DIR="${ParaMonteTest_SRC_DIR}/input"
ParaMonteTest_BIN_INPUT_DIR="${ParaMonteTest_BIN_DIR}/input"

# When generating code coverage, the tests MUST be performed.

if [ "${ParaMonteTest_RUN_ENABLED}" = "true" ]; then

    # first ensure the test input files exist in the test runtime bin folder

    if ! [ -d "${ParaMonteTest_BIN_INPUT_DIR}" ]; then
        mkdir "${ParaMonteTest_BIN_INPUT_DIR}"
    fi

    cp "${ParaMonteTest_SRC_INPUT_DIR}"/* "${ParaMonteTest_BIN_INPUT_DIR}"/

    # run the tests

    if [ "${MPI_ENABLED}" = "true" ]; then
        # first attempt to find an installation of the mpiexec on the system
        if [ -z ${MPIEXEC_PATH+x} ]; then
            MPIEXEC_PATH_RESET_ENABLED="true"
            MPIEXEC_PATH=$(command -v mpiexec)
            echo >&2
            echo >&2 "-- ${BUILD_NAME} - Inferred mpiexec path: $(command -v mpiexec)"
            echo >&2
        else
            MPIEXEC_PATH_RESET_ENABLED="false"
        fi
        if [ -f "${MPIEXEC_PATH}" ]; then
            echo >&2 "-- ${BUILD_NAME} - running command: ${MPIEXEC_PATH} -n ${FOR_COARRAY_NUM_IMAGES} ./testParaMonte"
            (cd ${ParaMonte_BLD_DIR}/test/bin && "${MPIEXEC_PATH}" -n ${FOR_COARRAY_NUM_IMAGES} ./testParaMonte)
            verify $? "test run"
            if [ "${MPIEXEC_PATH_RESET_ENABLED}" = "true" ]; then unset MPIEXEC_PATH; fi
        else
            echo >&2
            echo >&2 "-- ${BUILD_NAME} - ${warning}: No tests of the ParaMonte library will be performed."
            echo >&2 "-- ${BUILD_NAME} - ${warning}: The mpiexec executable could not be found on your system."
            echo >&2 "-- ${BUILD_NAME} - ${warning}: If you do not have an MPI library installed on your system,"
            echo >&2 "-- ${BUILD_NAME} - ${warning}: ParaMonte may be able to install one for you."
            echo >&2
        fi
    else
        if [ "${PMCS}" = "gnu" ]; then
            if [[ -f "${SETUP_FILE_PATH}" ]]; then
                cp ${SETUP_FILE_PATH} ${ParaMonte_BLD_DIR}/test/bin/
            fi
            if [ "${CAF_ENABLED}" = "true" ]; then
                (cd ${ParaMonte_BLD_DIR}/test/bin && source ${SETUP_FILE_PATH} && cafrun -np ${FOR_COARRAY_NUM_IMAGES} ./testParaMonte)
                verify $? "test run"
            else
                (cd ${ParaMonte_BLD_DIR}/test/bin && ./testParaMonte)
                verify $? "test run"
            fi
        fi
        if [ "${PMCS}" = "intel" ]; then
            if [ "${CAF_ENABLED}" = "true" ]; then
                (export FOR_COARRAY_NUM_IMAGES && cd ${ParaMonte_BLD_DIR}/test/bin && ./testParaMonte)
                verify $? "test run"
            else
                (cd ${ParaMonte_BLD_DIR}/test/bin && ./testParaMonte)
                verify $? "test run"
            fi
        fi
    fi
    #if [ "${CODECOV_ENABLED}" = "true" ]; then # stop on test failure when code coverage analysis is being performed.
    #    verify $? "test run"
    #fi
else
    echo >&2 "skipping ParaMonte library test run..."
fi

####################################################################################################################################
# get ParaMonte lib name
####################################################################################################################################

if [ "${LTYPE}" = "static" ]; then
    PMLIB_EXT=".a"
else
    if [ "${PLATFORM}" = "linux" ]; then PMLIB_EXT=".so"; fi
    if [ "${PLATFORM}" = "darwin" ]; then PMLIB_EXT=".dylib"; fi
    if [ "${PLATFORM}" = "mingw" ] || [ "${PLATFORM}" = "cygwin" ]; then PMLIB_EXT=".dll"; fi
fi
PMLIB_FULL_PATH="$(ls ${ParaMonte_LIB_DIR}/libparamonte_*_${BTYPE}_${LTYPE}_${MEMORY_ALLOCATION}*${PMLIB_EXT} | sort -V | tail -n1)"
PMLIB_FULL_NAME=${PMLIB_FULL_PATH##*/}
PMLIB_BASE_NAME=${PMLIB_FULL_NAME%.*}
export PMLIB_FULL_PATH
export PMLIB_FULL_NAME
export PMLIB_BASE_NAME

echo >&2
echo >&2 "-- ${BUILD_NAME} - ParaMonte installed-library full path: ${PMLIB_FULL_PATH}"
echo >&2 "-- ${BUILD_NAME} - ParaMonte installed-library full name: ${PMLIB_FULL_NAME}"
echo >&2 "-- ${BUILD_NAME} - ParaMonte installed-library base name: ${PMLIB_BASE_NAME}"
echo >&2

####################################################################################################################################
# build ParaMonte MATLAB test
####################################################################################################################################

####################################################################################################################################
# check MATLAB's existence
####################################################################################################################################

# it is imperative to nullify these MATLAB variables for all language builds at all times

unset MATLAB_ROOT_DIR
unset MATLAB_EXE_PATH
unset MATLAB_BIN_DIR
unset MATLAB_LIB_DIR
      MATLAB_INC_DIR="."
unset MATLAB_LIBMX_FILE
unset MATLAB_LIBMEX_FILE
unset MATLAB_LIBMAT_FILE
unset MATLAB_VERSION_FILE
# set "MATLAB_INC_DIR_FLAG="

if [ "${INTERFACE_LANGUAGE}" = "matlab" ] && [ "${LTYPE}" = "dynamic" ] && [ "${CFI_ENABLED}" = "true" ]; then

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #: check MATLAB's existence
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    echo >&2
    echo >&2 "-- ${BUILD_NAME}MATLAB - searching for a MATLAB installation on your system..."

    unset MATLAB_ROOT_DIR
    unset MATLAB_EXE_PATH
    unset MATLAB_BIN_DIR

    if command -v matlab >/dev/null 2>&1; then
        MATLAB_EXE_PATH=$(command -v matlab)
        MATLAB_BIN_DIR=$(dirname $MATLAB_EXE_PATH)
        MATLAB_ROOT_DIR=$(dirname $MATLAB_BIN_DIR)
        # echo >&2 "-- ${BUILD_NAME} - MATLAB detected at: ${MATLAB_EXE_PATH}"
    else
        echo >&2 "-- ${BUILD_NAME} - MATLAB could not be found in among the search paths."
        echo >&2 "-- ${BUILD_NAME} - searching for MATLAB in the default installation directories..."
    fi

    if [ -z ${MATLAB_EXE_PATH+x} ]; then
        if [ "${isMacOS}" = "true" ]; then
            INSTALL_LOC_LIST="/Applications/MATLAB_"
        else
            INSTALL_LOC_LIST="/usr/local/MATLAB"
        fi
        MATLAB_VERSION_LIST="R2025b R2025a R2024b R2024a R2023b R2023a R2022b R2022a R2021b R2021a R2020b R2020a R2019b R2019a R2018b R2018a R2017b R2017a"

        for INSTALL_LOC in $INSTALL_LOC_LIST; do
            for MATLAB_VERSION in $MATLAB_VERSION_LIST; do
                if [ "${isMacOS}" = "true" ]; then
                    MATLAB_ROOT_DIR="${INSTALL_LOC_LIST}${MATLAB_VERSION}.app"
                else
                    MATLAB_ROOT_DIR="${INSTALL_LOC_LIST}/${MATLAB_VERSION}"
                fi
                MATLAB_BIN_DIR="${MATLAB_ROOT_DIR}/bin"
                MATLAB_EXE_PATH="${MATLAB_BIN_DIR}/matlab"
                if [[ -f "${MATLAB_EXE_PATH}" ]]; then
                    MATLAB_ROOT_DIR="${MATLAB_ROOT_DIR}"
                    MATLAB_EXE_PATH="${MATLAB_EXE_PATH}"
                    MATLAB_BIN_DIR="${MATLAB_BIN_DIR}"
                    echo >&2 "-- ${BUILD_NAME}MATLAB - MATLAB ${MATLAB_VERSION} installation detected at: ${MATLAB_EXE_PATH}"
                    echo >&2
                    break 2
                fi
            done
            # if [ -z ${MATLAB_EXE_PATH+x} ]; then break; fi
        done
    fi

    if [ -z ${MATLAB_EXE_PATH+x} ]; then

        EXAMPLE_MATLAB_VERSION="R2020a"
        if [ "${isMacOS}" = "true" ]; then
            EXAMPLE_MATLAB_ROOT_DIR="/Applications/MATLAB_${EXAMPLE_MATLAB_VERSION}.app"
        else
            EXAMPLE_MATLAB_ROOT_DIR="/usr/local/MATLAB/${EXAMPLE_MATLAB_VERSION}"
        fi
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: Exhausted all possible search paths for a MATLAB installation, but failed to find MATLAB."
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: The ParaMonte MATLAB kernel will not be functional without building the required DLL libraries."
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: Please add MATLAB to your environmental variable PATH and rerun the install script."
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: For example, in your current terminal, try:"
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: "
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}:     export PATH=\"PATH_TO_MATLAB_BIN_DIR:$PATH\""
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: "
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: where PATH_TO_MATLAB_BIN_DIR must be replaced with path to the bin folder of the current"
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: installation of MATLAB on your system. Typical MATLAB bin installation path on a 64-bit ${OSNAME}"
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: Operating Systems is a string like the following:"
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: "
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}:     \"${EXAMPLE_MATLAB_ROOT_DIR}/bin\""
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: "
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: where ${EXAMPLE_MATLAB_VERSION} in the path points to the MATLAB ${EXAMPLE_MATLAB_VERSION} version installation on the system. You can also "
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: find the installation location of MATLAB by typing the following command in your MATLAB session:"
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: "
        echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}:     matlabroot"
        echo >&2
        if [ "${YES_TO_ALL_DISABLED}" = "true" ]; then
            answerNotGiven=true
            while [ "${answerNotGiven}" = "true" ]; do
                read -p "-- ${BUILD_NAME}MATLAB - Do you wish to continue with the rest of the installation process (y/n)? " answer
                if [[ $answer == [yY] || $answer == [yY][eE][sS] ]]; then
                    answer=y
                    answerNotGiven=false
                fi
                if [[ $answer == [nN] || $answer == [nN][oO] ]]; then
                    answer=n
                    answerNotGiven=false
                fi
                if [ "${answerNotGiven}" = "true" ]; then
                    echo >&2 "-- ${BUILD_NAME}MATLAB - please enter either y or n"
                fi
            done
        else
            echo >&2 "-- ${BUILD_NAME}MATLAB - Do you wish to continue with the rest of the installation process (y/n)? y"
            answer=y
        fi
        if [ "${answer}" = "y" ]; then
            echo >&2
            echo >&2 "-- ${BUILD_NAME}MATLAB - ${warning}: skipping the ParaMonte MATLAB dynamic library build..."
            echo >&2
        else
            echo >&2
            echo >&2 "-- ${BUILD_NAME}MATLAB - exiting the ParaMonte MATLAB library build."
            echo >&2
            exit 1
        fi

    else

        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        #: Build ParaMonte MATLAB
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        echo >&2
        echo >&2 "-- ${BUILD_NAME}MATLAB - MATLAB installation detected at: ${MATLAB_EXE_PATH}"
        echo >&2

        ParaMonteMATLAB_BLD_LIB_DIR="${ParaMonte_BLD_DIR}/lib"
        if [ -d "${MATLAB_BIN_DIR}" ] && [ "${DRYRUN_ENABLED}" != "true" ]; then
            # cd "${ParaMonteMATLAB_BLD_LIB_DIR}" && fname=$(find -type f -name 'libparamonte_*');
            # PMLIB_NAME_EXT=${fname:2} # removing first two characters './'
            # PMLIB_NAME=$(basename -- "$PMLIB_NAME_EXT")
            # PMLIB_EXT="${PMLIB_NAME##*.}"
            # PMLIB_NAME="${PMLIB_NAME%.*}"
            # PMLIB_MATLAB_NAME="${PMLIB_NAME/_c_/_}"
            PMLIB_MATLAB_NAME="${PMLIB_BASE_NAME/_matlab_/_}"
            # MATLAB_BUILD_FLAGS="-DDLL_ENABLED "
            # if ${BTYPE}==debug   set "MATLAB_BUILD_FLAGS=!MATLAB_BUILD_FLAGS!/Od /Z7"
            # if ${BTYPE}==testing set "MATLAB_BUILD_FLAGS=!MATLAB_BUILD_FLAGS!/O2"
            # if ${BTYPE}==release set "MATLAB_BUILD_FLAGS=!MATLAB_BUILD_FLAGS!/Od"
            # if ${BTYPE}==debug   set "MATLAB_BUILD_FLAGS=!MATLAB_BUILD_FLAGS!!INTEL_CPP_DEBUG_FLAGS!"
            # if ${BTYPE}==testing set "MATLAB_BUILD_FLAGS=!MATLAB_BUILD_FLAGS!!INTEL_CPP_TESTING_FLAGS!"
            # if ${BTYPE}==release set "MATLAB_BUILD_FLAGS=!MATLAB_BUILD_FLAGS!!INTEL_CPP_RELEASE_FLAGS!"
            MEX_FLAGS="-v -nojvm"
            CFLAGS="COMPFLAGS='-fPIC -shared -Wl,-rpath=.'"
            LINKFLAGS="LINKFLAGS='-fPIC -shared -Wl,-rpath=.'"
            if [ "${BTYPE}" = "debug" ]; then MEX_FLAGS="${MEX_FLAGS} -g"; fi
            if [ "${BTYPE}" = "release" ]; then MEX_FLAGS="${MEX_FLAGS} -O"; fi
            echo >&2 "-- ${BUILD_NAME}MATLAB - generating the ParaMonte MATLAB dynamic library: ${ParaMonteMATLAB_BLD_LIB_DIR}${PMLIB_MATLAB_NAME}"
            echo >&2 "-- ${BUILD_NAME}MATLAB - compiler options: ${MATLAB_BUILD_FLAGS}"
            echo >&2 "-- ${BUILD_NAME}MATLAB - compiler command: ${MATLAB_BIN_DIR}/mex ${MEX_FLAGS} ${CFLAGS} ${LINKFLAGS} ${ParaMonteKernel_SRC_DIR}/paramonte.m.c ${PMLIB_FULL_PATH} -output ${PMLIB_MATLAB_NAME}"
            # CC=icl COMPFLAGS="${MATLAB_BUILD_FLAGS}"
            cd "${ParaMonteMATLAB_BLD_LIB_DIR}"
            "${MATLAB_BIN_DIR}/mex" ${MEX_FLAGS} "${CFLAGS}" "${LINKFLAGS}" "${ParaMonteKernel_SRC_DIR}/paramonte.m.c" ${PMLIB_FULL_PATH} -output ${PMLIB_MATLAB_NAME}
            if [ $? -eq 0 ]; then
                echo >&2 "-- ${BUILD_NAME}MATLAB - ${BoldGreen}The ParaMonte MATLAB dynamic library build appears to have succeeded.${ColorReset}"
            else
                echo >&2
                echo >&2 "-- ${BUILD_NAME}MATLAB - ${BoldRed}Fatal Error${ColorReset}: The ParaMonte MATLAB library build failed."
                echo >&2 "-- ${BUILD_NAME}MATLAB - Please make sure you have the following components installed"
                echo >&2 "-- ${BUILD_NAME}MATLAB - on your system before rerunning the installation script:"
                echo >&2 "-- ${BUILD_NAME}MATLAB - "
                echo >&2 "-- ${BUILD_NAME}MATLAB -     -- MATLAB, including MATLAB compilers."
                echo >&2 "-- ${BUILD_NAME}MATLAB -     -- The GNU compiler collection 7 or newer while being compatible with MATLAB."
                echo >&2 "-- ${BUILD_NAME}MATLAB - "
                echo >&2 "-- ${BUILD_NAME}MATLAB - Note that MATLAB is compatible only with the GNU compiler collections on ${OSNAME} "
                echo >&2 "-- ${BUILD_NAME}MATLAB - and not with the Intel compilers. Once you are sure of the existence of these components in your "
                echo >&2 "-- ${BUILD_NAME}MATLAB - ${OSNAME} command line environment, run the following command:"
                echo >&2 "-- ${BUILD_NAME}MATLAB - "
                echo >&2 "-- ${BUILD_NAME}MATLAB -     \"${MATLAB_BIN_DIR}/mex\" -setup C"
                echo >&2 "-- ${BUILD_NAME}MATLAB - "
                echo >&2 "-- ${BUILD_NAME}MATLAB - You should then be able to see a message like the following, "
                echo >&2 "-- ${BUILD_NAME}MATLAB - "
                echo >&2 "-- ${BUILD_NAME}MATLAB -     MEX configured to use 'gcc' for C language compilation."
                echo >&2 "-- ${BUILD_NAME}MATLAB - "
                echo >&2 "-- ${BUILD_NAME}MATLAB - This ensures that your MATLAB is already set to use an appropriate compiler for building applications."
                echo >&2 "-- ${BUILD_NAME}MATLAB - If you cannot identify the cause of the failure, Please report this error at: "
                echo >&2 "-- ${BUILD_NAME}MATLAB - "
                echo >&2 "-- ${BUILD_NAME}MATLAB -     https://github.com/cdslaborg/paramonte/issues"
                echo >&2 "-- ${BUILD_NAME}MATLAB - "
                echo >&2 "-- ${BUILD_NAME}MATLAB - gracefully exiting The ParaMonte build script."
                echo >&2
                exit 1
            fi
            cd "${sourceFileDir}"
        fi

    fi

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #: Build ParaMonte MATLAB test
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    echo >&2
    echo >&2 "-- ${BUILD_NAME}MATLABTest - building ParaMonte MATLAB test..."
    ParaMonteMATLABTest_BLD_DIR=${ParaMonte_BLD_DIR}/test/MATLAB
    export ParaMonteMATLABTest_BLD_DIR

    MATLAB_TEST_FILENAME=testParaMonte_${BTYPE}
    if [ "${CAF_ENABLED}" = "true" ]; then MATLAB_TEST_FILENAME=${MATLAB_TEST_FILENAME}_${CAFTYPE}; fi
    if [ "${MPI_ENABLED}" = "true" ]; then MATLAB_TEST_FILENAME=${MATLAB_TEST_FILENAME}_mpi; fi
    if [ "${OMP_ENABLED}" = "true" ]; then MATLAB_TEST_FILENAME=${MATLAB_TEST_FILENAME}_omp; fi
    MATLAB_TEST_FILENAME=${MATLAB_TEST_FILENAME}.m

    if [ -d "${ParaMonteMATLABTest_BLD_DIR}" ]; then
        echo >&2 "-- ${BUILD_NAME}MATLABTest - ${ParaMonteMATLABTest_BLD_DIR} already exists. skipping..."
    else
        echo >&2 "-- ${BUILD_NAME}MATLABTest - generating MATLAB files directory: ${ParaMonteMATLABTest_BLD_DIR}"
        mkdir "${ParaMonteMATLABTest_BLD_DIR}"
    fi

    # copy necessary ParaMonte MATLAB library files in MATLAB's directory

    echo >&2 "-- ${BUILD_NAME}MATLABTest - copying the ParaMonte library files to the MATLAB directory"
    echo >&2 "-- ${BUILD_NAME}MATLABTest - from: ${ParaMonteInterfaceMATLAB_SRC_DIR}/paramonte"
    echo >&2 "-- ${BUILD_NAME}MATLABTest -   to: ${ParaMonteMATLABTest_BLD_DIR}/paramonte/"
    cp -R "${ParaMonteInterfaceMATLAB_SRC_DIR}"/paramonte "${ParaMonteMATLABTest_BLD_DIR}"/
    echo >&2

    # copy necessary ParaMonte library auxiliary files

    echo >&2 "-- ${BUILD_NAME}MATLABTest - copying the ParaMonte library auxiliary files"
    echo >&2 "-- ${BUILD_NAME}MATLABTest - from: ${ParaMonteInterface_SRC_DIR}/auxil"
    echo >&2 "-- ${BUILD_NAME}MATLABTest -   to: ${ParaMonteMATLABTest_BLD_DIR}/paramonte/"
    cp -R "${ParaMonteInterface_SRC_DIR}/auxil" "${ParaMonteMATLABTest_BLD_DIR}/paramonte/"
    echo >&2

    # copy necessary ParaMonte MATLAB dynamic library files in MATLAB's directory

    if [ -d "${ParaMonteMATLABTest_BLD_DIR}/paramonte/lib" ]; then
        echo >&2 "-- ${BUILD_NAME}MATLABTest - ${ParaMonteMATLABTest_BLD_DIR}/paramonte/lib already exists. skipping..."
    else
        echo >&2 "-- ${BUILD_NAME}MATLABTest - generating MATLAB files directory: ${ParaMonteMATLABTest_BLD_DIR}/paramonte/lib"
        mkdir "${ParaMonteMATLABTest_BLD_DIR}/paramonte/lib"
    fi

    echo >&2 "-- ${BUILD_NAME}MATLABTest - copying the ParaMonte library files..."
    echo >&2 "-- ${BUILD_NAME}MATLABTest - from: ${PMLIB_FULL_PATH}"
    echo >&2 "-- ${BUILD_NAME}MATLABTest -   to: ${ParaMonteMATLABTest_BLD_DIR}/paramonte/lib"
    cp -R "${ParaMonte_LIB_DIR}/"* "${ParaMonteMATLABTest_BLD_DIR}/paramonte/lib"
    echo >&2

    # copy necessary ParaMonte MATLAB library files in MATLAB's directory

    echo >&2 "-- ${BUILD_NAME}MATLABTest - copying the ParaMonte library test files to the MATLAB directory"
    echo >&2 "-- ${BUILD_NAME}MATLABTest - from: ${ParaMonteMATLABTest_SRC_DIR}/${MATLAB_TEST_FILENAME}"
    echo >&2 "-- ${BUILD_NAME}MATLABTest -   to: ${ParaMonteMATLABTest_BLD_DIR}/"
    cp "${ParaMonteMATLABTest_SRC_DIR}"/${MATLAB_TEST_FILENAME} "${ParaMonteMATLABTest_BLD_DIR}"/
    echo >&2

    # copy necessary input files in MATLAB's directory

    echo >&2 "-- ${BUILD_NAME}MATLABTest - copying input files to the MATLAB directory"
    echo >&2 "-- ${BUILD_NAME}MATLABTest - from: ${ParaMonteTest_SRC_DIR}/input"
    echo >&2 "-- ${BUILD_NAME}MATLABTest -   to: ${ParaMonteMATLABTest_BLD_DIR}/input/"
    cp -R "${ParaMonteTest_SRC_DIR}"/input "${ParaMonteMATLABTest_BLD_DIR}"/
    echo >&2

fi

####################################################################################################################################
# build ParaMonte Python test
####################################################################################################################################

if [ "${INTERFACE_LANGUAGE}" = "python" ] && [ "${LTYPE}" = "dynamic" ] && [ "${CFI_ENABLED}" = "true" ]; then

    echo >&2
    echo >&2 "-- ${BUILD_NAME}PythonTest - building ParaMonte Python test..."
    ParaMontePythonTest_BLD_DIR=${ParaMonte_BLD_DIR}/test/Python
    export ParaMontePythonTest_BLD_DIR

    PYTHON_TEST_FILENAME=testParaMonte_${BTYPE}
    if [ "${CAF_ENABLED}" = "true" ]; then PYTHON_TEST_FILENAME=${PYTHON_TEST_FILENAME}_${CAFTYPE}; fi
    if [ "${MPI_ENABLED}" = "true" ]; then PYTHON_TEST_FILENAME=${PYTHON_TEST_FILENAME}_mpi; fi
    if [ "${OMP_ENABLED}" = "true" ]; then PYTHON_TEST_FILENAME=${PYTHON_TEST_FILENAME}_omp; fi
    PYTHON_TEST_FILENAME=${PYTHON_TEST_FILENAME}.py

    if [[ -d "$ParaMontePythonTest_BLD_DIR" ]]
    then
        echo >&2 "-- ${BUILD_NAME}PythonTest - ${ParaMontePythonTest_BLD_DIR} already exists. skipping..."
    else
        echo >&2 "-- ${BUILD_NAME}PythonTest - generating Python files directory: ${ParaMontePythonTest_BLD_DIR}"
        mkdir "${ParaMontePythonTest_BLD_DIR}"
    fi

    # copy necessary ParaMonte Python library files in Python's directory

    echo >&2 "-- ${BUILD_NAME}PythonTest - copying the ParaMonte library files to the Python directory"
    echo >&2 "-- ${BUILD_NAME}PythonTest - from: ${ParaMonteInterfacePython_SRC_DIR}/paramonte"
    echo >&2 "-- ${BUILD_NAME}PythonTest -   to: ${ParaMontePythonTest_BLD_DIR}/paramonte/"
    cp -R "${ParaMonteInterfacePython_SRC_DIR}"/paramonte "${ParaMontePythonTest_BLD_DIR}"/
    echo >&2

    # copy necessary ParaMonte library auxiliary files

    echo >&2 "-- ${BUILD_NAME}PythonTest - copying the ParaMonte library auxiliary files"
    echo >&2 "-- ${BUILD_NAME}PythonTest - from: ${ParaMonteInterface_SRC_DIR}/auxil"
    echo >&2 "-- ${BUILD_NAME}PythonTest -   to: ${ParaMontePythonTest_BLD_DIR}/paramonte/"
    cp -R "${ParaMonteInterface_SRC_DIR}/auxil" "${ParaMontePythonTest_BLD_DIR}/paramonte/"
    echo >&2

    # copy necessary ParaMonte Python DLL files in Python's directory

    echo >&2 "-- ${BUILD_NAME}PythonTest - copying the ParaMonte shared library files to the Python directory"
    echo >&2 "-- ${BUILD_NAME}PythonTest - from: ${ParaMonte_LIB_DIR}/${PMLIB_FULL_NAME}"
    echo >&2 "-- ${BUILD_NAME}PythonTest -   to: ${ParaMontePythonTest_BLD_DIR}/paramonte/"
    cp "${ParaMonte_LIB_DIR}"/* "${ParaMontePythonTest_BLD_DIR}"/paramonte/
    echo >&2

    # copy necessary ParaMonte Python library files in Python's directory

    echo >&2 "-- ${BUILD_NAME}PythonTest - copying the ParaMonte library test files to the Python directory"
    echo >&2 "-- ${BUILD_NAME}PythonTest - from: ${ParaMontePythonTest_SRC_DIR}/${PYTHON_TEST_FILENAME}"
    echo >&2 "-- ${BUILD_NAME}PythonTest -   to: ${ParaMontePythonTest_BLD_DIR}/"
    cp "${ParaMontePythonTest_SRC_DIR}"/${PYTHON_TEST_FILENAME} "${ParaMontePythonTest_BLD_DIR}"/
    echo >&2

    # copy necessary input files in Python's directory

    echo >&2 "-- ${BUILD_NAME}PythonTest - copying input files to the Python directory"
    echo >&2 "-- ${BUILD_NAME}PythonTest - from: ${ParaMonteTest_SRC_DIR}/input"
    echo >&2 "-- ${BUILD_NAME}PythonTest -   to: ${ParaMontePythonTest_BLD_DIR}/input/"
    cp -R "${ParaMonteTest_SRC_DIR}"/input "${ParaMontePythonTest_BLD_DIR}"/
    echo >&2

fi

####################################################################################################################################
#### Generate *.gcda *.gcno codecov files for the kernel test source files
####################################################################################################################################

if [ "${CODECOV_ENABLED}" = "true" ]; then

    # This version extraction relies on "version " appearing before the lcov version number.
    # This must be adjusted to any potential future changes in the LCOV version string.
    if command -v lcov >/dev/null 2>&1; then
        lcovVersion=$(lcov -v | grep -Po '(?<=version )[^;]+')
    else
        unset lcovVersion
    fi

    unset htmlSubDir
    if [ "${MPI_ENABLED}" = "true" ]; then
        parallelismText="MPI Parallel"
        htmlSubDir="mpi"
    elif [ "${CAF_ENABLED}" = "true" ]; then
        parallelismText="Coarray Parallel"
        htmlSubDir="caf"
    else
        parallelismText="Serial"
        htmlSubDir="serial"
    fi
    htmlDir="${ParaMonte_ROOT_DIR}/codecov/kernel/${htmlSubDir}"
    #htmlDir="${ParaMonte_ROOT_DIR}/codecov/kernel/${ParaMonteVersion}/${PMLIB_BASE_NAME}"
    htmlTitleCodeCov="ParaMonte ${ParaMonteVersion} :: ${parallelismText} Kernel - Code Coverage Report"

    if [[ ${PMCS} == [gG][nN][uU] ]]; then

        if command -v gcov >/dev/null 2>&1; then

            gcovPath="$(command -v gcov)"
            echo >&2
            echo >&2 "-- ${BUILD_NAME}CodeCoverage - GNU gcov detected at: ${gcovPath}"
            echo >&2 "-- ${BUILD_NAME}CodeCoverage - Invoking gcov to generate coverage report..."
            echo >&2
            if [ "${HEAP_ARRAY_ENABLED}" = "true" ]; then
                MEMORY_ALLOCATION="heap"
            else
                MEMORY_ALLOCATION="stack"
            fi

            gcovKernelDataDir=$(find "${ParaMonte_OBJ_DIR}" -name ParaMonte_mod*.o)
            gcovKernelDataDir=$(dirname "${gcovKernelDataDir}")

            if [ -d "${gcovKernelDataDir}" ]; then

                gcovKernelDir="${ParaMonte_BLD_DIR}"/gcov
                if ! [ -d "${gcovKernelDir}" ]; then
                    mkdir -p "${gcovKernelDir}"
                fi
                cd "${gcovKernelDir}"

                ####################################################################################################################
                #### GCOV KERNEL SOURCE: Generate *.gcda *.gcno codecov files for the kernel source files
                ####################################################################################################################

                for srcFileName in "${ParaMonteKernel_SRC_DIR}"/*.f90; do
                    if ! [[ "${srcFileName}" =~ .*"ifport.f90".* ]]; then

                        srcFileNameBase=$(basename -- "${srcFileName}")
                        #if ! [[ "${srcFileName}" =~ .*".inc.f90".* ]]; then
                        #objFilePath=$(find "${gcovKernelDataDir}" -name ${srcFileNameBase}.o)
                        # The following assumes that cmake names the object files with full source file name (including file extension).
                        objFilePath="${gcovKernelDataDir}/${srcFileNameBase}.o"
                        echo >&2 gcov "${ParaMonteKernel_SRC_DIR}/${srcFileName}" -o "${objFilePath}"

                        gcov "${ParaMonteKernel_SRC_DIR}/${srcFileName}" -o "${objFilePath}" \
                        || {
                            echo >&2
                            echo >&2 "-- ${BUILD_NAME}CodeCoverage - Fatal Error: Code Coverage analysis via GNU gcov tool failed."
                            echo >&2
                            exit 1
                        }
                        #fi

                    fi
                done

                ####################################################################################################################
                #### GCOV TEST SOURCE: Generate *.gcda *.gcno codecov files for the kernel test source files
                ####################################################################################################################

                #### first attempt to infer the cmake object files' directory

                gcovKernelTestDataDir=$(find "${ParaMonteTest_OBJ_DIR}" -name Test_mod*.o)
                gcovKernelTestDataDir=$(dirname "${gcovKernelTestDataDir}")

                if [ -d "${gcovKernelTestDataDir}" ]; then

                    gcovKernelTestDir="${ParaMonteTest_BLD_DIR}"/gcov
                    rm -rf "${gcovKernelTestDir}"
                    if ! [ -d "${gcovKernelTestDir}" ]; then
                        mkdir -p "${gcovKernelTestDir}"
                    fi
                    cd "${gcovKernelTestDir}"

                    for srcFileName in "${ParaMonteTest_SRC_DIR}"/*.f90; do
                        if ! [[ "${srcFileName}" =~ .*"Test_mod.f90".* ]]; then

                            srcFileNameBase=$(basename -- "${srcFileName}")
                            objFilePath="${gcovKernelTestDataDir}/${srcFileNameBase}.o"
                            echo >&2 gcov "${ParaMonteTest_SRC_DIR}/${srcFileName}" -o "${objFilePath}"

                            gcov "${ParaMonteTest_SRC_DIR}/${srcFileName}" -o "${objFilePath}" \
                            || {
                                echo >&2
                                echo >&2 "-- ${BUILD_NAME}CodeCoverage - Fatal Error: The ParaMonte Code Coverage analysis of the test source files via GNU gcov tool failed."
                                echo >&2 "-- ${BUILD_NAME}CodeCoverage - Fatal Error: The ParaMonte test source file: ${ParaMonteTest_SRC_DIR}/${srcFileName}"
                                echo >&2 "-- ${BUILD_NAME}CodeCoverage - skipping..."
                                echo >&2
                            }
                            #fi
                        fi
                    done

                else

                    echo >&2
                    echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${warning}: the directory for the *.gcda *.gcno codecov data files does not exist."
                    echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${warning}: the expected directory path: ${gcovKernelTestDataDir}"
                    echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${warning}: skipping code coverage report generation for the test files..."
                    echo >&2

                fi

                ####################################################################################################################
                # LCOV: generate lcov summary report file
                ####################################################################################################################

                if command -v lcov >/dev/null 2>&1; then

                    #### generate the kernel code coverage report file

                    lcovKernelDir="${ParaMonte_BLD_DIR}"/lcov
                    rm -rf "${lcovKernelDir}"
                    mkdir -p "${lcovKernelDir}"
                    cd "${lcovKernelDir}"

                    lcovOutputKernelFilePath="${lcovKernelDir}/paramonte.kernel.coverage.info"

                    unset branchCoverageFlag
                    # Add the following flag to lcov to enable branch coverage:
                    # branchCoverageFlag="--rc lcov_branch_coverage=1"
                    # "${branchCoverageFlag}"

                    lcov --capture \
                    --directory "${gcovKernelDataDir}" \
                    --output-file "${lcovOutputKernelFilePath}" \
                    || {
                        echo >&2
                        echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${pmwarn} Code Coverage report generation via LCOV tool failed."
                        echo >&2
                        #exit 1
                    }

                    #### generate the kernel code coverage report file

                    lcovOutputCombinedFilePath="${lcovKernelDir}/paramonte.combined.coverage.info"

                    unset lcovOutputTestFilePath
                    if ls "${gcovKernelTestDir}"/*.gcov 1> /dev/null 2>&1; then

                        echo >&2
                        echo >&2 "-- ${BUILD_NAME}CodeCoverage - generating the code coverage report file for the ParaMonte test files..."
                        echo >&2

                        lcovOutputTestFilePath="${lcovKernelDir}/paramonte.test.coverage.info"

                        lcov --capture \
                        --directory "${gcovKernelTestDataDir}" \
                        --output-file "${lcovOutputTestFilePath}" \
                        && {

                            echo >&2
                            echo >&2 "-- ${BUILD_NAME}CodeCoverage - Combining all LCOV code coverage report files as a single final report file..."
                            echo >&2

                            lcov --add-tracefile "${lcovOutputKernelFilePath}" -a "${lcovOutputTestFilePath}" -o "${lcovOutputCombinedFilePath}"

                        } || {

                            echo >&2
                            echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${warning}: Code Coverage report generation for the ParaMonte test source files via lcov tool failed."
                            echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${warning}: skipping..."
                            echo >&2

                        }

                    else

                        echo >&2
                        echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${warning}: Failed to detect the *.gcda *.gcno codecov data files for the ParaMonte test source files."
                        echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${warning}: the expected directory path for the files: ${gcovKernelTestDir}"
                        echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${warning}: The coverage report for the ParaMonte test source file will not be included."
                        echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${warning}: skipping the lcov code coverage report generation for the test files..."
                        echo >&2

                    fi

                    ################################################################################################################
                    # HTML: convert the lcov summary file to the final html report files
                    ################################################################################################################

                    if command -v genhtml >/dev/null 2>&1; then

                        if [ -d "${htmlDir}" ]; then
                            rm -rf "${htmlDir}"
                        else
                            mkdir -p "${htmlDir}"
                        fi

                        if ! [ -f "${lcovOutputCombinedFilePath}" ]; then
                            cp "${lcovOutputKernelFilePath}" "${lcovOutputCombinedFilePath}" || {
                                echo >&2
                                echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${pmwarn} copy action failed:"
                                echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${pmwarn} from: ${lcovOutputKernelFilePath}"
                                echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${pmwarn}   to: ${lcovOutputCombinedFilePath}"
                                #echo >&2 "-- ${BUILD_NAME}CodeCoverage - "
                                #echo >&2 "-- ${BUILD_NAME}CodeCoverage - gracefully exiting The ParaMonte build script."
                                #echo >&2
                                #exit 1
                            }
                        fi

                        genhtml \
                        "${lcovOutputCombinedFilePath}" \
                        --output-directory "${htmlDir}" \
                        --legend \
                        --title "${htmlTitleCodeCov}" \
                        && {

                            echo >&2
                            echo >&2 "-- ${BUILD_NAME}CodeCoverage - The code coverage build files are stored at: ${ParaMonte_BLD_DIR}"
                            echo >&2 "-- ${BUILD_NAME}CodeCoverage - The code coverage report files are stored at: ${htmlDir}"
                            echo >&2

                            # postprocess the html files

                            pmlinkopen='<a href="https:\/\/www.cdslab.org\/paramonte\/" target="_blank">'
                            pmlinklogo='<img alt="The ParaMonte Documentation Website" src="https:\/\/www.cdslab.org\/paramonte\/notes\/api\/kernel\/logo.png"\/>'
                            pmlinkclose='<\/a>'

                            original='<tr><td class="title">LCOV - code coverage report<\/td><\/tr>'
                            modified='<tr><td class="title">'
                            modified+="${pmlinkopen}"
                            modified+="${pmlinklogo}"
                            modified+="${pmlinkclose}"
                            modified+='<\/td><\/tr>'

                            footer='<tr><td class="versionInfo">'
                            footer+='<a href="https:\/\/www.cdslab.org\/paramonte"><b>ParaMonte: Plain Powerful Parallel Monte Carlo Library<\/b><\/a>&nbsp;<br>'
                            footer+='<a href="https:\/\/www.cdslab.org" target="_blank"><b>The Computational Data Science Lab<\/b><\/a><br>'
                            footer+="&copy; Copyright 2012 - $(date +%Y)"
                            footer+='<\/td><\/tr>'

                            shopt -s globstar
                            for htmlFilePath in "${htmlDir}"/**/*.html; do # Whitespace-safe and recursive
                                sed -i "s/${original}/${modified}/g" "${htmlFilePath}"
                                sed -i "/<tr><td class=\"versionInfo\">/c\\${footer}" "${htmlFilePath}"
                            done

                            #scfile="${ParaMonte_ROOT_DIR}/auxil/sc.html"
                            #if [ -f "${scfile}" ]; then
                            #    echo >&2
                            #    echo >&2 "-- ${BUILD_NAME}CodeCoverage - processing the sc file contents..."
                            #    echo >&2
                            #    sccontents=`cat "${ParaMonte_ROOT_DIR}/auxil/sc.html"`
                            #    sed -e '/<\/body>/r${scfile}' "${htmlFilePath}"
                            #else
                            #    echo >&2
                            #    echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${warning} ${ParaMonte_ROOT_DIR}/auxil/sc.html is missing in your clone."
                            #    echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${warning} This is not critical, unless you are a ParaMonte developer and"
                            #    echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${warning} aim to publicly release this code coverage report. To obtain a "
                            #    echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${warning} copy of the file, contact the ParaMonte lead developer at"
                            #    echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${warning} "
                            #    echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${warning} shahmoradi@utexas.edu"
                            #    echo >&2
                            #fi

                        } || {
                            echo >&2
                            echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${pmwarn} Code Coverage report generation via genhtml failed."
                            echo >&2
                            #exit 1
                        }
                        # "${branchCoverageFlag}" \
                        #--title "<a href=\"https://github.com/cdslaborg/paramonte\" target=\"_blank\">ParaMonte::kernel</a> code coverage report" \

                        ## generate test files code coverage
                        #
                        #gcovKernelTestDataDir=$(find "${ParaMonte_OBJ_DIR}" -name Test_mod*.o)
                        #gcovKernelTestDataDir=$(dirname "${gcovKernelTestDataDir}")
                        #
                        #lcovKernelTestDir="${ParaMonte_BLD_DIR}"/test/lcov
                        #mkdir -p "${lcovKernelTestDir}"
                        #cd "${lcovKernelTestDir}"
                        #
                        #lcov --capture --directory "${gcovKernelTestDataDir}" --output-file ./paramonte.coverage.info
                        #
                        #genhtml paramonte.coverage.info --output-directory "${lcovKernelTestDir}/html"

                    else
                        echo >&2
                        echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${pmwarn} Failed to find the GENHTML test coverage summarizer."
                        echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${pmwarn} The genhtml program is required to generate the coverage report."
                        echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${pmwarn} If you believe genhtml is already installed on your system,"
                        echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${pmwarn} please make sure the path its directory is added to the"
                        echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${pmwarn} PATH environmental variable of your terminal."
                        echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${pmwarn} Once added, rerun the ParaMonte code coverage."
                        #echo >&2 "-- ${BUILD_NAME}CodeCoverage - "
                        #echo >&2 "-- ${BUILD_NAME}CodeCoverage - gracefully exiting The ParaMonte build script."
                        #echo >&2
                        #exit 1
                    fi

                else

                    echo >&2
                    echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${pmwarn} Failed to find the LCOV test coverage summarizer."
                    echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${pmwarn} The lcov program is required to generate the coverage report."
                    echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${pmwarn} If you believe lcov is already installed on your system,"
                    echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${pmwarn} please make sure the path its directory is added to the"
                    echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${pmwarn} PATH environmental variable of your terminal."
                    echo >&2 "-- ${BUILD_NAME}CodeCoverage - ${pmwarn} Once added, rerun the ParaMonte code coverage."
                    #echo >&2 "-- ${BUILD_NAME}CodeCoverage - "
                    #echo >&2 "-- ${BUILD_NAME}CodeCoverage - gracefully exiting The ParaMonte build script."
                    #echo >&2
                    #exit 1

                fi

                cd "${ParaMonte_ROOT_DIR}"

            else

                echo >&2
                echo >&2 "-- ${BUILD_NAME}CodeCoverage - Fatal Error: Failed to find the ParaMonte library objects directory."
                echo >&2 "-- ${BUILD_NAME}CodeCoverage - "
                echo >&2 "-- ${BUILD_NAME}CodeCoverage - gracefully exiting The ParaMonte build script."
                echo >&2
                exit 1

            fi

        else

            echo >&2
            echo >&2 "-- ${BUILD_NAME}CodeCoverage - Fatal Error: Failed to find the GNU gcov test coverage program."
            echo >&2 "-- ${BUILD_NAME}CodeCoverage - The gcov program is required to generate the coverage report."
            echo >&2 "-- ${BUILD_NAME}CodeCoverage - If you believe gcov is already installed on your system,"
            echo >&2 "-- ${BUILD_NAME}CodeCoverage - please make sure the path its directory is added to the"
            echo >&2 "-- ${BUILD_NAME}CodeCoverage - PATH environmental variable of your terminal."
            echo >&2 "-- ${BUILD_NAME}CodeCoverage - Once added, rerun the ParaMonte code coverage."
            echo >&2 "-- ${BUILD_NAME}CodeCoverage - "
            echo >&2 "-- ${BUILD_NAME}CodeCoverage - gracefully exiting The ParaMonte build script."
            echo >&2
            exit 1

        fi

    elif [[ ${PMCS} == [iI][nN][tT][eE][lL] ]]; then

            echo >&2
            echo >&2 "-- ${BUILD_NAME}CodeCoverage - Fatal Error: Intel code coverage not supported."
            echo >&2
            exit 1

    fi

    exit

fi

####################################################################################################################################
# build ParaMonte example
####################################################################################################################################

cd "${ParaMonte_ROOT_DIR}"
source "${ParaMonte_ROOT_DIR}/example/buildParaMonteExample.sh"

####################################################################################################################################
# copy ParaMonte binary / library files to the bin directory
####################################################################################################################################

echo >&2
echo >&2 "-- ${BUILD_NAME} - copying the ParaMonte binary/library files to the bin directory..."
ParaMonte_BIN_DIR=${ParaMonte_ROOT_DIR}/bin
export ParaMonte_BIN_DIR
mkdir -p ${ParaMonte_BIN_DIR}


####################################################################################################################################
# mission accomplished
####################################################################################################################################

ParaMonte_BIN_DIR_CURRENT="${ParaMonte_BIN_DIR}"
if [ "${INTERFACE_LANGUAGE}" = "python" ]; then
    ParaMonte_BIN_DIR_CURRENT="${ParaMonte_BIN_DIR_CURRENT}/libparamonte_Python"
elif [ "${INTERFACE_LANGUAGE}" = "matlab" ]; then
    ParaMonte_BIN_DIR_CURRENT="${ParaMonte_BIN_DIR_CURRENT}/libparamonte_MATLAB"
elif [ "${INTERFACE_LANGUAGE}" = "matdram" ]; then
    ParaMonte_BIN_DIR_CURRENT="${ParaMonte_BIN_DIR_CURRENT}/libparamonte_MatDRAM"
else
    ParaMonte_BIN_DIR_CURRENT="${ParaMonte_BIN_DIR_CURRENT}/${PMLIB_BASE_NAME}"
fi

echo >&2
echo >&2 "-- ${BUILD_NAME} - ParaMonte build directory: ${ParaMonte_BLD_DIR}"
echo >&2 "-- ${BUILD_NAME} - ParaMonte install directory: ${ParaMonte_BIN_DIR_CURRENT}"
echo >&2
echo >&2 "-- ${BUILD_NAME} - mission accomplished"
echo >&2

####################################################################################################################################
# clean up variables
####################################################################################################################################

if ! [ -z ${SETUP_FILE_PATH+x} ]; then
    rm "${SETUP_FILE_PATH}"
fi

if [ "${CLEAN}" = "true" ]; then

    for SUITE in $SUITE_LIST
    do
        for LANG in $LANG_LIST
        do
            suiteLangCompilerName="${SUITE}${LANG}CompilerName"
            suiteLangCompilerPath="${SUITE}${LANG}CompilerPath"
            suiteLangCompilerVersion="${SUITE}${LANG}CompilerVersion"
            eval unset ${suiteLangCompilerName}
            eval unset ${suiteLangCompilerPath}
            eval unset ${suiteLangCompilerVersion}
        done
    done

    unset PMCS
    unset BTYPE
    unset LTYPE
    unset CAFTYPE
    unset CFI_ENABLED
    unset MPI_ENABLED
    unset OMP_ENABLED
    unset HEAP_ARRAY_ENABLED
    unset SUITE_LIST
    unset LANG_LIST
    unset ParaMonteTest_RUN_ENABLED
    unset ParaMonteExample_RUN_ENABLED
    unset Fortran_COMPILER_PATH
    unset MPIEXEC_PATH

fi

echo
