#!/bin/bash
#**********************************************************************************************************************************
#**********************************************************************************************************************************
#
#  ParaMonte: plain powerful parallel Monte Carlo library.
#
#  Copyright (C) 2012-present, The Computational Data Science Lab
#
#  This file is part of ParaMonte library. 
#
#  ParaMonte is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, version 3 of the License.
#
#  ParaMonte is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
#
#**********************************************************************************************************************************
#**********************************************************************************************************************************

# NOTE: Do not change the contents of this file unless you know what the consequences are.
# This is the Bash sctipt file that builds objects, dynamic libraries, 
# as well as the test and example binaries of the ParaMonte library on non-Windows systems.
# Upon invocation of this file from a Bash command-line interface, 
# this file will first call the configuration file configParaMonte.bat to read the user's
# requested configuration for building the ParaMonte library.

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: set build type: release, debug, testing :: set library type: static, dynamic
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

FILE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

ParaMonte_ROOT_DIR=${FILE_DIR}
export ParaMonte_ROOT_DIR
#export ParaMonte_ROOT_DIR="${ParaMonte_ROOT_DIR:-${PWD%/}}"

if [[ ! -f "${ParaMonte_ROOT_DIR}/src/ParaMonte/ParaMonte_mod.f90" ]]; then
  echo >&2
  echo >&2 "-- ParaMonte - FATAL: build failed."
  echo >&2 "-- ParaMonte - FATAL: Please run this script inside the top-level ParaMonte library root directory."
  echo >&2 "-- ParaMonte - FATAL: This is the directory which contains this file in the GitHub repository of ParaMonte."
  echo >&2
  exit 1
fi

UNAME_PLATFORM="$(uname -s)"
case "${UNAME_PLATFORM}" in
    Linux*)     PLATFORM=linux;;
    Darwin*)    PLATFORM=mac;;
    CYGWIN*)    PLATFORM=cygwin;;
    MINGW*)     PLATFORM=mingw;;
    *)          PLATFORM="unknown:${UNAME_PLATFORM}"
esac
if [[ "$PLATFORM" =~ .*"unknown".* ]]; then
    echo >&2
    echo >&2 "-- ParaMonte - FATAL: build failed. unrecognized platform - ${PLATFORM}"
    echo >&2 "-- ParaMonte - supported platforms include: Linux, Darwin, CYGWIN, MINGW"
    echo >&2 "-- ParaMonte - ParaMonte build has been only tested on Linux and Darwin platforms."
    echo >&2
    echo >&2 "-- ParaMonte - gracefully exiting."
    echo >&2
else
    export PLATFORM
fi


ARCHITECTURE=$(uname -p)
if [[ "$ARCHITECTURE" =~ .*"64".* ]]; then
    ARCHITECTURE="x64"
else
    ARCHITECTURE=$(uname -m)
    if [[ "$ARCHITECTURE" =~ .*"64".* ]]; then ARCHITECTURE="x64"; fi
fi
export ARCHITECTURE

read -r ParaMonteVersion < .VERSION

#echo "$(cat ./auxil/ParaMonteBanner.txt)"

echo >&2 
echo >&2 "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
echo >&2 "::::                                                                                                                            ::::"
echo >&2 "                                      ParaMonte library version ${ParaMonteVersion} build on ${PLATFORM}                            "
echo >&2 "::::                                                                                                                            ::::"
echo >&2 "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
echo >&2

echo >&2
echo >&2 "-- ParaMonte - current directory: ${FILE_DIR}"
echo >&2 "-- ParaMonte - current system's platform: ${PLATFORM}"
echo >&2 "-- ParaMonte - current system's architecture: ${ARCHITECTURE}"

####################################################################################################################################
# Configure ParaMonte Build
####################################################################################################################################

chmod a+x ./configParaMonte.sh

echo >&2
echo >&2 "-- ParaMonte - configuring ParaMonte Build..."

echo >&2
echo >&2 "-- ParaMonte - default configuration: "

#output=$("./configParaMonte.sh")
#echo "${output}"

. ./configParaMonte.sh

usage()
{
cat << EndOfMessage

    usage:

        buildParaMonte.sh
        -s <compiler suite: intel/gnu>
        -b <build mode: release/testing/debug>
        -l <library type: static/dynamic>
        -c <coarray: none/single/shared/distributed>
        -m <mpi enabled: true/false>
        -i <C-Fortran interface enabled: true/false>
        -e <heap allocation enabled: true/false>
        -t <ParaMonte test run enabled: true/false>
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

        -s | --compiler_suite   : the ParaMonte library build compiler suite: intel, gnu
        -b | --build_type       : the ParaMonte library build type: release, testing, debug
        -l | --lib_type         : the ParaMonte library type: static, dynamic
        -c | --caf_type         : the ParaMonte library Coarray Fortran parallelism type: none, single, shared, distributed
        -m | --mpi_enabled      : the ParaMonte library MPI parallelism enabled?: true, false
        -i | --cfi_enabled      : the ParaMonte library C-Fortran interface enabled? must be true if the library is to be called from non-Fortran languages: true, false
        -e | --heap_enabled     : the ParaMonte library heap array allocation enabled?: true, false
        -t | --test_enabled     : the ParaMonte library test run enabled?: true, false
        -x | --exam_enabled     : the ParaMonte library examples run enabled?: true, false
        -f | --fortran          : path to Fortran compiler. If provided, the ParaMonte library will be built via the specified compiler.
        -M | --mpiexec          : path to mpiexec routine. If provided, it will be used to find the MPI library.
        -F | --fresh            : enables a fresh installation of all of the prerequisites of ParaMonte library. Applicable only to GNU compiler suite.
        -y | --yes-to-all       : if a fresh installation of all of the prerequisites is needed, automatically answer yes to all permission requests.
        -B | --bootstrap        : enables robust bootstrap build when building the required GCC version with an old GCC version. Applicable only to GNU compiler suite.
        -n | --num_images       : the default number of processes (coarray images) on which the ParaMonte examples/tests (if any) will be run: positive integer
        -a | --clean            : clean the environmental variables upon exit, if flag is provided.
        -h | --help             : help with the sctipt usage


EndOfMessage
}

unset PMCS
unset MPIEXEC_PATH
unset GCC_BOOTSTRAP
unset Fortran_COMPILER_PATH
FRESH_INSTALL_ENABLED=false
YES_TO_ALL_DISABLED=true
CLEAN=false

while [ "$1" != "" ]; do
    case $1 in
        -s | --compiler_suite ) shift
                                PMCS=$1
                                ;;
        -b | --build_type )     shift
                                BTYPE=$1
                                ;;
        -l | --lib_type )       shift
                                LTYPE=$1
                                ;;
        -c | --caf_type )       shift
                                CAFTYPE=$1
                                ;;
        -m | --mpi_enabled )    shift
                                MPI_ENABLED=$1
                                ;;
        -i | --cfi_enabled )    shift
                                CFI_ENABLED=$1
                                ;;
        -e | --heap_enabled )   shift
                                HEAP_ARRAY_ENABLED=$1
                                ;;
        -t | --test_enabled )   shift
                                ParaMonteTest_RUN_ENABLED=$1; export ParaMonteTest_RUN_ENABLED
                                ;;
        -x | --exam_enabled )   shift
                                ParaMonteExample_RUN_ENABLED=$1; export ParaMonteExample_RUN_ENABLED
                                ;;
        -f | --fortran )        shift
                                Fortran_COMPILER_PATH=$1; export Fortran_COMPILER_PATH
                                ;;
        -M | --mpiexec )        shift
                                MPIEXEC_PATH=$1; export MPIEXEC_PATH
                                ;;
        -n | --num_images )     shift
                                FOR_COARRAY_NUM_IMAGES=$1
                                ;;
        -F | --fresh )          FRESH_INSTALL_ENABLED=true; export FRESH_INSTALL_ENABLED
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
        * )                     echo >&2 "-- ParaMonte - FATAL: the input flag is not recognized: $1"
                                usage
                                exit 1
    esac
    shift
done

export PMCS
export BTYPE
export LTYPE
export CAFTYPE
export MPI_ENABLED
export CFI_ENABLED
export FOR_COARRAY_NUM_IMAGES

echo >&2
echo >&2 "-- ParaMonte - current requested configuration: "

#OUTPUT=$("./configParaMonte.sh")
#echo "${OUTPUT}"
. ./configParaMonte.sh

# check flag consistencies

if [ "${CFI_ENABLED}" = "true" ]; then
    if [ "${CAFTYPE}" = "shared" ] || [ "${CAFTYPE}" = "distributed" ]; then
        echo >&2
        echo >&2 "-- ParaMonte - FATAL: incompatible input flags specified by the user:"
        echo >&2 "-- ParaMonte - FATAL:     -i | --cfi_enabled : ${CFI_ENABLED}"
        echo >&2 "-- ParaMonte - FATAL:     -c | --caf_type : ${CAFTYPE}"
        echo >&2 "-- ParaMonte - FATAL: coarray parallelism is not available in C language."
        echo >&2
        echo >&2 "-- ParaMonte - gracefully exiting."
        echo >&2
        exit 1
    fi
fi

if [ "${MPI_ENABLED}" = "true" ]; then
    if [ "${CAFTYPE}" = "shared" ] || [ "${CAFTYPE}" = "distributed" ]; then
        echo >&2
        echo >&2 "-- ParaMonte - FATAL: incompatible input flags specified by the user:"
        echo >&2 "-- ParaMonte - FATAL:     -m | --mpi_enabled : ${MPI_ENABLED}"
        echo >&2 "-- ParaMonte - FATAL:     -c | --caf_type : ${CAFTYPE}"
        echo >&2 "-- ParaMonte - FATAL: coarray parallelism cannot be mixed with MPI in the current version of ParaMonte."
        echo >&2
        echo >&2 "-- ParaMonte - gracefully exiting."
        echo >&2
        exit 1
    fi
fi

if [ "${LTYPE}" = "dynamic" ]; then
    if [ "${CAFTYPE}" = "shared" ] || [ "${CAFTYPE}" = "distributed" ]; then
        echo >&2
        echo >&2 "-- ParaMonte - FATAL: incompatible input flags specified by the user:"
        echo >&2 "-- ParaMonte - FATAL:     -l | --lib_type : ${LTYPE}"
        echo >&2 "-- ParaMonte - FATAL:     -c | --caf_type : ${CAFTYPE}"
        echo >&2 "-- ParaMonte - FATAL: ParaMonte dynamic library build with coarray parallelism currently unsupported."
        echo >&2
        echo >&2 "-- ParaMonte - gracefully exiting."
        echo >&2
        exit 1
    fi
fi


####################################################################################################################################
# define utils
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
        echo >&2 "-- ParaMonte - ParaMonte $2 appears to have succeeded."
    else
        echo >&2 
        echo >&2 "    -- ParaMonte - FATAL: ParaMonte $2 appears to have failed."
        echo >&2 "    -- ParaMonte - FATAL: If the source of the error cannot be identified,"
        echo >&2 "    -- ParaMonte - FATAL: consider a fresh installation of ParaMonte's required compilers by calling"
        echo >&2 "    -- ParaMonte - FATAL: "
        echo >&2 "    -- ParaMonte - FATAL:     ./install --fresh"
        echo >&2 "    -- ParaMonte - FATAL: "
        echo >&2 "    -- ParaMonte - FATAL: If the error happens during the installation of ParaMonte prerequisites"
        echo >&2 "    -- ParaMonte - FATAL: it is possible that the current existing GCC compiler collection installed"
        echo >&2 "    -- ParaMonte - FATAL: on your system cannot compile the downloaded version of GCC that is required"
        echo >&2 "    -- ParaMonte - FATAL: for ParaMonte build. In such case, make sure you have a GCC compiler collection"
        echo >&2 "    -- ParaMonte - FATAL: version 7.1 or newer installed on your system, with an updated PATH environmental"
        echo >&2 "    -- ParaMonte - FATAL: variable, then reinstall ParaMonte."
        echo >&2 "    -- ParaMonte - FATAL: "
        echo >&2 "    -- ParaMonte - FATAL: If all ParaMonte installation attempts fail, please report this issue at"
        echo >&2 "    -- ParaMonte - FATAL: "
        echo >&2 "    -- ParaMonte - FATAL:     https://github.com/shahmoradi/paramonte/issues"
        echo >&2 "    -- ParaMonte - FATAL: "
        echo >&2 "    -- ParaMonte - FATAL: or by contacting the ParaMonte authors directly (e.g., shahmoradi@utexas.edu)."
        echo >&2 
        echo >&2 "    -- ParaMonte - gracefully exiting."
        echo >&2 
        exit 1
    fi
}

####################################################################################################################################
# check cmake version
####################################################################################################################################

cmakeVersion="$(cmake --version)"
cmakeVersion=${cmakeVersion:14:18}
#echo >&2 "cmake version: ${cmakeVersion}"
cmakeVersion=${cmakeVersion#c*[0-9]}
cmakeVersion=${cmakeVersion%\-*}
cmakeVersionRequired=3.14.0
echo "-- ParaMonte - cmake version: ${cmakeVersion}"
echo "-- ParaMonte - cmake version required: ${cmakeVersionRequired}"
cmakeInstallEnabled=false
compareVersions "${cmakeVersion}" "${cmakeVersionRequired}"
if [ "$?" = "2" ]; then
    cmakeInstallEnabled=true
    echo >&2 "-- ParaMonte - the current cmake installation is NOT ParaMonte compatible!"
else
    echo >&2 "-- ParaMonte - the current cmake installation is ParaMonte compatible!"
fi
export cmakeInstallEnabled

####################################################################################################################################
# set local dependencies
####################################################################################################################################

ParaMonte_REQ_DIR="${ParaMonte_ROOT_DIR}/build/prerequisites"
export ParaMonte_REQ_DIR
ParaMonte_REQ_INSTALL_DIR="${ParaMonte_REQ_DIR}/prerequisites/installations"
ParaMonte_GNU_ROOT_DIR="${ParaMonte_REQ_INSTALL_DIR}/gnu/8.3.0"
ParaMonte_MPI_ROOT_DIR="${ParaMonte_REQ_INSTALL_DIR}/mpich/3.2"
ParaMonte_CAF_ROOT_DIR="${ParaMonte_REQ_INSTALL_DIR}/opencoarrays/2.8.0"
ParaMonte_CMAKE_ROOT_DIR="${ParaMonte_REQ_INSTALL_DIR}/cmake/${cmakeVersionRequired}"

ParaMonte_GNU_BIN_DIR="${ParaMonte_GNU_ROOT_DIR}/bin"
ParaMonte_CAF_BIN_DIR="${ParaMonte_CAF_ROOT_DIR}/bin"
ParaMonte_MPI_BIN_DIR="${ParaMonte_MPI_ROOT_DIR}/bin"
ParaMonte_CMAKE_BIN_DIR="${ParaMonte_CMAKE_ROOT_DIR}/bin"

ParaMonte_GNU_LIB_DIR="${ParaMonte_GNU_ROOT_DIR}/lib64"
ParaMonte_MPI_LIB_DIR="${ParaMonte_MPI_ROOT_DIR}/lib"
ParaMonte_CAF_LIB_DIR="${ParaMonte_CAF_ROOT_DIR}/lib64"

ParaMonte_CAF_WRAPPER_PATH="${ParaMonte_CAF_BIN_DIR}/caf"
export ParaMonte_CAF_WRAPPER_PATH

ParaMonte_CAF_SETUP_PATH="${ParaMonte_CAF_ROOT_DIR}/setup.sh"
export ParaMonte_CAF_SETUP_PATH

ParaMonte_CMAKE_PATH="${ParaMonte_CMAKE_BIN_DIR}/cmake"
export ParaMonte_CMAKE_PATH

####################################################################################################################################
# set compiler suite
####################################################################################################################################

if [ -z ${PMCS+x} ]
then
    ############################################################
    SUITE_LIST="intel gnu"
    ############################################################
else
    if [[ ${PMCS} == [iI][nN][tT][eE][lL] ]]; then
        PMCS=intel
        SUITE_LIST=${PMCS}
    else
        if [[ ${PMCS} == [gG][nN][uU] ]]; then
            PMCS=gnu
            SUITE_LIST=${PMCS}
        else
            echo >&2 "-- ParaMonte - FATAL: the requested compiler suite ${PMCS} is unrecognized."
            echo >&2 "-- ParaMonte - FATAL: please choose either intel or gnu, or drop the option."
            echo >&2 "-- ParaMonte - FATAL: The installer will automatically find the proper compiler suite."
            echo >&2
            echo >&2 "-- ParaMonte - gracefully exiting."
            echo >&2
            exit 1
        fi
    fi
    
fi

####################################################################################################################################
# detect the compiler suites, C/Fortran compilers and CAF/MPI libraries
####################################################################################################################################

#prereqInstallEnabled=false

gnuCCompilerName=gcc
gnuFortranCompilerName=gfortran

intelCCompilerName=icc
intelFortranCompilerName=ifort

LANG_LIST="C Fortran"

for SUITE in $SUITE_LIST
do

    echo >&2
    echo >&2 "-- ParaMonteCompiler - checking for ${SUITE} compilers and libraries presence..."
    echo >&2

    for LANG in $LANG_LIST
    do

        suiteLangCompilerName="${SUITE}${LANG}CompilerName"
        if eval "command -v ${!suiteLangCompilerName} >/dev/null 2>&1"; then

            suiteLangCompilerPath="${SUITE}${LANG}CompilerPath"
            eval "unset ${suiteLangCompilerPath}"
            eval ${suiteLangCompilerPath}='$(command -v ${!suiteLangCompilerName})'

            echo >&2 "-- ParaMonteCompiler - ${SUITE} ${!suiteLangCompilerName} detected at: ${suiteLangCompilerPath}=${!suiteLangCompilerPath}"

            # get compiler version

            suiteLangCompilerVersion="${SUITE}${LANG}CompilerVersion"

            if [ "${LANG}" = "C" ]; then

                eval ${suiteLangCompilerVersion}="$(${!suiteLangCompilerName} -dumpversion)"
                echo >&2 "-- ParaMonteCompiler - ${SUITE} ${LANG} compiler version: ${suiteLangCompilerVersion}=${!suiteLangCompilerVersion}"

            fi

            if [ "${LANG}" = "Fortran" ]; then

                cd ./auxil/

                if ${!suiteLangCompilerPath} getCompilerVersion.f90 -o getCompilerVersion.exe; then

                    chmod +x getCompilerVersion.exe
                    ./getCompilerVersion.exe && {
                        eval ${suiteLangCompilerVersion}='$(head -n 1 getCompilerVersion.tmp)'
                        echo >&2 "-- ParaMonteCompiler - ${SUITE} ${LANG} compiler version: ${suiteLangCompilerVersion}=${!suiteLangCompilerVersion}"
                        isParaMonteCompatibleCompiler=$(head -n 1 isParaMonteCompatibleCompiler.tmp)
                        if [ "$isParaMonteCompatibleCompiler" = "true" ]; then
                            echo >&2 "-- ParaMonteCompiler - ${SUITE} ${LANG} compiler is ParaMonte compatible!"
                            eval "export $suiteLangCompilerPath"
                        else
                            echo >&2 "-- ParaMonteCompiler - ${SUITE} ${LANG} compiler is not ParaMonte compatible...skipping"
                            unset ${suiteLangCompilerPath}
                        fi
                        rm *.tmp *.exe
                        #cd ..
                        #break
                    } || {
                        echo >&2 "-- ParaMonteCompiler - failed to detect the ${SUITE} ${LANG} compiler version...skipping"
                        unset ${suiteLangCompilerPath}
                    }

                else

                    echo >&2 "-- ParaMonteCompiler - failed to detect the ${SUITE} ${LANG} compiler version...skipping"
                    unset ${suiteLangCompilerPath}

                fi

                cd ..

            fi

        else

            echo >&2 "-- ParaMonteCompiler - ${SUITE} ${!suiteLangCompilerName} not found."
            unset ${suiteLangCompilerPath}

        fi

    done

done

####################################################################################################################################
# detect MPI wrappers
####################################################################################################################################

gnuCMpiWrapperName="mpic++"
gnuFortranMpiWrapperName=mpifort

intelCMpiWrapperName=mpiicc
intelFortranMpiWrapperName=mpiifort

for SUITE in $SUITE_LIST
do

    echo >&2
    echo >&2 "-- ParaMonteMPI - checking for ${SUITE} MPI wrappers and libraries presence..."
    echo >&2

    for LANG in $LANG_LIST
    do

        suiteLangMpiWrapperName="${SUITE}${LANG}MpiWrapperName"
        suiteLangMpiWrapperPath="${SUITE}${LANG}MpiWrapperPath"
        if eval "command -v ${!suiteLangMpiWrapperName} >/dev/null 2>&1"; then

            eval "unset ${suiteLangMpiWrapperPath}"
            eval ${suiteLangMpiWrapperPath}='$(command -v ${!suiteLangMpiWrapperName})'

            echo >&2 "-- ParaMonteMPI - ${SUITE} ${!suiteLangMpiWrapperName} detected at: ${suiteLangMpiWrapperPath}=${!suiteLangMpiWrapperPath}"

        else

            echo >&2 "-- ParaMonteMPI - failed to detect the ${SUITE} ${LANG} MPI wrapper...skipping"
            unset ${suiteLangMpiWrapperPath}

        fi

    done

done

####################################################################################################################################
# detect CAF wrapper
####################################################################################################################################

echo >&2  

CAF_ENABLED=false
if [ "${CAFTYPE}" != "none" ]; then
    CAF_ENABLED=true
    #if [ -z ${intelFortranMpiWrapperPath+x} ]; then
        if command -v caf >/dev/null 2>&1; then
            cafCompilerPath=$(command -v caf)
            echo >&2 "-- ParaMonteCAF - OpenCoarrays Fortran compiler wrapper detected at: ${cafCompilerPath}"
            cafVersion="$(caf -dumpversion)"
            cafVersionRequired="7.3.0"
            echo >&2 "-- ParaMonteCAF - caf version: ${cafVersion}"
            echo >&2 "-- ParaMonteCAF - caf version required: ${cafVersionRequired}"
            compareVersions "$cafVersion" "$cafVersionRequired"
            if [ "$?" = "2" ]; then
                cafInstallEnabled=true
                mpiInstallEnabled=true
                gnuInstallEnabled=true
                #PMCS=caf
                #COMPILER_VERSION=unknownversion
            else
                #if [ "$(printf '%s\n' "$cafVersionRequired" "$currentver" | sort -V | head -n1)" = "$cafVersionRequired" ]; then
                echo >&2 "-- ParaMonteCAF - OpenCoarrays Fortran compiler wrapper is ParaMonte compatible!"
            fi
        else
            cafInstallEnabled=true
            mpiInstallEnabled=true
            gnuInstallEnabled=true
        fi
    #fi
    if [ "${cafInstallEnabled}" = "true" ]; then
        echo >&2 "-- ParaMonteCAF - NOTE: OpenCoarrays caf compiler wrapper could not be found on your system."
        echo >&2 
    fi
fi
export CAF_ENABLED

####################################################################################################################################
# set ParaMonte compiler suite
####################################################################################################################################

prereqInstallAllowed=false
if [ -z ${PMCS+x} ]; then prereqInstallAllowed=true; fi

if [ -z ${Fortran_COMPILER_PATH+x} ]; then

    if [ -z ${PMCS+x} ]; then

        # if no preference, then priority is with intel

        if ! ${intelFortranCompilerPath+false}; then
            if [ "${MPI_ENABLED}" = "true" ] || [ "${CAF_ENABLED}" = "true" ]; then
                if ! ${intelFortranMpiWrapperPath+false}; then
                    PMCS=intel
                    COMPILER_VERSION=${intelFortranCompilerVersion}
                    Fortran_COMPILER_PATH="${intelFortranCompilerPath}"
                    MPIEXEC_PATH=$(dirname "${intelFortranMpiWrapperPath}")/mpiexec
                    #MPIEXEC_PATH=$(command -v mpiexec)
                    prereqInstallAllowed=false
                else
                    PMCS=gnu
                fi
            else
                PMCS=intel
                COMPILER_VERSION=${intelFortranCompilerVersion}
                Fortran_COMPILER_PATH="${intelFortranCompilerPath}"
                prereqInstallAllowed=false
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
                        MPIEXEC_PATH=$(dirname "${intelFortranMpiWrapperPath}")/mpiexec
                    else
                        errorOccurred=true
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
            echo >&2 "-- ParaMonte - FATAL: the Fortran compiler and wrapper components of the"
            echo >&2 "-- ParaMonte - FATAL: requested compiler suite ${PMCS} could not be detected."
            echo >&2 "-- ParaMonte - FATAL: Please make sure the ${PMCS} compiler suite is installed on your system"
            echo >&2 "-- ParaMonte - FATAL: and avariable in the environmental variable PATH of your shell."
            echo >&2
            echo >&2 "-- ParaMonte - gracefully exiting."
            echo >&2
            exit 1
        fi

    fi

    if [ "${PMCS}" = "gnu" ]; then

        if ! ${gnuFortranCompilerPath+false}; then
            COMPILER_VERSION=${gnuFortranCompilerVersion}
            Fortran_COMPILER_PATH="${gnuFortranCompilerPath}"
            gnuInstallEnabled=false
        else
            gnuInstallEnabled=true
            if [ "${prereqInstallAllowed}" = "true" ]; then
                echo >&2
                echo >&2 "-- ParaMonte - WARNING: GNU Fortran compiler could not be found on your system."
                echo >&2 "-- ParaMonte - WARNING: If you do not have GNU compiler suite installed on your system."
                echo >&2 "-- ParaMonte - WARNING: ParaMonte may be able to install the compiler suite for you. To do so,"
                echo >&2 "-- ParaMonte - WARNING: drop the input argument -s or --compiler_suite when calling the script."
                echo >&2
            fi
        fi

        if [ "${MPI_ENABLED}" = "true" ]; then
            if [ -z ${MPIEXEC_PATH+x} ]; then
                mpiInstallEnabled=true
                gnuInstallEnabled=true
                if [ "${prereqInstallAllowed}" = "true" ]; then
                    echo >&2
                    echo >&2 "-- ParaMonte - WARNING: GNU Fortran compiler could not be found on your system."
                    echo >&2 "-- ParaMonte - WARNING: If you do not have GNU compiler suite installed on your system."
                    echo >&2 "-- ParaMonte - WARNING: ParaMonte may be able to install the compiler suite for you. To do so,"
                    echo >&2 "-- ParaMonte - WARNING: drop the input argument -s or --compiler_suite when calling the script."
                    echo >&2
                fi
            else
                mpiInstallEnabled=false
            fi
        else
            mpiInstallEnabled=false
        fi

        if [ "${CAF_ENABLED}" = "true" ]; then
            if ! ${cafCompilerPath+false}; then
                COMPILER_VERSION=unknownversion
                Fortran_COMPILER_PATH="${cafCompilerPath}"
                cafInstallEnabled=false
            else
                cafInstallEnabled=true
                mpiInstallEnabled=true
                gnuInstallEnabled=true
            fi
        fi

    fi

else # if fortran compiler path defined

    cafInstallEnabled=false
    mpiInstallEnabled=false
    gnuInstallEnabled=false

    PMCS=unknownsuite
    COMPILER_VERSION=unknownversion
    Fortran_COMPILER_NAME=${Fortran_COMPILER_PATH##*/}
    echo >&2
    echo >&2 "-- ParaMonteCompiler - user-requested compiler path: ${Fortran_COMPILER_PATH}"
    echo >&2 "-- ParaMonteCompiler - user-requested compiler name: ${Fortran_COMPILER_NAME}"
    if [ "${Fortran_COMPILER_NAME}" = "gfortran" ] || [ "${Fortran_COMPILER_NAME}" = "caf" ]; then
        PMCS=gnu
    fi
    if [ "${Fortran_COMPILER_NAME}" = "ifort" ]; then
        PMCS=intel
    fi

    if [ "${MPI_ENABLED}" = "true" ] && [ -z ${MPIEXEC_PATH+x} ]; then
        mpiInstallEnabled=true
        gnuInstallEnabled=true
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
#    echo >&2 "-- ParaMonteMPI - intel Fortran MPI compiler wrapper: ${intelFortranMpiWrapperPath}"
#    echo >&2 "-- ParaMonteMPI - intel Fortran compiler: ${intelFortranCompilerPath}"
#    echo >&2 
#    echo >&2 "-- ParaMonteMPI - gnu Fortran MPI compiler wrapper: ${gnuFortranMpiWrapperPath}"
#    echo >&2 "-- ParaMonteMPI - gnu Fortran compiler: ${gnuFortranCompilerPath}"
#    echo >&2 
#    if [ -z ${intelFortranMpiWrapperPath+x} ] || [ -z ${intelFortranCompilerPath+x} ]; then
#        #if [ -z ${gnuFortranMpiWrapperPath+x} ] || [ -z ${gnuFortranCompilerPath+x} ]; then
#            mpiInstallEnabled=true
#            gnuInstallEnabled=true
#        #fi
#    fi
#fi

####################################################################################################################################
# set up ParaMonte library prerequisites
####################################################################################################################################

if [ "${FRESH_INSTALL_ENABLED}" = "true" ]; then
    prereqInstallAllowed=true
    cmakeInstallEnabled=true
    cafInstallEnabled=true
    mpiInstallEnabled=true
    gnuInstallEnabled=true
fi

# chmod 777 -R "${ParaMonte_ROOT_DIR}/auxil/prerequisites"

if [ "${prereqInstallAllowed}" = "true" ]; then

    if [ "${cafInstallEnabled}" = "true" ] || [ "${mpiInstallEnabled}" = "true" ] || [ "${gnuInstallEnabled}" = "true" ] || [ "${cmakeInstallEnabled}" = "true" ]; then

        #ParaMonte_REQ_DIR="${ParaMonte_ROOT_DIR}/build/prerequisites"
        #ParaMonte_CAF_SETUP_PATH="${ParaMonte_REQ_DIR}/prerequisites/installations/opencoarrays/2.8.0/setup.sh"

        if [ "${FRESH_INSTALL_ENABLED}" = "true" ]; then
            rm -rf "${ParaMonte_REQ_DIR}"
        fi

        answer=y
        if [ ! -d "${ParaMonte_REQ_DIR}" ]; then

            echo >&2
            echo >&2 "-- ParaMonte - WARNING: ParaMonte build with the requested configuration requires the installations"
            echo >&2 "-- ParaMonte - WARNING: of either OpenCoarrays, MPICH MPI library, GNU compilers, or CMAKE on your system."
            echo >&2 "-- ParaMonte - WARNING: ParaMonte can install all the prerequisites on your system from the web, if needed."
            echo >&2 "-- ParaMonte - WARNING: Thre prerequisites (GNU compilers) may occupy up to 5Gb of your system's memory."
            echo >&2
            if [ "${YES_TO_ALL_DISABLED}" = "true" ]; then
                answerNotGiven=true
                while [ "${answerNotGiven}" = "true" ]; do
                    read -p "-- ParaMonte - Do you wish to continue with the installation of the prerequisites (y/n)? " answer
                    if [[ $answer == [yY] || $answer == [yY][eE][sS] ]]; then
                        answer=y
                        answerNotGiven=false
                    fi
                    if [[ $answer == [nN] || $answer == [nN][oO] ]]; then
                        answer=n
                        answerNotGiven=false
                    fi
                    if [ "${answerNotGiven}" = "true" ]; then
                        echo >&2 "-- ParaMonte - please enter either y or n"
                    fi
                done
            else
                echo >&2 "-- ParaMonte - Do you wish to continue with the installation of the prerequisites (y/n)? y"
                answer=y
            fi
            echo >&2

            if [ "${answer}" = "y" ]; then

                echo >&2 "-- ParaMonte - generating directory: ${ParaMonte_REQ_DIR}/"
                mkdir -p "${ParaMonte_REQ_DIR}/"

                # cp -rv "${ParaMonte_ROOT_DIR}/auxil/prerequisites" "${ParaMonte_REQ_DIR}/../"
                cp -rv "${ParaMonte_ROOT_DIR}/auxil/prerequisites.tar.gz" "${ParaMonte_REQ_DIR}/../"
                verify $? "installation setup of prerequisites"
                if ! [ -d "${ParaMonte_REQ_DIR}" ]; then
                    mkdir "${ParaMonte_REQ_DIR}"
                fi
                (cd "${ParaMonte_REQ_DIR}/../" && tar xvzf opencoarrays.tar.gz -C prerequisites --strip-components 1)
                verify $? "unpacking of prerequisites"
                # chmod +x -R "${ParaMonte_REQ_DIR}"

            else

                echo >&2 "-- ParaMonte - WARNING: ParaMonte installation will proceed with no guarantee of success."
                echo >&2
                #exit 1

            fi

        fi

        if [ "${answer}" = "y" ]; then

            # check cmake

            if [ "${cmakeInstallEnabled}" = "true" ]; then
                # ParaMonte_CMAKE_BIN_DIR="${ParaMonte_REQ_DIR}/prerequisites/installations/cmake/${cmakeVersionRequired}/bin"
                # ParaMonte_CMAKE_PATH="${ParaMonte_CMAKE_BIN_DIR}/cmake"
                if [[ -f "${ParaMonte_CMAKE_PATH}" ]]; then
                    echo >&2 "-- ParaMonte - cmake ${cmakeVersionRequired} detected."
                    #cmakeInstallEnabled=false
                else
                    echo >&2 "-- ParaMonte - cmake ${cmakeVersionRequired} missing."
                    echo >&2 "-- ParaMonte - installing the prerequisites...this can take a while."
                    (cd ${ParaMonte_REQ_DIR} && chmod +x ./install.sh && ./install.sh  --yes-to-all --package cmake --install-version ${cmakeVersionRequired} )
                    verify $? "installation of cmake"
                fi
            fi

            # check mpi

            CURRENT_PKG="MPICH library"
            if [ "${mpiInstallEnabled}" = "true" ]; then
                #ParaMonte_MPI_BIN_DIR="${ParaMonte_REQ_DIR}/prerequisites/installations/mpich/3.2/bin"
                MPIEXEC_PATH="${ParaMonte_MPI_BIN_DIR}/mpiexec"
                if [[ -f "${MPIEXEC_PATH}" ]]; then
                    echo >&2 "-- ParaMonte - ${CURRENT_PKG} detected."
                    #mpiInstallEnabled=false
                else
                    ##########################################################################
                    echo >&2 "-- ParaMonte - ${CURRENT_PKG} missing."
                    echo >&2 "-- ParaMonte - installing the prerequisites...this can take a while."
                    (cd ${ParaMonte_REQ_DIR} && chmod +x ./install.sh && ./install.sh --yes-to-all ${GCC_BOOTSTRAP} ) || 
                    {
                        if [ -z ${GCC_BOOTSTRAP+x} ]; then
                            echo >&2
                            read -p "-- ParaMonte - ${CURRENT_PKG} installation failed. Shall I retry with bootstrap (y/n)? " answer 
                            echo >&2
                            if [[ $answer == [yY] || $answer == [yY][eE][sS] ]]; then
                                GCC_BOOTSTRAP="--bootstrap"
                                (cd ${ParaMonte_REQ_DIR} && chmod +x ./install.sh && ./install.sh --yes-to-all ${GCC_BOOTSTRAP} )
                                verify $? "${CURRENT_PKG} installation"
                            else
                                verify 1 "${CURRENT_PKG} installation"
                            fi
                        else
                            verify 1 "${CURRENT_PKG} installation"
                        fi
                    }
                    ##########################################################################
                fi
            fi

            # check gnu

            CURRENT_PKG="GNU compiler collection"
            if [ "${gnuInstallEnabled}" = "true" ]; then
                #ParaMonte_GNU_BIN_DIR="${ParaMonte_REQ_DIR}/prerequisites/installations/gnu/8.3.0/bin"
                Fortran_COMPILER_PATH="${ParaMonte_GNU_BIN_DIR}/gfortran"
                if [[ -f "${Fortran_COMPILER_PATH}" ]]; then
                    echo >&2 "-- ParaMonte - ${CURRENT_PKG} detected."
                    #gnuInstallEnabled=false
                else
                    ##########################################################################
                    echo >&2 "-- ParaMonte - ${CURRENT_PKG} missing."
                    echo >&2 "-- ParaMonte - installing the prerequisites...this can take a while."
                    (cd ${ParaMonte_REQ_DIR} && chmod +x ./install.sh && ./install.sh --yes-to-all ${GCC_BOOTSTRAP} ) || 
                    {
                        if [ -z ${GCC_BOOTSTRAP+x} ]; then
                            echo >&2
                            read -p "-- ParaMonte - ${CURRENT_PKG} installation failed. Shall I retry with bootstrap (y/n)? " answer 
                            echo >&2
                            if [[ $answer == [yY] || $answer == [yY][eE][sS] ]]; then
                                GCC_BOOTSTRAP="--bootstrap"
                                (cd ${ParaMonte_REQ_DIR} && chmod +x ./install.sh && ./install.sh --yes-to-all ${GCC_BOOTSTRAP} )
                                verify $? "${CURRENT_PKG} installation"
                            else
                                verify 1 "${CURRENT_PKG} installation"
                            fi
                        else
                            verify 1 "${CURRENT_PKG} installation"
                        fi
                    }
                    ##########################################################################
                fi
            fi

            # # set up setup.sh file

            # SETUP_FILE_PATH="${ParaMonte_ROOT_DIR}/build/setup.sh"
            # export SETUP_FILE_PATH

            # echo "# ParaMonte runtime environment setup script." > ${SETUP_FILE_PATH}
            # echo "# Source this Bash script in your Bash environment like," >> ${SETUP_FILE_PATH}
            # echo "#     source ./setup.sh" >> ${SETUP_FILE_PATH}
            # echo "# before compiling your source files and linking with ParaMonte library." >> ${SETUP_FILE_PATH}
            # echo "" >> ${SETUP_FILE_PATH}

            # check caf

            CURRENT_PKG="OpenCoarrays compiler wrapper"
            if [ "${cafInstallEnabled}" = "true" ]; then
                #ParaMonte_CAF_WRAPPER_PATH="${ParaMonte_CAF_BIN_DIR}/caf"
                if [[ -f "${ParaMonte_CAF_WRAPPER_PATH}" ]]; then
                    echo >&2 "-- ParaMonte - ${CURRENT_PKG} detected."
                    #cafInstallEnabled=false
                else
                    ##########################################################################
                    echo >&2 "-- ParaMonte - ${CURRENT_PKG} missing."
                    echo >&2 "-- ParaMonte - installing the prerequisites...this can take a while."
                    (cd ${ParaMonte_REQ_DIR} && chmod +x ./install.sh && ./install.sh ${GCC_BOOTSTRAP} --yes-to-all) || 
                    {
                        if [ -z ${GCC_BOOTSTRAP+x} ]; then
                            echo >&2
                            read -p "-- ParaMonte - ${CURRENT_PKG} installation failed. Shall I retry with bootstrap (y/n)? " answer 
                            echo >&2
                            if [[ $answer == [yY] || $answer == [yY][eE][sS] ]]; then
                                GCC_BOOTSTRAP="--bootstrap"
                                (cd ${ParaMonte_REQ_DIR} && chmod +x ./install.sh && ./install.sh ${GCC_BOOTSTRAP} --yes-to-all)
                                verify $? "${CURRENT_PKG} installation"
                            else
                                verify 1 "${CURRENT_PKG} installation"
                            fi
                        else
                            verify 1 "${CURRENT_PKG} installation"
                        fi
                    }
                    ##########################################################################
                fi
                Fortran_COMPILER_PATH="${ParaMonte_CAF_WRAPPER_PATH}"
            fi

            if [ -f "${ParaMonte_CAF_SETUP_PATH}" ]; then
                source "${ParaMonte_CAF_SETUP_PATH}"
                # echo "" >> ${SETUP_FILE_PATH}
                # echo "source ${ParaMonte_CAF_SETUP_PATH}" >> ${SETUP_FILE_PATH}
                # echo "" >> ${SETUP_FILE_PATH}
            fi

        fi

    fi

fi

####################################################################################################################################
# set up PATH & LD_LIBRARY_PATH
####################################################################################################################################

if [ -d "${ParaMonte_CMAKE_BIN_DIR}" ]; then
    if [ -z ${PATH+x} ]; then
        PATH="${ParaMonte_CMAKE_BIN_DIR}"
    else
        if [[ ":$PATH:" != *":${ParaMonte_CMAKE_BIN_DIR}:"* ]]; then
            PATH="${ParaMonte_CMAKE_BIN_DIR}:${PATH}"
        fi
    fi
fi

####################################################################################################################################
# set up setup.sh file
####################################################################################################################################

if [ "${PMCS}" = "gnu" ] && [ "${prereqInstallAllowed}" = "true" ]; then

    SETUP_FILE_PATH="${ParaMonte_ROOT_DIR}/build/setup.sh"
    export SETUP_FILE_PATH
    {
    echo "# ParaMonte runtime environment setup script."
    echo "# Source this Bash script in your Bash environment like,"
    echo "#     source ./setup.sh"
    echo "# before compiling your source files and linking with ParaMonte library."
    echo ""
    } > ${SETUP_FILE_PATH}
    chmod +x ${SETUP_FILE_PATH}

    if [[ -f "${SETUP_FILE_PATH}" ]]; then
        if [[ -f "${ParaMonte_CAF_SETUP_PATH}" ]]; then
            {
            echo ""
            echo "source ${ParaMonte_CAF_SETUP_PATH}"
            echo ""
            } >> ${SETUP_FILE_PATH}
            if [ "${LTYPE}" = "dynamic" ]; then
                {
                echo "if [ -z \${LD_LIBRARY_PATH+x} ]; then"
                echo "    LD_LIBRARY_PATH=."
                echo "else"
                echo "    FILE_DIR=\"\$( cd \"\$( dirname \"\${BASH_SOURCE[0]}\" )\" >/dev/null 2>&1 && pwd )\""
                echo "    if [[ \":\$LD_LIBRARY_PATH:\" != *\":${FILE_DIR}:\"* ]]; then"
                echo "        LD_LIBRARY_PATH=\"${FILE_DIR}:\${LD_LIBRARY_PATH}\""
                echo "    fi"
                echo "fi"
                echo "export LD_LIBRARY_PATH"
                echo ""
                } >> ${SETUP_FILE_PATH}
            fi
        fi
    fi

    if [ -d "${ParaMonte_GNU_BIN_DIR}" ]; then
        if [ -z ${PATH+x} ]; then
            PATH="${ParaMonte_GNU_BIN_DIR}"
        else
            if [[ ":$PATH:" != *":${ParaMonte_GNU_BIN_DIR}:"* ]]; then
                PATH="${ParaMonte_GNU_BIN_DIR}:${PATH}"
            fi
        fi
        if [ -z ${LD_LIBRARY_PATH+x} ]; then
            LD_LIBRARY_PATH="${ParaMonte_GNU_LIB_DIR}"
        else
            if [[ ":$LD_LIBRARY_PATH:" != *":${ParaMonte_GNU_LIB_DIR}:"* ]]; then
                LD_LIBRARY_PATH="${ParaMonte_GNU_LIB_DIR}:${LD_LIBRARY_PATH}"
            fi
        fi
        if [[ -f "${SETUP_FILE_PATH}" ]]; then
            {
            echo "if [ -z \${PATH+x} ]; then"
            echo "    export PATH=\"${ParaMonte_GNU_BIN_DIR}\""
            echo "else"
            echo "    if [[ \":\$PATH:\" != *\":${ParaMonte_GNU_BIN_DIR}:\"* ]]; then"
            echo "        export PATH=\"${ParaMonte_GNU_BIN_DIR}:\${PATH}\""
            echo "    fi"
            echo "fi"
            echo "if [ -z \${LD_LIBRARY_PATH+x} ]; then"
            echo "    export LD_LIBRARY_PATH=\"${ParaMonte_GNU_LIB_DIR}\""
            echo "else"
            echo "    if [[ \":\$LD_LIBRARY_PATH:\" != *\":${ParaMonte_GNU_LIB_DIR}:\"* ]]; then"
            echo "        export LD_LIBRARY_PATH=\"${ParaMonte_GNU_LIB_DIR}:\${LD_LIBRARY_PATH}\""
            echo "    fi"
            echo "fi"
            } >> ${SETUP_FILE_PATH}
        fi
    fi

    if [ -d "${ParaMonte_MPI_BIN_DIR}" ]; then
        if [ -z ${PATH+x} ]; then
            PATH="${ParaMonte_MPI_BIN_DIR}"
        else
            if [[ ":$PATH:" != *":${ParaMonte_MPI_BIN_DIR}:"* ]]; then
                PATH="${ParaMonte_MPI_BIN_DIR}:${PATH}"
            fi
        fi
        if [ -z ${LD_LIBRARY_PATH+x} ]; then
            LD_LIBRARY_PATH="${ParaMonte_MPI_LIB_DIR}"
        else
            if [[ ":$LD_LIBRARY_PATH:" != *":${ParaMonte_MPI_LIB_DIR}:"* ]]; then
                LD_LIBRARY_PATH="${ParaMonte_MPI_LIB_DIR}:${LD_LIBRARY_PATH}"
            fi
        fi
        if [[ -f "${SETUP_FILE_PATH}" ]]; then
            {
            echo "if [ -z \${PATH+x} ]; then"
            echo "    export PATH=\"${ParaMonte_MPI_BIN_DIR}\""
            echo "else"
            echo "    if [[ \":\$PATH:\" != *\":${ParaMonte_MPI_BIN_DIR}:\"* ]]; then"
            echo "        export PATH=\"${ParaMonte_MPI_BIN_DIR}:\${PATH}\""
            echo "    fi"
            echo "fi"
            echo "if [ -z \${LD_LIBRARY_PATH+x} ]; then"
            echo "    export LD_LIBRARY_PATH=\"${ParaMonte_MPI_LIB_DIR}\""
            echo "else"
            echo "    if [[ \":\$LD_LIBRARY_PATH:\" != *\":${ParaMonte_MPI_LIB_DIR}:\"* ]]; then"
            echo "        export LD_LIBRARY_PATH=\"${ParaMonte_MPI_LIB_DIR}:\${LD_LIBRARY_PATH}\""
            echo "    fi"
            echo "fi"
            } >> ${SETUP_FILE_PATH}
        fi
    fi

    if [ -d "${ParaMonte_CAF_BIN_DIR}" ]; then
        if [ -z ${PATH+x} ]; then
            PATH="${ParaMonte_CAF_BIN_DIR}"
        else
            if [[ ":$PATH:" != *":${ParaMonte_CAF_BIN_DIR}:"* ]]; then
                PATH="${ParaMonte_CAF_BIN_DIR}:${PATH}"
            fi
        fi
        if [ -z ${LD_LIBRARY_PATH+x} ]; then
            LD_LIBRARY_PATH="${ParaMonte_CAF_LIB_DIR}"
        else
            if [[ ":$LD_LIBRARY_PATH:" != *":${ParaMonte_CAF_LIB_DIR}:"* ]]; then
                LD_LIBRARY_PATH="${ParaMonte_CAF_LIB_DIR}:${LD_LIBRARY_PATH}"
            fi
        fi
        if [[ -f "${SETUP_FILE_PATH}" ]]; then
            {
            echo "if [ -z \${PATH+x} ]; then"
            echo "    export PATH=\"${ParaMonte_CAF_BIN_DIR}\""
            echo "else"
            echo "    if [[ \":\$PATH:\" != *\":${ParaMonte_CAF_BIN_DIR}:\"* ]]; then"
            echo "        export PATH=\"${ParaMonte_CAF_BIN_DIR}:\${PATH}\""
            echo "    fi"
            echo "fi"
            echo "if [ -z \${LD_LIBRARY_PATH+x} ]; then"
            echo "    export LD_LIBRARY_PATH=\"${ParaMonte_CAF_LIB_DIR}\""
            echo "else"
            echo "    if [[ \":\$LD_LIBRARY_PATH:\" != *\":${ParaMonte_CAF_LIB_DIR}:\"* ]]; then"
            echo "        export LD_LIBRARY_PATH=\"${ParaMonte_CAF_LIB_DIR}:\${LD_LIBRARY_PATH}\""
            echo "    fi"
            echo "fi"
            } >> ${SETUP_FILE_PATH}
        fi
    fi

    export PATH
    export LD_LIBRARY_PATH

fi

####################################################################################################################################
# set Fortran compiler version one last time
####################################################################################################################################

echo >&2
echo >&2 "-- ParaMonteCompiler - Fortran compiler path: ${Fortran_COMPILER_PATH}"
echo >&2 "-- ParaMonteCompiler - MPI mpiexec path: ${MPIEXEC_PATH}"
echo >&2

if [ "${PMCS}" = "gnu" ] || [ "${COMPILER_VERSION}" = "unknownversion" ]; then

    cd ./auxil/

    LANG=Fortran
    isUnknownVersion=false
    if ${Fortran_COMPILER_PATH} getCompilerVersion.f90 -o getCompilerVersion.exe; then

        chmod +x getCompilerVersion.exe
        ./getCompilerVersion.exe && {
            COMPILER_VERSION=$(head -n 1 getCompilerVersion.tmp)
            echo >&2 "-- ParaMonteCompiler - ${PMCS} ${LANG} compiler version: ${COMPILER_VERSION}"
            isParaMonteCompatibleCompiler=$(head -n 1 isParaMonteCompatibleCompiler.tmp)
            if [ "$isParaMonteCompatibleCompiler" = "true" ]; then
                echo >&2 "-- ParaMonteCompiler - ${PMCS} ${LANG} compiler is ParaMonte compatible!"
            else
                echo >&2 "-- ParaMonteCompiler - ${SUITE} ${LANG} compiler is not ParaMonte compatible..."
                echo >&2 "-- ParaMonteCompiler - ParaMonte installation failed."
                echo >&2
                echo >&2 "-- ParaMonte - gracefully exiting."
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

    cd ..

    if [ "${isUnknownVersion}" = "true" ]; then
        echo >&2 "-- ParaMonteCompiler - failed to detect the ${PMCS} ${LANG} compiler version...skipping"
        COMPILER_VERSION=unknownversion
    fi

fi

####################################################################################################################################
# set ParaMonte build dir
####################################################################################################################################

export PMCS
export COMPILER_VERSION
echo >&2 "-- ParaMonte - selected compiler suite: ${PMCS}"
echo >&2 "-- ParaMonte - selected compiler version: ${COMPILER_VERSION}"

unset PARALLELIZATION_DIR
if [ "${OMP_ENABLED}" = "true" ]; then PARALLELIZATION_DIR=${PARALLELIZATION_DIR}omp; fi
if [ "${MPI_ENABLED}" = "true" ]; then PARALLELIZATION_DIR=${PARALLELIZATION_DIR}mpi; fi
if [ "${CAF_ENABLED}" = "true" ]; then PARALLELIZATION_DIR=${PARALLELIZATION_DIR}caf${CAFTYPE}; fi
if [ -z ${PARALLELIZATION_DIR+x} ]; then PARALLELIZATION_DIR=serial; fi
export PARALLELIZATION_DIR

unset MEMORY_ALLOCATION
if [ "${HEAP_ARRAY_ENABLED}" = "true" ]; then MEMORY_ALLOCATION="heap"; fi
if [ "${HEAP_ARRAY_ENABLED}" = "false" ]; then MEMORY_ALLOCATION="stack"; fi
if [ -z ${MEMORY_ALLOCATION+x} ]; then
    echo >&2 "-- ParaMonte - FATAL: HEAP_ARRAY_ENABLED must be set to either true or false."
    echo >&2 "-- ParaMonte - FATAL: you have provided HEAP_ARRAY_ENABLED=${HEAP_ARRAY_ENABLED}"
    echo >&2
    echo >&2 "-- ParaMonte - gracefully exiting."
    echo >&2
    exit 1
fi
export MEMORY_ALLOCATION

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set ParaMonte library build directories
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

echo >&2 
echo >&2 "-- ParaMonte - setting up build directories..."
echo >&2 

ParaMonte_BLD_DIR=${ParaMonte_ROOT_DIR}/build/${PLATFORM}${ARCHITECTURE}/${PMCS}/${COMPILER_VERSION}/${BTYPE}/${LTYPE}/${MEMORY_ALLOCATION}/${PARALLELIZATION_DIR}
if [ ${CFI_ENABLED} = "true" ]; then ParaMonte_BLD_DIR=${ParaMonte_BLD_DIR}/C; fi
if [ ${CFI_ENABLED} = "false" ]; then ParaMonte_BLD_DIR=${ParaMonte_BLD_DIR}/Fortran; fi
if [ -z ${CFI_ENABLED+x} ]; then
    echo >&2 "-- ParaMonte - FATAL: CFI_ENABLED must be set to either true or false."
    echo >&2 "-- ParaMonte - FATAL: you have provided CFI_ENABLED=${CFI_ENABLED}"
    echo >&2
    echo >&2 "-- ParaMonte - gracefully exiting."
    echo >&2
    exit 1
fi
export ParaMonte_BLD_DIR

echo >&2 "-- ParaMonte - ParaMonte build directory: ${ParaMonte_BLD_DIR}"
if [ -d "${ParaMonte_BLD_DIR}" ]; then
    echo >&2 "-- ParaMonte - ParaMonte build directory already exists. skipping..."
else
    echo >&2 "-- ParaMonte - generating ParaMonte build directory..."
    mkdir -p ${ParaMonte_BLD_DIR}
fi
echo >&2 "-- ParaMonte - all generated build files will be stored at: ${ParaMonte_BLD_DIR}"

# set object/module/lib files directories

ParaMonte_OBJ_DIR=${ParaMonte_BLD_DIR}/obj; export ParaMonte_OBJ_DIR
ParaMonte_MOD_DIR=${ParaMonte_BLD_DIR}/mod; export ParaMonte_MOD_DIR
ParaMonte_LIB_DIR=${ParaMonte_BLD_DIR}/lib; export ParaMonte_LIB_DIR

echo >&2 "-- ParaMonte - ParaMonte object files directory: ${ParaMonte_OBJ_DIR}"
echo >&2 "-- ParaMonte - ParaMonte module files directory: ${ParaMonte_MOD_DIR}"
echo >&2 "-- ParaMonte - ParaMonte library files directory: ${ParaMonte_LIB_DIR}"

# make bin directory

ParaMonte_BIN_DIR=${ParaMonte_ROOT_DIR}/bin
echo >&2 "-- ParaMonte - ParaMonte binaries directory: ${ParaMonte_BIN_DIR}"
if [[ -d "${ParaMonte_BIN_DIR}" ]]; then
    echo >&2 "-- ParaMonte - ParaMonte binaries directory already exists. skipping..."
else
    echo >&2 "-- ParaMonte - generating ParaMonte binaries directory..."
    mkdir "${ParaMonte_BIN_DIR}/"
fi
export ParaMonte_BIN_DIR

####################################################################################################################################
# set ParaMonte library source directories
####################################################################################################################################

           ParaMonteTest_SRC_DIR=${ParaMonte_ROOT_DIR}/src/test
        ParaMonteExample_SRC_DIR=${ParaMonte_ROOT_DIR}/example
      ParaMonteInterface_SRC_DIR=${ParaMonte_ROOT_DIR}/src/interface
     ParaMonteInterfaceC_SRC_DIR=${ParaMonteInterface_SRC_DIR}/C
ParaMonteInterfacePython_SRC_DIR=${ParaMonteInterface_SRC_DIR}/Python
     ParaMontePythonTest_SRC_DIR=${ParaMonteInterfacePython_SRC_DIR}/test

export ParaMonteTest_SRC_DIR
export ParaMonteExample_SRC_DIR
export ParaMontePythonTest_SRC_DIR
export ParaMonteInterface_SRC_DIR
export ParaMonteInterfaceC_SRC_DIR
export ParaMonteInterfacePython_SRC_DIR

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
    FC_OPTION="-DFC=${Fortran_COMPILER_PATH}"
fi
if [ -z ${MPIEXEC_PATH+x} ]; then
    MPIEXEC_OPTION=""
else
    MPIEXEC_OPTION="-DMPIEXEC_EXECUTABLE=${MPIEXEC_PATH}"
fi

echo >&2
echo >&2 "-- ParaMonte - CMAKE Fortran compiler option: ${FC_OPTION}"
echo >&2 "-- ParaMonte - CMAKE mpiexec option: ${MPIEXEC_OPTION}"
echo >&2

if [ "${gnuInstallEnabled}" = "true" ] || [ "${mpiInstallEnabled}" = "true" ] || [ "${cafInstallEnabled}" = "true" ]; then
    ParaMonte_CAF_SETUP_PATH_CMD="source ${ParaMonte_CAF_SETUP_PATH}"
else
    ParaMonte_CAF_SETUP_PATH_CMD=""
fi

(cd ${ParaMonte_BLD_DIR} && \
${ParaMonte_CAF_SETUP_PATH_CMD} && \
cmake  \
--verbose=1 \
${FC_OPTION} \
${MPIEXEC_OPTION} \
-DPMCS=${PMCS} \
-DMPI_ENABLED=${MPI_ENABLED} \
-DCAFTYPE=${CAFTYPE} \
-DBTYPE=${BTYPE} \
-DLTYPE=${LTYPE} \
-DHEAP_ARRAY_ENABLED=${HEAP_ARRAY_ENABLED} \
-DCFI_ENABLED=${CFI_ENABLED} \
-DOMP_ENABLED=${OMP_ENABLED} \
${ParaMonte_ROOT_DIR} \
)
verify $? "build with cmake"

(cd ${ParaMonte_BLD_DIR} && \
make \
)
verify $? "build with make"

(cd ${ParaMonte_BLD_DIR} && \
make install \
)
verify $? "installation"

LD_LIBRARY_PATH=${ParaMonte_BLD_DIR}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH
if [ "${ParaMonteTest_RUN_ENABLED}" = "true" ]; then
    if [ "${MPI_ENABLED}" = "true" ]; then
        (cd ${ParaMonte_BLD_DIR}/test/bin && \
        "${MPIEXEC_PATH}" -np ${FOR_COARRAY_NUM_IMAGES} ./testParaMonte \
        )
    else
        if [ "${PMCS}" = "gnu" ]; then
            if [[ -f "${SETUP_FILE_PATH}" ]]; then
                cp ${SETUP_FILE_PATH} ${ParaMonte_BLD_DIR}/test/bin/
            fi
            if [ "${CAF_ENABLED}" = "true" ]; then
                (cd ${ParaMonte_BLD_DIR}/test/bin && \
                source ${SETUP_FILE_PATH} &&\
                cafrun -np ${FOR_COARRAY_NUM_IMAGES} ./testParaMonte \
                )
            else
                (cd ${ParaMonte_BLD_DIR}/test/bin && \
                ./testParaMonte \
                )
            fi
        fi
        if [ "${PMCS}" = "intel" ]; then
            if [ "${CAF_ENABLED}" = "true" ]; then
                (export FOR_COARRAY_NUM_IMAGES && \
                cd ${ParaMonte_BLD_DIR}/test/bin && \
                ./testParaMonte \
                )
            else
                (cd ${ParaMonte_BLD_DIR}/test/bin && \
                ./testParaMonte \
                )
            fi
        fi
    fi
    verify $? "test run"
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
PMLIB_FULL_PATH="$(ls ${ParaMonte_LIB_DIR}/libparamonte_${LTYPE}_${MEMORY_ALLOCATION}_${BTYPE}_*${PMLIB_EXT} | sort -V | tail -n1)"
PMLIB_FULL_NAME=${PMLIB_FULL_PATH##*/}
PMLIB_BASE_NAME=${PMLIB_FULL_NAME%.*}

echo >&2 
echo >&2 "-- ParaMonte - ParaMonte installed-library full path: ${PMLIB_FULL_PATH}"
echo >&2 "-- ParaMonte - ParaMonte installed-library full name: ${PMLIB_FULL_NAME}"
echo >&2 "-- ParaMonte - ParaMonte installed-library base name: ${PMLIB_BASE_NAME}"
echo >&2 

####################################################################################################################################
# build ParaMonte Python test
####################################################################################################################################

if [ "${LTYPE}" = "dynamic" ] && [ "${CFI_ENABLED}" = "true" ]; then

    echo >&2 
    echo >&2 "-- ParaMontePythonTest - building ParaMonte Python test..."
    ParaMontePythonTest_BLD_DIR=${ParaMonte_BLD_DIR}/test/Python
    export ParaMontePythonTest_BLD_DIR
    
    PYTHON_TEST_FILENAME=testParaMonte_${BTYPE}
    if [ "${CAF_ENABLED}" = "true" ]; then PYTHON_TEST_FILENAME=${PYTHON_TEST_FILENAME}_${CAFTYPE}; fi
    if [ "${MPI_ENABLED}" = "true" ]; then PYTHON_TEST_FILENAME=${PYTHON_TEST_FILENAME}_mpi; fi
    if [ "${OMP_ENABLED}" = "true" ]; then PYTHON_TEST_FILENAME=${PYTHON_TEST_FILENAME}_omp; fi
    PYTHON_TEST_FILENAME=${PYTHON_TEST_FILENAME}.py

    if [[ -d "$ParaMontePythonTest_BLD_DIR" ]]
    then
        echo >&2 "-- ParaMontePythonTest - ${ParaMontePythonTest_BLD_DIR} already exists. skipping..."
    else
        echo >&2 "-- ParaMontePythonTest - generating Python files directory: ${ParaMontePythonTest_BLD_DIR}"
        mkdir "${ParaMontePythonTest_BLD_DIR}"
    fi

    # copy necessary ParaMonte Python library files in Python's directory

    echo >&2 "-- ParaMontePythonTest - copying ParaMonte library files to the Python directory"
    echo >&2 "-- ParaMontePythonTest - from: ${ParaMonteInterfacePython_SRC_DIR}/paramonte"
    echo >&2 "-- ParaMontePythonTest -   to: ${ParaMontePythonTest_BLD_DIR}/paramonte/"
    cp -R ${ParaMonteInterfacePython_SRC_DIR}/paramonte ${ParaMontePythonTest_BLD_DIR}/
    echo >&2 

    # copy necessary ParaMonte Python DLL files in Python's directory

    echo >&2 "-- ParaMontePythonTest - copying ParaMonte shared library files to the Python directory"
    echo >&2 "-- ParaMontePythonTest - from: ${ParaMonte_LIB_DIR}/${PMLIB_FULL_NAME}"
    echo >&2 "-- ParaMontePythonTest -   to: ${ParaMontePythonTest_BLD_DIR}/paramonte/"
    cp ${ParaMonte_LIB_DIR}/${PMLIB_FULL_NAME} ${ParaMontePythonTest_BLD_DIR}/paramonte/
    echo >&2 

    # copy necessary ParaMonte Python library files in Python's directory

    echo >&2 "-- ParaMontePythonTest - copying ParaMonte library test files to the Python directory"
    echo >&2 "-- ParaMontePythonTest - from: ${ParaMontePythonTest_SRC_DIR}/${PYTHON_TEST_FILENAME}"
    echo >&2 "-- ParaMontePythonTest -   to: ${ParaMontePythonTest_BLD_DIR}/"
    cp ${ParaMontePythonTest_SRC_DIR}/${PYTHON_TEST_FILENAME} ${ParaMontePythonTest_BLD_DIR}/
    echo >&2 

    # copy necessary input files in Python's directory

    echo >&2 "-- ParaMontePythonTest - copying input files to the Python directory"
    echo >&2 "-- ParaMontePythonTest - from: ${ParaMonteTest_SRC_DIR}/input"
    echo >&2 "-- ParaMontePythonTest -   to: ${ParaMontePythonTest_BLD_DIR}/input/"
    cp -R ${ParaMonteTest_SRC_DIR}/input ${ParaMontePythonTest_BLD_DIR}/
    echo >&2 

fi

####################################################################################################################################
# build ParaMonte example
####################################################################################################################################

source ./example/buildParaMonteExample.sh

####################################################################################################################################
# copy ParaMonte binary / library files to the bin directory
####################################################################################################################################

echo >&2 
echo >&2 "-- ParaMonte - copying ParaMonte binary/library files to the bin directory..."
ParaMonte_BIN_DIR=${ParaMonte_ROOT_DIR}/bin
export ParaMonte_BIN_DIR
mkdir -p ${ParaMonte_BIN_DIR}


####################################################################################################################################
# mission accomplished
####################################################################################################################################

echo >&2
echo >&2 "-- ParaMonte - ParaMonte binary/library directory: ${ParaMonte_BIN_DIR}"
echo >&2 "-- ParaMonte - ParaMonte build directory: ${ParaMonte_BLD_DIR}"
echo >&2
echo >&2 "-- ParaMonte - mission accomplished"
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
