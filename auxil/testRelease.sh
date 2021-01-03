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
####       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
####
####################################################################################################################################
####################################################################################################################################

BUILD_NAME="ReleaseTest"
workingDir="$(pwd)"

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
#### determine the platform
####################################################################################################################################

UNAME_PLATFORM="$(uname -s)"
case "${UNAME_PLATFORM}" in
    Linux*)     PLATFORM=linux;;
    Darwin*)    PLATFORM=darwin;;
    CYGWIN*)    PLATFORM=windows;;
    MINGW*)     PLATFORM=windows;;
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
#### fetch options
####################################################################################################################################

LANG_LIST="c cpp fortran python matlab"
BTYPE_LIST="debug release"
LTYPE_LIST="shared"
PARALLELISM_LIST="none impi mpich openmpi"
PMCS_LIST="gnu intel"
MEMORY_LIST="heap"
pmReleaseLink="https://github.com/cdslaborg/paramonte/releases/latest/download"
unset pmKernelVersionNumber

while [ "$1" != "" ]; do
    case $1 in
        -v | --version )        shift
                                pmKernelVersionNumber="$1"
                                pmReleaseLink="https://github.com/cdslaborg/paramonte/releases/download/v${pmKernelVersionNumber}"
                                ;;
        -L | --lang )           shift
                                LANG_LIST="$1"
                                ;;
        -s | --compiler_suite ) shift
                                PMCS_LIST="$1"
                                ;;
        -b | --build )          shift
                                BTYPE_LIST="$1"
                                ;;
        -l | --lib )            shift
                                LTYPE_LIST="$1"
                                ;;
        -p | --par )            shift
                                PARALLELISM_LIST="$1"
                                ;;
        -m | --mem )            shift
                                MEMORY_LIST="$1"
                                ;;
        -t | --test )           shift
                                TTYPE="$1"
                                ;;
        -x | --exam_enabled )   shift
                                ParaMonteExample_RUN_ENABLED="$1"
                                ;;
        -D | --deploy )         deploy_flag="--deploy"
                                ;;
        -f | --fortran )        shift
                                Fortran_COMPILER_PATH="$1"
                                ;;
        -M | --mpiexec )        shift
                                MPIEXEC_PATH="$1"
                                ;;
        -F | --fresh )          fresh_flag="--fresh"
                                ;;
        -O | --local )          local_flag="--local"
                                ;;
        -d | --dryrun )         dryrun_flag="--dryrun"
                                ;;
        -y | --yes-to-all )     yes_to_all_flag="--yes-to-all"
                                ;;
        -B | --bootstrap )      gcc_bootstrap_flag="--bootstrap"
                                ;;
        -a | --matdram )        MatDRAM_ENABLED="true"
                                ;;
        -c | --codecov )        codecov_flag="--codecov"
                                ;;
        -n | --nproc )          shift
                                FOR_COARRAY_NUM_IMAGES="$1"
                                ;;
        -h | --help )           usage
                                echo >&2 ""
                                echo >&2 ""
                                exit
                                ;;
        * )                     usage
                                echo >&2 ""
                                echo >&2 "-- ParaMonte - FATAL: The specified flag $1 does not exist."
                                echo >&2 ""
                                echo >&2 "-- ParaMonte - gracefully exiting."
                                echo >&2 ""
                                echo >&2 ""
                                exit 1
    esac
    shift
done

####################################################################################################################################
#### download libraries and run tests
####################################################################################################################################

for PMCS in $PMCS_LIST; do
    for LANG in $LANG_LIST; do
        for BTYPE in $BTYPE_LIST; do
            for LTYPE in $LTYPE_LIST; do
                for MEMORY in $MEMORY_LIST; do
                    for PARALLELISM in $PARALLELISM_LIST; do

                        if [ "${LANG}" = "matlab" ] || [ "${LANG}" = "python" ]; then
                            BENABLED=false
                        else
                            BENABLED=true
                        fi

                        if [ "${PARALLELISM}" = "none" ]; then
                            parSuffix=""
                        else
                            parSuffix="_${PARALLELISM}"
                            if ([ "${PMCS}" = "intel" ] && ([ "${PARALLELISM}" = "openmpi" ] || [ "${PARALLELISM}" = "mpich" ])) \
                            || ([ "${PMCS}" = "gnu" ] && [ "${PARALLELISM}" = "impi" ]) \
                            then
                                BENABLED=false
                            fi
                        fi

                        if [ "${BENABLED}" = "true" ]; then

                            pmLibName="libparamonte_${LANG}_${PLATFORM}_${ARCHITECTURE}_${PMCS}_${BTYPE}_${LTYPE}_${MEMORY}${parSuffix}"
                            compressedFileExt=".tar.gz"
                            untar=(tar xvzf "${pmLibName}${compressedFileExt}")

                            if [ "${isMacOS}" = "true" ]; then
                                fetch="curl -OL"
                            else
                                fetch="wget"
                            fi
                            fetch+=" ${pmReleaseLink}/${pmLibName}${compressedFileExt}"


                            tempDir=$(mktemp -d "${TMPDIR:-/tmp}/cversion.XXXXXXXXX")
                            echo >&2 "-- ${BUILD_NAME} - changing directory to: ${tempDir}"

                            cd "${tempDir}" \
                            && \
                            echo >&2 && echo >&2 "-- ${BUILD_NAME} - ${fetch}" && echo >&2 && $(${fetch}) \
                            && \
                            echo >&2 && echo >&2 "-- ${BUILD_NAME} - ${untar[@]}" && echo >&2 && "${untar[@]}" \
                            && \
                            cd ${pmLibName} \
                            && \
                            ./build.sh && ./run.sh \
                            || {
                                echo >&2
                                echo >&2 "-- ${BUILD_NAME} - test FAILED for ${pmLibName}."
                                echo >&2 "-- ${BUILD_NAME} - gracefully exiting."
                                echo >&2
                                cd "${workingDir}"
                                exit 1
                            }

                        fi

                    done
                done
            done
        done
    done
done

cd "${workingDir}"