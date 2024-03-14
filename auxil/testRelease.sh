#!/bin/bash
####################################################################################################################################
####################################################################################################################################
####                                                                                                                            ####
####    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           ####
####                                                                                                                            ####
####    Copyright (C) 2012-present, The Computational Data Science Lab                                                          ####
####                                                                                                                            ####
####    This file is part of the ParaMonte library.                                                                             ####
####                                                                                                                            ####
####    LICENSE                                                                                                                 ####
####                                                                                                                            ####
####       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          ####
####                                                                                                                            ####
####################################################################################################################################
####################################################################################################################################

build_name="ReleaseTest"
workingDir="$(pwd)"
srcFileDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

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
    echo >&2 "-- ${build_name} - FATAL: Build failed. unrecognized platform - ${PLATFORM}"
    echo >&2 "-- ${build_name} - FATAL: The supported platforms include: Linux, Darwin, CYGWIN, MINGW"
    echo >&2 "-- ${build_name} - FATAL: The ParaMonte build has been only tested on Linux and Darwin platforms."
    echo >&2
    echo >&2 "-- ${build_name} - gracefully exiting."
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

FOR_COARRAY_NUM_IMAGES=2
LANG_LIST="c cpp fortran python matlab"
BTYPE_LIST="debug release"
LTYPE_LIST="shared"
PARALLELISM_LIST="none impi mpich openmpi"
csid_LIST="gnu intel"
MEMORY_LIST="heap"
pmReleaseLink="https://github.com/cdslaborg/paramonte/releases/latest/download"
unset pmVersionNumber

while [ "$1" != "" ]; do
    case $1 in
        -v | --version )        shift
                                pmVersionNumber="$1"
                                pmReleaseLink="https://github.com/cdslaborg/paramonte/releases/download/v${pmVersionNumber}"
                                ;;
        -L | --lang )           shift
                                LANG_LIST="$1"
                                ;;
        -s | --compiler_suite ) shift
                                csid_LIST="$1"
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
        -D | --deps )           flag_deps="--deps"
                                ;;
        -f | --fortran )        shift
                                Fortran_COMPILER_PATH="$1"
                                ;;
        -M | --mpiexec )        shift
                                MPIEXEC_PATH="$1"
                                ;;
        -F | --fresh )          flag_fresh="--fresh"
                                ;;
        -O | --local )          local_flag="--local"
                                ;;
        -d | --dryrun )         dryrun_flag="--dryrun"
                                ;;
        -y | --yes-to-all )     yes_to_all_flag="--yes-to-all"
                                ;;
        -B | --bootstrap )      gcc_bootstrap_flag="--bootstrap"
                                ;;
        -c | --codecov )        flag_codecov="--codecov"
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

for csid in $csid_LIST; do
    for LANG in $LANG_LIST; do
        for BTYPE in $BTYPE_LIST; do
            for LTYPE in $LTYPE_LIST; do
                for MEMORY in $MEMORY_LIST; do
                    for PARALLELISM in $PARALLELISM_LIST; do

                        BENABLED=true

                        if [ "${isWindows}" = "true" ] && ([ "${csid}" = "gnu" ] || [ "${PARALLELISM}" = "mpich" ] || [ "${PARALLELISM}" = "openmpi" ]); then
                            BENABLED=false
                        fi

                        if [ "${LANG}" = "matlab" ] || [ "${LANG}" = "python" ]; then
                            BENABLED=false
                        fi

                        if [ "${PARALLELISM}" = "none" ]; then
                            nproc_flag=""
                            parSuffix=""
                        else
                            parSuffix="_${PARALLELISM}"
                            nproc_flag="--nproc ${FOR_COARRAY_NUM_IMAGES}"
                            if ([ "${csid}" = "intel" ] && ([ "${PARALLELISM}" = "openmpi" ] || [ "${PARALLELISM}" = "mpich" ])) \
                            || ([ "${csid}" = "intel" ] && [ "${PARALLELISM}" = "impi" ] && [ "${isMacOS}" = "true" ]) \
                            || ([ "${csid}" = "gnu" ] && [ "${PARALLELISM}" = "impi" ]) \
                            then
                                BENABLED=false
                            fi
                        fi

                        if [ "${BENABLED}" = "true" ]; then

                            tempDir=$(mktemp -d "${TMPDIR:-/tmp}/cversion.XXXXXXXXX")
                            echo >&2 "-- ${build_name} - changing directory to: ${tempDir}"
                            cd "${tempDir}"

                            pmLibName="libparamonte_${LANG}_${PLATFORM}_${ARCHITECTURE}_${csid}_${BTYPE}_${LTYPE}_${MEMORY}${parSuffix}"
                            compressedFileExt=".tar.gz"
                            untar=(tar xvzf "${pmLibName}${compressedFileExt}")

                            compressedFileName="${pmLibName}"
                            if [ "${isWindows}" = "true" ]; then
                                fetch="${srcFileDir}/wget"
                                compressedFileName+=".zip"
                                untar=("${srcFileDir}/7z" -Y e "${tempDir}/${compressedFileName}" -o"${tempDir}/${pmLibName}")
                            else
                                compressedFileName+=".tar.gz"
                                untar=(tar xvzf "${tempDir}/${compressedFileName}")
                                if [ "${isMacOS}" = "true" ]; then
                                    fetch="curl -OL"
                                else
                                    fetch="wget"
                                fi
                            fi
                            fetch+=" ${pmReleaseLink}/${compressedFileName}"

                            echo >&2 && echo >&2 "-- ${build_name} - ${fetch}" && echo >&2 && $(${fetch}) \
                            && \
                            echo >&2 && echo >&2 "-- ${build_name} - ${untar[@]}" && echo >&2 && "${untar[@]}" \
                            && \
                            cd ${pmLibName} \
                            && \
                            ./build.sh && ./run.sh ${nproc_flag} \
                            || {
                                echo >&2
                                echo >&2 "-- ${build_name} - test FAILED for ${pmLibName}."
                                echo >&2 "-- ${build_name} - gracefully exiting."
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

####################################################################################################################################
#### test matlab
####################################################################################################################################

for LANG in $LANG_LIST; do

    if [ "${LANG}" = "matlab" ] || [ "${LANG}" = "python" ]; then

        unset plexe
        if [ "${LANG}" = "matlab" ]; then
            plexe=matlab
            cmd=(${plexe} -batch main)
        elif [ "${LANG}" = "python" ]; then
            for plexe in python python3; do
                if command -v ${plexe} >/dev/null 2>&1; then
                    pythonPath="$(command -v ${plexe})"
                    if ! [[ "${pythonPath}" =~ .*"Microsoft".* ]]; then
                        pythonVersion="$(${plexe} --version || echo echo >&2 "error occurred while capturing ${plexe} version. skipping...")"
                        if [[ "${pythonVersion}" =~ .*"3.".* ]]; then
                            break
                        fi
                    fi
                    unset pythonPath
                fi
            done
            if ! [ -z ${pythonPath+x} ]; then
                cmd=(${plexe} main.py)
            fi
        fi

        if ! [ -z ${plexe+x} ]; then

            tempDir=$(mktemp -d "${TMPDIR:-/tmp}/cversion.XXXXXXXXX")
            echo >&2 "-- ${build_name} - changing directory to: ${tempDir}"
            cd "${tempDir}"

            pmLibName="libparamonte_${LANG}_${PLATFORM}_${ARCHITECTURE}"

            compressedFileName="${pmLibName}"
            if [ "${isWindows}" = "true" ]; then
                fetch="${srcFileDir}/wget"
                compressedFileName+=".zip"
                untar=("${srcFileDir}/7z" -Y e "${tempDir}/${compressedFileName}" -o"${tempDir}/${pmLibName}")
            else
                compressedFileName+=".tar.gz"
                untar=(tar xvzf "${tempDir}/${compressedFileName}")
                if [ "${isMacOS}" = "true" ]; then
                    fetch="curl -OL"
                else
                    fetch="wget"
                fi
            fi
            fetch+=" ${pmReleaseLink}/${compressedFileName}"

            echo >&2 && echo >&2 "-- ${build_name} - ${fetch}" && echo >&2 && $(${fetch}) \
            && \
            echo >&2 && echo >&2 "-- ${build_name} - ${untar[@]}" && echo >&2 && "${untar[@]}" \
            && \
            ls "${tempDir}/" && cd "${tempDir}/${pmLibName}" \
            && \
            echo >&2 && echo >&2 "-- ${build_name} - ${cmd[@]}" && echo >&2 && "${cmd[@]}" \
            && \
            {
            if [[ "${PARALLELISM}" =~ .*"impi".* ]] || [[ "${PARALLELISM}" =~ .*"mpich".* ]] || [[ "${PARALLELISM}" =~ .*"openmpi".* ]]; then
                cmd=(mpiexec)
                if [ "${PLATFORM}" = "windows" ]; then
                    cmd+=(-localonly)
                fi
                cmd+=(-n 2 ${plexe})
                if [ "${LANG}" = "matlab" ]; then
                    cmd+=(-batch main_mpi)
                elif [ "${LANG}" = "python" ]; then
                    cmd+=(main_mpi.py)
                fi
            fi
            echo >&2 && echo >&2 "-- ${build_name} - ${cmd[@]}" && echo >&2 && "${cmd[@]}"
            } || {
                echo >&2
                echo >&2 "-- ${build_name} - test FAILED for ${pmLibName}."
                echo >&2 "-- ${build_name} - gracefully exiting."
                echo >&2
                cd "${workingDir}"
                exit 1
            }

        fi

    fi

done

####################################################################################################################################

cd "${workingDir}"
