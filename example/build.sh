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

# In this script, the priority is always given to the Intel over GNU compilers if they are detected.

FILE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

usage()
{
cat << EndOfMessage

This is a simple standalone Bash script for building C/Fortran applications that use ParaMonte library on Unix-like OS.
The compiler name and flags are automatically inferred from the library name. Note this requires the existence of the 
required compiler and libraries on your system. These libraries should be automatically installed on your system 
when you build ParaMonte via the provided build-scripts at the root of the project.
You can override the compiler options and flags by providing the following optional arguments to this build-script.
However, the onus on you to ensure compiler/library availability and compatibility with the ParaMonte library.

    usage:

        build.sh -c <compiler name/path: gcc/g++/gfortran/icc/icpc/ifort/clang> -f <compiler_flags>

    example:

        build.sh -c /opt/apps/gcc/8.3.0/bin/gfortran -n 3

    flag definitions:

        -c | --compiler         : the compiler name to be used for building the ParaMonte example.
                                : possible values: gcc/g++/gfortran/icc/icpc/ifort/clang
                                : possible values: mpicc/mpic++/mpifort/mpiicc/mpiicpc/mpiifort/caf
                                : WARNING: if you provide the compiler name the onus is on YOU
                                : WARNING: to ensure that the compiler exists on your system
                                : WARNING: and that a proper compiler has been chosen.
                                : the default will be determined based on the runtime environment.
        -f | --compiler_flags   : compiler flags. If specified, it will override all default compiler flags.
                                : If not provided, appropriate values will be set for each flag.
                                : If multiple space-delimited flags are passed, enclose all with "".
        -n | --nproc            : the default number of processes (coarray images) on which the ParaMonte 
                                : examples/tests (if any) will be run: positive integer
                                : If not provided, the default is 3.
        -h | --help             : help with the script usage

NOTE: ALL FLAGS ARE OPTIONAL. If not provided, appropriate values will be set for each missing flag.
NOTE: 
NOTE: If the compiler flags, -f | --compiler_flags, option is specified, it will override all default compiler flags.
NOTE: If not provided, appropriate values will be set for each flag.
NOTE: 
NOTE: Upon finishing the build, the script will generate another Bash script named run.sh in 
NOTE: the same directory, which can be used to run the executable. Usage:
NOTE: 
NOTE:     ./run.sh
NOTE: 
NOTE: or, 
NOTE: 
NOTE:     source ./run.sh

EndOfMessage
}

unset FOR_COARRAY_NUM_IMAGES
unset USER_SELECTED_COMPILER
unset USER_SELECTED_COMPILER_FLAGS

while [ "$1" != "" ]; do
    case $1 in
        -c | --compiler )       shift
                                USER_SELECTED_COMPILER=$1
                                ;;
        -f | --compiler_flags ) shift
                                USER_SELECTED_COMPILER_FLAGS=$1
                                ;;
        -n | --nproc )          shift
                                FOR_COARRAY_NUM_IMAGES=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     echo >&2 "-- ParaMonteExampleBuildScript - FATAL: the input flag is not recognized: $1"
                                usage
                                exit 1
    esac
    shift
done

####################################################################################################################################
# get the platform name. required to define __WIN64__ for binding GFortran applications to Intel-prebuilt ParaMonte libraries on Windows.
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
    echo >&2 "-- ${BUILD_NAME} - FATAL: build failed. unrecognized platform - ${PLATFORM}"
    echo >&2 "-- ${BUILD_NAME} - supported platforms include: Linux, Darwin, CYGWIN, MINGW"
    echo >&2 "-- ${BUILD_NAME} - ParaMonte build has been only tested on Linux and Darwin platforms."
    echo >&2
    echo >&2 "-- ${BUILD_NAME} - gracefully exiting."
    echo >&2
    exit 1
else
    export PLATFORM
fi

pmLibExtList=""
isMinGW=false
isLinux=false
isDarwin=false
isCygwin=false

if [[ "${PLATFORM}" =~ .*"darwin".* ]]; then
    isDarwin=true
    OSNAME="macOS"
    pmLibExtList=".dylib .a"
elif [[ "${PLATFORM}" =~ .*"linux".* ]]; then
    isLinux=true
    OSNAME="Linux"
    pmLibExtList=".so .a"
elif [[ "${PLATFORM}" =~ .*"mingw".* ]]; then
    isMinGW=true
    OSNAME="MinGW"
    pmLibExtList=".dll .lib"
elif [[ "${PLATFORM}" =~ .*"cygwin".* ]]; then
    isCygwin=true
    OSNAME="Cygwin"
    pmLibExtList=".dll .lib"
fi

####################################################################################################################################
# build ParaMonteExample objects and executable
####################################################################################################################################

PM_EXAM_EXE_NAME="main.exe"
export PM_EXAM_EXE_NAME

if [ -z ${FOR_COARRAY_NUM_IMAGES+x} ]; then
    FOR_COARRAY_NUM_IMAGES=3
fi
export FOR_COARRAY_NUM_IMAGES

for pmLibExt in ${pmLibExtList}; do
    pmLibFullPath="$(ls -dp ${FILE_DIR}/*libparamonte_*${pmLibExt} | sort -V | tail -n1)"
    if [ -f "${pmLibFullPath}" ]; then
        break
    else
        unset pmLibFullPath
    fi
done

if [ -f "${pmLibFullPath}" ]; then
    pmLibFullName=${pmLibFullPath##*/}
    pmLibBaseName=${pmLibFullName%.*}
else
    echo >&2
    echo >&2 "-- ParaMonteExample - FATAL: pmLibFullPath=${pmLibFullPath}"
    echo >&2 "-- ParaMonteExample - FATAL: The ParaMonte library file could not be found in the current directory."
    echo >&2 "-- ParaMonteExample - FATAL: This is the file that is prefixed with libparamonte_ and suffixed with"
    echo >&2 "-- ParaMonteExample - FATAL: the library file extension which can be either: dll, so, dylib, lib, a"
    echo >&2 "-- ParaMonteExample - FATAL: Ensure the ParaMonte library file exists in the current directory, then rerun the script."
    echo >&2 "-- ParaMonteExample - FATAL: If the problem persists, please report this issue at:"
    echo >&2 "-- ParaMonteExample - FATAL: "
    echo >&2 "-- ParaMonteExample - FATAL:     https://github.com/cdslaborg/paramonte/issues"
    echo >&2 "-- ParaMonteExample - FATAL: "
    echo >&2 "-- ParaMonteExample - FATAL: Gracefully exiting..."
    echo >&2
    exit 1
fi

####################################################################################################################################
# determine the existing GNU and Intel compilers, generate a list of potential compilers, infer the example's language and src files
####################################################################################################################################

unset compilerList
iccIntelCompilerDetected=false
icpcIntelCompilerDetected=false
ifortIntelCompilerDetected=false

if [ "$isMinGW" = "true" ] || [ "$isCygwin" = "true" ]; then
    cIntelCompilerName="icl"
    cppIntelCompilerName="icl"
    cIntelMPIWrapperName="mpicl"
    cppIntelMPIWrapperName="mpicl"
else
    cIntelCompilerName="icc"
    cppIntelCompilerName="icpc"
    cIntelMPIWrapperName="mpiicc"
    cppIntelMPIWrapperName="mpiicpc"
fi

if [[ "${pmLibFullName}" =~ .*"_c_".* ]]; then

    pmSrcExt=c
    pmExamLang=C

    clist="gcc"

    if command -v ${cIntelCompilerName} >/dev/null 2>&1; then
        iccIntelCompilerDetected=true
        declare -a compilerList=("${cIntelCompilerName}" "${clist}")
    else
        declare -a compilerList=("${clist}")
    fi

elif [[ "${pmLibFullName}" =~ .*"_cpp_".* ]]; then

    pmSrcExt=cpp
    pmExamLang=C++

    cpattern="g++"
    clist=$(( IFS=:; for p in $PATH; do unset lsout; lsout=$(ls -dm "$p"/${cpattern}*); if ! [[ -z "${lsout// }" ]]; then echo "${lsout}, "; fi; done ) 2>/dev/null)

    if command -v ${cppIntelCompilerName} >/dev/null 2>&1; then
        icpcIntelCompilerDetected=true
        declare -a compilerList=("${cppIntelCompilerName}" "${clist}")
    else
        declare -a compilerList=("${clist}")
    fi

elif [[ "${pmLibFullName}" =~ .*"_fortran_".* ]]; then

    pmSrcExt=f90
    pmExamLang=Fortran
    pmSrcFiles="paramonte.${pmSrcExt}"

    cpattern="gfortran"
    clist=$(( IFS=:; for p in $PATH; do unset lsout; lsout=$(ls -dm "$p"/${cpattern}*); if ! [[ -z "${lsout// }" ]]; then echo "${lsout}, "; fi; done ) 2>/dev/null)

    if command -v ifort >/dev/null 2>&1; then
        ifortIntelCompilerDetected=true
        declare -a compilerList=("ifort" "${clist}")
    else
        declare -a compilerList=("${clist}")
    fi

else

    echo >&2
    echo >&2 "-- ParaMonteExample - FATAL: The ParaMonte library example's programming language is unrecognizable."
    echo >&2 "-- ParaMonteExample - FATAL: The supported languages are: C, C++, Fortran."
    echo >&2 "-- ParaMonteExample - FATAL: Ensure the ParaMonte library file exists in the current directory, then rerun the script."
    echo >&2 "-- ParaMonteExample - FATAL: If the problem persists, please report this issue at:"
    echo >&2 "-- ParaMonteExample - FATAL: "
    echo >&2 "-- ParaMonteExample - FATAL:     https://github.com/cdslaborg/paramonte/issues"
    echo >&2 "-- ParaMonteExample - FATAL: "
    echo >&2 "-- ParaMonteExample - FATAL: Gracefully exiting..."
    echo >&2
    exit 1

fi

pmSrcFiles="${pmSrcFiles} logfunc.${pmSrcExt} main.${pmSrcExt}"; export pmSrcFiles

echo >&2
echo >&2 "-- ParaMonteExample${pmExamLang} - The ParaMonte library's full path: ${pmLibFullPath}"
echo >&2 "-- ParaMonteExample${pmExamLang} - The ParaMonte library's full name: ${pmLibFullName}"
echo >&2 "-- ParaMonteExample${pmExamLang} - The ParaMonte library's base name: ${pmLibBaseName}"
echo >&2

####################################################################################################################################
# infer the ParaMonte compiler suite
####################################################################################################################################

unset pmCompilerSuite

if [[ "${pmLibFullName}" =~ .*"_intel_".* ]]; then
    pmCompilerSuite=intel
elif [[ "${pmLibFullName}" =~ .*"_gnu_".* ]]; then
    pmCompilerSuite=gnu
else
    echo >&2
    echo >&2 "-- ParaMonteExample${pmExamLang} - WARNING: The compiler suite with which the ParaMonte library was build is unrecognizable."
    echo >&2 "-- ParaMonteExample${pmExamLang} - WARNING: There is a strong likelihood that this build will fail."
    echo >&2
fi

####################################################################################################################################
# build setup
####################################################################################################################################

unset pmBuildType
unset GNU_C_COMPILER_FLAGS
unset INTEL_C_COMPILER_FLAGS

GNU_Fortran_COMPILER_FLAGS="-std=gnu -cpp"

if [ "$isMinGW" = "true" ] || [ "$isCygwin" = "true" ]; then

    GNU_Fortran_COMPILER_FLAGS="${GNU_Fortran_COMPILER_FLAGS} -D__WIN64__"

    INTEL_C_COMPILER_FLAGS_DEBUG="/debug:full /Zi /Od /Wall /traceback /Qcheck-pointers:rw"
    INTEL_C_COMPILER_FLAGS_TESTING="/Od"
    INTEL_C_COMPILER_FLAGS_RELEASE="/O3 /Qip /Qipo /Qunroll /Qunroll-aggressive /Ob2 /Qparallel /Qinline-dllimport"

    INTEL_Fortran_COMPILER_FLAGS_DEBUG="/debug:full /Zi /CB /Od"
    INTEL_Fortran_COMPILER_FLAGS_TESTING="/Od"
    INTEL_Fortran_COMPILER_FLAGS_RELEASE="/O3 /Qip /Qipo /Qunroll /Qunroll-aggressive"
    INTEL_Fortran_COMPILER_FLAGS="/standard-semantics /fpp /nologo /F0x1000000000"

else

    INTEL_Fortran_COMPILER_FLAGS_DEBUG="-O0 -debug full"
    INTEL_Fortran_COMPILER_FLAGS_TESTING="-O0"
    INTEL_Fortran_COMPILER_FLAGS_RELEASE="-O3 -ip -ipo -unroll -unroll-aggressive -finline-functions"
    INTEL_Fortran_COMPILER_FLAGS="-standard-semantics -fpp"

    INTEL_C_COMPILER_FLAGS_DEBUG="-O0 -debug full"
    INTEL_C_COMPILER_FLAGS_TESTING="-O0"
    INTEL_C_COMPILER_FLAGS_RELEASE="-O3"
    if [ "${isDarwin}" = "true" ]; then
        INTEL_C_COMPILER_FLAGS="-no-multibyte-chars"
    fi

fi

if [[ "${pmLibFullName}" =~ .*"release".* ]]; then

    pmBuildType=release

    GNU_C_COMPILER_FLAGS="${GNU_C_COMPILER_FLAGS} -O3"
    GNU_Fortran_COMPILER_FLAGS="${GNU_Fortran_COMPILER_FLAGS} -O3 -funroll-loops -finline-functions -ftree-vectorize"

    INTEL_C_COMPILER_FLAGS="${INTEL_C_COMPILER_FLAGS} ${INTEL_C_COMPILER_FLAGS_RELEASE}"
    INTEL_Fortran_COMPILER_FLAGS="${INTEL_Fortran_COMPILER_FLAGS} ${INTEL_Fortran_COMPILER_FLAGS_RELEASE}"

fi

if [[ "${pmLibFullName}" =~ .*"testing".* ]]; then

    pmBuildType=testing

    GNU_C_COMPILER_FLAGS="${GNU_C_COMPILER_FLAGS} -O0"
    GNU_Fortran_COMPILER_FLAGS="${GNU_Fortran_COMPILER_FLAGS} -O0"

    INTEL_C_COMPILER_FLAGS="${INTEL_C_COMPILER_FLAGS} ${INTEL_C_COMPILER_FLAGS_TESTING}"
    INTEL_Fortran_COMPILER_FLAGS="${INTEL_Fortran_COMPILER_FLAGS} ${INTEL_Fortran_COMPILER_FLAGS_TESTING}"

fi

if [[ "${pmLibFullName}" =~ .*"debug".* ]]; then

    pmBuildType=debug

    GNU_C_COMPILER_FLAGS="${GNU_C_COMPILER_FLAGS} -O0 -g"
    GNU_Fortran_COMPILER_FLAGS="${GNU_Fortran_COMPILER_FLAGS} -O0 -g"

    INTEL_C_COMPILER_FLAGS="${INTEL_C_COMPILER_FLAGS} ${INTEL_C_COMPILER_FLAGS_DEBUG}"
    INTEL_Fortran_COMPILER_FLAGS="${INTEL_Fortran_COMPILER_FLAGS} ${INTEL_Fortran_COMPILER_FLAGS_DEBUG}"

fi

if [ "${pmExamLang}" = "C++" ]; then
    GNU_C_COMPILER_FLAGS="${GNU_C_COMPILER_FLAGS} -std=gnu++11"
fi

echo >&2 "-- ParaMonteExample${pmExamLang} - ParaMonte build type: ${pmBuildType}"

####################################################################################################################################
# library type
####################################################################################################################################

if [[ "${pmLibFullName}" =~ .*"dynamic".* ]]; then
    pmLibType=dynamic
else
    pmLibType=static
fi
echo >&2 "-- ParaMonteExample${pmExamLang} - The ParaMonte library type: ${pmLibType}"

####################################################################################################################################
# MPI
####################################################################################################################################

MPI_ENABLED=false
if [[ "${pmLibFullName}" =~ .*"_mpi".* ]] || [[ "${pmLibFullName}" =~ .*"impi".* ]] || [[ "${pmLibFullName}" =~ .*"mpich".* ]] || [[ "${pmLibFullName}" =~ .*"openmpi".* ]]; then
    MPI_ENABLED=true
fi
echo >&2 "-- ParaMonteExample${pmExamLang} - MPI_ENABLED: ${MPI_ENABLED}"

if [ "${MPI_ENABLED}" = "true" ] && [ "${pmLibType}" = "static" ]; then
    compilerListTemp=()
    if [ "${pmExamLang}" = "Fortran" ]; then
        if command -v mpiifort >/dev/null 2>&1; then
            compilerListTemp+=("mpiifort")
        fi
        if command -v mpiifort >/dev/null 2>&1; then
            compilerListTemp+=("mpifort")
        fi
    fi
    if [ "${pmExamLang}" = "C" ]; then
        if command -v ${cIntelMPIWrapperName} >/dev/null 2>&1; then
            compilerListTemp+=("${cIntelMPIWrapperName}")
        fi
        if command -v mpicc >/dev/null 2>&1; then
            compilerListTemp+=("mpicc")
        fi
    fi
    if [ "${pmExamLang}" = "C++" ]; then
        if command -v ${cppIntelMPIWrapperName} >/dev/null 2>&1; then
            compilerListTemp+=("${cppIntelMPIWrapperName}")
        fi
        if command -v mpicxx >/dev/null 2>&1; then
            compilerListTemp+=("mpicxx")
        fi
        if command -v "mpic++" >/dev/null 2>&1; then
            compilerListTemp+=("mpic++")
        fi
        if command -v "mpicc" >/dev/null 2>&1; then
            compilerListTemp+=("mpicc")
        fi
    fi
    if [ ${#compilerListTemp[@]} -eq 0 ]; then
        compilerList=("${compilerListTemp[@]}")
    fi
fi

####################################################################################################################################
# CAF
####################################################################################################################################

CAFTYPE=none
CAF_ENABLED=false
if [[ "${pmLibFullName}" =~ .*"cafsingle".* ]]; then
    #-fcoarray=single
    GNU_Fortran_COMPILER_FLAGS="${GNU_Fortran_COMPILER_FLAGS}"
    INTEL_Fortran_COMPILER_FLAGS="${INTEL_Fortran_COMPILER_FLAGS} -coarray=single"
    CAF_ENABLED=true
fi
if [[ "${pmLibFullName}" =~ .*"cafshared".* ]]; then
    #-fcoarray=shared
    GNU_Fortran_COMPILER_FLAGS="${GNU_Fortran_COMPILER_FLAGS}"
    INTEL_Fortran_COMPILER_FLAGS="${INTEL_Fortran_COMPILER_FLAGS} -coarray=shared"
    CAF_ENABLED=true
fi
if [[ "${pmLibFullName}" =~ .*"cafdistributed".* ]]; then
    #-fcoarray=distributed
    GNU_Fortran_COMPILER_FLAGS="${GNU_Fortran_COMPILER_FLAGS}"
    INTEL_Fortran_COMPILER_FLAGS="${INTEL_Fortran_COMPILER_FLAGS} -coarray=distributed"
    CAF_ENABLED=true
fi
echo >&2 "-- ParaMonteExample${pmExamLang} - CAFTYPE: ${CAFTYPE}"

if [ "${CAF_ENABLED}" = "true" ]; then

    if [[ "${pmLibFullName}" =~ .*"_gnu_".* ]]; then
        declare -a compilerList=("caf")
    fi

    if [[ "${pmLibFullName}" =~ .*"_intel_".* ]]; then
        declare -a compilerList=("ifort")
    fi

    if [ "$isMinGW" = "true" ] || [ "$isCygwin" = "true" ]; then
        echo >&2
        echo >&2 "-- ParaMonteExample${pmExamLang} - WARNING: Building Coarray Fortran applications via GNU Compilers is unsupported."
        echo >&2 "-- ParaMonteExample${pmExamLang} - WARNING: This application build will likely fail."
        echo >&2
        #exit 1
    fi

fi

####################################################################################################################################
# report compiler choice and compiler flags
####################################################################################################################################

echo >&2 "-- ParaMonteExample${pmExamLang} - ParaMonte library's compiler suite: ${pmCompilerSuite}"
echo >&2 "-- ParaMonteExample${pmExamLang} - inferred compiler choice(s):"
echo >&2 "-- ParaMonteExample${pmExamLang} - "
compilerListLen=${#compilerList[@]}
compilerListLenMinusOne="$(($compilerListLen-1))"
for i in $(seq 0 $compilerListLenMinusOne)
do
    csvCompilerList="${compilerList[$i]}"
    for COMPILER in $(echo ${csvCompilerList} | sed "s/,/ /g")
    do
        echo >&2 "-- ParaMonteExample${pmExamLang} -     ${COMPILER}"
    done
done
echo >&2

if [ -z ${USER_SELECTED_COMPILER+x} ]; then
    echo >&2 "-- ParaMonteExample${pmExamLang} - user-selected compiler/linker: none"
else
    declare -a compilerList=("${USER_SELECTED_COMPILER}")
    echo >&2 "-- ParaMonteExample${pmExamLang} - user-selected compiler/linker: ${USER_SELECTED_COMPILER}"
fi

if ! [ -z ${USER_SELECTED_COMPILER_FLAGS+x} ]; then
    echo >&2 "-- ParaMonteExample${pmExamLang} - user-selected compiler/linker flags: ${USER_SELECTED_COMPILER_FLAGS}"
fi

####################################################################################################################################
# build the example
####################################################################################################################################

if [ -f "${FILE_DIR}/setup.sh" ]; then
    source "${FILE_DIR}/setup.sh"
fi

if [ -z ${LD_LIBRARY_PATH+x} ]; then
    LD_LIBRARY_PATH="${FILE_DIR}"
else
    if [[ ":${LD_LIBRARY_PATH}:" != *":${FILE_DIR}:"* ]]; then
        LD_LIBRARY_PATH="${FILE_DIR}:${LD_LIBRARY_PATH}"
    fi
fi
export LD_LIBRARY_PATH

BUILD_SUCCEEDED=false
RUN_FILE_NAME="run.sh"

#for COMPILER in ${compilerList}
compilerListLen=${#compilerList[@]}
compilerListLenMinusOne="$(($compilerListLen-1))"
for i in $(seq 0 $compilerListLenMinusOne)
do

    csvCompilerList="${compilerList[$i]}"

    for COMPILER in $(echo ${csvCompilerList} | sed "s/,/ /g")
    do

        #### Infer the compiler flags

        unset COMPILER_FLAGS

        if [ -z ${USER_SELECTED_COMPILER+x} ] && [ -z ${USER_SELECTED_COMPILER_FLAGS+x} ]; then

            if [ "${pmExamLang}" = "C" ] || [ "${pmExamLang}" = "C++" ]; then

                if [ "${COMPILER}" = "icl" ] || [ "${COMPILER}" = "icc" ] || [ "${COMPILER}" = "icpc" ] || [ "${COMPILER}" = "mpiicc" ] || [ "${COMPILER}" = "mpiicpc" ]; then
                    COMPILER_FLAGS="${INTEL_C_COMPILER_FLAGS}"
                else
                    COMPILER_FLAGS="${GNU_C_COMPILER_FLAGS}"
                fi

            elif [ "${pmExamLang}" = "Fortran" ]; then

                # -DIS_COMPATIBLE_COMPILER

                if [ "${COMPILER}" = "ifort" ] || [ "${COMPILER}" = "mpiifort" ]; then
                    COMPILER_FLAGS="${INTEL_Fortran_COMPILER_FLAGS}"
                else
                    COMPILER_FLAGS="${GNU_Fortran_COMPILER_FLAGS}"
                fi

            fi

        else

            COMPILER_FLAGS="${USER_SELECTED_COMPILER_FLAGS}"

        fi

        echo >&2
        echo >&2 "-- ParaMonteExample${pmExamLang} - The ParaMonte example's compilation command:"
        echo >&2 "-- ParaMonteExample${pmExamLang} -     ${COMPILER} ${COMPILER_FLAGS} ${pmSrcFiles} -c"

        # Intel compiler does not support -dumpversion

        #CPATH="$(command -v "${COMPILER}")"
        #if [[ -f "${CPATH}" && "${CPATH}" =~ .*"ifort".* && "${CPATH}" =~ .*"intel".* ]]; then
        #    CVERSION="${COMPILER} --version"
        #elif [[ "${COMPILER}" =~ .*"gfortran".* ]]; then
        #    CVERSION="${COMPILER} -dumpversion"
        #fi
        #${CVERSION} >/dev/null 2>&1 && \

        ${COMPILER} ${COMPILER_FLAGS} ${pmSrcFiles} -c >/dev/null 2>&1 && \
        {

            echo >&2 "-- ParaMonteExample${pmExamLang} - The example's source file compilation appears to have succeeded."

            csvLinkerList=${COMPILER}
            LINKER_FLAGS=
            if ([ "${pmExamLang}" = "C" ] || [ "${pmExamLang}" = "C++" ]) && [ "${pmLibType}" = "static" ]; then
                if [ "${pmCompilerSuite}" = "intel" ]; then
                    if [ "${pmLibType}" = "static" ] && ( [ "${MPI_ENABLED}" = "true" ] || [ "${CAF_ENABLED}" = "true" ] ); then
                        csvLinkerList="mpiifort"
                    else
                        csvLinkerList="ifort"
                    fi
                    LINKER_FLAGS="-nofor_main"
                    #if [ "${MPI_ENABLED}" = "true" ]; then
                    #    LINKER="mpiifort" # xxx point of weakness: assumes intel mpi to have been installed
                    #fi
                fi
                if [ "${pmCompilerSuite}" = "gnu" ]; then
                    if [ "${pmLibType}" = "static" ] && ( [ "${MPI_ENABLED}" = "true" ] || [ "${CAF_ENABLED}" = "true" ] ); then
                        csvLinkerList="mpifort"
                    else
                        csvLinkerList="gfortran"
                    fi
                    #csvLinkerList="${clist}"
                    #if [ "${MPI_ENABLED}" = "true" ]; then
                    #    LINKER="mpifort"
                    #fi
                fi
            fi

            for LINKER in $(echo ${csvLinkerList} | sed "s/,/ /g")
            do

                echo >&2
                echo >&2 "-- ParaMonteExample${pmExamLang} - The ParaMonte example's compilation command:"
                echo >&2 "-- ParaMonteExample${pmExamLang} -     ${LINKER} ${COMPILER_FLAGS} ${LINKER_FLAGS} ${pmSrcFiles//.$pmSrcExt/.o} ${pmLibFullName} -o ${PM_EXAM_EXE_NAME}"

                # Intel compiler does not support -dumpversion
                #${LINKER} -dumpversion >/dev/null 2>&1 && \

                if ([ "$isMinGW" = "true" ] || [ "$isCygwin" = "true" ]) \
                && ([ "${LINKER}" = "icl" ] || [ "${LINKER}" = "ifort" ] || [ "${LINKER}" = "mpicl" ] || [ "${LINKER}" = "mpiifort" ]); then
                    pmLibFullNameCurrent="${pmLibBaseName}.lib"
                    pmObjExt="obj"
                else
                    pmLibFullNameCurrent="${pmLibFullName}"
                    pmObjExt="o"
                fi

                ${LINKER} ${COMPILER_FLAGS} ${LINKER_FLAGS} ${pmSrcFiles//.$pmSrcExt/.$pmObjExt} "${pmLibFullNameCurrent}" -o ${PM_EXAM_EXE_NAME} >/dev/null 2>&1 && \
                {
                #if [ $? -eq 0 ]; then

                    BUILD_SUCCEEDED=true

                    echo >&2 "-- ParaMonteExample${pmExamLang} - The example's source file linking appears to have succeeded."

                    {
                    echo "#!/bin/bash"
                    echo "# ParaMonte example runtime setup script."
                    echo "# "
                    } > ${RUN_FILE_NAME}
                    if [ "${MPI_ENABLED}" = "true" ] || [ "${CAF_ENABLED}" = "true" ]; then
                        {
                        echo "# usage:"
                        echo "# "
                        echo "#     ./${RUN_FILE_NAME} -n number_of_processors"
                        echo "# "
                        echo "# or,"
                        echo "# "
                        echo "#     source ./${RUN_FILE_NAME} -n number_of_processors"
                        echo "# "
                        echo "# where number_of_processors is an integer representing the number"
                        echo "# of physical processors on which the example will be run."
                        echo ""
                        echo "while [ \"\$1\" != \"\" ]; do"
                        echo "    case \$1 in"
                        echo "        -n | --nproc )        shift"
                        echo "                              FOR_COARRAY_NUM_IMAGES=\$1"
                        echo "                              ;;"
                        echo "        * )                   echo >\&2 \"-- ParaMonteExampleRunScript - FATAL: the input flag is not recognized: \$1\""
                        echo "                              exit 1"
                        echo "    esac"
                        echo "    shift"
                        echo "done"
                        echo ""
                        } >> ${RUN_FILE_NAME}
                    else
                        {
                        echo "# usage:"
                        echo "# "
                        echo "#     ./${RUN_FILE_NAME}"
                        echo "# "
                        echo "# or,"
                        echo "# "
                        echo "#     source ./${RUN_FILE_NAME}"
                        echo ""
                        } >> ${RUN_FILE_NAME}
                    fi
                    {
                    echo ""
                    echo "FILE_DIR=\"\$( cd \"\$( dirname \"\${BASH_SOURCE[0]}\" )\" >/dev/null 2>&1 && pwd )\""
                    echo "if [ -z \${PATH+x} ]; then"
                    echo "    PATH=."
                    echo "else"
                    echo "    if [[ \":\$PATH:\" != *\":${FILE_DIR}:\"* ]]; then"
                    echo "        PATH=\"${FILE_DIR}:\${PATH}\""
                    echo "    fi"
                    echo "fi"
                    echo "export LD_LIBRARY_PATH"
                    echo "if [ -z \${LD_LIBRARY_PATH+x} ]; then"
                    echo "    LD_LIBRARY_PATH=${FILE_DIR}"
                    echo "else"
                    echo "    if [[ \":\$LD_LIBRARY_PATH:\" != *\":${FILE_DIR}:\"* ]]; then"
                    echo "        LD_LIBRARY_PATH=\"${FILE_DIR}:\${LD_LIBRARY_PATH}\""
                    echo "    fi"
                    echo "fi"
                    echo "export LD_LIBRARY_PATH"
                    echo "export PATH"
                    echo ""
                    echo ""
                    echo "if [ -z \${FOR_COARRAY_NUM_IMAGES+x} ]; then"
                    echo "    FOR_COARRAY_NUM_IMAGES=${FOR_COARRAY_NUM_IMAGES}"
                    echo "fi"
                    echo ""
                    } >> ${RUN_FILE_NAME}

                    if [ -f "${FILE_DIR}/setup.sh" ]; then
                        echo "source ${FILE_DIR}/setup.sh" >> ${RUN_FILE_NAME}
                        echo "" >> ${RUN_FILE_NAME}
                    fi

                    echo "# run ParaMonte example executable" >> ${RUN_FILE_NAME}
                    echo "" >> ${RUN_FILE_NAME}
                    echo "chmod +x ${PM_EXAM_EXE_NAME}" >> ${RUN_FILE_NAME}
                    if [ "${MPI_ENABLED}" = "true" ]; then
                        echo "mpiexec -n \${FOR_COARRAY_NUM_IMAGES} ./${PM_EXAM_EXE_NAME} || mpiexec --oversubscribe -n \${FOR_COARRAY_NUM_IMAGES} ./${PM_EXAM_EXE_NAME}" >> ${RUN_FILE_NAME}
                        echo "" >> ${RUN_FILE_NAME}
                    else
                        if [ "${CAF_ENABLED}" = "true" ]; then
                            if [ "${pmCompilerSuite}" = "intel" ]; then
                                echo "export FOR_COARRAY_NUM_IMAGES && ./${PM_EXAM_EXE_NAME}" >> ${RUN_FILE_NAME}
                            else
                                echo "cafrun -np \${FOR_COARRAY_NUM_IMAGES} ./${PM_EXAM_EXE_NAME}" >> ${RUN_FILE_NAME}
                            fi
                        else
                            echo "./${PM_EXAM_EXE_NAME}" >> ${RUN_FILE_NAME}
                        fi
                    fi

                    chmod +x ${RUN_FILE_NAME}

                    break

                } || {
                #else

                    echo >&2 "-- ParaMonteExample${pmExamLang} - The example's source files compilation and linking appear to have failed. skipping..."
                    echo >&2 "-- ParaMonteExample${pmExamLang} - If the compiler is missing or unidentified, you can pass the path to the compiler to the build script:"
                    echo >&2 "-- ParaMonteExample${pmExamLang} - For instructions, type on the command line:"
                    echo >&2 "-- ParaMonteExample${pmExamLang} - "
                    echo >&2 "-- ParaMonteExample${pmExamLang} -     cd $(pwd)"
                    echo >&2 "-- ParaMonteExample${pmExamLang} -     ./build.sh --help"
                    echo >&2

                #fi
                }

            done

        } || {

            echo >&2 "-- ParaMonteExample${pmExamLang} - The example's source files compilation and linking appear to have failed. skipping..."
            echo >&2 "-- ParaMonteExample${pmExamLang} - If the compiler is missing or unidentified, you can pass the path to the compiler to the build script."
            echo >&2 "-- ParaMonteExample${pmExamLang} - For instructions, type on the command line:"
            echo >&2 "-- ParaMonteExample${pmExamLang} - "
            echo >&2 "-- ParaMonteExample${pmExamLang} -     cd $(pwd)"
            echo >&2 "-- ParaMonteExample${pmExamLang} -     ./build.sh --help"
            echo >&2

        }

        if [ "${BUILD_SUCCEEDED}" = "true" ]; then break; fi

    done

    if [ "${BUILD_SUCCEEDED}" = "true" ]; then break; fi

done

echo >&2
if [ "${BUILD_SUCCEEDED}" = "true" ]; then
    echo >&2 "-- ParaMonteExample${pmExamLang} - To run the example's executable with the proper environmental setup, try:"
    echo >&2 "-- ParaMonteExample${pmExamLang} - "
    echo >&2 "-- ParaMonteExample${pmExamLang} -     ./${RUN_FILE_NAME}"
    echo >&2 "-- ParaMonteExample${pmExamLang} - "
    echo >&2 "-- ParaMonteExample${pmExamLang} - or,"
    echo >&2 "-- ParaMonteExample${pmExamLang} - "
    echo >&2 "-- ParaMonteExample${pmExamLang} -     source ./${RUN_FILE_NAME}"
    echo >&2 "-- ParaMonteExample${pmExamLang} - "
    echo >&2 "-- ParaMonteExample${pmExamLang} - Look at the contents of ${RUN_FILE_NAME} to change the runtime settings as you wish."
    echo >&2
    exit 0
else
    echo >&2 "-- ParaMonteExample${pmExamLang} - exhausted all possible compiler names to build and run ParaMonte example but failed."
    echo >&2 "-- ParaMonteExample${pmExamLang} - Please consider reporting the circumstances surrounding this issue at"
    echo >&2 "-- ParaMonteExample${pmExamLang} - "
    echo >&2 "-- ParaMonteExample${pmExamLang} -    https://github.com/cdslaborg/paramonte/issues"
    echo >&2 "-- ParaMonteExample${pmExamLang} - "
    echo >&2 "-- ParaMonteExample${pmExamLang} - or directly to ParaMonte authors (e.g., shahmoradi@utexas.edu)"
    echo >&2 "-- ParaMonteExample${pmExamLang} - "
    echo >&2 "-- ParaMonteExample${pmExamLang} - gracefully exiting ParaMonte"
    echo >&2
    exit 1
fi
