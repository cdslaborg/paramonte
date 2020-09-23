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
NOTE:     source ./run.sh

EndOfMessage
}

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
# build ParaMonteExample objects and executable
####################################################################################################################################

PM_EXAM_EXE_NAME="main.exe"
export PM_EXAM_EXE_NAME

if [ -z ${FOR_COARRAY_NUM_IMAGES+x} ]; then
    FOR_COARRAY_NUM_IMAGES=3
fi
export FOR_COARRAY_NUM_IMAGES

PMLIB_FULL_PATH="$(ls -d ${FILE_DIR}/*libparamonte_* | sort -V | tail -n1)"
PMLIB_FULL_NAME=${PMLIB_FULL_PATH##*/}
PMLIB_BASE_NAME=${PMLIB_FULL_NAME%.*}

####################################################################################################################################
# get ParaMonte's compiler suite
####################################################################################################################################

unset FILE_EXT
unset SRC_FILES
unset COMPILER_LIST
unset EXAMPLE_LANGUAGE
unset PM_COMPILER_SUITE

if [[ "$PMLIB_FULL_NAME" =~ .*"_fortran_".* ]]; then
    FILE_EXT=f90
    SRC_FILES="paramonte.${FILE_EXT}"; export SRC_FILES
    EXAMPLE_LANGUAGE=Fortran
    if [[ "$PMLIB_FULL_NAME" =~ .*"_intel_".* ]]; then
        PM_COMPILER_SUITE=intel
        COMPILER_LIST=ifort
    fi
    if [[ "$PMLIB_FULL_NAME" =~ .*"_gnu_".* ]]; then
        PM_COMPILER_SUITE=gnu
        COMPILER_LIST=gfortran
    fi
fi

if [[ "$PMLIB_FULL_NAME" =~ .*"_c_".* ]]; then
    FILE_EXT=c
    EXAMPLE_LANGUAGE=C
    if [[ "$PMLIB_FULL_NAME" =~ .*"_intel_".* ]]; then
        PM_COMPILER_SUITE=intel
        COMPILER_LIST=icc
    fi
    if [[ "$PMLIB_FULL_NAME" =~ .*"_gnu_".* ]]; then
        PM_COMPILER_SUITE=gnu
        COMPILER_LIST=gcc
    fi
fi

if [[ "$PMLIB_FULL_NAME" =~ .*"_cpp_".* ]]; then
    FILE_EXT=cpp
    EXAMPLE_LANGUAGE=C++
    if [[ "$PMLIB_FULL_NAME" =~ .*"_intel_".* ]]; then
        PM_COMPILER_SUITE=intel
        COMPILER_LIST=icpc
    fi
    if [[ "$PMLIB_FULL_NAME" =~ .*"_gnu_".* ]]; then
        PM_COMPILER_SUITE=gnu
        COMPILER_LIST=g++
    fi
fi

SRC_FILES="${SRC_FILES} logfunc.${FILE_EXT} main.${FILE_EXT}"

echo >&2
echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - ParaMonte library's full path: ${PMLIB_FULL_PATH}"
echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - ParaMonte library's full name: ${PMLIB_FULL_NAME}"
echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - ParaMonte library's base name: ${PMLIB_BASE_NAME}"

####################################################################################################################################
# build setup
####################################################################################################################################

unset PM_BLD_TYPE
unset GNU_C_COMPILER_FLAGS
if [[ "$PMLIB_FULL_NAME" =~ .*"darwin".* ]]; then
    INTEL_C_COMPILER_FLAGS="-no-multibyte-chars"
else
    unset INTEL_C_COMPILER_FLAGS
fi
GNU_Fortran_COMPILER_FLAGS="-std=gnu"
INTEL_Fortran_COMPILER_FLAGS="-standard-semantics"

if [[ "$PMLIB_FULL_NAME" =~ .*"release".* ]]; then
    PM_BLD_TYPE=release
    GNU_C_COMPILER_FLAGS="${GNU_C_COMPILER_FLAGS} -O3"
    GNU_Fortran_COMPILER_FLAGS="${GNU_Fortran_COMPILER_FLAGS} -O3 -funroll-loops -finline-functions -ftree-vectorize"
    INTEL_C_COMPILER_FLAGS="${INTEL_C_COMPILER_FLAGS} -O3"
    INTEL_Fortran_COMPILER_FLAGS="${INTEL_Fortran_COMPILER_FLAGS} -O3 -ip -ipo -unroll -unroll-aggressive -finline-functions"
fi
if [[ "$PMLIB_FULL_NAME" =~ .*"testing".* ]]; then
    PM_BLD_TYPE=testing
    GNU_C_COMPILER_FLAGS="${GNU_C_COMPILER_FLAGS} -O0"
    GNU_Fortran_COMPILER_FLAGS="${GNU_Fortran_COMPILER_FLAGS} -O0"
    INTEL_C_COMPILER_FLAGS="${INTEL_C_COMPILER_FLAGS} -O0"
    INTEL_Fortran_COMPILER_FLAGS="${INTEL_Fortran_COMPILER_FLAGS} -O0"
fi
if [[ "$PMLIB_FULL_NAME" =~ .*"debug".* ]]; then
    PM_BLD_TYPE=debug
    GNU_C_COMPILER_FLAGS="${GNU_C_COMPILER_FLAGS} -O0 -g"
    GNU_Fortran_COMPILER_FLAGS="${GNU_Fortran_COMPILER_FLAGS} -O0 -g"
    INTEL_C_COMPILER_FLAGS="${INTEL_C_COMPILER_FLAGS} -O0 -debug full"
    INTEL_Fortran_COMPILER_FLAGS="${INTEL_Fortran_COMPILER_FLAGS} -O0 -debug full"
fi
echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - ParaMonte build type: ${PM_BLD_TYPE}"

####################################################################################################################################
# library type
####################################################################################################################################

if [[ "${PMLIB_FULL_NAME}" =~ .*"dynamic".* ]]; then
    PM_LIB_TYPE=dynamic
else
    PM_LIB_TYPE=static
fi
echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - ParaMonte library type: ${PM_LIB_TYPE}"

####################################################################################################################################
# MPI
####################################################################################################################################

MPI_ENABLED=false
if [[ "$PMLIB_FULL_NAME" =~ .*"mpi".* ]]; then
    MPI_ENABLED=true
fi
echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - MPI_ENABLED: ${MPI_ENABLED}"

if [ "${MPI_ENABLED}" = "true" ]; then
    if [ "${EXAMPLE_LANGUAGE}" = "Fortran" ]; then
        COMPILER_LIST="mpiifort mpifort"
    fi
    if [ "${EXAMPLE_LANGUAGE}" = "C" ]; then
        COMPILER_LIST="mpiicc mpicc"
    fi
    if [ "${EXAMPLE_LANGUAGE}" = "C++" ]; then
        COMPILER_LIST="mpiicpc mpicxx mpic++ mpicc"
    fi
fi

####################################################################################################################################
# CAF
####################################################################################################################################

CAFTYPE=none
CAF_ENABLED=false
if [[ "$PMLIB_FULL_NAME" =~ .*"cafsingle".* ]]; then
    #-fcoarray=single
    GNU_Fortran_COMPILER_FLAGS="${GNU_Fortran_COMPILER_FLAGS}"
    INTEL_Fortran_COMPILER_FLAGS="${INTEL_Fortran_COMPILER_FLAGS} -coarray=single"
    CAF_ENABLED=true
fi
if [[ "$PMLIB_FULL_NAME" =~ .*"cafshared".* ]]; then
    #-fcoarray=shared
    GNU_Fortran_COMPILER_FLAGS="${GNU_Fortran_COMPILER_FLAGS}"
    INTEL_Fortran_COMPILER_FLAGS="${INTEL_Fortran_COMPILER_FLAGS} -coarray=shared"
    CAF_ENABLED=true
fi
if [[ "$PMLIB_FULL_NAME" =~ .*"cafdistributed".* ]]; then
    #-fcoarray=distributed
    GNU_Fortran_COMPILER_FLAGS="${GNU_Fortran_COMPILER_FLAGS}"
    INTEL_Fortran_COMPILER_FLAGS="${INTEL_Fortran_COMPILER_FLAGS} -coarray=distributed"
    CAF_ENABLED=true
fi
echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - CAFTYPE: ${CAFTYPE}"

if [ "${PM_COMPILER_SUITE}" = "gnu" ] && [ "${CAF_ENABLED}" = "true" ]; then
    COMPILER_LIST=caf
fi

####################################################################################################################################
# report compiler choice
####################################################################################################################################

echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - ParaMonte library's compiler suite: ${PM_COMPILER_SUITE}"
echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - inferred compiler choice(s): ${COMPILER_LIST}"

if [ -z ${USER_SELECTED_COMPILER+x} ]; then
    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - user-selected compiler/linker: none"
else
    COMPILER_LIST="${USER_SELECTED_COMPILER}"
    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - user-selected compiler/linker: ${USER_SELECTED_COMPILER}"
fi

if [ -z ${USER_SELECTED_COMPILER+x} ] && [ -z ${USER_SELECTED_COMPILER_FLAGS+x} ]; then
    if [ "${PM_COMPILER_SUITE}" = "intel" ]; then
        if [ "${EXAMPLE_LANGUAGE}" = "C" ] || [ "${EXAMPLE_LANGUAGE}" = "C++" ]; then COMPILER_FLAGS=${INTEL_C_COMPILER_FLAGS}; fi
        if [ "${EXAMPLE_LANGUAGE}" = "Fortran" ]; then COMPILER_FLAGS="${INTEL_Fortran_COMPILER_FLAGS} -fpp"; fi # -DIS_COMPATIBLE_COMPILER
    fi
    if [ "${PM_COMPILER_SUITE}" = "gnu" ]; then
        if [ "${EXAMPLE_LANGUAGE}" = "C" ] || [ "${EXAMPLE_LANGUAGE}" = "C++" ]; then COMPILER_FLAGS=${GNU_C_COMPILER_FLAGS}; fi
        if [ "${EXAMPLE_LANGUAGE}" = "Fortran" ]; then COMPILER_FLAGS="${GNU_Fortran_COMPILER_FLAGS} -cpp"; fi # -DIS_COMPATIBLE_COMPILER
    fi
    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - inferred compiler/linker flags(s): ${COMPILER_FLAGS}"
else
    COMPILER_FLAGS="${USER_SELECTED_COMPILER_FLAGS}"
    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - user-selected compiler/linker flags: ${USER_SELECTED_COMPILER_FLAGS}"
fi

####################################################################################################################################
# build example
####################################################################################################################################

BUILD_SUCCEEDED=false
RUN_FILE_NAME="run.sh"

for COMPILER in ${COMPILER_LIST}
do

    echo >&2
    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - compiling ParaMonte example with ${COMPILER}"
    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - ${COMPILER} ${COMPILER_FLAGS} ${SRC_FILES} ${PMLIB_FULL_PATH} -o ${PM_EXAM_EXE_NAME}"

    ${COMPILER} ${COMPILER_FLAGS} ${SRC_FILES} -c

    LINKER=${COMPILER}
    LINKER_FLAGS=
    if ([ "${EXAMPLE_LANGUAGE}" = "C" ] || [ "${EXAMPLE_LANGUAGE}" = "C++" ]) && [ "${PM_LIB_TYPE}" = "static" ]; then
        if [ "${PM_COMPILER_SUITE}" = "intel" ]; then
            LINKER_FLAGS="-nofor_main"
            if [ "${MPI_ENABLED}" = "true" ]; then
                LINKER="mpiifort" # xxx point of weakness: assumes intel mpi to have been installed
            else
                LINKER="ifort"
            fi
        fi
        if [ "${PM_COMPILER_SUITE}" = "gnu" ]; then
            if [ "${MPI_ENABLED}" = "true" ]; then
                LINKER="mpifort"
            else
                LINKER="gfortran"
            fi
        fi
    fi

    echo >&2
    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - linking ParaMonte example with ${LINKER}"
    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - ${LINKER} ${COMPILER_FLAGS} ${LINKER_FLAGS} logfunc.o main.o ${PMLIB_FULL_PATH} -o ${PM_EXAM_EXE_NAME}"

    ${LINKER} ${COMPILER_FLAGS} ${LINKER_FLAGS} logfunc.o main.o ${PMLIB_FULL_NAME} -o ${PM_EXAM_EXE_NAME}

    if [ $? -eq 0 ]; then

        BUILD_SUCCEEDED=true

        echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - example build appears to have succeeded."

        {
        echo "# ParaMonte example runtime setup script."
        echo "# "
        } > ${RUN_FILE_NAME}
        if [ "${MPI_ENABLED}" = "true" ] || [ "${CAF_ENABLED}" = "true" ]; then
            {
            echo "# usage:"
            echo "# "
            echo "#     source ./${RUN_FILE_NAME} -n number_of_processors"
            echo "# "
            echo "# where number_of_processors is an integer representing the number"
            echo "# of physical proceesors on which the example will be run."
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
        echo "    LD_LIBRARY_PATH=."
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

        if [ -f "./setup.sh" ]; then
            echo "source ./setup.sh" >> ${RUN_FILE_NAME}
            echo "" >> ${RUN_FILE_NAME}
        fi

        echo "# run ParaMonte example executable" >> ${RUN_FILE_NAME}
        echo "" >> ${RUN_FILE_NAME}
        echo "chmod +x ${PM_EXAM_EXE_NAME}" >> ${RUN_FILE_NAME}
        if [ "${MPI_ENABLED}" = "true" ]; then
            echo "mpiexec -n \${FOR_COARRAY_NUM_IMAGES} ./${PM_EXAM_EXE_NAME}" >> ${RUN_FILE_NAME}
            echo "" >> ${RUN_FILE_NAME}
        else
            if [ "${CAF_ENABLED}" = "true" ]; then
                if [ "${PM_COMPILER_SUITE}" = "intel" ]; then
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

    else

        echo >&2
        echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - example build appears to have failed. skipping..."

    fi

done

echo >&2

if [ "${BUILD_SUCCEEDED}" = "true" ]; then

    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - To run the example's executable with the proper environmental setup, try:"
    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - "
    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} -     source ./${RUN_FILE_NAME}"
    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - "
    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - Look at the contents of ${RUN_FILE_NAME} to change the runtime settings as you wish."
    echo >&2

else

    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - exhausted all possible compiler names to build and run ParaMonte example but failed."
    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - Please consider reporting the circumstances surrounding this issue at"
    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - "
    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} -    https://github.com/cdslaborg/paramonte/issues"
    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - "
    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - or directly to ParaMonte authors (e.g., shahmoradi@utexas.edu)"
    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - "
    echo >&2 "-- ParaMonteExample${EXAMPLE_LANGUAGE} - gracefully exiting ParaMonte"
    echo >&2

    exit 1

fi

exit 0
