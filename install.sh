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
# this script will parse the user-provided flags and their values 
# to build the ParaMonte library.
# to redirect output to the external file install.sh.out, try:
# install.sh >install.sh.out 2>&1

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

echo "$(cat ./auxil/ParaMonteBanner.txt)"

usage()
{
    echo "$(cat ${ParaMonte_ROOT_DIR}/install.sh.usage.txt)"
}

unset LANG_LIST
unset BTYPE_LIST
unset LTYPE_LIST
unset CAFTYPE_LIST
unset MPI_ENABLED_LIST
unset HEAP_ARRAY_ENABLED_LIST
unset ParaMonteTest_RUN_ENABLED
unset ParaMonteExample_RUN_ENABLED
unset Fortran_COMPILER_PATH
unset MPIEXEC_PATH

PMCS_LIST="none"
fresh_flag=""
yes_to_all_flag=""
gcc_bootstrap_flag=""
FOR_COARRAY_NUM_IMAGES=3

while [ "$1" != "" ]; do
    case $1 in
        -L | --lang )           shift
                                LANG_LIST="$1"
                                ;;
        -s | --compiler_suite ) shift
                                PMCS_LIST="$1"
                                ;;
        -b | --build_type )     shift
                                BTYPE_LIST="$1"
                                ;;
        -l | --lib_type )       shift
                                LTYPE_LIST="$1"
                                ;;
        -c | --caf_type )       shift
                                CAFTYPE_LIST="$1"
                                ;;
        -m | --mpi_enabled )    shift
                                MPI_ENABLED_LIST="$1"
                                ;;
        -e | --heap_enabled )   shift
                                HEAP_ARRAY_ENABLED_LIST="$1"
                                ;;
        -t | --test_enabled )   shift
                                ParaMonteTest_RUN_ENABLED="$1"
                                ;;
        -x | --exam_enabled )   shift
                                ParaMonteExample_RUN_ENABLED="$1"
                                ;;
        -f | --fortran )        shift
                                Fortran_COMPILER_PATH="$1"
                                ;;
        -M | --mpiexec )        shift
                                MPIEXEC_PATH="$1"
                                ;;
        -F | --fresh )          fresh_flag="--fresh"
                                ;;
        -y | --yes-to-all )     yes_to_all_flag="--yes-to-all"
                                ;;
        -B | --bootstrap )      shift
                                gcc_bootstrap_flag="--bootstrap"
                                ;;
        -n | --num_images )     shift
                                FOR_COARRAY_NUM_IMAGES="$1"
                                ;;
        -h | --help )           usage
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
# auxil
####################################################################################################################################

getLowerCaseChar()
{
    if [[ $1 =~ [A-Z] ]];then
        n=$(printf "%d" "'$1")
        n=$((n+32))
        printf \\$(printf "%o" "$n")
    else
        printf "%s" "$1"
    fi
}

getUpperCaseChar()
{
    if [[ $1 =~ [a-z] ]];then
        n=$(printf "%d" "'$1")
        n=$((n-32))
        printf \\$(printf "%o" "$n")
    else
        printf "%s" "$1"
    fi
}

getLowerCase() {
    word="$@"
    for((i=0;i<${#word};i++)); do
        ch="${word:$i:1}"
        getLowerCaseChar "$ch"
    done
}

getUpperCase() {
    word="$@"
    for((i=0;i<${#word};i++)); do
        ch="${word:$i:1}"
        getUpperCaseChar "$ch"
    done
}

isnumeric() {
    isNumericValue=true
    word="$@"
    for((i=0;i<${#word};i++)); do
        ch="${word:$i:1}"
        if ! [[ $ch =~ [0-9] ]];then
            isNumericValue=false
            break
        fi
    done
    echo $isNumericValue
}

####################################################################################################################################
# verify arguments
####################################################################################################################################

reportConflict()
{
    usage
    echo >&2 ""
    echo >&2 "-- ParaMonte - FATAL: conflicting flag values detected."
    echo >&2 "-- ParaMonte - FATAL: $1"
    echo >&2 ""
    echo >&2 "-- ParaMonte - gracefully exiting."
    echo >&2 ""
    echo >&2 ""
    exit 1
}

reportBadValue()
{
    usage
    echo >&2 ""
    echo >&2 "-- ParaMonte - FATAL: The requested input value $2 specified with "
    echo >&2 "-- ParaMonte - FATAL: the input flag $1 is not supported."
    if ! [ -z ${3+x} ]; then
    echo >&2 "-- ParaMonte - FATAL: $3"
    fi
    echo >&2 ""
    echo >&2 "-- ParaMonte - gracefully exiting."
    echo >&2 ""
    echo >&2 ""
    exit 1
}

if ! [ -z ${LANG_LIST+x} ]; then
    for LANG in $LANG_LIST; do
        if  [[ $LANG != [cC]
            && $LANG != [fF][oO][rR][tT][rR][aA][nN]
            && $LANG != [pP][yY][tT][hH][oO][nN] ]]; then
            reportBadValue "-L or --lang" $LANG
        fi
    done
    if [ "${LANG_LIST}" = "" ]; then
        unset LANG_LIST
    else
        LANG_LIST="$(getLowerCase $LANG_LIST)"
    fi
fi

if ! [ -z ${PMCS_LIST+x} ]; then
    for PMCS in $PMCS_LIST; do
        if  [[ $PMCS != [nN][oO][nN][eE] 
            && $PMCS != [iI][nN][tT][eE][lL]
            && $PMCS != [gG][nN][uU] ]]; then
            reportBadValue "-s or --compiler_suite" $PMCS
        fi
    done
    if [ "${PMCS_LIST}" = "" ]; then
        unset PMCS_LIST
    else
        PMCS_LIST="$(getLowerCase $PMCS_LIST)"
    fi
fi

if ! [ -z ${BTYPE_LIST+x} ]; then
    for BTYPE in $BTYPE_LIST; do
        if  [[ $BTYPE != [rR][eE][lL][eE][aA][sS][eE]
            && $BTYPE != [tT][eE][sS][tT][iI][nN][gG] 
            && $BTYPE != [dD][eE][bB][uU][gG] ]]; then
            reportBadValue "-b or --build_type" $BTYPE
        fi
    done
    if [ "${BTYPE_LIST}" = "" ]; then
        unset BTYPE_LIST
    else
        BTYPE_LIST="$(getLowerCase $BTYPE_LIST)"
    fi
fi

if ! [ -z ${LTYPE_LIST+x} ]; then
    for LTYPE in $LTYPE_LIST; do
        if  [[ $LTYPE != [dD][yY][nN][aA][mM][iI][cC] 
            && $LTYPE != [sS][tT][aA][tT][iI][cC] ]]; then
            reportBadValue "-l or --lib_type" $LTYPE
        fi
    done
    if [ "${LTYPE_LIST}" = "" ]; then
        unset LTYPE_LIST
    else
        LTYPE_LIST="$(getLowerCase $LTYPE_LIST)"
    fi
fi

if ! [ -z ${CAFTYPE_LIST+x} ]; then
    for CAFTYPE in $CAFTYPE_LIST; do
        if  [[ $CAFTYPE != [nN][oO][nN][eE] 
            && $CAFTYPE != [sS][iI][nN][gG][lL][eE] 
            && $CAFTYPE != [sS][hH][aA][rR][eE][dD] 
            && $CAFTYPE != [dD][iI][sS][tT][rR][iI][bB][uU][tT][eE][dD] ]]; then
            reportBadValue "-c or --caf_type" $CAFTYPE
        fi
    done
    if [ "${CAFTYPE_LIST}" = "" ]; then
        unset CAFTYPE_LIST
    else
        CAFTYPE_LIST="$(getLowerCase $CAFTYPE_LIST)"
    fi
fi

if ! [ -z ${MPI_ENABLED_LIST+x} ]; then
    for MPI_ENABLED in $MPI_ENABLED_LIST; do
        if  [[ $MPI_ENABLED != [tT][rR][uU][eE] 
            && $MPI_ENABLED != [fF][aA][lL][sS][eE] ]]; then
            reportBadValue "-m or --mpi_enabled" $MPI_ENABLED
        fi
    done
    if [ "${MPI_ENABLED_LIST}" = "" ]; then
        unset MPI_ENABLED_LIST
    else
        MPI_ENABLED_LIST="$(getLowerCase $MPI_ENABLED_LIST)"
    fi
fi

if ! [ -z ${HEAP_ARRAY_ENABLED_LIST+x} ]; then
    for HEAP_ARRAY_ENABLED in $HEAP_ARRAY_ENABLED_LIST; do
        if  [[ $HEAP_ARRAY_ENABLED != [tT][rR][uU][eE] 
            && $HEAP_ARRAY_ENABLED != [fF][aA][lL][sS][eE] ]]; then
            reportBadValue "-e or --heap_enabled" $HEAP_ARRAY_ENABLED
        fi
    done
    if [ "${HEAP_ARRAY_ENABLED_LIST}" = "" ]; then
        unset HEAP_ARRAY_ENABLED_LIST
    else
        HEAP_ARRAY_ENABLED_LIST="$(getLowerCase $HEAP_ARRAY_ENABLED_LIST)"
    fi
fi

if ! [ -z ${ParaMonteTest_RUN_ENABLED+x} ]; then
    if  [[ $ParaMonteTest_RUN_ENABLED != [tT][rR][uU][eE] 
        && $ParaMonteTest_RUN_ENABLED != [fF][aA][lL][sS][eE] ]]; then
        reportBadValue "-t or --test_enabled" $ParaMonteTest_RUN_ENABLED
    fi
    if [ "${ParaMonteTest_RUN_ENABLED}" = "" ]; then
        unset ParaMonteTest_RUN_ENABLED
    else
        ParaMonteTest_RUN_ENABLED="$(getLowerCase $ParaMonteTest_RUN_ENABLED)"
    fi
fi

if ! [ -z ${ParaMonteExample_RUN_ENABLED+x} ]; then
    if  [[ $ParaMonteExample_RUN_ENABLED != [tT][rR][uU][eE] 
        && $ParaMonteExample_RUN_ENABLED != [fF][aA][lL][sS][eE] ]]; then
        reportBadValue "-x or --exam_enabled" $ParaMonteExample_RUN_ENABLED
    fi
    if [ "${ParaMonteExample_RUN_ENABLED}" = "" ]; then
        unset ParaMonteExample_RUN_ENABLED
    else
        ParaMonteExample_RUN_ENABLED="$(getLowerCase $ParaMonteExample_RUN_ENABLED)"
    fi
fi

fortran_flag=""
if ! [ -z ${Fortran_COMPILER_PATH+x} ]; then
    if [[ -f "${Fortran_COMPILER_PATH}" ]]; then
        fortran_flag="--mpiexec ${Fortran_COMPILER_PATH}"
    else
        if [ "${Fortran_COMPILER_PATH}" = "" ]; then
            unset Fortran_COMPILER_PATH
        else
            reportBadValue "-f or --fortran" "${Fortran_COMPILER_PATH}" "The value specified must be the path to the Fortran compiler executable file."
        fi
    fi
fi

mpiexec_flag=""
if ! [ -z ${MPIEXEC_PATH+x} ]; then
    if [[ -f "${MPIEXEC_PATH}" ]]; then
        mpiexec_flag="--mpiexec ${MPIEXEC_PATH}"
    else
        if [ "${MPIEXEC_PATH}" = "" ]; then
            unset MPIEXEC_PATH
        else
            reportBadValue "-M or --mpiexec" "${MPIEXEC_PATH}" "The value specified must be the path to the mpiexec executable file."
        fi
    fi
fi

num_images_flag=""
if ! [ -z ${FOR_COARRAY_NUM_IMAGES+x} ]; then
    isNumericValue="$(isnumeric ${FOR_COARRAY_NUM_IMAGES})"
    if [ "${isNumericValue}" = "true" ]; then
        num_images_flag="--num_images ${FOR_COARRAY_NUM_IMAGES}"
    else
        reportBadValue "-n or --num_images" $FOR_COARRAY_NUM_IMAGES "The number of processors must be a positive integer."
    fi
fi

####################################################################################################################################
# verify arguments consistencies
####################################################################################################################################

if ! [ -z ${CAFTYPE_LIST+x} ]; then
    for CAFTYPE in $CAFTYPE_LIST; do
        if  [ "${CAFTYPE}" != "none" ]; then
            for LANG in $LANG_LIST; do
                if  [ "${LANG}" != "fortran" ]; then
                    reportConflict "Coarray Fortran parallelism cannot be used to build library for ${LANG} language."
                fi
            done
            for MPI_ENABLED in $MPI_ENABLED_LIST; do
                if  [ "${MPI_ENABLED}" = "true" ]; then
                    reportConflict "Coarray Fortran parallelism cannot be mixed with MPI."
                fi
            done
            for LTYPE in $LTYPE_LIST; do
                if  [ "${LTYPE}" = "dynamic" ]; then
                    reportConflict "Coarray Fortran parallelism cannot be used with dynamic library build option."
                fi
            done
        fi
    done
fi

if ! [ -z ${LANG_LIST+x} ]; then
    for LANG in $LANG_LIST; do
        if  [ "${LANG}" = "python" ]; then
            for LTYPE in $LTYPE_LIST; do
                if  [ "${LTYPE}" = "static" ]; then
                    reportConflict "ParaMonte static library build is not possible for usage from Python language."
                fi
            done
        fi
    done
fi

####################################################################################################################################
# configure build
####################################################################################################################################

if [ -z ${LANG_LIST+x} ]; then
    LANG_LIST="fortran c python"
fi
CFI_ENABLED_LIST=""
C_IS_MISSING=true
Fortran_IS_MISSING=true
for LANG in $LANG_LIST; do
    if [ "${LANG}" = "fortran" ]; then
        if [ "${Fortran_IS_MISSING}" = "true" ]; then
            CFI_ENABLED_LIST="${CFI_ENABLED_LIST} false"
            Fortran_IS_MISSING=false
        fi
    fi
    if [ "${LANG}" = "python" ] || [ "${LANG}" = "c" ]; then
        if [ "${C_IS_MISSING}" = "true" ]; then
            CFI_ENABLED_LIST="${CFI_ENABLED_LIST} true"
            C_IS_MISSING=false
        fi
    fi
done

if [ -z ${PMCS_LIST+x} ]; then
    PMCS_LIST="none"
fi
if [ -z ${BTYPE_LIST+x} ]; then
    BTYPE_LIST="release testing debug"
fi
if [ -z ${LTYPE_LIST+x} ]; then
    LTYPE_LIST="static dynamic"
fi
if [ -z ${CAFTYPE_LIST+x} ]; then
    CAFTYPE_LIST="none single shared distributed"
fi
if [ -z ${MPI_ENABLED_LIST+x} ]; then
    MPI_ENABLED_LIST="true false"
fi
if [ -z ${HEAP_ARRAY_ENABLED_LIST+x} ]; then
    HEAP_ARRAY_ENABLED_LIST="true false"
fi
if [ -z ${ParaMonteTest_RUN_ENABLED+x} ]; then
    ParaMonteTest_RUN_ENABLED="true"
fi
if [ -z ${ParaMonteExample_RUN_ENABLED+x} ]; then
    ParaMonteExample_RUN_ENABLED="true"
fi

if [ "${LANG_LIST}" = "python" ]; then
    LTYPE_LIST="dynamic"
    CAFTYPE_LIST="none"
    HEAP_ARRAY_ENABLED_LIST="true"
fi

for PMCS in $PMCS_LIST; do

    for CFI_ENABLED in $CFI_ENABLED_LIST; do

        for BTYPE in $BTYPE_LIST; do

            for LTYPE in $LTYPE_LIST; do

                for HEAP_ARRAY_ENABLED in $HEAP_ARRAY_ENABLED_LIST; do

                    for MPI_ENABLED in $MPI_ENABLED_LIST; do

                        for CAFTYPE in $CAFTYPE_LIST; do

                            BENABLED=true

                            test_enabled_flag="--test_enabled ${ParaMonteTest_RUN_ENABLED}"
                            exam_enabled_flag="--exam_enabled ${ParaMonteExample_RUN_ENABLED}"

                            if [ "${PMCS}" = "none" ]; then
                               compiler_suite_flag=""
                            else
                               compiler_suite_flag="--compiler_suite ${PMCS}"
                            fi

                            caftype_flag="--caf_type ${CAFTYPE}"
                            lib_type_flag="--lib_type ${LTYPE}"
                            build_type_flag="--build_type ${BTYPE}"
                            cfi_enabled_flag="--cfi_enabled ${CFI_ENABLED}"
                            mpi_enabled_flag="--mpi_enabled ${MPI_ENABLED}"
                            heap_enabled_flag="--heap_enabled ${HEAP_ARRAY_ENABLED}"

                            # verify no conflict

                            #if [ "${LTYPE}" = "dynamic" ] && [ "${HEAP_ARRAY_ENABLED}" != "true" ]; then
                            #    BENABLED=false
                            #fi

                            if ! [ "${CAFTYPE}" = "none" ]; then
                                if [ "${MPI_ENABLED}" = "true" ] || [ "${CFI_ENABLED}" = "true" ] || [ "${LTYPE}" = "dynamic" ]; then
                                    BENABLED=false
                                fi
                            fi

                            if [ "${LTYPE}" = "dynamic" ] && [ "${HEAP_ARRAY_ENABLED}" = "false" ]; then
                                BENABLED=false
                            fi

                            if [ "${BENABLED}" = "true" ]; then

                                echo >&2 ""
                                echo >&2 "************************************************************************************************************************************"
                                echo >&2 ""
                                echo >&2 "-- ParaMonte - invoking: "
                                echo >&2 ""
                                echo >&2 "                          buildParaMonte.sh \ "
                                if ! [ "${compiler_suite_flag}" = "" ]; then
                                echo >&2 "                          ${compiler_suite_flag} \ "
                                fi
                                echo >&2 "                          ${build_type_flag} \ "
                                echo >&2 "                          ${lib_type_flag} \ "
                                echo >&2 "                          ${cfi_enabled_flag} \ "
                                echo >&2 "                          ${heap_enabled_flag} \ "
                                echo >&2 "                          ${mpi_enabled_flag} \ "
                                echo >&2 "                          ${caftype_flag} \ "
                                echo >&2 "                          ${test_enabled_flag} \ "
                                echo >&2 "                          ${exam_enabled_flag} \ "
                                if ! [ "${yes_to_all_flag}" = "" ]; then
                                echo >&2 "                          ${yes_to_all_flag} \ "
                                fi
                                if ! [ "${fresh_flag}" = "" ]; then
                                echo >&2 "                          ${fresh_flag} \ "
                                fi
                                if ! [ "${gcc_bootstrap_flag}" = "" ]; then
                                echo >&2 "                          ${gcc_bootstrap_flag} \ "
                                fi
                                if ! [ "${fortran_flag}" = "" ]; then
                                echo >&2 "                          ${fortran_flag} \ "
                                fi
                                if ! [ "${mpiexec_flag}" = "" ]; then
                                echo >&2 "                          ${mpiexec_flag} \ "
                                fi
                                if ! [ "${num_images_flag}" = "" ]; then
                                echo >&2 "                          ${num_images_flag} \ "
                                fi
                                echo >&2 "                          --clean"
                                echo >&2 ""
                                echo >&2 "************************************************************************************************************************************"
                                echo >&2 ""

                                (cd ${ParaMonte_ROOT_DIR} && \
                                chmod +x ./buildParaMonte.sh && \
                                ./buildParaMonte.sh \
                                ${compiler_suite_flag} \
                                ${build_type_flag} \
                                ${lib_type_flag} \
                                ${cfi_enabled_flag} \
                                ${heap_enabled_flag} \
                                ${mpi_enabled_flag} \
                                ${caftype_flag} \
                                ${test_enabled_flag} \
                                ${exam_enabled_flag} \
                                ${yes_to_all_flag} \
                                ${fresh_flag} \
                                ${gcc_bootstrap_flag} \
                                ${fortran_flag} \
                                ${mpiexec_flag} \
                                ${num_images_flag} \
                                )

                                fresh_flag=""

                            fi

                        done

                    done

                done

            done

        done

    done

done

echo >&2 ""
echo >&2 "-- ParaMonte - all build files are stored at ${ParaMonte_ROOT_DIR}/build/"
echo >&2 "-- ParaMonte - the library files are ready to use at ${ParaMonte_ROOT_DIR}/bin/"
echo >&2 ""
echo >&2 "-- ParaMonte - mission accomplished."
echo >&2 ""
