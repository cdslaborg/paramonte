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
#   NOTE: Do not change the contents of this file unless you know what the consequences are.
#   This is the Bash script file that builds objects, dynamic libraries, 
#   as well as the test and example binaries of the ParaMonte library on non-Windows systems.
#   Upon invocation of this file from a Bash command-line interface, 
#   this script will parse the user-provided flags and their values 
#   to build the ParaMonte library.
#   to redirect output to the external file install.sh.out, try:
#
#       install.sh >install.sh.out 2>&1
#
#   to redirect output to the external file install.sh.out and run the installation in background, try:
#
#       install.sh >install.sh.out 2>&1 &
#       jobs; disown

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
unset PARALLELISM_LIST
unset MEMORY_LIST
unset ParaMonteTest_RUN_ENABLED
unset ParaMonteExample_RUN_ENABLED
unset Fortran_COMPILER_PATH
unset MPIEXEC_PATH

PMCS_LIST="none"
fresh_flag=""
yes_to_all_flag=""
gcc_bootstrap_flag=""
FOR_COARRAY_NUM_IMAGES=3
MatDRAM_ENABLED="false"
dryrun_flag=""

while [ "$1" != "" ]; do
    case $1 in
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
        -d | --dryrun )         dryrun_flag="--dryrun"
                                ;;
        -y | --yes-to-all )     yes_to_all_flag="--yes-to-all"
                                ;;
        -B | --bootstrap )      gcc_bootstrap_flag="--bootstrap"
                                ;;
#       -a | --matdram )        shift
#                               MatDRAM_ENABLED="true"
#                               ;;
        -n | --nproc )          shift
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
# determine whether to build MatDRAM or not. NOTE: If true, all other builds will be disabled. NOT IMPLEMENTED YET. NOT NEEDED.
####################################################################################################################################

#export MatDRAM_ENABLED
#if [ "${MatDRAM_ENABLED}" = "true" ]; then
#fi

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
    echo >&2 "-- ParaMonte - WARNING: conflicting flag values detected."
    echo >&2 "-- ParaMonte - WARNING: $1"
    echo >&2 "-- ParaMonte - WARNING: The requested build configuration will ignored."
    echo >&2 "-- ParaMonte - skipping..."
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
            && ($LANG != "c++" && $LANG != "C++")
            && $LANG != [fF][oO][rR][tT][rR][aA][nN]
            && $LANG != [mM][aA][tT][lL][aA][bB]
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
            reportBadValue "-b or --build" $BTYPE
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
            reportBadValue "-l or --lib" $LTYPE
        fi
    done
    if [ "${LTYPE_LIST}" = "" ]; then
        unset LTYPE_LIST
    else
        LTYPE_LIST="$(getLowerCase $LTYPE_LIST)"
    fi
fi

if ! [ -z ${PARALLELISM_LIST+x} ]; then
    for PARALLELISM in $PARALLELISM_LIST; do
        if  [[ $PARALLELISM != [nN][oO][nN][eE] 
            && $PARALLELISM != [cC][aA][fF][sS][iI][nN][gG][lL][eE] 
            && $PARALLELISM != [cC][aA][fF][sS][hH][aA][rR][eE][dD] 
            && $PARALLELISM != [cC][aA][fF][dD][iI][sS][tT][rR][iI][bB][uU][tT][eE][dD] 
            && $PARALLELISM != [mM][pP][iI] 
            ]]; then
            reportBadValue "-p or --parallelism" $PARALLELISM
        fi
    done
    if [ "${PARALLELISM_LIST}" = "" ]; then
        unset PARALLELISM_LIST
    else
        PARALLELISM_LIST="$(getLowerCase $PARALLELISM_LIST)"
    fi
fi

if ! [ -z ${MEMORY_LIST+x} ]; then
    for MEMORY in $MEMORY_LIST; do
        if  [[ $MEMORY != [hH][eE][aA][pP] 
            && $MEMORY != [sS][tT][aA][cC][kK] ]]; then
            reportBadValue "-m or --memory" $MEMORY
        fi
    done
    if [ "${MEMORY_LIST}" = "" ]; then
        unset MEMORY_LIST
    else
        MEMORY_LIST="$(getLowerCase $MEMORY_LIST)"
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

nproc_flag=""
if ! [ -z ${FOR_COARRAY_NUM_IMAGES+x} ]; then
    isNumericValue="$(isnumeric ${FOR_COARRAY_NUM_IMAGES})"
    if [ "${isNumericValue}" = "true" ]; then
        nproc_flag="--nproc ${FOR_COARRAY_NUM_IMAGES}"
    else
        reportBadValue "-n or --nproc" $FOR_COARRAY_NUM_IMAGES "The number of processors must be a positive integer."
    fi
fi

####################################################################################################################################
# verify arguments consistencies
####################################################################################################################################

if ! [ -z ${PARALLELISM_LIST+x} ]; then
    for PARALLELISM in $PARALLELISM_LIST; do
        if  [[ "${PARALLELISM}" =~ .*"caf".* ]]; then
            for LANG in $LANG_LIST; do
                if  [ "${LANG}" != "fortran" ]; then
                    reportConflict "Coarray Fortran parallelism cannot be used to build the ParaMonte library for the ${LANG} language."
                fi
            done
            if  [[ "${PARALLELISM}" =~ .*"mpi".* ]]; then
                reportConflict "Coarray Fortran parallelism cannot be mixed with MPI."
            fi
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
        if  [ "${LANG}" = "matlab" ] || [ "${LANG}" = "python" ]; then
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
    LANG_LIST="c c++ fortran matlab python"
fi
# # CFI_ENABLED_LIST=""
# C_IS_MISSING=true
# Fortran_IS_MISSING=true
# MATLAB_IS_MISSING=true
# Python_IS_MISSING=true
# for LANG in $LANG_LIST; do
    # if [ "${LANG}" = "c" ]; then
        # if [ "${C_IS_MISSING}" = "true" ]; then
            # # CFI_ENABLED_LIST="${CFI_ENABLED_LIST} true"
            # C_IS_MISSING=false
        # fi
    # fi
    # if [ "${LANG}" = "fortran" ]; then if [ "${Fortran_IS_MISSING}" = "true" ]; then Fortran_IS_MISSING=false; fi; fi
    # if [ "${LANG}" = "matlab" ]; then if [ "${MATLAB_IS_MISSING}" = "true" ]; then MATLAB_IS_MISSING=false; fi; fi
    # if [ "${LANG}" = "python" ]; then if [ "${Python_IS_MISSING}" = "true" ]; then Python_IS_MISSING=false; fi; fi
# done

if [ -z ${PMCS_LIST+x} ]; then
    PMCS_LIST="none"
fi
if [ -z ${BTYPE_LIST+x} ]; then
    BTYPE_LIST="release testing debug"
fi
if [ -z ${LTYPE_LIST+x} ]; then
    LTYPE_LIST="static dynamic"
fi
if [ -z ${PARALLELISM_LIST+x} ]; then
    PARALLELISM_LIST="none mpi cafsingle cafshared cafdistributed"
fi
if [ -z ${MEMORY_LIST+x} ]; then
    MEMORY_LIST="stack heap"
fi
if [ -z ${ParaMonteTest_RUN_ENABLED+x} ]; then
    ParaMonteTest_RUN_ENABLED="true"
fi
if [ -z ${ParaMonteExample_RUN_ENABLED+x} ]; then
    ParaMonteExample_RUN_ENABLED="true"
fi

if [ "${LANG_LIST}" = "matlab" ]; then
    MEMORY_LIST="heap"
    LTYPE_LIST="dynamic"
    if [ -z ${PARALLELISM_LIST+x} ]; then PARALLELISM_LIST="none mpi"; fi
fi

if [ "${LANG_LIST}" = "python" ]; then
    MEMORY_LIST="heap"
    LTYPE_LIST="dynamic"
    if [ -z ${PARALLELISM_LIST+x} ]; then PARALLELISM_LIST="none mpi"; fi
fi

for PMCS in $PMCS_LIST; do

    #for CFI_ENABLED in $CFI_ENABLED_LIST; do
    for INTERFACE_LANGUAGE in $LANG_LIST; do

        for BTYPE in $BTYPE_LIST; do

            for LTYPE in $LTYPE_LIST; do

                for MEMORY in $MEMORY_LIST; do

                    for PARALLELISM in $PARALLELISM_LIST; do

                        BENABLED=true

                        interface_language_flag="--lang ${INTERFACE_LANGUAGE}"
                        if [ "${INTERFACE_LANGUAGE}" = "fortran" ]; then
                            CFI_ENABLED="false"
                        else
                            CFI_ENABLED="true"
                        fi
                        cfi_enabled_flag="--cfi_enabled ${CFI_ENABLED}"

                        test_enabled_flag="--test_enabled ${ParaMonteTest_RUN_ENABLED}"
                        exam_enabled_flag="--exam_enabled ${ParaMonteExample_RUN_ENABLED}"

                        if [ "${PMCS}" = "none" ]; then
                           compiler_suite_flag=""
                        else
                           compiler_suite_flag="--compiler_suite ${PMCS}"
                        fi

                        caftype_flag="--caf none"
                        mpi_enabled_flag="--mpi_enabled false"
                        if [[ "${PARALLELISM}" =~ .*"caf".* ]]; then
                            caftype_flag="--caf ${PARALLELISM:3}"
                        elif [[ "${PARALLELISM}" =~ .*"mpi".* ]]; then
                            mpi_enabled_flag="--mpi_enabled true"
                        fi
                        lib_flag="--lib ${LTYPE}"
                        build_flag="--build ${BTYPE}"
                        heap_enabled_flag="--heap_enabled true"
                        if [ "${MEMORY}" = "stack" ]; then
                            heap_enabled_flag="--heap_enabled false"
                        fi

                        # verify no conflict

                        if [[ "${PARALLELISM}" =~ .*"caf".* ]]; then
                            if [ "${CFI_ENABLED}" = "true" ] || [ "${LTYPE}" = "dynamic" ]; then
                                BENABLED=false
                            fi
                        fi

                        if [ "${LTYPE}" = "dynamic" ]; then
                            if [ "${MEMORY}" = "stack" ]; then BENABLED=false; fi
                        else
                            if [ "${INTERFACE_LANGUAGE}" = "matlab" ] || [ "${INTERFACE_LANGUAGE}" = "python" ]; then BENABLED=false; fi
                        fi

                        if [ "${BENABLED}" = "true" ]; then

                            echo >&2 ""
                            echo >&2 "************************************************************************************************************************************"
                            echo >&2 ""
                            echo >&2 "-- ParaMonte - invoking: "
                            echo >&2 ""
                            echo >&2 "                          buildParaMonte.sh \ "
                            echo >&2 "                          ${interface_language_flag} \ "
                            if ! [ "${compiler_suite_flag}" = "" ]; then
                            echo >&2 "                          ${compiler_suite_flag} \ "
                            fi
                            echo >&2 "                          ${build_flag} \ "
                            echo >&2 "                          ${lib_flag} \ "
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
                            if ! [ "${dryrun_flag}" = "" ]; then
                            echo >&2 "                          ${dryrun_flag} \ "
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
                            if ! [ "${nproc_flag}" = "" ]; then
                            echo >&2 "                          ${nproc_flag} \ "
                            fi
                            echo >&2 "                          --clean"
                            echo >&2 ""
                            echo >&2 "************************************************************************************************************************************"
                            echo >&2 ""

                            (cd ${ParaMonte_ROOT_DIR} && \
                            chmod +x ./buildParaMonte.sh && \
                            ./buildParaMonte.sh \
                            ${interface_language_flag} \
                            ${compiler_suite_flag} \
                            ${build_flag} \
                            ${lib_flag} \
                            ${cfi_enabled_flag} \
                            ${heap_enabled_flag} \
                            ${mpi_enabled_flag} \
                            ${caftype_flag} \
                            ${test_enabled_flag} \
                            ${exam_enabled_flag} \
                            ${yes_to_all_flag} \
                            ${fresh_flag} \
                            ${dryrun_flag} \
                            ${gcc_bootstrap_flag} \
                            ${fortran_flag} \
                            ${mpiexec_flag} \
                            ${nproc_flag} \
                            ) || {
                            echo >&2 ""
                            echo >&2 "-- ParaMonte "
                            echo >&2 "-- ParaMonte - Fatal Error: The ParaMonte library build failed for the following configuration:"
                            echo >&2 "-- ParaMonte - "
                            echo >&2 "-- ParaMonte -               language: ${INTERFACE_LANGUAGE}"
                            echo >&2 "-- ParaMonte -             build type: ${BTYPE}"
                            echo >&2 "-- ParaMonte -           library type: ${LTYPE}"
                            echo >&2 "-- ParaMonte -      memory allocation: ${MEMORY}"
                            echo >&2 "-- ParaMonte -            parallelism: ${PARALLELISM}"
                            echo >&2 "-- ParaMonte - "
                            echo >&2 "-- ParaMonte - If you cannot identify the cause of the failure, please report this error at: "
                            echo >&2 "-- ParaMonte - "
                            echo >&2 "-- ParaMonte -     https://github.com/cdslaborg/paramonte/issues"
                            echo >&2 "-- ParaMonte - "
                            echo >&2 "-- ParaMonte - gracefully exiting..."
                            echo >&2 ""
                            cd "${ParaMonte_ROOT_DIR}"
                            #return
                            exit 1
                            }

                            fresh_flag=""

                        fi

                    done

                done

            done

        done

        # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        # :: if MATLAB, generate MatDRAM
        # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        #if [ "${INTERFACE_LANGUAGE}" = "matlab" ]; then
        #
        #    echo >&2 "-- ParaMonte - Generating the ParaMonte::MatDRAM library..."
        #    echo >&2 ""
        #
        #    MatDRAM_ORIGIN_PATH=./bin/MATLAB
        #    MatDRAM_DESTINATION_PATH=./bin/MatDRAM
        #    if ! [ -d "${MatDRAM_DESTINATION_PATH}" ]; then
        #        mkdir -p "${MatDRAM_DESTINATION_PATH}"
        #    fi
        #    echo >&2 "-- ParaMonte - copying the MatDRAM library files..."
        #    echo >&2 "-- ParaMonte - from: ${MatDRAM_ORIGIN_PATH}"
        #    echo >&2 "-- ParaMonte -   to: ${MatDRAM_DESTINATION_PATH}"
        #    cp -frp "${MatDRAM_ORIGIN_PATH}" -T "${MatDRAM_DESTINATION_PATH}"
        #
        #    # delete the binary files
        #
        #    rm -rf "${MatDRAM_DESTINATION_PATH}/paramonte/lib"
        #
        #    # delete the mpi example file
        #
        #    rm -rf "${MatDRAM_DESTINATION_PATH}/main_mpi.m"
        #
        #fi

    done

done

echo >&2 ""
echo >&2 "-- ParaMonte - all build files are stored at ${ParaMonte_ROOT_DIR}/build/"
echo >&2 "-- ParaMonte - the library files are ready to use at ${ParaMonte_ROOT_DIR}/bin/"
echo >&2 ""
echo >&2 "-- ParaMonte - mission accomplished."
echo >&2 ""
