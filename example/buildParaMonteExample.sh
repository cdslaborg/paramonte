#!/bin/sh
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

# NOTE: This is not a standalone build-script. It must only be called by buildParaMonte.sh script in the root directory of the project.

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# build ParaMonte library example objects and executable
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

echo >&2 
echo >&2 "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
echo >&2 "::::                                                                                                                            ::::"
echo >&2 "                                                  ParaMonte Library Examples Build                                                  "
echo >&2 "::::                                                                                                                            ::::"
echo >&2 "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
echo >&2 

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# setup examples' interface language
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

abort() {
    echo >&2
    echo >&2 "-- ParaMonteExample${LANG_NAME} - Copy action failed. Please resolve the error. Gracefully exiting..."
    echo >&2
    exit 1
}

unset LANG_NAME
LANG_IS_C=false
LANG_IS_CPP=false
LANG_IS_MATLAB=false
LANG_IS_Python=false
LANG_IS_Fortran=false
LANG_IS_DYNAMIC=false
LANG_IS_COMPILED=false

if [ "${INTERFACE_LANGUAGE}" = "python" ]; then
    LANG_IS_DYNAMIC=true
    LANG_IS_Python=true
    LANG_FILE_EXT=py
    LANG_NAME=Python
fi

if [ "${INTERFACE_LANGUAGE}" = "matlab" ]; then
    LANG_IS_DYNAMIC=true
    LANG_IS_MATLAB=true
    LANG_FILE_EXT=m
    LANG_NAME=MATLAB
fi

if [ "${INTERFACE_LANGUAGE}" = "fortran" ]; then
    LANG_IS_COMPILED=true
    LANG_IS_Fortran=true
    LANG_FILE_EXT=f90
    LANG_NAME=Fortran
fi

if [ "${INTERFACE_LANGUAGE}" = "c++" ]; then
    LANG_IS_COMPILED=true
    LANG_FILE_EXT=cpp
    LANG_IS_CPP=true
    LANG_NAME=C++
fi

if [ "${INTERFACE_LANGUAGE}" = "c" ]; then
    LANG_IS_COMPILED=true
    LANG_FILE_EXT=c
    LANG_IS_C=true
    LANG_NAME=C
fi

if [ -z ${LANG_NAME+x} ]; then
    echo >&2
    echo >&2 "-- ParaMonteExample - Fatal Error: unrecognized or no language specified. exiting..."
    echo >&2
    exit 1
else
    LANG_ABBR="${LANG_NAME//+/p}"
#    if [ "${INTERFACE_LANGUAGE}" = "c++" ]; then
#        LANG_ABBR="cpp"
#    else
#        LANG_ABBR="${INTERFACE_LANGUAGE}"
#    fi
fi

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# build examples
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# set and make example directories

ParaMonteExample_BLD_DIR="${ParaMonte_BLD_DIR}/example"
ParaMonteInterface_SRC_DIR_CURRENT="${ParaMonteInterface_SRC_DIR}/${LANG_NAME}"

# select examples to build

EXAM_LIST="mvn"
export EXAM_LIST

echo >&2
echo >&2 "-- ParaMonteExample${LANG_NAME} - The ParaMonte library example interface language: ${INTERFACE_LANGUAGE}"
echo >&2

# make example build directory

ParaMonteExample_BLD_DIR=${ParaMonte_BLD_DIR}/example
echo >&2 "-- ParaMonteExample${LANG_NAME} - ParaMonte examples root directory: ${ParaMonteExample_BLD_DIR}"
if [[ -d "${ParaMonteExample_BLD_DIR}" ]]; then
    echo >&2 "-- ParaMonteExample${LANG_NAME} - ParaMonte  examples root directory already exists. skipping..."
else
    echo >&2 "-- ParaMonteExample${LANG_NAME} - generating ParaMonte examples root directory..."
    mkdir "${ParaMonteExample_BLD_DIR}/"
fi
export ParaMonteExample_BLD_DIR

echo >&2 
echo >&2 "-- ParaMonteExample${LANG_NAME} - generating ParaMonte examples in ${LANG_NAME} language..."
echo >&2 "-- ParaMonteExample${LANG_NAME} - The ParaMonte ${LANG_NAME} examples directory: ${ParaMonteExample_BLD_DIR}"

for EXAM_NAME in $EXAM_LIST
do

    echo >&2
    echo >&2 "-- ParaMonteExample${LANG_NAME} - EXAM_NAME=${EXAM_NAME}"
    echo >&2

    ParaMonteExample_BLD_DIR_CURRENT="${ParaMonteExample_BLD_DIR}/${EXAM_NAME}"
    if [[ -d "${ParaMonteExample_BLD_DIR_CURRENT}" ]]; then
        echo >&2 "-- ParaMonteExample${LANG_NAME} - ParaMonte ${EXAM_NAME} example build directory already exists. deleting the old contents..."
        rm -rf "${ParaMonteExample_BLD_DIR_CURRENT}"
    fi
    echo >&2 "-- ParaMonteExample${LANG_NAME} - generating ParaMonte ${EXAM_NAME} example build directory..."
    mkdir -p "${ParaMonteExample_BLD_DIR_CURRENT}/"

    # The ParaMonte library kernel files

    ParaMonteExample_LIB_DIR_CURRENT="${ParaMonteExample_BLD_DIR_CURRENT}"
    if [ "${LANG_IS_Python}" = "true" ]; then ParaMonteExample_LIB_DIR_CURRENT="${ParaMonteExample_LIB_DIR_CURRENT}/paramonte"; fi
    if [ "${LANG_IS_MATLAB}" = "true" ]; then ParaMonteExample_LIB_DIR_CURRENT="${ParaMonteExample_LIB_DIR_CURRENT}/paramonte/lib"; fi

    if [[ -d "${ParaMonteExample_LIB_DIR_CURRENT}" ]]; then
        echo >&2 "-- ParaMonteExample${LANG_NAME} - ParaMonte ${EXAM_NAME} example library directory already exists. skipping..."
    else
        echo >&2 "-- ParaMonteExample${LANG_NAME} - generating ParaMonte ${EXAM_NAME} example library directory..."
        mkdir -p "${ParaMonteExample_LIB_DIR_CURRENT}/"
    fi

    echo >&2
    echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library files..."
    echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${PMLIB_FULL_PATH}"
    echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_LIB_DIR_CURRENT}/"
    cp -R "${ParaMonte_LIB_DIR}/"libparamonte_* "${ParaMonteExample_LIB_DIR_CURRENT}/" || abort

    # The ParaMonte library dll dependency files

    if ! [ "${PMCS}" = "intel" ]; then
        if ! [ -z ${Fortran_COMPILER_PATH+x} ]; then
            Fortran_COMPILER_DIR=$(dirname "${Fortran_COMPILER_PATH}")
            Fortran_COMPILER_ROOT_DIR="${Fortran_COMPILER_DIR}"/..
            for Fortran_COMPILER_LIB_SUBDIR in "lib64"
            do
                Fortran_COMPILER_LIB_DIR="${Fortran_COMPILER_ROOT_DIR}"/"${Fortran_COMPILER_LIB_SUBDIR}"
                if [ -d "${Fortran_COMPILER_LIB_DIR}" ]; then
                    echo >&2
                    echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library dll dependency files..."
                    echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${PMLIB_FULL_PATH}"
                    echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_LIB_DIR_CURRENT}/"
                    cp -R "${Fortran_COMPILER_LIB_DIR}/"libgfortran.so.* "${ParaMonteExample_LIB_DIR_CURRENT}/" || abort
                    break
                fi
            done
        else
            echo >&2
            echo >&2 "-- ParaMonteExample${LANG_NAME} - WARNING: the ParaMonte library dll dependency files could not be found."
            echo >&2
        fi
    fi

    # The ParaMonte library example required files

    echo >&2
    echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library ${EXAM_NAME} example build files in ${LANG_NAME} language..."

    # The ParaMonte library license file

    echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library license file..."
    echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonte_ROOT_DIR}/LICENSE.md"
    echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/LICENSE.md"
    cp "${ParaMonte_ROOT_DIR}/LICENSE.md" "${ParaMonteExample_BLD_DIR_CURRENT}/LICENSE.md" || abort

    # The ParaMonte library README.md file

    echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library README.md file..."
    echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterface_SRC_DIR_CURRENT}/README.md"
    echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/README.md"
    cp "${ParaMonteInterface_SRC_DIR_CURRENT}/README.md" "${ParaMonteExample_BLD_DIR_CURRENT}/README.md" || abort

    if [ "${LANG_IS_COMPILED}" = "true" ]; then

        # The ParaMonte library CHANGES.md file

        echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library CHANGES.md file..."
        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonte_ROOT_DIR}/CHANGES.md"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/CHANGES.md"
        cp "${ParaMonte_ROOT_DIR}/CHANGES.md" "${ParaMonteExample_BLD_DIR_CURRENT}/CHANGES.md" || abort

        # The ParaMonte library build script

        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteExample_SRC_DIR}"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}"
        cp "${ParaMonteExample_SRC_DIR}/build.sh" "${ParaMonteExample_BLD_DIR_CURRENT}/" || abort
        chmod +x ${ParaMonteExample_BLD_DIR_CURRENT}/build.sh
        if [[ -f "${SETUP_FILE_PATH}" ]]; then
            cp "${SETUP_FILE_PATH}" "${ParaMonteExample_BLD_DIR_CURRENT}/" || abort
        fi

        # The ParaMonte library example header/module files

        echo >&2
        echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library C header file paramonte.h..."
        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterface_SRC_DIR_CURRENT}"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/"
        cp ${ParaMonteInterface_SRC_DIR_CURRENT}/* ${ParaMonteExample_BLD_DIR_CURRENT}/ || abort

        if [ "${LANG_NAME}" = "Fortran" ]; then
            echo >&2
            echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library Fortran module file paramonte.mod..."
            echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonte_MOD_DIR}"
            echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}"
            cp ${ParaMonte_MOD_DIR}/paradram_mod.mod ${ParaMonteExample_BLD_DIR_CURRENT}/ || abort
        fi

    fi

    if [ "${LANG_IS_DYNAMIC}" = "true" ]; then

        # The ParaMonte library CHANGES.md file

        echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library CHANGES.md file..."
        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterface_SRC_DIR_CURRENT}/CHANGES.md"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/CHANGES.md"
        cp "${ParaMonteInterface_SRC_DIR_CURRENT}/CHANGES.md" "${ParaMonteExample_BLD_DIR_CURRENT}/CHANGES.md" || abort

        # The ParaMonte library interface files

        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterface_SRC_DIR_CURRENT}/paramonte"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/"
        cp -R "${ParaMonteInterface_SRC_DIR_CURRENT}/paramonte" "${ParaMonteExample_BLD_DIR_CURRENT}/" || abort

        # The ParaMonte library auxiliary file (this must be done to generate the auxil folder)

        echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library auxiliary files"
        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterface_SRC_DIR}/auxil"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/"
        cp -R "${ParaMonteInterface_SRC_DIR}/auxil" "${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/" || abort
        echo >&2

        # The ParaMonte kernel version file (this must appear only after the above)

        echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library kernel version file..."
        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonte_ROOT_DIR}/.VERSION"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/auxil/.VERSION_KERNEL"
        cp "${ParaMonte_ROOT_DIR}/.VERSION" "${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/auxil/.VERSION_KERNEL" || abort

        # The ParaMonte interface version file (this must appear only after the above)

        echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library interface version file..."
        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterface_SRC_DIR}/${LANG_NAME}/.VERSION"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/auxil/.VERSION_INTERFACE"
        cp "${ParaMonteInterface_SRC_DIR}/${LANG_NAME}/.VERSION" "${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/auxil/.VERSION_INTERFACE" || abort

        if [ "${LANG_IS_Python}" = "true" ]; then

            # PyPI build - The ParaMonte library Python setup files

            echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library Python setup files..."
            echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterface_SRC_DIR_CURRENT}/setup/*"
            echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/"
            cp ${ParaMonteInterface_SRC_DIR_CURRENT}/setup/* "${ParaMonteExample_BLD_DIR_CURRENT}/" || abort

        fi

    fi

    # The ParaMonte library example input files

    ParaMonteExample_INP_DIR_CURRENT="${ParaMonteExample_SRC_DIR}/${EXAM_NAME}/input"
    echo >&2
    echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library ${EXAM_NAME} example input files in ${LANG_NAME} language..."
    echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteExample_INP_DIR_CURRENT}"
    echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}"
    cp -R ${ParaMonteExample_INP_DIR_CURRENT}/* ${ParaMonteExample_BLD_DIR_CURRENT}/ || abort

    # The ParaMonte library example source files

    ParaMonteExample_SRC_DIR_CURRENT="${ParaMonteExample_SRC_DIR}/${EXAM_NAME}/${LANG_NAME}"

    echo >&2
    echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library ${EXAM_NAME} example source files in ${LANG_NAME} language..."

    # copy the example files

    mainFileName=main.${LANG_FILE_EXT}
    if [ "${LANG_IS_DYNAMIC}" = "true" ] && [ "${MPI_ENABLED}" = "true" ]; then mainFileName=main_mpi.${LANG_FILE_EXT}; fi

    echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteExample_SRC_DIR}/main.${LANG_FILE_EXT}"
    echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/"
    cp "${ParaMonteExample_SRC_DIR}/${mainFileName}" "${ParaMonteExample_BLD_DIR_CURRENT}/" || abort

    echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteExample_SRC_DIR_CURRENT}/*"
    echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}"
    cp "${ParaMonteExample_SRC_DIR_CURRENT}"/* "${ParaMonteExample_BLD_DIR_CURRENT}/" || abort

done

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: copy the first example to the bin directory
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

ParaMonteExample_BLD_DIR_CURRENT="${ParaMonteExample_BLD_DIR}/mvn"

if ! [ -d "${ParaMonte_BIN_DIR}" ]; then mkdir "${ParaMonte_BIN_DIR}/"; fi

ParaMonteExample_BIN_DIR_CURRENT="${ParaMonte_BIN_DIR}/${PMLIB_BASE_NAME}"
if [ "${LANG_IS_DYNAMIC}" = "true" ]; then ParaMonteExample_BIN_DIR_CURRENT="${ParaMonte_BIN_DIR}/libparamonte_${LANG_ABBR}"; fi
if ! [ -d "${ParaMonteExample_BIN_DIR_CURRENT}" ]; then mkdir "${ParaMonteExample_BIN_DIR_CURRENT}/"; fi

echo >&2 "-- ParaMonteExample${LANG_NAME} - The ParaMonte ${LANG_NAME} library binary directory: ${ParaMonteExample_BIN_DIR_CURRENT}"

if [ -d "${ParaMonteExample_BIN_DIR_CURRENT}" ]; then
    echo >&2 "-- ParaMonteExample${LANG_NAME} - ParaMonte binary/library directory already exists. skipping..."
else
    echo >&2 "-- ParaMonteExample${LANG_NAME} - generating ParaMonte binary/library directory..."
    mkdir -p "${ParaMonteExample_BIN_DIR_CURRENT}"
fi

echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library files to the bin folder..."
echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteExample_BLD_DIR_CURRENT}"
echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BIN_DIR_CURRENT}"
cp -R "${ParaMonteExample_BLD_DIR_CURRENT}"/* "${ParaMonteExample_BIN_DIR_CURRENT}" || abort

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# The ParaMonte library example build and run if requested
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

for EXAM_NAME in $EXAM_LIST
do

    echo >&2
    if [ "${ParaMonteExample_RUN_ENABLED}" = "true" ] && [ "${LANG_IS_COMPILED}" = "true" ]; then
        ParaMonteExample_BLD_DIR_CURRENT="${ParaMonteExample_BLD_DIR}/${EXAM_NAME}"
        (cd ${ParaMonteExample_BLD_DIR_CURRENT} && ./build.sh)
        if [ $? -eq 0 ]; then
            echo >&2 "-- ParaMonteExample - ParaMonte example build successful."
            echo >&2
        else
            echo >&2 "-- ParaMonteExample - Failed to build the example. exiting..."
            echo >&2
            exit 1
        fi
        (cd ${ParaMonteExample_BLD_DIR_CURRENT} && ./run.sh)
        if [ $? -eq 0 ]; then
            echo >&2 "-- ParaMonteExample - ParaMonte example run successful."
            echo >&2
        else
            echo >&2 "-- ParaMonteExample - Failed to run the example. exiting..."
            echo >&2
            exit 1
        fi
    fi

done
