#!/bin/sh
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

unset LANG_NAME
LANG_IS_C=false
LANG_IS_MATLAB=false
LANG_IS_Python=false
LANG_IS_Fortran=false
LANG_IS_DYNAMIC=false
LANG_IS_FortranC=false

if [ "${INTERFACE_LANGUAGE}" = "matlab" ]; then
    LANG_IS_DYNAMIC=true
    LANG_IS_MATLAB=true
    LANG_FILE_EXT=m
    LANG_NAME=MATLAB
fi

if [ "${INTERFACE_LANGUAGE}" = "python" ]; then
    LANG_IS_DYNAMIC=true
    LANG_IS_Python=true
    LANG_FILE_EXT=py
    LANG_NAME=Python
fi

if [ "${INTERFACE_LANGUAGE}" = "fortran" ]; then
    LANG_IS_FortranC=true
    LANG_IS_Fortran=true
    LANG_FILE_EXT=f90
    LANG_NAME=Fortran
fi

if [ "${INTERFACE_LANGUAGE}" = "c" ]; then
    LANG_IS_FortranC=true
    LANG_FILE_EXT=c
    LANG_IS_C=true
    LANG_NAME=C
fi

if [ -z ${LANG_NAME+x} ]; then
    echo >&2
    echo >&2 "-- ParaMonteExample - Fatal Error: unrecognized or no language specified. exiting..."
    echo >&2
    exit 1
fi

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# build examples
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# set and make example directories

set ParaMonteExample_BLD_DIR=${ParaMonte_BLD_DIR}/example

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
    cp -R "${ParaMonte_LIB_DIR}/"* "${ParaMonteExample_LIB_DIR_CURRENT}/"

    # The ParaMonte library example required files

    echo >&2
    echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library ${EXAM_NAME} example build files in ${LANG_NAME} language..."

    if [ "${LANG_IS_FortranC}" = "true" ]; then

        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteExample_SRC_DIR}"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}"
        cp "${ParaMonteExample_SRC_DIR}/build.sh" "${ParaMonteExample_BLD_DIR_CURRENT}/"
        chmod +x ${ParaMonteExample_BLD_DIR_CURRENT}/build.sh
        if [[ -f "${SETUP_FILE_PATH}" ]]; then
            cp "${SETUP_FILE_PATH}" "${ParaMonteExample_BLD_DIR_CURRENT}/"
        fi

        # The ParaMonte library example header/module files

        if [ "${LANG_NAME}" = "Fortran" ]; then
            echo >&2
            echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library Fortran module file paramonte.mod..."
            echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonte_MOD_DIR}"
            echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}"
            cp ${ParaMonte_MOD_DIR}/paradram_mod.mod ${ParaMonteExample_BLD_DIR_CURRENT}/
            echo >&2
            echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library Fortran module file paramonte.f90..."
            echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterfaceFortran_SRC_DIR}"
            echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/"
            cp ${ParaMonteInterfaceFortran_SRC_DIR}/* ${ParaMonteExample_BLD_DIR_CURRENT}/
        fi
        if [ "${LANG_NAME}" = "C" ]; then
            echo >&2
            echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library C header file paramonte.h..."
            echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterfaceC_SRC_DIR}"
            echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/"
            cp ${ParaMonteInterfaceC_SRC_DIR}/* ${ParaMonteExample_BLD_DIR_CURRENT}/
        fi

    fi

    if [ "${LANG_IS_DYNAMIC}" = "true" ]; then

        # The ParaMonte kernel version file

        echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library kernel version file..."
        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonte_ROOT_DIR}/.VERSION"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/.VERSION_KERNEL"
        cp "${ParaMonte_ROOT_DIR}/.VERSION" "${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/.VERSION_KERNEL"


        if [ "${LANG_IS_MATLAB}" = "true" ]; then

            # The ParaMonte MATLAB library files

            echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterfaceMATLAB_SRC_DIR}/paramonte"
            echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/"
            cp -R "${ParaMonteInterfaceMATLAB_SRC_DIR}/paramonte" "${ParaMonteExample_BLD_DIR_CURRENT}/"

            # The ParaMonte library license files

            echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library license file..."
            echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonte_ROOT_DIR}/LICENSE"
            echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/LICENSE"
            cp "${ParaMonte_ROOT_DIR}/LICENSE" "${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/LICENSE"

            # The ParaMonte library CHANGES.md files

            echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library CHANGES.md file..."
            echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterfaceMATLAB_SRC_DIR}/CHANGES.md"
            echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/CHANGES.md"
            cp "${ParaMonteInterfaceMATLAB_SRC_DIR}/CHANGES.md" "${ParaMonteExample_BLD_DIR_CURRENT}/CHANGES.md"

        fi

        if [ "${LANG_IS_Python}" = "true" ]; then

            # The ParaMonte Python library files

            echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterfacePython_SRC_DIR}/paramonte"
            echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/"
            cp -R "${ParaMonteInterfacePython_SRC_DIR}/paramonte" "${ParaMonteExample_BLD_DIR_CURRENT}/"

            # PyPI build - ParaMonte library Python setup files

            echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library Python setup files..."
            echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterfacePython_SRC_DIR}/setup/*"
            echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/"
            cp ${ParaMonteInterfacePython_SRC_DIR}/setup/* "${ParaMonteExample_BLD_DIR_CURRENT}/"

            # PyPI build - ParaMonte library license files

            echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library license file..."
            echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonte_ROOT_DIR}/LICENSE"
            echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/LICENSE"
            cp "${ParaMonte_ROOT_DIR}/LICENSE" "${ParaMonteExample_BLD_DIR_CURRENT}/LICENSE"

            # PyPI build - ParaMonte library CHANGES.md files

            echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library CHANGES.md file..."
            echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterfacePython_SRC_DIR}/setup/CHANGES.md"
            echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/CHANGES.md"
            cp "${ParaMonteInterfacePython_SRC_DIR}/setup/CHANGES.md" "${ParaMonteExample_BLD_DIR_CURRENT}/CHANGES.md"

        fi

    fi

    # The ParaMonte library example input files

    ParaMonteExample_INP_DIR_CURRENT="${ParaMonteExample_SRC_DIR}/${EXAM_NAME}/input"
    echo >&2
    echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library ${EXAM_NAME} example input files in ${LANG_NAME} language..."
    echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteExample_INP_DIR_CURRENT}"
    echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}"
    cp -R ${ParaMonteExample_INP_DIR_CURRENT}/* ${ParaMonteExample_BLD_DIR_CURRENT}/

    # The ParaMonte library example source files

    ParaMonteExample_SRC_DIR_CURRENT="${ParaMonteExample_SRC_DIR}/${EXAM_NAME}/${LANG_NAME}"

    echo >&2
    echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library ${EXAM_NAME} example source files in ${LANG_NAME} language..."

    if [ "${LANG_IS_FortranC}" = "true" ]; then

        # copy the example files

        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteExample_SRC_DIR_CURRENT}/*"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}"
        cp ${ParaMonteExample_SRC_DIR_CURRENT}/* ${ParaMonteExample_BLD_DIR_CURRENT}/

        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteExample_SRC_DIR}/main.${LANG_FILE_EXT}"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/"
        cp ${ParaMonteExample_SRC_DIR}/main.${LANG_FILE_EXT} ${ParaMonteExample_BLD_DIR_CURRENT}/

    fi

    if [ "${LANG_IS_MATLAB}" = "true" ] || [ "${LANG_IS_Python}" = "true" ]; then

        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteExample_SRC_DIR_CURRENT}/README.md"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/"
        cp ${ParaMonteExample_SRC_DIR_CURRENT}/README.md ${ParaMonteExample_BLD_DIR_CURRENT}/

        scriptFileName=main.${LANG_FILE_EXT}
        if [ "${MPI_ENABLED}" = "true" ]; then scriptFileName=main_mpi.${LANG_FILE_EXT}; fi
        cp ${ParaMonteExample_SRC_DIR_CURRENT}/${scriptFileName} ${ParaMonteExample_BLD_DIR_CURRENT}/

    fi

done

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: copy the first example to the bin directory
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

ParaMonteExample_BLD_DIR_CURRENT="${ParaMonteExample_BLD_DIR}/mvn"

if ! [ -d "${ParaMonte_BIN_DIR}" ]; then mkdir "${ParaMonte_BIN_DIR}/"; fi

ParaMonteExample_BIN_DIR_CURRENT="${ParaMonte_BIN_DIR}/${PMLIB_BASE_NAME}"
if [ "${LANG_IS_DYNAMIC}" = "true" ]; then ParaMonteExample_BIN_DIR_CURRENT="${ParaMonte_BIN_DIR}/${LANG_NAME}"; fi
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
cp -R "${ParaMonteExample_BLD_DIR_CURRENT}"/* "${ParaMonteExample_BIN_DIR_CURRENT}" || {
    echo >&2
    echo >&2 "-- ParaMonteExample${LANG_NAME} - FATAL: copy failed."
    echo >&2
    exit 1
}

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# The ParaMonte library example build and run if requested
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

for EXAM_NAME in $EXAM_LIST
do

    echo >&2
    if [ "${ParaMonteExample_RUN_ENABLED}" = "true" ] && [ "${LANG_IS_FortranC}" = "true" ]; then
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
