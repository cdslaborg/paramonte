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

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# build ParaMonte library example objects and executable
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

echo >&2 
echo >&2 "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
echo >&2 "::::                                                                                                                            ::::"
echo >&2 "                                                  ParaMonte Library Examples Build                                                  "
echo >&2 "::::                                                                                                                            ::::"
echo >&2 "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
echo >&2 

# set and make example directories

set ParaMonteExample_BLD_DIR=${ParaMonte_BLD_DIR}/example

# select examples to build

EXAM_LIST="mvn"
export EXAM_LIST

# select languages for which to build

set EXAM_LANG_LIST=
if [ "${CFI_ENABLED}" = "true" ]; then
    EXAM_LANG_LIST="C"
    if [ "${LTYPE}" = "dynamic" ]; then
        EXAM_LANG_LIST="${EXAM_LANG_LIST} Python"
    fi
else
    EXAM_LANG_LIST="Fortran"
fi
export EXAM_LANG_LIST

echo >&2
echo >&2 "-- ParaMonte - EXAM_LANG_LIST=${EXAM_LANG_LIST}"
echo >&2

# make example build directory

ParaMonteExample_BLD_DIR=${ParaMonte_BLD_DIR}/example
echo >&2 "-- ParaMonte - ParaMonte examples root directory: ${ParaMonteExample_BLD_DIR}"
if [[ -d "${ParaMonteExample_BLD_DIR}" ]]; then
    echo >&2 "-- ParaMonte - ParaMonte  examples root directory already exists. skipping..."
else
    echo >&2 "-- ParaMonte - generating ParaMonte examples root directory..."
    mkdir "${ParaMonteExample_BLD_DIR}/"
fi
export ParaMonteExample_BLD_DIR

# build examples

for EXAM_LANG in $EXAM_LANG_LIST
do

echo >&2
echo >&2 "-- ParaMonte - EXAM_LANG=${EXAM_LANG}"

    ParaMonteExample_LNG_DIR="${ParaMonteExample_BLD_DIR}/$EXAM_LANG"

    echo >&2 
    echo >&2 "-- ParaMonteExample${EXAM_LANG} - generating ParaMonte examples in ${EXAM_LANG} language..."
    echo >&2 "-- ParaMonteExample${EXAM_LANG} - ParaMonte ${EXAM_LANG} examples directory: ${ParaMonteExample_LNG_DIR}"

    IS_Python_LANG=false
    if [ "${EXAM_LANG}" = "Python" ]; then IS_Python_LANG=true; fi

    IS_FortranC_LANG=false
    if [ "${EXAM_LANG}" = "C" ]; then 
        IS_FortranC_LANG=true;
        LANG_FILE_EXT="c"
    fi
    if [ "${EXAM_LANG}" = "Fortran" ]; then
        IS_FortranC_LANG=true;
        LANG_FILE_EXT="f90"
    fi

    for EXAM in $EXAM_LIST
    do

        echo >&2
        echo >&2 "-- ParaMonte - EXAM=${EXAM}"

        # set the binaries directory

        # set ParaMonteExample_BIN_DIR_CURRENT=${ParaMonte_BIN_DIR}/${PMLIB_BASE_NAME}_${ParaMonteVersion}
        ParaMonteExample_BIN_DIR_CURRENT="${ParaMonte_BIN_DIR}/${PMLIB_BASE_NAME}"
        if [ "${IS_Python_LANG}" = "true" ]; then ParaMonteExample_BIN_DIR_CURRENT="${ParaMonte_BIN_DIR}/Python"; fi

        echo >&2
        echo >&2 "-- ParaMonteExample${EXAM_LANG} - current ParaMonte example binary/library directory: ${ParaMonteExample_BIN_DIR_CURRENT}"

        if [ -d "${ParaMonteExample_BIN_DIR_CURRENT}" ]; then
            echo >&2 "-- ParaMonte - ParaMonte binary/library directory already exists. skipping..."
        else
            echo >&2 "-- ParaMonte - generating ParaMonte binary/library directory..."
            mkdir -p "${ParaMonteExample_BIN_DIR_CURRENT}"
        fi

        echo >&2
        ParaMonteExample_BLD_DIR_CURRENT="${ParaMonteExample_LNG_DIR}/${EXAM}"
        if [[ -d "${ParaMonteExample_BLD_DIR_CURRENT}" ]]; then
            echo >&2 "-- ParaMonteExample${EXAM_LANG} - ParaMonte ${EXAM} example build directory already exists. skipping..."
        else
            echo >&2 "-- ParaMonteExample${EXAM_LANG} - generating ParaMonte ${EXAM} example build directory..."
            mkdir -p "${ParaMonteExample_BLD_DIR_CURRENT}/"
        fi

        # ParaMonte library files

        ParaMonteExample_LIB_DIR_CURRENT="${ParaMonteExample_BLD_DIR_CURRENT}"
        if [ "${IS_Python_LANG}" = "true" ]; then ParaMonteExample_LIB_DIR_CURRENT="${ParaMonteExample_LIB_DIR_CURRENT}/paramonte"; fi
        if [[ -d "${ParaMonteExample_LIB_DIR_CURRENT}" ]]; then
            echo >&2 "-- ParaMonteExample${EXAM_LANG} - ParaMonte ${EXAM} example library directory already exists. skipping..."
        else
            echo >&2 "-- ParaMonteExample${EXAM_LANG} - generating ParaMonte ${EXAM} example library directory..."
            mkdir -p "${ParaMonteExample_LIB_DIR_CURRENT}/"
        fi

        echo >&2
        echo >&2 "-- ParaMonteExample${EXAM_LANG} - copying the ParaMonte library files..."
        echo >&2 "-- ParaMonteExample${EXAM_LANG} - from: ${PMLIB_FULL_PATH}"
        echo >&2 "-- ParaMonteExample${EXAM_LANG} -   to: ${ParaMonteExample_LIB_DIR_CURRENT}/"
        cp "${PMLIB_FULL_PATH}" "${ParaMonteExample_LIB_DIR_CURRENT}/"

        if [ "${IS_Python_LANG}" = "true" ]; then
            mkdir -p "${ParaMonteExample_BIN_DIR_CURRENT}/paramonte/" && \
            cp "${PMLIB_FULL_PATH}" "${ParaMonteExample_BIN_DIR_CURRENT}/paramonte/"
        else
            cp "${PMLIB_FULL_PATH}" "${ParaMonteExample_BIN_DIR_CURRENT}/"
        fi

        # ParaMonte library example build files

        echo >&2
        echo >&2 "-- ParaMonteExample${EXAM_LANG} - copying the ParaMonte library ${EXAM} example build files in ${EXAM_LANG} language..."

        if [ "${IS_FortranC_LANG}" = "true" ]; then
            echo >&2 "-- ParaMonteExample${EXAM_LANG} - from: ${ParaMonteExample_SRC_DIR}"
            echo >&2 "-- ParaMonteExample${EXAM_LANG} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}"
            cp "${ParaMonteExample_SRC_DIR}/build.sh" "${ParaMonteExample_BLD_DIR_CURRENT}/"
            cp "${ParaMonteExample_SRC_DIR}/build.sh" "${ParaMonteExample_BIN_DIR_CURRENT}/"
            chmod +x ${ParaMonteExample_BLD_DIR_CURRENT}/build.sh
            chmod +x ${ParaMonteExample_BIN_DIR_CURRENT}/build.sh
            if [[ -f "${SETUP_FILE_PATH}" ]]; then
                cp "${SETUP_FILE_PATH}" "${ParaMonteExample_BLD_DIR_CURRENT}/"
                cp "${SETUP_FILE_PATH}" "${ParaMonteExample_BIN_DIR_CURRENT}/"
            fi
        fi

        if [ "${IS_Python_LANG}" = "true" ]; then

            echo >&2 "-- ParaMonteExample${EXAM_LANG} - from: ${ParaMonteInterfacePython_SRC_DIR}/paramonte"
            echo >&2 "-- ParaMonteExample${EXAM_LANG} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/"
            cp -R "${ParaMonteInterfacePython_SRC_DIR}/paramonte" "${ParaMonteExample_BLD_DIR_CURRENT}/"
            cp -R "${ParaMonteInterfacePython_SRC_DIR}/paramonte" "${ParaMonteExample_BIN_DIR_CURRENT}/"

            # PyPI build - ParaMonte library Python setup files

            echo >&2 "-- ParaMonteExample${EXAM_LANG} - copying the ParaMonte library Python setup files..."
            echo >&2 "-- ParaMonteExample${EXAM_LANG} - from: ${ParaMonteInterfacePython_SRC_DIR}/setup/*"
            echo >&2 "-- ParaMonteExample${EXAM_LANG} -   to: ${ParaMonteExample_BIN_DIR_CURRENT}/"
            cp ${ParaMonteInterfacePython_SRC_DIR}/setup/* "${ParaMonteExample_BIN_DIR_CURRENT}/"

            # PyPI build - ParaMonte library license files

            echo >&2 "-- ParaMonteExample${EXAM_LANG} - copying the ParaMonte library license file..."
            echo >&2 "-- ParaMonteExample${EXAM_LANG} - from: ${ParaMonte_ROOT_DIR}/LICENSE"
            echo >&2 "-- ParaMonteExample${EXAM_LANG} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/LICENSE"
            cp "${ParaMonte_ROOT_DIR}/LICENSE" "${ParaMonteExample_BIN_DIR_CURRENT}/LICENSE"

            # PyPI build - ParaMonte library CHANGES.md files

            echo >&2 "-- ParaMonteExample${EXAM_LANG} - copying the ParaMonte library CHANGES.md file..."
            echo >&2 "-- ParaMonteExample${EXAM_LANG} - from: ${ParaMonteInterfacePython_SRC_DIR}/setup/CHANGES.md"
            echo >&2 "-- ParaMonteExample${EXAM_LANG} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/CHANGES.md"
            cp "${ParaMonteInterfacePython_SRC_DIR}/setup/CHANGES.md" "${ParaMonteExample_BIN_DIR_CURRENT}/CHANGES.md"

        fi

        # ParaMonte library example input files

        ParaMonteExample_INP_DIR_CURRENT="${ParaMonteExample_SRC_DIR}/${EXAM}/input"
        echo >&2
        echo >&2 "-- ParaMonteExample${EXAM_LANG} - copying the ParaMonte library ${EXAM} example input files in ${EXAM_LANG} language..."
        echo >&2 "-- ParaMonteExample${EXAM_LANG} - from: ${ParaMonteExample_INP_DIR_CURRENT}"
        echo >&2 "-- ParaMonteExample${EXAM_LANG} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}"
        cp -R ${ParaMonteExample_INP_DIR_CURRENT}/* ${ParaMonteExample_BLD_DIR_CURRENT}/
        cp -R ${ParaMonteExample_INP_DIR_CURRENT}/* ${ParaMonteExample_BIN_DIR_CURRENT}/

        # ParaMonte library example source files

        ParaMonteExample_SRC_DIR_CURRENT="${ParaMonteExample_SRC_DIR}/${EXAM}/${EXAM_LANG}"

        echo >&2
        echo >&2 "-- ParaMonteExample${EXAM_LANG} - copying the ParaMonte library ${EXAM} example source files in ${EXAM_LANG} language..."

        if [ "${IS_FortranC_LANG}" = "true" ]; then
            echo >&2 "-- ParaMonteExample${EXAM_LANG} - from: ${ParaMonteExample_SRC_DIR_CURRENT}/*"
            echo >&2 "-- ParaMonteExample${EXAM_LANG} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}"
            cp ${ParaMonteExample_SRC_DIR_CURRENT}/* ${ParaMonteExample_BLD_DIR_CURRENT}/
            cp ${ParaMonteExample_SRC_DIR_CURRENT}/* ${ParaMonteExample_BIN_DIR_CURRENT}/
            echo >&2 "-- ParaMonteExample${EXAM_LANG} - from: ${ParaMonteExample_SRC_DIR}/main.${LANG_FILE_EXT}"
            echo >&2 "-- ParaMonteExample${EXAM_LANG} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/"
            cp ${ParaMonteExample_SRC_DIR}/main.${LANG_FILE_EXT} ${ParaMonteExample_BLD_DIR_CURRENT}/
            cp ${ParaMonteExample_SRC_DIR}/main.${LANG_FILE_EXT} ${ParaMonteExample_BIN_DIR_CURRENT}/
        fi

        if [ "${IS_Python_LANG}" = "true" ]; then

            echo >&2 "-- ParaMonteExample${EXAM_LANG} - from: ${ParaMonteExample_SRC_DIR_CURRENT}/README.md"
            echo >&2 "-- ParaMonteExample${EXAM_LANG} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/"
            cp ${ParaMonteExample_SRC_DIR_CURRENT}/README.md ${ParaMonteExample_BLD_DIR_CURRENT}/

            PythonScriptFileName=main.py
            if [ "${MPI_ENABLED}" = "true" ]; then PythonScriptFileName=main_mpi.py; fi
            cp ${ParaMonteExample_SRC_DIR_CURRENT}/${PythonScriptFileName} ${ParaMonteExample_BLD_DIR_CURRENT}/
            cp ${ParaMonteExample_SRC_DIR_CURRENT}/${PythonScriptFileName} ${ParaMonteExample_BIN_DIR_CURRENT}/

        fi

        # ParaMonte library example header/module files

        if [ "${EXAM_LANG}" = "Fortran" ]; then
            echo >&2
            echo >&2 "-- ParaMonteExample${EXAM_LANG} - copying the ParaMonte library Fortran module file paramonte.mod..."
            echo >&2 "-- ParaMonteExample${EXAM_LANG} - from: ${ParaMonte_MOD_DIR}"
            echo >&2 "-- ParaMonteExample${EXAM_LANG} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}"
            cp ${ParaMonte_MOD_DIR}/paradram_mod.mod ${ParaMonteExample_BLD_DIR_CURRENT}/
            cp ${ParaMonte_MOD_DIR}/paradram_mod.mod ${ParaMonteExample_BIN_DIR_CURRENT}/
            echo >&2
            echo >&2 "-- ParaMonteExample${EXAM_LANG} - copying the ParaMonte library Fortran module file paramonte.f90..."
            echo >&2 "-- ParaMonteExample${EXAM_LANG} - from: ${ParaMonteInterfaceFortran_SRC_DIR}"
            echo >&2 "-- ParaMonteExample${EXAM_LANG} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/"
            cp ${ParaMonteInterfaceFortran_SRC_DIR}/* ${ParaMonteExample_BLD_DIR_CURRENT}/
            cp ${ParaMonteInterfaceFortran_SRC_DIR}/* ${ParaMonteExample_BIN_DIR_CURRENT}/
        fi
        if [ "${EXAM_LANG}" = "C" ]; then
            echo >&2
            echo >&2 "-- ParaMonteExample${EXAM_LANG} - copying the ParaMonte library C header file paramonte.h..."
            echo >&2 "-- ParaMonteExample${EXAM_LANG} - from: ${ParaMonteInterfaceC_SRC_DIR}"
            echo >&2 "-- ParaMonteExample${EXAM_LANG} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/"
            cp ${ParaMonteInterfaceC_SRC_DIR}/* ${ParaMonteExample_BLD_DIR_CURRENT}/
            cp ${ParaMonteInterfaceC_SRC_DIR}/* ${ParaMonteExample_BIN_DIR_CURRENT}/
        fi

        # ParaMonte library example build and run if requested

        echo >&2
        if [ "${ParaMonteExample_RUN_ENABLED}" = "true" ] && [ "${IS_FortranC_LANG}" = "true" ]; then
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

    # echo >&2 

done
