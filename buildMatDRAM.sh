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
#
# NOTE: This is not a standalone file. This file must be called via install.sh script file. 
# NOTE: This is the Bash script file that builds the ParaMonte MatDRAM library.
# NOTE: Do not change the contents of this file unless you know what the consequences are.

BUILD_NAME="MatDRAM"; export BUILD_NAME
verify() {
    if [ $1 -eq 0 ]; then
        echo >&2 "-- ${BUILD_NAME} - ParaMonte::MatDRAM $2 appears to have succeeded."
    else
        echo >&2
        echo >&2 "    -- ${BUILD_NAME} - FATAL: ParaMonte::MatDRAM $2 appears to have failed."
        echo >&2 "    -- ${BUILD_NAME} - FATAL: If all ParaMonte::MatDRAM installation attempts fail, please report this issue at"
        echo >&2 "    -- ${BUILD_NAME} - FATAL: "
        echo >&2 "    -- ${BUILD_NAME} - FATAL:     https://github.com/shahmoradi/paramonte/issues"
        echo >&2 "    -- ${BUILD_NAME} - FATAL: "
        echo >&2 "    -- ${BUILD_NAME} - FATAL: or by contacting the ParaMonte authors directly (e.g., shahmoradi@utexas.edu)."
        echo >&2
        echo >&2 "    -- ${BUILD_NAME} - gracefully exiting."
        echo >&2
        exit 1
    fi
}

printCopyFailMsg() {
    echo >&2
    echo >&2 "-- ParaMonteExample${LANG_NAME} - Copy action failed. Please resolve the error. Gracefully exiting..."
    echo >&2
    exit 1
}

FILE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

ParaMonte_ROOT_DIR="${FILE_DIR}"
export ParaMonte_ROOT_DIR

ParaMonteExample_SRC_DIR="${ParaMonte_ROOT_DIR}/example"
ParaMonteInterface_SRC_DIR="${ParaMonte_ROOT_DIR}/src/interface"

echo >&2 "-- ${BUILD_NAME} - project root directory: ${ParaMonte_ROOT_DIR}"

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: build MatDRAM library example objects and executable
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

echo >&2
echo >&2 "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
echo >&2 "::::                                                                                                                            ::::"
echo >&2 "                                                  MatDRAM MATLAB Library Build                                                      "
echo >&2 "::::                                                                                                                            ::::"
echo >&2 "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
echo >&2

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: make bin directory
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

ParaMonte_BIN_DIR="${ParaMonte_ROOT_DIR}/bin"
echo >&2 -- "${BUILD_NAME} - The MatDRAM binaries directory: ${ParaMonte_BIN_DIR}"
if ! [ -d "${ParaMonte_BIN_DIR}" ]; then
    mkdir "${ParaMonte_BIN_DIR}"
    verify $? "binary directory creation"
fi

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: setup examples' interface language
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

LANG_FILE_EXT="m"
LANG_NAME="MATLAB"

if [ -z ${LANG_NAME+x} ]; then
    echo >&2 
    echo >&2 "-- ${BUILD_NAME}Example - Fatal Error: unrecognized or no language specified. exiting..."
    echo >&2 
    exit 1
fi

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: build examples
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# select examples to build

EXAM_LIST="mvn"

# set and make example directories

MatDRAM_BLD_DIR="${ParaMonte_ROOT_DIR}/build/libparamonte_matdram"
MatDRAM_EXP_DIR="${MatDRAM_BLD_DIR}/example"

ParaMonteInterface_SRC_DIR_CURRENT="${ParaMonteInterface_SRC_DIR}/${LANG_NAME}"

echo >&2 
echo >&2 "-- ${BUILD_NAME} - generating the MatDRAM library examples in ${LANG_NAME} language..."
echo >&2 "-- ${BUILD_NAME} - The MatDRAM ${LANG_NAME} examples directory: ${MatDRAM_EXP_DIR}"

for EXAM_NAME in $EXAM_LIST
do

    MatDRAM_BLD_DIR_CURRENT="${MatDRAM_EXP_DIR}/${EXAM_NAME}"

    echo >&2 "-- ${BUILD_NAME} - The MatDRAM library ${EXAM_NAME} example directory: ${MatDRAM_BLD_DIR_CURRENT}"
    if [ -d "${MatDRAM_BLD_DIR_CURRENT}" ]; then
        echo >&2 "-- ${BUILD_NAME} - previous example build detected. deleting the old contents..."
        rm -rf "${MatDRAM_BLD_DIR_CURRENT}"
        verify $? "deletion of the old files"
        echo >&2 "-- ${BUILD_NAME} - regenerating the MatDRAM library ${EXAM_NAME} example directory: ${MatDRAM_BLD_DIR_CURRENT}"
    fi
    mkdir -p "${MatDRAM_BLD_DIR_CURRENT}/"
    verify $? "ParaMonte::MatDRAM example directory creation"

    # The MatDRAM library example required files

    echo >&2 "-- ${BUILD_NAME} - copying the MatDRAM library ${EXAM_NAME} example required files in ${LANG_NAME} language..."

    # The MatDRAM library README.md file

    # echo >&2 "-- ${BUILD_NAME} - copying the MatDRAM library README.md file..."
    # echo >&2 "-- ${BUILD_NAME} - from: ${ParaMonteInterface_SRC_DIR_CURRENT}/README.md"
    # echo >&2 "-- ${BUILD_NAME} -   to: ${MatDRAM_BLD_DIR_CURRENT}/paramonte/README.md"
    # cp "${ParaMonteInterface_SRC_DIR_CURRENT}/README.md" "${MatDRAM_BLD_DIR_CURRENT}/README.md" || printCopyFailMsg

    # The MatDRAM library CHANGES.md file

    # echo >&2 "-- ${BUILD_NAME} - copying the MatDRAM library CHANGES.md file..."
    # echo >&2 "-- ${BUILD_NAME} - from: ${ParaMonteInterface_SRC_DIR_CURRENT}/CHANGES.md"
    # echo >&2 "-- ${BUILD_NAME} -   to: ${MatDRAM_BLD_DIR_CURRENT}/paramonte/CHANGES.md"
    # copy "${ParaMonteInterface_SRC_DIR_CURRENT}/CHANGES.md" "${MatDRAM_BLD_DIR_CURRENT}/CHANGES.md" || printCopyFailMsg

    # The MatDRAM library license file

    echo >&2
    echo >&2 "-- ${BUILD_NAME} - copying the MatDRAM library license file..."
    echo >&2 "-- ${BUILD_NAME} - from: ${ParaMonte_ROOT_DIR}/LICENSE.md"
    echo >&2 "-- ${BUILD_NAME} -   to: ${MatDRAM_BLD_DIR_CURRENT}/LICENSE.md"
    cp "${ParaMonte_ROOT_DIR}/LICENSE.md" "${MatDRAM_BLD_DIR_CURRENT}/LICENSE.md" || printCopyFailMsg

    # The MatDRAM library README file

    echo >&2
    echo >&2 "-- ${BUILD_NAME} - copying the MatDRAM library license file..."
    echo >&2 "-- ${BUILD_NAME} - from: ${ParaMonteInterface_SRC_DIR_CURRENT}/README.MatDRAM.md"
    echo >&2 "-- ${BUILD_NAME} -   to: ${MatDRAM_BLD_DIR_CURRENT}/README.MatDRAM.md"
    cp "${ParaMonteInterface_SRC_DIR_CURRENT}/README.MatDRAM.md" "${MatDRAM_BLD_DIR_CURRENT}/README.md" || printCopyFailMsg

    # The MatDRAM library interface files

    echo >&2
    echo >&2 "-- ${BUILD_NAME} - from: ${ParaMonteInterface_SRC_DIR_CURRENT}/libparamonte_matdram"
    echo >&2 "-- ${BUILD_NAME} -   to: ${MatDRAM_BLD_DIR_CURRENT}/paramonte/"
    cp -R "${ParaMonteInterface_SRC_DIR_CURRENT}/paramonte" "${MatDRAM_BLD_DIR_CURRENT}/paramonte/" || printCopyFailMsg

    # The MatDRAM library ParaDRAM class method files, required for postprocess reading of output of data

    COPY_PATH_SOURCE="${ParaMonteInterface_SRC_DIR_CURRENT}/paramonte/interface/@ParaDRAM/readMarkovChain.m"
    COPY_PATH_DESTIN="${MatDRAM_BLD_DIR_CURRENT}/paramonte/kernel/@ParaDRAM_class/"
    echo >&2
    echo >&2 "-- ${BUILD_NAME} - from: ${COPY_PATH_SOURCE}"
    echo >&2 "-- ${BUILD_NAME} -   to: ${COPY_PATH_DESTIN}"
    cp "${COPY_PATH_SOURCE}" "${COPY_PATH_DESTIN}" || printCopyFailMsg

    COPY_PATH_SOURCE="${ParaMonteInterface_SRC_DIR_CURRENT}/paramonte/interface/@ParaMonteSampler/readChain.m"
    COPY_PATH_DESTIN="${MatDRAM_BLD_DIR_CURRENT}/paramonte/kernel/@ParaDRAM_class/"
    echo >&2
    echo >&2 "-- ${BUILD_NAME} - from: ${COPY_PATH_SOURCE}"
    echo >&2 "-- ${BUILD_NAME} -   to: ${COPY_PATH_DESTIN}"
    cp "${COPY_PATH_SOURCE}" "${COPY_PATH_DESTIN}" || printCopyFailMsg

    COPY_PATH_SOURCE="${ParaMonteInterface_SRC_DIR_CURRENT}/paramonte/interface/@ParaMonteSampler/readSample.m"
    COPY_PATH_DESTIN="${MatDRAM_BLD_DIR_CURRENT}/paramonte/kernel/@ParaDRAM_class/"
    echo >&2
    echo >&2 "-- ${BUILD_NAME} - from: ${COPY_PATH_SOURCE}"
    echo >&2 "-- ${BUILD_NAME} -   to: ${COPY_PATH_DESTIN}"
    cp "${COPY_PATH_SOURCE}" "${COPY_PATH_DESTIN}" || printCopyFailMsg

    COPY_PATH_SOURCE="${ParaMonteInterface_SRC_DIR_CURRENT}/paramonte/interface/@ParaMonteSampler/readRestart.m"
    COPY_PATH_DESTIN="${MatDRAM_BLD_DIR_CURRENT}/paramonte/kernel/@ParaDRAM_class/"
    echo >&2
    echo >&2 "-- ${BUILD_NAME} - from: ${COPY_PATH_SOURCE}"
    echo >&2 "-- ${BUILD_NAME} -   to: ${COPY_PATH_DESTIN}"
    cp "${COPY_PATH_SOURCE}" "${COPY_PATH_DESTIN}" || printCopyFailMsg

    COPY_PATH_SOURCE="${ParaMonteInterface_SRC_DIR_CURRENT}/paramonte/interface/@ParaMonteSampler/readReport.m"
    COPY_PATH_DESTIN="${MatDRAM_BLD_DIR_CURRENT}/paramonte/kernel/@ParaDRAM_class/"
    echo >&2
    echo >&2 "-- ${BUILD_NAME} - from: ${COPY_PATH_SOURCE}"
    echo >&2 "-- ${BUILD_NAME} -   to: ${COPY_PATH_DESTIN}"
    cp "${COPY_PATH_SOURCE}" "${COPY_PATH_DESTIN}" || printCopyFailMsg

    COPY_PATH_SOURCE="${ParaMonteInterface_SRC_DIR_CURRENT}/paramonte/interface/@ParaMonteSampler/readTabular.m"
    COPY_PATH_DESTIN="${MatDRAM_BLD_DIR_CURRENT}/paramonte/kernel/@ParaDRAM_class/"
    echo >&2
    echo >&2 "-- ${BUILD_NAME} - from: ${COPY_PATH_SOURCE}"
    echo >&2 "-- ${BUILD_NAME} -   to: ${COPY_PATH_DESTIN}"
    cp "${COPY_PATH_SOURCE}" "${COPY_PATH_DESTIN}" || printCopyFailMsg

    COPY_PATH_SOURCE="${ParaMonteInterface_SRC_DIR_CURRENT}/paramonte/interface/@ParaMonteSampler/getFilePathList.m"
    COPY_PATH_DESTIN="${MatDRAM_BLD_DIR_CURRENT}/paramonte/kernel/@ParaDRAM_class/"

    echo >&2
    echo >&2 "-- ${BUILD_NAME} - from: ${COPY_PATH_SOURCE}"
    echo >&2 "-- ${BUILD_NAME} -   to: ${COPY_PATH_DESTIN}"
    cp "${COPY_PATH_SOURCE}" "${COPY_PATH_DESTIN}" || printCopyFailMsg

    # The MatDRAM library banner file

    echo >&2
    echo >&2 "-- ${BUILD_NAME} - copying the MatDRAM library auxiliary files"
    echo >&2 "-- ${BUILD_NAME} - from: ${ParaMonteInterface_SRC_DIR}/auxil/"
    echo >&2 "-- ${BUILD_NAME} -   to: ${MatDRAM_BLD_DIR_CURRENT}/paramonte/auxil/"
    cp -R "${ParaMonteInterface_SRC_DIR}/auxil" "${MatDRAM_BLD_DIR_CURRENT}/paramonte/" || printCopyFailMsg

    # The MatDRAM library kernel version file (must appear only after the above)

    echo >&2
    echo >&2 "-- ${BUILD_NAME} - copying the MatDRAM library kernel version file..."
    echo >&2 "-- ${BUILD_NAME} - from: ${ParaMonte_ROOT_DIR}/.VERSION"
    echo >&2 "-- ${BUILD_NAME} -   to: ${MatDRAM_BLD_DIR_CURRENT}/paramonte/auxil/.VERSION_KERNEL"
    cp "${ParaMonte_ROOT_DIR}/.VERSION" "${MatDRAM_BLD_DIR_CURRENT}/paramonte/auxil/.VERSION_KERNEL" || printCopyFailMsg

    # The MatDRAM library interface version file (must appear only after the above)

    echo >&2
    echo >&2 "-- ${BUILD_NAME} - copying the MatDRAM library interface version file..."
    echo >&2 "-- ${BUILD_NAME} - from: ${ParaMonteInterface_SRC_DIR_CURRENT}/.VERSION"
    echo >&2 "-- ${BUILD_NAME} -   to: ${MatDRAM_BLD_DIR_CURRENT}/paramonte/auxil/.VERSION_INTERFACE"
    cp "${ParaMonteInterface_SRC_DIR_CURRENT}/.VERSION" "${MatDRAM_BLD_DIR_CURRENT}/paramonte/auxil/.VERSION_INTERFACE" || printCopyFailMsg

    # The MatDRAM library example input files

    ParaMonteExample_SRC_DIR_CURRENT="${ParaMonteExample_SRC_DIR}/${EXAM_NAME}/${LANG_NAME}"
    mainFileName=main.${LANG_FILE_EXT}
    echo >&2
    echo >&2 "-- ${BUILD_NAME} - copying the MatDRAM library ${EXAM_NAME} example source files in ${LANG_NAME} language..."
    cp "${ParaMonteExample_SRC_DIR}/${mainFileName}" "${MatDRAM_BLD_DIR_CURRENT}/" || printCopyFailMsg

    echo >&2
    echo >&2 "-- ${BUILD_NAME} - from: ${ParaMonteExample_SRC_DIR_CURRENT}/*"
    echo >&2 "-- ${BUILD_NAME} -   to: ${MatDRAM_BLD_DIR_CURRENT}"
    cp -R "${ParaMonteExample_SRC_DIR_CURRENT}"/* "${MatDRAM_BLD_DIR_CURRENT}/" || printCopyFailMsg

done

echo >&2 

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: copy the first example to the bin directory
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

MatDRAM_BLD_DIR_CURRENT="${MatDRAM_EXP_DIR}/mvn"

MatDRAMExample_BIN_DIR_CURRENT="${ParaMonte_BIN_DIR}/libparamonte_matdram"

echo >&2 "-- ${BUILD_NAME} - The MatDRAM ${LANG_NAME} library binary directory: ${MatDRAMExample_BIN_DIR_CURRENT}"

if ! [ -d ${MatDRAMExample_BIN_DIR_CURRENT} ]; then
    mkdir "${MatDRAMExample_BIN_DIR_CURRENT}/"
    verify $? "directory creation"
fi

echo >&2
echo >&2 "-- ${BUILD_NAME} - copying the MatDRAM library files to the bin folder..."
echo >&2 "-- ${BUILD_NAME} - from: ${MatDRAM_BLD_DIR_CURRENT}"
echo >&2 "-- ${BUILD_NAME} -   to: ${MatDRAMExample_BIN_DIR_CURRENT}"
cp -R "${MatDRAM_BLD_DIR_CURRENT}"/* "${MatDRAMExample_BIN_DIR_CURRENT}" || printCopyFailMsg

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: build test
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

MatDRAMTest_BLD_DIR="${MatDRAM_BLD_DIR}/test"
MatDRAMTest_SRC_DIR="${ParaMonteInterface_SRC_DIR_CURRENT}/test"

echo >&2
echo >&2 "-- ${BUILD_NAME} - copying the MatDRAM library files to the test folder..."
echo >&2 "-- ${BUILD_NAME} - from: ${MatDRAM_BLD_DIR_CURRENT}"
echo >&2 "-- ${BUILD_NAME} -   to: ${MatDRAMTest_BLD_DIR}"
cp -R "${MatDRAM_BLD_DIR_CURRENT}"/* "${MatDRAMTest_BLD_DIR}" || printCopyFailMsg

cd ${MatDRAMTest_BLD_DIR}
for file in main.m logfunc.m LICENSE.md; do
    if [ -f "${MatDRAMTest_BLD_DIR}/${file}" ]; then
        rm "${MatDRAMTest_BLD_DIR}/${file}" || verify $? "MatDRAM test ${file} file removal"
    fi
done

echo >&2
cp "${MatDRAMTest_SRC_DIR}/testParaMonte.m" "${MatDRAMTest_BLD_DIR}/testParaMonte.m" || printCopyFailMsg
cp "${MatDRAMTest_SRC_DIR}/getLogFunc.m" "${MatDRAMTest_BLD_DIR}/getLogFunc.m" || printCopyFailMsg

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: copy the first example to the bin directory
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if [ "${ParaMonteExample_RUN_ENABLED}" = "true" ]; then

    for EXAM_NAME in ${EXAM_LIST}; do

        echo >&2 "-- ${BUILD_NAME} - Running the MatDRAM library's ${EXAM_NAME} example."

        # The MatDRAM library example build and run if requested

        MatDRAM_BLD_DIR_CURRENT="${MatDRAM_EXP_DIR}/${EXAM_NAME}"
        cd "${MatDRAM_BLD_DIR_CURRENT}" && matlab -batch "main"
        if [ $? -eq 0 ]; then
            echo >&2 "-- ${BUILD_NAME} - The MatDRAM library example build/run appears to have succeeded."
            echo >&2
        else
            echo >&2 "-- ${BUILD_NAME} - Fatal Error: The MatDRAM library example build/run failed. exiting..."
            echo >&2
            exit 1
        fi

    done

fi

echo >&2 
echo >&2 "-- ${BUILD_NAME} -      the MatDRAM build files are located at: ${MatDRAM_BLD_DIR_CURRENT}"
echo >&2 "-- ${BUILD_NAME} -    the MatDRAM library files are located at: ${MatDRAMExample_BIN_DIR_CURRENT}"
echo >&2 
