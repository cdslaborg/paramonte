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
####       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
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

BUILD_NAME="buildParaMonteExample"
verify() {
    if [ $1 -eq 0 ]; then
        echo >&2 "-- ${pmattn} ${BoldGreen}The ParaMonte library $2 appears to have succeeded.${ColorReset}"
    else
        echo >&2
        echo >&2 "    -- ${BUILD_NAME} - FATAL: ParaMonte $2 appears to have failed."
        echo >&2 "    -- ${BUILD_NAME} - FATAL: If the source of the error cannot be identified,"
        echo >&2 "    -- ${BUILD_NAME} - FATAL: consider a fresh installation of ParaMonte's required compilers by calling"
        echo >&2 "    -- ${BUILD_NAME} - FATAL: "
        echo >&2 "    -- ${BUILD_NAME} - FATAL:     ./install --fresh"
        echo >&2 "    -- ${BUILD_NAME} - FATAL: "
        echo >&2 "    -- ${BUILD_NAME} - FATAL: If the error happens during the installation of ParaMonte prerequisites"
        echo >&2 "    -- ${BUILD_NAME} - FATAL: it is possible that the current existing GCC compiler collection installed"
        echo >&2 "    -- ${BUILD_NAME} - FATAL: on your system cannot compile the downloaded version of GCC that is required"
        echo >&2 "    -- ${BUILD_NAME} - FATAL: for ParaMonte build. In such case, make sure you have a GCC compiler collection"
        echo >&2 "    -- ${BUILD_NAME} - FATAL: version 7.1 or newer installed on your system, with an updated PATH environmental"
        echo >&2 "    -- ${BUILD_NAME} - FATAL: variable, then reinstall ParaMonte."
        echo >&2 "    -- ${BUILD_NAME} - FATAL: "
        echo >&2 "    -- ${BUILD_NAME} - FATAL: If all ParaMonte installation attempts fail, please report this issue at"
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

if [ "${isMacOS}" = "true" ]; then
    sharedFileExt="dylib"
else
    sharedFileExt="so"
fi

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

ParaMonteInterface_SRC_DIR_CURRENT="${ParaMonteInterface_SRC_DIR}/${LANG_NAME}"

# select examples to build

EXAM_LIST="mvn"
export EXAM_LIST

echo >&2
echo >&2 "-- ParaMonteExample${LANG_NAME} - The ParaMonte library example interface language: ${INTERFACE_LANGUAGE}"
echo >&2

# make example build directory

ParaMonteExample_BLD_DIR="${ParaMonte_BLD_DIR}"/example
echo >&2 "-- ParaMonteExample${LANG_NAME} - ParaMonte examples root directory: ${ParaMonteExample_BLD_DIR}"
if [[ -d "${ParaMonteExample_BLD_DIR}" ]]; then
    echo >&2 "-- ParaMonteExample${LANG_NAME} - ParaMonte  examples root directory already exists. skipping..."
else
    echo >&2 "-- ParaMonteExample${LANG_NAME} - generating ParaMonte examples root directory..."
    mkdir "${ParaMonteExample_BLD_DIR}/"
    verify $? "directory creation"
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
        verify $? "deletion of the old files"
    fi
    echo >&2 "-- ParaMonteExample${LANG_NAME} - generating ParaMonte ${EXAM_NAME} example build directory..."
    mkdir -p "${ParaMonteExample_BLD_DIR_CURRENT}/"
    verify $? "recursive directory creation"

    # The ParaMonte library kernel files

    ParaMonteExample_LIB_DIR_CURRENT="${ParaMonteExample_BLD_DIR_CURRENT}"
    if [ "${LANG_IS_Python}" = "true" ]; then ParaMonteExample_LIB_DIR_CURRENT="${ParaMonteExample_LIB_DIR_CURRENT}/paramonte"; fi
    if [ "${LANG_IS_MATLAB}" = "true" ]; then ParaMonteExample_LIB_DIR_CURRENT="${ParaMonteExample_LIB_DIR_CURRENT}/paramonte/lib"; fi

    if [[ -d "${ParaMonteExample_LIB_DIR_CURRENT}" ]]; then
        echo >&2 "-- ParaMonteExample${LANG_NAME} - ParaMonte ${EXAM_NAME} example library directory already exists. skipping..."
    else
        echo >&2 "-- ParaMonteExample${LANG_NAME} - generating ParaMonte ${EXAM_NAME} example library directory..."
        mkdir -p "${ParaMonteExample_LIB_DIR_CURRENT}/"
        verify $? "recursive directory creation"
    fi

    echo >&2
    echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library files..."
    echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${PMLIB_FULL_PATH}"
    echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_LIB_DIR_CURRENT}/"
    cp -R "${ParaMonte_LIB_DIR}/"libparamonte_* "${ParaMonteExample_LIB_DIR_CURRENT}/" || printCopyFailMsg

    # Copy the ParaMonte library dll dependency files

    if [ "${deploy_enabled}" = "true" ]; then

        #### gnu

        if [ "${PMCS}" = "gnu" ] && ! [ "${CAF_ENABLED}" = "true" ]; then

            #### first create the dependencies list

            unset inspector
            if [ "${isMacOS}" = "true" ] && command -v otool >/dev/null 2>&1 && command -v install_name_tool >/dev/null 2>&1; then
                inspector="otool -L"
            elif ! [ "${isMacOS}" = "true" ] && command -v ldd >/dev/null 2>&1; then # other unix-like
                inspector="ldd"
            fi

            if [ -z ${inspector+x} ]; then

                echo >&2
                echo >&2 "-- ParaMonteExample${LANG_NAME} - WARNING: The ldd (Linux) or otool/install_name_tool (macOS) could not be found."
                echo >&2 "-- ParaMonteExample${LANG_NAME} - WARNING: Skipping the shared library file copying..."
                echo >&2

            else

                #### loop over all existing shared files, and their dependencies recursively

                ishared=0
                sharedFilePathList=($(ls "${ParaMonteExample_LIB_DIR_CURRENT}"/libparamonte_*.${sharedFileExt}))
                sharedFilePathListLen=${#sharedFilePathList[@]}

                while [ "$ishared" -lt "${sharedFilePathListLen}" ]; do

                    sharedFilePath="${sharedFilePathList[$ishared]}"

                    echo >&2
                    echo >&2 "-- ParaMonteExample${LANG_NAME} - checking the dependencies of ${sharedFilePath}"

                    # @todo: point of weakness: The following for-loop assumes no white space in the dependency paths.

                    dependencyList=()
                    for dependencyFilePath in $(${inspector} "${sharedFilePath}"); do

                        dependencyFileName="$(basename "${dependencyFilePath}")"

                        # copy dependencyFilePath only if they are GNU shared files and are not MPI-related.

                        if  [ -f "${dependencyFilePath}" ] && ! [ -f "${ParaMonteExample_LIB_DIR_CURRENT}/${dependencyFileName}" ] && \
                          ( [[ "${dependencyFilePath}" =~ .*"gnu".* ]] || \
                            [[ "${dependencyFilePath}" =~ .*"gcc".* ]] || \
                            [[ "${dependencyFilePath}" =~ .*"gfortran".* ]] ) && \
                        ! ( [[ "${dependencyFilePath}" =~ .*"mpich".* ]] || \
                            [[ "${dependencyFilePath}" =~ .*"open-mpi".* ]] || 
                            [[ "${dependencyFilePath}" =~ .*"openmpi".* ]] 
                            ); then
                            #if ! [ "${isMacOS}" = "true" ] || ( [ "${isMacOS}" = "true" ] && ! [[ "${dependencyFilePath}" =~ .*"/usr/lib/".* ]] ); then
                                echo >&2 "-- ParaMonteExample${LANG_NAME} - dependency detected: ${dependencyFilePath}"
                                dependencyList+=("${dependencyFilePath}")
                            #fi
                        fi

                    done

                    dependencyListLen=${#dependencyList[@]}
                    dependencyListLenMinusOne="$(($dependencyListLen-1))"

                    if [ ${dependencyListLen} -eq 0 ]; then
                        echo >&2
                        echo >&2 "-- ParaMonteExample${LANG_NAME} - NOTE: No shared file dependencies were detected in the shared library file: ${sharedFilePath}"
                        echo >&2 "-- ParaMonteExample${LANG_NAME} - NOTE: Skipping the shared library file copying..."
                        echo >&2
                    else
                        echo >&2
                        echo >&2 "-- ParaMonteExample${LANG_NAME} - ${dependencyListLen} shared library file dependencies were detected."
                        echo >&2
                        for idep in $(seq 0 $dependencyListLenMinusOne); do
                            dependencyName=$(basename "${dependencyList[idep]}")
                            dependencyPathDestin="${ParaMonteExample_LIB_DIR_CURRENT}/${dependencyName}"
                            echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library dependency shared file..."
                            echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${dependencyList[idep]}"
                            echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${dependencyPathDestin}"
                            (yes | \cp -rf "${dependencyList[idep]}" "${dependencyPathDestin}") >/dev/null 2>&1 && {
                                echo >&2 "-- ParaMonteExample${LANG_NAME} - appending the shared file list with: ${dependencyPathDestin}"
                                sharedFilePathList+=("${dependencyPathDestin}")
                                if [ "${isMacOS}" = "true" ]; then
                                    echo >&2 "-- ParaMonteExample${LANG_NAME} - changing the install_name to @rpath for the dependency file..."
                                    install_name_tool -change \
                                    "${dependencyList[idep]}" \
                                    "@rpath/${dependencyName}" \
                                    "${sharedFilePath}" || {
                                        echo >&2
                                        echo >&2 "-- ParaMonteExample${LANG_NAME} - FATAL: Changing the install_name of the dependency file to @rpath failed."
                                        echo >&2
                                        exit 1
                                        #if [ "$BASH_SOURCE" == "$0" ]; then exit 30; else return 88; fi # return with an error message
                                    }
                                fi
                            } || {
                                echo >&2
                                echo >&2 "-- ParaMonteExample${LANG_NAME} - FATAL: The dependency file copy attempt failed at: ${dependencyList[idep]}"
                                echo >&2
                                exit 1
                                #if [ "$BASH_SOURCE" == "$0" ]; then exit 30; else return 88; fi # return with an error message
                            }
                            echo >&2
                        done
                    fi
                    sharedFilePathListLen=${#sharedFilePathList[@]}
                    ishared=$((ishared+1))

                    echo >&2
                    echo >&2 "-- ParaMonteExample${LANG_NAME} - The updated shared file list:"
                    printf   '%s\n' "${sharedFilePathList[@]}"
                    echo >&2 "-- ParaMonteExample${LANG_NAME} - The length of the shared file list: ${sharedFilePathListLen}"
                    echo >&2 "-- ParaMonteExample${LANG_NAME} - The current index through the list: ${ishared}"
                    echo >&2

                done # with while loop

            fi

        fi

        #### no need to execute the following anymore. The above block takes care of the shared files in a more robust way.

        if [ "false" = "true" ] && [ "${PMCS}" = "gnu" ] && ! [ "${CAF_ENABLED}" = "true" ] && ! [ "${isMacOS}" = "true" ]; then # caf does not have lib dependency

            ########################################################################################################################

            # On macOS, simply copying the dependencies does not work.
            # @todo A better strategy is to ask the user to install the relevant GFortran version on their system on their behalf via the build script.

            if ! [ -z ${Fortran_COMPILER_PATH+x} ]; then

                if [ "${isMacOS}" = "true" ]; then
                    Fortran_COMPILER_LIB_DIR_LIST="/usr/local/gfortran/lib"
                    sharedFileExt="dylib"
                else
                    sharedFileExt="so"
                    Fortran_COMPILER_DIR=$(dirname "${Fortran_COMPILER_PATH}")
                    FortranCompilerVersion=$("${Fortran_COMPILER_PATH}" -dumpversion)
                    FortranCompilerMajorVersion="$(cut -d '.' -f 1 <<< "$FortranCompilerVersion")"
                    Fortran_COMPILER_ROOT_DIR="${Fortran_COMPILER_DIR}"/..
                    #copySucceeded=false
                    Fortran_COMPILER_LIB_DIR_LIST="${Fortran_COMPILER_ROOT_DIR}/lib64:/usr/lib/gcc/x86_64-linux-gnu/${FortranCompilerMajorVersion}"
                fi
                FILE_LIST="libgfortran:libquadmath"

                for FILE in ${FILE_LIST//:/ }
                do
                    for Fortran_COMPILER_LIB_DIR in ${Fortran_COMPILER_LIB_DIR_LIST//:/ }
                    do
                        if [ -d "${Fortran_COMPILER_LIB_DIR}" ]; then
                            flist=$(( IFS=:; unset lsout; lsout=$(ls -dm "${Fortran_COMPILER_LIB_DIR}/${FILE}"*.${sharedFileExt}*); if ! [[ -z "${lsout// }" ]]; then echo "${lsout}, "; fi) 2>/dev/null)
                            for fpath in $(echo $flist | sed "s/,/ /g"); do
                                echo >&2
                                echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library dll dependency file..."
                                echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${fpath}"
                                echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_LIB_DIR_CURRENT}/"
                                (yes | \cp -rf "${fpath}" "${ParaMonteExample_LIB_DIR_CURRENT}/") >/dev/null 2>&1 || {
                                    echo >&2
                                    echo >&2 "-- ParaMonteExample${LANG_NAME} - WARNING: A ParaMonte library dll dependency file copy attempt failed at: ${fpath}"
                                    echo >&2
                                    continue
                                }
                            done
                        fi
                    done
                done

            else

                echo >&2
                echo >&2 "-- ParaMonteExample${LANG_NAME} - WARNING: The ParaMonte shared library dependency files could not be found."
                echo >&2

            fi

            #### copy MPI shared library files. UPDATE: Not a good idea to copy MPI installation files. It only creates a mess in the binary folders and rarely works.

            if [ "false" = "true" ] && ! [ -z ${MPIEXEC_PATH+x} ] && ! [ "${MPILIB_NAME}" = "openmpi" ] && [ "${MPI_ENABLED}" = "true" ]; then

                # copy files only if it is not openMPI. OpenMPI files interfer for local installations of the library on the user system. Example:
                # mpiexec: Error: unknown option "-np"
                # mpiexec: Error: unknown option "--oversubscribe"

                #echo >&2
                #echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the the mpiexec executable file..."
                #echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${MPIEXEC_PATH}*"
                #echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_LIB_DIR_CURRENT}/"
                #yes | \cp -rf "${MPIEXEC_PATH}"* "${ParaMonteExample_LIB_DIR_CURRENT}/"
                MPI_BIN_DIR=$(dirname "${MPIEXEC_PATH}")
                MPI_ROOT_DIR="${MPI_BIN_DIR}"/..
                MPI_LIB_DIR_LIST="${MPI_ROOT_DIR}/lib:${MPI_ROOT_DIR}/lib64"
                # copy all so files
                for MPI_LIB_DIR in ${MPI_LIB_DIR_LIST//:/ }
                do
                    echo >&2
                    echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library dll dependency file..."
                    echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${MPI_LIB_DIR}/*.so*"
                    echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_LIB_DIR_CURRENT}/"
                    if [ -d "${MPI_LIB_DIR}" ]; then
                        cp -rf "${MPI_LIB_DIR}/"*.so* "${ParaMonteExample_LIB_DIR_CURRENT}/" && break || {
                            echo >&2 "-- ParaMonteExample${LANG_NAME} - copy failed. skipping..."
                        }
                    fi
                done
                #keyList="libmpi:libmpifort"
                #for key in ${keyList//:/ }
                #do
                #    for MPI_LIB_DIR in ${MPI_LIB_DIR_LIST//:/ }
                #    do
                #        if [ -d "${MPI_LIB_DIR}" ]; then
                #            find "${MPI_LIB_DIR}" -name "${key}.so"* -print0 |
                #            while IFS= read -r -d '' libMpiPath; do
                #                echo >&2
                #                echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library dll dependency file..."
                #                echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${libMpiPath}"
                #                echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_LIB_DIR_CURRENT}/"
                #                #yes | \cp -rf "${libMpiPath}" "${ParaMonteExample_LIB_DIR_CURRENT}/" && {
                #                yes | \cp -rf "${libMpiPath}" "${ParaMonteExample_LIB_DIR_CURRENT}/"
                #            done
                #        fi
                #    done
                #done

                ##### copy numa shared files
                #
                #FILE_LIST="libnuma.so:libpciaccess.so"
                #FILE_DIR_LIST="/usr/lib64"
                #for FILE in ${FILE_LIST//:/ }
                #do
                #    for FILE_DIR in ${FILE_DIR_LIST//:/ }
                #    do
                #        if [ -d "${FILE_DIR}" ]; then
                #            flist=$(( IFS=:; unset lsout; lsout=$(ls -dm "${FILE_DIR}/${FILE}"*); if ! [[ -z "${lsout// }" ]]; then echo "${lsout}, "; fi) 2>/dev/null)
                #            for fpath in $(echo $flist | sed "s/,/ /g"); do
                #                echo >&2
                #                echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library dll dependency file..."
                #                echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${fpath}"
                #                echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_LIB_DIR_CURRENT}/"
                #                (yes | \cp -rf "${fpath}" "${ParaMonteExample_LIB_DIR_CURRENT}/") >/dev/null 2>&1 || {
                #                    echo >&2
                #                    echo >&2 "-- ParaMonteExample${LANG_NAME} - WARNING: A ParaMonte library dll dependency file copy attempt failed at: ${fpath}"
                #                    echo >&2
                #                    continue
                #                }
                #            done
                #        fi
                #    done
                #done

            fi

        fi

        #### intel - UPDATE: Not a good idea to copy MPI installation files. It only creates a mess in the binary folders and rarely works.

        if [ "false" = "true" ] && [ "${PMCS}" = "intel" ] && ! [ -z ${MPIEXEC_PATH+x} ] && [ "${MPI_ENABLED}" = "true" ]; then
            #echo >&2
            #echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the the mpiexec executable file..."
            #echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${MPIEXEC_PATH}*"
            #echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_LIB_DIR_CURRENT}/"
            #yes | \cp -rf "${MPIEXEC_PATH}"* "${ParaMonteExample_LIB_DIR_CURRENT}/"
            MPI_BIN_DIR=$(dirname "${MPIEXEC_PATH}")
            MPI_ROOT_DIR="${MPI_BIN_DIR}"/..
            MPI_LIB_DIR_LIST="${MPI_ROOT_DIR}/lib:${MPI_ROOT_DIR}/lib64"
            # copy all so files
            for MPI_LIB_DIR in ${MPI_LIB_DIR_LIST//:/ }
            do
                echo >&2
                echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library dll dependency file..."
                echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${MPI_LIB_DIR}/libmpifort.so*"
                echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${MPI_LIB_DIR}/libmpi.so*"
                echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_LIB_DIR_CURRENT}/"
                thisIntelBuild="${BTYPE}"; if [ "${BTYPE}" = "testing" ]; then thisIntelBuild="release"; fi
                if [ -d "${MPI_LIB_DIR}" ]; then
                    cp -rf "${MPI_LIB_DIR}/"libmpifort.so* "${ParaMonteExample_LIB_DIR_CURRENT}/" && \
                    cp -rf "${MPI_LIB_DIR}/${thisIntelBuild}_mt/"libmpi.so* "${ParaMonteExample_LIB_DIR_CURRENT}/" && \
                    break || {
                        echo >&2 "-- ParaMonteExample${LANG_NAME} - copy failed. skipping..."
                    }
                fi
            done
        fi

    fi

    # The ParaMonte library example required files

    echo >&2
    echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library ${EXAM_NAME} example build files in ${LANG_NAME} language..."

    # The ParaMonte library license file

    echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library license file..."
    echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonte_ROOT_DIR}/LICENSE.md"
    echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/LICENSE.md"
    cp "${ParaMonte_ROOT_DIR}/LICENSE.md" "${ParaMonteExample_BLD_DIR_CURRENT}/LICENSE.md" || printCopyFailMsg

    # The ParaMonte library README.md file

    echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library README.md file..."
    echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterface_SRC_DIR_CURRENT}/README.md"
    echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/README.md"
    cp "${ParaMonteInterface_SRC_DIR_CURRENT}/README.md" "${ParaMonteExample_BLD_DIR_CURRENT}/README.md" || printCopyFailMsg

    if [ "${LANG_IS_COMPILED}" = "true" ]; then

        # The ParaMonte library CHANGES.md file

        echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library CHANGES.md file..."
        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonte_ROOT_DIR}/CHANGES.md"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/CHANGES.md"
        cp "${ParaMonte_ROOT_DIR}/CHANGES.md" "${ParaMonteExample_BLD_DIR_CURRENT}/CHANGES.md" || printCopyFailMsg

        # The ParaMonte library build script

        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteExample_SRC_DIR}/build.sh"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}"
        cp "${ParaMonteExample_SRC_DIR}/build.sh" "${ParaMonteExample_BLD_DIR_CURRENT}/" || printCopyFailMsg
        chmod +x "${ParaMonteExample_BLD_DIR_CURRENT}"/build.sh
        if [[ -f "${SETUP_FILE_PATH}" ]]; then
            cp "${SETUP_FILE_PATH}" "${ParaMonteExample_BLD_DIR_CURRENT}/" || printCopyFailMsg
        fi

        # The ParaMonte library example header/module files

        echo >&2
        echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library C header file paramonte.h..."
        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterface_SRC_DIR_CURRENT}"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/"
        cp "${ParaMonteInterface_SRC_DIR_CURRENT}"/* "${ParaMonteExample_BLD_DIR_CURRENT}/" || printCopyFailMsg

        if [ "${LANG_NAME}" = "Fortran" ]; then
            echo >&2
            echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library Fortran module file paramonte.mod..."
            echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonte_MOD_DIR}"
            echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}"
            cp "${ParaMonte_MOD_DIR}/paradram_mod.mod" "${ParaMonteExample_BLD_DIR_CURRENT}/" || printCopyFailMsg
        fi

    fi

    if [ "${LANG_IS_DYNAMIC}" = "true" ]; then

        # The ParaMonte library CHANGES.md file

        echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library CHANGES.md file..."
        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterface_SRC_DIR_CURRENT}/CHANGES.md"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/CHANGES.md"
        cp "${ParaMonteInterface_SRC_DIR_CURRENT}/CHANGES.md" "${ParaMonteExample_BLD_DIR_CURRENT}/CHANGES.md" || printCopyFailMsg

        # The ParaMonte library interface files

        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterface_SRC_DIR_CURRENT}/paramonte"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/"
        cp -R "${ParaMonteInterface_SRC_DIR_CURRENT}/paramonte" "${ParaMonteExample_BLD_DIR_CURRENT}/" || printCopyFailMsg

        # The ParaMonte library auxiliary file (this must be done to generate the auxil folder)

        echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library auxiliary files"
        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterface_SRC_DIR}/auxil"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/"
        cp -R "${ParaMonteInterface_SRC_DIR}/auxil" "${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/" || printCopyFailMsg
        echo >&2

        # The ParaMonte kernel version file (this must appear only after the above)

        echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library kernel version file..."
        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonte_ROOT_DIR}/.VERSION"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/auxil/.VERSION_KERNEL"
        cp "${ParaMonte_ROOT_DIR}/.VERSION" "${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/auxil/.VERSION_KERNEL" || printCopyFailMsg

        # The ParaMonte interface version file (this must appear only after the above)

        echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library interface version file..."
        echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterface_SRC_DIR}/${LANG_NAME}/.VERSION"
        echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/auxil/.VERSION_INTERFACE"
        cp "${ParaMonteInterface_SRC_DIR}/${LANG_NAME}/.VERSION" "${ParaMonteExample_BLD_DIR_CURRENT}/paramonte/auxil/.VERSION_INTERFACE" || printCopyFailMsg

        if [ "${LANG_IS_Python}" = "true" ]; then

            # PyPI build - The ParaMonte library Python setup files

            echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library Python setup files..."
            echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteInterface_SRC_DIR_CURRENT}/setup/*"
            echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/"
            cp "${ParaMonteInterface_SRC_DIR_CURRENT}/setup/"* "${ParaMonteExample_BLD_DIR_CURRENT}/" || printCopyFailMsg

        fi

    fi

    # The ParaMonte library example input files

    ParaMonteExample_INP_DIR_CURRENT="${ParaMonteExample_SRC_DIR}/${EXAM_NAME}/input"
    echo >&2
    echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library ${EXAM_NAME} example input files in ${LANG_NAME} language..."
    echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteExample_INP_DIR_CURRENT}"
    echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}"
    cp -R "${ParaMonteExample_INP_DIR_CURRENT}"/* "${ParaMonteExample_BLD_DIR_CURRENT}"/ || printCopyFailMsg

    # The ParaMonte library example source files

    ParaMonteExample_SRC_DIR_CURRENT="${ParaMonteExample_SRC_DIR}/${EXAM_NAME}/${LANG_NAME}"

    echo >&2
    echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library ${EXAM_NAME} example source files in ${LANG_NAME} language..."

    # copy the example files

    mainFileName=main.${LANG_FILE_EXT}
    if [ "${LANG_IS_DYNAMIC}" = "true" ] && [ "${MPI_ENABLED}" = "true" ]; then mainFileName=main_mpi.${LANG_FILE_EXT}; fi

    echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteExample_SRC_DIR}/main.${LANG_FILE_EXT}"
    echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}/"
    cp "${ParaMonteExample_SRC_DIR}/${mainFileName}" "${ParaMonteExample_BLD_DIR_CURRENT}/" || printCopyFailMsg

    echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteExample_SRC_DIR_CURRENT}/*"
    echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BLD_DIR_CURRENT}"
    cp "${ParaMonteExample_SRC_DIR_CURRENT}"/* "${ParaMonteExample_BLD_DIR_CURRENT}/" || printCopyFailMsg

done

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: copy the first example to the bin directory
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

ParaMonteExample_BLD_DIR_CURRENT="${ParaMonteExample_BLD_DIR}/mvn"

if ! [ -d "${ParaMonte_BIN_DIR}" ]; then
    mkdir "${ParaMonte_BIN_DIR}/";
    verify $? "directory creation"
fi

ParaMonteExample_BIN_DIR_CURRENT="${ParaMonte_BIN_DIR}/${PMLIB_BASE_NAME}"
if [ "${LANG_IS_DYNAMIC}" = "true" ]; then ParaMonteExample_BIN_DIR_CURRENT="${ParaMonte_BIN_DIR}/libparamonte_${LANG_ABBR}"; fi
if ! [ -d "${ParaMonteExample_BIN_DIR_CURRENT}" ]; then
    mkdir "${ParaMonteExample_BIN_DIR_CURRENT}/";
    verify $? "recursive directory creation"
fi

echo >&2 "-- ParaMonteExample${LANG_NAME} - The ParaMonte ${LANG_NAME} library install directory: ${ParaMonteExample_BIN_DIR_CURRENT}"

if [ -d "${ParaMonteExample_BIN_DIR_CURRENT}" ]; then
    echo >&2 "-- ParaMonteExample${LANG_NAME} - The ParaMonte install directory already exists. skipping..."
else
    echo >&2 "-- ParaMonteExample${LANG_NAME} - generating ParaMonte install directory..."
    mkdir -p "${ParaMonteExample_BIN_DIR_CURRENT}"
    verify $? "recursive directory creation"
fi

echo >&2 "-- ParaMonteExample${LANG_NAME} - copying the ParaMonte library files to the install directory..."
echo >&2 "-- ParaMonteExample${LANG_NAME} - from: ${ParaMonteExample_BLD_DIR_CURRENT}"
echo >&2 "-- ParaMonteExample${LANG_NAME} -   to: ${ParaMonteExample_BIN_DIR_CURRENT}"
cp -rf "${ParaMonteExample_BLD_DIR_CURRENT}"/* "${ParaMonteExample_BIN_DIR_CURRENT}" || printCopyFailMsg

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# The ParaMonte library example build and run if requested
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

for EXAM_NAME in $EXAM_LIST
do

    echo >&2
    if [ "${ParaMonteExample_RUN_ENABLED}" = "true" ] && [ "${LANG_IS_COMPILED}" = "true" ]; then
        ParaMonteExample_BLD_DIR_CURRENT="${ParaMonteExample_BLD_DIR}/${EXAM_NAME}"
        (cd "${ParaMonteExample_BLD_DIR_CURRENT}" && ./build.sh)
        if [ $? -eq 0 ]; then
            echo >&2 "-- ParaMonteExample - ParaMonte example build successful."
            echo >&2
        else
            echo >&2 "-- ParaMonteExample - Failed to build the example. exiting..."
            echo >&2
            exit 1
        fi
        (cd "${ParaMonteExample_BLD_DIR_CURRENT}" && ./run.sh)
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
