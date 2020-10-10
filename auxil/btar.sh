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
# Usage:
#
#   ./auxil/btar.sh --dir ./bin/
#
####################################################################################################################################
# parse arguments
####################################################################################################################################

AUXIL_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
TARGET_DIR="${AUXIL_DIR}/../bin/"

while [ "$1" != "" ]; do
    case $1 in
        -d | --dir )    shift
                        TARGET_DIR=$1
                        ;;
        -h | --help )   #usage
                        echo >&2 "-- ParaMonte - Sadly, there is no help here, because you should not be here if you do not know why you are here."
                        exit
                        ;;
        * )             echo >&2 "-- ParaMonte - FATAL: the input flag is not recognized: $1"
                        #usage
                        exit 1
    esac
    shift
done

if [ -d "${TARGET_DIR}" ]; then
    echo >&2
    echo >&2 "-- ParaMonte - compressing all subdirectories in the directory: ${TARGET_DIR}"
    echo >&2
    cd "${TARGET_DIR}"
    for subdir in ./*; do
        if [ -d "${subdir}" ]; then
            tarfile="${subdir}.tar.gz"
            #cd "${subdir}"
            if [ -f "${tarfile}" ]; then
                echo >&2 "-- ParaMonte - WARNING: compressed subdirectory already exists: ${tarfile}"
                echo >&2 "-- ParaMonte - WARNING: skipping..."
            else
                echo >&2 "-- ParaMonte - compressing subdirectory: ${subdir}"
                tar -zcvf ${tarfile} --exclude="${subdir}/setup.sh" ${subdir}
                #cd ..
            fi
        else
            echo >&2 "-- ParaMonte - WARNING: non-directory object detected: ${subdir}"
        fi
    done
else
    echo >&2
    echo >&2 "-- ParaMonte - FATAL: The requested input target directory ${TARGET_DIR} specified" 
    echo >&2 "-- ParaMonte - FATAL: with the input flag --dir does not exist."
    echo >&2
    exit 1
fi
