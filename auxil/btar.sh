#!/bin/bash
####################################################################################################################################
####################################################################################################################################
##
##   ParaMonte: plain powerful parallel Monte Carlo library.
##
##   Copyright (C) 2012-present, The Computational Data Science Lab
##
##   This file is part of the ParaMonte library.
##
##   ParaMonte is free software: you can redistribute it and/or modify it
##   under the terms of the GNU Lesser General Public License as published
##   by the Free Software Foundation, version 3 of the License.
##
##   ParaMonte is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
##   GNU Lesser General Public License for more details.
##
##   You should have received a copy of the GNU Lesser General Public License
##   along with the ParaMonte library. If not, see,
##
##       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
##
##   ACKNOWLEDGMENT
##
##   As per the ParaMonte library license agreement terms,
##   if you use any parts of this library for any purposes,
##   we ask you to acknowledge the use of the ParaMonte library
##   in your work (education/research/industry/development/...)
##   by citing the ParaMonte library as described on this page:
##
##       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
##
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
    echo >&2 "-- ParaMonte - FATAL: The requested input target directory ${TARGET_DIR} specified" 
    echo >&2 "-- ParaMonte - FATAL: with the input flag --dir does not exist."
    exit 1
fi
