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

# Usage: make dclean
# or
# Usage: make distclean
# This CMake script will delete build directories and files to bring the
# package back to it's distribution state

# We want to start from the top of the source dir, so if we are in build
# we want to start one directory up

GET_FILENAME_COMPONENT(BASEDIR ${CMAKE_SOURCE_DIR} NAME)
IF(${BASEDIR} STREQUAL "build")
    SET(TOPDIR "${CMAKE_SOURCE_DIR}/..")
ELSE()
    SET(TOPDIR "${CMAKE_SOURCE_DIR}")
ENDIF()

MACRO(GET_PARENT_DIRECTORIES search_string return_list grandparents)
    FILE(GLOB_RECURSE new_list ${search_string})
    SET(dir_list "")
    FOREACH(file_path ${new_list})
        GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
        # Remove an extra directory component to return grandparent
        IF(${grandparents})
            # Tack on a fake extension to trick CMake into removing a second
            # path component
            SET(dir_path "${dir_path}.tmp")
            GET_FILENAME_COMPONENT(dir_path ${dir_path} PATH)
        ENDIF(${grandparents})
        SET(dir_list ${dir_list} ${dir_path})
    ENDFOREACH()
    LIST(REMOVE_DUPLICATES dir_list)
    SET(${return_list} ${dir_list})
ENDMACRO()

# Find directories and files that we will want to remove
FILE(GLOB_RECURSE CMAKECACHE "${TOPDIR}/*CMakeCache.txt")
FILE(GLOB_RECURSE CMAKECACHE "${CMAKE_BINARY_DIR}/*CMakeCache.txt")
FILE(GLOB_RECURSE CMAKEINSTALL "${TOPDIR}/*cmake_install.cmake" "${TOPDIR}/*install_manifest.txt")
FILE(GLOB_RECURSE CMAKEINSTALL "${CMAKE_BINARY_DIR}/*cmake_install.cmake" "${CMAKE_BINARY_DIR}/*install_manifest.txt")
FILE(GLOB_RECURSE MAKEFILE "${TOPDIR}/*Makefile")
FILE(GLOB_RECURSE MAKEFILE "${CMAKE_BINARY_DIR}/*Makefile")
FILE(GLOB_RECURSE CMAKETESTFILES "${TOPDIR}/*CTestTestfile.cmake")
FILE(GLOB_RECURSE CMAKETESTFILES "${CMAKE_BINARY_DIR}/*CTestTestfile.cmake")
SET(TOPDIRECTORIES "${CMAKE_BINARY_DIR}/lib" 
                   "${CMAKE_BINARY_DIR}/obj"
                   "${CMAKE_BINARY_DIR}/mod"
                   "${CMAKE_BINARY_DIR}/bin"
                   "${CMAKE_BINARY_DIR}/test"
                   "${CMAKE_BINARY_DIR}/example"
)

# CMake has trouble finding directories recursively, so locate these
# files and then save the parent directory of the files
GET_PARENT_DIRECTORIES(Makefile.cmake CMAKEFILES 0)
GET_PARENT_DIRECTORIES(LastTest.log CMAKETESTING 1)

# Place these files and directories into a list
SET(DEL 
    ${MAKEFILE}
    ${CMAKECACHE}
    ${CMAKEFILES}
    ${CMAKEINSTALL}
    ${CMAKETESTING}
    ${CMAKETESTFILES}
    ${TOPDIRECTORIES}
)

# If we are not in the build dir, delete that as well
#IF(NOT (${BASEDIR} STREQUAL "build"))
#    FILE(GLOB BUILD "${TOPDIR}/build")
#    SET(DEL ${DEL} ${BUILD})
#ENDIF()

# Loop over the directories and delete each one
FOREACH(D ${DEL})
    message( STATUS " Removing${D}")
    IF(EXISTS ${D})
        FILE(REMOVE_RECURSE ${D})
    ENDIF()
ENDFOREACH()
