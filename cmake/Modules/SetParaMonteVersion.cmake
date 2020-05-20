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

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Name project and specify source languages. Parse version from .VERSION file so that more info can be added and easier to get from scripts
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

file( STRINGS ".VERSION" first_line LIMIT_COUNT 1 )

string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+(-rc[0-9]+)?" ParaMonteVersion "${first_line}")

if( (NOT (ParaMonteVersion MATCHES "[0-9]+\\.[0-9]+\\.[0-9]+(-rc[0-9]+)?")) AND (EXISTS "${CMAKE_SOURCE_DIR}/.git"))
    find_package(Git QUIET)
    if(GIT_FOUND)
        execute_process ( COMMAND "${GIT_EXECUTABLE}" describe --abbrev=0
                          WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
                          RESULT_VARIABLE git_status
                          OUTPUT_VARIABLE git_output
                          OUTPUT_STRIP_TRAILING_WHITESPACE
                        )
        if((git_status STREQUAL "0") AND (git_output MATCHES "[0-9]+\\.[0-9]+\\.[0-9]+(-rc[0-9]+)?"))
            set(ParaMonteVersion "${git_output}")
        endif()
        execute_process ( COMMAND "${GIT_EXECUTABLE}" describe --always
                          WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
                          RESULT_VARIABLE git_status
                          OUTPUT_VARIABLE full_git_describe
                          OUTPUT_STRIP_TRAILING_WHITESPACE
                        )
        if(NOT (git_status STREQUAL "0"))
            set(full_git_describe NOTFOUND)
        endif()
    else()
        message ( WARNING 
                " \n"
                " ${pmwarn}\n"
                " ${pmattn} Could not find git executable!"
                )
    endif()
endif()

if(NOT (ParaMonteVersion MATCHES "[0-9]+\\.[0-9]+\\.[0-9]+(-rc[0-9]+)?"))
    message ( WARNING 
            " \n"
            " ${pmwarn}\n"
            " ${pmattn} Could not extract version from git, falling back on .VERSION, line 3."
            )
    file( STRINGS ".VERSION" ParaMonteVersion REGEX "[0-9]+\\.[0-9]+\\.[0-9]+(-rc[0-9]+)?" )
endif()

if(NOT full_git_describe)
    set(full_git_describe ${ParaMonteVersion})
endif()

string(REGEX REPLACE "-rc[0-9]+$" ".0" PARAMONTE_CMAKE_PROJECT_VERSION "${ParaMonteVersion}" )

#file(READ "${CMAKE_CURRENT_LIST_DIR}/cmake/ParaMonteBannerBatch.txt" ParaMonte_BANNER)
##file(READ "${CMAKE_CURRENT_LIST_DIR}/cmake/ParaMonteBannerCmake.txt" ParaMonte_BANNER)
## replace build and version place holders with the correct values
##string(CONFIGURE "${ParaMonte_BANNER}" ParaMonte_BANNER @ONLY)
##string(STRIP "${ParaMonte_BANNER}" ParaMonte_BANNER)
#message("\n${ParaMonte_BANNER}")

if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/.git")
   message( STATUS "${pmattn} Build from git repository detected" )
endif()

message( STATUS "${pmattn} Running with CMake from: ${CMAKE_COMMAND}" )
message( STATUS "${pmattn} ParaMonte top source dir: ${CMAKE_CURRENT_SOURCE_DIR}" )

