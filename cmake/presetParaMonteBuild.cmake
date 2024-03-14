####################################################################################################################################
####################################################################################################################################
####                                                                                                                            ####
####    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           ####
####                                                                                                                            ####
####    Copyright (C) 2012-present, The Computational Data Science Lab                                                          ####
####                                                                                                                            ####
####    This file is part of the ParaMonte library.                                                                             ####
####                                                                                                                            ####
####    LICENSE                                                                                                                 ####
####                                                                                                                            ####
####       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          ####
####                                                                                                                            ####
####################################################################################################################################
####################################################################################################################################

#   This script sets CMAKE_BUILD_TYPE and CMAKE_CONFIGURATION_TYPES
#   Note that these variables must be set before project() or enable_language()
#   CMAKE commands are called in a project when a new build tree is first created.

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set the desired compiler suite
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


set(CMAKE_CONFIGURATION_TYPES "debug; testing; release; ipo; tuned; native" FORCE) # "RelWithDebInfo"

if (DEFINED build)

    string(TOLOWER ${build} build)

    message(NOTICE "${pmattn} User-specified `build` detected.")

    #   Warn user if conflicting values are specified (or cached) independently by two separate variables: `build` `CMAKE_BUILD_TYPE`

    if (DEFINED CMAKE_BUILD_TYPE AND NOT "${CMAKE_BUILD_TYPE}" STREQUAL "")
        string(TOLOWER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_LOWER)
        if (NOT "${build}" STREQUAL "${CMAKE_BUILD_TYPE_LOWER}")
            message(NOTICE
                    "${pmwarn} The CMake variable `CMAKE_BUILD_TYPE` and the CMake argument \n"
                    "${pmwarn} `build(=${build})` are both simultaneously defined but different.\n"
                    "${pmwarn} The variable `CMAKE_BUILD_TYPE(=${CMAKE_BUILD_TYPE})` will be overwritten\n"
                    "${pmwarn} with the value of `build`."
                    )
        endif()
    endif()

    #   Preset internal variable based on the user choice.

    if (# release
        #"${build}" MATCHES "[Ii][Pp][Oo]" OR
        #"${build}" MATCHES "[Tt][Uu][Nn][Ee][Dd]" OR
        #"${build}" MATCHES "[Rr][Ee][Ll][Ee][Aa][Ss][Ee]"
         "${build}" STREQUAL "release" OR "${build}" STREQUAL "ipo" OR "${build}" STREQUAL "tuned" OR "${build}" STREQUAL "native"
        )
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING "The ParaMonte library build type." FORCE)
    elseif (
           #"${build}" MATCHES "[Tt][Ee][Ss][Tt][Ii][Nn][Gg]" OR
           #"${build}" MATCHES "[Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo]"
            "${build}" STREQUAL "testing" OR "${build}" STREQUAL "relwithdebinfo"
            )
        set(build "testing")
        set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "The ParaMonte library build type." FORCE)
    elseif (
           #"${build}" MATCHES "[Dd][Ee][Bb][Uu][Gg]" OR
            "${build}" STREQUAL "debug"
            )
        set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "The ParaMonte library build type." FORCE)
    else()
        printUsage()
        message(FATAL_ERROR
                "${pmfatal} The user-specified build is unsupported. build=${build}\n"
                "${pmfatal} See the above usage notes for the possible choices.\n"
                )
    endif()

elseif (DEFINED CMAKE_BUILD_TYPE)

    message(NOTICE "${pmattn} Preset CMAKE_BUILD_TYPE detected. CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
    set(build "${CMAKE_BUILD_TYPE}")
    string(TOLOWER "${build}" build)

else()

    set(build "release")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "ParaMonte library build type." FORCE)
    message(NOTICE "${pmattn} Defaulting to `release` build...")

endif()

message(NOTICE "${pmattn} CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
message(NOTICE "${pmattn} build=${build}")

#set (build ${build} CACHE STRING "Select CMAKE configuration to build" FORCE)
#set_property(CACHE build PROPERTY STRINGS ${CMAKE_CONFIGURATION_TYPES})
#unset(build)
