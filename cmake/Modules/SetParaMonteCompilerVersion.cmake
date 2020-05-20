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
# Set CMAKE_Fortran_COMPILER_VERSION if CMake doesn't do it for us
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#unset(CMAKE_Fortran_COMPILER_VERSION)
if ( NOT CMAKE_Fortran_COMPILER_VERSION )
    if ( NOT (CMAKE_VERSION VERSION_LESS 3.3.1) )
        message ( AUTHOR_WARNING 
                " \n"
                " ${pmwarn}\n"
                " ${pmattn} CMake ${CMAKE_VERSION} should know about Fortran compiler versions\n"
                " ${pmattn} but is missing CMAKE_Fortran_COMPILER_VERSION variable."
                )
    endif()
    # No CMAKE_Fortran_COMPILER_VERSION set, build our own
    # Try extracting it directly from ISO_Fortran_ENV's compiler_version
    # Write program for introspection
    file( WRITE "${CMAKE_BINARY_DIR}/getCompilerVersion.f90"
        "
            program main
                use iso_Fortran_env, only: compiler_version, output_unit
                write(output_unit,'(a)') compiler_version()
            end program
        "
        )
    try_run( PROG_RAN COMPILE_SUCCESS "${CMAKE_BINARY_DIR}" "${CMAKE_BINARY_DIR}/getCompilerVersion.f90" RUN_OUTPUT_VARIABLE VER_STRING )
    if ( COMPILE_SUCCESS )
        string( REGEX MATCH "[0-9]+\\.[0-9]+(\\.[0-9]+)?" DETECTED_VER "${VER_STRING}" )
        message( STATUS "${pmattn} Detected Fortran compiler as ${VER_STRING}" )
        message( STATUS "${pmattn} Extracted version number: ${DETECTED_VER}" )
    endif()
    if( ( NOT COMPILE_SUCCESS ) OR ( NOT DETECTED_VER ) )
        message ( WARNING
                " \n"
                " ${pmwarn}\n"
                " ${pmattn} Could not reliably detect Fortran compiler version.\n"
                " ${pmattn} It will be inferred it from the C compiler if it matches the Fortran compiler ID."
                )
    endif()
    if( "${CMAKE_C_COMPILER_ID}" MATCHES "${CMAKE_Fortran_COMPILER_ID}" )
        set( DETECTED_VER "${CMAKE_C_COMPILER_VERSION}" )
    else()
        message ( FATAL_ERROR
                " \n"
                "${pmfatal}\n"
                "   Exhausted all possible means of detecting the Fortran compiler version, all in vain.\n"
                "   Please install the latest GNU or Intel Fortran compilers."
                )
    endif()
    set( CMAKE_Fortran_COMPILER_VERSION "${DETECTED_VER}" )
endif()

# if(CMAKE_BUILD_TYPE MATCHES "Debug|DEBUG|debug")
#     add_definitions(-DDBG_ENABLED)
# endif()

# We have populated CMAKE_Fortran_COMPILER_VERSION if it was missing

if(intel_compiler AND (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 18.0.0))
    set(ParaMonte_aware_compiler true)
elseif(gnu_compiler AND (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 7.0.0))
    set(ParaMonte_aware_compiler true)
else()
    set(ParaMonte_aware_compiler false)
endif()

if (NOT ParaMonte_aware_compiler)
    message ( FATAL_ERROR 
            " \n"
            "${pmfatal}\n"
            "   Building ParaMonte minimally requires either:\n"
            "       1. Intel Fortran compiler version 18.0.0 or newer, or\n"
            "       2. GNU Fortran compiler version 7.0.0 or newer.\n"
            "   Please install one of the above two comiler suites and rebuild."
            )
endif()

#if(gnu_compiler AND (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 5.4.0))
#  message( STATUS "Disabling optimization flags due to GCC < 5.4 bug")
#  set(CMAKE_Fortran_FLAGS_RELEASE -O0
#    CACHE STRING "Flags used by the compiler during release builds." FORCE)
#  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -DNDEBUG -O0"
#    CACHE STRING "Flags used by the compiler during release builds with debug info" FORCE)
#  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O0")
#endif()
