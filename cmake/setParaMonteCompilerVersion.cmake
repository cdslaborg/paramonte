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

#   This script defines the following variables:
#   `CMAKE_Fortran_COMPILER_VERSION` : (if undefined already) the full version of the compiler in `major[.minor[.patch[.tweak]]]` format.

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set CMAKE_Fortran_COMPILER_VERSION if CMake doesn't do it for us. This must happen after calling project().
# UPDATE: Set the CMAKE_Fortran_COMPILER_VERSION even if CMake sets it, because stupid CMake sets it to
# a version number different from the output of `compiler_version()` or even `ifort --version`.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (CMAKE_Fortran_COMPILER_VERSION)
    if (NOT (CMAKE_VERSION VERSION_LESS 3.3.1))
        message(NOTICE
                "${pmwarn} Preset CMAKE_Fortran_COMPILER_VERSION detected: ${CMAKE_Fortran_COMPILER_VERSION}\n"
                "${pmwarn} The value of CMAKE_Fortran_COMPILER_VERSION will be overwritten with the version inferred from `compiler_options()`."
                )
    endif()
endif()

#   No CMAKE_Fortran_COMPILER_VERSION set, build our own
#   Try extracting it directly from ISO_Fortran_ENV's `compiler_version()`.

file(WRITE "${CMAKE_BINARY_DIR}/getCompilerVersion.F90"
    "
    use iso_Fortran_env, only: compiler_version, output_unit
    write(output_unit,'(a)') compiler_version()
    end
    "
    )

try_run(PROGRAM_RUN_NOTICE
        PROGRAM_COMPILE_SUCCEEDED
        "${CMAKE_BINARY_DIR}" # Binary dir
        "${CMAKE_BINARY_DIR}/getCompilerVersion.F90" # source file
        RUN_OUTPUT_VARIABLE COMPILER_VERSION_STRING
        )

if (PROGRAM_RUN_NOTICE STREQUAL FAILED_TO_RUN OR NOT PROGRAM_COMPILE_SUCCEEDED)
    message(NOTICE
            "\n"
            "${pmwarn} \n"
            "${pmwarn} Could not reliably detect Fortran compiler version.\n"
            "${pmwarn} It will be inferred it from the C compiler if it matches the Fortran compiler ID."
            )
    if("${CMAKE_C_COMPILER_ID}" MATCHES "${CMAKE_Fortran_COMPILER_ID}")
        set(COMPILER_VERSION_NUMBER "${CMAKE_C_COMPILER_VERSION}")
    else()
        set(COMPILER_VERSION_NUMBER "")
        message(NOTICE
                "\n"
                "${pmwarn} \n"
                "${pmwarn} Exhausted all possible means of detecting the Fortran compiler version.\n"
                "${pmwarn} Please install of the compilers supported by the library.\n"
                "${pmwarn} Proceeding with no guarantee of build success...\n"
                )
    endif()
else()
    string(REPLACE "\n" "" COMPILER_VERSION_STRING "${COMPILER_VERSION_STRING}")
    string(REGEX MATCH "[0-9]+\\.[0-9]+(\\.[0-9]+)?" COMPILER_VERSION_NUMBER "${COMPILER_VERSION_STRING}" )
    message(NOTICE "${pmattn} The full Fortran compiler version string: ${COMPILER_VERSION_STRING}" )
    message(NOTICE "${pmattn} The Fortran compiler version number: ${COMPILER_VERSION_NUMBER}" )
endif()

set(CMAKE_Fortran_COMPILER_VERSION "${COMPILER_VERSION_NUMBER}")

