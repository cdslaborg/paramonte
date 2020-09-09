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

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set the desired compiler suite
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set(COMPILER_IS_PREDEFINED FALSE)
if (DEFINED FC)
    set (CMAKE_Fortran_COMPILER "${FC}")
    set (COMPILER_IS_PREDEFINED TRUE)
endif()

# if (DEFINED MPI_Fortran_COMPILER)
    # set (MPI_Fortran_COMPILER "${MPI_Fortran_COMPILER}")
    # set (COMPILER_IS_PREDEFINED TRUE)
# endif()

if (DEFINED MPIEXEC_EXECUTABLE)
    set (MPIEXEC_EXECUTABLE "${MPIEXEC_EXECUTABLE}")
    set (COMPILER_IS_PREDEFINED TRUE)
endif()

if (COMPILER_IS_PREDEFINED)
    return()
endif()

set(PMCS "${PMCS}" CACHE STRING "The requested ParaMonte build compiler suite" FORCE)
unset (requested_compiler_suite)
if( DEFINED PMCS AND (NOT ${PMCS} STREQUAL "") )
    unset (CMAKE_C_COMPILER)
    unset (CMAKE_CXX_COMPILER)
    unset (CMAKE_Fortran_COMPILER)
    if("${PMCS}" MATCHES "[Ii][Nn][Tt][Ee][Ll]")
        set (requested_compiler_suite "Intel")
        message ( STATUS "${pmnote} requested_compiler_suite=${requested_compiler_suite}" )
    elseif("${PMCS}" MATCHES "[Gg][Nn][Uu]")
        set (requested_compiler_suite "GNU")
        message ( STATUS "${pmnote} requested_compiler_suite=${requested_compiler_suite}" )
    else()
        message ( FATAL_ERROR 
                " \n"
                "${pmfatal}\n"
                "   The requested compiler suite for ParaMonte library build as specified\n"
                "   by the environment variable PMCS=${PMCS} is not any of the two supported\n"
                "   comipler suites: Intel or GNU. Please ensure one of these compiler suites\n"
                "   is installed on your system, then provide its name to cmake by passing\n"
                " \n"
                "       -DPMCS=Intel\n"
                " \n"
                "       or\n"
                " \n"
                "       -DPMCS=GNU\n"
                " \n"
                "   when invoking cmake."
                )
    endif()
endif()


if ( (NOT DEFINED requested_compiler_suite) OR (requested_compiler_suite STREQUAL "Intel") )

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # detect intel compilers
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # Intel icc

    execute_process(COMMAND which icc OUTPUT_VARIABLE icc_path RESULT_VARIABLE icc_errstat ERROR_QUIET)
    string(REGEX REPLACE "\n$" "" icc_path "${icc_path}") # strip new line char from the end
    if (icc_errstat STREQUAL "0")

        message( STATUS "${pmattn} Intel C compiler detected at: ${icc_path}" )
        set(CMAKE_C_COMPILER "${icc_path}" CACHE FILEPATH "Intel C compiler" FORCE)

        execute_process(COMMAND which mpiicc OUTPUT_VARIABLE mpiicc_path RESULT_VARIABLE mpiicc_errstat ERROR_QUIET)
        string(REGEX REPLACE "\n$" "" mpiicc_path "${mpiicc_path}") # strip new line char from the end
        if (mpiicc_errstat STREQUAL "0")
            message( STATUS "${pmattn} Intel MPI C wrapper detected at: ${mpiicc_path}" )
            set(MPI_C_COMPILER "${mpiicc_path}" CACHE FILEPATH "Intel MPI C wrapper" FORCE)
        else()
            message ( WARNING 
                    " \n"
                    " ${pmwarn}\n"
                    " ${pmattn} Intel MPI C wrapper could not be found."
                    )
        endif()

    else()

        message ( WARNING
                " \n"
                " ${pmwarn}\n"
                " ${pmattn} Intel C compiler could not be found."
                )

    endif()

    # Intel icpc

    execute_process(COMMAND which icpc OUTPUT_VARIABLE icpc_path RESULT_VARIABLE icpc_errstat ERROR_QUIET)
    string(REGEX REPLACE "\n$" "" icpc_path "${icpc_path}") # strip new line char from the end
    if (icpc_errstat STREQUAL "0")

        message( STATUS "${pmattn} Intel C++ compiler detected at: ${icpc_path}" )
        set(CMAKE_CXX_COMPILER "${icpc_path}" CACHE FILEPATH "Intel C++ compiler" FORCE)

        execute_process(COMMAND which mpiicpc OUTPUT_VARIABLE mpiicpc_path RESULT_VARIABLE mpiicpc_errstat ERROR_QUIET)
        string(REGEX REPLACE "\n$" "" mpiicpc_path "${mpiicpc_path}") # strip new line char from the end
        if (mpiicpc_errstat STREQUAL "0")
            message( STATUS "${pmattn} Intel MPI C++ wrapper detected at: ${mpiicpc_path}" )
            set(MPI_CXX_COMPILER "${mpiicpc_path}" CACHE FILEPATH "Intel MPI C++ wrapper" FORCE)
        else()
            message ( WARNING 
                    " \n"
                    " ${pmwarn}\n"
                    " ${pmattn} Intel MPI C++ wrapper could not be found."
                    )
        endif()

    else()

        message ( WARNING
                " \n"
                " ${pmwarn}\n"
                " ${pmattn} Intel C++ compiler could not be found."
                )

    endif()

    # Intel ifort

    execute_process(COMMAND which ifort OUTPUT_VARIABLE ifort_path RESULT_VARIABLE ifort_errstat ERROR_QUIET)
    string(REGEX REPLACE "\n$" "" ifort_path "${ifort_path}") # strip new line char from the end
    if (ifort_errstat STREQUAL "0")

        message( STATUS "${pmattn} Intel Fortran compiler detected at: ${ifort_path}" )
        set(CMAKE_Fortran_COMPILER "${ifort_path}" CACHE FILEPATH "Intel Fortran compiler" FORCE)

        execute_process(COMMAND which mpiifort OUTPUT_VARIABLE mpiifort_path RESULT_VARIABLE mpiifort_errstat ERROR_QUIET)
        string(REGEX REPLACE "\n$" "" mpiifort_path "${mpiifort_path}") # strip new line char from the end
        if (mpiifort_errstat STREQUAL "0")
            message( STATUS "${pmattn} Intel MPI Fortran wrapper detected at: ${mpiifort_path}" )
            set(MPI_Fortran_COMPILER "${mpiifort_path}" CACHE FILEPATH "Intel MPI Fortran wrapper" FORCE)
        else()
            message ( WARNING
                    " \n"
                    " ${pmwarn}\n"
                    " ${pmattn} Intel MPI Fortran wrapper could not be found."
                    )
        endif()

    else()

        message ( WARNING
                " \n"
                " ${pmwarn}\n"
                " ${pmattn} Intel Fortran compiler could not be found."
                )

    endif()

    # check if all components have been detected

    unset( intel_compiler )
    if (icc_errstat         STREQUAL "0" AND 
        icpc_errstat        STREQUAL "0" AND 
        ifort_errstat       STREQUAL "0" AND 
        mpiicc_errstat      STREQUAL "0" AND 
        mpiicpc_errstat     STREQUAL "0" AND 
        mpiifort_errstat    STREQUAL "0" )
        set( COMPILER_SUITE "Intel" )
        set( CMAKE_Fortran_COMPILER_ID "Intel" )
        set( intel_compiler TRUE)
    elseif (requested_compiler_suite STREQUAL "Intel")
        message ( WARNING
                " \n"
                " ${pmwarn}\n"
                " ${pmattn} Failed to detect all required components of the Intel Compiler Collection.\n"
                " ${pmattn} ParaMonte build will proceed, but there is no guarantee of a\n"
                " ${pmattn} successful build via Intel compiler suite."
                )
    endif()

elseif ( (NOT DEFINED requested_compiler_suite) OR (requested_compiler_suite STREQUAL "GNU") )

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # detect GNU compilers
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # GNU C compiler

    execute_process(COMMAND which gcc OUTPUT_VARIABLE gcc_path RESULT_VARIABLE gcc_errstat ERROR_QUIET)
    string(REGEX REPLACE "\n$" "" gcc_path "${gcc_path}") # strip new line char from the end
    if (gcc_errstat STREQUAL "0")

        message( STATUS "${pmattn} GNU C compiler detected at: ${gcc_path}" )
        set(CMAKE_C_COMPILER "${gcc_path}" CACHE FILEPATH "GNU C compiler" FORCE)

        execute_process(COMMAND which mpicc OUTPUT_VARIABLE mpicc_path RESULT_VARIABLE mpicc_errstat ERROR_QUIET)
        string(REGEX REPLACE "\n$" "" mpicc_path "${mpicc_path}") # strip new line char from the end
        if (mpicc_errstat STREQUAL "0")
            message( STATUS "${pmattn} GNU MPI C wrapper detected at: ${mpicc_path}" )
            set(MPI_C_COMPILER "${mpicc_path}" CACHE FILEPATH "GNU MPI C wrapper" FORCE)
        else()
            message ( WARNING 
                    " \n"
                    " ${pmwarn}\n"
                    " ${pmattn} GNU MPI C wrapper could not be found."
                    )
        endif()

    else()

        message ( WARNING
                " \n"
                " ${pmwarn}\n"
                " ${pmattn} GNU C compiler could not be found."
                )

    endif()

    # GNU C++ compiler

    execute_process(COMMAND which g++ OUTPUT_VARIABLE gcxx_path RESULT_VARIABLE gcxx_errstat ERROR_QUIET)
    string(REGEX REPLACE "\n$" "" gcxx_path "${gcxx_path}") # strip new line char from the end
    if (gcxx_errstat STREQUAL "0")

        message( STATUS "${pmattn} GNU C++ compiler detected at: ${gcxx_path}" )
        set(CMAKE_CXX_COMPILER "${gcxx_path}" CACHE FILEPATH "GNU C++ compiler" FORCE)

        execute_process(COMMAND which mpicxx OUTPUT_VARIABLE mpicxx_path RESULT_VARIABLE mpicxx_errstat ERROR_QUIET)
        string(REGEX REPLACE "\n$" "" mpicxx_path "${mpicxx_path}") # strip new line char from the end
        if (mpicxx_errstat STREQUAL "0")
            message( STATUS "${pmattn} GNU MPI C++ wrapper detected at: ${mpicxx_path}" )
            set(MPI_CXX_COMPILER "${mpicxx_path}" CACHE FILEPATH "GNU MPI C++ wrapper" FORCE)
        else()
            message ( WARNING 
                    " \n"
                    " ${pmwarn}\n"
                    " ${pmattn} GNU MPI C++ wrapper could not be found."
                    )
        endif()

    else()

        message ( WARNING
                " \n"
                " ${pmwarn}\n"
                " ${pmattn} GNU C++ compiler could not be found."
                )

    endif()

    # GNU Fortran compiler

    if( DEFINED CAFTYPE AND (
        "${CAFTYPE}" MATCHES "single" OR
        "${CAFTYPE}" MATCHES "shared" OR
        "${CAFTYPE}" MATCHES "distributed"
        ))
        set(CAF_ENABLED ON CACHE BOOL "Enable Coarray Fortran parallelism" FORCE)
        execute_process(COMMAND which caf OUTPUT_VARIABLE gfort_path RESULT_VARIABLE gfort_errstat ERROR_QUIET)
        if (NOT gfort_errstat STREQUAL "0")
        message ( WARNING
                " \n"
                " ${pmwarn}\n"
                " ${pmattn} Failed to detect Coarray-aware Fortran compiler wrapper caf.\n"
                " ${pmattn} Please ensure you have OpenCoarrays installed on your system properly.\n"
                " ${pmattn} ParaMonte may not be built correctly."
                )
        execute_process(COMMAND which gfortran OUTPUT_VARIABLE gfort_path RESULT_VARIABLE gfort_errstat ERROR_QUIET)
        endif()
    else()
        execute_process(COMMAND which gfortran OUTPUT_VARIABLE gfort_path RESULT_VARIABLE gfort_errstat ERROR_QUIET)
    endif()
    string(REGEX REPLACE "\n$" "" gfort_path "${gfort_path}") # strip new line char from the end
    if (gfort_errstat STREQUAL "0")

        message( STATUS "${pmattn} GNU Fortran compiler detected at: ${gfort_path}" )
        set(CMAKE_Fortran_COMPILER "${gfort_path}" CACHE FILEPATH "GNU Fortran compiler" FORCE)

        execute_process(COMMAND which mpifort OUTPUT_VARIABLE mpifort_path RESULT_VARIABLE mpifort_errstat ERROR_QUIET)
        string(REGEX REPLACE "\n$" "" mpifort_path "${mpifort_path}") # strip new line char from the end
        if (mpifort_errstat STREQUAL "0")
            message( STATUS "${pmattn} GNU MPI Fortran wrapper detected at: ${mpifort_path}" )
            set(MPI_Fortran_COMPILER "${mpifort_path}" CACHE FILEPATH "GNU MPI Fortran wrapper" FORCE)
        else()
            message ( WARNING
                    " \n"
                    " ${pmwarn}\n"
                    " ${pmattn} GNU MPI Fortran wrapper could not be found."
                    )
        endif()

    else()

        message ( WARNING
                " \n"
                " ${pmwarn}\n"
                " ${pmattn} GNU Fortran compiler could not be found."
                )

    endif()

    # check if all components have been detected

    unset( gnu_compiler )
    if (gcc_errstat         STREQUAL "0" AND 
        gcxx_errstat        STREQUAL "0" AND 
        gfort_errstat       STREQUAL "0" AND 
        mpicc_errstat       STREQUAL "0" AND 
        mpicxx_errstat      STREQUAL "0" AND 
        mpifort_errstat     STREQUAL "0" )
        set( COMPILER_SUITE "GNU" )
        set( CMAKE_Fortran_COMPILER_ID "GNU" )
        set( gnu_compiler TRUE)
    elseif (requested_compiler_suite STREQUAL "GNU")
        message ( WARNING
                " \n"
                " ${pmwarn}\n"
                " ${pmattn} Failed to detect all required components of the GNU Compiler Collection.\n"
                " ${pmattn} ParaMonte build will proceed, but there is no guarantee of a\n"
                " ${pmattn} successful build via GNU compiler suite."
                )
    endif()

endif()
