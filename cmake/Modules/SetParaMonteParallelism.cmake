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
# set Coarray Fortran (CAF) parallelization model
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#        none:  No Coarray Parallelization will be performed. Serial library will be built.   
#      single:  Coarray Parallelism will be invoked. However, only one image will perform the tasks, and there is no parallel communications.
#      shared:  Coarray Parallelism will be invoked. It causes the underlying Intel® Message Passing Interface (MPI) parallelization
#               to run on one node with multiple cores or processors with shared memory.
# distributed:  This option requires a special license to be installed for Intel compiler suite, and is only available on Linux systems, although you can specify it here.
#               On the Linux systems, it causes the underlying Intel® MPI Library parallelization to run in a multi-node environment (multiple CPUs with distributed memory).

set ( CAFTYPE "none" CACHE STRING "Coarray Fortran configuration to build (none/single/shared/distributed)" )

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set number of Fortran Coarray (MPI) images (that will run in parallel, if Coarray/MPI parallelism is requested by the user)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set ( FOR_COARRAY_NUM_IMAGES 3 CACHE STRING "Number of Fortran Coarray/MPI images/processes" )

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set other parallelization model flags
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#set(MPI_ENABLED FALSE CACHE BOOL "Enable build with MPI-parallelism")
#set(OMP_ENABLED FALSE CACHE BOOL "Enable build with MPI-parallelism" FORCE)
option(MPI_ENABLED "Enable build with MPI-parallelism" OFF)
option(OMP_ENABLED "Enable build with OpenMP-parallelism (unsupported)" OFF)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set Fortran/C compiler version: compiler version matters when optimizations are enabled as object files will be different across different versions
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#set ( COMPILER_VERSION "" CACHE STRING "Fortran/C compiler version" )

include(SetParaMonteFlags)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

if (CAF_ENABLED OR MPI_ENABLED)

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Find MPI and set some flags so that FC and CC can point to gfortran and gcc
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # If the user passes FC=mpifort etc. check and prefer that location

    get_filename_component( FTN_COMPILER_NAME "${CMAKE_Fortran_COMPILER}" NAME )
    get_filename_component( FTN_COMPILER_DIR "${CMAKE_Fortran_COMPILER}" REALPATH )
    if (FTN_COMPILER_NAME MATCHES "^[mM][pP][iI]")
        message( DEPRECATION "${pmattn} Setting the Fortran compiler to an MPI wrapper script is deprecated!")
        set (MPI_Fortran_COMPILER "${CMAKE_Fortran_COMPILER}")
    endif()

   #get_filename_component( C_COMPILER_NAME "${CMAKE_C_COMPILER}" NAME )
   #get_filename_component( C_COMPILER_DIR "${CMAKE_C_COMPILER}" REALPATH )
   #if (C_COMPILER_NAME MATCHES "^[mM][pP][iI]")
   #    message( DEPRECATION "${pmattn} Setting the Fortran compiler to an MPI wrapper script is deprecated!")
   #    set (MPI_C_COMPILER "${CMAKE_C_COMPILER}")
   #endif()

    find_package( MPI )

    if ( (NOT MPI_Fortran_FOUND) OR (NOT MPIEXEC) )
    #if ( (NOT MPI_C_FOUND) OR (NOT MPI_Fortran_FOUND) OR (NOT MPIEXEC) )
        # Get default install location of MPICH from install.sh
        message ( WARNING 
                " \n"
                " ${pmwarn}\n"
                " ${pmattn} Could not find all MPI components!\n"
                " ${pmattn} MPI_C_FOUND         = ${MPI_C_FOUND}\n"
                " ${pmattn} MPI_Fortran_FOUND   = ${MPI_Fortran_FOUND}\n"
                " ${pmattn} MPIEXEC             = ${MPIEXEC}\n" 
                )
        #execute_process (   COMMAND "./install.sh" -P mpich
        #                    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        #                    OUTPUT_VARIABLE DEFAULT_MPICH_INSTALL_LOC
        #                    OUTPUT_QUIET
        #                    OUTPUT_STRIP_TRAILING_WHITESPACE
        #                )
        #find_program    (   MY_MPI_EXEC NAMES mpirun mpiexec lamexec srun
        #                    PATHS "${DEFAULT_MPICH_INSTALL_LOC}" ENV PATH
        #                    HINTS "${FTN_COMPILER_DIR}" "${C_COMPILER_DIR}"
        #                    PATH_SUFFIXES bin
        #                )
        #set ( MPI_HOME "${MPI_HOME}" "${MY_MPI_EXEC}" "${MY_MPI_EXEC}/.." )
        #find_package( MPI REQUIRED )
    endif()
    list(REMOVE_DUPLICATES MPI_Fortran_INCLUDE_PATH)

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Test for consistent MPI environment
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if (NOT MPIEXEC)
        message ( WARNING 
                " \n"
                " ${pmwarn}\n"
                " ${pmattn} CMake failed to find `mpiexec` or similar. If building with `./install.sh`\n"
                " ${pmattn} please report this bug to the ParaMonte developers at\n"
                " ${pmattn} https://github.com/cdslaborg/ParaMonte/issues, \n"
                " ${pmattn} otherwise point CMake to the desired MPI runtime.\n"
                )
        #message ( WARNING 
        #        " \n"
        #        " ${pmwarn}\n"
        #        " ${pmattn} Building ParaMonte in serial mode..." 
        #        )
    #else()
    #  add_definitions(-DHAVE_MPI)
    endif()

    message ( STATUS "")
    message ( STATUS "${pmattn} MPI_Fortran_FOUND: ${MPI_Fortran_FOUND}")
    message ( STATUS "${pmattn} MPI_Fortran_COMPILER: ${MPI_Fortran_COMPILER}")
    message ( STATUS "${pmattn} MPI_Fortran_COMPILE_FLAGS: ${MPI_Fortran_COMPILE_FLAGS}")
    message ( STATUS "${pmattn} MPI_Fortran_LINK_FLAGS: ${MPI_Fortran_LINK_FLAGS}")
    message ( STATUS "${pmattn} MPI_Fortran_LIBRARIES: ${MPI_Fortran_LIBRARIES}")
    message ( STATUS "${pmattn} MPI_C_FOUND: ${MPI_C_FOUND}")
    message ( STATUS "${pmattn} MPI_C_COMPILER: ${MPI_C_COMPILER}")
    message ( STATUS "${pmattn} MPI_C_COMPILE_FLAGS: ${MPI_C_COMPILE_FLAGS}")
    message ( STATUS "${pmattn} MPI_C_LINK_FLAGS: ${MPI_C_LINK_FLAGS}")
    message ( STATUS "${pmattn} MPI_C_LIBRARIES: ${MPI_C_LIBRARIES}")
    message ( STATUS "${pmattn} MPIEXEC: ${MPIEXEC}")
    message ( STATUS "")

    get_filename_component(MPIEXEC_DIR "${MPIEXEC}" DIRECTORY)
    get_filename_component(MPIFC_DIR "${MPI_Fortran_COMPILER}" DIRECTORY)
   #get_filename_component(MPICC_DIR "${MPI_C_COMPILER}" DIRECTORY)

    #get_filename_component(MPIEXEC_RELATIVE_LOC "${MPIEXEC}" PROGRAM)
    #get_filename_component(MPIEXEC_ABS_LOC "${MPIEXEC_RELATIVE_LOC}" REALPATH)
    #get_filename_component(MPIEXEC_DIR "${MPIEXEC_ABS_LOC}" DIRECTORY)
    #get_filename_component(MPICC_RELATIVE_LOC "${MPI_C_COMPILER}" PROGRAM)
    #get_filename_component(MPICC_ABS_LOC "${MPICC_RELATIVE_LOC}" REALPATH)
    #get_filename_component(MPICC_DIR "${MPICC_ABS_LOC}" DIRECTORY)
    #get_filename_component(MPIFC_RELATIVE_LOC "${MPI_Fortran_COMPILER}" PROGRAM)
    #get_filename_component(MPIFC_ABS_LOC "${MPIFC_RELATIVE_LOC}" REALPATH)
    #get_filename_component(MPIFC_DIR "${MPIFC_ABS_LOC}" DIRECTORY)

    if (MPIEXEC_DIR STREQUAL MPIFC_DIR)
    #if ((MPIEXEC_DIR STREQUAL MPICC_DIR) AND (MPIEXEC_DIR STREQUAL MPIFC_DIR))
        message ( STATUS "${pmattn} MPI runtime and compile time environments appear to be consistent")
    else()
        message ( WARNING 
                " \n"
                " ${pmwarn}\n"
                " ${pmattn} MPIEXEC is in ${MPIEXEC_DIR},\n"
                " ${pmattn} which differs from the location of MPICC and/or MPIFC which are in\n"
                " ${pmattn} ${MPICC_DIR} and \n"
                " ${pmattn} ${MPIFC_DIR}, respectively.\n"
                " ${pmattn} This is likely indicative of a problem. If building with `./install.sh` please report\n"
                " ${pmattn} this to the ParaMonte developers by filing a new issue at:\n"
                " ${pmattn} https://github.com/cdslaborg/ParaMonte/issues/new"
                )
    endif()

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Work around bug #317 present on fedora systems
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if( (MPI_C_LINK_FLAGS MATCHES "noexecstack") OR (MPI_Fortran_LINK_FLAGS MATCHES "noexecstack") )
        message ( WARNING
                " \n"
                " ${pmwarn}\n"
                " ${pmattn} The `noexecstack` linker flag was found in the MPI_<lang>_LINK_FLAGS variable. This is\n"
                " ${pmattn} known to cause segmentation faults for some Fortran codes. See, e.g.,\n"
                " ${pmattn} https://gcc.gnu.org/bugzilla/show_bug.cgi?id=71729 or\n"
                " ${pmattn} https://github.com/cdslaborg/ParaMonte/issues/317.\n"
                " ${pmattn} `noexecstack` is being replaced with `execstack`\n"
                )
        string(REPLACE "noexecstack" "execstack" MPI_C_LINK_FLAGS_FIXED ${MPI_C_LINK_FLAGS})
        string(REPLACE "noexecstack" "execstack" MPI_Fortran_LINK_FLAGS_FIXED ${MPI_Fortran_LINK_FLAGS})
        set(MPI_C_LINK_FLAGS "${MPI_C_LINK_FLAGS_FIXED}" CACHE STRING "MPI C linking flags" FORCE)
        set(MPI_Fortran_LINK_FLAGS "${MPI_Fortran_LINK_FLAGS_FIXED}" CACHE STRING "MPI Fortran linking flags" FORCE)
    endif()

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Make sure a simple "hello world" C mpi program compiles
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#    set(OLD_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
#    set(CMAKE_REQUIRED_FLAGS ${MPI_C_COMPILE_FLAGS} ${MPI_C_LINK_FLAGS})
#    set(OLD_INCLUDES ${CMAKE_REQUIRED_INCLUDES})
#    set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
#    set(OLD_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
#    set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
#    include (CheckCSourceCompiles)
#    CHECK_C_SOURCE_COMPILES("
#    #include <mpi.h>
#    #include <stdio.h>
#    int main(int argc, char** argv) {
#        MPI_Init(NULL, NULL);
#        int world_size;
#        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#        int world_rank;
#        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#        char processor_name[MPI_MAX_PROCESSOR_NAME];
#        int name_len;
#        MPI_Get_processor_name(processor_name, &name_len);
#        printf('Hello world from processor %s, rank %d out of %d processors',
#                processor_name, world_rank, world_size);
#        MPI_Finalize();
#    }
#    "
#    MPI_C_COMPILES )
#
#    set(CMAKE_REQUIRED_FLAGS ${OLD_REQUIRED_FLAGS})
#    set(CMAKE_REQUIRED_INCLUDES ${OLD_INCLUDES})
#    set(CMAKE_REQUIRED_LIBRARIES ${OLD_LIBRARIES})
#    unset(OLD_REQUIRED_FLAGS)
#    unset(OLD_INCLUDES)
#    unset(OLD_LIBRARIES)
#
#    if (NOT MPI_C_COMPILES)
#        message (FATAL_ERROR 
#                " \n"
#                "${pmfatal}\n"
#                "   MPI_C is missing!\n"
#                "   Try setting MPI_C_COMPILER to the appropriate C compiler wrapper script and reconfigure.\n"
#                "   i.e., `cmake -DMPI_C_COMPILER=/path/to/mpicc ..` or set it by editing the cache using\n"
#                "   cmake-gui or ccmake.\n"
#                )
#    endif()

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Make sure a simple "hello world" Fortran mpi program compiles.
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    set(OLD_REQUIRED_FLAGS ${CMAKE_Fortran_FLAGS})
    if (gnu_compiler)
        set(CMAKE_Fortran_FLAGS "-ffree-form ${MPI_Fortran_COMPILE_FLAGS} ${MPI_Fortran_LINK_FLAGS}" )
    elseif(intel_compiler)
        set(CMAKE_Fortran_FLAGS "-free ${MPI_Fortran_COMPILE_FLAGS} ${MPI_Fortran_LINK_FLAGS}" )
    endif()
    string(REPLACE ";" " " CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS})
    set(OLD_INCLUDES ${CMAKE_REQUIRED_INCLUDES})
    set(CMAKE_REQUIRED_INCLUDES ${MPI_Fortran_INCLUDE_PATH})
    set(OLD_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
    set(CMAKE_REQUIRED_LIBRARIES ${MPI_Fortran_LIBRARIES})

    include (CheckFortranSourceCompiles)
    CHECK_Fortran_SOURCE_COMPILES(
    "
    program mpi_hello
        use mpi
        implicit none
        integer :: ierr, mpi_world_size, mpi_world_rank, res_len
        character(len=MPI_MAX_PROCESSOR_NAME) :: proc
        call mpi_init(ierr)
        call mpi_comm_size(MPI_COMM_WORLD,mpi_world_size,ierr)
        call mpi_comm_rank(MPI_COMM_WORLD,mpi_world_rank,ierr)
        call mpi_get_processor_name(proc,res_len,ierr)
        write(*,*) 'Hello from processor ', trim(proc), ' rank ', mpi_world_rank, ' out of ', mpi_world_size, '.'
        call mpi_finalize(ierr)
    end program
    "
    MPI_Fortran_MODULE_COMPILES
    SRC_EXT F90
    )
    set(CMAKE_Fortran_FLAGS ${OLD_REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_INCLUDES ${OLD_INCLUDES})
    set(CMAKE_REQUIRED_LIBRARIES ${OLD_LIBRARIES})
    unset(OLD_REQUIRED_FLAGS)
    unset(OLD_INCLUDES)
    unset(OLD_LIBRARIES)

    if ( NOT MPI_Fortran_MODULE_COMPILES )
    #    add_definitions(-DMPI_WORKING_MODULE)
    #else()
        message ( WARNING 
                " \n"
                " ${pmwarn}\n"
                " ${pmattn} It appears that the Fortran MPI compiler is not working.\n"
                " ${pmattn} For ParaMonte-compatible compilers, this may be irrelavent.\n"
                " ${pmattn} It is possible that the build will succeed, despite this fishy circumstance.\n"
                )
    endif()

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Setup MPI flags
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if(NOT CMAKE_C_COMPILE_FLAGS)
        set(CMAKE_C_COMPILE_FLAGS "")
    endif()
    if(NOT CMAKE_C_LINK_FLAGS)
        set(CMAKE_C_LINK_FLAGS "")
    endif()
    if(NOT CMAKE_Fortran_COMPILE_FLAGS)
        set(CMAKE_Fortran_COMPILE_FLAGS "")
    endif()
    if(NOT CMAKE_Fortran_LINK_FLAGS)
        set(CMAKE_Fortran_LINK_FLAGS "")
    endif()
    set(CMAKE_C_COMPILE_FLAGS "${CMAKE_C_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS}")
    set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} ${MPI_C_LINK_FLAGS}")
    set(CMAKE_Fortran_COMPILE_FLAGS "${CMAKE_Fortran_COMPILE_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
    set(CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} ${MPI_Fortran_LINK_FLAGS}")

    include_directories(SYSTEM ${MPI_INCLUDE_PATH})

endif() # CAF_ENABLED OR MPI_ENABLED
