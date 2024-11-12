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

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Preset the desired MPI library. This must be done before calling project() or enable_language().
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#   This script predefines `MPIEXEC_EXECUTABLE` <b>if and only if</b> the optional CMake argument `me` is specified by the user.
#   This is crucial for enforcing the use of the specified MPI library over the default choice of CMake and must happen prior to project() initiation.

if (DEFINED me)
    if (DEFINED MPIEXEC_EXECUTABLE AND NOT "${me}" STREQUAL "${MPIEXEC_EXECUTABLE}")
        message(WARNING
                "\n"
                "${pmwarn} The CMake variable `MPIEXEC_EXECUTABLE` and the ParaMonte CMake argument `me` are both simultaneously defined but different.\n"
                "${pmwarn} The CMake variable `MPIEXEC_EXECUTABLE` will be overwritten with the value of `me`.\n"
                "\n"
                "${pmwarn}     MPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}\n"
                "${pmwarn}     me=${me}\n"
                "\n"
                )
    endif()
    set(MPIEXEC_EXECUTABLE "${me}" CACHE FILEPATH "CMake mpiexec-executable path" FORCE)
    message(NOTICE "${pmattn} User-specified mpiexec-executable choice detected. me=${MPIEXEC_EXECUTABLE}")
elseif (DEFINED MPIEXEC_EXECUTABLE)
    set(MPIEXEC_EXECUTABLE "${MPIEXEC_EXECUTABLE}" CACHE FILEPATH "CMake mpiexec executable" FORCE)
    message(NOTICE "${pmattn} Preset mpiexec-executable choice detected. MPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}")
endif()

####################################################################################################################################

if ("${CAF_ENABLED}" OR "${MPI_ENABLED}")

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Find MPI and set some flags so that FC and CC can point to gfortran and gcc
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # If the user passes FC=mpifort etc. check and prefer that location

    get_filename_component(FTN_COMPILER_NAME "${CMAKE_Fortran_COMPILER}" NAME)
    get_filename_component(FTN_COMPILER_DIR "${CMAKE_Fortran_COMPILER}" REALPATH)
    if (FTN_COMPILER_NAME MATCHES "^[mM][pP][iI]")
        message(DEPRECATION "${pmwarn} Setting the Fortran compiler to an MPI wrapper script is deprecated!")
        set (MPI_Fortran_COMPILER "${CMAKE_Fortran_COMPILER}")
    endif()

   #get_filename_component(C_COMPILER_NAME "${CMAKE_C_COMPILER}" NAME)
   #get_filename_component(C_COMPILER_DIR "${CMAKE_C_COMPILER}" REALPATH)
   #if (C_COMPILER_NAME MATCHES "^[mM][pP][iI]")
   #    message( DEPRECATION "${pmattn} Setting the Fortran compiler to an MPI wrapper script is deprecated!")
   #    set (MPI_C_COMPILER "${CMAKE_C_COMPILER}")
   #endif()

    #if (NOT MPI_Fortran_FOUND)
        find_package(MPI)
    #endif()

    if (NOT (MPI_Fortran_FOUND AND MPIEXEC_EXECUTABLE))
    #if ((NOT MPI_C_FOUND) OR (NOT MPI_Fortran_FOUND) OR (NOT MPIEXEC_EXECUTABLE) )
        # Get default install location of MPICH from install.sh
        message(NOTICE
                "${pmwarn} Could not find all MPI components!\n"
                "${pmwarn} MPI_C_FOUND          = ${MPI_C_FOUND}\n"
                "${pmwarn} MPI_Fortran_FOUND    = ${MPI_Fortran_FOUND}\n"
                "${pmwarn} MPIEXEC_EXECUTABLE   = ${MPIEXEC_EXECUTABLE}"
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

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Test for consistent MPI environment
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if (NOT MPIEXEC_EXECUTABLE)
        message(NOTICE
                "${pmwarn} CMake failed to find `mpiexec` or similar. If building with `./install.sh`\n"
                "${pmwarn} please report this bug to the ParaMonte developers at\n"
                "${pmwarn} \n"
                "${pmwarn}     https://github.com/cdslaborg/paramonte/issues, \n"
                "${pmwarn} \n"
                "${pmwarn} otherwise point CMake to the desired MPI runtime."
                )
    #else()
    #  add_definitions(-DHAVE_MPI)
    endif()

    message (NOTICE "${pmattn} MPI_Fortran_FOUND: ${MPI_Fortran_FOUND}")
    message (NOTICE "${pmattn} MPI_Fortran_COMPILER: ${MPI_Fortran_COMPILER}")
    message (NOTICE "${pmattn} MPI_Fortran_COMPILE_FLAGS: ${MPI_Fortran_COMPILE_FLAGS}")
    message (NOTICE "${pmattn} MPI_Fortran_INCLUDE_DIRS=${MPI_Fortran_INCLUDE_DIRS}")
    message (NOTICE "${pmattn} MPI_Fortran_LINK_FLAGS: ${MPI_Fortran_LINK_FLAGS}")
    message (NOTICE "${pmattn} MPI_Fortran_LIBRARIES: ${MPI_Fortran_LIBRARIES}")
    message (NOTICE "${pmattn} MPI_C_FOUND: ${MPI_C_FOUND}")
    message (NOTICE "${pmattn} MPI_C_COMPILER: ${MPI_C_COMPILER}")
    message (NOTICE "${pmattn} MPI_C_COMPILE_FLAGS: ${MPI_C_COMPILE_FLAGS}")
    message (NOTICE "${pmattn} MPI_C_LINK_FLAGS: ${MPI_C_LINK_FLAGS}")
    message (NOTICE "${pmattn} MPI_C_LIBRARIES: ${MPI_C_LIBRARIES}")
    message (NOTICE "${pmattn} MPIEXEC_EXECUTABLE: ${MPIEXEC_EXECUTABLE}")
    message (NOTICE "${pmattn} MPI_Fortran_VERSION: ${MPI_Fortran_VERSION}")
    message (NOTICE "${pmattn} MPI_C_VERSION: ${MPI_C_VERSION}")
    message (NOTICE "${pmattn} MPI_VERSION: ${MPI_C_VERSION}")

    get_filename_component(MPIEXEC_DIR "${MPIEXEC_EXECUTABLE}" DIRECTORY)
    get_filename_component(MPIFC_DIR "${MPI_Fortran_COMPILER}" DIRECTORY)
    #get_filename_component(MPICC_DIR "${MPI_C_COMPILER}" DIRECTORY)

    #get_filename_component(MPIEXEC_RELATIVE_LOC "${MPIEXEC_EXECUTABLE}" PROGRAM)
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
        message (NOTICE "${pmattn} MPI runtime and compile time environments appear to be consistent.")
    else()
        message (NOTICE
                "${pmwarn} mpiexec executable path is ${MPIEXEC_DIR},\n"
                "${pmwarn} which differs from the location of MPICC and/or MPIFC which are respectively in\n"
                "${pmwarn} \"${MPICC_DIR}\"\n"
                "${pmwarn} \"${MPIFC_DIR}\"\n"
                "${pmwarn} This is likely indicative of a problem. If building with `./install.sh` please report\n"
                "${pmwarn} this to the ParaMonte developers by filing a new issue at:\n"
                "${pmwarn} https://github.com/cdslaborg/paramonte/issues/new"
                )
    endif()

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Work around bug #317 present on fedora systems
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if( (MPI_C_LINK_FLAGS MATCHES "noexecstack") OR (MPI_Fortran_LINK_FLAGS MATCHES "noexecstack") )
        message (NOTICE
                "\n"
                "${pmattn} The `noexecstack` linker flag was found in the MPI_<lang>_LINK_FLAGS variable. This is\n"
                "${pmattn} known to cause segmentation faults for some Fortran codes. See, e.g.,\n"
                "${pmwarn} \n"
                "${pmattn} https://gcc.gnu.org/bugzilla/show_bug.cgi?id=71729\n"
                "${pmattn} https://github.com/cdslaborg/paramonte/issues/317\n"
                "${pmwarn} \n"
                "${pmattn} `noexecstack` is being replaced with `execstack`\n"
                "\n"
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
    if (${csid_is_gnu})
        set(CMAKE_Fortran_FLAGS "-ffree-form ${MPI_Fortran_COMPILE_FLAGS} ${MPI_Fortran_LINK_FLAGS}" )
    elseif(${csid_is_intel})
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
        use mpi !mpi_f08
        implicit none
        integer :: ierr, mpi_world_size, mpi_world_rank, res_len
        character(len=MPI_MAX_PROCESSOR_NAME) :: proc
        call mpi_init(ierr)
        call mpi_comm_size(MPI_COMM_WORLD, mpi_world_size, ierr)
        call mpi_comm_rank(MPI_COMM_WORLD, mpi_world_rank, ierr)
        call mpi_get_processor_name(proc, res_len, ierr)
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

    if (NOT MPI_Fortran_MODULE_COMPILES )
    #    add_definitions(-DMPI_WORKING_MODULE)
    #else()
        message(NOTICE
                "${pmwarn} The Fortran MPI compiler wrapper appears to not be working.\n"
                "${pmwarn} For ParaMonte-compatible compilers, this may be irrelevant.\n"
                "${pmwarn} It is possible that the build will succeed, despite this fishy circumstance."
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
    set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} ${MPI_C_LINK_FLAGS}")
    set(CMAKE_C_COMPILE_FLAGS "${CMAKE_C_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS}")
    set(CMAKE_Fortran_COMPILE_FLAGS "${CMAKE_Fortran_COMPILE_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
    set(CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} ${MPI_Fortran_LINK_FLAGS}")

    #include_directories(SYSTEM ${MPI_Fortran_INCLUDE_DIRS})
    #set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}" "${MPI_Fortran_INCLUDE_DIRS}") This is now added to the paramonte target.
    set(MPIEXEC_EXECUTABLE "${MPIEXEC_EXECUTABLE}" CACHE BOOL "mpiexec path" FORCE)
    set(MPIEXEC_EXECUTABLE "${MPIEXEC_EXECUTABLE}" CACHE BOOL "mpiexec path" FORCE)
    set(MPI_Fortran_FOUND "${MPI_Fortran_FOUND}" CACHE BOOL "MPI found." FORCE)

    #if (${csid_is_intel} AND ${MPI_ENABLED})
    #    set(MPI_LINK_FLAGS "${MPI_LINK_FLAGS}" "-static_mpi")
    #    set(MPI_COMPILER_FLAGS "${MPI_COMPILER_FLAGS}" "-static_mpi")
    #    set(MPI_Fortran_LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}" "-static_mpi")
    #    set(MPI_Fortran_COMPILE_OPTIONS "${MPI_Fortran_COMPILE_OPTIONS}" "-static_mpi")
    #endif()

endif() # CAF_ENABLED OR MPI_ENABLED
