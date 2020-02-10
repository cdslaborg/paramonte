cmake_minimum_required(VERSION 3.10)

cmake_policy(VERSION 3.10...3.14)

if(DEFINED ENV{SOURCE_DATE_EPOCH})
  set(SOURCE_DATE_EPOCH "$ENV{SOURCE_DATE_EPOCH}")
endif()

set(CMAKE_SKIP_RPATH ON
  CACHE BOOL "Don't add a build-dir rpath")
set(CMAKE_BUILD_WITH_INSTALL_RPATH ON
  CACHE BOOL "Build using the install rpath")
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON
  CACHE BOOL "Use relative rpaths")
# Set the type/configuration of build to perform
set ( CMAKE_CONFIGURATION_TYPES "Debug" "Release" "MinSizeRel" "RelWithDebInfo" "CodeCoverage" )
set ( CMAKE_BUILD_TYPE "Release"
  CACHE STRING "Select which configuration to build." )
set_property ( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${CMAKE_CONFIGURATION_TYPES} )

# Add option and check environment to determine if developer tests should be run
set( _TF OFF )
if( NOT DEFINED ENV{OPENCOARRAYS_DEVELOPER})
  set ( ENV{OPENCOARRAYS_DEVELOPER} FALSE )
endif()
if($ENV{OPENCOARRAYS_DEVELOPER})
  set( _TF ON )
endif()
option(CAF_RUN_DEVELOPER_TESTS "Run tests intended only for developers" ${_TF})
mark_as_advanced(CAF_RUN_DEVELOPER_TESTS)
option( CAF_ENABLE_ISO_Fortran_BINDING
  "Build and install F2018 C interop iso_fortran_binding.h header; experimental!"
  ${_TF} )


# Name project and specify source languages
# Parse version from .VERSION file so that more info can be added and easier to get from scripts
file(STRINGS ".VERSION" first_line
  LIMIT_COUNT 1
  )

string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+(-rc[0-9]+)?"
  OpenCoarraysVersion "${first_line}")

if((NOT (OpenCoarraysVersion MATCHES "[0-9]+\\.[0-9]+\\.[0-9]+(-rc[0-9]+)?")) AND (EXISTS "${CMAKE_SOURCE_DIR}/.git"))
  find_package(Git QUIET)
  if(GIT_FOUND)
    execute_process(COMMAND "${GIT_EXECUTABLE}" describe --abbrev=0
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
      RESULT_VARIABLE git_status
      OUTPUT_VARIABLE git_output
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if((git_status STREQUAL "0") AND (git_output MATCHES "[0-9]+\\.[0-9]+\\.[0-9]+(-rc[0-9]+)?"))
      set(OpenCoarraysVersion "${git_output}")
    endif()
    execute_process(COMMAND "${GIT_EXECUTABLE}" describe --always
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
      RESULT_VARIABLE git_status
      OUTPUT_VARIABLE full_git_describe
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT (git_status STREQUAL "0"))
      set(full_git_describe NOTFOUND)
    endif()
    # Create a source distribution target using git archive
    # e.g., `make dist` will package a release using current git state
    add_custom_target(dist
      COMMAND "${CMAKE_COMMAND}" -P "${CMAKE_SOURCE_DIR}/cmake/makeDist.cmake" "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}"
      VERBATIM)
  else()
    message( WARNING "Could not find git executable!")
  endif()
endif()

if(NOT (OpenCoarraysVersion MATCHES "[0-9]+\\.[0-9]+\\.[0-9]+(-rc[0-9]+)?"))
  message( WARNING "Could not extract version from git, falling back on .VERSION, line 3.")
  file(STRINGS ".VERSION" OpenCoarraysVersion
  REGEX "[0-9]+\\.[0-9]+\\.[0-9]+(-rc[0-9]+)?"
  )
endif()

if(NOT full_git_describe)
  set(full_git_describe ${OpenCoarraysVersion})
endif()

string(REGEX REPLACE "-rc[0-9]+$"
       ".0" OPENCOARRAYS_CMAKE_PROJECT_VERSION
       "${OpenCoarraysVersion}")

file(READ "${CMAKE_CURRENT_LIST_DIR}/cmake/Banner.txt.in" OC_BANNER)
string(CONFIGURE "${OC_BANNER}" OC_BANNER @ONLY)
string(STRIP "${OC_BANNER}" OC_BANNER)
message(STATUS "${OC_BANNER}")

if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/.git")
  message( STATUS "Build from git repository detected")
endif()

message( STATUS "Running with CMake from: ${CMAKE_COMMAND}" )
message( STATUS "Current source dir (for OpenCoarrays): ${CMAKE_CURRENT_SOURCE_DIR}")

project(opencoarrays VERSION "${OPENCOARRAYS_CMAKE_PROJECT_VERSION}" LANGUAGES C Fortran)

if(DEFINED CMAKE_BUILD_TYPE AND "${CMAKE_BUILD_TYPE}" MATCHES "[Rr][Ee][Ll]")
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-working-directory")
  endif()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdebug-prefix-map=${CMAKE_SOURCE_DIR}=.")
  if(CMAKE_C_COMPILER_ID MATCHES "GNU")
    set(CMAKE_C_FLAGS       "${CMAKE_C_FLAGS} -fno-working-directory")
  endif()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdebug-prefix-map=${CMAKE_SOURCE_DIR}=.")
endif()

message( STATUS "Building OpenCoarrays version: ${full_git_describe}" )
set(OpenCoarrays_dist_string "OpenCoarrays-${full_git_describe}")
message( STATUS "Building for target architecture: ${CMAKE_SYSTEM_PROCESSOR}" )

#Print an error message on an attempt to build inside the source directory tree:
if ("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
  message(FATAL_ERROR "ERROR! "
    "CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
    " == CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}"
    "\nThis archive does not support in-source builds:\n"
    "You must now delete the CMakeCache.txt file and the CMakeFiles/ directory under "
    "the 'src' source directory or you will not be able to configure correctly!"
    "\nYou must now run something like:\n"
    "  $ rm -r CMakeCache.txt CMakeFiles/"
    "\n"
    "Please create a directory outside the opencoarrays source tree and build under that outside directory "
    "in a manner such as\n"
    "  $ mkdir build-opencarrays\n"
    "  $ cd build-opencoarrays\n"
    "  $ CC=gcc FC=gfortran cmake -DBUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/install/dir /path/to/opencoarrays/src/dir \n"
    "\nsubstituting the appropriate syntax for your shell (the above line assumes the bash shell)."
    )
endif()

# Silence warnings about dereferencing unset variables
if(NOT CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_FLAGS "")
endif()
if(NOT CMAKE_REQUIRED_LIBRARIES)
  set(CMAKE_REQUIRED_LIBRARIES "")
endif()
if(NOT CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_INCLUDES "")
endif()

#Report untested Fortran compiler unless explicitly directed to build all examples.
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" )
  set(gfortran_compiler true)
  set ( CMAKE_C_FLAGS_CODECOVERAGE "-fprofile-arcs -ftest-coverage -O0"
    CACHE STRING "Code coverage C compiler flags")
  set ( CMAKE_Fortran_FLAGS_CODECOVERAGE "-fprofile-arcs -ftest-coverage -O0"
    CACHE STRING "Code coverage Fortran compiler flags")
else()
  message(WARNING
    "\n"
    "Attempting to build with untested Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}. "
    "Please report any failures to opencoarrays@googlegroups.com\n\n"
  )
endif()

#-----------------------------------------------------------------
# Set CMAKE_Fortran_COMPILER_VERSION if CMake doesn't do it for us
#-----------------------------------------------------------------
if (  NOT CMAKE_Fortran_COMPILER_VERSION )
  if ( NOT (CMAKE_VERSION VERSION_LESS 3.3.1) )
    message( AUTHOR_WARNING
     "CMake ${CMAKE_VERSION} should know about Fortran compiler versions but is missing CMAKE_Fortran_COMPILER_VERSION variable."
    )
  endif()
  # No CMAKE_Fortran_COMPILER_VERSION set, build our own
  # Try extracting it directly from ISO_FORTRAN_ENV's compiler_version
  # Write program for introspection
  file( WRITE "${CMAKE_BINARY_DIR}/get_compiler_ver.f90"
    "program main
     use iso_fortran_env, only: compiler_version, output_unit
     write(output_unit,'(a)') compiler_version()
end program"
  )
  try_run( PROG_RAN COMPILE_SUCCESS
    "${CMAKE_BINARY_DIR}" "${CMAKE_BINARY_DIR}/get_compiler_ver.f90"
    RUN_OUTPUT_VARIABLE VER_STRING
  )
  if ( COMPILE_SUCCESS )
    string( REGEX MATCH "[0-9]+\\.[0-9]+(\\.[0-9]+)?"
      DETECTED_VER "${VER_STRING}"
    )
    message( STATUS "Detected Fortran compiler as ${VER_STRING}" )
    message( STATUS "Extracted version number: ${DETECTED_VER}" )
  endif()
  if( ( NOT COMPILE_SUCCESS ) OR ( NOT DETECTED_VER ) )
    message( WARNING "Could not reliably detect Fortran compiler version. We'll infer it from
the C compiler if it matches the Fortran compiler ID." )
  endif()
  if( "${CMAKE_C_COMPILER_ID}" MATCHES "${CMAKE_Fortran_COMPILER_ID}" )
    set( DETECTED_VER "${CMAKE_C_COMPILER_VERSION}" )
  else()
    message( FATAL_ERROR "Exhausted all possible means of detecting the Fortran compiler version, cannot proceed!" )
  endif()
  set( CMAKE_Fortran_COMPILER_VERSION "${DETECTED_VER}" )
endif()

if(CMAKE_BUILD_TYPE MATCHES "Debug|DEBUG|debug")
  add_definitions(-DEXTRA_DEBUG_OUTPUT)
endif()

  # We have populated CMAKE_Fortran_COMPILER_VERSION if it was missing
  if(gfortran_compiler AND (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 5.0.0))
    set(opencoarrays_aware_compiler true)
    add_definitions(-DPREFIX_NAME=_gfortran_caf_)
  else()
    set(opencoarrays_aware_compiler false)
    add_definitions(-DPREFIX_NAME=_caf_extensions_)
  endif()
  if(gfortran_compiler AND (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 5.4.0))
    # GCC patch to fix issue accepted for 5.4 release
    # See https://github.com/sourceryinstitute/opencoarrays/issues/28 and
    # https://groups.google.com/forum/#!msg/opencoarrays/RZOwwYTqG80/46S9eL696dgJ
    message( STATUS "Disabling optimization flags due to GCC < 5.4 bug")
    set(CMAKE_Fortran_FLAGS_RELEASE -O0
      CACHE STRING "Flags used by the compiler during release builds." FORCE)
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -DNDEBUG -O0"
      CACHE STRING "Flags used by the compiler during release builds with debug info" FORCE)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O0")
  endif()
  if ( gfortran_compiler AND ( NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.0.0 ) )
    add_definitions(-DGCC_GE_7) # Tell library to build against GFortran 7.x bindings b/c we might be using clang for C
  endif()
  if ( gfortran_compiler AND ( NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 8.0.0 ) )
    add_definitions(-DGCC_GE_8) # Tell library to build against GFortran 8.x bindings w/ descriptor change
  endif()

if(gfortran_compiler)
  set(OLD_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
  set(CMAKE_REQUIRED_FLAGS $<$<COMPILE_LANGUAGE:Fortran>:-fcoarray=single -ffree-form>)
endif()
include(CheckFortranSourceCompiles)
CHECK_Fortran_SOURCE_COMPILES("
  program main
    implicit none
    integer :: i
    i = this_image()
  end program
" Check_Simple_Coarray_Fortran_Source_Compiles)
if(gfortran_compiler)
  set (CMAKE_REQUIRED_FLAGS ${OLD_REQUIRED_FLAGS})
  unset(OLD_REQUIRED_FLAGS)
endif()


#----------------------------------------------------------------------------
# Find MPI and set some flags so that FC and CC can point to gfortran and gcc
#----------------------------------------------------------------------------

# If the user passes FC=mpifort etc. check and prefer that location
get_filename_component( FTN_COMPILER_NAME "${CMAKE_Fortran_COMPILER}"
  NAME )
get_filename_component( C_COMPILER_NAME "${CMAKE_C_COMPILER}"
  NAME )
get_filename_component( FTN_COMPILER_DIR "${CMAKE_Fortran_COMPILER}"
  REALPATH )
get_filename_component( C_COMPILER_DIR "${CMAKE_C_COMPILER}"
  REALPATH )

if (FTN_COMPILER_NAME MATCHES "^[mM][pP][iI]")
  message( DEPRECATION "Setting the Fortran compiler to an MPI wrapper script is deprecated!")
  set (MPI_Fortran_COMPILER "${CMAKE_Fortran_COMPILER}")
endif()
if (C_COMPILER_NAME MATCHES "^[mM][pP][iI]")
  message( DEPRECATION "Setting the Fortran compiler to an MPI wrapper script is deprecated!")
  set (MPI_C_COMPILER "${CMAKE_C_COMPILER}")
endif()

find_package( MPI )

if ( (NOT MPI_C_FOUND) OR (NOT MPI_Fortran_FOUND) OR (NOT MPIEXEC))
  # Get default install location of MPICH from install.sh
  message(WARNING "Could not find all MPI components!")
  message(WARNING "
MPI_C_FOUND        = ${MPI_C_FOUND}
MPI_Fortran_FOUND  = ${MPI_Fortran_FOUND}
MPIEXEC = ${MPIEXEC}
")
  execute_process( COMMAND "./install.sh" -P mpich
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    OUTPUT_VARIABLE DEFAULT_MPICH_INSTALL_LOC
    OUTPUT_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  find_program (MY_MPI_EXEC NAMES mpirun mpiexec lamexec srun
    PATHS "${DEFAULT_MPICH_INSTALL_LOC}" ENV PATH
    HINTS "${FTN_COMPILER_DIR}" "${C_COMPILER_DIR}"
    PATH_SUFFIXES bin)
  set ( MPI_HOME "${MPI_HOME}" "${MY_MPI_EXEC}" "${MY_MPI_EXEC}/.." )
  find_package( MPI REQUIRED )
endif()
list(REMOVE_DUPLICATES MPI_Fortran_INCLUDE_PATH)

# Test for consistent MPI environment
if (NOT MPIEXEC)
  message ( ERROR "CMake failed to find `mpiexec` or similar. If building with `./install.sh` please
report this bug to the OpenCoarrays developers at
https://github.com/sourceryinstitute/opencoarrays/issues, otherwise point CMake
to the desired MPI runtime.")
else()
  add_definitions(-DHAVE_MPI)
endif()

get_filename_component(MPIEXEC_RELATIVE_LOC "${MPIEXEC}"
  PROGRAM)
get_filename_component(MPIEXEC_ABS_LOC "${MPIEXEC_RELATIVE_LOC}"
  REALPATH)
get_filename_component(MPIEXEC_DIR "${MPIEXEC_ABS_LOC}"
  DIRECTORY)

get_filename_component(MPICC_RELATIVE_LOC "${MPI_C_COMPILER}"
  PROGRAM)
get_filename_component(MPICC_ABS_LOC "${MPICC_RELATIVE_LOC}"
  REALPATH)
get_filename_component(MPICC_DIR "${MPICC_ABS_LOC}"
  DIRECTORY)

get_filename_component(MPIFC_RELATIVE_LOC "${MPI_Fortran_COMPILER}"
  PROGRAM)
get_filename_component(MPIFC_ABS_LOC "${MPIFC_RELATIVE_LOC}"
  REALPATH)
get_filename_component(MPIFC_DIR "${MPIFC_ABS_LOC}"
  DIRECTORY)

if ((MPIEXEC_DIR STREQUAL MPICC_DIR) AND (MPIEXEC_DIR STREQUAL MPIFC_DIR))
  message ( STATUS "MPI runtime and compile time environments appear to be consistent")
else()
  message ( WARNING "MPIEXEC is in \"${MPIEXEC_DIR},\"
which differs from the location of MPICC and/or MPIFC which are in
\"${MPICC_DIR}\" and \"${MPIFC_DIR},\" respectively.
This is likely indicative of a problem. If building with `./install.sh` please report
this to the OpenCoarrays developers by filing a new issue at:
https://github.com/sourceryinstitute/OpenCoarrays/issues/new")
endif()

#-----------------------------------------------
# Work around bug #317 present on fedora systems
#-----------------------------------------------
if( (MPI_C_LINK_FLAGS MATCHES "noexecstack") OR (MPI_Fortran_LINK_FLAGS MATCHES "noexecstack") )
  message ( WARNING
"The `noexecstack` linker flag was found in the MPI_<lang>_LINK_FLAGS variable. This is
known to cause segmentation faults for some Fortran codes. See, e.g.,
https://gcc.gnu.org/bugzilla/show_bug.cgi?id=71729 or
https://github.com/sourceryinstitute/OpenCoarrays/issues/317.

`noexecstack` is being replaced with `execstack`"
    )
  string(REPLACE "noexecstack"
    "execstack" MPI_C_LINK_FLAGS_FIXED ${MPI_C_LINK_FLAGS})
  string(REPLACE "noexecstack"
    "execstack" MPI_Fortran_LINK_FLAGS_FIXED ${MPI_Fortran_LINK_FLAGS})
  set(MPI_C_LINK_FLAGS "${MPI_C_LINK_FLAGS_FIXED}" CACHE STRING
    "MPI C linking flags" FORCE)
  set(MPI_Fortran_LINK_FLAGS "${MPI_Fortran_LINK_FLAGS_FIXED}" CACHE STRING
    "MPI Fortran linking flags" FORCE)
endif()

#--------------------------------------------------------
# Make sure a simple "hello world" C mpi program compiles
#--------------------------------------------------------
set(OLD_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
set(CMAKE_REQUIRED_FLAGS ${MPI_C_COMPILE_FLAGS} ${MPI_C_LINK_FLAGS})
set(OLD_INCLUDES ${CMAKE_REQUIRED_INCLUDES})
set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
set(OLD_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
include (CheckCSourceCompiles)
CHECK_C_SOURCE_COMPILES("
#include <mpi.h>
#include <stdio.h>
int main(int argc, char** argv) {
  MPI_Init(NULL, NULL);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  printf('Hello world from processor %s, rank %d out of %d processors',
         processor_name, world_rank, world_size);
  MPI_Finalize();
}"
MPI_C_COMPILES)
set(CMAKE_REQUIRED_FLAGS ${OLD_REQUIRED_FLAGS})
set(CMAKE_REQUIRED_INCLUDES ${OLD_INCLUDES})
set(CMAKE_REQUIRED_LIBRARIES ${OLD_LIBRARIES})
unset(OLD_REQUIRED_FLAGS)
unset(OLD_INCLUDES)
unset(OLD_LIBRARIES)

if (NOT MPI_C_COMPILES)
  message(FATAL_ERROR "MPI_C is missing! "
    "Try setting MPI_C_COMPILER to the appropriate C compiler wrapper script and reconfigure. "
    "i.e., `cmake -DMPI_C_COMPILER=/path/to/mpicc ..` or set it by editing the cache using "
    "cmake-gui or ccmake."
    )
endif()

#--------------------------------------------------------------
# Make sure a simple "hello world" Fortran mpi program compiles
# Try using mpi.mod first then fall back on includ 'mpif.h'
#--------------------------------------------------------------
set(OLD_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
set(CMAKE_REQUIRED_FLAGS "-ffree-form" ${MPI_Fortran_COMPILE_FLAGS} ${MPI_Fortran_LINK_FLAGS})
set(OLD_INCLUDES ${CMAKE_REQUIRED_INCLUDES})
set(CMAKE_REQUIRED_INCLUDES ${MPI_Fortran_INCLUDE_PATH})
set(OLD_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_Fortran_LIBRARIES})
include (CheckFortranSourceCompiles)
CHECK_Fortran_SOURCE_COMPILES("
program mpi_hello
use mpi
implicit none
integer :: ierr, mpi_world_size, mpi_world_rank, res_len
character*(MPI_MAX_PROCESSOR_NAME) :: proc
call mpi_init(ierr)
call mpi_comm_size(MPI_COMM_WORLD,mpi_world_size,ierr)
call mpi_comm_rank(MPI_COMM_WORLD,mpi_world_rank,ierr)
call mpi_get_processor_name(proc,res_len,ierr)
write(*,*) 'Hello from processor ', trim(proc), ' rank ', mpi_world_rank, ' out of ', mpi_world_size, '.'
call mpi_finalize(ierr)
end program
"
MPI_Fortran_MODULE_COMPILES)
set(CMAKE_REQUIRED_FLAGS ${OLD_REQUIRED_FLAGS})
set(CMAKE_REQUIRED_INCLUDES ${OLD_INCLUDES})
set(CMAKE_REQUIRED_LIBRARIES ${OLD_LIBRARIES})
unset(OLD_REQUIRED_FLAGS)
unset(OLD_INCLUDES)
unset(OLD_LIBRARIES)

#--------------------------------
# If that failed try using mpif.h
#--------------------------------
set(OLD_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
set(CMAKE_REQUIRED_FLAGS "-ffree-form" ${MPI_Fortran_COMPILE_FLAGS} ${MPI_Fortran_LINK_FLAGS})
set(OLD_INCLUDES ${CMAKE_REQUIRED_INCLUDES})
set(CMAKE_REQUIRED_INCLUDES ${MPI_Fortran_INCLUDE_PATH})
set(OLD_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_Fortran_LIBRARIES})
include (CheckFortranSourceCompiles)
CHECK_Fortran_SOURCE_COMPILES("
program mpi_hello
implicit none
include 'mpif.h'
integer :: ierr, mpi_world_size, mpi_world_rank, res_len
character*(MPI_MAX_PROCESSOR_NAME) :: proc
call mpi_init(ierr)
call mpi_comm_size(MPI_COMM_WORLD,mpi_world_size,ierr)
call mpi_comm_rank(MPI_COMM_WORLD,mpi_world_rank,ierr)
call mpi_get_processor_name(proc,res_len,ierr)
write(*,*) 'Hello from processor ', trim(proc), ' rank ', mpi_world_rank, ' out of ', mpi_world_size, '.'
call mpi_finalize(ierr)
end program
"
  MPI_Fortran_INCLUDE_COMPILES)
set(CMAKE_REQUIRED_FLAGS ${OLD_REQUIRED_FLAGS})
set(CMAKE_REQUIRED_INCLUDES ${OLD_INCLUDES})
set(CMAKE_REQUIRED_LIBRARIES ${OLD_LIBRARIES})
unset(OLD_REQUIRED_FLAGS)
unset(OLD_INCLUDES)
unset(OLD_LIBRARIES)

if ( (NOT MPI_Fortran_MODULE_COMPILES) AND (NOT MPI_Fortran_INCLUDE_COMPILES) )
  message ( WARNING "It appears that the Fortran MPI compiler is not working. "
    "For OpenCoarrays Aware compilers, this may be irrelavent: "
    "  The src/extensions/opencoarrays.F90 module will be disabled, but it is "
    "  possible that the build will succeed, despite this fishy circumstance."
    )
endif()

if ( MPI_Fortran_MODULE_COMPILES )
  add_definitions(-DMPI_WORKING_MODULE)
else()
  message ( WARNING "It appears that MPI was built with a different Fortran compiler. "
    "It is possible that this may cause unpredictable behavior. The build will continue "
    "using `mpif.h` BUT please report any suspicious behavior to the OpenCoarrays "
    "developers."
    )
endif()

#----------------
# Setup MPI flags
#----------------
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

#---------------------------------------------------
# Use standardized GNU install directory conventions
#---------------------------------------------------
include(GNUInstallDirs)
set(mod_dir_tail "${OpenCoarrays_dist_string}_${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/${OpenCoarrays_dist_string}-tests")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

#-----------------
# Install manpages
#-----------------
install(FILES "${CMAKE_SOURCE_DIR}/doc/man/man1/caf.1" "${CMAKE_SOURCE_DIR}/doc/man/man1/cafrun.1"
  DESTINATION "${CMAKE_INSTALL_MANDIR}/man1"
  PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
  COMPONENT documentation)

#---------------------------------------------------
# Define macro for compiling with caf wrapper script
#---------------------------------------------------
if(APPLE)
  find_program(CODESIGN codesign)
endif()

function(caf_compile_executable target main_depend)
  set(includes "")
  foreach(includedir ${MPI_Fortran_INCLUDE_PATH})
    list(APPEND includes "-I${includedir}")
  endforeach()
  string(TOUPPER "${CMAKE_BUILD_TYPE}" build_type)
  separate_arguments(config_Fortran_flags UNIX_COMMAND "${CMAKE_Fortran_FLAGS_${build_type}}")
  get_directory_property( DirDefs DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" COMPILE_DEFINITIONS )
  set(localDefs "")
  foreach(d ${DirDefs})
    list(APPEND localDefs "-D${d}")
  endforeach()
    add_custom_command(OUTPUT "${target}"
      COMMAND "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/caf"
                 ${includes} ${localDefs} ${config_Fortran_flags}
                 -o "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${target}"
                 "${CMAKE_CURRENT_SOURCE_DIR}/${main_depend}"
		 "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libopencoarrays_test_utilities.a"
		 "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libopencoarrays_mod.a"
		 ${ARGN}
      DEPENDS "${main_depend}" ${ARGN} caf_mpi_static opencoarrays_test_utilities opencoarrays_mod
      VERBATIM
      )
  add_custom_target("build_${target}" ALL
    DEPENDS "${target}")
endfunction(caf_compile_executable)

enable_testing()

#--------------------------------------------------------
# Setup shellcheck if present for testing/linting scripts
#--------------------------------------------------------
find_program(SHELLCHECK_EXE shellcheck
  DOC "Path to shellcheck executable for linting scripts"
  )
if (CAF_RUN_DEVELOPER_TESTS OR $ENV{OPENCOARRAYS_DEVELOPER})
  if(NOT SHELLCHECK_EXE)
    message( AUTHOR_WARNING "OpenCoarrays developers should install shellcheck to test/lint all shell scripts.
    See https://github.com/koalaman/shellcheck#installing for info on obtaining shellcheck.")
  endif()
endif()

function(lint_script script_dir script_name)
  if (SHELLCHECK_EXE)
    add_test(NAME "shellcheck:${script_name}"
      COMMAND ${SHELLCHECK_EXE} -x "${script_dir}/${script_name}"
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
  elseif (CAF_RUN_DEVELOPER_TESTS OR $ENV{OPENCOARRAYS_DEVELOPER})
    message( AUTHOR_WARNING "test: shellcheck:${script_name} not run because shellcheck not installed." )
  endif()
endfunction()

#-----------------------------------------------
# Setup script style testing & enforcement macro
#-----------------------------------------------

find_program(style_pl style.pl "${CMAKE_SOURCE_DIR}/developer-scripts/")
function(check_script_style script_full_path)
  if(style_pl)
    add_test(NAME "style:${script_full_path}"
      COMMAND "${style_pl}" "${script_full_path}")
  endif()
endfunction()

#------------------------------------------------------------------------------
# Add custom properties on targets for controling number of images during tests
#------------------------------------------------------------------------------
define_property(TARGET
  PROPERTY MIN_IMAGES
  BRIEF_DOCS "Minimum allowable images for the test <integer>"
  FULL_DOCS "Property to mark executable targets run as tests that they require at least <MIN_IMAGES> images to run"
  )

define_property(TARGET
  PROPERTY POWER_2_IMGS
  BRIEF_DOCS "True if test must be run with a power of 2 images (T/F)"
  FULL_DOCS "Property to mark executable targets run as tests that they require 2^n images."
  )


#-------------------------------
# Recurse into the src directory
#-------------------------------
add_subdirectory(src)

#-----------------------------------------------------
# Publicize installed location to other CMake projects
#-----------------------------------------------------
install(EXPORT OpenCoarraysTargets
  NAMESPACE
    OpenCoarrays::
  DESTINATION
    "${CMAKE_INSTALL_LIBDIR}/cmake/opencoarrays"
)
include(CMakePackageConfigHelpers) # standard CMake module
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/OpenCoarraysConfigVersion.cmake"
  VERSION "${opencoarrays_VERSION}"
  COMPATIBILITY AnyNewerVersion
)
configure_file("${CMAKE_SOURCE_DIR}/cmake/pkg/OpenCoarraysConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/OpenCoarraysConfig.cmake" @ONLY)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/OpenCoarraysConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/OpenCoarraysConfigVersion.cmake"
  DESTINATION
    "${CMAKE_INSTALL_LIBDIR}/cmake/opencoarrays"
)

#------------------------------------------
# Add portable unistall command to makefile
#------------------------------------------
# Adapted from the CMake Wiki FAQ
configure_file ( "${CMAKE_SOURCE_DIR}/cmake/uninstall.cmake.in" "${CMAKE_BINARY_DIR}/uninstall.cmake"
    @ONLY)

add_custom_target ( uninstall
  COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/uninstall.cmake" )

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure)
# See JSON-Fortran's CMakeLists.txt file to find out how to get the check target to depend
# on the test executables

#------------------------------------
# Add command to hash installed files
#------------------------------------
configure_file( "${CMAKE_SOURCE_DIR}/cmake/hash-installed.cmake.in" "${CMAKE_BINARY_DIR}/hash-installed.cmake"
  @ONLY)
add_custom_target( hash_installed
  COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/hash-installed.cmake")

#---------------------------------------------------------------------------------------
# Define macro for adding CAF tests, and ensuring proper flags are passed to MPI runtime
#---------------------------------------------------------------------------------------
# Test for OpenMPI moved to src/mpi/CMakeLists.txt and passed back up through set(... PARENT_SCOPE)
function(add_caf_test name num_caf_img test_target)
  # Function to add MPI tests.
  if(TARGET ${test_target})
    get_target_property(min_test_imgs ${test_target} MIN_IMAGES)
  elseif(TARGET build_${test_target})
    get_target_property(min_test_imgs build_${test_target} MIN_IMAGES)
  endif()
  if(min_test_imgs)
    if(num_caf_img LESS min_test_imgs)
      message( FATAL_ERROR "Test ${name} requires ${min_test_imgs} but was only given ${num_caf_images}" )
    endif()
  endif()
  # Add a host file for OMPI
  if ( openmpi )
    set(test_parameters --hostfile ${CMAKE_BINARY_DIR}/hostfile)
  endif()
  if ( ((N_CPU LESS num_caf_img) OR (N_CPU EQUAL 0)) )
    message(STATUS "Test ${name} is oversubscribed: ${num_caf_img} CAF images requested with ${N_CPU} system processor available.")
    if ( openmpi )
      if (min_test_imgs)
	set( num_caf_img ${min_test_imgs} )
      elseif ( N_CPU LESS 2 )
        set( num_caf_img 2 )
      endif()
      set (test_parameters ${test_parameters} --oversubscribe)
      message( STATUS "Open-MPI back end detected, passing --oversubscribe for oversubscribed test, ${name}, with ${num_caf_img} ranks/images." )
    endif()
  endif()
  if( openmpi)
    if(NOT DEFINED ENV{TRAVIS})
      set (test_parameters ${test_parameters} --allow-run-as-root)
    endif()
  endif()
  set(test_parameters -np ${num_caf_img} ${test_parameters})
  if(DEFINED ARGN)
    add_test(NAME ${name} COMMAND "bash" "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/cafrun" ${test_parameters} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test_target}" ${ARGN})
  else()
    add_test(NAME ${name} COMMAND "bash" "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/cafrun" ${test_parameters} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test_target}")
  endif()
  set_property(TEST ${name} PROPERTY PASS_REGULAR_EXPRESSION "Test passed.")
endfunction(add_caf_test)

#--------------
# Add OCA tests
#--------------

if(opencoarrays_aware_compiler)
  if (CAF_RUN_DEVELOPER_TESTS OR $ENV{OPENCOARRAYS_DEVELOPER})
    message ( STATUS "Running Developer tests is enabled. Some tests may fail for open issues." )
  endif()
  # Unit tests targeting each libcaf_mpi function, argument, and branch of code
  add_caf_test(initialize_mpi 2 initialize_mpi)
  add_caf_test(register 2 register)
  add_caf_test(register_vector 2 register_vector)
  add_caf_test(register_alloc_vector 2 register_alloc_vector)
  add_caf_test(allocate_as_barrier 2 allocate_as_barrier)
  if(gfortran_compiler AND (NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.0.0) OR (CAF_RUN_DEVELOPER_TESTS OR $ENV{OPENCOARRAYS_DEVELOPER}))
    if( CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.0.0 )
      message( AUTHOR_WARNING "Developer tests requested and GFortran < 7: test allocate_as_barrier_proc may fail" )
    endif()
    add_caf_test(allocate_as_barrier_proc 8 allocate_as_barrier_proc)
  endif()
  if (gfortran_compiler AND (NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.0.0) OR (CAF_RUN_DEVELOPER_TESTS OR $ENV{OPENCOARRAYS_DEVELOPER}))
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.0.0)
      message( AUTHOR_WARNING "Allocatable components of coarray derived type developer tests enabled, despite lack of support in GFortran < 7\n These tests should fail." )
    endif()
    add_caf_test(register_alloc_comp_1 2 register_alloc_comp_1)
    add_caf_test(register_alloc_comp_2 2 register_alloc_comp_2)
    add_caf_test(register_alloc_comp_3 2 register_alloc_comp_3)
    add_caf_test(async_comp_alloc 6 async_comp_alloc)
    add_caf_test(async_comp_alloc_2 2 async_comp_alloc_2)
    add_caf_test(comp_allocated_1 2 comp_allocated_1)
    add_caf_test(comp_allocated_2 2 comp_allocated_2)
    add_caf_test(alloc_comp_get_convert_nums 2 alloc_comp_get_convert_nums)
    if(NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 8)
      add_caf_test(team_number 8 team_number)
      add_caf_test(teams_subset 3 teams_subset)
      add_caf_test(get_communicator 3  get_communicator)
      add_caf_test(teams_coarray_get 5 teams_coarray_get)
      add_caf_test(teams_coarray_get_by_ref 5 teams_coarray_get_by_ref)
      add_caf_test(teams_coarray_send 5 teams_coarray_send)
      add_caf_test(teams_coarray_send_by_ref 5 teams_coarray_send_by_ref)
      add_caf_test(teams_coarray_sendget 5 teams_coarray_sendget)
      add_caf_test(sync_team 8 sync_team)
      add_caf_test(alloc_comp_multidim_shape 2 alloc_comp_multidim_shape)
    endif()
  endif()

  if (gfortran_compiler)
    if((NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.0.0) OR (CAF_RUN_DEVELOPER_TESTS OR $ENV{OPENCOARRAYS_DEVELOPER}))
      add_caf_test(send_convert_nums 2 send_convert_nums)
      add_caf_test(sendget_convert_nums 3 sendget_convert_nums)
      add_caf_test(sendget_convert_char_array 3 sendget_convert_char_array)
    else()
      message( AUTHOR_WARNING "Skipping the following tests due to GFortran < 7.0.0 lack of compatibility:
        sendget_convert_nums.f90
        sendget_convert_char_array.f90
        send_convert_nums.f90")
    endif()
    if((NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.3.0) OR (CAF_RUN_DEVELOPER_TESTS OR $ENV{OPENCOARRAYS_DEVELOPER}))
      add_caf_test(send_convert_char_array 2 send_convert_char_array)
      add_caf_test(alloc_comp_send_convert_nums 2 alloc_comp_send_convert_nums)
    elseif((NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.4.0) OR (CAF_RUN_DEVELOPER_TESTS OR $ENV{OPENCOARRAYS_DEVELOPER}))
      add_caf_test(send-strided-self 2 send-strided-self)
    endif()
    if((CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.3.0) AND (CAF_RUN_DEVELOPER_TESTS OR $ENV{OPENCOARRAYS_DEVELOPER}))
      message( AUTHOR_WARNING "Skipping the following tests to GFortran < 7.3.0 lack of compatibility:
        send_convert_char_array.f90
        alloc_comp_send_convert_nums.f90")
    endif()
    if((CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.4.0) AND (CAF_RUN_DEVELOPER_TESTS OR $ENV{OPENCOARRAYS_DEVELOPER}))
      message( AUTHOR_WARNING "Skipping the following test to GFortran < 7.4.0 lack of compatibility:
        send-strided-self.f90")
    endif()
  endif()

  # Pure get tests
  add_caf_test(get_array 2 get_array)
  add_caf_test(get_self 2 get_self)
  add_caf_test(get_convert_nums 2 get_convert_nums)
  add_caf_test(get_convert_char_array 2 get_convert_char_array)
  add_caf_test(get_with_offset_1d 2 get_with_offset_1d)
  add_caf_test(whole_get_array 2 whole_get_array)
  add_caf_test(strided_get 2 strided_get)

  # Pure send tests
  add_caf_test(send_array 2 send_array)
  add_caf_test(convert-before-put 3 convert-before-put)
  add_caf_test(send_with_vector_index 2 send_with_vector_index)

  # Pure sendget tests
  add_caf_test(strided_sendget 3 strided_sendget)
  add_caf_test(get_with_vector_index 4 get_with_vector_index)

  # Collective subroutine tests
  add_caf_test(co_sum 4 co_sum_test)
  add_caf_test(co_broadcast 4 co_broadcast_test)
  add_caf_test(co_broadcast_derived_type 4 co_broadcast_derived_type_test)
  if((gfortran_compiler AND (NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10.0.0)) OR (CAF_RUN_DEVELOPER_TESTS OR $ENV{OPENCOARRAYS_DEVELOPER}))
    add_caf_test(co_broadcast_allocatable_components 4 co_broadcast_allocatable_components_test)
  endif()
  add_caf_test(co_min 4 co_min_test)
  add_caf_test(co_max 4 co_max_test)
  add_caf_test(co_reduce 4 co_reduce_test)
  add_caf_test(co_reduce_res_im 4 co_reduce_res_im)
  add_caf_test(co_reduce_string 4 co_reduce_string)
  add_caf_test(syncimages_status 8 syncimages_status)
  add_caf_test(sync_ring_abort_np3 3 sync_image_ring_abort_on_stopped_image)
  add_caf_test(sync_ring_abort_np7 7 sync_image_ring_abort_on_stopped_image)
  add_caf_test(simpleatomics 8 atomics)

  # Synchronization tests
  add_caf_test(syncall 8 syncall)
  add_caf_test(syncimages 8 syncimages)
  add_caf_test(syncimages2 8 syncimages2)
  add_caf_test(duplicate_syncimages 8 duplicate_syncimages)

  # possible logic error in the following test
#  add_caf_test(increment_my_neighbor 32 increment_my_neighbor)

  # Integration tests verifying the use of libcaf_mpi in applications
  add_caf_test(hello_multiverse 2 hello_multiverse)
  add_caf_test(coarray_burgers_pde 2 coarray_burgers_pde)
  add_caf_test(co_heat 2 co_heat)
  add_caf_test(asynchronous_hello_world 3 asynchronous_hello_world)

  # Regression tests based on reported issues
  if((gfortran_compiler AND (NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.0.0)) OR (CAF_RUN_DEVELOPER_TESTS OR $ENV{OPENCOARRAYS_DEVELOPER}))
    if( CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.0.0 )
      message( AUTHOR_WARNING "Developer tests requested and GFortran < 7: test source-alloc-no-sync may fail" )
    endif()
    # GFortran PR 78505 only fixed on trunk/gcc 7, issue #243
    add_caf_test(source-alloc-no-sync 8 source-alloc-sync)
  endif()

  # Fixed on GCC 6 > 6.4, 7 > 7.3 and 8.x
  if (FALSE) # FIXME!
    add_caf_test(non-coarray-vector-indexed-lhs issue-322-non-coarray-vector-idx-lhs)
  endif()

  # Open GCC 7 regressions
  if ((CAF_RUN_DEVELOPER_TESTS OR $ENV{OPENCOARRAYS_DEVELOPER}) OR
      (gfortran_compiler AND (
	  (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.0.0) OR
	  (NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.3.0)
	  )))
    add_caf_test(put-allocatable-coarray-comp 2 issue-422-send)
    add_caf_test(get-put-allocatable-comp 3 issue-422-send-get)
  else()
    message( AUTHOR_WARNING
      "Skipping regressions in GFortran 7.0:
       put-allocatable-coarray-comp (issue-422)
       get-put-allocatable-comp     (issue-422)")
  endif()

  add_caf_test(allocatable_p2p_event_post 4 allocatable_p2p_event_post)
  # Fixed GCC 7 regressions, should run on GCC 6 and 7
  add_caf_test(static_event_post_issue_293 3 static_event_post_issue_293)


  # These co_reduce (#172, fixed by PR #332, addl discussion in PR
  # #331) tests are for bugs not regressions. Should be fixed in all
  # version of GCC, I beleive
  add_caf_test(co_reduce-factorial 4 co_reduce-factorial)
  add_caf_test(co_reduce-factorial-int8 4 co_reduce-factorial-int8)
  add_caf_test(co_reduce-factorial-int64 4 co_reduce-factorial-int64)

  # issues reported by @neok-m4700
  add_caf_test(issue-493-coindex-slice            8 issue-493-coindex-slice)
  add_caf_test(issue-488-multi-dim-cobounds-true  8 issue-488-multi-dim-cobounds true)
  add_caf_test(issue-488-multi-dim-cobounds-false 8 issue-488-multi-dim-cobounds false)
  add_caf_test(issue-503-multidim-array-broadcast 2 issue-503-multidim-array-broadcast)
  add_caf_test(issue-503-non-contig-red-ndarray   2 issue-503-non-contig-red-ndarray)
  add_caf_test(issue-552-send_by_ref-singleton    2 issue-552-send_by_ref-singleton)
  if((gfortran_compiler AND (NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.0.0)) OR (CAF_RUN_DEVELOPER_TESTS OR $ENV{OPENCOARRAYS_DEVELOPER}))
    add_caf_test(issue-515-mimic-mpi-gatherv 2 issue-515-mimic-mpi-gatherv)
  endif()

  # IMAGE FAIL tests
  if(NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.0.0)
    if(CAF_ENABLE_FAILED_IMAGES)
      add_caf_test(image_status_test_1 4 image_status_test_1)
      if ((NOT DEFINED ENV{TRAVIS}) OR (CAF_RUN_DEVELOPER_TESTS OR $ENV{OPENCOARRAYS_DEVELOPER}))
	add_caf_test(image_fail_test_1 4 image_fail_test_1)
	set_property(TEST image_fail_test_1 PROPERTY FAIL_REGULAR_EXPRESSION "Test failed.")
	# No other way to check that image_fail_test_1 passes.
	add_caf_test(image_fail_and_sync_test_1 4 image_fail_and_sync_test_1)
        add_caf_test(image_fail_and_sync_test_2 4 image_fail_and_sync_test_2)

	add_caf_test(image_fail_and_sync_test_3 4 image_fail_and_sync_test_3)
	add_caf_test(image_fail_and_status_test_1 4 image_fail_and_status_test_1)
	add_caf_test(image_fail_and_failed_images_test_1 4 image_fail_and_failed_images_test_1)
	add_caf_test(image_fail_and_stopped_images_test_1 4 image_fail_and_stopped_images_test_1)
	add_caf_test(image_fail_and_get_test_1 4 image_fail_and_get_test_1)
      endif()
    endif()
  endif()
else()
  add_test(co_sum_extension test-co_sum-extension.sh)
  set_property(TEST co_sum_extension PROPERTY PASS_REGULAR_EXPRESSION "Test passed.")

  add_test(co_broadcast_extension test-co_broadcast-extension.sh)
  set_property(TEST co_broadcast_extension PROPERTY PASS_REGULAR_EXPRESSION "Test passed.")

  add_test(co_min_extension test-co_min-extension.sh)
  set_property(TEST co_min_extension PROPERTY PASS_REGULAR_EXPRESSION "Test passed.")

  add_test(co_max_extension test-co_max-extension.sh)
  set_property(TEST co_max_extension PROPERTY PASS_REGULAR_EXPRESSION "Test passed.")

  add_test(co_reduce_extension test-co_reduce-extension.sh)
  set_property(TEST co_reduce_extension PROPERTY PASS_REGULAR_EXPRESSION "Test passed.")
endif()

# Test bash installation scripts
include(cmake/AddInstallationScriptTest.cmake )
add_installation_script_test(installation-scripts.sh src/tests/installation/)

if( CAF_ENABLE_FAILED_IMAGES AND ( NOT
      (	"${CMAKE_C_COMPILER_ID}" STREQUAL "GNU" ) AND ( "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64" ) )
    )
  message ( WARNING
    "ISO_Fortran_BINDING has only been tested with recent GCC C compilers on x86_64 architectures!
     Some people have reported problems on other architectures or using other C compilers (Clang).
     There are likely bugs present; use at your own risk!")
endif()

# Test ISO_Fortran_binding library
if( CAF_ENABLE_ISO_Fortran_BINDING OR ( CAF_RUN_DEVELOPER_TESTS OR $ENV{OPENCOARRAYS_DEVELOPER} ) )
  # See https://github.com/sourceryinstitute/OpenCoarrays/issues/523#issuecomment-401613209
  add_test(NAME ISO_Fortran_binding_tests COMMAND "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ISO_Fortran_binding_tests")
endif()

# Lint the Travis-CI scripts
set(TRAVIS_SCRIPTS
  install.linux.sh
  install.osx.sh
  test-script.InstallScript.sh
  test-script.cmake.sh)
foreach(SCRIPT ${TRAVIS_SCRIPTS})
  if(EXISTS "${CMAKE_SOURCE_DIR}/developer-scripts/travis/${SCRIPT}")
    lint_script("${CMAKE_SOURCE_DIR}/developer-scripts/travis" ${SCRIPT})
  endif()
endforeach()

# Create and install a manifest of installed files, and another containing their sha256 hashes
foreach(VAR CMAKE_BINARY_DIR CMAKE_INSTALL_DATADIR PROJECT_NAME SOURCE_DATE_EPOCH)
  install(CODE "set(${VAR} \"${${VAR}}\")")
endforeach()
install(SCRIPT "${CMAKE_CURRENT_LIST_DIR}/cmake/InstallFileManifest.cmake")
