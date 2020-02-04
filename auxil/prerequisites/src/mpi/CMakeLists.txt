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

if (NOT MPI_C_FOUND)
  find_package(MPI REQUIRED)
endif()

#-----------------------------------------------------------
# Add the actual library targets, then adjust them as needed
#-----------------------------------------------------------
add_library(opencoarrays_mod STATIC ../extensions/opencoarrays.F90)
set_target_properties(opencoarrays_mod
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${mod_dir_tail}"
  POSITION_INDEPENDENT_CODE TRUE)
target_compile_options(opencoarrays_mod
  PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${MPI_Fortran_COMPILE_FLAGS}>)
target_link_libraries(opencoarrays_mod
  PUBLIC ${MPI_Fortran_LINK_FLAGS}
  PUBLIC ${MPI_Fortran_LIBRARIES})
target_link_libraries(opencoarrays_mod
  PUBLIC caf_mpi_static)
target_include_directories(opencoarrays_mod PUBLIC
  $<$<COMPILE_LANGUAGE:Fortran>:${MPI_Fortran_INCLUDE_PATH}>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${mod_dir_tail}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${mod_dir_tail}>)

add_library(caf_mpi SHARED mpi_caf.c ../common/caf_auxiliary.c)
add_library(caf_mpi_static STATIC  mpi_caf.c ../common/caf_auxiliary.c)
target_link_libraries(caf_mpi
  PUBLIC ${MPI_C_LINK_FLAGS}
  PUBLIC ${MPI_C_LIBRARIES})
target_link_libraries(caf_mpi_static
  PUBLIC ${MPI_C_LINK_FLAGS}
  PUBLIC ${MPI_C_LIBRARIES})
set_target_properties(caf_mpi_static
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE)
target_include_directories(caf_mpi PUBLIC
  $<$<COMPILE_LANGUAGE:C>:${MPI_C_INCLUDE_PATH}>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_include_directories(caf_mpi_static PUBLIC
  $<$<COMPILE_LANGUAGE:C>:${MPI_C_INCLUDE_PATH}>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_options(caf_mpi
  PUBLIC $<$<COMPILE_LANGUAGE:C>:${MPI_C_COMPILE_FLAGS}>)
target_compile_options(caf_mpi_static
  PUBLIC $<$<COMPILE_LANGUAGE:C>:${MPI_C_COMPILE_FLAGS}>)

set(CAF_SO_VERSION 0)
if(gfortran_compiler)
  if(NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 8.0.0)
    set(CAF_SO_VERSION 3)
  elseif(NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.0.0)
    set(CAF_SO_VERSION 2)
  elseif(NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 6.0.0)
    set(CAF_SO_VERSION 1)
  endif()
endif()

set_target_properties ( caf_mpi
  PROPERTIES
  SOVERSION ${CAF_SO_VERSION}
#  VERSION ${PROJECT_VERSION}
)

# Create a symlink in the include dir
if(UNIX)
  add_custom_command(TARGET caf_mpi
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink "./${mod_dir_tail}/opencoarrays.mod" "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/opencoarrays.mod"
    COMMENT "Creating symlink ${CMAKE_INSTALL_INCLUDEDIR}/opencoarrays.mod --> ${CMAKE_INSTALL_INCLUDEDIR}/${mod_dir_tail}/opencoarrays.mod")
endif()

install(DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  FILES_MATCHING PATTERN "opencoarrays.mod")

set_target_properties( caf_mpi_static
  PROPERTIES
  SOVERSION ${CAF_SO_VERSION}
#  VERSION ${PROJECT_VERSION}
)

if (gfortran_compiler)
  target_compile_options(caf_mpi INTERFACE $<$<COMPILE_LANGUAGE:Fortran>:-fcoarray=lib>)
  target_compile_options(caf_mpi_static INTERFACE $<$<COMPILE_LANGUAGE:Fortran>:-fcoarray=lib>)
endif()

install(TARGETS opencoarrays_mod caf_mpi caf_mpi_static EXPORT OpenCoarraysTargets
  DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

#----------------------------------
# Determine if we're using Open MPI
#----------------------------------
cmake_host_system_information(RESULT N_CPU QUERY NUMBER_OF_LOGICAL_CORES)
set(N_CPU ${N_CPU} PARENT_SCOPE)
cmake_host_system_information(RESULT HOST_NAME QUERY HOSTNAME)
set(HOST_NAME ${HOST_NAME} PARENT_SCOPE)
execute_process(COMMAND ${MPIEXEC} --version
  OUTPUT_VARIABLE mpi_version_out)
if (mpi_version_out MATCHES "[Oo]pen[ -][Mm][Pp][Ii]")
  message( STATUS "OpenMPI detected")
  set ( openmpi true PARENT_SCOPE)
  # Write out a host file because OMPI's mpiexec is dumb
  file(WRITE ${CMAKE_BINARY_DIR}/hostfile "${HOST_NAME} slots=${N_CPU}\n")
  message( STATUS "hostfile written to: ${CMAKE_BINARY_DIR}/hostfile")
  if(NOT DEFINED ENV{TRAVIS})
    message( STATUS "Open-MPI back end detected, passing --allow-run-as-root to allow tests to pass when run with sudo or as root." )
  endif()
endif ()

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  set(gfortran_compiler true)
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Cray")
  set(cray_compiler true)
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "PGI")
  set(portland_group_compiler true)
endif()

if(gfortran_compiler AND (NOT opencoarrays_aware_compiler))
  target_compile_definitions(opencoarrays_mod
    PUBLIC -DCOMPILER_SUPPORTS_CAF_INTRINSICS)
endif()

option(CAF_EXPOSE_INIT_FINALIZE "Expose caf_init and caf_finalize in opencoarrays module" FALSE)
if(CAF_EXPOSE_INIT_FINALIZE)
  target_compile_definitions(opencoarrays_mod
    PRIVATE -DEXPOSE_INIT_FINALIZE)
endif()

include(CheckIncludeFile)
CHECK_INCLUDE_FILE("alloca.h" HAVE_ALLOCA)
if(NOT HAVE_ALLOCA)
  target_compile_definitions(caf_mpi
    PRIVATE -DALLOCA_MISSING)
  target_compile_definitions(caf_mpi_static
    PRIVATE -DALLOCA_MISSING)
  message(WARNING "Could not find <alloca.h>. Assuming functionality is provided elsewhere.")
endif()

#----------------------------------------------------------------------
# Test if MPI implementation provides features needed for failed images
#----------------------------------------------------------------------
set(MPI_HAS_FAULT_TOL_EXT YES)

CHECK_INCLUDE_FILE("signal.h" HAVE_SIGNAL_H)
if (NOT HAVE_SIGNAL_H)
  set(MPI_HAS_FAULT_TOL_EXT NO)
  message( FATAL_ERROR "Currently, OpenCoarrays cannot build without signal.h")
endif()

include(CheckSymbolExists)
CHECK_SYMBOL_EXISTS(SIGKILL "signal.h" HAVE_SIGKILL)
if(NOT HAVE_SIGKILL) # try -D_POSIX, needed for mingw-w64, maybe others, see #435
                     # https://github.com/sourceryinstitute/OpenCoarrays/issues/435#issuecomment-323592433
  list( APPEND CMAKE_REQUIRED_DEFINITIONS -D_POSIX)
  CHECK_SYMBOL_EXISTS(SIGKILL "signal.h" HAVE_SIGKILL2)
  if(HAVE_SIGKILL2)
    set(HAVE_SIGKILL ${HAVE_SIGKILL2})
    foreach(lib caf_mpi caf_mpi_static)
      target_compile_definitions(${lib}
    	PUBLIC -D_POSIX)
    endforeach()
  endif()
endif()

if (NOT HAVE_SIGKILL)
  set(MPI_HAS_FAULT_TOL_EXT NO)
  message (FATAL_ERROR "Currently, OpenCoarrays cannot build without SIGKILL from signal.h")
endif()

set(NEEDED_SYMBOLS MPIX_ERR_PROC_FAILED;MPIX_ERR_REVOKED;MPIX_Comm_failure_ack;MPIX_Comm_failure_get_acked;MPIX_Comm_shrink;MPIX_Comm_agree)

set(old_cmake_required_includes "${CMAKE_REQUIRED_INCLUDES}")
if(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${MPI_C_INCLUDE_PATH})
else()
  set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
endif()
set(old_cmake_required_flags "${CMAKE_REQUIRED_FLAGS}")
if(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS} ${MPI_C_COMPILE_FLAGS} ${MPI_C_LINK_FLAGS})
else()
  set(CMAKE_REQUIRED_FLAGS ${MPI_C_COMPILE_FLAGS} ${MPI_C_LINK_FLAGS})
endif()
set(old_cmake_required_libraries "${CMAKE_REQUIRED_LIBRARIES}")
if(CMAKE_REQUIRED_LIBRARIES)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${MPI_C_LIBRARIES})
else()
  set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
endif()

set(MPI_HEADERS mpi.h)
include(CheckIncludeFiles)
CHECK_INCLUDE_FILES("mpi.h;mpi-ext.h" HAVE_MPI_EXT)
if(HAVE_MPI_EXT)
    foreach(lib caf_mpi caf_mpi_static)
      target_compile_definitions(${lib}
    	PRIVATE -DHAVE_MPI_EXT_H)
    endforeach()
  set(MPI_HEADERS ${MPI_HEADERS};mpi-ext.h)
endif()

foreach(symbol ${NEEDED_SYMBOLS})
  CHECK_SYMBOL_EXISTS(${symbol} "${MPI_HEADERS}" HAVE_${symbol})
  if(NOT HAVE_${symbol})
    message( STATUS "\${HAVE_${symbol}} = ${HAVE_${symbol}}")
    message( STATUS
      "Note: Failed Images not supported by the current MPI implementation! (Needs MPIX experimental features--as of MPI3)")
    set(MPI_HAS_FAULT_TOL_EXT NO)
    break() # no need to keep looking
  endif()
endforeach(symbol)
set(CMAKE_REQUIRED_INCLUDES ${old_cmake_required_includes})
set(CMAKE_REQUIRED_FLAGS ${old_cmake_required_flags})
set(CMAKE_REQUIRED_LIBRARIES ${old_cmake_required_libraries})

if(MPI_HAS_FAULT_TOL_EXT) # AND (NOT openmpi))
  option(CAF_ENABLE_FAILED_IMAGES "Enable failed images support" FALSE)
  message(STATUS "The MPI implementation appears to have the experimental features for ULFM
that will allow you to build OpenCoarrays with failed images support. However,
ULFM support did not make it into the MPI-4 standard, and it is known to
trigger some bugs. Because of this we have disabled it by default. You may
add the `-DCAF_ENABLE_FAILED_IMAGES:BOOL=ON` CMake flag or edit the value
with ccmake or cmake-gui if you would like to experiment with failed images.")
else()
  set(CAF_ENABLE_FAILED_IMAGES FALSE
    CACHE
    BOOL
    "Enable failed images support (no support in the selected MPI implementation)"
    FORCE)
endif()

if(CAF_ENABLE_FAILED_IMAGES)
  foreach(lib caf_mpi caf_mpi_static)
    target_compile_definitions(${lib}
      PUBLIC -DUSE_FAILED_IMAGES)
  endforeach()
endif()

#---------------------------------------------------
# Windows Intel MPI compatibility, see GH issue #435
#---------------------------------------------------
set(old_cmake_required_includes "${CMAKE_REQUIRED_INCLUDES}")
if(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${MPI_C_INCLUDE_PATH})
else()
  set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
endif()
CHECK_INCLUDE_FILES("mpi.h" HAVE_MPI_H)
CHECK_SYMBOL_EXISTS(I_MPI_VERSION "mpi.h" HAVE_Intel_MPI)
if(HAVE_Intel_MPI AND WIN32)
  foreach(lib caf_mpi caf_mpi_static)
    target_compile_definitions(${lib}
      PUBLIC -DUSE_GCC)
  endforeach()
endif()
set(CMAKE_REQUIRED_INCLUDES ${old_cmake_required_includes})

##############################################
# Configure `caf` and `cafrun` wrapper scripts
##############################################
# List of caf.in variables needing configuration:
#
# @CAF_VERSION@ @opencoarrays_aware_compiler@ @Fortran_COMPILER@ @CAF_MODDIR@
# @CAF_MPI_Fortran_LINK_FLAGS@ @CAF_MPI_Fortran_COMPILE_FLAGS@
# @CAF_LIBS@ @CAF_MPI_LIBS@
#

set(CAF_VERSION "${full_git_describe}")
set(Fortran_COMPILER "${CMAKE_Fortran_COMPILER}")
set(CAF_MODDIR "${CMAKE_INSTALL_INCLUDEDIR}/${mod_dir_tail}")
set(MOD_DIR_FLAG "${CMAKE_Fortran_MODDIR_FLAG}")
set(CAF_MPI_LIBS "")
foreach( lib IN LISTS MPI_Fortran_LIBRARIES)
  set(CAF_MPI_LIBS "${CAF_MPI_LIBS} \"${lib}\"")
endforeach()
string(STRIP "${CAF_MPI_LIBS}" CAF_MPI_LIBS)
set(CAF_MPI_Fortran_LINK_FLAGS "")
foreach( lflag IN LISTS MPI_Fortran_LINK_FLAGS)
  set(CAF_MPI_Fortran_LINK_FLAGS "${CAF_MPI_Fortran_LINK_FLAGS} ${lflag}" )
endforeach()
string(STRIP "${CAF_MPI_Fortran_LINK_FLAGS}" CAF_MPI_Fortran_LINK_FLAGS)
set(CAF_MPI_Fortran_COMPILE_FLAGS "")
foreach( fcflag IN LISTS MPI_Fortran_COMPILE_FLAGS)
  set(CAF_MPI_Fortran_COMPILE_FLAGS "${CAF_MPI_Fortran_COMPILE_FLAGS} ${fcflag}" )
endforeach()
string(STRIP "${CAF_MPI_Fortran_COMPILE_FLAGS}" CAF_MPI_Fortran_COMPILE_FLAGS)
set_target_properties(caf_mpi_static
  PROPERTIES OUTPUT_NAME caf_mpi)
get_target_property(libcaf_static caf_mpi_static OUTPUT_NAME)
set(CAF_LIBS "${CMAKE_INSTALL_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}${libcaf_static}${CMAKE_STATIC_LIBRARY_SUFFIX}")

configure_file("${CMAKE_SOURCE_DIR}/src/extensions/caf.in" "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/caf"
  @ONLY)


# List of carrun.in variables needing configuration:
#
# @CAF_VERSION@ @MPIEXEC@ @MPIEXEC_NUMPROC_FLAG@ @MPIEXEC_PREFLAGS@ @MPIEXEC_POSTFLAGS@
# @HAVE_FAILED_IMG@
#

if(CAF_ENABLE_FAILED_IMAGES)
  set(HAVE_FAILED_IMG "true")
else()
  set(HAVE_FAILED_IMG "false")
endif()

configure_file("${CMAKE_SOURCE_DIR}/src/extensions/cafrun.in" "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/cafrun"
  @ONLY)

install(PROGRAMS "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/cafrun"
  "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/caf"
  DESTINATION "${CMAKE_INSTALL_FULL_BINDIR}")

lint_script("${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}" caf)
check_script_style("${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/caf")
lint_script("${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}" cafrun)
check_script_style("${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/cafrun")
