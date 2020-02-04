set(include_directory ${CMAKE_CURRENT_SOURCE_DIR}/include-files)
set(library_directory ${CMAKE_CURRENT_SOURCE_DIR}/library)
set(config_directory  ${CMAKE_CURRENT_BINARY_DIR}/library)

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Cray")
  configure_file(${include_directory}/cray_capabilities.txt ${config_directory}/compiler_capabilities.txt COPYONLY)
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  configure_file(${include_directory}/intel_capabilities.txt ${config_directory}/compiler_capabilities.txt COPYONLY)
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  configure_file(${include_directory}/gfortran_capabilities.txt ${config_directory}/compiler_capabilities.txt COPYONLY)
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL")
  configure_file(${include_directory}/ibm_capabilities.txt ${config_directory}/compiler_capabilities.txt COPYONLY)
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "PGI")
  configure_file(${include_directory}/portlandgroup_capabilities.txt ${config_directory}/compiler_capabilities.txt COPYONLY)
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "NAG")
  configure_file(${include_directory}/nag_capabilities.txt ${config_directory}/compiler_capabilities.txt COPYONLY)
else()
  message ("Unknown Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif()

set(CMAKE_Fortran_COMPILER "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/caf")
add_executable(coarray_burgers_pde
  main.F90
  global_field.F90
  local_field.F90
  ${library_directory}/ForTrilinos_assertion_utility.F90
  ${library_directory}/object_interface.F90
  ${library_directory}/co_object_interface.F90
)
add_dependencies(coarray_burgers_pde caf_mpi_static)
target_include_directories(coarray_burgers_pde PRIVATE ${config_directory})
