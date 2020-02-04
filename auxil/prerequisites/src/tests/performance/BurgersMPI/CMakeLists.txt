set(include_directory ${CMAKE_CURRENT_SOURCE_DIR}/../../integration/pde_solvers/include-files)
set(library_directory ${CMAKE_CURRENT_SOURCE_DIR}/../../integration/pde_solvers/library)
set(config_directory  ${CMAKE_CURRENT_BINARY_DIR}/../../integration/pde_solvers/library)

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

add_executable( mpi_burgers_pde
  main.F90
  ${library_directory}/object_interface.F90
  ${library_directory}/ForTrilinos_assertion_utility.F90
  ${library_directory}/ForTrilinos_error.F90
  kind_parameters.F90
  shared.F90
  mpi_module.F90
  mpi_share.F90
  input_file.F90
  periodic_2nd_order.F90
)
target_include_directories(mpi_burgers_pde PRIVATE ${config_directory})
