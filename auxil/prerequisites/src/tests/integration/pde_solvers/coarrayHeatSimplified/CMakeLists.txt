set(CMAKE_Fortran_COMPILER "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/caf")
add_library(local_field OBJECT local_field.f90)
add_library(global_field OBJECT global_field.f90)
add_dependencies(local_field caf_mpi_static)
add_dependencies(global_field local_field caf_mpi_static)
add_executable(co_heat
  main.f90
  $<TARGET_OBJECTS:local_field>
  $<TARGET_OBJECTS:global_field>
)
add_dependencies(co_heat caf_mpi_static)
