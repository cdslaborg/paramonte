IF (("${CMAKE_SYSTEM_NAME}" MATCHES "Linux" ) AND (${HIGH_RESOLUTION_TIMER}))
  message("\nUsing assembly-language timers that tick once per clock cycle.")
  message("WARNING: walltime.o is designed for an x86 CPU with a 3.6 GHz clock.")
  message("Because modern processors can throttle their clock spees, do not ")
  message("rely on this timer for an absolute value, but it can be useful for ")
  message("comparison or relative execution execution times on the same platform.\n")
  add_executable(mpi_distributed_transpose
    mpi_distributed_transpose.F90
    walltime.o
    )
  target_compile_definitions(mpi_distributed_transpose
    PRIVATE -DHAVE_WALLTIME)
else()
  add_executable(mpi_distributed_transpose
    mpi_distributed_transpose.F90
  )
endif()
