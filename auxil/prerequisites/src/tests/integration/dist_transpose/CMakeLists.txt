set(DIST_TRANSP_SRCS coarray_distributed_transpose.F90)
if (HIGH_RESOLUTION_TIMER)
  LIST(APPEND DIST_TRANSP_SRCS -DHAVE_WALLTIME walltime.o)
endif()
caf_compile_executable(coarray_distributed_transpose
  ${DIST_TRANSP_SRCS} )
