caf_compile_executable(initialize_mpi initialize_mpi.F90)
caf_compile_executable(register register.f90)
set_target_properties(build_register
  PROPERTIES MIN_IMAGES 2)
caf_compile_executable(register_vector register_vector.f90)
caf_compile_executable(register_alloc_vector register_alloc_vector.f90)
caf_compile_executable(allocate_as_barrier allocate_as_barrier.f90)
caf_compile_executable(allocate_as_barrier_proc allocate_as_barrier_proc.f90)

if (gfortran_compiler AND (NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.0.0))
  caf_compile_executable(async_comp_alloc async_comp_alloc.f90)
  set_target_properties(build_async_comp_alloc
    PROPERTIES MIN_IMAGES 2)
  caf_compile_executable(async_comp_alloc_2 async_comp_alloc_2.f90)
  caf_compile_executable(register_alloc_comp_1 register_alloc_comp_1.f90)
  caf_compile_executable(register_alloc_comp_2 register_alloc_comp_2.f90)
  caf_compile_executable(register_alloc_comp_3 register_alloc_comp_3.f90)
  caf_compile_executable(comp_allocated_1 comp_allocated_1.f90)
  caf_compile_executable(comp_allocated_2 comp_allocated_2.f90)
elseif((CAF_RUN_DEVELOPER_TESTS OR $ENV{OPENCOARRAYS_DEVELOPER}))
  message( AUTHOR_WARNING "Skipping building the following tests due to GFortran < 7.x lack of compatibility:
    async_comp_alloc.f90
    async_comp_alloc_2.f90
    register_alloc_comp_1.f90
    register_alloc_comp_2.f90
    register_alloc_comp_3.f90
    comp_allocated_1.f90
    comp_allocated_2.f90" )
endif()
