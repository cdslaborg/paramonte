# Fortran tests

caf_compile_executable(increment_my_neighbor increment_neighbor.f90)
caf_compile_executable(atomics testAtomics.f90)

# C tests
#include(CMakeForceCompiler)
#CMAKE_FORCE_C_COMPILER(mpicc GNU)
#add_executable(C_sync_images syncimages2.c ../../../mpi/mpi_caf.c ../../../common/caf_auxiliary.c)
