set(directories_to_build mpi tests)
if( CAF_ENABLE_ISO_Fortran_BINDING OR ( CAF_RUN_DEVELOPER_TESTS OR $ENV{OPENCOARRAYS_DEVELOPER} ) )
  list( APPEND directories_to_build iso-fortran-binding )
endif()
foreach(directory ${directories_to_build})
  add_subdirectory(${directory})
endforeach()

set(N_CPU ${N_CPU} PARENT_SCOPE)
set(HOST_NAME ${HOST_NAME} PARENT_SCOPE)

if(openmpi)
  set(openmpi ${openmpi} PARENT_SCOPE)
endif()

install( FILES libcaf.h libcaf-gfortran-descriptor.h libcaf-version-def.h
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
