if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  set(gfortran_compiler true)
endif()

function(generate_test_script base_name num_images)
  set(source ${base_name}.F90)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/../collectives/${source}
    ${CMAKE_CURRENT_BINARY_DIR}/${source}
    COPYONLY
  )
  # Now we write the script that compiles and runs the test
  set(harness "${CMAKE_BINARY_DIR}/bin_staging/test-${base_name}-extension.sh")
  file(WRITE  "${harness}" "#!/bin/bash\n")
  file(APPEND "${harness}" "cd ${CMAKE_CURRENT_BINARY_DIR}\n")
  set(executable "${base_name}_extension")
  if (opencoarrays_aware_compiler)
    # Explicitly include the directory containing the .mod file because the caf script won't pass it
    # automatically when the script wraps an OpenCoarrays-aware compiler
    file(APPEND
      "${harness}"
      "FC=mpifort ${CMAKE_INSTALL_PREFIX}/bin/caf ${source} -o ${executable} -I${CMAKE_BINARY_DIR}/mod -DUSE_EXTENSIONS\n"
    )
  else()
    file(APPEND "${harness}" "FC=mpifort ${CMAKE_INSTALL_PREFIX}/bin/caf ${source} -o ${executable} -DUSE_EXTENSIONS\n")
  endif()
  file(APPEND "${harness}" "${CMAKE_INSTALL_PREFIX}/bin/cafrun -np ${num_images} ./${executable}\n")
  file( INSTALL "${harness}"
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )

endfunction(generate_test_script)

generate_test_script(co_sum 4)
generate_test_script(co_broadcast 4)
generate_test_script(co_min 4)
generate_test_script(co_max 4)
generate_test_script(co_reduce 4)
