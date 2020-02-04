# checking whether the machine is of type 64-bit before proceeding further
if ( ("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "x86_64") AND ("${CMAKE_SYSTEM_NAME}" MATCHES "Linux") )
  # Default to older SSE-instruction-based FFT library
  if (NOT (DEFINED ENV{TRAVIS}) AND (NOT SKIP_ASSEMBLY_DEPS) )
    if (LEGACY_ARCHITECTURE OR (NOT DEFINED(LEGACY_ARCHITECTURE)))
      set(fft_library ${CMAKE_CURRENT_SOURCE_DIR}/libfft_sse.a )
    else()
      message(WARNING "Attempting to use libfft_avx.a, which only works on recent x86 architectures.")
      set(fft_library ${CMAKE_CURRENT_SOURCE_DIR}/libfft_avx.a )
    endif()
    caf_compile_executable( coarray_navier_stokes
      coarray-shear_coll.F90
      ${walltime_o}
      ${fft_library}
    )
  endif()
else()
  # Skip Navier-Stokes build until an appropriate FFT has been found.
endif()
