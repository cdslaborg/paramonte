add_executable(caf_psnap cafpsnap.f90 timemeasure.c)
target_link_libraries(caf_psnap OpenCoarrays)
