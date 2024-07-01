####################################################################################################################################
####################################################################################################################################
####                                                                                                                            ####
####    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           ####
####                                                                                                                            ####
####    Copyright (C) 2012-present, The Computational Data Science Lab                                                          ####
####                                                                                                                            ####
####    This file is part of the ParaMonte library.                                                                             ####
####                                                                                                                            ####
####    LICENSE                                                                                                                 ####
####                                                                                                                            ####
####       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          ####
####                                                                                                                            ####
####################################################################################################################################
####################################################################################################################################

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copy the ParaMonte library benchmark build scripts to the deployment benchmark subdirectories.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set(collection "benchmark")
set(${collection}BR "${bench}")
set(${collection}PP "${benchpp}")
set(origin "${paramonte_benchmark_dir}/${lang}/") # The ending forward slash is essential for copying the contents of the ${collection} language folder.
set(destin "${paramonte_bld_pkg_benchmark_dir}")

if (${lang_is_dynamic})
    include(setParaMonteCollectionDynamic)
elseif ("${lang}" STREQUAL "c" OR "${lang}" STREQUAL "cpp" OR "${lang}" STREQUAL "fortran")
    include(setParaMonteCollectionStatic)
else()
    message(FATAL_ERROR "${pmfatal} Internal ParaMonte-CMake error: Unrecognized compiled language: lang=${lang}")
endif()