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
# Copy the ParaMonte library collection build scripts to the collection package subdirectories.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set(collection "example")
set(${collection}BR "${exam}")
set(${collection}PP "${exampp}")
set(origin "${paramonte_example_dir}/${lang}/") # The ending forward slash is essential for copying the contents of the ${collection} language folder.
set(destin "${paramonte_bld_pkg_example_dir}")

if (${lang_is_dynamic})
    include(setParaMonteCollectionDynamic)
elseif ("${lang}" STREQUAL "c" OR "${lang}" STREQUAL "cpp" OR "${lang}" STREQUAL "fortran")
    include(setParaMonteCollectionStatic)
else()
    message(FATAL_ERROR "${pmfatal} Internal ParaMonte-CMake error: Unrecognized compiled language: lang=${lang}")
endif()