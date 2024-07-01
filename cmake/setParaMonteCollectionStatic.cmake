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

if ("${collection}" STREQUAL "example")
    set(isCollectionExam TRUE)
else()
    set(isCollectionExam FALSE)
endif()

if ("${collection}" STREQUAL "benchmark")
    set(isCollectionBench TRUE)
else()
    set(isCollectionBench FALSE)
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Infer the runtime shell type: windows / posix
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set(OS "$ENV{OS}")
string(TOLOWER "${OS}" OS)
if ("${OS}" MATCHES ".*windows.*")
    set(shellexec "CMD /C")
    set(paramonte_collection_build_script "build.bat")
    message(NOTICE "${pmattn} Assuming Windows-compliant runtime shell as requested for build ${collection}...")
else()
    set(shellexec "sh -c")
    set(paramonte_collection_build_script "build.sh")
    message(NOTICE "${pmattn} Assuming POSIX-compliant runtime shell as requested for build ${collection}...")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# The rule for the ParaMonte library ${collection} build. Build ${collection} only if requested.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# This phony target is only relevant for building and running ${collection}.
add_custom_target(${collection})

message(NOTICE "${pmattn} Setting the rules for copying the ParaMonte library ${collection} in ${lang} language to the installation directory...")
message(NOTICE "${pmattn} ${collection}BR=${${collection}BR}")
message(NOTICE "${pmattn} ${collection}PP=${${collection}PP}")

if (EXISTS "${origin}")

    #string(TOLOWER "${${collection}BR}" ${collection}BR_lower)
    #string(TOLOWER "${${collection}PP}" ${collection}PP_lower)
    if (NOT "${${collection}BR}" STREQUAL "" AND NOT "${${collection}BR}" MATCHES "[Nn][Oo][Nn][Ee]")
        set(${collection}BREnabled ON)
    else()
        set(${collection}BREnabled OFF)
    endif()
    if (NOT "${${collection}BR}" STREQUAL "" AND NOT "${${collection}PP}" MATCHES "[Nn][Oo][Nn][Ee]")
        set(${collection}PPEnabled ON)
    else()
        set(${collection}PPEnabled OFF)
    endif()

    if (${collection}BREnabled)
        message(NOTICE "${pmattn} Enabling the ParaMonte library ${collection} build as requested...")
    else()
        message(NOTICE "${pmwarn} Skipping the ParaMonte library ${collection} build as requested...")
    endif()

    if (${collection}PPEnabled)
        message(NOTICE "${pmattn} Enabling the ParaMonte library ${collection} post-processing build as requested...")
    else()
        message(NOTICE "${pmwarn} Skipping the ParaMonte library ${collection} post-processing build as requested...")
    endif()

    # Adding a custom target runs the rule only at compile time,
    # which subsequently does not allow the addition of the build scripts to the ${collection} folders,
    # because the ${collection} folders have not been added yet. The above approach has the limitation that
    # it always requires build generation via CMake to update the ${collection} files and the build scripts.
    #add_custom_target(${collection} ALL COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}")

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Copy the ${collection} files from the ${origin} to the ${destin} directory.
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    message(NOTICE  "${pmattn} Copying the paramonte::${lang} library ${collection} files to the installation directory...\n"
                    "${pmattn}      - from: ${origin}\n"
                    "${pmattn}      -   to: ${destin}")
    #add_custom_target(${collection} COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}")
    file(COPY "${origin}" DESTINATION "${destin}")

    if (NOT ${lang_is_dynamic})
        if (EXISTS "${paramonte_${collection}_dir}/generic/")
            message(NOTICE  "${pmattn} Copying the paramonte::${lang} library generic ${collection} files to the installation directory...\n"
                            "${pmattn}      - from: ${paramonte_example_dir}/generic\n"
                            "${pmattn}      -   to: ${destin}")
            file(COPY "${paramonte_${collection}_dir}/generic/" DESTINATION "${destin}")
        else()
            message(NOTICE  "${pmwarn} The paramonte library generic ${collection} folder does not exists. skipping...")
        endif()
    endif()

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Generate the contents of the build scripts.
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #   For C/C++ examples:
    #
    #       We generate a `CMakeLists.txt` and the bash and batch scripts to call Cmake and run the example.
    #       This approach is essential because the C/C++ compiler name and flags can be other than those of gnu/Intel compilers.
    #
    #   For Fortran examples:
    #
    #       We generate only the bash and batch build scripts without CMake files because only two Fortran compiler vendors are supported.
    #       This approach is currently pragmatic and clean, but will have to be extended to CMake generic build files with more compilers.
    #
    #   For dynamic-language examples:
    #
    #       We generate only the bash and batch execution scripts.
    #
    unset(collection_cmakelists_path)
    unset(collection_cmakelists_contents)
    unset(collection_cmakelists_bash_path)
    unset(collection_cmakelists_batch_path)
    unset(collection_cmakelists_bash_contents)
    unset(collection_cmakelists_batch_contents)
    set(collection_bld_bash_path "${CMAKE_CURRENT_BINARY_DIR}/build.sh")
    set(collection_bld_batch_path "${CMAKE_CURRENT_BINARY_DIR}/build.bat")

    #### Build the contents of the collection CMakeList.txt.

    if ("${lang}" STREQUAL "c")
        set(binary_lang "C")
    elseif ("${lang}" STREQUAL "cpp")
        set(binary_lang "CXX")
    elseif ("${lang}" STREQUAL "fortran")
        set(binary_lang "Fortran")
    else()
        message(FATAL_ERROR "${pmfatal} Internal ParaMonte-CMake error: Unrecognized compiled language: lang=${lang}")
    endif()

    #### get the compiler name.

    #get_filename_component(binary_compiler "${CMAKE_${binary_lang}_COMPILER}" NAME)
    set(binary_compiler "${CMAKE_${binary_lang}_COMPILER}")

    #### Build the contents of CMakeLists.

    # WARNING
    # The first paramonet library path search option below in `find_library` is
    # essential for correct (original) library identification in code coverage builds.

    string(CONCAT collection_cmakelists_contents "${collection_cmakelists_contents}"
        "cmake_minimum_required(VERSION 3.14)\n"
        "set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE})\n"
        "if (\"\${CMAKE_Fortran_COMPILER}\" STREQUAL \"\")\n"
        "    set(CMAKE_Fortran_COMPILER \"${binary_compiler}\")\n"
        "endif()\n"
        "project(main LANGUAGES ${binary_lang})\n"
        "set(CMAKE_MACOSX_RPATH ON)\n"
        "set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)\n"
        "add_executable(binary main${lang_ext})\n"
        "set_target_properties(binary PROPERTIES OUTPUT_NAME \"main\" SUFFIX \".exe\")\n"
        "target_include_directories(binary PUBLIC \"\${CMAKE_CURRENT_SOURCE_DIR}/../../../inc\")\n"
        "find_library(pmlib NAMES ${libname} ${libname}.a ${libname}.dll ${libname}.dylib ${libname}.lib ${libname}.so PATHS\n"
        "            \"\${CMAKE_CURRENT_SOURCE_DIR}/../../../../lib\" # library directory for code coverage.\n"
        "            \"\${CMAKE_CURRENT_SOURCE_DIR}/../../../lib\" # library directory for examples/benchmarks.\n"
        "            )\n"
        "target_link_libraries(binary PUBLIC \"\${pmlib}\")\n"
        "target_link_libraries(binary PUBLIC \"\${pmlib}\")\n"
        "if (APPLE)\n"
        "    set(rpath_prop \"@loader_path\")\n"
        "elseif(UNIX)\n"
        "    set(rpath_prop \"$ORIGIN\")\n"
        "endif()\n"
        "set_target_properties(binary PROPERTIES INSTALL_RPATH \"${rpath_prop}\")\n"
    )
    #if (${csid_is_gnu} AND ${codecov_enabled})
    #    string(CONCAT collection_cmakelists_contents "${collection_cmakelists_contents}"
    #    "find_library(pmlib NAMES ${libname} ${libname}.a ${libname}.dll ${libname}.dylib ${libname}.lib ${libname}.so PATHS \"\${CMAKE_CURRENT_SOURCE_DIR}/../../../../lib\" \"\${CMAKE_CURRENT_SOURCE_DIR}/../../../lib\")\n"
    #    )
    #else()
    #    string(CONCAT collection_cmakelists_contents "${collection_cmakelists_contents}"
    #    "find_library(pmlib NAMES ${libname} ${libname}.a ${libname}.dll ${libname}.dylib ${libname}.lib ${libname}.so PATHS \"\${CMAKE_CURRENT_SOURCE_DIR}/../../../lib\")\n"
    #    )
    #endif()

        #set(CMAKE_MACOSX_RPATH 1)
        #set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

    unset(mkl_flag_bash)
    unset(mkl_flag_batch)
    unset(binary_ipo_options)
    unset(binary_link_options)
    unset(binary_compile_options)
    unset(binary_compile_definitions)
    #if (${csid_is_gnu} OR ${csid_is_intel})
        #set(binary_link_options "-Wl,-rpath,../../../lib")
    #endif()

    if ("${lang}" STREQUAL "c" OR "${lang}" STREQUAL "cpp")

        string( CONCAT
                collection_cmakelists_contents
                "${collection_cmakelists_contents}"
                "if (\"\${CMAKE_${binary_lang}_COMPILER_ID}\" STREQUAL \"GNU\")\n"
                "    target_link_libraries(binary PUBLIC m)\n"
                "endif()\n"
                )

    elseif ("${lang}" STREQUAL "fortran")

        if (${csid_is_gnu})

            set(binary_ipo_options -flto=3)
            set(binary_compile_options "${binary_compile_options}" -cpp -ffree-line-length-none)
            if ("${build}" STREQUAL "testing")
                set(binary_compile_options "${binary_compile_options}" -O2)
            elseif ("${build}" STREQUAL "debug")
                set(binary_compile_options "${binary_compile_options}" -O0 -g -fcheck=all -fbacktrace)
            else()
                set(binary_compile_options "${binary_compile_options}" -O3)
            endif()

            if (${codecov_enabled})
                set(binary_link_options "${binary_link_options}" "--coverage")
                set(binary_compile_options "${binary_compile_options}" "--coverage")
            endif()

        elseif (${csid_is_intel})

            if (WIN32)
                set(binary_ipo_options /Qipo)
                set(mkl_compile_options /Qmkl:sequential)
                set(mkl_compile_definitions /DMKL_ENABLED=1)
                #set(binary_link_options "${binary_link_options}" /F0x1000000000) The intel /F flag does not seem to work anymore.
                set(binary_compile_options "${binary_compile_options}" /standard-semantics /fpp) # /F0x1000000000 /heap-arrays We use editbin software now to increase stack size.
                if ("${build}" STREQUAL "testing")
                    set(binary_compile_options "${binary_compile_options}" /O2)
                elseif ("${build}" STREQUAL "debug")
                    set(binary_compile_options "${binary_compile_options}" /Od /debug:full /CB /Qinit:snan,arrays /warn:all /gen-interfaces /traceback /check:all /fpe-all:0 /Qtrapuv)
                else()
                    set(binary_compile_options "${binary_compile_options}" /O3)
                endif()
            else()
                set(binary_ipo_options -ipo)
                set(mkl_compile_options -qmkl=sequential)
                set(mkl_compile_definitions -DMKL_ENABLED=1)
                set(binary_compile_options "${binary_compile_options}" -fpp -standard-semantics)
                if ("${build}" STREQUAL "testing")
                    set(binary_compile_options "${binary_compile_options}" -O2)
                elseif ("${build}" STREQUAL "debug")
                    set(binary_compile_options "${binary_compile_options}" -O0 -g3 -CB -debug full -traceback -check all -fpe0)
                else()
                    set(binary_compile_options "${binary_compile_options}" -O3)
                endif()
            endif()

        endif()

    endif()

    if (DEFINED binary_ipo_options AND ("${build}" STREQUAL "ipo" OR "${build}" STREQUAL "tuned" OR "${build}" STREQUAL "native"))
        set(binary_compile_options "${binary_compile_options}" "${binary_ipo_options}")
        set(binary_link_options "${binary_link_options}" "${binary_ipo_options}")
    endif()

    if (OpenBLAS_ENABLED)
        # Warning: The following openblas library name will have to updated to `openblas_64`
        # if the compiler options for long integer kinds for array size are enabled.
        string( CONCAT
                collection_cmakelists_contents
                "${collection_cmakelists_contents}"
                "find_library(oblib NAMES openblas libopenblas PATHS \"\${CMAKE_CURRENT_SOURCE_DIR}/../../../lib\")\n"
                "set_property(TARGET binary APPEND PROPERTY LINK_LIBRARIES \"${oblib}\")\n"
                )
    elseif (BLAS_ENABLED OR LAPACK_ENABLED)
        if ("${lang}" STREQUAL "fortran" AND ${csid_is_intel})
            set(binary_compile_options "${binary_compile_options} ${mkl_compile_options}")
            set(binary_compile_definitions "${binary_compile_definitions} ${mkl_compile_definitions}")
        else()
            if (BLAS_ENABLED AND BLAS_FOUND AND DEFINED BLAS_LIBRARIES)
                string( CONCAT
                        collection_cmakelists_contents
                        "${collection_cmakelists_contents}"
                        "set_property(TARGET binary APPEND PROPERTY LINK_LIBRARIES \"${BLAS_LIBRARIES}\")\n"
                        )
            elseif (LAPACK_ENABLED AND LAPACK_FOUND AND DEFINED LAPACK_LIBRARIES)
                string( CONCAT
                        collection_cmakelists_contents
                        "${collection_cmakelists_contents}"
                        "set_property(TARGET binary APPEND PROPERTY LINK_LIBRARIES \"${LAPACK_LIBRARIES}\")\n"
                        )
            endif()
        endif()
    endif()

    if (DEFINED binary_compile_definitions)
        string( CONCAT
                collection_cmakelists_contents
                "${collection_cmakelists_contents}"
                "set_property(TARGET binary APPEND PROPERTY COMPILE_DEFINITIONS \"${binary_compile_definitions}\")\n"
                )
    endif()

    if (DEFINED binary_compile_options)
        string( CONCAT
                collection_cmakelists_contents
                "${collection_cmakelists_contents}"
                "set_property(TARGET binary APPEND PROPERTY COMPILE_OPTIONS \"${binary_compile_options}\")\n"
                )
    endif()

    if (DEFINED binary_link_options)
        string( CONCAT
                collection_cmakelists_contents
                "${collection_cmakelists_contents}"
                "set_property(TARGET binary APPEND PROPERTY LINK_OPTIONS \"${binary_link_options}\")\n"
                )
    endif()

    #### Add parallelism instructions and libraries to the CMakeLists file if needed.

    if(WIN32)
        set(bincmd "main.exe")
    else()
        set(bincmd "./main.exe")
    endif()
    if(MPI_ENABLED)
        string( CONCAT
                collection_cmakelists_contents
                "${collection_cmakelists_contents}"
                "find_package(MPI)\n"
                "if(MPI_${binary_lang}_FOUND)\n"
                "    target_link_libraries(binary PUBLIC MPI::MPI_${binary_lang})\n"
                "else()\n"
                "    message(WARNING \"CMake could not find an MPI library in the current environment.\")\n"
                "    message(WARNING \"The example build & run may fail.\")\n"
                "endif()\n"
                )
        if (WIN32 AND ("${CMAKE_${binary_lang}_COMPILER_ID}" MATCHES ".*Intel.*"))
            set(bincmd "mpiexec -localonly -n ${nproc} ${bincmd}")
        else()
            set(bincmd "mpiexec -n ${nproc} ${bincmd}")
        endif()
    elseif(OMP_ENABLED)
        string( CONCAT
                collection_cmakelists_contents
                "${collection_cmakelists_contents}"
                "find_package(OpenMP)\n"
                "if(OpenMP_${binary_lang}_FOUND)\n"
                "    target_link_libraries(binary PUBLIC OpenMP::OpenMP_${binary_lang})\n"
                "else()\n"
                "    message(WARNING \"CMake could not find an OpenMP library in the current environment.\")\n"
                "    message(WARNING \"The example build & run may fail.\")\n"
                "endif()\n"
                )
    endif()

    #### Add the example execution make rule.

    string(CONCAT
        collection_cmakelists_contents
        "${collection_cmakelists_contents}"
        "add_custom_target(run COMMAND ${bincmd} DEPENDS binary WORKING_DIRECTORY \"\$\{CMAKE_CURRENT_SOURCE_DIR\}\")\n"
        )
    set(collection_cmakelists_bash_name "build.sh") # cmake.sh
    set(collection_cmakelists_batch_name "build.bat") # cmake.bat
    set(collection_cmakelists_path "${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt")
    set(collection_cmakelists_bash_path "${CMAKE_CURRENT_BINARY_DIR}/${collection_cmakelists_bash_name}")
    set(collection_cmakelists_batch_path "${CMAKE_CURRENT_BINARY_DIR}/${collection_cmakelists_batch_name}")

    set(binary_cmake_build "cmake . -G \"${CMAKE_GENERATOR}\" && ${CMAKE_MAKE_PROGRAM}")

    set(collection_cmakelists_bash_contents "#!/bin/bash\n")
    if (WIN32)
        string( CONCAT
                collection_cmakelists_bash_contents
                "${collection_cmakelists_bash_contents}"
                "# Add the library path to the\n"
                "# PATH environment variable.\n"
                "export PATH=../../../lib:$PATH\n"
                )
    endif()
    string( CONCAT
            collection_cmakelists_bash_contents
            "${collection_cmakelists_bash_contents}"
            "${binary_cmake_build}\n"
            "${CMAKE_MAKE_PROGRAM} run"
            )

    string( CONCAT
            collection_cmakelists_batch_contents
            "REM Add the library path to the\n"
            "REM PATH environment variable.\n"
            "setlocal EnableDelayedExpansion\n"
            "set \"PATH=..\\..\\..\\lib;%PATH%\"\n"
            "REM Enlarge the library stack size.\n"
            "editbin ..\\..\\..\\lib\\*.dll /stack:1000000000\n"
            "REM Generate build scripts.\n"
            "${binary_cmake_build}\n"
            "REM Enlarge the binary stack size.\n"
            "editbin *.exe /stack:1000000000\n"
            "REM Run the binary.\n"
            "${CMAKE_MAKE_PROGRAM} run"
            )

    #### Write the files.

    if (DEFINED collection_cmakelists_path AND DEFINED collection_cmakelists_contents)
        file(WRITE "${collection_cmakelists_path}" "${collection_cmakelists_contents}")
        message(NOTICE "${pmattn} ${collection} CMakeLists.txt script:")
        message(NOTICE "${collection_cmakelists_contents}")
    endif()
    if (DEFINED collection_cmakelists_bash_path AND DEFINED collection_cmakelists_bash_contents)
        file(WRITE "${collection_cmakelists_bash_path}" "${collection_cmakelists_bash_contents}")
        message(NOTICE "${pmattn} ${collection} CMakeLists Bash build script:")
        message(NOTICE "${collection_cmakelists_bash_contents}")
    endif()
    if (WIN32 AND DEFINED collection_cmakelists_batch_path AND DEFINED collection_cmakelists_batch_contents)
        file(WRITE "${collection_cmakelists_batch_path}" "${collection_cmakelists_batch_contents}")
        message(NOTICE "${pmattn} ${collection} CMakeLists Batch build script:")
        message(NOTICE "${collection_cmakelists_batch_contents}")
    endif()

    unset(collection_cmakelists_batch_contents)
    unset(collection_cmakelists_bash_contents)
    unset(collection_cmakelists_contents)

    #if ("${lang}" STREQUAL "c" OR "${lang}" STREQUAL "cpp")
    #
    #    #set(collection_bld_bash_path "${CMAKE_CURRENT_BINARY_DIR}/${collection_cmakelists_bash_name}")
    #    #set(collection_bld_batch_path "${CMAKE_CURRENT_BINARY_DIR}/${collection_cmakelists_batch_name}")
    #
    #elseif ("${lang}" STREQUAL "fortran")
    #
    #    if (${csid_is_gnu} OR ${csid_is_intel})
    #
    #        # Rather than copying the build script files from the ${collection}/auxil folder, generate the files on the fly.
    #
    #        # Set the first line of the build scripts.
    #        set(collection_bld_bash_pre "#!/bin/bash")
    #        set(collection_bld_batch_pre "set PATH=..\\..\\..\\lib;%PATH%")
    #
    #        # Set the library name for ${collection} build scripts.
    #        set(collection_bld_bash_compile "-o main.exe main${lang_ext} ../../../lib/${libname}*")
    #        if (OpenBLAS_ENABLED)
    #            set(collection_bld_bash_compile "${collection_bld_bash_compile} ../../../lib/openblas*")
    #        endif()
    #        if (WIN32 AND "${lib}" STREQUAL "shared" AND ${csid_is_intel})
    #            set(collection_bld_batch_compile "/exe:main.exe main${lang_ext} ..\\..\\..\\lib\\${libname}*.lib")
    #        else()
    #            set(collection_bld_batch_compile "/exe:main.exe main${lang_ext} ..\\..\\..\\lib\\${libname}*")
    #        endif()
    #        unset(libname)
    #
    #        # Generate the generic build scripts.
    #
    #        unset(mkl_flag_bash)
    #        unset(mkl_flag_batch)
    #
    #        # Set the last line of the build scripts.
    #
    #        if (${MPI_ENABLED})
    #            set(collection_bld_bash_post "mpiexec -n ${nproc} ./main.exe")
    #            if (${csid_is_intel})
    #                set(collection_bld_batch_post "mpiexec -localonly -n ${nproc} ./main.exe") #${MPIEXEC_EXECUTABLE}
    #                set(fcname "mpiifort")
    #            else()
    #                set(collection_bld_batch_post "mpiexec -n ${nproc} main.exe")
    #                set(fcname "mpifort")
    #            endif()
    #        else()
    #            get_filename_component(fcname "${CMAKE_Fortran_COMPILER}" NAME)
    #            set(collection_bld_bash_post "./main.exe")
    #            set(collection_bld_batch_post "main.exe")
    #        endif()
    #
    #        if (isCollectionBench)
    #            set(collection_bld_bash_compile "-DBLAS_ENABLED=${BLAS_ENABLED} -DLAPACK_ENABLED=${LAPACK_ENABLED} ${collection_bld_bash_compile}")
    #            set(collection_bld_batch_compile "/DBLAS_ENABLED=${BLAS_ENABLED} /DLAPACK_ENABLED=${LAPACK_ENABLED} ${collection_bld_batch_compile}")
    #        endif()
    #
    #        if (${csid_is_gnu})
    #
    #            set(collection_bld_bash_compile "-cpp -ffree-line-length-none -Wl,-rpath,../../../lib -I../../../inc ${collection_bld_bash_compile}")
    #            if ("${build}" STREQUAL "debug")
    #                set(collection_bld_bash_compile "${fcname} -O0 -g -fcheck=all -fbacktrace ${collection_bld_bash_compile}")
    #            elseif ("${build}" STREQUAL "testing")
    #                set(collection_bld_bash_compile "${fcname} -O2 ${collection_bld_bash_compile}")
    #            else()
    #                set(collection_bld_bash_compile "${fcname} -O3 -flto=3 ${collection_bld_bash_compile}")
    #            endif()
    #
    #        elseif (${csid_is_intel})
    #
    #            # \todo
    #            # The windows case needs to be tested and configured on Windows system.
    #            set(mkl_flag_bash " -qmkl=sequential -DMKL_ENABLED=1") # The initial white space is essential.
    #            set(mkl_flag_batch " /Qmkl:sequential /DMKL_ENABLED=1") # The initial white space is essential.
    #
    #            set(collection_bld_bash_compile "-standard-semantics -fpp -Wl,-rpath,../../../lib -I../../../inc ${collection_bld_bash_compile}")
    #            if ("${build}" STREQUAL "debug")
    #                set(collection_bld_bash_compile "${fcname} -O0 -g3 -CB -debug full -traceback -check all -fpe0 ${collection_bld_bash_compile}")
    #            elseif ("${build}" STREQUAL "testing")
    #                set(collection_bld_bash_compile "${fcname} -O2 ${collection_bld_bash_compile}")
    #            else()
    #                set(collection_bld_bash_compile "${fcname} -O3 -ipo ${collection_bld_bash_compile}")
    #            endif()
    #
    #            set(collection_bld_batch_compile "/standard-semantics /fpp /I:..\\..\\..\\inc  ${collection_bld_batch_compile}")
    #            if ("${build}" STREQUAL "debug")
    #                set(collection_bld_batch_compile "${fcname} /Od /debug:full /CB /Qinit:snan,arrays /warn:all /gen-interfaces /traceback /check:all /fpe-all:0 /Qtrapuv ${collection_bld_batch_compile}")
    #            elseif ("${build}" STREQUAL "testing")
    #                set(collection_bld_batch_compile "${fcname} /O2 ${collection_bld_batch_compile}")
    #            else()
    #                set(collection_bld_batch_compile "${fcname} /O3 /Qipo ${collection_bld_batch_compile}")
    #            endif()
    #
    #        endif()
    #        unset(fcname)
    #
    #        # Generate the full ${collection} build scripts.
    #
    #        if (isCollectionBench)
    #            # Generate the full ${collection} build scripts with BLAS/LAPACK linking.
    #            file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/BLAPACK")
    #            set(collection_bld_bash_contents "${collection_bld_bash_compile}")
    #            if (OpenBLAS_ENABLED)
    #                set(collection_bld_bash_contents "${collection_bld_bash_pre}\n${collection_bld_bash_compile}\n${collection_bld_bash_post}")
    #            elseif (BLAS_ENABLED AND LAPACK_ENABLED)
    #                set(collection_bld_bash_contents "${collection_bld_bash_pre}\n${collection_bld_bash_compile} ${mkl_flag_bash} -llapack -lblas\n${collection_bld_bash_post}")
    #            elseif (BLAS_ENABLED)
    #                set(collection_bld_bash_contents "${collection_bld_bash_pre}\n${collection_bld_bash_compile} ${mkl_flag_bash} -lblas\n${collection_bld_bash_post}")
    #            else()
    #                set(collection_bld_bash_contents "${collection_bld_bash_pre}\n${collection_bld_bash_compile}\n${collection_bld_bash_post}")
    #            endif()
    #        else()
    #            set(collection_bld_bash_contents "${collection_bld_bash_pre}\n${collection_bld_bash_compile}\n${collection_bld_bash_post}")
    #        endif()
    #        set(collection_bld_batch_contents "${collection_bld_batch_pre}\n${collection_bld_batch_compile}\n${collection_bld_batch_post}")
    #
    #    else()
    #
    #        # \todo XXX this needs improvement.
    #        message(WARNING "${pmwarn} Unrecognized unsupported compiler detected. ParaMonte ${collection} build scripts will not function.")
    #
    #    endif()
    #
    #endif()

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Copy the build scripts to the ${collection} folders.
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # \todo
    # The following mess is due to the lack of support for CMake write/copy with custom file attributes.
    # Ideally, instead of having file(WRITE...) followed by many file(COPY...) to each example folder,
    # we could simply have many file(WRITE...) to all example folders.

    if (DEFINED collection_bld_batch_path AND DEFINED collection_bld_batch_contents)
        file(WRITE "${collection_bld_batch_path}" "${collection_bld_batch_contents}")
        message(NOTICE "${pmattn} ${collection} Batch build script:")
        message(NOTICE "${collection_bld_batch_contents}")
    endif()
    if (DEFINED collection_bld_bash_path AND DEFINED collection_bld_bash_contents)
        file(WRITE "${collection_bld_bash_path}" "${collection_bld_bash_contents}")
        message(NOTICE "${pmattn} ${collection} Bash build script:")
        message(NOTICE "${collection_bld_bash_contents}")
    endif()
    unset(collection_bld_batch_contents)
    unset(collection_bld_bash_contents)

    # Loop over all directories to add the build scripts.

    setSubDirList(subDirList "${destin}")
    foreach(subDir ${subDirList})
        setSubDirList(subSubDirList "${destin}/${subDir}")
        foreach(subSubDir ${subSubDirList})
            set(collection_current_dir "${destin}/${subDir}/${subSubDir}")
            message(NOTICE "${pmattn} Adding the build scripts to ${collection} at: ${collection_current_dir}")
            if (EXISTS "${collection_cmakelists_path}")
                file(COPY "${collection_cmakelists_path}" DESTINATION "${collection_current_dir}/"
                FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
            endif()
            if (EXISTS "${collection_cmakelists_bash_path}")
                file(COPY "${collection_cmakelists_bash_path}" DESTINATION "${collection_current_dir}/"
                FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
            endif()
            if (EXISTS "${collection_cmakelists_batch_path}")
                file(COPY "${collection_cmakelists_batch_path}" DESTINATION "${collection_current_dir}/"
                FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
            endif()
            if (WIN32 AND EXISTS "${collection_bld_batch_path}")
                file(COPY "${collection_bld_batch_path}" DESTINATION "${collection_current_dir}/"
                FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
            endif()
            if (EXISTS "${collection_bld_bash_path}")
                file(COPY "${collection_bld_bash_path}" DESTINATION "${collection_current_dir}/"
                FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
            endif()
        endforeach()
    endforeach()

    unset(collection_cmakelists_batch_path)
    unset(collection_cmakelists_bash_path)
    unset(collection_cmakelists_path)
    unset(collection_bld_batch_path)
    unset(collection_bld_bash_path)
    unset(collection_current_dir)

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Build and run the ParaMonte ${collection}.
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if(${collection}BREnabled)

        message(NOTICE "${pmattn} Enabling the selected paramonte::${lang} library ${collection} builds...")

        #### Check if an input list of ${collection} folders is specified.

        unset(paramonte_${collection}BR_file)
        # Assume an absolute path specified.
        set(pathtemp "${${collection}BR}")
        if (EXISTS "${pathtemp}" AND NOT IS_DIRECTORY "${pathtemp}")
            set(paramonte_${collection}BR_file "${${collection}BR}")
        else()
            # Assume a relative path with respect to root directory.
            set(pathtemp "${paramonte_dir}/${${collection}BR}")
            if (EXISTS "${pathtemp}" AND NOT IS_DIRECTORY "${pathtemp}")
                set(paramonte_${collection}BR_file "${pathtemp}")
            else()
                # Assume a relative path with respect to ${collection}/${lang} root directory.
                set(pathtemp "${paramonte_dir}/${collection}/${lang}/${${collection}BR}")
                if (EXISTS "${pathtemp}" AND NOT IS_DIRECTORY "${pathtemp}")
                    set(paramonte_${collection}BR_file "${pathtemp}")
                endif()
            endif()
        endif()

        #### Create the collection items list.

        unset(paramonte_${collection}BR_dir_list)
        if (DEFINED paramonte_${collection}BR_file)
            message(NOTICE "${pmattn} The specified ${collection} file list detected at: ${paramonte_${collection}BR_file}")
            file(STRINGS "${paramonte_${collection}BR_file}" paramonte_${collection}BR_dir_list)
        else()
            string(REPLACE "." ";" ${collection}List "${${collection}BR}")
            string(REPLACE "," ";" ${collection}List "${${collection}List}")
            string(REPLACE " " ";" ${collection}List "${${collection}List}")
            message(NOTICE "${pmattn} Enabling the requested paramonte::${lang} library ${collection} as outlined by the flag value \"${${collection}BR}\"")
            setSubDirList(paramonte_collection_mod_list "${destin}")
            if(lang_is_dynamic)
                file(GLOB_RECURSE paramonte_${collection}BR_dir_list
                    LIST_DIRECTORIES true
                    #RELATIVE "${destin}"
                    *main${lang_ext}
                    )
            else()
                foreach(paramonte_collection_inc_dir ${paramonte_collection_mod_list})
                    set(paramonte_collection_inc_dir_abs "${destin}/${paramonte_collection_inc_dir}")
                    setSubDirList(paramonte_collection_sub_dir_list "${paramonte_collection_inc_dir_abs}")
                    foreach(paramonte_collection_sub_dir ${paramonte_collection_sub_dir_list})
                        set(paramonte_collection_current ${paramonte_collection_inc_dir}/${paramonte_collection_sub_dir})
                        foreach(bitem ${${collection}List}) # no quotations for ${collection}List!
                            #message(NOTICE "${pmwarn} ${paramonte_collection_inc_dir} ${paramonte_collection_sub_dir} ${bitem}")
                            string(TOLOWER "${bitem}" bitem_lower)
                            if ("${bitem_lower}" STREQUAL "all"
                                OR
                                "${paramonte_collection_inc_dir}" MATCHES ".*${bitem}.*" # OR "${bitem}" MATCHES ".*${paramonte_collection_inc_dir}.*"
                                OR
                                "${paramonte_collection_sub_dir}" MATCHES ".*${bitem}.*" # OR "${bitem}" MATCHES ".*${paramonte_collection_sub_dir}.*"
                                )
                                list(APPEND paramonte_${collection}BR_dir_list ${paramonte_collection_current})
                            endif()
                        endforeach()
                    endforeach()
                endforeach()
            endif()
        endif()

        #### Loop over all collection items to add the custom commands.

        if (DEFINED paramonte_${collection}BR_dir_list)
            #set(counter 1)
            foreach(paramonte_${collection}BR_dir_rel ${paramonte_${collection}BR_dir_list})
                set(paramonte_${collection}BR_dir_abs "${destin}/${paramonte_${collection}BR_dir_rel}")
                set(paramonte_collection_build_script_path "${paramonte_${collection}BR_dir_abs}/${paramonte_collection_build_script}")
                if (EXISTS "${paramonte_collection_build_script_path}" AND NOT IS_DIRECTORY "${paramonte_collection_build_script_path}")
                    message(NOTICE "${pmattn} Enabling the ${collection} build at: ${paramonte_${collection}BR_dir_abs}")
                    file(REMOVE_RECURSE "${paramonte_${collection}BR_dir_abs}/*.exe")
                    file(REMOVE_RECURSE "${paramonte_${collection}BR_dir_abs}/*.out")
                    set(mod_proc_name "${paramonte_${collection}BR_dir_rel}${collection}")
                    string(REPLACE "/" "_" mod_proc_name "${mod_proc_name}")
                    string(REPLACE "." "_" mod_proc_name "${mod_proc_name}")
                    string(REPLACE "-" "_" mod_proc_name "${mod_proc_name}")
                    # ensure the target name is unique.
                    # This target counting is crucial when the user-specified collection keywords overlap.
                    if (TARGET "${mod_proc_name}BR")
                        message(NOTICE "${pmwarn} Target \"${mod_proc_name}BR\" has been already added. skipping...")
                    else()
                        add_custom_target(  "${mod_proc_name}BR" #ALL
                                            WORKING_DIRECTORY "${paramonte_${collection}BR_dir_abs}"
                                            COMMAND "${paramonte_collection_build_script_path}"
                                            COMMENT "${pmattn} Building and running the ${collection} scripts: ${BoldCyan} ${paramonte_collection_build_script_path} ${ColorReset}"
                                            USES_TERMINAL
                                            )
                        add_dependencies(${collection} "${mod_proc_name}BR")
                    endif()
                    #math(EXPR counter "${counter}+1")
                endif()
            endforeach()
        else()
            printUsage()
            message(FATAL_ERROR
                    "\n"
                    "${pmfatal} Unsupported user-specified ParaMonte library ${collection} detected. ${collection}BR=${${collection}BR}\n"
                    "${pmfatal} Follow the guidelines above to appropriately specify the testing mode or drop the option.\n"
                    "\n"
                    )
        endif()

    endif()

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Build and run the ParaMonte ${collection} postprocessing (PP).
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if(${collection}PPEnabled)

        if(Python3_Interpreter_FOUND)

            message(NOTICE "${pmattn} Enabling the selected paramonte::${lang} library ${collection} builds...")

            unset(paramonte_${collection}PP_dir_list)

            if ("${${collection}PP}" STREQUAL "")
                if (DEFINED paramonte_${collection}BR_dir_list)
                    set(paramonte_${collection}PP_dir_list "${paramonte_${collection}BR_dir_list}")
                endif()
            endif()

            if (NOT DEFINED paramonte_${collection}PP_dir_list)

                # Check if an input list of ${collection}PP folders is specified.

                unset(paramonte_${collection}PP_file)
                # Assume an absolute path specified.
                set(pathtemp "${${collection}PP}")
                if (EXISTS "${pathtemp}" AND NOT IS_DIRECTORY "${pathtemp}")
                    set(paramonte_${collection}PP_file "${${collection}PP}")
                else()
                    # Assume a relative path with respect to root directory.
                    set(pathtemp "${paramonte_dir}/${${collection}PP}")
                    if (EXISTS "${pathtemp}" AND NOT IS_DIRECTORY "${pathtemp}")
                        set(paramonte_${collection}PP_file "${pathtemp}")
                    else()
                        # Assume a relative path with respect to ${collection}/${lang} root directory.
                        set(pathtemp "${paramonte_dir}/${collection}/${lang}/${${collection}PP}")
                        if (EXISTS "${pathtemp}" AND NOT IS_DIRECTORY "${pathtemp}")
                            set(paramonte_${collection}PP_file "${pathtemp}")
                        endif()
                    endif()
                endif()

                unset(paramonte_${collection}PP_file)
                if (EXISTS "${${collection}PP}")
                    set(paramonte_${collection}PP_file "${${collection}PP}")
                else()
                    set(paramonte_${collection}PP_file "${paramonte_dir}/${${collection}PP}")
                    #file(TO_NATIVE_PATH "${paramonte_${collection}PP_file}" paramonte_${collection}PP_file)
                    if (NOT EXISTS "${paramonte_${collection}PP_file}")
                        unset(paramonte_${collection}PP_file)
                    endif()
                endif()

                # If no input list of ${collection}PP folders is specified, create a list of ${collection}PP directories.

                if (DEFINED paramonte_${collection}PP_file)
                    message(NOTICE "${pmattn} The specified ${collection} postprocessing file list detected at: ${paramonte_${collection}PP_file}")
                    file(STRINGS "${paramonte_${collection}PP_file}" paramonte_${collection}PP_dir_list)
                else() #if (NOT "${${collection}PP_lower}" STREQUAL "none")
                    # Build a list of ${collection} postprocessing.
                    string(REPLACE "." ";" ${collection}PPList "${${collection}PP}")
                    string(REPLACE "," ";" ${collection}PPList "${${collection}PPList}")
                    string(REPLACE " " ";" ${collection}PPList "${${collection}PPList}")
                    message(NOTICE "${pmattn} Enabling the requested paramonte::${lang} library ${collection} postprocessing as outlined by the flag value ${${collection}PP}")
                    setSubDirList(paramonte_${collection}PP_mod_list "${destin}")
                    foreach(paramonte_${collection}PP_inc_dir ${paramonte_${collection}PP_mod_list})
                        set(paramonte_${collection}PP_inc_dir_abs "${destin}/${paramonte_${collection}PP_inc_dir}")
                        setSubDirList(paramonte_${collection}PP_sub_dir_list "${paramonte_${collection}PP_inc_dir_abs}")
                        foreach(paramonte_${collection}PP_sub_dir ${paramonte_${collection}PP_sub_dir_list})
                            set(paramonte_${collection}PP_current ${paramonte_${collection}PP_inc_dir}/${paramonte_${collection}PP_sub_dir})
                            foreach(bitem ${${collection}PPList}) # no quotations for ${collection}PPList!
                                string(TOLOWER "${bitem}" bitem_lower)
                                if ("${bitem_lower}" STREQUAL "all"
                                    OR
                                    "${paramonte_${collection}PP_inc_dir}" MATCHES ".*${bitem}.*" #OR "${bitem}" MATCHES ".*${paramonte_${collection}PP_inc_dir}.*"
                                    OR
                                    "${paramonte_${collection}PP_sub_dir}" MATCHES ".*${bitem}.*" #OR "${bitem}" MATCHES ".*${paramonte_${collection}PP_sub_dir}.*"
                                    )
                                    list(APPEND paramonte_${collection}PP_dir_list ${paramonte_${collection}PP_current})
                                endif()
                            endforeach()
                        endforeach()
                    endforeach()
                endif()

            endif()

            # Run the ${collection}PP list.

            if (DEFINED paramonte_${collection}PP_dir_list)
                #set(counter 1)
                foreach(paramonte_${collection}PP_dir_rel ${paramonte_${collection}PP_dir_list})
                    set(paramonte_${collection}PP_dir_abs "${destin}/${paramonte_${collection}PP_dir_rel}")
                    set(paramonte_collection_build_script_path "${paramonte_${collection}PP_dir_abs}/${paramonte_collection_build_script}")
                    if (EXISTS "${paramonte_${collection}PP_dir_abs}/main.py")
                        file(GLOB old_image_files LIST_DIRECTORIES TRUE "${paramonte_${collection}PP_dir_abs}/*.png")
                        foreach(old_image_file ${old_image_files})
                            file(REMOVE "${old_image_file}")
                        endforeach()
                        set(mod_proc_name "${paramonte_${collection}PP_dir_rel}${collection}")
                        string(REPLACE "/" "_" mod_proc_name "${mod_proc_name}")
                        string(REPLACE "." "_" mod_proc_name "${mod_proc_name}")
                        string(REPLACE "-" "_" mod_proc_name "${mod_proc_name}")
                        if (TARGET "${mod_proc_name}BR")
                       #if (TARGET ${collection}BR${counter})
                           #add_custom_command( TARGET ${collection}BR${counter} POST_BUILD
                            add_custom_command( TARGET "${mod_proc_name}BR" POST_BUILD
                                                WORKING_DIRECTORY "${paramonte_${collection}PP_dir_abs}"
                                                COMMAND "${Python3_EXECUTABLE}" "${paramonte_${collection}PP_dir_abs}/main.py"
                                                COMMENT "${pmattn} Running the ${collection} post-processing scripts: ${Cyan} ${paramonte_${collection}PP_dir_abs}/main.py ${ColorReset}"
                                                USES_TERMINAL
                                                )
                           ##add_custom_target(  ${collection}PP${counter} #ALL
                           # add_custom_target(  "${mod_proc_name}PP" #ALL
                           #                     WORKING_DIRECTORY "${paramonte_${collection}PP_dir_abs}"
                           #                     COMMAND "${Python3_EXECUTABLE}" "${paramonte_${collection}PP_dir_abs}/main.py"
                           #                     COMMENT "${pmattn} Running the ${collection} post-processing scripts: ${Cyan} ${paramonte_${collection}PP_dir_abs}/main.py ${ColorReset}"
                           #                     USES_TERMINAL
                           #                     )
                           # add_dependencies("${mod_proc_name}PP" "${mod_proc_name}BR")
                           ##add_dependencies(${collection}PP${counter} ${collection}BR${counter})
                        elseif (TARGET "${mod_proc_name}PP")
                            message(NOTICE "${pmwarn} Target \"${mod_proc_name}PP\" has been already added. skipping...")
                        else()
                           #add_custom_target(  ${collection}PP${counter} #ALL
                            add_custom_target(  "${mod_proc_name}PP" #ALL
                                                WORKING_DIRECTORY "${paramonte_${collection}PP_dir_abs}"
                                                COMMAND "${Python3_EXECUTABLE}" "${paramonte_${collection}PP_dir_abs}/main.py"
                                                COMMENT "${pmattn} Running the ${collection} post-processing scripts: ${Cyan} ${paramonte_${collection}PP_dir_abs}/main.py ${ColorReset}"
                                                USES_TERMINAL
                                                )
                            add_dependencies(${collection} "${mod_proc_name}PP")
                           #add_dependencies(${collection} ${collection}PP${counter})
                        endif()
                    endif()
                    #math(EXPR counter "${counter}+1")
                endforeach()
            else()
                printUsage()
                message(FATAL_ERROR
                        "\n"
                        "${pmfatal} Unsupported user-specified ParaMonte library ${collection} detected. ${collection}PP=${${collection}PP}\n"
                        "${pmfatal} Follow the guidelines above to appropriately specify the testing mode or drop the option.\n"
                        "\n"
                        )
            endif()

        else()

            message(WARNING
                    "${pmwarn} The ParaMonte::${lang} ${collection} postprocessing could not be set up\n"
                    "${pmwarn} because a Python3 interpreter appears to be missing on your system.\n"
                    "${pmwarn} Skipping the ParaMonte::${lang} library ${collection} postprocessing...\n"
                    )

        endif()

    endif()

else()

    message(WARNING
            "\n"
            "${pmwarn} The ParaMonte::${lang} ${collection} directory in the root directory of the ParaMonte library ${collection} is missing.\n"
            "${pmwarn} paramonte_collection_${lang}_dir=${origin}\n"
            "${pmwarn} Skipping the ParaMonte::${lang} library ${collection} build...\n"
            "\n"
            )

endif()