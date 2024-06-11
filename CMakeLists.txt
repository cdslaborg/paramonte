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

#   CMake project file for the ParaMonte library:
#
#   WARNING:
#   ALL MESSAGE LINES MUST BEGIN WITH A WHITESPACE CHARACTER in CMake.
#   This is essential for proper printing of all message.
#
#   NOTE:
#   A variable reference has the form ${variable_name} and is evaluated inside a Quoted Argument or an Unquoted Argument.
#   A variable reference is replaced by the value of the variable, or by the empty string if the variable is not set.
#
#   For usage guidelines, see the contents of `CMakeLists.md` file in the root directory of the project.
#
#   Commands to build ParaMonte on TACC:
#
#       idev -N 1 -n 68 -m 15 # on KNL
#       idev -p skx-dev -N 1 -n 48 -t 00:15:00 # on SKX
#       cmake .. -DMPI_ENABLED=TRUE -DCFI_ENABLED=TRUE -DDLL_ENABLED=TRUE -DCMAKE_BUILD_TYPE=debug
#       cmake ../../.. -Dbuild=release -Dpar=cafdist
#       make
#       make install
#       make deploy
#       make test
#       make example
#       make benchmark


cmake_minimum_required(VERSION 3.18)

#cmake_policy(VERSION 3.10...3.23)
cmake_policy(VERSION 3.10...3.28)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Add our local modules to the module path
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
#unset(cmakeFileList)
#file(GLOB cmakeFiles "${CMAKE_MODULE_PATH}/string/*.cmake")
#set(cmakeFileList "${cmakeFileList}" "${cmakeFiles}")
#file(GLOB cmakeFiles "${CMAKE_MODULE_PATH}/string/encoding/*.cmake")
#set(cmakeFileList "${cmakeFileList}" "${cmakeFiles}")
#foreach(cmakeFile ${cmakeFileList})
#    message("INCLUDE ${cmakeFile}")
#    include(${cmakeFile})
#endforeach(cmakeFile)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Define CMake colors
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)

    string(ASCII 27 ESC)
    set(ColorReset  "${ESC}[0m")
    set(ColorBold   "${ESC}[1m")
    set(Red         "${ESC}[31m")
    set(Green       "${ESC}[32m")
    set(Yellow      "${ESC}[33m")
    set(Blue        "${ESC}[34m")
    set(Magenta     "${ESC}[35m")
    set(Cyan        "${ESC}[36m")
    set(White       "${ESC}[37m")
    set(BoldRed     "${ESC}[1;31m")
    set(BoldGreen   "${ESC}[1;32m")
    set(BoldYellow  "${ESC}[1;33m")
    set(BoldBlue    "${ESC}[1;34m")
    set(BoldMagenta "${ESC}[1;35m")
    set(BoldCyan    "${ESC}[1;36m")
    set(BoldWhite   "${ESC}[1;37m")

    set(pmcolor "${BoldCyan}")
    set(pmnote  " -- ${pmcolor}ParaMonte CMake -${ColorReset}") # This should not be really used anywhere, except in the following.
    set(pmattn  "${pmnote} ${BoldYellow}NOTE:${ColorReset}")
    set(pmwarn  "${pmnote} ${BoldMagenta}WARNING:${ColorReset}")
    set(pmfatal "${pmnote} ${BoldRed}FATAL ERROR:${ColorReset}")
    set(passed  "${BoldGreen}PASSED${ColorReset}")
    set(failed  "${BoldRed}FAILED${ColorReset}")
    set(found   "${BoldGreen}FOUND${ColorReset}")
    set(missing "${BoldMagenta}MISSING${ColorReset}")

    function(printUsage)
        file(READ "${CMAKE_SOURCE_DIR}/CMakeLists.md" paramonte_cmake_usage)
        message("\n${pmnote} CMake Usage Guidelines:${ColorReset}\n")
        message("${paramonte_cmake_usage}")
    endfunction()

    #function(message)
    #   list(GET ARGV 0 MessageType)
    #   if(MessageType STREQUAL FATAL_ERROR OR MessageType STREQUAL SEND_ERROR)
    #       list(REMOVE_AT ARGV 0)
    #       _message(${MessageType} "${BoldRed}${ARGV}${ColourReset}")
    #   elseif(MessageType STREQUAL WARNING)
    #       list(REMOVE_AT ARGV 0)
    #       _message(${MessageType} "${BoldYellow}${ARGV}${ColourReset}")
    #   elseif(MessageType STREQUAL AUTHOR_WARNING)
    #       list(REMOVE_AT ARGV 0)
    #       _message(${MessageType} "${BoldCyan}${ARGV}${ColourReset}")
    #   elseif(MessageType STREQUAL STATUS)
    #       list(REMOVE_AT ARGV 0)
    #       _message(${MessageType} "${Green}${ARGV}${ColourReset}")
    #   else()
    #       _message("${ARGV}")
    #   endif()
    #endfunction()

    macro(setSubDirList result current_dir)
        file(GLOB children RELATIVE ${current_dir} ${current_dir}/*)
        set(dirlist "")
        foreach(child ${children})
            if(IS_DIRECTORY ${current_dir}/${child})
                list(APPEND dirlist ${child})
            endif()
        endforeach()
        set(${result} ${dirlist})
    endmacro()

    # Assign logical variable to the result of the input expression.
    macro(setBool var)
         if(${ARGN})
             set(${var} ON)
         else()
             set(${var} OFF)
         endif()
    endmacro()

    #   Convert file contents into a CMake list (where each element in the list is one line of the file).
    #
    #   WARNING
    #
    #       there are two caveats wherein this macro does not work as intended:
    #
    #       (1) It puts each line of the (presumed to be text) file into an element of a
    #       cmake list including empty lines as unique elements... However, if you are
    #       familiar with cmake lists at all, you'll instantly say to yourself -- "ah
    #       ha, so blank lines will be skipped in a FOREACH loop." True : caveat #1 :
    #       blank lines are undetectable in a cmake FOREACH structure. But the non-blank
    #       lines are usually the important ones... so skipping the blank ones is
    #       probably ok... depends on your task, I suppose.
    #
    #       (2) It does not handle trailing backslashes at the end of a line
    #       "correctly." It will end up with consecutive lines with trailing backslashes
    #       all put together as one line in the cmake list. caveat #2 : you'll need to
    #       do more work if you want to detect trailing backslashes and get lines
    #       including those stuffed into a cmake list properly
    #
    #       Other than that, it's quite useful and works with CMake 2.4.6 (probably
    #       earlier versions, too, but most of my experience using this technique has
    #       been with 2.4.6.)
    macro(setLineList filepath lineList)
      file(READ "${filepath}" ${lineList})
      string(REGEX REPLACE ";" "\\\\;" ${lineList} "${${lineList}}")
      string(REGEX REPLACE "\n" ";" ${lineList} "${${lineList}}")
    endmacro()

endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set CMake verbosity level.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (NOT DEFINED VERBOSE)
    set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "CMake Verbosity level." FORCE)
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Print an error message on an attempt to build inside the source directory tree:
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if ("${lib}" STREQUAL "")
    set(lib "shared" CACHE STRING "ParaMonte library type: shared, static" FORCE)
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Print an error message on an attempt to build inside the source directory tree:
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    if ("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
        printUsage()
        if (WIN32)
            set(delfiles "del /F CMakeCache.txt; rmdir /Q /s CMakeFiles")
        else()
            set(delfiles "rm -rf CMakeCache.txt CMakeFiles/")
        endif()
        message(FATAL_ERROR
                "\n"
                "${pmfatal} CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}\n"
                "${pmfatal} CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}\n"
                "\n"
                "${pmfatal} This archive does not support in-source builds.\n"
                "\n"
                "${pmfatal} You must now delete the CMakeCache.txt file and the CMakeFiles/ directory under\n"
                "${pmfatal} in the ParaMonte library's root directory, otherwise you will not be able to configure correctly.\n"
                "${pmfatal} You must now run something like:\n"
                "\n"
                "${pmfatal}     ${delfiles}\n"
                "\n"
                "${pmfatal} Please create a directory wither inside or outside the ParaMonte library source tree, \n"
                "${pmfatal} change the directory to it, and build in this directory, in a manner such as the following,\n"
                "\n"
                )
    endif()
    if(DEFINED ENV{SOURCE_DATE_EPOCH})
        set(SOURCE_DATE_EPOCH "$ENV{SOURCE_DATE_EPOCH}")
    endif()
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Display ParaMonte Banner
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

file(READ "${CMAKE_CURRENT_LIST_DIR}/src/auxil/.paramonte.banner" paramonte_banner)
message("\n${pmcolor}${paramonte_banner}${ColorReset}\n")
unset(paramonte_banner)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: set the dependency packaging requirement. This is not used anymore as of version 2.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set(DEPS_ENABLED 0)
if(DEFINED deps)
    message(NOTICE "${pmattn} deps=${deps}")
    string(TOLOWER "${deps}" deps)
    if("${deps}" STREQUAL "all")
        set(DEPS_ENABLED 1)
    endif()
endif()
message(NOTICE "${pmattn} DEPS_ENABLED=${DEPS_ENABLED}")
unset(deps)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: set the required root directories.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

message(NOTICE "${pmattn} setting up the ParaMonte directory hierarchy...")
message(NOTICE "${pmattn} project root directory: ${CMAKE_SOURCE_DIR}")

# define paramonte_dir as CMAKE_SOURCE_DIR

set(paramonte_dir "${CMAKE_SOURCE_DIR}")
set(paramonte_ext_dir "${paramonte_dir}/external")
set(paramonte_ext_obl_dir "${paramonte_dir}/external/OpenBLAS")
set(paramonte_example_dir "${paramonte_dir}/example" CACHE PATH "ParaMonte example source files directory")
set(paramonte_benchmark_dir "${paramonte_dir}/benchmark" CACHE PATH "ParaMonte benchmark source files directory")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set the ParaMonte source files directory
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set(paramonte_src_dir "${CMAKE_SOURCE_DIR}/src" CACHE PATH "ParaMonte interface directory")#/interface

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Set the ParaMonte other-language source files directory and loop over them to ensure their existence.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    set(list_lang c cpp fortran go java julia mathematica matlab python r)
    if (EXISTS "${paramonte_src_dir}/")
        message(NOTICE "${pmattn} Interface source files directory: ${paramonte_src_dir}")
        foreach(item ${list_lang}) # WARNING: lang cannot be used as loop variable because it is already used as an input variable to CMake.
            set(paramonte_src_${item}_dir "${paramonte_src_dir}/${item}" CACHE PATH "ParaMonte ${item} source files directory.")
            if (EXISTS "${paramonte_src_${item}_dir}/")
                message(NOTICE "${pmattn} Checking the existence of ${paramonte_src_${item}_dir} ... ${found}")
                set(paramonte_src_${item}_main_dir "${paramonte_src_${item}_dir}/main" CACHE PATH "ParaMonte ${item} main source files directory.")
                set(paramonte_src_${item}_test_dir "${paramonte_src_${item}_dir}/test" CACHE PATH "ParaMonte ${item} test source files directory.")
                set(paramonte_src_${item}_test_input_dir "${paramonte_src_${item}_test_dir}/input" CACHE PATH "ParaMonte ${item} test input source files directory.")
            else()
                message(NOTICE "${pmattn} Checking the existence of ${paramonte_src_${item}_dir} ... ${missing}")
            endif()
            set(collections "benchmark" "example")
            foreach(collection in ${collections})
                set(paramonte_${collection}_${item}_dir "${paramonte_${collection}_dir}/${item}" CACHE PATH "ParaMonte ${item} ${collection} files directory.")
                if (EXISTS "${paramonte_${collection}_${item}_dir}/")
                    message(NOTICE "${pmattn} Checking the existence of ${paramonte_${collection}_${item}_dir} ... ${found}")
                else()
                    message(NOTICE "${pmattn} Checking the existence of ${paramonte_${collection}_${item}_dir} ... ${missing}")
                endif()
            endforeach()
        endforeach()
    else()
        message(FATAL_ERROR
                "\n"
                "${pmfatal} The ParaMonte interface source directories does not exist.\n"
                "${pmfatal} The integrity of the ParaMonte library is compromised.\n"
                "${pmfatal} Clone the ParaMonte library from GitHub and retry:\n"
                "\n"
                "   https://github.com/cdslaborg/paramonte\n"
                "\n"
                )
    endif()
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Set and generate the required build and binary installation directories.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)

    if (DEFINED fresh)
        string(TOLOWER "${fresh}" fresh)
    endif()
    message(NOTICE "${pmattn} Setting up the ParaMonte build directories...")

    #if (DEFINED bdir)
    #    message(NOTICE "${pmattn} User-specified ParaMonte build directory detected. bdir=\"${bdir}\"")
    #else()
    #    set(bdir "${CMAKE_SOURCE_DIR}/bld/${os}/${arch}/${csid}/${CMAKE_Fortran_COMPILER_VERSION}/${build}/${lib}/${mem}/${par}/${lang}")
    #    message(NOTICE "${pmattn} Defining and generating the ParaMonte build directory: bdir=\"${bdir}\"")
    #    set(ENV{bdir} "${bdir}")
    #    #set(paramonte_bld_dir "${CMAKE_BINARY_DIR}/${csid}/${CMAKE_Fortran_COMPILER_VERSION}/${CMAKE_BUILD_TYPE}/${LTYPE}/${PARALLELISM_DIR}")
    #    #set(paramonte_bld_dir "${CMAKE_BINARY_DIR}")
    #endif()

    set(bdir "${CMAKE_BINARY_DIR}")
    set(paramonte_bld_dir "${bdir}")
    set(ENV{paramonte_bld_dir} "${paramonte_bld_dir}")
    set(CMAKE_INSTALL_PREFIX "${paramonte_bld_dir}" CACHE BOOL "ParaMonte/OpenBLAS/dependencies install directory prefix." FORCE)
    message(NOTICE "${pmattn} All build files will be stored at paramonte_bld_dir=\"${paramonte_bld_dir}\"")

    #### Purge files or folders within the build folder if requested.

    file(GLOB bld_dir_items LIST_DIRECTORIES true RELATIVE "${paramonte_bld_dir}" "${paramonte_bld_dir}/*")
    message(NOTICE "${pmattn} The current contents of the ParaMonte build directory: \"${bld_dir_items}\"")
    message(NOTICE "${pmattn} Purging any requested contents prior to the ParaMonte library build...")
    if (NOT "${fresh}" STREQUAL "none" AND NOT "${fresh}" STREQUAL "")
        foreach(fitem ${fresh})
            set (purge_failed 1)
            if ("${fitem}" MATCHES ".*[\.][\.].*")
                message(NOTICE "${pmwarn} Purging parent-level or higher folders and objects are impossible: \"${fitem}\"")
            elseif ("${fitem}" STREQUAL "bdir" OR "${fitem}" STREQUAL "all")
                message(NOTICE "${pmattn} Purging the entire ParaMonte library build directory: \"${paramonte_bld_dir}/\"")
                foreach(bld_dir_item ${bld_dir_items})
                    file(REMOVE_RECURSE "${paramonte_bld_dir}/${bld_dir_item}")
                endforeach()
                set (purge_failed 0)
            elseif (NOT "${fitem}" STREQUAL "ddir") # This will be taken care of later.
                unset(dirsep)
                if ("${fitem}" MATCHES ".*[/].*")
                    set(dirsep "/")
                elseif ("${fitem}" MATCHES ".*[\\].*")
                    set(dirsep "\\")
                endif()
                if (DEFINED dirsep)
                    list(GET fitems 0 fitems1)
                    list(GET fitems 1 fitems2)
                    if(IS_DIRECTORY "${paramonte_bld_dir}/${fitems1}")
                        set(purge_path "${paramonte_bld_dir}/${fitems1}/${fitems2}")
                        message(NOTICE "${pmattn} Purging the build sub-subdirectory: \"${purge_path}/\"")
                        file(REMOVE_RECURSE "${purge_path}")
                        set (purge_failed 0)
                    endif()
                else()
                    set(purge_path "${paramonte_bld_dir}/${fitem}")
                    message(NOTICE "${pmattn} Purging the build subdirectory: \"${purge_path}/\"")
                    file(REMOVE_RECURSE "${purge_path}")
                    set (purge_failed 0)
                endif()
            endif()
            if (purge_failed)
                message(NOTICE "${pmwarn} Nonexistent or unrecognized pattern to remove from the current build: \"${fitem}\"")
            endif()
        endforeach()
    endif()

    set(subdir_list benchmark example inc lib obj pkg)
    foreach(subdir ${subdir_list})
        set(paramonte_bld_${subdir}_dir "${paramonte_bld_dir}/${subdir}")
        message(NOTICE "${pmattn} paramonte_bld_${subdir}_dir: ${paramonte_bld_${subdir}_dir}")
        if (NOT EXISTS "${paramonte_bld_${subdir}_dir}/")
            file(MAKE_DIRECTORY "${paramonte_bld_${subdir}_dir}")
        elseif ("${fresh}" MATCHES ".*${subdir}.*")
            file(REMOVE_RECURSE "${paramonte_bld_${subdir}_dir}/")
            file(MAKE_DIRECTORY "${paramonte_bld_${subdir}_dir}")
        endif()
    endforeach()

    # paramonte_bld_pkg_dir must be created before calling the library target file.
    set(paramonte_bld_pkg_benchmark_dir "${paramonte_bld_pkg_dir}/benchmark")
    set(paramonte_bld_pkg_example_dir "${paramonte_bld_pkg_dir}/example")
    set(paramonte_bld_pkg_fpp_dir "${paramonte_bld_pkg_dir}/fpp")
    set(paramonte_bld_pkg_lib_dir "${paramonte_bld_pkg_dir}/lib")
    set(paramonte_bld_pkg_inc_dir "${paramonte_bld_pkg_dir}/inc")

endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Set the default performance profiling option.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set(perfprof_enabled FALSE)
if ("${perfprof}" MATCHES "[Aa][Ll][Ll]")
    set(perfprof_enabled TRUE)
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Set the code coverage preprocessor flag.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (DEFINED test)
    string(TOLOWER "${test}" test)
endif()

if ("${test}" STREQUAL "all")
    if (NOT "${mod}" STREQUAL "" AND NOT "${mod}" STREQUAL "all")
        message(WARNING
                "\n"
                "${pmwarn} Specifying `all` for the value of build configuration option `test`\n"
                "${pmwarn} requires setting the build configuration option `mod` also to `all`.\n"
                "${pmwarn} You have specified `-Dmod=${mod}`.\n"
                "${pmwarn} Resetting `mod` to `all`...\n"
                "\n"
                )
        set(mod "all" CACHE "modules to compiler." FORCE)
    endif()
endif()

set(codecov_enabled FALSE)
if("${codecov}" MATCHES "[Aa][Ll][Ll]")
    add_compile_definitions("CODECOV_ENABLED=1")
    set(codecov_enabled TRUE)
    if (NOT DEFINED test)# OR "${test}" STREQUAL "none"
        set(test "all")
    endif()
    if (NOT DEFINED checking)
        set(checking "all")
    endif()
endif()

if (DEFINED fpp)
    string(TOLOWER "${fpp}" fpp)
endif()

if (codecov_enabled)
    if ("${fpp}" STREQUAL "all")
        if (NOT "${fpp}" STREQUAL "" AND NOT "${fpp}" STREQUAL "all")
            message(WARNING
                    "\n"
                    "${pmwarn} Specifying `all` for the value of build configuration option `test`\n"
                    "${pmwarn} requires setting the build configuration option `fpp` also to `all`.\n"
                    "${pmwarn} You have specified `-Dfpp=${fpp}`.\n"
                    "${pmwarn} Resetting `fpp` to `all`...\n"
                    "\n"
                    )
            set(fpp "all" CACHE "Fortran preprocessing rule." FORCE)
        endif()
    endif()
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Predefine the build mode. This must be done before calling project() or enable_language().
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

include(presetParaMonteBuild) # Predefines `CMAKE_BUILD_TYPE` and `CMAKE_CONFIGURATION_TYPES` based on the user choice `build`.

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Preset the desired compiler suite. This must be done before calling project() or enable_language().
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

include(presetParaMonteCompiler) # Predefines `CMAKE_Fortran_COMPILER` based on the user choice `fc`.

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Specify source languages and standards.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/.git")
        message(NOTICE "${pmattn} Build from git repository detected.")
    endif()
    message(NOTICE "${pmattn} CMake executable path: ${CMAKE_COMMAND}")
    message(NOTICE "${pmattn} ParaMonte root directory: ${CMAKE_CURRENT_SOURCE_DIR}")
    #   Languages are enabled by the project() command.
    #   Language-specific built-in variables, such as CMAKE_CXX_COMPILER, CMAKE_CXX_COMPILER_ID etc are set by invoking `project()`.
    #   An implicit `project()` command will be automatically generated if no project command is in the top-level CMakeLists file.
    set(CMAKE_Fortran_STANDARD 2008)
    set(CMAKE_SHARED_LIBRARY_PREFIX_Fortran "") # The default is `lib`.
    set(CMAKE_STATIC_LIBRARY_PREFIX_Fortran "") # The default is `lib`.
    enable_language(Fortran) # Fortran C CXX # redundant
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Reset the compiler suite name variable `CMAKE_Fortran_COMPILER_ID` if undefined. This must be done after project initiation.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)

    if (NOT DEFINED CMAKE_Fortran_COMPILER_ID OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "")
        #set(csid "${CMAKE_Fortran_COMPILER_ID}" CACHE STRING "CMAKE_Fortran_COMPILER_ID suite" FORCE)
        message(WARNING
                "\n"
                "${pmwarn} The CMake Compiler Suite variable is undefined. CMAKE_Fortran_COMPILER_ID=${CMAKE_Fortran_COMPILER_ID}\n"
                "${pmwarn} This name is required for naming the output binary. skipping...\n"
                "\n"
                )
        if ("${CMAKE_Fortran_COMPILER}" MATCHES "ifort" OR "${CMAKE_Fortran_COMPILER_ID}" MATCHES "intel")
            set(CMAKE_Fortran_COMPILER_ID "Intel")
        elseif("${CMAKE_Fortran_COMPILER}" MATCHES "gfortran" OR "${CMAKE_Fortran_COMPILER_ID}" MATCHES "gnu")
            set(CMAKE_Fortran_COMPILER_ID "GNU")
        endif()
    endif()

    set(csid_is_gnu FALSE)
    set(csid_is_intel FALSE)

    if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "[Ii][Nn][Tt][Ee][Ll]")
        set(csid_is_intel TRUE)
    elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "[Gg][Nn][Uu]")
        set(csid_is_gnu TRUE)
    else()
        message(WARNING
                "\n"
                "${pmwarn} \n"
                "${pmwarn} The compiler choice to build the ParaMonte library is unsupported and untested:\n"
                "${pmwarn}     CMAKE_Fortran_COMPILER_ID=${CMAKE_Fortran_COMPILER_ID}\n"
                "${pmwarn}     CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}\n"
                "${pmwarn} Currently, only the Intel and the GNU compiler suites are officially supported\n"
                "${pmwarn} for building the ParaMonte library. You can specify either of these compiler suites\n"
                "${pmwarn} by providing their names or path to their executables to cmake as arguments\n"
                "${pmwarn} \n"
                "${pmwarn}       -Dfc=ifort\n"
                "${pmwarn} \n"
                "${pmwarn}    or\n"
                "${pmwarn} \n"
                "${pmwarn}        -Dfc=ifort\n"
                "${pmwarn} \n"
                "${pmwarn}    Proceeding without any guarantee of build success. You are now in an uncharted territory.\n"
                "${pmwarn}    Please report any failures to: https://github.com/cdslaborg/paramonte/issues/new/choose\n"
                "\n"
                )
    endif()

    if (${csid_is_gnu} OR ${csid_is_intel})
        #set(CMAKE_C_FLAGS_INIT "")
        #set(CMAKE_C_FLAGS_DEBUG_INIT "")
        #set(CMAKE_C_FLAGS_RELEASE_INIT "")
        #set(CMAKE_CXX_FLAGS_INIT "")
        #set(CMAKE_CXX_FLAGS_DEBUG_INIT "")
        #set(CMAKE_CXX_FLAGS_RELEASE_INIT "")
        set(CMAKE_Fortran_FLAGS "")
        set(CMAKE_Fortran_FLAGS_DEBUG "")
        set(CMAKE_Fortran_FLAGS_RELEASE "")
        set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "")
        set(CMAKE_Fortran_FLAGS_INIT "")
        set(CMAKE_Fortran_FLAGS_DEBUG_INIT "")
        set(CMAKE_Fortran_FLAGS_RELEASE_INIT "")
        set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO_INIT "")
    endif()

    set(csid "${CMAKE_Fortran_COMPILER_ID}")
    string(TOLOWER "${csid}" csid)

    message(NOTICE "${pmattn} CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}")
    message(NOTICE "${pmattn} CMAKE_Fortran_COMPILER_ID=${CMAKE_Fortran_COMPILER_ID}")
    message(NOTICE "${pmattn} CMAKE_Fortran_COMPILER_VERSION=${CMAKE_Fortran_COMPILER_VERSION}")

endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set CMAKE_Fortran_COMPILER_VERSION if CMake doesn't do it for us. This must happen after calling project().
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    #   This script defines the following variables:
    #  `CMAKE_Fortran_COMPILER_VERSION` (if undefined already): the full version of the compiler in `major[.minor[.patch[.tweak]]]` format.
    include(setParaMonteCompilerVersion)
    #   We have populated CMAKE_Fortran_COMPILER_VERSION if it was missing.
    #   Now extract the compiler major version.
    if (
        (${csid_is_gnu} AND (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10.3.0)) OR
        (${csid_is_intel} AND (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 21.0.0))
        )
        message(FATAL_ERROR
                "\n"
                "${pmfatal} Building ParaMonte minimally requires either:\n"
                "${pmfatal} \n"
                "${pmfatal}     1. GNU Fortran compiler version 11.0.0 or newer.\n"
                "${pmfatal}     2. Intel Fortran compiler version 21.0.0 or newer.\n"
                "${pmfatal} \n"
                "${pmfatal} Please install one of the above two comiler suites and rebuild.\n"
                "\n"
                )
    else()
        if (NOT (${csid_is_gnu} OR ${csid_is_intel}))
            message(WARNING
                    "\n"
                    "${pmwarn} Untested Fortran compiler version.\n"
                    "${pmwarn} The ParaMonte library is currently only tested with the following compilers:\n"
                    "${pmwarn} \n"
                    "${pmwarn}      1. GNU Fortran compiler version 11.0.0 or newer.\n"
                    "${pmwarn}      2. Intel Fortran compiler version 21.0.0 or newer\n"
                    "${pmwarn} \n"
                    "${pmwarn} The build will proceed with no guarantee of success...\n"
                    "${pmwarn} \n"
                    "\n"
                    )
        endif()
        string(REPLACE "." ";" COMPILER_VERSION_LIST ${CMAKE_Fortran_COMPILER_VERSION})
        list(GET COMPILER_VERSION_LIST 0 COMPILER_VERSION_MAJOR)
        unset(COMPILER_VERSION_LIST)
        message(NOTICE "${pmattn} COMPILER_VERSION_MAJOR=${COMPILER_VERSION_MAJOR}")
    endif()
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Test modern/Coarray Fortran compiler awareness.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

include(testParaMonteCompiler)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set the processor architecture.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    #   CMAKE_HOST_SYSTEM_PROCESSOR:
    #
    #   On Windows, the value of the PROCESSOR_ARCHITECTURE environment variable is used.
    #   The options are: AMD64, IA64, ARM64, EM64T, X86. Source: this SuperUser answer.
    #
    #   On macOS, the value of uname -m is used by default.
    #   However, since this might vary based on whether you're using x86 or ARM CMake, version 3.19.2+ will use the value of CMAKE_APPLE_SILICON_PROCESSOR instead, if it is set.
    #   It also normalizes Power Macintosh to powerpc. The possible values here are x86_64, arm64, and powerpc.
    #
    #   Note that CMake will use the value of CMAKE_HOST_SYSTEM_PROCESSOR (and CMAKE_SYSTEM_NAME) for a number of tasks (one should not change its value):
    #
    #       -#  The default value of CPACK_SYSTEM_NAME is ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}.
    #       -#  The language modules Modules/CMake<LANG>Information.cmake all optionally include platform modules suffixed with -${CMAKE_SYSTEM_PROCESSOR}.cmake.
    #           In most cases, no such modules exist. On Android, very many platform modules exist with processor-specific variants.
    #       -#  The FindJNI module uses it to guide its search.
    #       -#  It is passed to ARMClang via --mcpu (compile) and --cpu (link).
    #
    if (
        "${CMAKE_HOST_SYSTEM_PROCESSOR}" STREQUAL "x64" OR
        "${CMAKE_HOST_SYSTEM_PROCESSOR}" STREQUAL "amd64" OR
        "${CMAKE_HOST_SYSTEM_PROCESSOR}" STREQUAL "AMD64" OR
        "${CMAKE_HOST_SYSTEM_PROCESSOR}" MATCHES "[Xx]86[_-]64" #   [_-] must have - as the last character.
        )
        set(arch "amd64")
    else()
        set(arch "${CMAKE_HOST_SYSTEM_PROCESSOR}")
    endif()
    string(TOLOWER ${arch} arch)
    message(NOTICE "${pmattn} Building for target architecture: ${arch}")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Silence warnings about dereferencing unset variables.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    if(NOT CMAKE_REQUIRED_FLAGS)
      set(CMAKE_REQUIRED_FLAGS "")
    endif()
    if(NOT CMAKE_REQUIRED_LIBRARIES)
      set(CMAKE_REQUIRED_LIBRARIES "")
    endif()
    if(NOT CMAKE_REQUIRED_INCLUDES)
      set(CMAKE_REQUIRED_INCLUDES "")
    endif()
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Predefine the blas/lapack mode. This must be done before calling project() or enable_language().
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# WARNING
# No generic-target compile definition must be added above this line.
# Though generally harmless, the definitions will leak to OpenBLAS build.
if (TRUE)
    # Set up BLAS / LAPACK requirements.
    set(BLAS_REQUIRED 0)
    if (DEFINED blas AND NOT "${blas}" MATCHES "[nN][oO][nN][eE]" AND NOT "${blas}" MATCHES "[Dd][Ee][Ff][Aa][Uu][Ll][Tt]")
        set(BLAS_REQUIRED 1)
    endif()
    set(LAPACK_REQUIRED 0)
    if (DEFINED lapack AND NOT "${lapack}" MATCHES "[nN][oO][nN][eE]" AND NOT "${lapack}" MATCHES "[Dd][Ee][Ff][Aa][Uu][Ll][Tt]")
        set(LAPACK_REQUIRED 1)
    endif()
    # Build the OpenBLAS source library if it exists and no other blas library is enabled.
    set(BLAS_ENABLED 0)
    set(LAPACK_ENABLED 0)
    set(OpenBLAS_ENABLED 0)

    if (EXISTS "${paramonte_ext_obl_dir}/") # This is source directory of OpenBLAS # "${blas}" MATCHES "[Oo][Pp][eE][nN][Bb][Ll][Aa][Ss]" OR "${lapack}" MATCHES "[Oo][Pp][eE][nN][Bb][Ll][Aa][Ss]"
        if (BLAS_REQUIRED OR LAPACK_REQUIRED)
            message(NOTICE "${pmwarn} A requested BLAS library is already enabled, while OpenBLAS source is detected in the project root directory. Skipping OpenBLAS source build...")
            message(NOTICE "${pmwarn} lapack=${lapack}")
            message(NOTICE "${pmwarn} blas=${blas}")
        elseif ("${blas}" MATCHES "[Dd][Ee][Ff][Aa][Uu][Ll][Tt]") # default to local OpenBLAS project build.
            # Full list of options at: https://github.com/xianyi/OpenBLAS/blob/develop/Makefile.rule
            if (FALSE)
                # Always build OpenBLAS as static position-independent library.
                # This is likely what we want to do for non-compiled languages.
                set(BUILD_STATIC_LIBS 1 CACHE BOOL "OpenBLAS no static library." FORCE)
                set(BUILD_SHARED_LIBS 0 CACHE BOOL "OpenBLAS no shared library." FORCE)
            elseif ("${lib}" MATCHES "[Dd][Yy][Nn][Aa][Mm][Ii][Cc]" OR "${lib}" MATCHES "[sS][hH][aA][rR][eE][dD]") # DLL_ENABLED
                    set(BUILD_STATIC_LIBS 0 CACHE BOOL "OpenBLAS no static library." FORCE)
                    set(BUILD_SHARED_LIBS 1 CACHE BOOL "OpenBLAS no shared library." FORCE)
            else()
                set(BUILD_STATIC_LIBS 1 CACHE BOOL "OpenBLAS no static library." FORCE)
                set(BUILD_SHARED_LIBS 0 CACHE BOOL "OpenBLAS no shared library." FORCE)
            endif()
            set(NO_CBLAS 1 CACHE BOOL "OpenBLAS No CBLAS." FORCE)
            set(NO_LAPACKE 1 CACHE BOOL "OpenBLAS no C-LAPACK." FORCE)
            if (NOT DEFINED DYNAMIC_ARCH)
                set(DYNAMIC_ARCH 0 CACHE BOOL "OpenBLAS multiple dispatch." FORCE)
            endif()
            #set(USE_THREAD ${OMP_ENABLED} CACHE BOOL "OpenBLAS threading." FORCE)
            set(USE_THREAD 1 CACHE BOOL "OpenBLAS threading." FORCE)
            # The following install prefix is now set globally for the entire project.
            #set(CMAKE_INSTALL_PREFIX "${paramonte_bld_dir}" CACHE BOOL "OpenBLAS install directory." FORCE)
            # setup openblas files directory
            set(paramonte_bld_obl_dir "${paramonte_bld_dir}/obl")
            message(NOTICE "${pmattn} paramonte_bld_obl_dir: ${paramonte_bld_obl_dir}")
            if (NOT EXISTS "${paramonte_bld_obl_dir}/")
                file(MAKE_DIRECTORY "${paramonte_bld_obl_dir}")
            elseif ("${fresh}" MATCHES ".*obl.*" OR "${fresh}" MATCHES ".*openblas.*")
                file(REMOVE_RECURSE "${paramonte_bld_obl_dir}/")
            endif()
            add_subdirectory("${paramonte_ext_obl_dir}" "${paramonte_bld_obl_dir}")
            set(OpenBLAS_ENABLED 1)
            set(LAPACK_ENABLED 1)
            set(BLAS_ENABLED 1)
        endif()
    endif()
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set the target language.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    set(C_ENABLED 0)
    set(CPP_ENABLED 0)
    set(FORTRAN_ENABLED 0)
    set(GO_ENABLED 0)
    set(JAVA_ENABLED 0)
    set(JULIA_ENABLED 0)
    set(MATHEMATICA_ENABLED 0)
    set(MATLAB_ENABLED 0)
    set(PYTHON_ENABLED 0)
    set(R_ENABLED 0)
    set(lang_is_dynamic 1)
    string(TOLOWER "${lang}" lang)
    if ("${lang}" STREQUAL "c")
        set(lang_is_dynamic 0)
        set(lang_ext ".c")
        set(C_ENABLED 1)
    elseif ("${lang}" MATCHES "[Cc][XxPp+][XxPp+]")
        set(lang_is_dynamic 0)
        set(lang_ext ".cpp")
        set(CPP_ENABLED 1)
        set(lang "cpp")
    elseif ("${lang}" MATCHES "[Cc]#" OR "${lang}" MATCHES "[Cc][sS][hH][aA][rR][pP]")
        set(lang_is_dynamic 0)
        set(CSHARP_ENABLED 1)
        set(lang_ext ".cs")
        set(lang "csharp")
    elseif ("${lang}" MATCHES "[Ff][Oo][Rr][Tt][Rr][Aa][Nn]")
        set(lang_is_dynamic 0)
        set(FORTRAN_ENABLED 1)
        set(lang_ext ".F90")
    elseif ("${lang}" MATCHES "[Gg][Oo]")
        set(lang_ext ".go")
        set(GO_ENABLED 1)
    elseif ("${lang}" MATCHES "[Jj][Aa][Vv][Aa]")
        set(lang_ext ".java")
        set(JAVA_ENABLED 1)
    elseif ("${lang}" MATCHES "[Jj][Uu][Ll][Ii][Aa]")
        set(JULIA_ENABLED 1)
        set(lang_ext ".jl")
    elseif ("${lang}" MATCHES "[mM][aA][tT][hH][eE][mM][aA][tT][iI][cC][aA]")
        set(MATHEMATICA_ENABLED 1)
        set(lang_ext ".wls")
    elseif ("${lang}" MATCHES "[mM][aA][tT][lL][aA][bB]")
        set(MATLAB_ENABLED 1)
        set(lang_ext ".m")
    elseif ("${lang}" MATCHES "[pP][yY][tT][hH][oO][nN]")
        set(PYTHON_ENABLED 1)
        set(lang_ext ".py")
    elseif ("${lang}" STREQUAL "R" OR "${lang}" STREQUAL "r")
        set(lang_ext ".R")
        set(R_ENABLED 1)
    elseif ("${lang}" MATCHES "[rR][uU][sS][tT]")
        set(RUST_ENABLED 1)
        set(lang_ext ".rs")
    elseif ("${lang}" MATCHES "[sS][aA][sS]")
        set(lang_ext ".sas")
        set(SAS_ENABLED 1)
    elseif ("${lang}" MATCHES "[sS][wW][iI][fF][tT]")
        set(lang_ext ".swift")
        set(SWIFT_ENABLED 1)
    elseif (NOT DEFINED lang OR "${lang}" STREQUAL "")
        set(FORTRAN_ENABLED 1)
        set(lang_ext ".F90")
        set(lang "fortran")
        message(NOTICE
                "${pmwarn} Unspecified target language. You can set the target language via -Dlang=target_lang\n"
                "${pmwarn} where `target_lang` is the name of the target language as written most commonly.\n"
                "${pmwarn} Defaulting to `${lang}` as the target language of choice."
                )
    else()
        printUsage()
        message(FATAL_ERROR
                "\n"
                "${pmfatal} The user-specified target language is unsupported. lang=${lang}\n"
                "${pmfatal} See the above usage notes for the possible choices.\n"
                "\n"
                )
    endif()
    add_compile_definitions("C_ENABLED=${C_ENABLED}")
    add_compile_definitions("CPP_ENABLED=${CPP_ENABLED}")
    add_compile_definitions("FORTRAN_ENABLED=${FORTRAN_ENABLED}")
    add_compile_definitions("GO_ENABLED=${GO_ENABLED}")
    add_compile_definitions("JAVA_ENABLED=${JAVA_ENABLED}")
    add_compile_definitions("JULIA_ENABLED=${JULIA_ENABLED}")
    add_compile_definitions("MATHEMATICA_ENABLED=${MATHEMATICA_ENABLED}")
    add_compile_definitions("MATLAB_ENABLED=${MATLAB_ENABLED}")
    add_compile_definitions("PYTHON_ENABLED=${PYTHON_ENABLED}")
    add_compile_definitions("R_ENABLED=${R_ENABLED}")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Enable or disable the Parameterized Derived Type interface. This is important as GNU gfortran as of 13.2 does not support PDTs.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    if (DEFINED pdt) # anything other than 0, N, NO, or FALSE
        if ("${pdt}" MATCHES "[Aa][Ll][Ll]")
            set(PDT_ENABLED 1)
        else()
            set(PDT_ENABLED 0)
        endif()
    elseif (${csid_is_gnu})
        set(PDT_ENABLED 0)
    elseif (${csid_is_intel})
        set(PDT_ENABLED 0)
    else()
        set(PDT_ENABLED 0)
    endif()
    add_compile_definitions("PDT_ENABLED=${PDT_ENABLED}")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set the C-Fortran interoperability. Enable C-Fortran interoperability (required to call ParaMonte from non-Fortran languages)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    if (DEFINED cfi)
        message(NOTICE "${pmattn} cfi=${cfi}")
        string(TOLOWER ${cfi} cfi)
    endif()
    if ("${cfi}" STREQUAL "all")
        set(CFI_ENABLED 1)
    elseif ("${cfi}" STREQUAL "none") # The only target language can be Fortran.
        if (NOT FORTRAN_ENABLED)
            printUsage()
            message(FATAL_ERROR
                    "\n"
                    "${pmfatal} The user-specified C-Fortran Interoperability (cfi=${cfi})\n"
                    "${pmfatal} evaluates to FALSE which is inconsistent with the target language.\n"
                    "${pmfatal} Set this option according to the guidelines above or drop it.\n"
                    "\n"
                    )
        endif()
        set(CFI_ENABLED 0)
    elseif (DEFINED cfi)
            printUsage()
            message(FATAL_ERROR
                    "\n"
                    "${pmfatal} Unrecognized specified value for the `cfi` C-Fortran Interoperability flag (cfi=${cfi})\n"
                    "${pmfatal} See the above documentation for possible usage.\n"
                    "\n"
                    )
    elseif (FORTRAN_ENABLED) # cfi is undefined.
        set(CFI_ENABLED 0)
    else() # for all other target languages.
        set(CFI_ENABLED 1)
    endif()
    add_compile_definitions("CFI_ENABLED=${CFI_ENABLED}")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Define the library linking type. This must happen after the choice of target language.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    #  shared:      Use this flag when you have R/Python/MATLAB/Julia code to which you need to link the ParaMonte library dynamically, using DLL files.
    #  static:      Use this flag when you have Fortran/C/C++ code to which you want to link the ParaMonte library statically.
    #               You can also link dynamically your Fortran/C/C++ codes using DLL files by specifying LTYPE=shared flag instead.
    unset(DLL_ENABLED) # If not set below, it will be set automatically based on the
    if ("${lib}" MATCHES "[Dd][Yy][Nn][Aa][Mm][Ii][Cc]" OR "${lib}" MATCHES "[sS][hH][aA][rR][eE][dD]")
        set(DLL_ENABLED 1)
    elseif ("${lib}" MATCHES "[sS][tT][aA][tT][iI][cC]")
        if (NOT (C_ENABLED OR CPP_ENABLED OR FORTRAN_ENABLED))
            #printUsage()
            message(NOTICE
                    "${pmwarn} The user-specified library build (lib=${lib})\n"
                    "${pmwarn} can be specified only with C, C++, and Fortran programming languages.\n"
                    "${pmwarn} The library build for all other (dynamic) languages has to be `shared`."
                    )
        endif()
        set(DLL_ENABLED 0)
    elseif (DEFINED lib)
        printUsage()
        message(FATAL_ERROR
                "\n"
                "${pmfatal} The user-specified library build is unsupported. lib=${lib}\n"
                "${pmfatal} See the above usage notes for the possible choices.\n"
                "\n"
                )
    elseif (C_ENABLED OR CPP_ENABLED OR FORTRAN_ENABLED)
        set(DLL_ENABLED 0)
    else()
        set(DLL_ENABLED 1)
    endif()
    add_compile_definitions("DLL_ENABLED=${DLL_ENABLED}")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set runtime checking flag. By this time, `build` must have become lowercase in setParaMonteBuild.cmake
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    if (DEFINED checking)
        string(TOLOWER ${checking} checking)
        if ("${checking}" STREQUAL "nocheck")
            set(CHECK_ENABLED 0)
        elseif ("${checking}" STREQUAL "checked")
            set(CHECK_ENABLED 1)
        else()
            printUsage()
            message(FATAL_ERROR
                    "\n"
                    "${pmfatal} The user-specified `checking` option is unsupported. checking=${checking}\n"
                    "${pmfatal} See the above usage notes for the possible choices.\n"
                    "\n"
                    )
        endif()
    elseif ("${build}" STREQUAL "debug" OR "${build}" STREQUAL "testing")
        set(CHECK_ENABLED 1)
    elseif ("${build}" STREQUAL "release" OR "${build}" STREQUAL "ipo" OR "${build}" STREQUAL "tuned" OR "${build}" STREQUAL "native")
        set(CHECK_ENABLED 0)
    else()
        printUsage()
        message(FATAL_ERROR
                "\n"
                "${pmfatal} Internal CMake file error occurred. Unrecognized `build`. build=${build}\n"
                "${pmfatal} This CMake file is compromised. Download a fresh clone of the ParaMonte library from\n"
                "\n"
                "${pmfatal}     https://github.com/cdslaborg/paramonte\n"
                "\n"
                )
    endif()
    add_compile_definitions("CHECK_ENABLED=${CHECK_ENABLED}")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set runtime procedure purity flag.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    if(DEFINED purity)
        message(NOTICE "${pmattn} purity=${purity}")
        string(TOLOWER "${purity}" purity)
        if("${purity}" STREQUAL "purity")
            add_compile_definitions("PURE=impure")
        elseif("${purity}" STREQUAL "pure")
            add_compile_definitions("PURE=pure")
        else()
            printUsage()
            message(FATAL_ERROR
                    "\n"
                    "${pmfatal} Unrecognized value for `purity` (purity=${purity}).\n"
                    "${pmfatal} See the documentation above for the usage.\n"
                    "\n"
                    )
        endif()
    elseif (${CHECK_ENABLED})# OR "${build}" STREQUAL "debug" OR "${build}" STREQUAL "testing")
        add_compile_definitions("PURE=impure")
    else()
        add_compile_definitions("PURE=pure")
    endif()

    if (${CHECK_ENABLED} AND NOT ${purity})
        printUsage()
        message(FATAL_ERROR
                "\n"
                "${pmfatal} Inconsistent user-specified configuration flags detected `checking`, `purity`: \"${checking}\", \"${purity}\"\n"
                "${pmfatal} The `purity` flag must be set to `impure` when checks are enabled (e.g., when `checking=checked`) \n"
                "${pmfatal} See the documentation above for the usage.\n"
                "\n"
                )
    endif()
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set ParaMonte OS variable and preprocessor macro.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    if(WIN32)
        add_compile_definitions("WINDOWS_ENABLED=1")
        set(os "windows")
    elseif(APPLE)
        add_compile_definitions("DARWIN_ENABLED=1")
        set(os "darwin")
    elseif(UNIX AND NOT APPLE)
        set(os "linux")
        add_compile_definitions("LINUX_ENABLED=1")
        # Check if the os is WSL.
        execute_process(COMMAND sh "-c" "uname -a"
                        RESULT_VARIABLE errorOccurred
                        OUTPUT_VARIABLE infos
                        WORKING_DIRECTORY "."
                        OUTPUT_STRIP_TRAILING_WHITESPACE
                        #COMMAND_ECHO STDOUT
                        )
        if(errorOccurred)
            message(WARNING
                    "${pmwarn} Failed to check for the Microsoft WSL. Here is the error message:\n"
                    "${errorOccurred}\n"
                    "${infos}\n"
                    "${pmwarn} You can ignore this warning if the platform is WSL > 1 or other Linux OS.\n"
                    "${pmwarn} CMake failed to fetch information about the Operating System. This is fine.\n"
                    "${pmwarn} However, if the current platform is Microsoft Windows Subsystem for Linux,\n"
                    "${pmwarn} The preprocessor macro `WSL_ENABLED=1` will not be set.\n"
                    "${pmwarn} This can cause compile errors if the platform is WSL 1.\n"
                    )
        elseif ("${infos}" MATCHES [Mm][Ii][Cc][Rr][Oo][Ss][Oo][Ff][Tt])
            add_compile_definitions("WSL_ENABLED=1")
        endif()
    endif()
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set the memory allocation resource.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    # 0: All automatic and stack arrays will be allocated on the stack. This could lead to stack overflow when calling DLL from other languages.
    # 1: All automatic and stack arrays will be allocated on the heap. Use this when calling ParaMonte from other languages.
    if ("${mem}" STREQUAL "" OR NOT DEFINED mem)
        if ("${lang}" STREQUAL "c" OR "${lang}" STREQUAL "cpp" OR "${lang}" STREQUAL "fortran" AND NOT "${lib}" STREQUAL "shared")
            set(mem "stack")
        else()
            set(mem "heap")
        endif()
    endif()
    string(TOLOWER "${mem}" mem)
    if ("${mem}" STREQUAL "heap") # MATCHES "[Hh][Ee][Aa][Pp]")
        set(HEAP_ARRAY_ENABLED 1)
    elseif ("${mem}" STREQUAL "stack") # MATCHES "[Ss][Tt][Aa][Cc][Kk]")
        set(HEAP_ARRAY_ENABLED 0)
    else()
        printUsage()
        message(FATAL_ERROR
                "\n"
                "${pmfatal} The user-specified library build is invalid. mem=${mem}\n"
                "${pmfatal} See the above usage notes for the possible choices.\n"
                "\n"
                )
    endif()
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Flags configuring the ParaMonte library build.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#option(DLL_ENABLED                             "Define DLL_ENABLED preprocessor flag" OFF)
#option(ParaMonte_OBJ_ENABLED                   "Build ParaMonte library object files" ON)
#option(ParaMonteTest_OBJ_ENABLED               "Build ParaMonte library test object files, intended only for developers" ON)
#option(ParaMonteTest_EXE_ENABLED               "Build ParaMonte library test executables, intended only for developers" OFF)
#option(ParaMonteTest_RUN_ENABLED               "Run ParaMonte library developer tests, intended only for developers" OFF)
#option(ParaMonteExampleMVNFortran_OBJ_ENABLED  "Build ParaMonte library examples' object files for Fortran " OFF)
#option(ParaMonteExampleMVNFortran_EXE_ENABLED  "Build ParaMonte library examples' executables for Fortran" OFF)
#option(ParaMonteExampleMVNFortran_RUN_ENABLED  "Run ParaMonte library examples' executables for Fortran" OFF)
#option(ParaMonteExampleMVNPython_BLD_ENABLED   "Build ParaMonte library examples' sctips for Python" OFF)
#option(ParaMonteExampleMVNPython_RUN_ENABLED   "Run ParaMonte library examples' sctips for Python" OFF)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Create the ParaMonte Version Preprocessor Flag. This is to avoid recompilation cascades caused by changing the value of the version macro
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    include(setParaMonteVersion)
    set(ParaMonteVersionFortranSourceFile "${paramonte_src_fortran_main_dir}/pm_paramonte@version.inc.F90")
    # do not regenerate the file if it already exists and the version is correct.
    # Otherwise, a file rewrite will lead to a recompilation which can be lengthy with interprocedural optimization.
    set(fileGenEnabled ON)
    if (EXISTS "${ParaMonteVersionFortranSourceFile}")
        file(READ  "${ParaMonteVersionFortranSourceFile}" contents)
        if ("${contents}" MATCHES "PARAMONTE_VERSION = \"${ParaMonteVersion}\"")
            set(fileGenEnabled OFF)
        endif()
    endif()
    if (fileGenEnabled)
        file(WRITE  "${ParaMonteVersionFortranSourceFile}"
                    "! WARNING - DO NOT CHANGE THE CONTENTS OF THIS FILE MANUALLY.\n"
                    "! WARNING - This file is auto-generated and overwritten by the ParaMonte CMake build scripts.\n"
                    "character(*, SK), parameter :: PARAMONTE_VERSION = \"${ParaMonteVersion}\""
                    )
    endif()
    #   Alternatively, we can also define the preprocessor macro PARAMONTE_VERSION.
    #   But changing the value of this macro (which is frequently done in production and release process(
    #   will lead to complete rebuild of the entire library which extremely lengthy and unnecessary.
    #   set(FPP_FLAGS "${FPP_FLAGS}" -DPARAMONTE_VERSION=\"'${ParaMonteVersion}'\") ! deprecated
    #   add_compile_definitions("PARAMONTE_VERSION=${ParaMonteVersion}")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Create the ParaMonte build timestamp for unsupported compilers. The supported compilers use the intrinsic `__TIMESTAMP__`.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (NOT (${csid_is_gnu} OR ${csid_is_intel}))
    string(TIMESTAMP PARAMONTE_BUILD_TIMESTAMP "%a %b %d %H:%M:%S %Y")
    add_compile_definitions("PARAMONTE_BUILD_TIMESTAMP=${PARAMONTE_BUILD_TIMESTAMP}")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Set preprocessor build flags.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#if ("${build}" STREQUAL "debug")
#    add_compile_definitions("DEBUG_ENABLED=1")
#elseif ("${build}" STREQUAL "testing")
#    add_compile_definitions("TESTING_ENABLD=1")
#elseif (NOT ("${build}" STREQUAL "release" OR "${build}" STREQUAL "ipo" OR "${build}" STREQUAL "tuned" OR "${build}" STREQUAL "native"))
#    message (FATAL_ERROR
#            "\n"
#            "${pmfatal} Internal ParaMonte CMake error occurred.\n"
#            "${pmfatal} The integrity of the ParaMonte CMake scripts is compromised.\n"
#            "${pmfatal} Unrecognized `build` option: build=${build}\n"
#            "${pmfatal} Clone the ParaMonte library from GitHub and retry:\n"
#            "\n"
#            "   https://github.com/cdslaborg/paramonte\n"
#            "\n"
#            )
#endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Define the ParaMonte parallelism preprocessor macros.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    unset(CAF_TYPE)
    set(CAF_ENABLED 0)
    set(MPI_ENABLED 0)
    set(OMP_ENABLED 0)
    if ("${par}" STREQUAL "none" OR "${par}" STREQUAL "") # MATCHES "[Nn][Oo][Nn][Ee]" AND NOT "${par}" STREQUAL "")
        set(par "serial")
    endif()
    string(TOLOWER ${par} par)
    if ("${par}" STREQUAL "cafsingle") # MATCHES "[Cc][Aa][Ff][Ss][Ii][Nn][Gg][Ll][Ee]")
        set(CAF_TYPE "single")
        set(CAF_ENABLED 1)
    elseif ("${par}" STREQUAL "cafshared") # MATCHES "[Cc][Aa][Ff][Ss][Hh][Aa][Rr][Ee][Dd]")
        set(CAF_TYPE "shared")
        set(CAF_ENABLED 1)
    elseif ("${par}" STREQUAL "cafdist") # MATCHES "[Cc][Aa][Ff][Dd][Ii][Ss][Tt]")
        set(CAF_TYPE "distributed")
        set(CAF_ENABLED 1)
    elseif ("${par}" STREQUAL "mpi") # MATCHES "[Mm][Pp][Ii]")
        set(MPI_ENABLED 1)
    elseif ("${par}" STREQUAL "openmp" OR "${par}" STREQUAL "omp") # MATCHES "[Oo][Mm][Pp]")
        set(OMP_ENABLED 1)
    elseif (NOT "${par}" STREQUAL "serial")
        printUsage()
        message(FATAL_ERROR
                "\n"
                "${pmfatal} The user-specified parallelism is unsupported. par=${par}\n"
                "${pmfatal} See the above usage notes for the possible choices.\n"
                "\n"
                )
    endif()
    add_compile_definitions("CAF_ENABLED=${CAF_ENABLED}")
    add_compile_definitions("MPI_ENABLED=${MPI_ENABLED}")
    add_compile_definitions("OMP_ENABLED=${OMP_ENABLED}")
endif()
if (${csid_is_intel} AND ${MPI_ENABLED})
    set(MPI_Fortran_COMPILER_FLAGS "-static_mpi")
    set(MPI_Fortran_LINK_FLAGS "-static_mpi")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Determine the supported types kinds by the specified Fortran compiler.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

include(setParaMonteKinds) # Defines all relevant <type>K_ENABLED macros to be lated included as the compiler FPP flags.

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Report build spec.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

message(NOTICE "${pmattn} operating system / platform: ${CMAKE_SYSTEM_NAME} / ${CMAKE_SYSTEM_PROCESSOR}")
message(NOTICE "${pmattn} C-Fortran compiler suite / version: ${CMAKE_Fortran_COMPILER_ID} / ${CMAKE_Fortran_COMPILER_VERSION}")
message(NOTICE "${pmattn} C-Fortran interoperation enabled - CFI_ENABLED: ${CFI_ENABLED}")
message(NOTICE "${pmattn} Coarray parallelization enabled - CAF_ENABLED: ${CAF_ENABLED}")
message(NOTICE "${pmattn} MPI parallelization enabled - MPI_ENABLED: ${MPI_ENABLED}")
message(NOTICE "${pmattn} library test enabled - TEST_RUN_ENABLED: ${ParaMonteTest_RUN_ENABLED}")
message(NOTICE "${pmattn} library build type: ${CMAKE_BUILD_TYPE}")
message(NOTICE "${pmattn} library type: ${lib}")
message(NOTICE "${pmattn} DLL_ENABLED: ${DLL_ENABLED}")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set the compiler flags (if the compiler is supported, otherwise, leave the default as it).
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    message(NOTICE "${pmattn} CMake default CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")
    message(NOTICE "${pmattn} CMake default CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
    message(NOTICE "${pmattn} CMake default CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
    message(NOTICE "${pmattn} CMake default CMAKE_Fortran_FLAGS_DEBUG: ${CMAKE_Fortran_FLAGS_DEBUG}")
    message(NOTICE "${pmattn} CMake default CMAKE_Fortran_FLAGS_RELEASE: ${CMAKE_Fortran_FLAGS_RELEASE}")
    message(NOTICE "${pmattn} CMake default CMAKE_Fortran_FLAGS_RELWITHDEBINFO: ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}")
    message(NOTICE "${pmattn} User-specified additional compiler flags: fcf=${fcf}")
    message(NOTICE "${pmattn} User-specified additional linker flags: flf=${flf}")
    set(FPP_ENABLED FALSE)
    set(INTEL_ENABLED 0)
    set(GNU_ENABLED 0)
    if (${csid_is_gnu} OR ${csid_is_intel})
        # Set the default mode for generating preprocessed files.
        message(NOTICE "${pmattn} fpp flag fpp=${fpp}")
        if ("${fpp}" STREQUAL "" OR NOT DEFINED fpp)
            if ("${lang}" STREQUAL "fortran")
                set(fpp "default")
            else()
                set(fpp "none")
            endif()
        endif()
        if ("${fpp}" MATCHES "[Dd][Ee][Ff][Aa][Uu][Ll][Tt]")
            set(FPP_ENABLED TRUE)
            if (${csid_is_gnu})
                set(GNU_ENABLED 1)
            else(${csid_is_intel})
                set(INTEL_ENABLED 1)
            endif()
        elseif ("${fpp}" MATCHES "[Gg][Ee][Nn][Ee][Rr][Ii][Cc]")
            set(FPP_ENABLED TRUE)
        elseif (NOT "${fpp}" MATCHES "[Nn][Oo][Nn][Ee]")
            message(FATAL_ERROR
                    "\n"
                    "${pmfatal} Unrecognized fpp option:fpp=${fpp}\n"
                    "${pmfatal} Allowed values:\n"
                    "\n"
                    "${pmfatal} none: Source files will not be preprocessed.\n"
                    "${pmfatal} gnu: Source files will be preprocessed according to the GNU Compiler syntax and toolset.\n"
                    "${pmfatal} intel: Source files will be preprocessed according to the Intel Compiler syntax and toolset.\n"
                    "${pmfatal} generic: Source files will be preprocessed assuming a generic Compiler syntax and toolset.\n"
                    "${pmfatal} default: Source files will be preprocessed assuming the default Compiler syntax and toolset.\n"
                    "\n"
                    )
        endif()
        unset(FC_FLAGS)
        unset(FL_FLAGS)
        unset(FCL_FLAGS)
        unset(FCL_FLAGS_DEBUG)
        unset(FCL_FLAGS_TESTING)
        unset(FCL_FLAGS_RELEASE)
        #unset(CMAKE_Fortran_FLAGS)
        option(MT_ENABLED "Enabled multithreaded library linking" ON)
        #   On exit, the following compiler-specific scripts append (but do not reset) the following CMake variables,
        #
        #   FL_FLAGS            :   The default link         flags for the Intel Fortran compiler.
        #   FC_FLAGS            :   The default compile      flags for the Intel Fortran compiler.
        #   FCL_FLAGS           :   The default compile/link flags for the Intel Fortran compiler.
        #   FCL_FLAGS_DEBUG     :   The debug   compile/link flags for the Intel Fortran compiler.
        #   FCL_FLAGS_TESTING   :   The testing compile/link flags for the Intel Fortran compiler.
        #   FCL_FLAGS_RELEASE   :   The release compile/link flags for the Intel Fortran compiler.
        #
        #   The scripts below are compiler specific and set the compile/link flags per Intel or GNU compilers.
        #   The compiler flags for all other unsupported compilers are currently left to be automatically set by CMake.
        if (${csid_is_gnu})
            include(setParaMonteFlagsGNU)
        elseif (${csid_is_intel})
            include(setParaMonteFlagsIntel)
        endif()
        #   User-specified flags will be added directly to the target library.
    #set(FC_FLAGS "${fcf}" "${FC_FLAGS}")
    #set(FL_FLAGS "${flf}" "${FC_FLAGS}")
        if ("${build}" STREQUAL "debug")
            set(FCL_FLAGS "${FCL_FLAGS}" "${FCL_FLAGS_DEBUG}")
        elseif ("${build}" STREQUAL "testing")
            set(FCL_FLAGS "${FCL_FLAGS}" "${FCL_FLAGS_TESTING}")
        elseif ("${build}" STREQUAL "release" OR "${build}" STREQUAL "ipo" OR "${build}" STREQUAL "tuned" OR "${build}" STREQUAL "native")
            set(FCL_FLAGS "${FCL_FLAGS}" "${FCL_FLAGS_RELEASE}")
        else()
            message(FATAL_ERROR
                    "\n"
                    "${pmfatal} Internal error occurred. Unrecognized `build`: ${build}\n"
                    "${pmfatal} The contents of the CMake file has been compromised.\n"
                    "\n"
                    )
        endif()
        message(NOTICE "${pmattn} Fortran compiler/linker flags - FCL_FLAGS: ${FCL_FLAGS}")
        message(NOTICE "${pmattn} Fortran compiler flags - FC_FLAGS: ${FC_FLAGS}")
        message(NOTICE "${pmattn} Fortran linker flags - FL_FLAGS: ${FL_FLAGS}")
    endif()
    #set(CMAKE_Fortran_FLAGS "${FPP_FLAGS}" CACHE STRING "Fortran compiler/linker flags" FORCE)
    #string(REPLACE ";" " " CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set number of Fortran Coarray (MPI) images (that will run in parallel, if Coarray/MPI parallelism is requested by the user).
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (NOT DEFINED nproc OR "${nproc}" STREQUAL "")
    set(nproc 3 CACHE STRING "Number of Fortran Coarray/MPI/OpenMP images/processes/threads." FORCE)
endif()

if (NOT DEFINED FOR_COARRAY_NUM_IMAGES)
    set(FOR_COARRAY_NUM_IMAGES "${nrpoc}" CACHE STRING "Number of Fortran Coarray images." FORCE)
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set ParaMonte parallelism: This MUST be included only after the above path setup. Kernel version MACRO setup depends on the paths.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

include(setParaMonteParallelism)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Set MPI library name.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    set(parname "${par}")
    unset(errorOccurred)
    if ("${par}" STREQUAL "mpi" AND NOT "${MPIEXEC_EXECUTABLE}" STREQUAL "")
        get_filename_component(mpiexec_dir "${MPIEXEC_EXECUTABLE}" DIRECTORY)
        get_filename_component(mpiexec_name "${MPIEXEC_EXECUTABLE}" NAME)
        execute_process(COMMAND sh "-c" "${mpiexec_name} --version"
                        RESULT_VARIABLE errorOccurred
                        OUTPUT_VARIABLE mpiVersionInfo
                        WORKING_DIRECTORY "${mpiexec_dir}"
                        OUTPUT_STRIP_TRAILING_WHITESPACE
                        )
        if(errorOccurred)
            set(mpiVersionInfo "${MPIEXEC_EXECUTABLE}")
        endif()
        if ("${mpiVersionInfo}" MATCHES [Ii][Nn][Tt][Ee][Ll])
            set(parname "impi")
            set(errorOccurred 0)
        elseif ("${mpiVersionInfo}" MATCHES [Mm][Pp][Ii][Cc][Hh])
            set(parname "mpich")
            set(errorOccurred 0)
        elseif ("${mpiVersionInfo}" MATCHES [oO][pP][eE][nN][rR][tT][eE] OR "${mpiVersionInfo}" MATCHES [oO][pP][eE][nN]-?[mM][pP][iI])
            set(parname "openmpi")
            set(errorOccurred 0)
        else() # Check for the presence of the mpichversion program. Assume mpich mpi implementation if it exists.
            set(mpichversion_path "${mpiexec_dir}/mpichversion")
            execute_process(COMMAND "${mpichversion_path}"
                            RESULT_VARIABLE errorOccurred
                            OUTPUT_VARIABLE mpiVersionInfo
                            WORKING_DIRECTORY "${CMAKE_ARGV3}"
                            OUTPUT_STRIP_TRAILING_WHITESPACE
                            )
            if(errorOccurred)
                set(errorOccurred 1)
            elseif ("${mpiVersionInfo}" MATCHES [Mm][Pp][Ii][Cc][Hh])
                set(parname "mpich")
                set(errorOccurred 0)
            else()
                set(errorOccurred 1)
            endif()
        endif()
    elseif ("${par}" STREQUAL "openmp" OR "${par}" STREQUAL "omp")
        find_package(OpenMP)
        set(parname "openmp")
    endif()
    if ("${errorOccurred}")
        message(WARNING
                "${pmwarn} Failed to fetch the mpi version information.\n"
                "${pmwarn} The library naming convention will default to `mpi`.\n"
                )
    endif()
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Find MATLAB if requested.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (MATLAB_ENABLED)
    if ("${matlabdir}" STREQUAL "")
        message(NOTICE "${pmattn} Searching for MATLAB...")
    else()
        set(Matlab_ROOT_DIR "${matlabdir}")
        message(NOTICE "${pmattn} Searching for MATLAB in the specified installation root directory: ${matlabdir}")
    endif()
    message(NOTICE "${pmattn} Sit back and relax...this search may take a few minutes...")
    find_package(Matlab)
    if (NOT Matlab_FOUND)
        message(WARNING
                "${pmwarn} Failed to find MATLAB on the system.\n"
                "${pmwarn} MATLAB is required for building the ParaMonte library MATLAB Mex files.\n"
                "${pmwarn} without MATLAB components, the ParaMonte MATLAB library will not function.\n"
                "${pmwarn} If you believe MATLAB exists on your system, specify the MATLAB directory containing\n"
                "${pmwarn} the target version of MATLAB to the ParaMonte build script via the command line argument:\n"
                "${pmwarn} \n"
                "${pmwarn}      install.sh --lang matlab --matlabdir path_to_matlab_root_directory\n"
                "${pmwarn} \n"
                "${pmwarn} or if you are calling CMake command directly, specify the path to the directory\n"
                "${pmwarn} by defining either `matlabdir` or `Matlab_ROOT_DIR` on the CMake command line:\n"
                "${pmwarn} \n"
                "${pmwarn}      cmake -Dlang=matlab -Dmatlabdir=path_to_matlab_root_directory\n"
                "${pmwarn} \n"
                "${pmwarn} The MATLAB root directory is the directory that is returned by the MATLAB command `matlabroot`.\n"
                "${pmwarn} An example of the MATLAB root directory (on Windows) is 'C:\\Program Files\\MATLAB\\R2017b'.\n"
                "${pmwarn} PROCEEDING WITHOUT COMPILING MEX FILES AND MATLAB DEPENDENCIES..."
                )
    endif()
    #matlab_add_mex(NAME mex_file_name SRC source_file.cpp)
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Find Python3 for the postprocessing of examples and benchmarks, if requested.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (NOT Python3_FOUND) # (NOT Python3_interpreter_FOUND)
    find_package(Python3 COMPONENTS Interpreter)
    if (Python3_Interpreter_FOUND)
        #set(Python3_ENABLED 1)
        message(NOTICE "${pmattn} Python3 Interpreter for postprocessing examples and benchmarks found.")
    else()
        message(NOTICE "${pmwarn} Python3 Interpreter for postprocessing examples and benchmarks was not found.")
    endif()
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Generate ParaMonte library build directories, object files, and shared libraries.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#include_directories(BEFORE ${paramonte_bld_lib_dir})

# message(NOTICE " ${pmattn} CMAKE_EXECUTABLE_SUFFIX_Fortran: ${CMAKE_EXECUTABLE_SUFFIX_Fortran}")

# Define the library name
if (DEFINED libname)
    message(NOTICE "${pmattn} User-specified ParaMonte library name detected. libname=\"${libname}\"")
else()
    #set(libname "libparamonte_${lang}_${os}_${arch}_${csid}${COMPILER_VERSION_MAJOR}_${build}_${lib}_${mem}_${parname}_${checking}")
    set(libname "libparamonte")
    #libname="paramonte_${lang}_${os}_${arch}_${csid}${CMAKE_Fortran_COMPILER_VERSION}_${build}_${lib}_${mem}${parname}"
endif()
set(libspec "${os}_${arch}_${csid}${COMPILER_VERSION_MAJOR}_${build}_${lib}_${mem}_${parname}_${checking}")
set(libpath "${os}/${arch}/${csid}${COMPILER_VERSION_MAJOR}/${build}/${lib}/${mem}/${parname}/${checking}")
if ("${lang}" STREQUAL "c" OR "${lang}" STREQUAL "cpp" OR "${lang}" STREQUAL "fortran")
    set(binname "${libname}_${lang}_${libspec}")
else()
    set(binname "${libname}_${lang}_${os}_${arch}")
endif()
message(NOTICE "${pmattn} The ParaMonte library name: libname=\"${libname}\"")
message(NOTICE "${pmattn} The ParaMonte library spec: libspec=\"${libspec}\"")
message(NOTICE "${pmattn} The ParaMonte library binary folder name: binname=\"${binname}\"")

# Have the .mod files placed in the lib folder
# set(CMAKE_Fortran_MODULE_DIRECTORY "${paramonte_bld_inc_dir}")

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
################################################## Initialize CMake project ########################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

project(ParaMonte VERSION "${ParaMonteVersion}" LANGUAGES C Fortran DESCRIPTION "A high-performance parallel Machine Learning library")
set(CMAKE_Fortran_PREPROCESS ON)
set(CMAKE_MACOSX_RPATH ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Find the BLAS library if needed.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    if (BLAS_REQUIRED)
        message(NOTICE "${pmattn} Searching for the requested BLAS library: ${blas}")
        if ("${DLL_ENABLED}" STREQUAL "1")
            set(BLA_STATIC OFF)
        else()
            set(BLA_STATIC ON)
        endif()
        if (NOT "${blas}" MATCHES "[aA][nN][yY]")
            set(BLA_VENDOR "${blas}")
        endif()
        find_package(BLAS REQUIRED)
        if (BLAS_FOUND)
            set(BLAS_ENABLED 1)
            message(NOTICE "${pmattn} The requested BLAS library found.")
            message(NOTICE "${pmattn} BLAS_LIBRARIES: ${BLAS_LIBRARIES}")
            message(NOTICE "${pmattn} BLAS_LINKER_FLAGS: ${BLAS_LINKER_FLAGS}")
        else()
            message(NOTICE "${pmwarn} BLAS library could not be found as requested. Skipping linking with BLAS...")
        endif()
    endif()
    add_compile_definitions("BLAS_ENABLED=${BLAS_ENABLED}")
    message(NOTICE "${pmattn} BLAS_ENABLED=${BLAS_ENABLED}")
    message(NOTICE "${pmattn} BLAS_REQUIRED=${BLAS_REQUIRED}")
    message(NOTICE "${pmattn} OpenBLAS_ENABLED=${OpenBLAS_ENABLED}")
    message(NOTICE "${pmattn} User-specified blas=${blas}")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Find the LAPACK library if needed.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)
    if (LAPACK_REQUIRED)
        # AND DEFINED bench AND NOT "${bench}" MATCHES "[nN][oO][nN][eE]")
        message(NOTICE "${pmattn} Searching for the requested LAPACK library: ${lapack}")
        if ("${DLL_ENABLED}" STREQUAL "1")
            set(BLA_STATIC OFF)
        else()
            set(BLA_STATIC ON)
        endif()
        if (NOT "${lapack}" MATCHES "[aA][nN][yY]")
            set(BLA_VENDOR "${lapack}")
        endif()
        find_package(LAPACK
                    #REQUIRED
                    )
        if (LAPACK_FOUND)
            set(LAPACK_ENABLED 1)
            message(NOTICE "${pmattn} The requested LAPACK library found.")
            message(NOTICE "${pmattn} BLAS_LIBRARIES: ${BLAS_LIBRARIES}")
            message(NOTICE "${pmattn} BLAS_LINKER_FLAGS: ${BLAS_LINKER_FLAGS}")
            message(NOTICE "${pmattn} LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}")
            message(NOTICE "${pmattn} LAPACK_LINKER_FLAGS: ${LAPACK_LINKER_FLAGS}")
        else()
            message(NOTICE "${pmwarn} LAPACK library could not be found as requested. Skipping linking with LAPACK...")
        endif()
    endif()
    add_compile_definitions("LAPACK_ENABLED=${LAPACK_ENABLED}")
    message(NOTICE "${pmattn} LAPACK_ENABLED=${LAPACK_ENABLED}")
    message(NOTICE "${pmattn} LAPACK_REQUIRED=${LAPACK_REQUIRED}")
    message(NOTICE "${pmattn} User-specified lapack vendor: ${lapack}")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set up the ParaMonte Fortran build.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (NOT ${csid_is_gnu} AND NOT ${csid_is_intel})
    add_compile_definitions("__FILE__=\"UNKNOWN_FILE\"")
    add_compile_definitions("__LINE__=0")
endif()
add_compile_definitions("GNU_ENABLED=${GNU_ENABLED}")
add_compile_definitions("INTEL_ENABLED=${INTEL_ENABLED}")
message(NOTICE "${pmattn} fpp flag GNU_ENABLED=${GNU_ENABLED}")
message(NOTICE "${pmattn} fpp flag INTEL_ENABLED=${INTEL_ENABLED}")

# This phony target must be defined before adding the main fortran source subdirectory
# because the fortran fpp processing depends on it.
add_custom_target(package)

# The source for the ParaMonte library and its binaries in the paramonte_bld_obj_dir folder
add_subdirectory("${paramonte_src_fortran_main_dir}" "${paramonte_bld_obj_dir}")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set up the ParaMonte MATLAB MEX build if needed.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (MATLAB_ENABLED)
    set(paramonte_bld_xrc_dir "${paramonte_bld_dir}/xrc")
    message(NOTICE "${pmattn} paramonte_bld_xrc_dir: ${paramonte_bld_xrc_dir}")
    if (NOT EXISTS "${paramonte_bld_xrc_dir}/")
        file(MAKE_DIRECTORY "${paramonte_bld_xrc_dir}")
    endif()
    set(paramonte_src_matlab_xrc_dir "${paramonte_src_matlab_dir}/xrc")
    message(NOTICE "${pmattn} paramonte_src_matlab_xrc_dir: ${paramonte_src_matlab_xrc_dir}")
    add_subdirectory("${paramonte_src_matlab_xrc_dir}" "${paramonte_bld_xrc_dir}")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copy the auxiliary files to the package directory.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)

    # This is now done in the paramonte_src_matlab_xrc_dir subdirectory above.
    #if (FPP_ENABLED)
    #    include(setParaMonteFPP)
    #endif()

    #   INCLUDE/interface files.

    if ("${lang}" STREQUAL "c" OR "${lang}" STREQUAL "cpp" OR "${lang}" STREQUAL "fortran")
        file(MAKE_DIRECTORY "${paramonte_bld_pkg_inc_dir}")
    endif()

    if ("${lang}" STREQUAL "c")

        set(origin "${paramonte_src_${lang}_main_dir}/inc")
        set(destin "${paramonte_bld_pkg_inc_dir}")
        message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library include files to the package directory...\n"
                        "${pmattn}      - from: ${origin}\n"
                        "${pmattn}      -   to: ${destin}")
        add_custom_target(copyDirInclude COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}" )
        add_dependencies(package copyDirInclude)

    elseif ("${lang}" STREQUAL "cpp")

        set(origin "${paramonte_src_${lang}_main_dir}/inc")
        set(destin "${paramonte_bld_pkg_inc_dir}")
        message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library include files to the package directory...\n"
                        "${pmattn}      - from: ${origin}\n"
                        "${pmattn}      -   to: ${destin}")
        add_custom_target(copyDirInclude COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}")
        add_dependencies(package copyDirInclude)

    elseif ("${lang}" STREQUAL "fortran")

        set(origin "${paramonte_bld_inc_dir}")
        set(destin "${paramonte_bld_pkg_inc_dir}")
        message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library include files to the package directory...\n"
                        "${pmattn}      - from: ${origin}\n"
                        "${pmattn}      -   to: ${destin}")
        add_custom_target(copyDirInclude COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}")
        add_dependencies(package copyDirInclude)

    elseif ("${lang}" STREQUAL "julia" OR "${lang}" STREQUAL "matlab" OR "${lang}" STREQUAL "python" OR "${lang}" STREQUAL "r")

        if ("${lang}" STREQUAL "matlab")
            set(pkgname "+pm")
        else()
            set(pkgname "paramonte")
        endif()

        set(paramonte_bld_pkg_paramonte_dir "${paramonte_bld_pkg_dir}/${pkgname}")

        # Copy library files

        set(origin "${paramonte_src_${lang}_main_dir}/${pkgname}")
        set(destin "${paramonte_bld_pkg_paramonte_dir}")
        message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library include files to the package directory...\n"
                        "${pmattn}      - from: ${origin}\n"
                        "${pmattn}      -   to: ${destin}")
        add_custom_target(copyDirInclude COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}")
        add_dependencies(package copyDirInclude)

        # Copy auxiliary files

        set(origin "${paramonte_src_dir}/auxil")
        set(destin "${paramonte_bld_pkg_paramonte_dir}/auxil")
        message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library auxiliary files to the package directory...\n"
                        "${pmattn}      - from: ${origin}\n"
                        "${pmattn}      -   to: ${destin}")
        add_custom_target(copyDirAuxil COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}")
        add_dependencies(package copyDirAuxil)

        if ("${lang}" STREQUAL "python")

            set(origin "${paramonte_src_${lang}_dir}/setup")
            set(destin "${paramonte_bld_pkg_dir}/setup")
            message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library setup files to the package directory...\n"
                            "${pmattn}      - from: ${origin}\n"
                            "${pmattn}      -   to: ${destin}")
            add_custom_target(copyDirSetup COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}")
            add_dependencies(package copyDirSetup)

        endif()

    endif()

    #   pmlib files

    if (lang_is_dynamic)
        set(paramonte_bld_pkg_lib_dir "${paramonte_bld_pkg_paramonte_dir}/lib/${libpath}")
    endif()

    set(origin "${paramonte_bld_lib_dir}")
    set(destin "${paramonte_bld_pkg_lib_dir}")
    message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library files to the package directory...\n"
                    "${pmattn}      - from: ${origin}\n"
                    "${pmattn}      -   to: ${destin}")
    add_custom_target(copyFileLib COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}")
    add_dependencies(package copyFileLib)

    #   LICENSE

    set(origin "${paramonte_dir}/LICENSE.md")
    set(destin "${paramonte_bld_pkg_dir}/LICENSE.md")
    message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library LICENSE.md file to the package directory...\n"
                    "${pmattn}      - from: ${origin}\n"
                    "${pmattn}      -   to: ${destin}")
    add_custom_target(copyFileLICENSE COMMAND ${CMAKE_COMMAND} -E copy "${origin}" "${destin}")
    add_dependencies(package copyFileLICENSE)

    #   README

    set(origin "${paramonte_src_${lang}_dir}/README.md")
    set(destin "${paramonte_bld_pkg_dir}/README.md")
    message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library README.md file to the package directory...\n"
                    "${pmattn}      - from: ${origin}\n"
                    "${pmattn}      -   to: ${destin}")
    add_custom_target(copyFileREADME COMMAND ${CMAKE_COMMAND} -E copy "${origin}" "${destin}")
    add_dependencies(package copyFileREADME)

    #   CHANGES

    set(origin "${paramonte_src_${lang}_dir}/CHANGES.md")
    set(destin "${paramonte_bld_pkg_dir}/CHANGES.md")
    message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library CHANGES.md file to the package directory...\n"
                    "${pmattn}      - from: ${origin}\n"
                    "${pmattn}      -   to: ${destin}")
    add_custom_target(copyFileCHANGES COMMAND ${CMAKE_COMMAND} -E copy "${origin}" "${destin}")
    add_dependencies(package copyFileCHANGES)

    #   VERSION

    set(origin "${paramonte_src_${lang}_dir}/VERSION.md")
    set(destin "${paramonte_bld_pkg_dir}/VERSION.md")
    message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library VERSION.md file to the package directory...\n"
                    "${pmattn}      - from: ${origin}\n"
                    "${pmattn}      -   to: ${destin}")
    add_custom_target(copyFileVERSION COMMAND ${CMAKE_COMMAND} -E copy "${origin}" "${destin}")
    add_dependencies(package copyFileVERSION)

    if ("${lang}" STREQUAL "julia" OR "${lang}" STREQUAL "matlab" OR "${lang}" STREQUAL "python" OR "${lang}" STREQUAL "r")
        set(origin "${paramonte_dir}/VERSION.md")
        set(destin "${paramonte_bld_pkg_paramonte_dir}/auxil/VERSION.md")
        message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library VERSION.md file to the package directory...\n"
                        "${pmattn}      - from: ${origin}\n"
                        "${pmattn}      -   to: ${destin}")
        add_custom_target(copyFileAuxilVERSION COMMAND ${CMAKE_COMMAND} -E copy "${origin}" "${destin}")
        add_dependencies(package copyFileAuxilVERSION)
    endif()

    # Copy examples/benchmarks folder to the package folder.

    #set(collections "benchmark" "example")
    #foreach(collection in ${collections})
    #    set(origin "${paramonte_${collection}_${lang}_dir}/") # The ending forward slash is essential for copying the contents of the ${collection} language folder.
    #    set(destin "${paramonte_bld_pkg_${collection}_dir}")
    #    message(NOTICE  "${pmattn} Copying the paramonte::${lang} library ${collection} files to the package directory...\n"
    #                    "${pmattn}      - from: ${origin}\n"
    #                    "${pmattn}      -   to: ${destin}")
    #    if (EXISTS "${origin}")
    #        file(COPY "${origin}" DESTINATION "${destin}")
    #    else()
    #        message(NOTICE  "${pmwarn} The paramonte library ${lang} ${collection} folder does not exists. skipping...")
    #    endif()
    #    if (EXISTS "${paramonte_${collection}_dir}/generic")
    #        message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library generic ${collection} files to the package directory...\n"
    #                        "${pmattn}      - from: ${paramonte_${collection}_dir}/generic\n"
    #                        "${pmattn}      -   to: ${destin}")
    #        file(COPY "${paramonte_${collection}_dir}/generic/" DESTINATION "${destin}")
    #    else()
    #        message(NOTICE  "${pmwarn} The paramonte library generic ${collection} folder does not exists. skipping...")
    #    endif()
    #endforeach()

endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Define the binary installation directory and copy binaries.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#add_custom_target(deploy)
#add_dependencies(deploy package)
if (DEFINED ddir)
    set(paramonte_bin_dir "${ddir}")
else()
    set(paramonte_bin_dir "${paramonte_dir}/bin")
endif()

foreach(fitem ${fresh})
    if ("${fitem}" STREQUAL "ddir" OR "${fitem}" STREQUAL "all")
        message(NOTICE "${pmattn} Purging the entire ParaMonte library deploy directory: \"${paramonte_bin_dir}/\"")
        file(REMOVE_RECURSE "${paramonte_bin_dir}/")
        break()
    endif()
endforeach()

# Copy the package to the binary installation directory.

set(origin "${paramonte_bld_pkg_dir}/")
set(destin "${paramonte_bin_dir}/${binname}")
foreach(fitem ${fresh})
    if ("${fitem}" STREQUAL "deploy" OR "${fitem}" STREQUAL "all")
        message(NOTICE "${pmattn} Purging the current ParaMonte library deploy directory: \"${destin}/\"")
        file(REMOVE_RECURSE "${destin}/")
        break()
    endif()
endforeach()

if (NOT ${codecov_enabled})
    message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library package to the final binary installation directory...\n"
                    "${pmattn}      - from: ${origin}\n"
                    "${pmattn}      -   to: ${destin}")
    add_custom_target(deploy COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}")
    message(NOTICE  "${pmattn} Adding the rule for creating the paramonte::${lang} library binary installation directory...")
    if (NOT EXISTS "${paramonte_bin_dir}")
        #file(MAKE_DIRECTORY "${paramonte_bin_dir}")
        add_custom_target(makeBinDir ALL COMMAND ${CMAKE_COMMAND} -E make_directory "${paramonte_bin_dir}")
        add_dependencies(deploy makeBinDir) # ensure the folder is created before deployment.
    endif()
else()
    message(NOTICE  "${pmattn} Skipping the rule for copying the paramonte::${lang} library package to the final binary installation directory...")
    add_custom_target(deploy)
endif()
add_dependencies(deploy package)

if (FALSE)
    #### There is no need for this check anymore as the
    #### CMake variable CMAKE_MAKE_PROGRAM can achieve the goal.
    if (C_ENABLED OR CPP_ENABLED OR FORTRAN_ENABLED)
        if (NOT DEFINED "${bldbin}" AND NOT "${CMAKE_GENERATOR}" STREQUAL "Unix Makefiles")
            message(NOTICE  "${pmwarn} The CMake variable `bldbin` was not set at the time of invoking CMake.")
            message(NOTICE  "${pmwarn} You can set it to the name of the build software binary to be used")
            message(NOTICE  "${pmwarn} by CMake for creating CMakeLists for library examples and benchmarks.")
            message(NOTICE  "${pmwarn} In MinGW environments, the `make` build software has the default name `mingw32-make` or `mingw32-make.exe`.")
            message(NOTICE  "${pmwarn} In MSYS environments, the `make` build software has the default name `make`.")
            message(NOTICE  "${pmwarn} Depending on your environment, you can set `-Dbldbin=\"BLDBIN\"`")
            message(NOTICE  "${pmwarn} where BLDBIN should be replaced with the actual `make` binary name.")
            message(NOTICE  "${pmwarn} When unset, the ParaMonte CMake script will attempt to infer the most")
            message(NOTICE  "${pmwarn} likely name for the build software based on the specified CMAKE_GENERATOR.")
            message(NOTICE  "${pmwarn} If the build software is incorrectly identified, the library examples and")
            message(NOTICE  "${pmwarn} benchmarks build scripts may not work. It is better to be set by the user.")
            if ("${CMAKE_GENERATOR}" STREQUAL "NMake Makefiles")
                set(bldbin "nmake")
            elseif ("${CMAKE_GENERATOR}" STREQUAL "Watcom WMake")
                set(bldbin "wmake")
            elseif ("${CMAKE_GENERATOR}" MATCHES ".*Makefiles.*")
                find_program(bldbin make DOC "Location of Make build software executable.")
                if ("${bldbin}" STREQUAL "bldbin-NOTFOUND")
                    message(NOTICE "${pmwarn} Failed to find the `make` build binary.")
                    if (MSYS OR MINGW)
                        message(NOTICE  "${pmwarn} Searching instead for the `mingw32-make` build binary.")
                        find_program(bldbin mingw32-make DOC "Location of Make build software executable.")
                        if ("${bldbin}" STREQUAL "bldbin-NOTFOUND")
                            find_program(bldbin mingw32-make.exe DOC "Location of Make build software executable.")
                        endif()
                        if ("${bldbin}" STREQUAL "bldbin-NOTFOUND")
                            unset(bldbin)
                        endif()
                    else()
                        unset(bldbin)
                    endif()
                endif()
            elseif ("${CMAKE_GENERATOR}" MATCHES ".*Ninja.*")
                set(bldbin "ninja")
            endif()
            if (DEFINED bldbin)
                message(NOTICE "${pmwarn} The identified build software binary: ${bldbin}.")
            else()
                message(NOTICE "${pmwarn} Failed to identify the build software binary.")
                message(NOTICE "${pmwarn} Proceeding with `make` as the default build software binary...")
                set(bldbin make)
            endif()
        elseif("${CMAKE_GENERATOR}" STREQUAL "Unix Makefiles")
            set(bldbin "make")
        endif()
    endif()
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# The source for the ParaMonte library examples.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (EXISTS "${paramonte_example_dir}")
    add_subdirectory("${paramonte_example_dir}" "${paramonte_bld_example_dir}")
else()
    message(WARNING
            "\n"
            "${pmwarn} The example directory in the root directory of the ParaMonte library is missing.\n"
            "${pmwarn} paramonte_example_dir=${paramonte_example_dir}\n"
            "${pmwarn} Skipping the ParaMonte library example build...\n"
            "\n"
            )
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# The source for the ParaMonte library benchmarks.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (EXISTS "${paramonte_benchmark_dir}")
    add_subdirectory("${paramonte_benchmark_dir}" "${paramonte_bld_benchmark_dir}")
else()
    message(WARNING
            "\n"
            "${pmwarn} The benchmark directory in the root directory of the ParaMonte library is missing.\n"
            "${pmwarn} paramonte_benchmark_dir=${paramonte_benchmark_dir}\n"
            "${pmwarn} Skipping the ParaMonte library benchmark build...\n"
            "\n"
            )
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# The source for the ParaMonte library tests and its binaries in the paramonte_bld_dir folder. Build tests only if requested.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (TRUE)

    if ("${test}" STREQUAL "" OR NOT DEFINED test)
        set(test "none")
    endif()
    string(TOLOWER "${test}" test)
    if ("${test}" STREQUAL "all" OR "${test}" STREQUAL "basic" OR "${test}" STREQUAL "sampler")

        message(NOTICE "${pmattn} Enabling the ParaMonte library ${test} tests...")

        # set up test build directory

        set(paramonte_bld_test_dir "${paramonte_bld_dir}/test")
        message(NOTICE "${pmattn} paramonte_bld_test_dir: ${paramonte_bld_test_dir}")
        if (NOT EXISTS "${paramonte_bld_test_dir}/")
            file(MAKE_DIRECTORY "${paramonte_bld_test_dir}")
        endif()

        # set up test object directory

        set(paramonte_bld_test_obj_dir "${paramonte_bld_test_dir}/obj")
        message(NOTICE "${pmattn} paramonte_bld_test_obj_dir: ${paramonte_bld_test_obj_dir}")
        if (NOT EXISTS "${paramonte_bld_test_obj_dir}/")
            file(MAKE_DIRECTORY "${paramonte_bld_test_obj_dir}")
        endif()

        # set up test module directory

        set(paramonte_bld_test_inc_dir "${paramonte_bld_test_dir}/inc")
        message(NOTICE "${pmattn} paramonte_bld_test_inc_dir: ${paramonte_bld_test_inc_dir}")
        if (NOT EXISTS "${paramonte_bld_test_inc_dir}/")
            file(MAKE_DIRECTORY "${paramonte_bld_test_inc_dir}")
        endif()

        # set up test bin directory

        #set(paramonte_bld_test_bin_dir "${paramonte_bld_test_dir}")
        set(paramonte_bld_test_bin_dir "${paramonte_bld_test_dir}/bin")
        message(NOTICE "${pmattn} paramonte_bld_test_bin_dir: ${paramonte_bld_test_bin_dir}")
        file(REMOVE_RECURSE "${paramonte_bld_test_bin_dir}/")
        if (NOT EXISTS "${paramonte_bld_test_bin_dir}/")
            file(MAKE_DIRECTORY "${paramonte_bld_test_bin_dir}")
        endif()

        # set up test bin directory

        set(paramonte_bld_test_bin_input_dir "${paramonte_bld_test_bin_dir}/input")
        message(NOTICE "${pmattn} paramonte_bld_test_bin_input_dir: ${paramonte_bld_test_bin_input_dir}")
        if (NOT EXISTS "${paramonte_bld_test_bin_input_dir}/")
            file(MAKE_DIRECTORY "${paramonte_bld_test_bin_input_dir}")
        endif()

        add_subdirectory("${paramonte_src_fortran_test_dir}" "${paramonte_bld_test_obj_dir}")

    elseif ("${test}" STREQUAL "none" OR NOT DEFINED test)

        message(NOTICE "${pmwarn} Skipping the ParaMonte library tests build as requested...")
        add_custom_target(test)

    else()

        printUsage()
        message(FATAL_ERROR
                "\n"
                "${pmfatal} Unsupported user-specified ParaMonte library test mode. test=${test}\n"
                "${pmfatal} Follow the guidelines above to appropriately specify the testing mode or drop the option.\n"
                "\n"
                )

    endif()

endif()

#message(NOTICE "CMAKE_Fortran_LINK_EXECUTABLE: ${CMAKE_Fortran_LINK_EXECUTABLE}")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Add a distclean target to the Makefile
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#add_custom_target(distclean COMMAND ${CMAKE_COMMAND} -P ${CMAKE_MODULE_PATH}/distclean.cmake)
#add_custom_target(dclean COMMAND ${CMAKE_COMMAND} -P ${CMAKE_MODULE_PATH}/distclean.cmake)