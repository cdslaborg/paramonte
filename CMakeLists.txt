####################################################################################################################################
####################################################################################################################################
####
####   MIT License
####
####   ParaMonte: plain powerful parallel Monte Carlo library.
####
####   Copyright (C) 2012-present, The Computational Data Science Lab
####
####   This file is part of the ParaMonte library.
####
####   Permission is hereby granted, free of charge, to any person obtaining a 
####   copy of this software and associated documentation files (the "Software"), 
####   to deal in the Software without restriction, including without limitation 
####   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
####   and/or sell copies of the Software, and to permit persons to whom the 
####   Software is furnished to do so, subject to the following conditions:
####
####   The above copyright notice and this permission notice shall be 
####   included in all copies or substantial portions of the Software.
####
####   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
####   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
####   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
####   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
####   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
####   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
####   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
####
####   ACKNOWLEDGMENT
####
####   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
####   As per the ParaMonte library license agreement terms, if you use any parts of 
####   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
####   work (education/research/industry/development/...) by citing the ParaMonte 
####   library as described on this page:
####
####       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
####
####################################################################################################################################
####################################################################################################################################

# CMake project file for ParaMonte
# idev -N 1 -n 68 -m 15 # on KNL
# idev -p skx-dev -N 1 -n 48 -t 00:15:00 # on SKX
# cmake .. -DPMCS=gnu -DMPI_ENABLED=TRUE -DCFI_ENABLED=TRUE -DDLL_ENABLED=TRUE -DCMAKE_BUILD_TYPE=debug
# cmake ../../.. -DPMCS=gnu -DBTYPE=release -DCAFTYPE=distributed
# make
# make install

 
cmake_minimum_required(VERSION 3.10)

cmake_policy(VERSION 3.10...3.14)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Define CMake colors
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if(WIN32)
    message( FATAL_ERROR 
    "-- ParaMonte - FATAL ERROR:\n"
    "   Do not build with CMAKE on Windows OS. Instead follow the instructions below:\n"
    "   1. Install Microsoft Visual Studio (MSVC).\n"
    "      You can obtain the MSVS Cummunity edition for free from the MS website.\n"
    "   2. Install Intel Parallel Studio for Windows. If you are a student, teacher,\n"
    "      or open-source develper, you can get a lisence for free from Intel.\n"
    "   3. Press the Windows key and find the Intel's dedicated Windows Terminal\n"
    "      by searching for the keyword 'intel', then open the terminal.\n"
    "   4. On the command line, go to the root directory of ParaMonte library.\n"
    "   5. Use the Windows Batch file 'configParaMonte.bat' to configure your build.\n"
    "      If you don't know how to do it, leave this file untouched\n"
    "      (the default values will be used).\n"
    "   6. On the command line run the Windows Batch file 'buildParaMonte.bat'.\n\n"
    )
else()
    string(ASCII 27 Esc)
    set(ColorReset  "${Esc}[m")
    set(ColorBold   "${Esc}[1m")
    set(Red         "${Esc}[31m")
    set(Green       "${Esc}[32m")
    set(Yellow      "${Esc}[33m")
    set(Blue        "${Esc}[34m")
    set(Magenta     "${Esc}[35m")
    set(Cyan        "${Esc}[36m")
    set(White       "${Esc}[37m")
    set(BoldRed     "${Esc}[1;31m")
    set(BoldGreen   "${Esc}[1;32m")
    set(BoldYellow  "${Esc}[1;33m")
    set(BoldBlue    "${Esc}[1;34m")
    set(BoldMagenta "${Esc}[1;35m")
    set(BoldCyan    "${Esc}[1;36m")
    set(BoldWhite   "${Esc}[1;37m")
endif()

# function(message)
    # list(GET ARGV 0 MessageType)
    # if(MessageType STREQUAL FATAL_ERROR OR MessageType STREQUAL SEND_ERROR)
        # list(REMOVE_AT ARGV 0)
        # _message(${MessageType} "${BoldRed}${ARGV}${ColourReset}")
    # elseif(MessageType STREQUAL WARNING)
        # list(REMOVE_AT ARGV 0)
        # _message(${MessageType} "${BoldYellow}${ARGV}${ColourReset}")
    # elseif(MessageType STREQUAL AUTHOR_WARNING)
        # list(REMOVE_AT ARGV 0)
        # _message(${MessageType} "${BoldCyan}${ARGV}${ColourReset}")
    # elseif(MessageType STREQUAL STATUS)
        # list(REMOVE_AT ARGV 0)
        # _message(${MessageType} "${Green}${ARGV}${ColourReset}")
    # else()
        # _message("${ARGV}")
    # endif()
# endfunction()

set(pmcolor "${BoldCyan}")
set(pmattn  " ${pmcolor}ParaMonte -${ColorReset}")
set(pmnote  "${pmattn} ${BoldMagenta}NOTE:${ColorReset}")
set(pmwarn  "${pmattn} ${BoldYellow}WARNING:${ColorReset}")
set(pmfatal "${pmattn} ${BoldRed}FATAL ERROR:${ColorReset}")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Print an error message on an attempt to build inside the source directory tree:
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if ("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
    message ( FATAL_ERROR 
            " \n"
            "${pmfatal}\n"
            "   CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}\n"
            "   CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}\n\n"
            "   This archive does not support in-source builds.\n\n"
            "   You must now delete the CMakeCache.txt file and the CMakeFiles/ directory under\n"
            "   in the ParaMonte library's root directory, otherwise you will not be able to configure correctly.\n"
            "   You must now run something like:\n\n"
            "       $ rm -r CMakeCache.txt CMakeFiles/\n\n"
            "   Please create a directory outside the ParaMonte library source tree or\n"
            "   create a directory inside the ParaMonte library source tree, then change the directory to it and\n"
            "   build under that outside directory, in a manner such as the following,\n\n"
            "       $ buildDir=$(uname -s)$(uname -m) # set build folder to OS name + architecture\n"
            "       $ mkdir $\{buildDir\}\n"
            "       $ cd $\{buildDir\}\n"
            "       $ CC=gcc FC=gfortran cmake -DBUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/install/dir /path/to/ParaMonte/src/dir\n\n"
            "   substituting the appropriate syntax for your shell (the above line assumes the bash shell)."
            )
endif()

if(DEFINED ENV{SOURCE_DATE_EPOCH})
    set(SOURCE_DATE_EPOCH "$ENV{SOURCE_DATE_EPOCH}")
endif()

set(CMAKE_SKIP_RPATH ON CACHE BOOL "Skip adding a build-dir rpath")
set(CMAKE_BUILD_WITH_INSTALL_RPATH ON CACHE BOOL "Build using the install rpath")
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON CACHE BOOL "Use relative rpaths")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Add our local modules to the module path
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_USER_MAKE_RULES_OVERRIDE "${CMAKE_SOURCE_DIR}/cmake/Modules/overrideDefault.cmake")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# banner
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

file(READ "${CMAKE_CURRENT_LIST_DIR}/auxil/.ParaMonteBanner" ParaMonte_BANNER)
message("${pmcolor}\n${ParaMonte_BANNER}${ColorReset}")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# define the build mode
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set ( CMAKE_CONFIGURATION_TYPES "Debug" "testing" "Release" ) # "RelWithDebInfo"

if (NOT DEFINED BTYPE)
    set ( BTYPE "Release" CACHE STRING "ParaMonte library build type." )
endif()
set (CMAKE_BUILD_TYPE ${BTYPE} CACHE STRING "Select CMAKE configuration to build" FORCE)

string(TOLOWER ${BTYPE} BTYPE)

set_property ( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${CMAKE_CONFIGURATION_TYPES} )

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Add option and check environment to determine if developer tests should be run
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# set( _TF OFF )
# if( NOT DEFINED ENV{ParaMonte_DEVELOPER})
  # set ( ENV{ParaMonte_DEVELOPER} FALSE )
# endif()
# if($ENV{ParaMonte_DEVELOPER})
  # set( _TF ON )
# endif()

# option(ParaMonteTest_RUN_ENABLED "Run tests intended only for developers" ${_TF})

# mark_as_advanced(ParaMonteTest_RUN_ENABLED)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set the desired compiler suite
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

include(SetParaMonteCompilerSuite)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Name project and specify source languages. Parse version from .VERSION file so that more info can be added and easier to get from scripts
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

include(SetParaMonteVersion)

project(ParaMonte   VERSION "${PARAMONTE_CMAKE_PROJECT_VERSION}"
                    #LANGUAGES Fortran C
                    LANGUAGES Fortran
                    DESCRIPTION "A high-performance parallel Monte Carlo library"
                    )
set(CMAKE_Fortran_STANDARD 2008)
#enable_language(Fortran C)
enable_language(Fortran)


if (NOT DEFINED PMCS)
    set(PMCS "${CMAKE_Fortran_COMPILER_ID}" CACHE STRING "The requested ParaMonte build compiler suite" FORCE)
endif()
string(TOLOWER ${PMCS} PMCS)

if ("${PMCS}" MATCHES "[Ii][Nn][Tt][Ee][Ll]")
    set(requested_compiler_suite "Intel")
    set(COMPILER_SUITE "Intel")
    set(intel_compiler TRUE)
    message ( STATUS "${pmnote} requested_compiler_suite=${requested_compiler_suite}" )
elseif("${PMCS}" MATCHES "[Gg][Nn][Uu]")
    set (requested_compiler_suite "GNU")
    message ( STATUS "${pmnote} requested_compiler_suite=${requested_compiler_suite}" )
    set(COMPILER_SUITE "GNU")
    set(gnu_compiler TRUE)
else()
    message ( FATAL_ERROR 
            " \n"
            "${pmfatal}\n"
            "   The compiler suite for building ParaMonte library as specified by \n"
            "   CMAKE_Fortran_COMPILER_ID=${PMCS} does not correspond to any of the two supported\n"
            "   comipler suites: Intel or GNU. Please ensure one of these compiler suites\n"
            "   is installed on your system, then provide its name to cmake by passing\n"
            " \n"
            "       -DPMCS=Intel\n"
            " \n"
            "       or\n"
            " \n"
            "       -DPMCS=GNU\n"
            " \n"
            "   when invoking cmake."
            )
endif()
string(TOLOWER ${COMPILER_SUITE} COMPILER_SUITE)

message(STATUS "CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}")

#if(DEFINED CMAKE_BUILD_TYPE AND "${CMAKE_BUILD_TYPE}" MATCHES "[Rr][Ee][Ll]") # release mode requested
#    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-working-directory")
#        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdebug-prefix-map=${CMAKE_SOURCE_DIR}=.")
#    endif()
#    if(CMAKE_C_COMPILER_ID MATCHES "GNU")
#        set(CMAKE_C_FLAGS       "${CMAKE_C_FLAGS} -fno-working-directory")
#        set(CMAKE_C_FLAGS       "${CMAKE_C_FLAGS} -fdebug-prefix-map=${CMAKE_SOURCE_DIR}=.")
#    endif()
#endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Test modern/Coarray Fortran compiler awareness
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

include(TestFortranCompiler)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ParaMonte version
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


message( STATUS "${pmattn} Building ParaMonte version: ${full_git_describe}" )
set(ParaMonte_dist_string "ParaMonte-${full_git_describe}")
message( STATUS "${pmattn} Building for target architecture: ${CMAKE_SYSTEM_PROCESSOR}" )

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Silence warnings about dereferencing unset variables
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if(NOT CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_FLAGS "")
endif()
if(NOT CMAKE_REQUIRED_LIBRARIES)
  set(CMAKE_REQUIRED_LIBRARIES "")
endif()
if(NOT CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_INCLUDES "")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set Fortran/C compiler suite: currently only intel is supported on Windows systems
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#set ( COMPILER_SUITE "intel" CACHE STRING "Fortran/C compiler suite (intel/gcc)" )
#message( STATUS "${pmattn} Building via the ${COMPILER_SUITE} compiler suit..." )
#if (COMPILER_SUITE MATCHES "intel|Intel|INTEL")
#    set (CMAKE_Fortran_COMPILER_ID "Intel" CACHE STRING "Fortran compiler ID" FORCE )
#    set (CMAKE_CXX_COMPILER_ID "Intel" CACHE STRING "CXX compiler ID" FORCE )
#    set (CMAKE_C_COMPILER_ID "Intel" CACHE STRING "C compiler ID" FORCE )
#elseif(COMPILER_SUITE MATCHES "gcc|GCC")
#    set (CMAKE_Fortran_COMPILER_ID "GNU" CACHE STRING "Fortran compiler ID" FORCE )
#    set (CMAKE_CXX_COMPILER_ID "GNU" CACHE STRING "CXX compiler ID" FORCE )
#    set (CMAKE_C_COMPILER_ID "GNU" CACHE STRING "C compiler ID" FORCE )
#endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Report untested Fortran compiler unless explicitly directed to build all examples.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" )
    set(gnu_compiler TRUE)
    set ( CMAKE_C_FLAGS_CODECOVERAGE "-fprofile-arcs -ftest-coverage -O0" CACHE STRING "Code coverage C compiler flags" )
    set ( CMAKE_Fortran_FLAGS_CODECOVERAGE "-fprofile-arcs -ftest-coverage -O0" CACHE STRING "Code coverage Fortran compiler flags" )
elseif( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel" )
    set(intel_compiler TRUE)
else()
    message ( WARNING
            " \n"
            " ${pmwarn}\n"
            " ${pmattn} Attempting to build with untested Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}.\n"
            " ${pmattn} Please report any failures to shahmoradi@utexas.edu"
            )
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set CMAKE_Fortran_COMPILER_VERSION if CMake doesn't do it for us
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

include(SetParaMonteCompilerVersion)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Check Simple Coarray Fortran Source Compiles
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#if(intel_compiler)
#    set(OLD_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
#    set(CMAKE_REQUIRED_FLAGS "-coarray=single")
#    #set(CMAKE_REQUIRED_FLAGS $<$<COMPILE_LANGUAGE:Fortran>:-coarray=single>)
#elseif(gnu_compiler)
#    set(OLD_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
#    set(CMAKE_REQUIRED_FLAGS $<$<COMPILE_LANGUAGE:Fortran>:-ffree-form -fcoarray=single>)
#endif()
#include(CheckFortranSourceCompiles)
#CHECK_Fortran_SOURCE_COMPILES( "program main; implicit none; write(*,*) this_image(); end program" Check_Simple_Coarray_Fortran_Source_Compiles )
#include(CheckFortranSourceRuns)
#check_Fortran_source_runs( "program main; implicit none; write(*,*) this_image(); end program" Check_Simple_Coarray_Fortran_Source_Compiles )
#message ("Check_Simple_Coarray_Fortran_Source_Compiles: ${Check_Simple_Coarray_Fortran_Source_Compiles}")
#if(gnu_compiler OR intel_compiler)
#    set (CMAKE_REQUIRED_FLAGS ${OLD_REQUIRED_FLAGS})
#    unset(OLD_REQUIRED_FLAGS)
#endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# flags configuring the ParaMonte library build
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

option(DLL_ENABLED                             "Define DLL_ENABLED preprocessor flag" OFF)
option(ParaMonte_OBJ_ENABLED                   "Build ParaMonte library object files" ON)
option(ParaMonteTest_OBJ_ENABLED               "Build ParaMonte library test object files, intended only for developers" ON)
option(ParaMonteTest_EXE_ENABLED               "Build ParaMonte library test executables, intended only for developers" OFF)
option(ParaMonteTest_RUN_ENABLED               "Run ParaMonte library developer tests, intended only for developers" OFF)
option(ParaMonteExampleMVNFortran_OBJ_ENABLED  "Build ParaMonte library examples' object files for Fortran " OFF)
option(ParaMonteExampleMVNFortran_EXE_ENABLED  "Build ParaMonte library examples' executables for Fortran" OFF)
option(ParaMonteExampleMVNFortran_RUN_ENABLED  "Run ParaMonte library examples' executables for Fortran" OFF)
option(ParaMonteExampleMVNPython_BLD_ENABLED   "Build ParaMonte library examples' sctips for Python" OFF)
option(ParaMonteExampleMVNPython_RUN_ENABLED   "Run ParaMonte library examples' sctips for Python" OFF)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# define library interface language (all must be lower-case)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set ( INTERFACE_LANGUAGE "c" CACHE STRING "Select the build interface language (c/c++/fortran/matlab/python: all must be lower case)" )

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# define library linking type
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# dynamic:      Use this flag when you have R/Python/MATLAB/Julia code to which you need to link the ParaMonte library dynamically, using DLL files. 
#  static:      Use this flag when you have Fortran/C/C++ code to which you want to link the ParaMonte library statically.
#               You can also link dynamically your Fortran/C/C++ codes using DLL files by specifying LTYPE=dynamic flag instead.

set ( LTYPE "static" CACHE STRING "Select which configuration to build (static/dynamic)" )

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set Fortran/C array allocation resource
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#  ON:          All automatic and stack arrays will be allocated on the heap. Use this when calling ParaMonte from other languages.
# OFF:          All automatic and stack arrays will be allocated on the stack. This could lead to stack overflow when calling DLL from other languages.

option( HEAP_ARRAY_ENABLED "Allocate all arrays and automatic objects on the Heap memory" OFF )

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set interoperability mode
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#  ON:          When you are calling ParaMonte library from any language other than Fortran.
# OFF:          When your objective function to be sampled by ParaMonte library is written in Fortran (as opposed to C).
#               It is also fine to set this flag to TRUE for Fortran application. However, the syntax of the Fortran objective function that
#               is passed to ParaMonte library will have to conform to the rules of C-Fortran interoperability standard,
#               as given in the abstract interface in the ParaMonte source file: ParaMonteLogFunc_mod.f90

option( CFI_ENABLED "Enable C-Fortran interoperability (required when calling ParaMonte from languages other than Fortran)" ON )

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set ParaMonte parallelism
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

include(SetParaMonteParallelism)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Use standardized GNU install directory conventions
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# include(GNUInstallDirs)
# set(mod_dir_tail "${ParaMonte_dist_string}_${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}")
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/${ParaMonte_dist_string}-tests")
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: set required root directories
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

message( STATUS "${pmattn} setting up the root directory..." )

# define ParaMonte_ROOT_DIR as CMAKE_SOURCE_DIR 

message( STATUS "${pmattn} project root directory: ${CMAKE_SOURCE_DIR}" )

# set the ParaMonte source files directory

set(ParaMonte_SRC_DIR
    "${CMAKE_SOURCE_DIR}/src/ParaMonte"
    CACHE PATH "ParaMonte source files directory" )

set(ParaMonteTest_SRC_DIR
    "${CMAKE_SOURCE_DIR}/src/test"
    CACHE PATH "ParaMonte test source files directory" )

# set the ParaMonte interface source files directory and loop over them to ensure their existence

set(ParaMonteInterface_SRC_DIR
    "${CMAKE_SOURCE_DIR}/src/interface"
    CACHE PATH "ParaMonte-interface source files directory" )
set(ParaMonteInterfaceC_SRC_DIR
    "${ParaMonteInterface_SRC_DIR}/C"
    CACHE PATH "ParaMonte-interface C source files directory" )
set(ParaMonteInterfacePython_SRC_DIR
    "${ParaMonteInterface_SRC_DIR}/Python"
    CACHE PATH "ParaMonte-interface Python source files directory" )

message( STATUS "${pmattn} interface source files directory: ${ParaMonteInterface_SRC_DIR}" )
foreach(X IN LISTS  ParaMonte_SRC_DIR 
                    ParaMonteInterface_SRC_DIR 
                    ParaMonteInterfaceC_SRC_DIR 
                    ParaMonteInterfacePython_SRC_DIR )
    if (NOT EXISTS "${X}/")
        message (FATAL_ERROR 
                " \n"
                "${pmfatal}\n"
                "   The source files directory does not exist:\n"
                "   ${X}"
                )
    endif()
endforeach()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: set required build directories
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

message( STATUS "${pmattn} setting up required build directories..." )

unset(PARALLELIZATION_DIR)
if (OMP_ENABLED)
    set(PARALLELIZATION_DIR "${PARALLELIZATION_DIR}omp")
endif()
if (MPI_ENABLED)
    set(PARALLELIZATION_DIR "${PARALLELIZATION_DIR}mpi")
endif()
if (CAF_ENABLED)
    set(PARALLELIZATION_DIR "${PARALLELIZATION_DIR}caf${CAFTYPE}")
endif()
if (NOT DEFINED PARALLELIZATION_DIR)
    set(PARALLELIZATION_DIR "serial")
endif()

# generate build directory

#set(ParaMonte_BLD_DIR "${CMAKE_BINARY_DIR}/${PMCS}/${CMAKE_Fortran_COMPILER_VERSION}/${CMAKE_BUILD_TYPE}/${LTYPE}/${PARALLELIZATION_DIR}")
set(ParaMonte_BLD_DIR "${CMAKE_BINARY_DIR}")
message( STATUS "${pmattn} ParaMonte_BLD_DIR: ${ParaMonte_BLD_DIR}" )
if (NOT EXISTS "${ParaMonte_BLD_DIR}/")
    file(MAKE_DIRECTORY "${ParaMonte_BLD_DIR}")
endif()
message( STATUS "${pmnote} all generated build files will be stored at: ${ParaMonte_BLD_DIR}" )

# generate object files directory

set(ParaMonte_OBJ_DIR "${ParaMonte_BLD_DIR}/obj")
message( STATUS "${pmattn} ParaMonte_OBJ_DIR: ${ParaMonte_OBJ_DIR}" )
if (NOT EXISTS "${ParaMonte_OBJ_DIR}/")
    file(MAKE_DIRECTORY "${ParaMonte_OBJ_DIR}")
endif()

# generate modules files directory

set(ParaMonte_MOD_DIR "${ParaMonte_BLD_DIR}/mod")
message( STATUS "${pmattn} ParaMonte_MOD_DIR: ${ParaMonte_MOD_DIR}" )
if (NOT EXISTS "${ParaMonte_MOD_DIR}/")
    file(MAKE_DIRECTORY "${ParaMonte_MOD_DIR}")
endif()

# generate lib files directory

set(ParaMonte_LIB_DIR "${ParaMonte_BLD_DIR}/lib")
message( STATUS "${pmattn} ParaMonte_LIB_DIR: ${ParaMonte_LIB_DIR}" )
if (NOT EXISTS "${ParaMonte_LIB_DIR}/")
    file(MAKE_DIRECTORY "${ParaMonte_LIB_DIR}")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: set required test build directories
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# set up test build directory

set(ParaMonteTest_BLD_DIR "${ParaMonte_BLD_DIR}/test")
message( STATUS "${pmattn} ParaMonteTest_BLD_DIR: ${ParaMonteTest_BLD_DIR}" )
if (NOT EXISTS "${ParaMonteTest_BLD_DIR}/")
    file(MAKE_DIRECTORY "${ParaMonteTest_BLD_DIR}")
endif()

# set up test object directory

set(ParaMonteTest_OBJ_DIR "${ParaMonteTest_BLD_DIR}/obj")
message( STATUS "${pmattn} ParaMonteTest_OBJ_DIR: ${ParaMonteTest_OBJ_DIR}" )
if (NOT EXISTS "${ParaMonteTest_OBJ_DIR}/")
    file(MAKE_DIRECTORY "${ParaMonteTest_OBJ_DIR}")
endif()

# set up test module directory

set(ParaMonteTest_MOD_DIR "${ParaMonteTest_BLD_DIR}/mod")
message( STATUS "${pmattn} ParaMonteTest_MOD_DIR: ${ParaMonteTest_MOD_DIR}" )
if (NOT EXISTS "${ParaMonteTest_MOD_DIR}/")
    file(MAKE_DIRECTORY "${ParaMonteTest_MOD_DIR}")
endif()

# set up test bin directory

#set(ParaMonteTest_BIN_DIR "${ParaMonteTest_BLD_DIR}")
set(ParaMonteTest_BIN_DIR "${ParaMonteTest_BLD_DIR}/bin")
message( STATUS "${pmattn} ParaMonteTest_BIN_DIR: ${ParaMonteTest_BIN_DIR}" )
if (NOT EXISTS "${ParaMonteTest_BIN_DIR}/")
    file(MAKE_DIRECTORY "${ParaMonteTest_BIN_DIR}")
endif()

# set up test bin directory

set(ParaMonteTest_INPUT_DIR "${ParaMonteTest_BIN_DIR}/input")
message( STATUS "${pmattn} ParaMonteTest_INPUT_DIR: ${ParaMonteTest_INPUT_DIR}" )
if (NOT EXISTS "${ParaMonteTest_INPUT_DIR}/")
    file(MAKE_DIRECTORY "${ParaMonteTest_INPUT_DIR}")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: generate ParaMonte library build directories, object files, and dynamic libraries
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#include_directories(BEFORE ${ParaMonte_LIB_DIR})

# message( STATUS " ${pmattn} CMAKE_EXECUTABLE_SUFFIX_Fortran: ${CMAKE_EXECUTABLE_SUFFIX_Fortran}" )

# Define the executable name
set(ParaMonteTest_EXE_NAME testParaMonte)

# Define the library name
if (NOT DEFINED PMLIB)

    if (HEAP_ARRAY_ENABLED)
        set(MEMORY_ALLOCATION heap)
    else()
        set(MEMORY_ALLOCATION stack)
    endif()

    if (WIN32)
        set(ARCHITECTURE ${PLATFORM})
    else()
        EXECUTE_PROCESS( COMMAND uname -m COMMAND tr -d "\n" OUTPUT_VARIABLE ARCHITECTURE )
        string(REPLACE "x86_64" "x64" ARCHITECTURE ${ARCHITECTURE})
    endif()
    string(TOLOWER "${CMAKE_SYSTEM_NAME}" PLATFORM)

    if (${INTERFACE_LANGUAGE} STREQUAL "c++" OR ${INTERFACE_LANGUAGE} STREQUAL "C++")
        set(LANG_ABBR "cpp")
    else()
        set(LANG_ABBR ${INTERFACE_LANGUAGE})
    endif()
    set(PMLIB paramonte_${LANG_ABBR}_${PLATFORM}_${ARCHITECTURE}_${PMCS}_${BTYPE}_${LTYPE}_${MEMORY_ALLOCATION} )

    if (CAF_ENABLED)
        set(PMLIB ${PMLIB}_caf${CAFTYPE})
    endif()

    if (MPI_ENABLED)
        set(PMLIB ${PMLIB}_mpi)
    endif()

    if (OMP_ENABLED)
        set(PMLIB ${PMLIB}_omp)
    endif()

    #if (MT_ENABLED)
    #    set(PMLIB ${PMLIB}_mt)
    #endif()

endif()

# Have the .mod files placed in the lib folder
# set(CMAKE_Fortran_MODULE_DIRECTORY "${ParaMonte_MOD_DIR}")

# The source for the ParaMonte library and its binaries in the ParaMonte_OBJ_DIR folder
#add_subdirectory(${ParaMonte_SRC_DIR} ${ParaMonte_OBJ_DIR})
add_subdirectory(${ParaMonte_SRC_DIR} "${ParaMonte_OBJ_DIR}")

# The source for the ParaMonte library tests and its binaries in the ParaMonte_BIN_DIR folder
add_subdirectory(${ParaMonteTest_SRC_DIR} ${ParaMonteTest_OBJ_DIR})

#message(STATUS "CMAKE_Fortran_LINK_EXECUTABLE: ${CMAKE_Fortran_LINK_EXECUTABLE}")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Add a distclean target to the Makefile
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

ADD_CUSTOM_TARGET(dclean COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake)
ADD_CUSTOM_TARGET(distclean COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake)
