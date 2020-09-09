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

if(CMAKE_Fortran_COMPILER_FORCED)
    # The compiler configuration was forced by the user.
    # Assume the user has configured all compiler information.
    set(CMAKE_Fortran_COMPILER_WORKS TRUE)
    return()
endif()

include(CMakeTestCompilerCommon)

# Remove any cached result from an older CMake version.
# We now store this in CMakeFortranCompiler.cmake.
unset(CMAKE_Fortran_COMPILER_WORKS CACHE)

# This file is used by EnableLanguage in cmGlobalGenerator to
# determine that the selected Fortran compiler can actually compile
# and link the most basic of programs.   If not, a fatal error
# is set and cmake stops processing commands and will not generate
# any makefiles or projects.
if(NOT CMAKE_Fortran_COMPILER_WORKS)
    PrintTestCompilerStatus("Fortran" "")
    file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompiler.f "
            PROGRAM TESTFortran
            PRINT *, 'Hello'
            END
        ")
    try_compile(CMAKE_Fortran_COMPILER_WORKS ${CMAKE_BINARY_DIR}
        ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompiler.f
        OUTPUT_VARIABLE OUTPUT)
    # Move result from cache to normal variable.
    set(CMAKE_Fortran_COMPILER_WORKS ${CMAKE_Fortran_COMPILER_WORKS})
    unset(CMAKE_Fortran_COMPILER_WORKS CACHE)
    set(FORTRAN_TEST_WAS_RUN 1)
endif()

if(NOT CMAKE_Fortran_COMPILER_WORKS)
    PrintTestCompilerStatus("Fortran" "  -- broken")
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
        "Determining if the Fortran compiler works/fails with the following output:\n${OUTPUT}\n\n")
    string(REPLACE "\n" "\n  " _output "${OUTPUT}")
    message(FATAL_ERROR "The Fortran compiler\n  \"${CMAKE_Fortran_COMPILER}\"\n"
                        "is not able to compile a simple test program.\nIt fails "
                        "with the following output:\n  ${_output}\n\n"
                        "CMake will not be able to correctly generate this project.")
else()
    if(FORTRAN_TEST_WAS_RUN)
        PrintTestCompilerStatus("Fortran" "  -- works")
        file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
            "Determining if the Fortran compiler works passed with "
            "the following output:\n${OUTPUT}\n\n")
    endif()

    # Try to identify the ABI and configure it into CMakeFortranCompiler.cmake
    include(${CMAKE_ROOT}/Modules/CMakeDetermineCompilerABI.cmake)
    CMAKE_DETERMINE_COMPILER_ABI(Fortran ${CMAKE_ROOT}/Modules/CMakeFortranCompilerABI.F)

    #####################################################################################
    # Test for Fortran 90 support by using an F2008-specific construct.
    if(NOT DEFINED CMAKE_Fortran_COMPILER_SUPPORTS_F2008)
        message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports Fortran 2008")
        file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompilerF2008.f90
        "
            PROGRAM TESTFortran2008
            block; integer stop ; stop = 1 ; do while ( stop == 0 ) ; end do; end block;
            END PROGRAM TESTFortran2008
        ")
        try_compile(CMAKE_Fortran_COMPILER_SUPPORTS_F2008 ${CMAKE_BINARY_DIR}
                    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompilerF2008.f90
                    OUTPUT_VARIABLE OUTPUT)
        if(CMAKE_Fortran_COMPILER_SUPPORTS_F2008)
            message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports Fortran 2008 -- yes")
            file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
                "Determining if the Fortran compiler supports Fortran 2008 passed with "
                "the following output:\n${OUTPUT}\n\n")
            set(CMAKE_Fortran_COMPILER_SUPPORTS_F2008 1)
        else()
            message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports Fortran 2008 -- no")
            file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
                "Determining if the Fortran compiler supports Fortran 2008 failed with "
                "the following output:\n${OUTPUT}\n\n")
            set(CMAKE_Fortran_COMPILER_SUPPORTS_F2008 0)
        endif()
        unset(CMAKE_Fortran_COMPILER_SUPPORTS_F2008 CACHE)
    endif()

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Test for Coarray Fortran (CAF) 2008 support by using an Coarray-specific construct.
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if( CAF_ENABLED OR intel_compiler )

        if(intel_compiler)
           set(OLD_REQUIRED_FLAGS ${CMAKE_Fortran_FLAGS})
           set(CMAKE_Fortran_FLAGS "-coarray=single")
           #set(CMAKE_Fortran_FLAGS $<$<COMPILE_LANGUAGE:Fortran>:-coarray=single>)
        elseif(gnu_compiler)
           set(OLD_REQUIRED_FLAGS ${CMAKE_Fortran_FLAGS})
           set(CMAKE_Fortran_FLAGS "-ffree-form -fcoarray=single")
        endif()

        message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports Coarray Fortran (CAF) 2008")
        file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompilerCAF.f90
            "
            program testFortranCompilerCAF
                block
                    implicit none
                    integer :: imageID, imageCount
                    imageID     = this_image()
                    imageCount  = num_images()
                end block
            end program testFortranCompilerCAF
            ")
        try_compile(CMAKE_Fortran_COMPILER_SUPPORTS_CAF ${CMAKE_BINARY_DIR}
                    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompilerCAF.f90
                    OUTPUT_VARIABLE OUTPUT)

        if(gnu_compiler OR intel_compiler)
           set (CMAKE_Fortran_FLAGS ${OLD_REQUIRED_FLAGS})
           unset(OLD_REQUIRED_FLAGS)
        endif()

        if(CMAKE_Fortran_COMPILER_SUPPORTS_CAF)
            message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports Coarray Fortran (CAF) 2008 -- yes")
            file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
                "Determining if the Fortran compiler supports Coarray Fortran 2008 passed with "
                "the following output:\n${OUTPUT}\n\n")
            set(CMAKE_Fortran_COMPILER_SUPPORTS_CAF 1)
        else()
            message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports Coarray Fortran (CAF) 2008 -- no")
            file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
                "Determining if the Fortran compiler supports Coarray Fortran (CAF) 2008 failed with "
                "the following output:\n${OUTPUT}\n\n")
            set(CMAKE_Fortran_COMPILER_SUPPORTS_CAF 0)
        endif()

        unset(CMAKE_Fortran_COMPILER_SUPPORTS_CAF CACHE)

    endif()

    # Re-configure to save learned information.
    configure_file(
        ${CMAKE_ROOT}/Modules/CMakeFortranCompiler.cmake.in
        ${CMAKE_PLATFORM_INFO_DIR}/CMakeFortranCompiler.cmake
        @ONLY
        )
    include(${CMAKE_PLATFORM_INFO_DIR}/CMakeFortranCompiler.cmake)

    if(CMAKE_Fortran_SIZEOF_DATA_PTR)
        foreach(f ${CMAKE_Fortran_ABI_FILES})
            include(${f})
        endforeach()
        unset(CMAKE_Fortran_ABI_FILES)
    endif()
endif()
