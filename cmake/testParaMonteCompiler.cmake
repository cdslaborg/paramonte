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

if(CMAKE_Fortran_COMPILER_FORCED)
    #   The compiler configuration was forced by the user.
    #   Assume the user has configured all compiler information.
    set(CMAKE_Fortran_COMPILER_WORKS TRUE)
    return()
endif()

include(CMakeTestCompilerCommon)

#   Remove any cached result from an older CMake version.
#   We now store this in CMakeFortranCompiler.cmake.
unset(CMAKE_Fortran_COMPILER_WORKS CACHE)

#   This file is used by EnableLanguage in cmGlobalGenerator to
#   determine that the selected Fortran compiler can actually compile
#   and link the most basic of programs.   If not, a fatal error
#   is set and cmake stops processing commands and will not generate
#   any makefiles or projects.

if(NOT CMAKE_Fortran_COMPILER_WORKS)
    PrintTestCompilerStatus("Fortran" "")
    file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompiler.f
        "
            PROGRAM TESTFortran
            PRINT *, 'Hello'
            END
        ")
    try_compile(
                CMAKE_Fortran_COMPILER_WORKS
                ${CMAKE_BINARY_DIR}
                ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompiler.f
                OUTPUT_VARIABLE program_output
                )
    # Move result from cache to normal variable.
    string(REPLACE "\n" "\n " program_output "${program_output}")
    set(CMAKE_Fortran_COMPILER_WORKS ${CMAKE_Fortran_COMPILER_WORKS})
    unset(CMAKE_Fortran_COMPILER_WORKS CACHE)
    set(FORTRAN_TEST_WAS_RUN 1)
endif()

if(NOT CMAKE_Fortran_COMPILER_WORKS)

    PrintTestCompilerStatus("Fortran" " -- broken")
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
        "Determining if the Fortran compiler yielded the following output:\n${program_output}\n\n")
    message(FATAL_ERROR
            "\n"
            "${pmfatal} \n"
            "${pmfatal} \n"
            "${pmfatal} ${pmfatal} The Fortran compiler\n"
            "${pmfatal} \n"
            "${pmfatal}        \"${CMAKE_Fortran_COMPILER}\"\n"
            "${pmfatal} \n"
            "${pmfatal}    is not able to compile a simple test program.\n"
            "${pmfatal}    It fails with the following output:\n"
            "${pmfatal} \n"
            "               ${program_output} "
            "${pmfatal} \n"
            "${pmfatal}     CMake will not be able to correctly generate this project."
            "\n"
            )
else()

    if(FORTRAN_TEST_WAS_RUN)
        PrintTestCompilerStatus("Fortran" " -- works")
        file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
            " Determining if the Fortran compiler works passed with the following output\n"
            " \n"
            " ${program_output}\n"
            " \n"
            )
    endif()

    # Try to identify the ABI and configure it into CMakeFortranCompiler.cmake
    include(${CMAKE_ROOT}/Modules/CMakeDetermineCompilerABI.cmake)
    CMAKE_DETERMINE_COMPILER_ABI(Fortran ${CMAKE_ROOT}/Modules/CMakeFortranCompilerABI.F)

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Test for Fortran 90 support by using an F2008-specific construct.
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if (DEFINED CMAKE_Fortran_COMPILER_SUPPORTS_F2008 AND CMAKE_Fortran_COMPILER_SUPPORTS_F2008)
        message(NOTICE "${pmattn} Does \"${CMAKE_Fortran_COMPILER}\" support Fortran 2008? ${CMAKE_Fortran_COMPILER_SUPPORTS_F2008} -- ${BoldGreen}PASSED${ResetColor}")
    else()
        set(FortranExamplePath "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompliance2008.F90")
        message(NOTICE "${pmattn} Generating compiler test file: \"${FortranExamplePath}\"")
        file(WRITE "${FortranExamplePath}"
        "
            module test_mod
                integer :: mv_test = 0
                interface exec
                module subroutine exec()
                end subroutine
                end interface
            end module
            submodule (test_mod) test_smod
                integer :: smv_test = 0
            contains
                module procedure exec
                end procedure
            end submodule
            program testFortranCompliance2008
                use test_mod
                block
                    integer :: i, A(5) = 0
                    integer :: B(size(A))
                    do concurrent(i = 1:size(A))
                        B(i) = A(i)
                    end do
                end block
            end program testFortranCompliance2008
        ")

        try_compile(CMAKE_Fortran_COMPILER_SUPPORTS_F2008 ${CMAKE_BINARY_DIR}
                    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompliance2008.F90
                    OUTPUT_VARIABLE program_output
                    )

        set(CMAKE_Fortran_COMPILER_SUPPORTS_F2008 "${CMAKE_Fortran_COMPILER_SUPPORTS_F2008}" CACHE STRING "CMAKE_Fortran_COMPILER_SUPPORTS_F2008?" FORCE)
        if(CMAKE_Fortran_COMPILER_SUPPORTS_F2008)
            message(NOTICE "${pmattn} Does \"${CMAKE_Fortran_COMPILER}\" support Fortran 2008? ${CMAKE_Fortran_COMPILER_SUPPORTS_F2008} -- ${BoldGreen}PASSED${ResetColor}")
            file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
                "Determining if the Fortran compiler supports Fortran 2008 passed with "
                "the following output:\n${program_output}\n\n")
        else()
            message(NOTICE "${pmattn} Does \"${CMAKE_Fortran_COMPILER}\" support Fortran 2008? ${CMAKE_Fortran_COMPILER_SUPPORTS_F2008} -- ${BoldRed}FAILED${ResetColor}")
            file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
                "Determining if the Fortran compiler supports Fortran 2008 failed with "
                "the following output:\n${program_output}\n\n")
            message(FATAL_ERROR "${program_output}")
        endif()
    endif()

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Test for Coarray Fortran (CAF) 2008 support by using an Coarray-specific construct.
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if(CAF_ENABLED)

        if(${csid_is_intel})
           set(OLD_REQUIRED_FLAGS ${CMAKE_Fortran_FLAGS})
           set(CMAKE_Fortran_FLAGS "-coarray=single")
        elseif(${csid_is_gnu})
           set(OLD_REQUIRED_FLAGS ${CMAKE_Fortran_FLAGS})
           set(CMAKE_Fortran_FLAGS "-ffree-form -fcoarray=single")
        endif()

        message(NOTICE "${pmattn} Checking whether ${CMAKE_Fortran_COMPILER} supports Coarray Fortran (CAF) 2008")

        file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompilerCAF.F90
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
                    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompilerCAF.F90
                    OUTPUT_VARIABLE program_output
                    )

        if(${csid_is_gnu} OR ${csid_is_intel})
           set (CMAKE_Fortran_FLAGS ${OLD_REQUIRED_FLAGS})
           unset(OLD_REQUIRED_FLAGS)
        endif()

        if(CMAKE_Fortran_COMPILER_SUPPORTS_CAF)
            message(NOTICE "${pmattn} Checking whether ${CMAKE_Fortran_COMPILER} supports Coarray Fortran (CAF) 2008 -- yes")
            file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
                "Determining if the Fortran compiler supports Coarray Fortran 2008 passed with "
                "the following output:\n${program_output}\n\n"
                )
            set(CMAKE_Fortran_COMPILER_SUPPORTS_CAF 1)
        else()
            message(NOTICE "${pmattn} Checking whether ${CMAKE_Fortran_COMPILER} supports Coarray Fortran (CAF) 2008 -- no")
            file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
                "Determining if the Fortran compiler supports Coarray Fortran (CAF) 2008 failed with "
                "the following output:\n${program_output}\n\n"
                )
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
