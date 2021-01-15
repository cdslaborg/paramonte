!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!   MIT License
!!!!
!!!!   ParaMonte: plain powerful parallel Monte Carlo library.
!!!!
!!!!   Copyright (C) 2012-present, The Computational Data Science Lab
!!!!
!!!!   This file is part of the ParaMonte library.
!!!!
!!!!   Permission is hereby granted, free of charge, to any person obtaining a
!!!!   copy of this software and associated documentation files (the "Software"),
!!!!   to deal in the Software without restriction, including without limitation
!!!!   the rights to use, copy, modify, merge, publish, distribute, sublicense,
!!!!   and/or sell copies of the Software, and to permit persons to whom the
!!!!   Software is furnished to do so, subject to the following conditions:
!!!!
!!!!   The above copyright notice and this permission notice shall be
!!!!   included in all copies or substantial portions of the Software.
!!!!
!!!!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!!!!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!!!!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!!!!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
!!!!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
!!!!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
!!!!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!!!!
!!!!   ACKNOWLEDGMENT
!!!!
!!!!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!!!!   As per the ParaMonte library license agreement terms, if you use any parts of
!!!!   this library for any purposes, kindly acknowledge the use of ParaMonte in your
!!!!   work (education/research/industry/development/...) by citing the ParaMonte
!!!!   library as described on this page:
!!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This module contains tests of the module 
!> [ParaDRAM_RefinedChain_mod](@ref paradram_refinedchain_mod),
!> [ParaDISE_RefinedChain_mod](@ref paradise_refinedchain_mod),
!> [ParaNest_RefinedChain_mod](@ref paranest_refinedchain_mod).
!>  \author Amir Shahmoradi

#if defined PARADRAM

#define ParaXXXX_RefinedChain_mod ParaDRAM_RefinedChain_mod

#elif defined PARADISE

#define ParaXXXX_RefinedChain_mod ParaDISE_RefinedChain_mod

#elif defined PARANEST

#define ParaXXXX_RefinedChain_mod ParaNest_RefinedChain_mod

#endif

    use Test_mod, only: Test_type, getLogFuncMVN
    use ParaXXXX_RefinedChain_mod
    implicit none

    private
    public :: test_RefinedChain

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_RefinedChain()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_getSkip4NewSampleSize_1, "test_getSkip4NewSampleSize_1")
        call Test%run(test_getSkip4NewSampleSize_2, "test_getSkip4NewSampleSize_2")
        call Test%run(test_getSkip4NewSampleSize_3, "test_getSkip4NewSampleSize_3")
        call Test%run(test_getSkip4NewSampleSize_4, "test_getSkip4NewSampleSize_4")
        call Test%finalize()
    end subroutine test_RefinedChain

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getSkip4NewSampleSize_1() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        integer(IK), parameter      :: oldSampleSize = 10_IK
        integer(IK), parameter      :: newSampleSize = 3_IK
        integer(IK), parameter      :: skip4NewSampleSize_ref = 4_IK
        integer(IK)                 :: skip4NewSampleSize
        skip4NewSampleSize = getSkip4NewSampleSize(oldSampleSize,newSampleSize)
        assertion = skip4NewSampleSize == skip4NewSampleSize_ref
        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(*,"(10(g0,:,', '))")
            write(*,"(10(g0,:,', '))") "skip4NewSampleSize_ref  ", skip4NewSampleSize_ref
            write(*,"(10(g0,:,', '))") "skip4NewSampleSize      ", skip4NewSampleSize
            write(*,"(10(g0,:,', '))")
            return
        end if
        ! LCOV_EXCL_STOP
    end function test_getSkip4NewSampleSize_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getSkip4NewSampleSize_2() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        integer(IK), parameter      :: oldSampleSize = 6_IK
        integer(IK), parameter      :: newSampleSize = 3_IK
        integer(IK), parameter      :: skip4NewSampleSize_ref = 2_IK
        integer(IK)                 :: skip4NewSampleSize
        skip4NewSampleSize = getSkip4NewSampleSize(oldSampleSize,newSampleSize)
        assertion = skip4NewSampleSize == skip4NewSampleSize_ref
        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(*,"(10(g0,:,', '))")
            write(*,"(10(g0,:,', '))") "skip4NewSampleSize_ref  ", skip4NewSampleSize_ref
            write(*,"(10(g0,:,', '))") "skip4NewSampleSize      ", skip4NewSampleSize
            write(*,"(10(g0,:,', '))")
            return
        end if
        ! LCOV_EXCL_STOP
    end function test_getSkip4NewSampleSize_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getSkip4NewSampleSize_3() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        integer(IK), parameter      :: oldSampleSize = 3_IK
        integer(IK), parameter      :: newSampleSize = 3_IK
        integer(IK), parameter      :: skip4NewSampleSize_ref = 1_IK
        integer(IK)                 :: skip4NewSampleSize
        skip4NewSampleSize = getSkip4NewSampleSize(oldSampleSize,newSampleSize)
        assertion = skip4NewSampleSize == skip4NewSampleSize_ref
        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(*,"(10(g0,:,', '))")
            write(*,"(10(g0,:,', '))") "skip4NewSampleSize_ref  ", skip4NewSampleSize_ref
            write(*,"(10(g0,:,', '))") "skip4NewSampleSize      ", skip4NewSampleSize
            write(*,"(10(g0,:,', '))")
            return
        end if
        ! LCOV_EXCL_STOP
    end function test_getSkip4NewSampleSize_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getSkip4NewSampleSize_4() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        integer(IK), parameter      :: oldSampleSize = 2_IK
        integer(IK), parameter      :: newSampleSize = 3_IK
        integer(IK)                 :: skip4NewSampleSize
        skip4NewSampleSize = getSkip4NewSampleSize(oldSampleSize,newSampleSize)
        assertion = skip4NewSampleSize < 0_IK
        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(*,"(10(g0,:,', '))")
            write(*,"(10(g0,:,', '))") "skip4NewSampleSize      ", skip4NewSampleSize
            write(*,"(10(g0,:,', '))")
            return
        end if
        ! LCOV_EXCL_STOP
    end function test_getSkip4NewSampleSize_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ParaXXXX_RefinedChain_mod
