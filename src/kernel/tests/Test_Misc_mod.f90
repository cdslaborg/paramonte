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
!!!!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This module contains tests of the module [Misc_mod](@ref misc_mod).
!>  \author Amir Shahmoradi

module Test_Misc_mod

    use Misc_mod
    use Err_mod, only: Err_type
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_Misc

    type(Test_type) :: Test

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_Misc()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_arth_IK_1, "test_arth_IK_1")
        call Test%run(test_arth_RK_1, "test_arth_RK_1")
        call Test%run(test_arth_RK_2, "test_arth_RK_2")
        call Test%run(test_swap_IK_1, "test_swap_IK_1")
        call Test%run(test_swap_RK_1, "test_swap_RK_1")
        call Test%run(test_swap_CK_1, "test_swap_RK_1")
        call Test%run(test_swap_SPI_1, "test_swap_SPI_1")
        call Test%run(test_swap_DPI_1, "test_swap_DPI_1")
        call Test%run(test_swap_SPR_1, "test_swap_SPR_1")
        call Test%run(test_swap_DPR_1, "test_swap_DPR_1")
        call Test%run(test_swap_SPC_1, "test_swap_SPC_1")
        call Test%run(test_swap_DPC_1, "test_swap_DPC_1")
        call Test%run(test_zroots_unity_1, "test_zroots_unity_1")
        call Test%run(test_copyArray_IK_1, "test_copyArray_IK_1")
        call Test%run(test_copyArray_IK_2, "test_copyArray_IK_2")
        call Test%run(test_copyArray_RK_1, "test_copyArray_RK_1")
        call Test%run(test_copyArray_RK_2, "test_copyArray_RK_2")
        call Test%run(test_resizeVector_RK_1, "test_resizeVector_RK_1")
        call Test%run(test_masked_swap_SPR_1, "test_masked_swap_SPR_1")
        call Test%run(test_masked_swap_SPR_2, "test_masked_swap_SPR_2")
        call Test%run(test_masked_swap_SPRV_1, "test_masked_swap_SPRV_1")
        call Test%run(test_masked_swap_SPRM_1, "test_masked_swap_SPRM_1")
        call Test%run(test_findUnique_1, "test_findUnique_1")
        call Test%finalize()
    end subroutine test_Misc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_swap_IK_1() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                 :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: vecLen = 3_IK
        integer(IK) , parameter :: Vector1_ref(vecLen) = [(+i,i=1,vecLen)]
        integer(IK) , parameter :: Vector2_ref(vecLen) = [(-i,i=1,vecLen)]
        integer(IK)             :: Vector1(vecLen)
        integer(IK)             :: Vector2(vecLen)
        Vector1 = Vector1_ref
        Vector2 = Vector2_ref
        call swap_IK(Vector1,Vector2)
        assertion = all(Vector1==Vector2_ref) .and. all(Vector2==Vector1_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1_ref =", Vector1_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2     =", Vector2
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2_ref =", Vector2_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1     =", Vector1
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_swap_IK_1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_swap_SPI_1() result(assertion)
        use Constants_mod, only: SPI
        implicit none
        logical                  :: assertion
        integer(SPI)             :: i
        integer(SPI) , parameter :: vecLen = 3_SPI
        integer(SPI) , parameter :: Vector1_ref(vecLen) = [(+i,i=1,vecLen)]
        integer(SPI) , parameter :: Vector2_ref(vecLen) = [(-i,i=1,vecLen)]
        integer(SPI)             :: Vector1(vecLen)
        integer(SPI)             :: Vector2(vecLen)
        Vector1 = Vector1_ref
        Vector2 = Vector2_ref
        call swap_SPI(Vector1,Vector2)
        assertion = all(Vector1==Vector2_ref) .and. all(Vector2==Vector1_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1_ref =", Vector1_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2     =", Vector2
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2_ref =", Vector2_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1     =", Vector1
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_swap_SPI_1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_swap_DPI_1() result(assertion)
        use Constants_mod, only: DPI
        implicit none
        logical                  :: assertion
        integer(DPI)             :: i
        integer(DPI) , parameter :: vecLen = 3
        integer(DPI) , parameter :: Vector1_ref(vecLen) = [(+i,i=1,vecLen)]
        integer(DPI) , parameter :: Vector2_ref(vecLen) = [(-i,i=1,vecLen)]
        integer(DPI)             :: Vector1(vecLen)
        integer(DPI)             :: Vector2(vecLen)
        Vector1 = Vector1_ref
        Vector2 = Vector2_ref
        call swap_DPI(Vector1,Vector2)
        assertion = all(Vector1==Vector2_ref) .and. all(Vector2==Vector1_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1_ref =", Vector1_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2     =", Vector2
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2_ref =", Vector2_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1     =", Vector1
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_swap_DPI_1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_swap_RK_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: vecLen = 3_IK
        real(RK)    , parameter :: Vector1_ref(vecLen) = [(real(+i,RK),i=1,vecLen)]
        real(RK)    , parameter :: Vector2_ref(vecLen) = [(real(-i,RK),i=1,vecLen)]
        real(RK)                :: Vector1(vecLen)
        real(RK)                :: Vector2(vecLen)
        Vector1 = Vector1_ref
        Vector2 = Vector2_ref
        call swap_RK(Vector1,Vector2)
        assertion = all(Vector1==Vector2_ref) .and. all(Vector2==Vector1_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1_ref =", Vector1_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2     =", Vector2
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2_ref =", Vector2_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1     =", Vector1
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_swap_RK_1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_swap_SPR_1() result(assertion)
        use Constants_mod, only: IK, SPR
        implicit none
        logical                 :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: vecLen = 3_IK
        real(SPR)   , parameter :: Vector1_ref(vecLen) = [(real(+i,SPR),i=1,vecLen)]
        real(SPR)   , parameter :: Vector2_ref(vecLen) = [(real(-i,SPR),i=1,vecLen)]
        real(SPR)               :: Vector1(vecLen)
        real(SPR)               :: Vector2(vecLen)
        Vector1 = Vector1_ref
        Vector2 = Vector2_ref
        call swap_SPR(Vector1,Vector2)
        assertion = all(Vector1==Vector2_ref) .and. all(Vector2==Vector1_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1_ref =", Vector1_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2     =", Vector2
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2_ref =", Vector2_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1     =", Vector1
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_swap_SPR_1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_swap_DPR_1() result(assertion)
        use Constants_mod, only: IK, DPR
        implicit none
        logical                 :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: vecLen = 3_IK
        real(DPR)   , parameter :: Vector1_ref(vecLen) = [(real(+i,DPR),i=1,vecLen)]
        real(DPR)   , parameter :: Vector2_ref(vecLen) = [(real(-i,DPR),i=1,vecLen)]
        real(DPR)               :: Vector1(vecLen)
        real(DPR)               :: Vector2(vecLen)
        Vector1 = Vector1_ref
        Vector2 = Vector2_ref
        call swap_DPR(Vector1,Vector2)
        assertion = all(Vector1==Vector2_ref) .and. all(Vector2==Vector1_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1_ref =", Vector1_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2     =", Vector2
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2_ref =", Vector2_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1     =", Vector1
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_swap_DPR_1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_swap_CK_1() result(assertion)
        use Constants_mod, only: IK, RK, CK
        implicit none
        logical                 :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: vecLen = 3_IK
        complex(CK) , parameter :: Vector1_ref(vecLen) = [(cmplx(+i,0.,RK),i=1,vecLen)]
        complex(CK) , parameter :: Vector2_ref(vecLen) = [(cmplx(-i,0.,RK),i=1,vecLen)]
        complex(CK)             :: Vector1(vecLen)
        complex(CK)             :: Vector2(vecLen)
        Vector1 = Vector1_ref
        Vector2 = Vector2_ref
        call swap_CK(Vector1,Vector2)
        assertion = all(Vector1==Vector2_ref) .and. all(Vector2==Vector1_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1_ref =", Vector1_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2     =", Vector2
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2_ref =", Vector2_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1     =", Vector1
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_swap_CK_1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_swap_SPC_1() result(assertion)
        use Constants_mod, only: IK, SPR, SPC
        implicit none
        logical                 :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: vecLen = 3_IK
        complex(SPC), parameter :: Vector1_ref(vecLen) = [(cmplx(+i,0.,SPR),i=1,vecLen)]
        complex(SPC), parameter :: Vector2_ref(vecLen) = [(cmplx(-i,0.,SPR),i=1,vecLen)]
        complex(SPC)            :: Vector1(vecLen)
        complex(SPC)            :: Vector2(vecLen)
        Vector1 = Vector1_ref
        Vector2 = Vector2_ref
        call swap_SPC(Vector1,Vector2)
        assertion = all(Vector1==Vector2_ref) .and. all(Vector2==Vector1_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1_ref =", Vector1_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2     =", Vector2
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2_ref =", Vector2_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1     =", Vector1
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_swap_SPC_1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_swap_DPC_1() result(assertion)
        use Constants_mod, only: IK, DPR, DPC
        implicit none
        logical                 :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: vecLen = 3_IK
        complex(DPC), parameter :: Vector1_ref(vecLen) = [(cmplx(+i,0.,DPR),i=1,vecLen)]
        complex(DPC), parameter :: Vector2_ref(vecLen) = [(cmplx(-i,0.,DPR),i=1,vecLen)]
        complex(DPC)            :: Vector1(vecLen)
        complex(DPC)            :: Vector2(vecLen)
        Vector1 = Vector1_ref
        Vector2 = Vector2_ref
        call swap_DPC(Vector1,Vector2)
        assertion = all(Vector1==Vector2_ref) .and. all(Vector2==Vector1_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1_ref =", Vector1_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2     =", Vector2
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2_ref =", Vector2_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1     =", Vector1
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_swap_DPC_1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_masked_swap_SPR_1() result(assertion)
        use Constants_mod, only: IK, SPR
        implicit none
        logical                 :: assertion
        logical     , parameter :: mask = .true.
        real(SPR)   , parameter :: scalar1_ref = 1._SPR
        real(SPR)   , parameter :: scalar2_ref = 2._SPR
        real(SPR)               :: scalar1
        real(SPR)               :: scalar2
        scalar1 = scalar1_ref
        scalar2 = scalar2_ref
        call masked_swap_SPR(scalar1,scalar2,mask)
        assertion = scalar1==scalar2_ref .and. scalar2==scalar1_ref
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "scalar1_ref =", scalar1_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "scalar2     =", scalar2
            write(Test%outputUnit,"(*(g0,:,' '))") "scalar2_ref =", scalar2_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "scalar1     =", scalar1
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_masked_swap_SPR_1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_masked_swap_SPR_2() result(assertion)
        use Constants_mod, only: IK, SPR
        implicit none
        logical                 :: assertion
        logical     , parameter :: mask = .false.
        real(SPR)   , parameter :: scalar1_ref = 1._SPR
        real(SPR)   , parameter :: scalar2_ref = 2._SPR
        real(SPR)               :: scalar1
        real(SPR)               :: scalar2
        scalar1 = scalar1_ref
        scalar2 = scalar2_ref
        call masked_swap_SPR(scalar1,scalar2,mask)
        assertion = scalar1==scalar1_ref .and. scalar2==scalar2_ref
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "scalar1_ref =", scalar1_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "scalar1     =", scalar1
            write(Test%outputUnit,"(*(g0,:,' '))") "scalar2_ref =", scalar2_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "scalar2     =", scalar2
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_masked_swap_SPR_2

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_masked_swap_SPRV_1() result(assertion)
        use Constants_mod, only: IK, SPR
        implicit none
        logical                 :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: vecLen = 3_IK
        logical     , parameter :: Mask(vecLen) = [ .false., .true., .false. ]
        real(SPR)   , parameter :: Vector2_ref(vecLen) = [(real(-i,SPR),i=1,vecLen)]
        real(SPR)   , parameter :: Vector1_ref(vecLen) = [(real(+i,SPR),i=1,vecLen)]
        real(SPR)               :: Vector1(vecLen)
        real(SPR)               :: Vector2(vecLen)
        Vector1 = Vector1_ref
        Vector2 = Vector2_ref
        call masked_swap_SPRV(Vector1,Vector2,Mask)
        assertion = all((Vector1==Vector2_ref) .eqv. Mask) .and. all( (Vector2==Vector1_ref) .eqv. Mask)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "Mask        =", Mask
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1_ref =", Vector1_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2     =", Vector2
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector2_ref =", Vector2_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector1     =", Vector1
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_masked_swap_SPRV_1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_masked_swap_SPRM_1() result(assertion)
        use Constants_mod, only: IK, SPR
        implicit none
        logical                 :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: nrow = 3_IK
        integer(IK) , parameter :: ncol = 2_IK
        logical     , parameter :: Mask(nrow,ncol) = reshape([ .false., .true., .false., .false., .true., .false. ], shape = shape(Mask))
        real(SPR)   , parameter :: Matrix2_ref(nrow,ncol) = reshape([(real(-i,SPR),i=1,nrow*ncol)], shape = shape(Mask))
        real(SPR)   , parameter :: Matrix1_ref(nrow,ncol) = reshape([(real(+i,SPR),i=1,nrow*ncol)], shape = shape(Mask))
        real(SPR)               :: Matrix1(nrow,ncol)
        real(SPR)               :: Matrix2(nrow,ncol)
        Matrix1 = Matrix1_ref
        Matrix2 = Matrix2_ref
        call masked_swap_SPRM(Matrix1,Matrix2,mask)
        assertion = all((Matrix1==Matrix2_ref) .eqv. mask) .and. all( (Matrix2==Matrix1_ref) .eqv. mask)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "Mask        =", Mask
            write(Test%outputUnit,"(*(g0,:,' '))") "Matrix1_ref =", Matrix1_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Matrix2     =", Matrix2
            write(Test%outputUnit,"(*(g0,:,' '))") "Matrix2_ref =", Matrix2_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Matrix1     =", Matrix1
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_masked_swap_SPRM_1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_arth_IK_1() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        integer(IK)                 :: i
        integer(IK) , parameter     :: first = 13, increment = 5, n = 10
        integer(IK) , parameter     :: ArithmeticProgression_ref(n) = [(i*increment+first,i=0,n-1)]
        integer(IK), allocatable    :: ArithmeticProgression(:)
        ArithmeticProgression = arth_IK(first = first,increment = increment, n = n)
        assertion = all(ArithmeticProgression == ArithmeticProgression_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "ArithmeticProgression_ref   =", ArithmeticProgression_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "ArithmeticProgression       =", ArithmeticProgression
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_arth_IK_1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_arth_RK_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                     :: assertion
        integer(IK)                 :: i
        real(RK)    , parameter     :: first = 13._RK, increment = 5._RK
        integer(IK) , parameter     :: n = 10_IK
        real(RK)    , parameter     :: ArithmeticProgression_ref(n) = [(real(i*increment+first,RK),i=0,n-1)]
        real(RK)    , allocatable   :: ArithmeticProgression(:)
        ArithmeticProgression = arth_RK(first = first, increment = increment, n = n)
        assertion = all(ArithmeticProgression == ArithmeticProgression_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "ArithmeticProgression_ref   =", ArithmeticProgression_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "ArithmeticProgression       =", ArithmeticProgression
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_arth_RK_1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_arth_RK_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                     :: assertion
        integer(IK)                 :: i
        real(RK)    , parameter     :: first = 13._RK, increment = 5._RK
        integer(IK) , parameter     :: n = 20_IK
        real(RK)    , parameter     :: ArithmeticProgression_ref(n) = [(real(i*increment+first,RK),i=0,n-1)]
        real(RK)    , allocatable   :: ArithmeticProgression(:)
        ArithmeticProgression = arth_RK(first = first, increment = increment, n = n)
        assertion = all(ArithmeticProgression == ArithmeticProgression_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "ArithmeticProgression_ref   =", ArithmeticProgression_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "ArithmeticProgression       =", ArithmeticProgression
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_arth_RK_2

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_zroots_unity_1() result(assertion)
        use Constants_mod, only: IK, RK, CK
        implicit none
        logical                     :: assertion
        complex(CK) , allocatable   :: Zroots(:)
        integer(IK) , parameter     :: n = 5, nn = 10
        real(RK)    , parameter     :: tolerance = 1.e-10_RK
        complex(CK) , parameter     :: Zroots_ref(nn) = [ (1.000000000000000_RK, 0.000000000000000_RK) &
                                                        , (0.3090169943749475_RK, 0.9510565162951535_RK) &
                                                        , (-0.8090169943749473_RK, 0.5877852522924732_RK) &
                                                        , (-0.8090169943749476_RK, -0.5877852522924730_RK) &
                                                        , (0.3090169943749472_RK, -0.9510565162951536_RK) &
                                                        , (1.000000000000000_RK, -0.2775557561562891E-15_RK) &
                                                        , (0.3090169943749478_RK, 0.9510565162951534_RK) &
                                                        , (-0.8090169943749472_RK, 0.5877852522924734_RK) &
                                                        , (-0.8090169943749477_RK, -0.5877852522924728_RK) &
                                                        , (0.3090169943749470_RK, -0.9510565162951536_RK) ]
        Zroots = zroots_unity(n,nn)
        assertion = all(abs(real(Zroots) - real(Zroots_ref)) < tolerance) .and. all(abs(aimag(Zroots) - aimag(Zroots_ref)) < tolerance)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.16,' '))")
            write(Test%outputUnit,"(*(g0.16,' '))") "Zroots_ref   =", Zroots_ref
            write(Test%outputUnit,"(*(g0.16,' '))") "Zroots       =", Zroots
            write(Test%outputUnit,"(*(g0.16,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_zroots_unity_1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_copyArray_IK_1() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        integer(IK)                 :: i
        integer(IK) , parameter     :: numCopied_ref = 5, numNotCopied_ref = 0
        integer(IK) , parameter     :: lenSource = 5, lenDestin = 10
        integer(IK) , parameter     :: Source(lenSource) = [(int(i,IK),i=1,lenSource)]
        integer(IK) , parameter     :: Destin_ref(lenDestin) = [1_IK, 2_IK, 3_IK, 4_IK, 5_IK, 0_IK, 0_IK, 0_IK, 0_IK, 0_IK]
        integer(IK)                 :: Destin(lenDestin)
        integer(IK)                 :: numCopied, numNotCopied
        Destin = 0_IK
        call copyArray_IK(Source = Source, Destination = Destin, numCopied = numCopied, numNotCopied = numNotCopied)
        assertion = all(Destin == Destin_ref) .and. (numCopied == numCopied_ref) .and. (numNotCopied == numNotCopied_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "numNotCopied_ref    =", numNotCopied_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "numNotCopied        =", numNotCopied
            write(Test%outputUnit,"(*(g0,:,' '))") "numCopied_ref       =", numCopied_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "numCopied           =", numCopied
            write(Test%outputUnit,"(*(g0,:,' '))") "Destin_ref          =", Destin_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Source              =", Source
            write(Test%outputUnit,"(*(g0,:,' '))") "Destin              =", Destin
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_copyArray_IK_1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_copyArray_IK_2() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        integer(IK)                 :: i
        integer(IK) , parameter     :: numCopied_ref = 5, numNotCopied_ref = 5
        integer(IK) , parameter     :: lenSource = 10, lenDestin = 5
        integer(IK) , parameter     :: Source(lenSource) = [(int(i,IK),i=1,lenSource)]
        integer(IK) , parameter     :: Destin_ref(lenDestin) = [1_IK, 2_IK, 3_IK, 4_IK, 5_IK]
        integer(IK)                 :: Destin(lenDestin)
        integer(IK)                 :: numCopied, numNotCopied
        Destin = 0_IK
        call copyArray_IK(Source = Source, Destination = Destin, numCopied = numCopied, numNotCopied = numNotCopied)
        assertion = all(Destin == Destin_ref) .and. (numCopied == numCopied_ref) .and. (numNotCopied == numNotCopied_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "numNotCopied_ref    =", numNotCopied_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "numNotCopied        =", numNotCopied
            write(Test%outputUnit,"(*(g0,:,' '))") "numCopied_ref       =", numCopied_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "numCopied           =", numCopied
            write(Test%outputUnit,"(*(g0,:,' '))") "Destin_ref          =", Destin_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Source              =", Source
            write(Test%outputUnit,"(*(g0,:,' '))") "Destin              =", Destin
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_copyArray_IK_2

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_copyArray_RK_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                     :: assertion
        integer(IK)                 :: i
        integer(IK) , parameter     :: numCopied_ref = 5, numNotCopied_ref = 0
        integer(IK) , parameter     :: lenSource = 5, lenDestin = 10
        real(RK)    , parameter     :: Source(lenSource) = [(real(i,RK),i=1,lenSource)]
        real(RK)    , parameter     :: Destin_ref(lenDestin) = [1._RK, 2._RK, 3._RK, 4._RK, 5._RK, 0._RK, 0._RK, 0._RK, 0._RK, 0._RK]
        real(RK)                    :: Destin(lenDestin)
        integer(IK)                 :: numCopied, numNotCopied
        Destin = 0._RK
        call copyArray_RK(Source = Source, Destination = Destin, numCopied = numCopied, numNotCopied = numNotCopied)
        assertion = all(Destin == Destin_ref) .and. (numCopied == numCopied_ref) .and. (numNotCopied == numNotCopied_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "numNotCopied_ref    =", numNotCopied_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "numNotCopied        =", numNotCopied
            write(Test%outputUnit,"(*(g0,:,' '))") "numCopied_ref       =", numCopied_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "numCopied           =", numCopied
            write(Test%outputUnit,"(*(g0,:,' '))") "Destin_ref          =", Destin_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Source              =", Source
            write(Test%outputUnit,"(*(g0,:,' '))") "Destin              =", Destin
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_copyArray_RK_1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_copyArray_RK_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                     :: assertion
        integer(IK)                 :: i
        integer(IK) , parameter     :: numCopied_ref = 5, numNotCopied_ref = 5
        integer(IK) , parameter     :: lenSource = 10, lenDestin = 5
        real(RK)    , parameter     :: Source(lenSource) = [(real(i,RK),i=1,lenSource)]
        real(RK)    , parameter     :: Destin_ref(lenDestin) = [1._RK, 2._RK, 3._RK, 4._RK, 5._RK]
        real(RK)                    :: Destin(lenDestin)
        integer(IK)                 :: numCopied, numNotCopied
        Destin = 0._RK
        call copyArray_RK(Source = Source, Destination = Destin, numCopied = numCopied, numNotCopied = numNotCopied)
        assertion = all(Destin == Destin_ref) .and. (numCopied == numCopied_ref) .and. (numNotCopied == numNotCopied_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "numNotCopied_ref    =", numNotCopied_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "numNotCopied        =", numNotCopied
            write(Test%outputUnit,"(*(g0,:,' '))") "numCopied_ref       =", numCopied_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "numCopied           =", numCopied
            write(Test%outputUnit,"(*(g0,:,' '))") "Destin_ref          =", Destin_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "Source              =", Source
            write(Test%outputUnit,"(*(g0,:,' '))") "Destin              =", Destin
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_copyArray_RK_2

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_resizeVector_RK_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                     :: assertion
        integer(IK)                 :: i
        integer(IK) , parameter     :: lenVector = 5, lenVectorNew = 10
        real(RK)    , parameter     :: Vector_ref(lenVector) = [(real(i,RK),i=1,lenVector)]
        real(RK)    , allocatable   :: Vector(:)
        Vector = Vector_ref
        call resizeVector_RK(Vector = Vector, from = lenVector, to = lenVectorNew)
        assertion = size(Vector) == lenVectorNew .and. all(Vector(1:lenVector) == Vector_ref(1:lenVector))
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector_ref      =", Vector_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "lenVectorNew    =", lenVectorNew
            write(Test%outputUnit,"(*(g0,:,' '))") "size(Vector)    =", size(Vector)
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_resizeVector_RK_1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_findUnique_1() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none
        logical                     :: assertion
        integer(IK), parameter      :: VECTOR(*) = [1,2,1,3,5,5,2]
        integer(IK), parameter      :: LEN_VECTOR = size(VECTOR)
        integer(IK), parameter      :: UNIQUE_VALUE(*) = [1,2,3,5]
        integer(IK), parameter      :: UNIQUE_COUNT(*) = [2,2,1,2]
        integer(IK), allocatable    :: UniqueValue(:), UniqueCount(:), ZeroLenVector(:)
        integer(IK)                 :: lenUnique

        call findUnique ( lenVector = LEN_VECTOR &
                        , Vector = VECTOR &
                        , UniqueValue = UniqueValue &
                        , UniqueCount = UniqueCount &
                        , lenUnique = lenUnique &
                        )

        assertion = all(UniqueValue==UNIQUE_VALUE) .and. all(UniqueCount==UNIQUE_COUNT)

        ! test with empty input vector

        allocate(ZeroLenVector(0))
        call findUnique ( lenVector = 0_IK & ! LCOV_EXCL_LINE
                        , Vector = ZeroLenVector & ! LCOV_EXCL_LINE
                        , UniqueValue = UniqueValue & ! LCOV_EXCL_LINE
                        , UniqueCount = UniqueCount & ! LCOV_EXCL_LINE
                        , lenUnique = lenUnique & ! LCOV_EXCL_LINE
                        )

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "VECTOR", VECTOR
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "UNIQUE_VALUE", UNIQUE_VALUE
            write(Test%outputUnit,"(*(g0,:,', '))") "UniqueValue ", UniqueValue
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "UNIQUE_COUNT", UNIQUE_COUNT
            write(Test%outputUnit,"(*(g0,:,', '))") "UniqueCount ", UniqueCount
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "lenUnique", lenUnique
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "VECTOR", ZeroLenVector
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "UniqueValue ", UniqueValue
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "UniqueCount ", UniqueCount
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "lenUnique", lenUnique
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_findUnique_1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_Misc_mod ! LCOV_EXCL_LINE