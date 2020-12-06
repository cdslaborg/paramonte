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

!>  \brief This module contains tests of the module [Parallelism_mod](@ref parallelism_mod).
!>  \author Amir Shahmoradi

module Test_Parallelism_mod

    use Parallelism_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_Parallelism

    type(Test_type) :: Test

    type :: TestFuncRosenBrock1D_type
        ! MATLAB reference code:
        ! function f = objectivefcn1(x)
        !     f = exp(-x^2) * cos(x) *sin(2*x);
        ! end
        ! CMD: x0 = [0.25];
        ! CMD: options = optimset('tolx',1.e-10)
        ! CMD: [x, f] = fminsearch(@objectivefcn1,x0,options)
        real(RK) :: XMIN_REF = -0.471354350447655_RK
        real(RK) :: FMIN_REF = -0.577293243101421_RK
    contains
        procedure, nopass :: get => getTestFuncRosenBrock1D
    end type TestFuncRosenBrock1D_type

    type :: TestFuncRosenBrock2D_type
        ! MATLAB reference code:
        ! function f = objectivefcn1(x)
        !     f = exp(-(x(1)-x(2))^2 - 2*x(1)^2)*cos(x(2))*sin(2*x(2));
        ! end
        ! CMD: x0 = [0.25,-0.25];
        ! CMD: options = optimset('tolx',1.e-10)
        ! CMD: [x, f] = fminsearch(@objectivefcn1,x0,options)
        real(RK) :: XMIN_REF(2) = [-0.169552635123540_RK, -0.508657910611664_RK]
        real(RK) :: FMIN_REF = -0.625285429865811_RK
    contains
        procedure, nopass :: get => getTestFuncRosenBrock2D
    end type TestFuncRosenBrock2D_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_Parallelism()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_BrentMinimum_type_1, "test_BrentMinimum_type_1")
        call Test%run(test_BrentMinimum_type_2, "test_BrentMinimum_type_2")
        call Test%run(test_PowellMinimum_type_1, "test_PowellMinimum_type_1") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call Test%finalize()
    end subroutine test_Parallelism

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_BrentMinimum_type_1() result(assertion)

        use Constants_mod, only: RK, IK
        implicit none

        logical                         :: assertion
        type(BrentMinimum_type)         :: BrentMinimum
        type(TestFuncRosenBrock1D_type) :: TestFuncRosenBrock1D

        assertion = .true.

        BrentMinimum = minimize( getFunc = getTestFuncRosenBrock1D, xtol = 1.e-8_RK )

        assertion = .not. BrentMinimum%Err%occurred
        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") BrentMinimum%Err%msg
                write(Test%outputUnit,"(*(g0,:,', '))")
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        assertion = assertion .and. abs(BrentMinimum%xmin-TestFuncRosenBrock1D%XMIN_REF) / abs(TestFuncRosenBrock1D%XMIN_REF) < 1.e-6_RK
        assertion = assertion .and. abs(BrentMinimum%fmin-TestFuncRosenBrock1D%FMIN_REF) / abs(TestFuncRosenBrock1D%FMIN_REF) < 1.e-6_RK

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "BrentMinimum%xmin", BrentMinimum%xmin
            write(Test%outputUnit,"(*(g0,:,', '))") "XMIN_REF         ", TestFuncRosenBrock1D%XMIN_REF
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "BrentMinimum%fmin", BrentMinimum%fmin
            write(Test%outputUnit,"(*(g0,:,', '))") "FMIN_REF         ", TestFuncRosenBrock1D%FMIN_REF
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "BrentMinimum%niter", BrentMinimum%niter
            write(Test%outputUnit,"(*(g0,:,', '))") "BrentMinimum%Bracket", BrentMinimum%Bracket
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_BrentMinimum_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_BrentMinimum_type_2() result(assertion)

        use Constants_mod, only: RK, IK
        implicit none

        logical                         :: assertion
        type(BrentMinimum_type)         :: BrentMinimum
        type(TestFuncRosenBrock1D_type) :: TestFuncRosenBrock1D

        assertion = .true.

        BrentMinimum = BrentMinimum_type( getFunc = getTestFuncRosenBrock1D &
                                        , xtol = 1.e-8_RK &
                                        , x0 = -1._RK &
                                        , x1 = 0._RK &
                                        , x2 = .5_RK &
                                        )

        assertion = .not. BrentMinimum%Err%occurred
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") BrentMinimum%Err%msg
                write(Test%outputUnit,"(*(g0,:,', '))")
            end if
            return
            ! LCOV_EXCL_STOP
        end if

        assertion = abs(BrentMinimum%xmin-TestFuncRosenBrock1D%XMIN_REF) / abs(TestFuncRosenBrock1D%XMIN_REF) < 1.e-6_RK
        assertion = abs(BrentMinimum%fmin-TestFuncRosenBrock1D%FMIN_REF) / abs(TestFuncRosenBrock1D%FMIN_REF) < 1.e-6_RK

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "BrentMinimum%xmin", BrentMinimum%xmin
            write(Test%outputUnit,"(*(g0,:,', '))") "XMIN_REF         ", TestFuncRosenBrock1D%XMIN_REF
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "BrentMinimum%fmin", BrentMinimum%fmin
            write(Test%outputUnit,"(*(g0,:,', '))") "FMIN_REF         ", TestFuncRosenBrock1D%FMIN_REF
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "BrentMinimum%niter", BrentMinimum%niter
            write(Test%outputUnit,"(*(g0,:,', '))") "BrentMinimum%Bracket", BrentMinimum%Bracket
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_BrentMinimum_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_PowellMinimum_type_1() result(assertion)

        use Constants_mod, only: RK, IK
        implicit none

        logical                         :: assertion
        type(PowellMinimum_type)        :: PowellMinimum
        type(TestFuncRosenBrock2D_type) :: TestFuncRosenBrock2D
        real(RK), allocatable           :: StartVec(:)

        assertion = .true.

        StartVec = [.25_RK,-.25_RK]

        PowellMinimum = PowellMinimum_type  ( ndim = 2_IK &
                                            , getFuncMD = getTestFuncRosenBrock2D &
                                            , StartVec = StartVec & ! LCOV_EXCL_LINE
                                            !, DirMat = reshape([1._RK, 0._RK, 0._RK, 1._RK], shape = [2,2])
                                            !, ftol = 1.e-8_RK &
                                            )

        assertion = .not. PowellMinimum%Err%occurred
        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") PowellMinimum%Err%msg
                write(Test%outputUnit,"(*(g0,:,', '))")
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        assertion = assertion .and. any(abs(PowellMinimum%xmin-TestFuncRosenBrock2D%XMIN_REF) / abs(TestFuncRosenBrock2D%XMIN_REF) < 1.e-6_RK)
        assertion = assertion .and. abs(PowellMinimum%fmin-TestFuncRosenBrock2D%FMIN_REF) / abs(TestFuncRosenBrock2D%FMIN_REF) < 1.e-6_RK

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "PowellMinimum%xmin", PowellMinimum%xmin
            write(Test%outputUnit,"(*(g0,:,', '))") "XMIN_REF          ", TestFuncRosenBrock2D%XMIN_REF
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "PowellMinimum%fmin", PowellMinimum%fmin
            write(Test%outputUnit,"(*(g0,:,', '))") "FMIN_REF          ", TestFuncRosenBrock2D%FMIN_REF
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "PowellMinimum%niter", PowellMinimum%niter
            write(Test%outputUnit,"(*(g0,:,', '))") "PowellMinimum%ftol", PowellMinimum%ftol
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "PowellMinimum%DirMat", PowellMinimum%DirMat
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_PowellMinimum_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getTestFuncRosenBrock1D(x) result(testFuncVal)
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: x
        real(RK)                :: testFuncVal
        testFuncVal = exp(-x**2) * cos(x) * sin(2*x)
    end function getTestFuncRosenBrock1D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getTestFuncRosenBrock2D(ndim,Point) result(testFuncVal)
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK) , intent(in)    :: ndim
        real(RK)    , intent(in)    :: Point(ndim)
        real(RK)                    :: testFuncVal
        testFuncVal = exp(-(Point(1)-Point(2))**2 - 2*Point(1)**2) * cos(Point(2)) * sin(2*Point(2))
    end function getTestFuncRosenBrock2D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_Parallelism_mod ! LCOV_EXCL_LINE