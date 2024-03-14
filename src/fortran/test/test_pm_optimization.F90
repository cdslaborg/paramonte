!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This module contains tests of the module [pm_optimization](@ref pm_optimization).
!>  \author Amir Shahmoradi

module test_pm_optimization

    use pm_optimization
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

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

    subroutine setTest()

        test = test_type(MODULE_NAME)
        call test%run(test_BrentMinimum_type_1, SK_"test_BrentMinimum_type_1")
        call test%run(test_BrentMinimum_type_2, SK_"test_BrentMinimum_type_2")
        call test%run(test_powellMin_type_1, SK_"test_powellMin_type_1") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_BrentMinimum_type_1() result(assertion)

        use pm_kind, only: RK, IK
        implicit none

        logical(LK)                     :: assertion
        type(brentMinimum_type)         :: BrentMinimum
        type(TestFuncRosenBrock1D_type) :: TestFuncRosenBrock1D

        assertion = .true._LK

        BrentMinimum = minimize( getFunc = getTestFuncRosenBrock1D, xtol = 1.e-8_RK )

        assertion = .not. BrentMinimum%err%occurred
        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (test%traceable) then
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") BrentMinimum%err%msg
                write(test%disp%unit,"(*(g0,:,', '))")
            end if
        end if
        ! LCOV_EXCL_STOP
        call test%assert(assertion)

        assertion = assertion .and. abs(BrentMinimum%xmin-TestFuncRosenBrock1D%XMIN_REF) / abs(TestFuncRosenBrock1D%XMIN_REF) < 1.e-6_RK
        call test%assert(assertion)
        assertion = assertion .and. abs(BrentMinimum%fmin-TestFuncRosenBrock1D%FMIN_REF) / abs(TestFuncRosenBrock1D%FMIN_REF) < 1.e-6_RK
        call test%assert(assertion)

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "BrentMinimum%xmin", BrentMinimum%xmin
            write(test%disp%unit,"(*(g0,:,', '))") "XMIN_REF         ", TestFuncRosenBrock1D%XMIN_REF
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "BrentMinimum%fmin", BrentMinimum%fmin
            write(test%disp%unit,"(*(g0,:,', '))") "FMIN_REF         ", TestFuncRosenBrock1D%FMIN_REF
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "BrentMinimum%niter", BrentMinimum%niter
            write(test%disp%unit,"(*(g0,:,', '))") "BrentMinimum%Bracket", BrentMinimum%Bracket
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_BrentMinimum_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_BrentMinimum_type_2() result(assertion)

        use pm_kind, only: RK, IK
        implicit none

        logical(LK)                     :: assertion
        type(brentMinimum_type)         :: BrentMinimum
        type(TestFuncRosenBrock1D_type) :: TestFuncRosenBrock1D

        assertion = .true._LK

        BrentMinimum = brentMinimum_type( getFunc = getTestFuncRosenBrock1D &
                                        , xtol = 1.e-8_RK &
                                        , x0 = -1._RK &
                                        , x1 = 0._RK &
                                        , x2 = .5_RK &
                                        )

        assertion = .not. BrentMinimum%err%occurred
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable) then
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") BrentMinimum%err%msg
                write(test%disp%unit,"(*(g0,:,', '))")
            end if
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

        assertion = abs(BrentMinimum%xmin-TestFuncRosenBrock1D%XMIN_REF) / abs(TestFuncRosenBrock1D%XMIN_REF) < 1.e-6_RK
        call test%assert(assertion)
        assertion = abs(BrentMinimum%fmin-TestFuncRosenBrock1D%FMIN_REF) / abs(TestFuncRosenBrock1D%FMIN_REF) < 1.e-6_RK
        call test%assert(assertion)

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "BrentMinimum%xmin", BrentMinimum%xmin
            write(test%disp%unit,"(*(g0,:,', '))") "XMIN_REF         ", TestFuncRosenBrock1D%XMIN_REF
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "BrentMinimum%fmin", BrentMinimum%fmin
            write(test%disp%unit,"(*(g0,:,', '))") "FMIN_REF         ", TestFuncRosenBrock1D%FMIN_REF
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "BrentMinimum%niter", BrentMinimum%niter
            write(test%disp%unit,"(*(g0,:,', '))") "BrentMinimum%Bracket", BrentMinimum%Bracket
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_BrentMinimum_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_powellMin_type_1() result(assertion)

        use pm_kind, only: RK, IK
        implicit none

        logical(LK)                     :: assertion
        type(powellMin_type)        :: powellMin
        type(TestFuncRosenBrock2D_type) :: TestFuncRosenBrock2D
        real(RK), allocatable           :: StartVec(:)

        assertion = .true._LK

        StartVec = [.25_RK,-.25_RK]

        powellMin = powellMin_type  ( ndim = 2_IK &
                                            , getFuncMD = getTestFuncRosenBrock2D &
                                            , StartVec = StartVec & ! LCOV_EXCL_LINE
                                            !, DirMat = reshape([1._RK, 0._RK, 0._RK, 1._RK], shape = [2,2])
                                            !, ftol = 1.e-8_RK &
                                            )

        assertion = .not. powellMin%err%occurred
        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (test%traceable) then
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") powellMin%err%msg
                write(test%disp%unit,"(*(g0,:,', '))")
            end if
        end if
        ! LCOV_EXCL_STOP
        call test%assert(assertion)

        assertion = assertion .and. any(abs(powellMin%xmin-TestFuncRosenBrock2D%XMIN_REF) / abs(TestFuncRosenBrock2D%XMIN_REF) < 1.e-6_RK)
        call test%assert(assertion)
        assertion = assertion .and. abs(powellMin%fmin-TestFuncRosenBrock2D%FMIN_REF) / abs(TestFuncRosenBrock2D%FMIN_REF) < 1.e-6_RK
        call test%assert(assertion)

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "powellMin%xmin", powellMin%xmin
            write(test%disp%unit,"(*(g0,:,', '))") "XMIN_REF          ", TestFuncRosenBrock2D%XMIN_REF
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "powellMin%fmin", powellMin%fmin
            write(test%disp%unit,"(*(g0,:,', '))") "FMIN_REF          ", TestFuncRosenBrock2D%FMIN_REF
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "powellMin%niter", powellMin%niter
            write(test%disp%unit,"(*(g0,:,', '))") "powellMin%ftol", powellMin%ftol
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "powellMin%DirMat", powellMin%DirMat
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_powellMin_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getTestFuncRosenBrock1D(x) result(testFuncVal)
        use pm_kind, only: RK
        implicit none
        real(RK), intent(in)    :: x
        real(RK)                :: testFuncVal
        testFuncVal = exp(-x**2) * cos(x) * sin(2*x)
    end function getTestFuncRosenBrock1D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getTestFuncRosenBrock2D(ndim,Point) result(testFuncVal)
        use pm_kind, only: RK, IK
        implicit none
        integer(IK) , intent(in)    :: ndim
        real(RK)    , intent(in)    :: Point(ndim)
        real(RK)                    :: testFuncVal
        testFuncVal = exp(-(Point(1)-Point(2))**2 - 2*Point(1)**2) * cos(Point(2)) * sin(2*Point(2))
    end function getTestFuncRosenBrock2D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_optimization ! LCOV_EXCL_LINE