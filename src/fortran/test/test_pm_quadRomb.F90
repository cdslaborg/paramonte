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

!>  \brief This module contains tests of the module [pm_quadRomb](@ref pm_quadRomb).
!>  \author Amir Shahmoradi

module test_pm_quadRomb

    use pm_quadRomb
    !use QuadPackDouble_pmod, only: dqagi
    use pm_except, only: POSBIG_RK
    use pm_kind, only: IK, LK, RK
    use pm_err, only: err_type
    use pm_test, only: test_type, LK

    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)
        call test%run(test_doQuadRombOpen_1, SK_"test_doQuadRombOpen_1")
        call test%run(test_doQuadRombOpen_2, SK_"test_doQuadRombOpen_2")
        call test%run(test_doQuadRombOpen_3, SK_"test_doQuadRombOpen_3")
        call test%run(test_doQuadRombOpen_4, SK_"test_doQuadRombOpen_4")
        call test%run(test_getQuadRomb_1, SK_"test_getQuadRomb_1")
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! integrate the Gaussian PDF via the midexp routine on the open interval \f$[0, +\infty)\f$.
    function test_doQuadRombOpen_1() result(assertion)

        use pm_quadRomb

        implicit none

        logical(LK)             :: assertion
        real(RK)                :: integral, relativeError, difference 
        real(RK), parameter     :: integral_ref = 0.5_RK
        real(RK), parameter     :: tolerance = 1.e-10_RK
        integer(IK)             :: numFuncEval, ierr

        call doQuadRombOpen ( getFunc = getTestFuncOpenInterval_1 &
                            , integrate = midexp &
                            , lowerLim = 0._RK &
                            , upperLim = POSBIG_RK & ! do not set this to huge() as GNU Fortran in debug mode crashes
                            , maxRelativeError = 0.1*tolerance &
                            , nRefinement = 10_IK &
                            , integral = integral &
                            , relativeError = relativeError &
                            , numFuncEval = numFuncEval &
                            , ierr = ierr &
                            )

        assertion = ierr == 0_IK

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable .and. .not. assertion) then
                write(test%disp%unit,"(*(g0))") "ierr = ", ierr, " /= 0"
            end if
            return
            ! LCOV_EXCL_STOP
        end if

        difference = abs(integral - integral_ref) / abs(integral_ref)
        assertion = difference < tolerance
        assertion = assertion .and. relativeError <= tolerance

        ! LCOV_EXCL_START
        if (test%traceable .and. .not. assertion) then
            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))") "integral_ref  = ", integral_ref
            write(test%disp%unit,"(*(g0))") "integral      = ", integral
            write(test%disp%unit,"(*(g0))") "difference    = ", difference
            write(test%disp%unit,"(*(g0))") "relativeError = ", relativeError
            write(test%disp%unit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_doQuadRombOpen_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_doQuadRombOpen_2() result(assertion)

        use pm_quadRomb

        implicit none

        logical(LK)             :: assertion
        real(RK)                :: relativeError, difference
        real(RK)                :: integral
        real(RK), parameter     :: integral_ref = 0.078649603525143_RK
        real(RK), parameter     :: tolerance = 1.e-10_RK
        integer(IK)             :: numFuncEval, ierr

        call doQuadRombOpen ( getFunc = getTestFuncOpenInterval_1 &
                            , integrate = midinf &
                            , lowerLim = 1._RK &
                            , upperLim = POSBIG_RK & ! do not set this to huge() as GNU Fortran in debug mode crashes
                            , maxRelativeError = 0.1*tolerance &
                            , nRefinement = 10_IK &
                            , integral = integral &
                            , relativeError = relativeError &
                            , numFuncEval = numFuncEval &
                            , ierr = ierr &
                            )

        assertion = ierr == 0_IK

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable .and. .not. assertion) then
                write(test%disp%unit,"(*(g0))") "ierr = ", ierr, " /= 0"
            end if
            return
            ! LCOV_EXCL_STOP
        end if

        difference = abs(integral - integral_ref) / abs(integral_ref)
        assertion = difference < tolerance
        assertion = assertion .and. relativeError <= tolerance

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))") "integral_ref  = ", integral_ref
            write(test%disp%unit,"(*(g0))") "integral      = ", integral
            write(test%disp%unit,"(*(g0))") "difference    = ", difference
            write(test%disp%unit,"(*(g0))") "relativeError = ", relativeError
            write(test%disp%unit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_doQuadRombOpen_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_doQuadRombOpen_3() result(assertion)

        use pm_quadRomb

        implicit none

        logical(LK)             :: assertion
        real(RK)                :: relativeError, difference
        real(RK)                :: integral
        real(RK), parameter     :: integral_ref = 0.078649603525143_RK
        real(RK), parameter     :: tolerance = 1.e-4_RK
        integer(IK)             :: numFuncEval, ierr

        call doQuadRombOpen ( getFunc = getTestFuncOpenInterval_1 &
                            , integrate = midinf &
                            , lowerLim = POSBIG_RK & ! do not set this to huge() as GNU Fortran in debug mode crashes
                            , upperLim = -1._RK &
                            , maxRelativeError = 1.e-10_RK &
                            , nRefinement = 10_IK &
                            , integral = integral &
                            , relativeError = relativeError &
                            , numFuncEval = numFuncEval &
                            , ierr = ierr &
                            )

        assertion = ierr == 0_IK
        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable .and. .not. assertion) then
                write(test%disp%unit,"(*(g0))") "ierr = ", ierr, " /= 0"
            end if
            return
            ! LCOV_EXCL_STOP
        end if

        difference = abs(integral - integral_ref) / abs(integral_ref)
        assertion = difference < tolerance
        assertion = assertion .and. relativeError <= tolerance

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))") "integral_ref  = ", integral_ref
            write(test%disp%unit,"(*(g0))") "integral      = ", integral
            write(test%disp%unit,"(*(g0))") "difference    = ", difference
            write(test%disp%unit,"(*(g0))") "relativeError = ", relativeError
            write(test%disp%unit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_doQuadRombOpen_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_doQuadRombOpen_4() result(assertion)

        use pm_quadRomb

        implicit none

        logical(LK)             :: assertion
        real(RK)                :: relativeError, difference
        real(RK)                :: integral
        real(RK), parameter     :: integral_ref = 1._RK
        real(RK), parameter     :: tolerance = 1.e-10_RK
        integer(IK)             :: numFuncEval, ierr

        call doQuadRombOpen ( getFunc = getTestFuncOpenInterval_1 &
                            , integrate = midpnt &
                            , lowerLim = -5._RK &
                            , upperLim = +5._RK &
                            , maxRelativeError = 1.e-10_RK &
                            , nRefinement = 10_IK &
                            , integral = integral &
                            , relativeError = relativeError &
                            , numFuncEval = numFuncEval &
                            , ierr = ierr &
                            )

        assertion = ierr == 0_IK
        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable .and. .not. assertion) then
                write(test%disp%unit,"(*(g0))") "ierr = ", ierr, " /= 0"
            end if
            return
            ! LCOV_EXCL_STOP
        end if

        difference = abs(integral - integral_ref) / abs(integral_ref)
        assertion = difference < tolerance
        assertion = assertion .and. relativeError <= tolerance

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))") "integral_ref  = ", integral_ref
            write(test%disp%unit,"(*(g0))") "integral      = ", integral
            write(test%disp%unit,"(*(g0))") "difference    = ", difference
            write(test%disp%unit,"(*(g0))") "relativeError = ", relativeError
            write(test%disp%unit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_doQuadRombOpen_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getQuadRomb_1() result(assertion)

        use pm_quadRomb

        implicit none

        logical(LK)             :: assertion
        real(RK)                :: relativeError, difference
        real(RK)                :: integral
        real(RK), parameter     :: integral_ref = 1._RK
        real(RK), parameter     :: tolerance = 1.e-10_RK
        integer(IK)             :: numFuncEval, ierr

        integral = getQuadRomb  ( getFunc = getTestFuncOpenInterval_1 &
                                , lowerLim = -5._RK &
                                , upperLim = +5._RK &
                                , maxRelativeError = 0.1*tolerance &
                                , nRefinement = 10_IK &
                                , relativeError = relativeError &
                                , numFuncEval = numFuncEval &
                                , ierr = ierr &
                                )
        assertion = ierr == 0_IK
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            if (test%traceable .and. .not. assertion) then
                write(test%disp%unit,"(*(g0))") "ierr = ", ierr, " /= 0"
            end if
            return
        ! LCOV_EXCL_STOP
        end if

        difference = abs(integral - integral_ref) / abs(integral_ref)
        assertion = difference < tolerance
        assertion = assertion .and. relativeError <= tolerance

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))") "integral_ref  = ", integral_ref
            write(test%disp%unit,"(*(g0))") "integral      = ", integral
            write(test%disp%unit,"(*(g0))") "difference    = ", difference
            write(test%disp%unit,"(*(g0))") "relativeError = ", relativeError
            write(test%disp%unit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getQuadRomb_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! The Gaussian PDF with sigma = 1/sqrt(2)
    function getTestFuncOpenInterval_1(x) result(funcVal)
        real(RK)    , parameter     :: LOG_INVERSE_SQRT_PI = log(sqrt(1._RK / acos(-1._RK)))
        real(RK)    , parameter     :: LOG_TINY = log(tiny(0._RK))
        real(RK)    , intent(in)    :: x
        real(RK)                    :: funcVal
        funcVal = LOG_INVERSE_SQRT_PI - x**2
        if (funcVal < LOG_TINY) then ! This takes care of the GNU Fortran 9.1 test crash in debug mode.
            funcVal = 0._RK
        else
            funcVal = exp(funcVal)
        end if
    end function getTestFuncOpenInterval_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! requires `midsql`.
    ! The MATLAB test function: https://www.mathworks.com/help/matlab/ref/integral.html
    ! LCOV_EXCL_START
    function getTestFuncOpenIntervalMATLAB(x) result(funcVal)
        real(RK), intent(in)    :: x
        real(RK)                :: funcVal
        funcVal = exp(-x**2) * log(x)**2
    end function getTestFuncOpenIntervalMATLAB
    ! LCOV_EXCL_STOP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_quadRomb ! LCOV_EXCL_LINE