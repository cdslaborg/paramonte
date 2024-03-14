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

!>  \brief This module contains tests of the module [pm_mathLogSumExp](@ref pm_mathLogSumExp).
!>  \author Amir Shahmoradi

module test_pm_mathLogSumExp

    use pm_mathLogSumExp
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    use pm_kind, only: LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED
        module function test_getLogSumExp_CK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if CK2_ENABLED
        module function test_getLogSumExp_CK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if CK1_ENABLED
        module function test_getLogSumExp_CK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED
        module function test_getLogSumExp_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
        module function test_getLogSumExp_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
        module function test_getLogSumExp_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED
        call test%run(test_getLogSumExp_CK3_1, SK_"test_getLogSumExp_CK3_1")
#endif
#if CK2_ENABLED
        call test%run(test_getLogSumExp_CK2_1, SK_"test_getLogSumExp_CK2_1")
#endif
#if CK1_ENABLED
        call test%run(test_getLogSumExp_CK1_1, SK_"test_getLogSumExp_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED
        call test%run(test_getLogSumExp_RK3_1, SK_"test_getLogSumExp_RK3_1")
#endif
#if RK2_ENABLED
        call test%run(test_getLogSumExp_RK2_1, SK_"test_getLogSumExp_RK2_1")
#endif
#if RK1_ENABLED
        call test%run(test_getLogSumExp_RK1_1, SK_"test_getLogSumExp_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogSumExp_RK_1() result(assertion)
        use pm_kind, only: RK, IK
        implicit none
        logical(LK)             :: assertion
        real(RK), parameter     :: LogValue(*) = [ log(0.5*huge(1._RK)), log(0.9*huge(1._RK)), log(0.1*huge(1._RK)) ]
        real(RK), parameter     :: tolerance = 1.e-10_RK
        real(RK), parameter     :: logSumExp_ref = 710.1881779865910_RK
        real(RK)                :: logSumExp
        real(RK)                :: difference
        logSumExp = getLogSumExp(LogValue, maxval(LogValue))
        difference = abs(logSumExp - logSumExp_ref)
        assertion  = difference < tolerance
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") "logSumExp_ref   = ", logSumExp_ref
            write(test%disp%unit,"(*(g0,:,' '))") "logSumExp       = ", logSumExp
            write(test%disp%unit,"(*(g0,:,' '))") "difference      = ", difference
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogSumExp_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogSumExp_RK_2() result(assertion)
        use pm_kind, only: RK, IK
        implicit none
        logical(LK)             :: assertion
        real(RK), parameter     :: LogValue(*) = [ log(0.5*huge(1._RK)), log(0.9*huge(1._RK)), log(0.1*huge(1._RK)) ]
        real(RK), parameter     :: tolerance = 1.e-10_RK
        real(RK), parameter     :: logSumExp_ref = 710.1881779865910_RK
        real(RK)                :: logSumExp
        real(RK)                :: difference
        logSumExp = getLogSumExp(LogValue, maxval(LogValue))
        difference = abs(logSumExp - logSumExp_ref)
        assertion  = difference < tolerance
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") "logSumExp_ref   = ", logSumExp_ref
            write(test%disp%unit,"(*(g0,:,' '))") "logSumExp       = ", logSumExp
            write(test%disp%unit,"(*(g0,:,' '))") "difference      = ", difference
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogSumExp_RK_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogSumExp_CK_1() result(assertion)
        use pm_kind, only: RK, IK
        implicit none
        logical(LK)             :: assertion
        complex(RK), parameter  :: LogValue(*) =    [ cmplx( log(0.5*huge(1._RK)), 0._RK, kind = RK ) &
                                                    , cmplx( log(0.9*huge(1._RK)), 0._RK, kind = RK ) &
                                                    , cmplx( log(0.1*huge(1._RK)), 0._RK, kind = RK ) &
                                                    ]
        real(RK), parameter     :: logSumExp_ref = cmplx(710.1881779865910_RK, 0._RK, RK)
        complex(RK), parameter  :: tolerance = cmplx(1.e-10_RK, 0._RK, RK)
        complex(RK)             :: logSumExp
        complex(RK)             :: difference
        logSumExp = getLogSumExp(LogValue, cmplx( maxval(real(LogValue,kind=RK)), kind = RK ))
        difference = abs(logSumExp - logSumExp_ref)
        assertion  = real(difference,RK) < real(tolerance,RK)
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") "logSumExp_ref   = ", logSumExp_ref
            write(test%disp%unit,"(*(g0,:,' '))") "logSumExp       = ", logSumExp
            write(test%disp%unit,"(*(g0,:,' '))") "difference      = ", difference
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogSumExp_CK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogSumExp_CK_2() result(assertion)
        use pm_kind, only: RK, IK
        implicit none
        logical(LK)             :: assertion
        complex(RK), parameter  :: LogValue(*) =    [ cmplx( log(0.5*huge(1._RK)), 0._RK, kind = RK ) &
                                                    , cmplx( log(0.9*huge(1._RK)), 0._RK, kind = RK ) &
                                                    , cmplx( log(0.1*huge(1._RK)), 0._RK, kind = RK ) &
                                                    ]
        real(RK), parameter     :: logSumExp_ref = cmplx(710.1881779865910_RK, 0._RK, RK)
        complex(RK), parameter  :: tolerance = cmplx(1.e-10_RK, 0._RK, RK)
        complex(RK)             :: logSumExp
        complex(RK)             :: difference
        logSumExp = getLogSumExp(LogValue, cmplx( maxval(real(LogValue,kind=RK)), kind = RK ))
        difference = abs(logSumExp - logSumExp_ref)
        assertion  = real(difference,RK) < real(tolerance,RK)
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") "logSumExp_ref   = ", logSumExp_ref
            write(test%disp%unit,"(*(g0,:,' '))") "logSumExp       = ", logSumExp
            write(test%disp%unit,"(*(g0,:,' '))") "difference      = ", difference
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogSumExp_CK_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_mathLogSumExp ! LCOV_EXCL_LINE