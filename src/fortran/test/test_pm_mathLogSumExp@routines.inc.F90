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

!>  \brief
!>  This include file contains procedure implementations of the tests of [getLogSumExp](@ref pm_mathLogSumExp::getLogSumExp).
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 2:06 AM, September 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) , parameter     :: NP = 10_IK
#if test_getLogSumExp_CK_ENABLED
        complex(CK) , parameter     :: TOL = cmplx(epsilon(1._CK),epsilon(1._CK),CK) * 10_IK
        complex(CK)                 :: logSumExp_ref
        complex(CK)                 :: logSumExp
        complex(CK)                 :: array(NP)
        complex(CK)                 :: maxArray
        complex(CK)                 :: diff(NP)
        real(CK)                    :: Temp(NP)
#elif test_getLogSumExp_RK_ENABLED
        real(RK)    , parameter     :: TOL = epsilon(1._RK) * 10_IK
        real(RK)                    :: logSumExp_ref
        real(RK)                    :: logSumExp
        real(RK)                    :: array(NP)
        real(RK)                    :: maxArray
        real(RK)                    :: diff(NP)
#else
#error "Unrecognized interface."
#endif

        assertion = .true._LK

        call runTest()
        call runTest(selection)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTest(control)

            type(selection_type), intent(in), optional :: control

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Compute the logSumExp for non-overflowing numbers
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if test_getLogSumExp_CK_ENABLED
            call random_number(Temp); array%re = Temp
            call random_number(Temp); array%im = Temp
            maxArray%re = maxval(real(array,kind=CK))
            maxArray%im = maxval(aimag(array))
#elif test_getLogSumExp_RK_ENABLED
            call random_number(array)
            maxArray = maxval(array)
#endif

            logSumExp_ref = maxArray + log(sum(exp(array - maxArray)))
            logSumExp = getLogSumExp(array, maxArray, cenabled)
            diff = abs(logSumExp - logSumExp_ref)
#if test_getLogSumExp_CK_ENABLED
            assertion = assertion .and. all(diff%re < TOL%re) .and. all(diff%im < TOL%im)
#elif test_getLogSumExp_RK_ENABLED
            assertion = assertion .and. all(diff < TOL)
#endif
            call report()
            call test%assert(assertion, desc = "The logSumExp must be computed correctly for two non-overflowing numbers.")

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Compute the logSumExp for non-overflowing numbers
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if test_getLogSumExp_CK_ENABLED
            array = huge(0._CK)**0.01_CK * array
            maxArray%re = maxval(real(array,kind=CK))
            maxArray%im = maxval(aimag(array))
#elif test_getLogSumExp_RK_ENABLED
            array = huge(0._RK)**0.01_RK * array
            maxArray = maxval(array)
#endif
            logSumExp_ref = maxArray + log(sum(exp(array - maxArray)))
            logSumExp = getLogSumExp(array, maxArray, cenabled)
            diff = abs(logSumExp - logSumExp_ref)
#if test_getLogSumExp_CK_ENABLED
            assertion = assertion .and. all(diff%re < TOL%re) .and. all(diff%im < TOL%im)
#elif test_getLogSumExp_RK_ENABLED
            assertion = assertion .and. all(diff < TOL)
#endif
            call report()
            call test%assert(assertion, desc = "The logSumExp must be computed correctly for two overflowing numbers.")
        
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "maxArray       ", maxArray
                write(test%disp%unit,"(*(g0,:,', '))") "logSumExp_ref  ", logSumExp_ref
                write(test%disp%unit,"(*(g0,:,', '))") "logSumExp      ", logSumExp
                write(test%disp%unit,"(*(g0,:,', '))") "diff           ", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL            ", TOL
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
