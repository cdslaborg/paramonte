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
!>  This include file contains procedure implementations of the tests of [pm_mathLogSubExp](@ref pm_mathLogSubExp).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_mathMinMax, only: setMinMax

#if     test_getLogSubExp_CK_ENABLED
        complex(CKC), parameter     :: TOL = epsilon(1._CKC) * 10_IK
        complex(CKC)                :: logSubExp_ref
        complex(CKC)                :: logSubExp
        complex(CKC)                :: smaller
        complex(CKC)                :: larger
        complex(CKC)                :: diff
        real(CKC)                   :: tmp
#elif   test_getLogSubExp_RK_ENABLED
        real(RKC)   , parameter     :: TOL = epsilon(1._RKC) * 10_IK
        real(RKC)                   :: logSubExp_ref
        real(RKC)                   :: logSubExp
        real(RKC)                   :: smaller
        real(RKC)                   :: larger
        real(RKC)                   :: diff
#else
#error  "Unrecognized interface."
#endif

        assertion = .true._LK

        call runTest(cenabled = .true._LK)
        call runTest(cenabled = .false._LK)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTest(cenabled)

            logical(LK), intent(in) :: cenabled

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Compute the logSubExp for non-overflowing numbers
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         test_getLogSubExp_CK_ENABLED
            call setUnifRand(larger)
            call setUnifRand(smaller)
            call setMinMax(smaller, larger)
            logSubExp_ref = larger + log(1._CK + exp(smaller - larger))
            logSubExp = getLogSubExp(smaller = smaller, larger = larger, cenabled = cenabled)
            diff = logSubExp - logSubExp_ref
            assertion = diff%re < TOL%re
#elif       test_getLogSubExp_RK_ENABLED
            call random_number(smaller)
            call random_number(larger)
            call setMinMax(smaller, larger)
            logSubExp_ref = larger + log(1._RK + exp(smaller - larger))
            logSubExp = getLogSubExp(smaller = smaller, larger = larger, cenabled = cenabled)
            diff = logSubExp - logSubExp_ref
            assertion = diff < TOL
#endif
            call report()
            call test%assert(assertion, desc = "The logSubExp must be computed correctly for two non-overflowing numbers.")

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Compute the logSubExp for non-overflowing numbers
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            smaller = -1.e30 * (2*smaller - 1)
            larger = +1.e30 * (2*larger - 1)
            call setMinMax(smaller, larger)

#if         test_getLogSubExp_CK_ENABLED
            logSubExp_ref = larger + log(1._CK + exp(smaller - larger))
            logSubExp = getLogSubExp(smaller = smaller, larger = larger, cenabled = cenabled)
            diff = logSubExp%re - logSubExp_ref%re
            assertion = diff%re < TOL%re
#elif       test_getLogSubExp_RK_ENABLED
            logSubExp_ref = larger + log(1._RK + exp(smaller - larger))
            logSubExp = getLogSubExp(smaller = smaller, larger = larger, cenabled = cenabled)
            diff = logSubExp - logSubExp_ref
            assertion = diff < TOL
#endif
            call report()
            call test%assert(assertion, desc = "The logSubExp must be computed correctly for two overflowing numbers.")

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "larger         ", larger
                write(test%disp%unit,"(*(g0,:,', '))") "smaller        ", smaller
                write(test%disp%unit,"(*(g0,:,', '))") "logSubExp_ref  ", logSubExp_ref
                write(test%disp%unit,"(*(g0,:,', '))") "logSubExp      ", logSubExp
                write(test%disp%unit,"(*(g0,:,', '))") "diff           ", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL            ", TOL
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine
