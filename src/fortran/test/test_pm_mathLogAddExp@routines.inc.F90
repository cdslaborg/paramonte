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
!>  This include file contains procedure implementations of the tests of [pm_mathLogAddExp](@ref pm_mathLogAddExp).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK_ENABLED
#define TYPEC complex
        integer, parameter :: TKG = CKG
#elif   RK_ENABLED
#define TYPEC real
        integer, parameter :: TKG = RKG
#else
#error  "Unrecognized interface."
#endif
        TYPEC(TKG)  , parameter :: TOL = epsilon(1._TKG) * 10_IK
        TYPEC(TKG)              :: logAddExp_ref
        TYPEC(TKG)              :: logAddExp
        TYPEC(TKG)              :: smaller
        TYPEC(TKG)              :: larger
        TYPEC(TKG)              :: diff

        assertion = .true._LK
        call runTest(cenabled = .true._LK)
        call runTest(cenabled = .false._LK)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTest(cenabled)

            logical(LK), intent(in) :: cenabled

            call setUnifRand(larger, -huge(smaller) / 10, huge(larger) / 10)
            call setUnifRand(smaller, -huge(smaller) / 10, huge(smaller) / 10)
            call setMinMax(smaller, larger)

            logAddExp_ref = larger + log(1._TKG + exp(smaller - larger))

            if (present(cenabled)) then
                logAddExp = getLogAddExp(smaller, larger, cenabled)
            else
                logAddExp = getLogAddExp(smaller, larger)
            end if

            diff = logAddExp - logAddExp_ref
            assertion = real(diff, TKG) < real(TOL, TKG)

            call report()
            call test%assert(assertion, desc = "The logAddExp must be computed correctly for two non-overflowing numbers.")

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "larger         ", larger
                write(test%disp%unit,"(*(g0,:,', '))") "smaller        ", smaller
                write(test%disp%unit,"(*(g0,:,', '))") "logAddExp_ref  ", logAddExp_ref
                write(test%disp%unit,"(*(g0,:,', '))") "logAddExp      ", logAddExp
                write(test%disp%unit,"(*(g0,:,', '))") "diff           ", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL            ", TOL
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

#undef  TYPEC