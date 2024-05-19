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
!>  This include file contains procedure implementations of the tests of [pm_cosmicRate](@ref pm_cosmicRate).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%
#if     getLogRateDensity_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG)   , parameter :: TOL = epsilon(0._RKG)*100
        real(RKG)               :: zplus1, logzplus1
        real(RKG)               :: logRateDensity_ref
        real(RKG)               :: logRateDensity
        real(RKG)               :: diff
        integer(IK) :: i
        assertion = .true._LK
        do i = 1_IK, 500_IK
            logzplus1 = getUnifRand(0._RKG, 10._RKG)
            zplus1 = exp(logzplus1)
#if         H06_ENABLED
            logRateDensity = getLogRateDensityH06(logzplus1)
#elif       L08_ENABLED
            logRateDensity = getLogRateDensityL08(logzplus1)
#elif       B10_ENABLED
            logRateDensity = getLogRateDensityB10(logzplus1)
#elif       P15_ENABLED
            logRateDensity = getLogRateDensityP15(logzplus1)
#elif       M14_ENABLED
            logRateDensity = getLogRateDensityM14(zplus1, logzplus1)
#elif       M17_ENABLED
            logRateDensity = getLogRateDensityM17(zplus1, logzplus1)
#elif       F18_ENABLED
            logRateDensity = getLogRateDensityF18(zplus1, logzplus1)
#else
#error      "Unrecognized interface."
#endif
            call report(int(__LINE__, IK))
        end do

    contains

#if     H06_ENABLED || L08_ENABLED || B10_ENABLED || P15_ENABLED

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pure elemental function getLogRateDensity_ref(logzplus1) result(logRateDensity)
            real(RKG), intent(in)   :: logzplus1
            real(RKG)               :: logRateDensity
#if         H06_ENABLED
            real(RKG), parameter :: G0 = +3.4_RKG
            real(RKG), parameter :: G1 = -0.3_RKG
            real(RKG), parameter :: G2 = -7.8_RKG
            real(RKG), parameter :: LOGZ0PLUS1 = log(1._RKG + 0.97_RKG)
            real(RKG), parameter :: LOGZ1PLUS1 = log(1._RKG + 4.50_RKG)
            real(RKG), parameter :: LOG_NORM_FAC_1 = LOGZ0PLUS1 * (G0 - G1)
            real(RKG), parameter :: LOG_NORM_FAC_2 = LOGZ1PLUS1 * (G1 - G2) + LOG_NORM_FAC_1
#elif       L08_ENABLED
            real(RKG), parameter :: G0 = +3.3000_RKG
            real(RKG), parameter :: G1 = +0.0549_RKG
            real(RKG), parameter :: G2 = -4.4600_RKG
            real(RKG), parameter :: LOGZ0PLUS1 = log(1._RKG + 0.993_RKG)
            real(RKG), parameter :: LOGZ1PLUS1 = log(1._RKG + 3.800_RKG)
            real(RKG), parameter :: LOG_NORM_FAC_1 = LOGZ0PLUS1 * (G0 - G1)
            real(RKG), parameter :: LOG_NORM_FAC_2 = LOGZ1PLUS1 * (G1 - G2) + LOG_NORM_FAC_1
#elif       B10_ENABLED
            real(RKG), parameter :: G0 = +3.14_RKG
            real(RKG), parameter :: G1 = +1.36_RKG
            real(RKG), parameter :: G2 = -2.92_RKG
            real(RKG), parameter :: LOGZ0PLUS1 = log(1._RKG + 0.97_RKG)
            real(RKG), parameter :: LOGZ1PLUS1 = log(1._RKG + 4.00_RKG)
            real(RKG), parameter :: LOG_NORM_FAC_1 = LOGZ0PLUS1 * (G0 - G1)
            real(RKG), parameter :: LOG_NORM_FAC_2 = LOGZ1PLUS1 * (G1 - G2) + LOG_NORM_FAC_1
#elif       P15_ENABLED
            real(RKG), parameter :: EXPONENT_HIGH_Z = -7.8_RKG
            real(RKG), parameter :: LOGZ1PLUS1 = log(1._RKG + 4.5_RKG)
            real(RKG), parameter :: LOG_NORM_FAC_2 = -EXPONENT_HIGH_Z * LOGZ1PLUS1
#else
#error      "Unrecognized interface."
#endif
#if         H06_ENABLED || L08_ENABLED || B10_ENABLED
            if (logzplus1 < LOGZ0PLUS1) then
                logRateDensity = logzplus1 * G0
            elseif (logzplus1 < LOGZ1PLUS1) then
                logRateDensity = logzplus1 * G1 + LOG_NORM_FAC_1
            else
                logRateDensity = logzplus1 * G2 + LOG_NORM_FAC_2
            end if
#elif       P15_ENABLED
            if (logzplus1 < LOGZ1PLUS1) then
                logRateDensity = 0._RKG
            else
                logRateDensity = logzplus1 * EXPONENT_HIGH_Z + LOG_NORM_FAC_2
            end if
#endif
        end function
#elif   M14_ENABLED || M17_ENABLED || F18_ENABLED
        pure elemental function getLogRateDensity_ref(zplus1, logzplus1) result(logRateDensity)
            real(RKG), intent(in)   :: zplus1, logzplus1
            real(RKG)               :: logRateDensity
#if         M14_ENABLED
            real(RKG), parameter :: LOG_AMPLITUDE = log(0.015_RKG)
            real(RKG), parameter :: EXPLONENT_LOWER = 2.7_RKG
            real(RKG), parameter :: EXPLONENT_UPPER = 5.6_RKG
            real(RKG), parameter :: ZPLUS1_BREAK = 2.9_RKG
            real(RKG), parameter :: ZPLUS1_COEFF = 1._RKG / (ZPLUS1_BREAK**EXPLONENT_UPPER)
#elif       M17_ENABLED
            real(RKG), parameter :: LOG_AMPLITUDE = log(0.01_RKG)
            real(RKG), parameter :: EXPLONENT_LOWER = 2.6_RKG
            real(RKG), parameter :: EXPLONENT_UPPER = 6.2_RKG
            real(RKG), parameter :: ZPLUS1_BREAK = 3.2_RKG
            real(RKG), parameter :: ZPLUS1_COEFF = 1._RKG / (ZPLUS1_BREAK**EXPLONENT_UPPER)
#elif       F18_ENABLED
            real(RKG), parameter :: LOG_AMPLITUDE = log(0.013_RKG)
            real(RKG), parameter :: EXPLONENT_LOWER = 2.99_RKG
            real(RKG), parameter :: EXPLONENT_UPPER = 6.19_RKG
            real(RKG), parameter :: ZPLUS1_BREAK = 2.63_RKG
            real(RKG), parameter :: ZPLUS1_COEFF = 1._RKG / (ZPLUS1_BREAK**EXPLONENT_UPPER)
#endif
            logRateDensity = LOG_AMPLITUDE + EXPLONENT_LOWER * logzplus1 - log(1._RKG + ZPLUS1_COEFF * zplus1**EXPLONENT_UPPER)
        end function
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line)
            integer(IK), intent(in) :: line
#if         H06_ENABLED || L08_ENABLED || B10_ENABLED || P15_ENABLED
            logRateDensity_ref = getLogRateDensity_ref(logzplus1)
#elif       M14_ENABLED || M17_ENABLED || F18_ENABLED
            logRateDensity_ref = getLogRateDensity_ref(zplus1, logzplus1)
#else
#error      "Unrecognized interface."
#endif
            diff = abs(logRateDensity - logRateDensity_ref)
            assertion = assertion .and. logical(diff <= TOL * abs(logRateDensity_ref), LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "logRateDensity_ref ", logRateDensity_ref
                write(test%disp%unit,"(*(g0,:,', '))") "logRateDensity     ", logRateDensity
                write(test%disp%unit,"(*(g0,:,', '))") "logzplus1          ", logzplus1
                write(test%disp%unit,"(*(g0,:,', '))") "zplus1             ", zplus1
                write(test%disp%unit,"(*(g0,:,', '))") "diff               ", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL                ", TOL
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, desc = "The log rate density must be computed correctly for the specified logzplus1.", line = line)
        end subroutine

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif