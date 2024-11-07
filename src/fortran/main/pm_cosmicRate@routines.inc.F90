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
!>  This file contains implementations of procedures [pm_cosmicRate](@ref pm_cosmicRate).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the array indexing rules.
#if     D1_ENABLED
        integer(IK) :: i
#define GET_ELEMENT(Array) Array(i)
#elif   D0_ENABLED
#define GET_ELEMENT(Array) Array
#define ALL
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getLogRateDensity_ENABLED && (H06_ENABLED || L08_ENABLED || B10_ENABLED || P15_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     H06_ENABLED
        real(RKG), parameter :: G0 = +3.4_RKG
        real(RKG), parameter :: G1 = -0.3_RKG
        real(RKG), parameter :: G2 = -7.8_RKG
        real(RKG), parameter :: LOGZ0PLUS1 = log(1._RKG + 0.97_RKG)
        real(RKG), parameter :: LOGZ1PLUS1 = log(1._RKG + 4.50_RKG)
        real(RKG), parameter :: LOG_NORM_FAC_1 = LOGZ0PLUS1 * (G0 - G1)
        real(RKG), parameter :: LOG_NORM_FAC_2 = LOGZ1PLUS1 * (G1 - G2) + LOG_NORM_FAC_1
#elif   L08_ENABLED
        real(RKG), parameter :: G0 = +3.3000_RKG
        real(RKG), parameter :: G1 = +0.0549_RKG
        real(RKG), parameter :: G2 = -4.4600_RKG
        real(RKG), parameter :: LOGZ0PLUS1 = log(1._RKG + 0.993_RKG)
        real(RKG), parameter :: LOGZ1PLUS1 = log(1._RKG + 3.800_RKG)
        real(RKG), parameter :: LOG_NORM_FAC_1 = LOGZ0PLUS1 * (G0 - G1)
        real(RKG), parameter :: LOG_NORM_FAC_2 = LOGZ1PLUS1 * (G1 - G2) + LOG_NORM_FAC_1
#elif   B10_ENABLED
        real(RKG), parameter :: G0 = +3.14_RKG
        real(RKG), parameter :: G1 = +1.36_RKG
        real(RKG), parameter :: G2 = -2.92_RKG
        real(RKG), parameter :: LOGZ0PLUS1 = log(1._RKG + 0.97_RKG)
        real(RKG), parameter :: LOGZ1PLUS1 = log(1._RKG + 4.00_RKG)
        real(RKG), parameter :: LOG_NORM_FAC_1 = LOGZ0PLUS1 * (G0 - G1)
        real(RKG), parameter :: LOG_NORM_FAC_2 = LOGZ1PLUS1 * (G1 - G2) + LOG_NORM_FAC_1
#elif   P15_ENABLED
        real(RKG), parameter :: EXPONENT_HIGH_Z = -7.8_RKG
        real(RKG), parameter :: LOGZ1PLUS1 = log(1._RKG + 4.5_RKG)
        real(RKG), parameter :: LOG_NORM_FAC_2 = -EXPONENT_HIGH_Z * LOGZ1PLUS1
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, all([logzplus1 >= 0._RKG]), SK_": The condition `logzplus1 >= 0.` must hold. logzplus1 = "//getStr(logzplus1)) ! fpp
#if     D1_ENABLED
        do concurrent(i = 1_IK : size(logzplus1, kind = IK))
#endif
#if         H06_ENABLED || L08_ENABLED || B10_ENABLED
            if (GET_ELEMENT(logzplus1) < LOGZ0PLUS1) then
                GET_ELEMENT(logRateDensity) = GET_ELEMENT(logzplus1) * G0
            elseif (GET_ELEMENT(logzplus1) < LOGZ1PLUS1) then
                GET_ELEMENT(logRateDensity) = GET_ELEMENT(logzplus1) * G1 + LOG_NORM_FAC_1
            else
                GET_ELEMENT(logRateDensity) = GET_ELEMENT(logzplus1) * G2 + LOG_NORM_FAC_2
            end if
#elif       P15_ENABLED
            if (GET_ELEMENT(logzplus1) < LOGZ1PLUS1) then
                GET_ELEMENT(logRateDensity) = 0._RKG
            else
                GET_ELEMENT(logRateDensity) = GET_ELEMENT(logzplus1) * EXPONENT_HIGH_Z + LOG_NORM_FAC_2
            end if
#endif
#if     D1_ENABLED
        end do
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getLogRateDensity_ENABLED && (M14_ENABLED || M17_ENABLED || F18_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     M14_ENABLED
        real(RKG), parameter :: LOG_AMPLITUDE = log(0.015_RKG)
        real(RKG), parameter :: EXPLONENT_LOWER = 2.7_RKG
        real(RKG), parameter :: EXPLONENT_UPPER = 5.6_RKG
        real(RKG), parameter :: ZPLUS1_BREAK = 2.9_RKG
        real(RKG), parameter :: ZPLUS1_COEFF = 1._RKG / (ZPLUS1_BREAK**EXPLONENT_UPPER)
#elif   M17_ENABLED
        real(RKG), parameter :: LOG_AMPLITUDE = log(0.01_RKG)
        real(RKG), parameter :: EXPLONENT_LOWER = 2.6_RKG
        real(RKG), parameter :: EXPLONENT_UPPER = 6.2_RKG
        real(RKG), parameter :: ZPLUS1_BREAK = 3.2_RKG
        real(RKG), parameter :: ZPLUS1_COEFF = 1._RKG / (ZPLUS1_BREAK**EXPLONENT_UPPER)
#elif   F18_ENABLED
        real(RKG), parameter :: LOG_AMPLITUDE = log(0.013_RKG)
        real(RKG), parameter :: EXPLONENT_LOWER = 2.99_RKG
        real(RKG), parameter :: EXPLONENT_UPPER = 6.19_RKG
        real(RKG), parameter :: ZPLUS1_BREAK = 2.63_RKG
        real(RKG), parameter :: ZPLUS1_COEFF = 1._RKG / (ZPLUS1_BREAK**EXPLONENT_UPPER)
#endif
        CHECK_ASSERTION(__LINE__, zplus1 >= 1._RKG, SK_": The condition `zplus1 >= 1.` must hold. logzplus1 = "//getStr(zplus1)) ! fpp
        CHECK_ASSERTION(__LINE__, logzplus1 >= 0._RKG, SK_": The condition `logzplus1 >= 0.` must hold. logzplus1 = "//getStr(logzplus1)) ! fpp
        logRateDensity = LOG_AMPLITUDE + EXPLONENT_LOWER * logzplus1 - log(1._RKG + ZPLUS1_COEFF * zplus1**EXPLONENT_UPPER)
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef GET_ELEMENT