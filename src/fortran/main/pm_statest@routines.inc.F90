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
!>  This include file contains the implementation of procedures in [pm_statest](@ref pm_statest).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getProbKS_ENABLED || setProbKS_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKC) :: ess ! effective sample size.
#if     getProbKS_ENABLED
        real(TKC) :: quanKS
#endif
#if     WIX_ENABLED
        CHECK_ASSERTION(__LINE__, 0 < weisum1, SK_"@getProbKS(): The condition `0 < weisum1` must hold. weisum1 = "//getStr(weisum1))
        ess = sqrt(real(weisum1, TKC))
#elif   WRX_ENABLED
        CHECK_ASSERTION(__LINE__, 0 < weisum1, SK_"@getProbKS(): The condition `0 < weisum1` must hold. weisum1 = "//getStr(weisum1))
        CHECK_ASSERTION(__LINE__, 0 < wsqsum1, SK_"@getProbKS(): The condition `0 < wsqsum1` must hold. wsqsum1 = "//getStr(wsqsum1))
        ess = weisum1 / sqrt(wsqsum1)
#elif   WII_ENABLED
        real(TKC) :: ess1, ess2
        CHECK_ASSERTION(__LINE__, 0 < weisum1, SK_"@getProbKS(): The condition `0 < weisum1` must hold. weisum1 = "//getStr(weisum1))
        CHECK_ASSERTION(__LINE__, 0 < weisum2, SK_"@getProbKS(): The condition `0 < weisum2` must hold. weisum2 = "//getStr(weisum2))
        ess1 = real(weisum1, TKC)
        ess2 = real(weisum2, TKC)
        ess = sqrt(ess1 * ess2 / (ess1 + ess2))
#elif   WRI_ENABLED
        real(TKC) :: ess1, ess2
        CHECK_ASSERTION(__LINE__, 0 < weisum1, SK_"@getProbKS(): The condition `0 < weisum1` must hold. weisum1 = "//getStr(weisum1))
        CHECK_ASSERTION(__LINE__, 0 < wsqsum1, SK_"@getProbKS(): The condition `0 < wsqsum1` must hold. wsqsum1 = "//getStr(wsqsum1))
        CHECK_ASSERTION(__LINE__, 0 < weisum2, SK_"@getProbKS(): The condition `0 < weisum2` must hold. weisum2 = "//getStr(weisum2))
        ess1 = weisum1**2 / wsqsum1
        ess2 = real(weisum2, TKC)
        ess = sqrt(ess1 * ess2 / (ess1 + ess2))
#elif   WRR_ENABLED
        real(TKC) :: ess1, ess2
        CHECK_ASSERTION(__LINE__, 0 < weisum1, SK_"@getProbKS(): The condition `0 < weisum1` must hold. weisum1 = "//getStr(weisum1))
        CHECK_ASSERTION(__LINE__, 0 < wsqsum1, SK_"@getProbKS(): The condition `0 < wsqsum1` must hold. wsqsum1 = "//getStr(wsqsum1))
        CHECK_ASSERTION(__LINE__, 0 < weisum2, SK_"@getProbKS(): The condition `0 < weisum2` must hold. weisum2 = "//getStr(weisum2))
        CHECK_ASSERTION(__LINE__, 0 < wsqsum2, SK_"@getProbKS(): The condition `0 < wsqsum2` must hold. wsqsum2 = "//getStr(wsqsum2))
        ess1 = weisum1**2 / wsqsum1
        ess2 = weisum2**2 / wsqsum2
        ess = sqrt(ess1 * ess2 / (ess1 + ess2))
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, 0 <= statKS, SK_"@getProbKS(): The condition `0 <= statKS` must hold. statKS = "//getStr(statKS))
        quanKS = statKS * (ess + 0.12_TKC + 0.11_TKC / ess)
        call setKolmCDF(probKS, quanKS)
        probKS = 1._TKC - probKS
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif