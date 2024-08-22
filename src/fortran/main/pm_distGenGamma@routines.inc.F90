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
!>  This include file contains the implementation of procedures in [pm_distGenGamma](@ref pm_distGenGamma).
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 12:20 PM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getGenGammaLogPDFNF_ENABLED && KDD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, kappa > 0._RKG, SK_"@getGenGammaLogPDFNF(): The condition `kappa > 0.` must hold. kappa = "//getStr(kappa)) ! fpp
        logPDFNF = -log_gamma(kappa)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getGenGammaLogPDFNF_ENABLED && KOD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, invOmega > 0._RKG, SK_"@getGenGammaLogPDFNF(): The condition `invOmega > 0.` must hold. invOmega = "//getStr(invOmega)) ! fpp
        logPDFNF = getGenGammaLogPDFNF(kappa)
        if (invOmega /= 1._RKG) logPDFNF = logPDFNF + log(invOmega)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getGenGammaLogPDFNF_ENABLED && KOS_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, invSigma > 0._RKG, SK_"@getGenGammaLogPDFNF(): The condition `invSigma > 0.` must hold. invSigma = "//getStr(invSigma)) ! fpp
        logPDFNF = getGenGammaLogPDFNF(kappa, invOmega)
        if (invSigma /= 1._RKG) logPDFNF = logPDFNF + log(invSigma)

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getGenGammaLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: kappa_def, invOmega_def, invSigma_def
        kappa_def = 1._RKG; if (present(kappa)) kappa_def = kappa
        invOmega_def = 1._RKG; if (present(invOmega)) invOmega_def = invOmega
        invSigma_def = 1._RKG; if (present(invSigma)) invSigma_def = invSigma
        call setGenGammaLogPDF(logPDF, x, getGenGammaLogPDFNF(kappa_def, invOmega_def, invSigma_def), kappa_def, invOmega_def, invSigma_def)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGenGammaLogPDF_ENABLED && DDDD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, x > 0._RKG, SK_"@getGenGammaLogPDFNF(): The condition `x > 0.` must hold. x = "//getStr(x)) ! fpp
        logPDF = -x

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGenGammaLogPDF_ENABLED && NKDD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, x > 0._RKG, SK_"@getGenGammaLogPDFNF(): The condition `x > 0.` must hold. x = "//getStr(x)) ! fpp
        CHECK_ASSERTION(__LINE__, kappa > 0._RKG, SK_"@getGenGammaLogPDFNF(): The condition `kappa > 0.` must hold. kappa = "//getStr(kappa)) ! fpp
        CHECK_ASSERTION(__LINE__, abs(getGenGammaLogPDFNF(kappa) - logPDFNF) <= 100 * epsilon(0._RKG), \
        SK_"@setGenGammaLogPDF(): The condition `abs(getGenGammaLogPDFNF(kappa) - logPDFNF) <= 100 * epsilon(0._RKG)` must hold. getGenGammaLogPDFNF(kappa), logPDFNF = " \
        //getStr([getGenGammaLogPDFNF(kappa), logPDFNF])) ! fpp
        logPDF = logPDFNF + log(x) * (kappa - 1._RKG) - x

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGenGammaLogPDF_ENABLED && NKOD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, x > 0._RKG, SK_"@getGenGammaLogPDFNF(): The condition `x > 0.` must hold. x = "//getStr(x)) ! fpp
        CHECK_ASSERTION(__LINE__, kappa > 0._RKG, SK_"@getGenGammaLogPDFNF(): The condition `kappa > 0.` must hold. kappa = "//getStr(kappa)) ! fpp
        CHECK_ASSERTION(__LINE__, invOmega > 0._RKG, SK_"@setGenGammaLogPDF(): The condition `invOmega > 0.` must hold. invOmega = "//getStr(invOmega)) ! fpp
        CHECK_ASSERTION(__LINE__, abs(getGenGammaLogPDFNF(kappa, invOmega) - logPDFNF) <= 100 * epsilon(0._RKG), \
        SK_"@setGenGammaLogPDF(): The condition `abs(getGenGammaLogPDFNF(kappa, invOmega) - logPDFNF) <= 100 * epsilon(0._RKG)` must hold. getGenGammaLogPDFNF(kappa, invOmega), logPDFNF = " \
        //getStr([getGenGammaLogPDFNF(kappa, invOmega), logPDFNF])) ! fpp
        logPDF = logPDFNF + log(x) * (kappa * invOmega - 1._RKG) - x**invOmega

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGenGammaLogPDF_ENABLED && NKOS_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: xscaled
        CHECK_ASSERTION(__LINE__, x > 0._RKG, SK_"@getGenGammaLogPDFNF(): The condition `x > 0.` must hold. x = "//getStr(x)) ! fpp
        CHECK_ASSERTION(__LINE__, kappa > 0._RKG, SK_"@getGenGammaLogPDFNF(): The condition `kappa > 0.` must hold. kappa = "//getStr(kappa)) ! fpp
        CHECK_ASSERTION(__LINE__, invOmega > 0._RKG, SK_"@setGenGammaLogPDF(): The condition `invOmega > 0.` must hold. invOmega = "//getStr(invOmega)) ! fpp
        CHECK_ASSERTION(__LINE__, invSigma > 0._RKG, SK_"@setGenGammaLogPDF(): The condition `invSigma > 0.` must hold. invSigma = "//getStr(invSigma)) ! fpp
        CHECK_ASSERTION(__LINE__, abs(getGenGammaLogPDFNF(kappa, invOmega, invSigma) - logPDFNF) <= 100 * epsilon(0._RKG), \
        SK_"@setGenGammaLogPDF(): The condition `abs(getGenGammaLogPDFNF(kappa, invOmega, invSigma) - logPDFNF) <= 100 * epsilon(0._RKG)` must hold. getGenGammaLogPDFNF(kappa, invOmega, invSigma), logPDFNF = " \
        //getStr([getGenGammaLogPDFNF(kappa, invOmega, invSigma), logPDFNF])) ! fpp
        xscaled = x * invSigma
        logPDF = logPDFNF + log(xscaled) * (kappa * invOmega - 1._RKG) - xscaled**invOmega

        !%%%%%%%%%%%%%%%%%%%%%
#elif   getGenGammaCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: info
        real(RKG) :: xnormed
        if (present(invSigma)) then
            xnormed = x * invSigma
        else
            xnormed = x
        end if
        if (present(invSigma)) xnormed = xnormed ** invOmega
        if (present(kappa)) then
            call setGenGammaCDF(cdf, xnormed, kappa, info)
        else
            call setGenGammaCDF(cdf, xnormed, info)
        end if
        if (info < 0_IK) error stop MODULE_NAME//SK_"@getGenGammaCDF(): The computation of the regularized Lower Incomplete Gamma function failed. This can happen if `kappa` is too large."

        !%%%%%%%%%%%%%%%%%%%%%
#elif   setGenGammaCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: gamincupp
#if     DDD_ENABLED
        real(RKG), parameter :: kappa = 1._RKG!, logGammaKappa = log_gamma(kappa)
        CHECK_ASSERTION(__LINE__, x > 0._RKG, SK_"@setGenGammaCDF(): The condition `x > 0.` must hold. x = "//getStr(x)) ! fpp
        !call setGammaIncLow(cdf, x, logGammaKappa, kappa, info, tol)
        call setGammaInc(cdf, gamincupp, x, kappa, info)
#else
        CHECK_ASSERTION(__LINE__, x > 0._RKG, SK_"@setGenGammaCDF(): The condition `x > 0.` must hold. x = "//getStr(x)) ! fpp
        CHECK_ASSERTION(__LINE__, kappa > 0._RKG, SK_"@setGenGammaCDF(): The condition `kappa > 0.` must hold. kappa = "//getStr(kappa)) ! fpp
        !check_assertion(__LINE__, abs(log_gamma(kappa) - logGammaKappa) < 100 * epsilon(0._RKG), SK_"@setGenGammaCDF(): The condition `abs(log_gamma(kappa) - logGammaKappa) < 100 * epsilon(0._RKG)` must hold. log_gamma(kappa), logGammaKappa = "//getStr([log_gamma(kappa), logGammaKappa])) ! fpp
#if     KDD_ENABLED
        !call setGammaIncLow(cdf, x, logGammaKappa, kappa, info, tol)
        call setGammaInc(cdf, gamincupp, x, kappa, info)
#else
        CHECK_ASSERTION(__LINE__, invOmega > 0._RKG, SK_"@setGenGammaCDF(): The condition `invOmega > 0.` must hold. invOmega = "//getStr(invOmega)) ! fpp
#if     KOD_ENABLED
        !call setGammaIncLow(cdf, x**invOmega, logGammaKappa, kappa, info, tol)
        call setGammaInc(cdf, gamincupp, x**invOmega, kappa, info)
#elif   KOS_ENABLED
        CHECK_ASSERTION(__LINE__, invSigma > 0._RKG, SK_"@setGenGammaCDF(): The condition `invSigma > 0.` must hold. invSigma = "//getStr(invSigma)) ! fpp
        !call setGammaIncLow(cdf, (x * invSigma)**invOmega, logGammaKappa, kappa, info, tol)
        call setGammaInc(cdf, gamincupp, (x * invSigma)**invOmega, kappa, info)
#else
#error  "Unrecognized interface."
#endif
#endif
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGenGammaRand_ENABLED && (RNGD_ENABLED || RNGF_ENABLED || RNGX_ENABLED) && KR_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Set the URNG.
#if     RNGD_ENABLED
#define RNG
#elif   RNGF_ENABLED || RNGX_ENABLED
#define RNG rng,
#elif   setGenGammaRand_ENABLED
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, omega > 0, SK_"@setGenGammaCDF(): The condition `omega > 0` must hold. omega = "//getStr(omega)) ! fpp
        call setGammaRand(RNG rand, kappa = kappa, sigma = sigma)
        rand = rand**omega
#undef  RNG

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif