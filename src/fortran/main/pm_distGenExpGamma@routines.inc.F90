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
!>  This include file contains the implementation of procedures in [pm_distGenExpGamma](@ref pm_distGenExpGamma).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 12:20 PM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getGenExpGammaLogPDFNF_ENABLED && KD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, kappa > 0._RKG, SK_"@getGenExpGammaLogPDFNF(): The condition `kappa > 0.` must hold. kappa = "//getStr(kappa)) ! fpp
        logPDFNF = -log_gamma(kappa)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getGenExpGammaLogPDFNF_ENABLED && KO_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        CHECK_ASSERTION(__LINE__, invOmega > 0._RKG, SK_"@getGenExpGammaLogPDFNF(): The condition `invOmega > 0.` must hold. invOmega = "//getStr(invOmega)) ! fpp
        logPDFNF = getGenExpGammaLogPDFNF(kappa)
        if (invOmega /= 1._RKG) logPDFNF = logPDFNF + log(invOmega)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getGenExpGammaLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: kappa_def, invOmega_def, logSigma_def
        kappa_def = 1._RKG; if (present(kappa)) kappa_def = kappa
        invOmega_def = 1._RKG; if (present(invOmega)) invOmega_def = invOmega
        logSigma_def = 0._RKG; if (present(logSigma)) logSigma_def = logSigma
        call setGenExpGammaLogPDF(logPDF, x, getGenExpGammaLogPDFNF(kappa_def, invOmega_def), kappa_def, invOmega_def, logSigma_def)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGenExpGammaLogPDF_ENABLED && DDDD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        logPDF = x - exp(x)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGenExpGammaLogPDF_ENABLED && NKDD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, kappa > 0._RKG, SK_"@getGenExpGammaLogPDFNF(): The condition `kappa > 0.` must hold. kappa = "//getStr(kappa)) ! fpp
        CHECK_ASSERTION(__LINE__, abs(getGenExpGammaLogPDFNF(kappa) - logPDFNF) <= 100 * epsilon(0._RKG), \
        SK_"@setGenExpGammaLogPDF(): The condition `abs(getGenExpGammaLogPDFNF(kappa) - logPDFNF) <= 100 * epsilon(0._RKG)` must hold. getGenExpGammaLogPDFNF(kappa), logPDFNF = " \
        //getStr([getGenExpGammaLogPDFNF(kappa), logPDFNF])) ! fpp
        logPDF = logPDFNF + kappa * x - exp(x)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGenExpGammaLogPDF_ENABLED && NKOD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: xscaled
        real(RKG), parameter :: LOG_SQRT_HUGE = log(sqrt(huge(0._RKG)))
        xscaled = x * invOmega
        CHECK_ASSERTION(__LINE__, kappa > 0._RKG, SK_"@getGenExpGammaLogPDFNF(): The condition `kappa > 0.` must hold. kappa = "//getStr(kappa)) ! fpp
        CHECK_ASSERTION(__LINE__, invOmega > 0._RKG, SK_"@setGenExpGammaLogPDF(): The condition `invOmega > 0.` must hold. invOmega = "//getStr(invOmega)) ! fpp
        CHECK_ASSERTION(__LINE__, abs(getGenExpGammaLogPDFNF(kappa, invOmega) - logPDFNF) <= 100 * epsilon(0._RKG), \
        SK_"@setGenExpGammaLogPDF(): The condition `abs(getGenExpGammaLogPDFNF(kappa, invOmega) - logPDFNF) <= 100 * epsilon(0._RKG)` must hold. getGenExpGammaLogPDFNF(kappa, invOmega), logPDFNF = " \
        //getStr([getGenExpGammaLogPDFNF(kappa, invOmega), logPDFNF])) ! fpp
        if (xscaled < LOG_SQRT_HUGE) then
            logPDF = logPDFNF + kappa * xscaled - exp(xscaled)
        else
            logPDF = -LOG_SQRT_HUGE
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGenExpGammaLogPDF_ENABLED && NKOS_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, abs(getGenExpGammaLogPDFNF(kappa, invOmega) - logPDFNF) <= 100 * epsilon(0._RKG), \
        SK_"@setGenExpGammaLogPDF(): The condition `abs(getGenExpGammaLogPDFNF(kappa, invOmega) - logPDFNF) <= 100 * epsilon(0._RKG)` must hold. getGenExpGammaLogPDFNF(kappa, invOmega), logPDFNF = " \
        //getStr([getGenExpGammaLogPDFNF(kappa, invOmega), logPDFNF])) ! fpp
        call setGenExpGammaLogPDF(logPDF, x - logSigma, logPDFNF, kappa, invOmega)

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getGenExpGammaCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: info
        real(RKG) :: xnormed
        if (present(logSigma)) then
            xnormed = x - logSigma
        else
            xnormed = x
        end if
        if (present(invOmega)) xnormed = xnormed * invOmega
        if (present(kappa)) then
            call setGenExpGammaCDF(cdf, xnormed, kappa, info)
        else
            call setGenExpGammaCDF(cdf, xnormed, info)
        end if
        if (info < 0_IK) error stop MODULE_NAME//SK_"@getGenExpGammaCDF(): The computation of the regularized Lower Incomplete Gamma function failed. This can happen if `kappa` is too large."

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGenExpGammaCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: gamincupp
#if     DDD_ENABLED
        real(RKG), parameter :: kappa = 1._RKG!, logGammaKappa = log_gamma(kappa)
        !call setGammaIncLow(cdf, exp(x), logGammaKappa, kappa, info)
        call setGammaInc(cdf, gamincupp, exp(x), kappa, info)
#else
        CHECK_ASSERTION(__LINE__, kappa > 0._RKG, SK_"@setGenExpGammaCDF(): The condition `kappa > 0.` must hold. kappa = "//getStr(kappa)) ! fpp
        !check_assertion(__LINE__, abs(log_gamma(kappa) - logGammaKappa) < 100 * epsilon(0._RKG), SK_"@setGenExpGammaCDF(): The condition `abs(log_gamma(kappa) - logGammaKappa) < 100 * epsilon(0._RKG)` must hold. log_gamma(kappa), logGammaKappa = "//getStr([log_gamma(kappa), logGammaKappa])) ! fpp
#if     KDD_ENABLED
        !call setGammaIncLow(cdf, exp(x), logGammaKappa, kappa, info)
        call setGammaInc(cdf, gamincupp, exp(x), kappa, info)
#else
        CHECK_ASSERTION(__LINE__, invOmega > 0._RKG, SK_"@setGenExpGammaCDF(): The condition `invOmega > 0.` must hold. invOmega = "//getStr(invOmega)) ! fpp
#if     KOD_ENABLED
        !call setGammaIncLow(cdf, exp(x * invOmega), logGammaKappa, kappa, info)
        call setGammaInc(cdf, gamincupp, exp(x * invOmega), kappa, info)
#elif   KOS_ENABLED
        !call setGammaIncLow(cdf, exp((x - logSigma) * invOmega), logGammaKappa, kappa, info)
        call setGammaInc(cdf, gamincupp, exp((x - logSigma) * invOmega), kappa, info)
#else
#error  "Unrecognized interface."
#endif
#endif
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif