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
!>  This include file contains the implementation of procedures in [pm_distExpGamma](@ref pm_distExpGamma).
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 12:20 PM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getExpGammaLogPDFNF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, kappa > 0._RKG, SK_"@getExpGammaLogPDFNF(): The condition `kappa > 0.` must hold. kappa = "//getStr(kappa))
        logPDFNF = -log_gamma(kappa)

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getExpGammaLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: kappa_def, logSigma_def
        kappa_def = 1._RKG; if (present(kappa)) kappa_def = kappa
        logSigma_def = 0._RKG; if (present(logSigma)) logSigma_def = logSigma
        call setExpGammaLogPDF(logPDF, x, getExpGammaLogPDFNF(kappa_def), kappa_def, logSigma_def)

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setExpGammaLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

#if     DDD_ENABLED
        logPDF = x - exp(x)
#elif   NKD_ENABLED
        CHECK_ASSERTION(__LINE__, 0._RKG < kappa, SK_"@getExpGammaLogPDFNF(): The condition `0. < kappa` must hold. kappa = "//getStr(kappa))
        logPDF = logPDFNF + kappa * x - exp(x)
#elif   NKS_ENABLED
        call setExpGammaLogPDF(logPDF, x - logSigma, logPDFNF, kappa)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%
#elif   getExpGammaCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: info
        real(RKG)   :: xnormed
        if (present(logSigma)) then
            xnormed = x - logSigma
        else
            xnormed = x
        end if
        if (present(kappa)) then
            call setExpGammaCDF(cdf, xnormed, kappa, info)
        else
            call setExpGammaCDF(cdf, xnormed, info)
        end if
        if (info < 0_IK) error stop MODULE_NAME//SK_"@getExpGammaCDF(): The computation of the regularized Lower Incomplete Gamma function failed. This can happen if `kappa` is too large."

        !%%%%%%%%%%%%%%%%%%%%%
#elif   setExpGammaCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: gamincupp
#if     DD_ENABLED
        real(RKG), parameter :: kappa = 1._RKG!, logGammaKappa = log_gamma(kappa)
        !call setGammaIncLow(cdf, gamincupp, exp(x), logGammaKappa, kappa, info)
        call setGammaInc(cdf, gamincupp, exp(x), kappa, info)
#else
        CHECK_ASSERTION(__LINE__, kappa > 0._RKG, SK_"@setExpGammaCDF(): The condition `kappa > 0.` must hold. kappa = "//getStr(kappa)) ! fpp
        !check_assertion(__LINE__, abs(log_gamma(kappa) - logGammaKappa) < 100 * epsilon(0._RKG), \
        !SK_"@setExpGammaCDF(): The condition `abs(log_gamma(kappa) - logGammaKappa) < 100 * epsilon(0._RKG)` must hold. log_gamma(kappa), logGammaKappa = "//\
        !getStr([log_gamma(kappa), logGammaKappa])) ! fpp
#if     KD_ENABLED
        !call setGammaIncLow(cdf, exp(x), logGammaKappa, kappa, info)
        call setGammaInc(cdf, gamincupp, exp(x), kappa, info)
#elif   KS_ENABLED
        !call setGammaIncLow(cdf, gamincupp, exp(x - logSigma), logGammaKappa, kappa, info)
        call setGammaInc(cdf, gamincupp, exp(x - logSigma), kappa, info)
#else
#error  "Unrecognized interface."
#endif
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif