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
!>  This include file contains the implementation of procedures in [pm_distNegExp](@ref pm_distNegExp).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%
#if     getNegExpLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        if (present(mu)) then
            if (present(invSigma)) then
                call setNegExpLogPDF(logPDF, x, mu, invSigma, log(invSigma))
            else
                call setNegExpLogPDF(logPDF, x, mu)
            end if
        else
            if (present(invSigma)) then
                call setNegExpLogPDF(logPDF, x, invSigma, log(invSigma))
            else
                call setNegExpLogPDF(logPDF, x)
            end if
        end if

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   setNegExpLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        ! Validate the input.

#if     DDD_ENABLED || DIL_ENABLED
        CHECK_ASSERTION(__LINE__, x <= 0._RKC, SK_"@setNegExpPDF(): The condition `x <= 0._RKC` must hold. x = "//getStr(x)) ! fpp
#elif   MDD_ENABLED || MIL_ENABLED
        CHECK_ASSERTION(__LINE__, x <= mu, SK_"@setNegExpPDF(): The condition `x <= mu` must hold. x = "//getStr([x, mu])) ! fpp
#else
#error  "Unrecognized interface."
#endif

#if     DIL_ENABLED || MIL_ENABLED
        CHECK_ASSERTION(__LINE__, 0._RKC < invSigma, SK_"@setNegExpPDF(): The condition `0._RKC < invSigma` must hold. invSigma = "//getStr(invSigma)) ! fpp
        CHECK_ASSERTION(__LINE__, abs(log(invSigma) - logInvSigma) < epsilon(0._RKC) * 10, SK_"@setNegExpPDF(): The condition `log(invSigma) == logInvSigma` must hold. invSigma, logInvSigma = "//getStr([invSigma, logInvSigma])) ! fpp
#elif   !(DDD_ENABLED || MDD_ENABLED)
#error  "Unrecognized interface."
#endif

        ! Compute the PDF.

#if     DDD_ENABLED
        logPDF = x
#elif   MDD_ENABLED
        logPDF = x - mu
#elif   DIL_ENABLED
        logPDF = logInvSigma + x * invSigma
#elif   MIL_ENABLED
        logPDF = logInvSigma + (x - mu) * invSigma
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%
#elif   getNegExpCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        if (present(mu)) then
            if (present(invSigma)) then
                call setNegExpCDF(cdf, x, mu, invSigma)
            else
                call setNegExpCDF(cdf, x - mu)
            end if
        else
            if (present(invSigma)) then
                call setNegExpCDF(cdf, x, invSigma)
            else
                call setNegExpCDF(cdf, x)
            end if
        end if

        !%%%%%%%%%%%%%%%%%%%
#elif   setNegExpCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        ! Validate the input.

#if     XDD_ENABLED || XDI_ENABLED
        CHECK_ASSERTION(__LINE__, x <= 0._RKC, SK_"@setNegExpCDF(): The condition `x <= 0._RKC` must hold. x = "//getStr(x)) ! fpp
#elif   XMI_ENABLED
        CHECK_ASSERTION(__LINE__, x <= mu, SK_"@setNegExpCDF(): The condition `x <= mu` must hold. x = "//getStr([x, mu])) ! fpp
#else
#error  "Unrecognized interface."
#endif

#if     XDI_ENABLED || XMI_ENABLED
        CHECK_ASSERTION(__LINE__, 0._RKC < invSigma, SK_"@setNegExpCDF(): The condition `0._RKC < invSigma` must hold. invSigma = "//getStr(invSigma)) ! fpp
#elif   !XDD_ENABLED
#error  "Unrecognized interface."
#endif

        ! Compute the CDF.

#if     XDD_ENABLED
        cdf = exp(x)
#elif   XDI_ENABLED
        cdf = exp(x * invSigma)
#elif   XMI_ENABLED
        cdf = exp((x - mu) * invSigma)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%
#elif   getNegExpRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        use pm_distUnif, only: getUnifRand
        call setNegExpRand(rand, getUnifRand(0._RKC, 1._RKC), sigma)
        if (present(mu)) rand = rand + mu

        !%%%%%%%%%%%%%%%%%%%%
#elif   setNegExpRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0._RKC < urand, SK_"@setNegExpRand(): The condition `0._RKC < urand` must hold. urand = "//getStr(urand)) ! fpp
        CHECK_ASSERTION(__LINE__, urand <= 1._RKC, SK_"@setNegExpRand(): The condition `urand <= 1._RKC` must hold. urand = "//getStr(urand)) ! fpp
#if     USD_ENABLED || USM_ENABLED
        CHECK_ASSERTION(__LINE__, 0._RKC < sigma, SK_"@setNegExpRand(): The condition `0._RKC < sigma` must hold. sigma = "//getStr(sigma)) ! fpp
#endif
#if     UDD_ENABLED
        rand = log(1._RKC - urand)
#elif   USD_ENABLED
        rand = log(1._RKC - urand) * sigma
#elif   USM_ENABLED
        rand = log(1._RKC - urand) * sigma + mu
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif