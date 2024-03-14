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
!>  This include file contains the implementation of procedures in [pm_distExpPDF](@ref pm_distExpPDF).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%
#if     getExpLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        if (present(mu)) then
            if (present(invSigma)) then
                call setExpLogPDF(logPDF, x, mu, invSigma, log(invSigma))
            else
                call setExpLogPDF(logPDF, x, mu)
            end if
        else
            if (present(invSigma)) then
                call setExpLogPDF(logPDF, x, invSigma, log(invSigma))
            else
                call setExpLogPDF(logPDF, x)
            end if
        end if

        !%%%%%%%%%%%%%%%%%%%
#elif   setExpLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        ! Validate the input.

#if     XDD_ENABLED || XDI_ENABLED
        CHECK_ASSERTION(__LINE__, 0 <= x, SK_"@setExpPDF(): The condition `0 <= x` must hold. x = "//getStr(x)) ! fpp
#elif   XMD_ENABLED || XMI_ENABLED
        CHECK_ASSERTION(__LINE__, mu <= x, SK_"@setExpPDF(): The condition `mu <= x` must hold. mu, x = "//getStr([mu, x])) ! fpp
#else
#error  "Unrecognized interface."
#endif

#if     XDI_ENABLED || XMI_ENABLED
        CHECK_ASSERTION(__LINE__, 0 < invSigma, SK_"@setExpPDF(): The condition `0 < invSigma` must hold. invSigma = "//getStr(invSigma)) ! fpp
        CHECK_ASSERTION(__LINE__, abs(log(invSigma) - logInvSigma) < epsilon(0._TKC) * 10, SK_"@setExpPDF(): The condition `log(invSigma) == logInvSigma` must hold. invSigma, logInvSigma = "//getStr([invSigma, logInvSigma])) ! fpp
#elif   !(XDD_ENABLED || XMD_ENABLED)
#error  "Unrecognized interface."
#endif

        ! Compute the PDF.

#if     XDD_ENABLED
        logPDF = -x
#elif   XMD_ENABLED
        logPDF = mu - x
#elif   XDI_ENABLED
        logPDF = logInvSigma - x * invSigma
#elif   XMI_ENABLED
        logPDF = logInvSigma - (x - mu) * invSigma
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%
#elif   getExpCDF_ENABLED
        !%%%%%%%%%%%%%%%%

        if (present(mu)) then
            if (present(invSigma)) then
                call setExpCDF(cdf, x, mu, invSigma)
            else
                call setExpCDF(cdf, x - mu)
            end if
        else
            if (present(invSigma)) then
                call setExpCDF(cdf, x, invSigma)
            else
                call setExpCDF(cdf, x)
            end if
        end if

        !%%%%%%%%%%%%%%%%
#elif   setExpCDF_ENABLED
        !%%%%%%%%%%%%%%%%

        ! Validate the input.

#if     XDD_ENABLED || XDI_ENABLED
        CHECK_ASSERTION(__LINE__, 0 <= x, SK_"@setExpCDF(): The condition `0 <= x` must hold. x = "//getStr(x)) ! fpp
#elif   XMI_ENABLED
        CHECK_ASSERTION(__LINE__, mu <= x, SK_"@setExpCDF(): The condition `mu <= x` must hold. mu, x = "//getStr([mu, x])) ! fpp
#else
#error  "Unrecognized interface."
#endif

#if     XDI_ENABLED || XMI_ENABLED
        CHECK_ASSERTION(__LINE__, 0 < invSigma, SK_"@setExpCDF(): The condition `0 < invSigma` must hold. invSigma = "//getStr(invSigma)) ! fpp
#elif   !XDD_ENABLED
#error  "Unrecognized interface."
#endif

        ! Compute the CDF.

#if     XDD_ENABLED
        cdf = 1._TKC - exp(-x)
#elif   XDI_ENABLED
        cdf = 1._TKC - exp(-x * invSigma)
#elif   XMI_ENABLED
        cdf = 1._TKC - exp(-(x - mu) * invSigma)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%
#elif   getExpRand_ENABLED
        !%%%%%%%%%%%%%%%%%

        call setUnifRand(rand)
        rand = 1._TKC - rand
        call setExpRand(rand, sigma)
        if (present(mu)) rand = rand + mu

        !%%%%%%%%%%%%%%%%%
#elif   setExpRand_ENABLED
        !%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0 < rand, SK_"@setExpRand(): The condition `0 < rand` must hold. rand = "//getStr(rand)) ! fpp
        CHECK_ASSERTION(__LINE__, rand <= 1, SK_"@setExpRand(): The condition `rand <= 1` must hold. rand = "//getStr(rand)) ! fpp
#if     SD_ENABLED || SM_ENABLED
        CHECK_ASSERTION(__LINE__, 0 < sigma, SK_"@setExpRand(): The condition `0 < sigma` must hold. sigma = "//getStr(sigma)) ! fpp
#endif
#if     DD_ENABLED
        rand = -log(rand)
#elif   SD_ENABLED
        rand = -log(rand) * sigma
#elif   SM_ENABLED
        rand = -log(rand) * sigma + mu
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif