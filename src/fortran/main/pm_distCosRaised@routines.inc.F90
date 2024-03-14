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
!>  This include file contains the implementation of procedures in [pm_distCosRaised](@ref pm_distCosRaised).
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%
#if     getCosRaisedPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        if (present(mu) .and. present(sigma)) then
            call setCosRaisedPDF(pdf, x, mu, 1._TKC / sigma)
        elseif (present(sigma)) then
            call setCosRaisedPDF(pdf, x, 0._TKC, 1._TKC / sigma)
        elseif (present(mu)) then
            call setCosRaisedPDF(pdf, x, mu)
        else
            call setCosRaisedPDF(pdf, x)
        end if

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   setCosRaisedPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        real(TKC), parameter :: PI = acos(-1._TKC)
#if     XMD_ENABLED || XDD_ENABLED
        real(TKC), parameter :: invSigma = 1._TKC
#if     XDD_ENABLED
        real(TKC), parameter :: mu = 0._TKC
#endif
#elif   XMI_ENABLED
        CHECK_ASSERTION(__LINE__, 0._TKC < invSigma, SK_"@setCosRaisedPDF(): The condition `0. < invSigma` must hold. invSigma = "//getStr(invSigma))
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, x <= mu + 1._TKC / invSigma, SK_"@getCosRaisedPDF(): The condition `x <= mu + 1. / invSigma` must hold. x, mu, invSigma = "//getStr([x, mu, invSigma]))
        CHECK_ASSERTION(__LINE__, x >= mu - 1._TKC / invSigma, SK_"@getCosRaisedPDF(): The condition `x >= mu - 1. / invSigma` must hold. x, mu, invSigma = "//getStr([x, mu, invSigma]))
        pdf = 0.5_TKC * invSigma * (1._TKC + cos(PI * invSigma * (x - mu)))

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   getCosRaisedCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        if (present(mu) .and. present(sigma)) then
            call setCosRaisedCDF(cdf, x, mu, 1._TKC / sigma)
        elseif (present(sigma)) then
            call setCosRaisedCDF(cdf, x, 0._TKC, 1._TKC / sigma)
        elseif (present(mu)) then
            call setCosRaisedCDF(cdf, x, mu)
        else
            call setCosRaisedCDF(cdf, x)
        end if

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   setCosRaisedCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        real(TKC) :: stanx
        real(TKC), parameter :: PI = acos(-1._TKC)
        real(TKC), parameter :: INV_PI = 1._TKC / PI
#if     XMD_ENABLED || XDD_ENABLED
        real(TKC), parameter :: invSigma = 1._TKC
#if     XDD_ENABLED
        real(TKC), parameter :: mu = 0._TKC
#endif
#elif   XMI_ENABLED
        CHECK_ASSERTION(__LINE__, 0._TKC < invSigma, SK_"@setCosRaisedCDF(): The condition `0. < invSigma` must hold. invSigma = "//getStr(invSigma))
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, x <= mu + 1._TKC / invSigma, SK_"@getCosRaisedCDF(): The condition `x <= mu + 1. / invSigma` must hold. x, mu, invSigma = "//getStr([x, mu, invSigma]))
        CHECK_ASSERTION(__LINE__, x >= mu - 1._TKC / invSigma, SK_"@getCosRaisedCDF(): The condition `x >= mu - 1. / invSigma` must hold. x, mu, invSigma = "//getStr([x, mu, invSigma]))
        stanx = (x - mu) * invSigma
        cdf = 0.5_TKC * (1._TKC + stanx + sin(stanx * PI) * INV_PI)
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif