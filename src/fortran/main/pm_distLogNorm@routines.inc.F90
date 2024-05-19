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
!>  This include file contains procedure implementations of [pm_distLogNorm](@ref pm_distLogNorm).
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%
#if     getLogNormLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: mu_def, invSigma
        CHECK_ASSERTION(__LINE__, x > 0._RKG, SK_"@getLogNormLogPDF(): The condition `x > 0._RKG` must hold. x = "//getStr(x))
        if (present(mu)) then
            mu_def = mu
        else
            mu_def = 0._RKG
        end if
        if (present(sigma)) then
            invSigma = 1._RKG / sigma
        else
            invSigma = 1._RKG
        end if
        call setLogNormLogPDF(logPDF, log(x), mu_def, invSigma, log(invSigma))

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   setLogNormLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG)   , parameter :: PI = acos(-1._RKG)
        real(RKG)   , parameter :: INVERSE_SQRT_TWO_PI = 1._RKG / sqrt(2._RKG * PI)
        real(RKG)   , parameter :: LOG_INVERSE_SQRT_TWO_PI = log(INVERSE_SQRT_TWO_PI)
#if     DS_ENABLED || MS_ENABLED
        CHECK_ASSERTION(__LINE__, invSigma > 0._RKG, SK_"@setLogNormLogPDF(): The condition `invSigma > 0._RKG` must hold. invSigma = "//getStr(invSigma))
        CHECK_ASSERTION(__LINE__, abs(logInvSigma - log(invSigma)) <= epsilon(0._RKG)*100, SK_"@setLogNormLogPDF(): The condition `abs(logInvSigma - log(invSigma)) <= epsilon(0._RKG)*100` must hold. logInvSigma, log(invSigma) = "//getStr([logInvSigma, log(invSigma)]))
#endif
#if     DD_ENABLED
        logPDF = LOG_INVERSE_SQRT_TWO_PI - 0.5_RKG * logx**2 - logx
#elif   DS_ENABLED
        logPDF = logInvSigma + LOG_INVERSE_SQRT_TWO_PI - 0.5_RKG * (logx * invSigma)**2 - logx
#elif   MD_ENABLED
        logPDF = LOG_INVERSE_SQRT_TWO_PI - 0.5_RKG * (logx - mu)**2 - logx
#elif   MS_ENABLED
        logPDF = logInvSigma + LOG_INVERSE_SQRT_TWO_PI - 0.5_RKG * ((logx - mu) * invSigma)**2 - logx
#else
#error  "Unrecognized interface"
#endif

        !%%%%%%%%%%%%%%%%%%%%
#elif   getLogNormCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, x >= 0._RKG, SK_"@getLogNormCDF(): The condition `x >= 0._RKG` must hold. x = "//getStr(x))
        if (x > 0._RKG) then
            if (present(mu) .and. present(sigma)) then
                call setLogNormCDF(cdf, log(x), mu, 1._RKG / sigma)
            elseif (present(sigma)) then
                call setLogNormCDF(cdf, log(x), 0._RKG, 1._RKG / sigma)
            elseif (present(mu)) then
                call setLogNormCDF(cdf, log(x), mu)
            else
                call setLogNormCDF(cdf, log(x))
            end if
        else
            cdf = 0._RKG
        end if

        !%%%%%%%%%%%%%%%%%%%%
#elif   setLogNormCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        real(RKG)   , parameter :: INVERSE_SQRT_TWO = 1._RKG / sqrt(2._RKG)
#if     DD_ENABLED
        cdf = 0.5_RKG * (1._RKG + erf(logx * INVERSE_SQRT_TWO))
#elif   MD_ENABLED
        cdf = 0.5_RKG * (1._RKG + erf((logx - mu) * INVERSE_SQRT_TWO))
#elif   MS_ENABLED
        CHECK_ASSERTION(__LINE__, invSigma > 0._RKG, SK_"@setLogNormCDF(): The condition `invSigma > 0._RKG` must hold. invSigma = "//getStr(invSigma))
        cdf = 0.5_RKG * (1._RKG + erf((logx - mu) * invSigma * INVERSE_SQRT_TWO))
#else
#error "Unrecognized interface"
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface"
        !%%%%%%%%%%%%%%%%%%%%%%%
#endif