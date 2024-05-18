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
!>  This include file contains implementations of the procedures in module [pm_mathGamma](@ref pm_mathGamma).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Sunday 11:23 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%
#if     getGammaInc_ENABLED
        !%%%%%%%%%%%%%%%%%%

        integer(IK) :: info
#if     Low_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = SK_"@getGammaIncLow()"
        call setGammaIncLow(gammaIncLow, x, log_gamma(kappa), kappa, info, tol)
#elif   Upp_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = SK_"@getGammaIncUpp()"
        call setGammaIncUpp(gammaIncUpp, x, log_gamma(kappa), kappa, info, tol)
#else
#error  "Unrecognized interface."
#endif
        if (info < 0_IK) error stop SK_"@file::"//__FILE__//SK_"@line::"//getStr(__LINE__)//&! LCOV_EXCL_LINE
        MODULE_NAME//PROCEDURE_NAME//SK_": The Incomplete Gamma function failed to converge." ! LCOV_EXCL_LINE

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGammaInc_ENABLED && Def_ENABLED && Low_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (x < kappa + 1._RKC) then
            call setGammaIncLowSeries(gammaIncLow, x, logGammaKappa, kappa, info, tol)
        else
            call setGammaIncUppContFrac(gammaIncLow, x, logGammaKappa, kappa, info, tol)
            gammaIncLow = 1._RKC - gammaIncLow
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGammaInc_ENABLED && Def_ENABLED && Upp_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (x < kappa + 1._RKC) then
            call setGammaIncLowSeries(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
            gammaIncUpp = 1._RKC - gammaIncUpp
        else
            call setGammaIncUppContFrac(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (Series_ENABLED && Low_ENABLED) || (ContFrac_ENABLED && Upp_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK)     , parameter :: ITMAX = 300_IK               !< The maximum allowed number of iterations.
        real(RKC)       , parameter :: EPS = 10 * epsilon(x)        !< The default relative accuracy.
        real(RKC)                   :: tol_def, delta
#if     Series_ENABLED && Low_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = SK_"@setGammaIncLowSeries()"
        real(RKC)                   :: kappaIncremented, summ
#elif   ContFrac_ENABLED && Upp_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = SK_"@setGammaIncUppContFrac()"
        real(RKC)       , parameter :: FPMIN_DEF = tiny(x) / EPS    !< A number near the smallest representable floating-point number.
        real(RKC)                   :: an, b, c, d, h
        real(RKC)                   :: fpmin
#else
#error  "Unrecognized interface."
#endif
        ! Compute the Gamma function using the fastest method decided at runtime.
        CHECK_ASSERTION(__LINE__, 0._RKC <= x, PROCEDURE_NAME//SK_": The input `x` must be positive. x = "//getStr(x)) ! fpp
        CHECK_ASSERTION(__LINE__, 0._RKC < getOption(EPS, tol), PROCEDURE_NAME//SK_": The condition `0. < tol` must hold. tol = "//getStr(getOption(EPS, tol))) ! fpp
        CHECK_ASSERTION(__LINE__, kappa > 0._RKC, PROCEDURE_NAME//SK_": The input `kappa` must be positive. kappa = "//getStr(kappa)) ! fpp
        CHECK_ASSERTION(__LINE__, abs(log_gamma(kappa) - logGammaKappa) <= 100 * epsilon(kappa), \
        PROCEDURE_NAME//SK_": The input `logGammaKappa` must equal `log_gamma(kappa)`. kappa, log_gamma(kappa), logGammaKappa = "//\
        getStr([kappa, log_gamma(kappa), logGammaKappa])) ! fpp
#if     Series_ENABLED && Low_ENABLED
        if (x > 0._RKC) then
            if (present(tol)) then
                tol_def = tol
            else
                tol_def = EPS
            end if
            kappaIncremented = kappa
            summ = 1._RKC / kappa
            delta = summ
            do info = 1, ITMAX
                kappaIncremented = kappaIncremented + 1._RKC
                delta = delta * x / kappaIncremented
                summ = summ + delta
                if (abs(summ) * tol_def < abs(delta)) cycle
                gammaIncLow = summ * exp(kappa * log(x) - x - logGammaKappa)
                return
            end do
            info = -info ! LCOV_EXCL_LINE
        else
            info = 0_IK
            gammaIncLow = 0._RKC
        end if
#elif   ContFrac_ENABLED && Upp_ENABLED
        if (x > 0._RKC) then
            if (present(tol)) then
                tol_def = tol
                fpmin = tiny(x) / tol_def
            else
                tol_def = EPS
                fpmin = FPMIN_DEF
            end if
            b = x + 1._RKC - kappa
            c = 1._RKC / fpmin
            d = 1._RKC / b
            h = d
            do info = 1, ITMAX
                an = -info * (info - kappa)
                b = b + 2._RKC
                d = an * d + b
                if (abs(d) < fpmin) d = fpmin
                c = b + an / c
                if (abs(c) < fpmin) c = fpmin
                d = 1._RKC / d
                delta = d * c
                h = h * delta
                if (tol_def < abs(delta - 1._RKC)) cycle
                gammaIncUpp = exp(kappa * log(x) - x - logGammaKappa) * h
                return
            end do
            info = -info ! LCOV_EXCL_LINE
        else
            info = 0_IK
            gammaIncUpp = 1._RKC
        end if
#endif
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif