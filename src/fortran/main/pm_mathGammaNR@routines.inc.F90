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
!>  This include file contains implementations of the procedures in module [pm_mathGammaNR](@ref pm_mathGammaNR).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Sunday 11:23 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getGammaIncLowNR_ENABLED || getGammaIncUppNR_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: info
#if     getGammaIncLowNR_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = SK_"@getGammaIncLowNR()"
        call setGammaIncLowNR(gammaIncLow, x, log_gamma(kappa), kappa, info, tol)
#elif   getGammaIncUppNR_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = SK_"@getGammaIncUppNR()"
        call setGammaIncUppNR(gammaIncUpp, x, log_gamma(kappa), kappa, info, tol)
#else
#error  "Unrecognized interface."
#endif
        if (info < 0_IK) error stop SK_"@file::"//__FILE__//SK_"@line::"//getStr(__LINE__)//&! LCOV_EXCL_LINE
        MODULE_NAME//PROCEDURE_NAME//SK_": The Incomplete Gamma function failed to converge." ! LCOV_EXCL_LINE

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGammaIncLowNR_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        if (x < kappa + 1._RKG) then
            call setGammaIncLowSeriesNR(gammaIncLow, x, logGammaKappa, kappa, info, tol)
        else
            call setGammaIncUppContFracNR(gammaIncLow, x, logGammaKappa, kappa, info, tol)
            gammaIncLow = 1._RKG - gammaIncLow
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGammaIncUppNR_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        if (x < kappa + 1._RKG) then
            call setGammaIncLowSeriesNR(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
            gammaIncUpp = 1._RKG - gammaIncUpp
        else
            call setGammaIncUppContFracNR(gammaIncUpp, x, logGammaKappa, kappa, info, tol)
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGammaIncLowSeriesNR_ENABLED || setGammaIncUppContFracNR_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG)                   :: tol_def, delta
        real(RKG)       , parameter :: EPS = 10 * epsilon(x) !< The default relative accuracy.
        integer(IK)     , parameter :: ITMAX = 10**6_IK !< The maximum allowed number of iterations.
#if     setGammaIncLowSeriesNR_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = SK_"@setGammaIncLowSeriesNR()"
        real(RKG)                   :: summ!, kappaIncremented
#elif   setGammaIncUppContFracNR_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = SK_"@setGammaIncUppContFracNR()"
        real(RKG)       , parameter :: FPMIN_DEF = tiny(x) / EPS    !< A number near the smallest representable floating-point number.
        real(RKG)                   :: an, b, c, d, h
        real(RKG)                   :: fpmin
#else
#error  "Unrecognized interface."
#endif
        ! Compute the Gamma function using the fastest method decided at runtime.
        CHECK_ASSERTION(__LINE__, 0._RKG <= x, PROCEDURE_NAME//SK_": The input `x` must be positive. x = "//getStr(x)) ! fpp
        CHECK_ASSERTION(__LINE__, 0._RKG < getOption(EPS, tol), PROCEDURE_NAME//SK_": The condition `0. < tol` must hold. tol = "//getStr(getOption(EPS, tol))) ! fpp
        CHECK_ASSERTION(__LINE__, kappa > 0._RKG, PROCEDURE_NAME//SK_": The input `kappa` must be positive. kappa = "//getStr(kappa)) ! fpp
        CHECK_ASSERTION(__LINE__, abs(log_gamma(kappa) - logGammaKappa) <= 100 * epsilon(kappa), \
        PROCEDURE_NAME//SK_": The input `logGammaKappa` must equal `log_gamma(kappa)`. kappa, log_gamma(kappa), logGammaKappa = "//\
        getStr([kappa, log_gamma(kappa), logGammaKappa])) ! fpp
#if     setGammaIncLowSeriesNR_ENABLED && 0
        if (x > 0._RKG) then
            if (present(tol)) then
                tol_def = tol
            else
                tol_def = EPS
            end if
            !kappaIncremented = kappa
            summ = 1._RKG / kappa
            delta = summ
            do info = 1, ITMAX
                !kappaIncremented = kappaIncremented + 1._RKG
                delta = delta * x / (kappa + info)
                summ = summ + delta
                if (abs(summ) * tol_def < abs(delta)) cycle
                gammaIncLow = summ * exp(kappa * log(x) - x - logGammaKappa)
                return
            end do
            info = -info ! LCOV_EXCL_LINE
        else
            info = 0_IK
            gammaIncLow = 0._RKG
        end if
#elif   setGammaIncLowSeriesNR_ENABLED && 1
        ! Use compensated summation to improve summation accuracy for large x ~ kappa cases.
        block
            real(RKG) :: err, temp, newval
            if (x > 0._RKG) then
                if (present(tol)) then
                    tol_def = tol
                else
                    tol_def = EPS
                end if
                !kappaIncremented = kappa
                summ = 1._RKG / kappa
                delta = summ
                err = 0._RKG
                do info = 1, ITMAX
                    !kappaIncremented = kappaIncremented + 1._RKG
                    !delta = delta * x / kappaIncremented
                    !!invdelta = (kappa + info) / (x * delta)
                    delta = delta * (x / (kappa + info))
                    temp = summ
                    newval = delta + err
                    summ = temp + newval
                    err = (temp - summ) + newval
                    !summ = summ + delta
                    if (abs(summ) * tol_def < abs(delta)) cycle
                    gammaIncLow = summ * exp(kappa * log(x) - x - logGammaKappa)
                    return
                end do
                info = -info ! LCOV_EXCL_LINE
            else
                info = 0_IK
                gammaIncLow = 0._RKG
            end if
        end block
#elif   setGammaIncUppContFracNR_ENABLED
        if (x > 0._RKG) then
            if (present(tol)) then
                tol_def = tol
                fpmin = tiny(x) / tol_def
            else
                tol_def = EPS
                fpmin = FPMIN_DEF
            end if
            b = x + 1._RKG - kappa
            c = 1._RKG / fpmin
            d = 1._RKG / b
            h = d
            do info = 1, ITMAX
                an = -info * (info - kappa)
                b = b + 2._RKG
                d = an * d + b
                if (abs(d) < fpmin) d = fpmin
                c = b + an / c
                if (abs(c) < fpmin) c = fpmin
                d = 1._RKG / d
                delta = d * c
                h = h * delta
                if (tol_def < abs(delta - 1._RKG)) cycle
                gammaIncUpp = exp(kappa * log(x) - x - logGammaKappa) * h
                return
            end do
            info = -info ! LCOV_EXCL_LINE
        else
            info = 0_IK
            gammaIncUpp = 1._RKG
        end if
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif