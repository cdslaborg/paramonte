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
!>  This include file contains implementations of the procedures in module [pm_mathGammaAM](@ref pm_mathGammaAM).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, July 22, 2024, 11:45 AM, NASA Goddard Space Flight Center, Washington, D.C.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getGammaIncLowAM_ENABLED || getGammaIncUppAM_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: info
#if     getGammaIncLowAM_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = SK_"@getGammaIncLowAM()"
        call setGammaIncLowAM(gammaIncLow, x, log_gamma(kappa), kappa, info)
#elif   getGammaIncUppAM_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = SK_"@getGammaIncUppAM()"
        call setGammaIncUppAM(gammaIncUpp, x, log_gamma(kappa), kappa, info)
#else
#error  "Unrecognized interface."
#endif
        if (info < 0_IK) error stop SK_"@file::"//__FILE__//SK_"@line::"//getStr(__LINE__)//&! LCOV_EXCL_LINE
        MODULE_NAME//PROCEDURE_NAME//SK_": The Incomplete Gamma function failed to converge." ! LCOV_EXCL_LINE

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGammaIncLowAM_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        if (x /= 0._RKG) then
            call setGammaIncAM(gammaIncLow, x, kappa, info)
            gammaIncLow = gammaIncLow * exp((kappa * log(x) - x) - logGammaKappa)
            if (kappa < x) gammaIncLow = 1._RKG - gammaIncLow
        else
            gammaIncLow = 0._RKG
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGammaIncUppAM_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        if (x /= 0._RKG) then
            call setGammaIncAM(gammaIncUpp, x, kappa, info)
            gammaIncUpp = gammaIncUpp * exp((kappa * log(x) - x) - logGammaKappa)
            if (x <= kappa) gammaIncUpp = 1._RKG - gammaIncUpp
        else
            gammaIncUpp = 0._RKG
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGammaIncAM_ENABLED && Def_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK), parameter :: ITMAX = 10_IK**6 !< The maximum allowed number of iterations.
        real(RKG), parameter :: EPS = 10 * epsilon(x) ! * 10 !< The default relative accuracy.
        real(RKG), parameter :: FPMIN = tiny(x) / EPS !< A number near the smallest representable floating-point number.
        real(RKG) :: kappa0, numer, denom, frac, coef, delta
        integer(IK) :: infoHalf

        CHECK_ASSERTION(__LINE__, 0 < kappa, SK_"@setGammaIncAM(): The condition `0 <= kappa` must hold. kappa = "//getStr(kappa)) ! fpp

        if(x < -9.0_RKG) then
            kappa0 = 5.0_RKG * sqrt(abs(x)) - 5.0_RKG
        elseif (0._RKG < x) then
            kappa0 = x
        else
            kappa0 = 0._RKG
        endif

        if (kappa0 <= kappa) then

            ! modified Lentz method of continued fraction using Eqn (15) and Algorithm 1 of Abergel, 2020.

            denom = kappa
            numer = 1._RKG
            gammaInc = numer / denom
            coef = 1._RKG / denom ! d
            frac = numer / FPMIN ! c
            do  info = 2, ITMAX
                ! Update the numerator continuant (a).
                infoHalf = info / 2
                if (info == infoHalf * 2_IK) then
                    numer = (1._RKG - kappa - infoHalf) * x
                else
                    numer = infoHalf * x
                end if
                ! Update the denominator continuant (b).
                denom = kappa + (info - 1_IK)
                ! Continue the rest of the Lentz method.
                coef = denom + numer * coef ! d
                if (coef == 0._RKG) coef = FPMIN !if (abs(coef) < FPMIN) coef = sign(FPMIN, coef)
                frac = denom + numer / frac ! c
                if (frac == 0._RKG) frac = FPMIN !if (abs(frac) < FPMIN) frac = sign(FPMIN, frac)
                coef = 1._RKG / coef
                delta = frac * coef
                gammaInc = gammaInc * delta
                if (abs(delta - 1._RKG) < EPS) return
            end do
            info = -info ! LCOV_EXCL_LINE

        elseif (0._RKG <= x) then

            ! modified Lentz method of continued fraction using Eqn (16) and Algorithm 1 of Abergel, 2020.

            numer = 1._RKG
            denom = (x - kappa) + 1._RKG
            gammaInc = numer / denom
            coef = 1._RKG / denom ! d
            frac = numer / FPMIN ! c
            do  info = 2, ITMAX
                ! Update the numerator continuant (a).
                numer = (info - 1_IK) * ((1_IK - info) + kappa)
                ! Update the denominator continuant (b).
                denom = (2_IK * info - 1_IK) + (x - kappa)
                ! Continue the rest of the Lentz method.
                coef = coef * numer + denom ! d
                if (coef == 0._RKG) coef = FPMIN
                frac = denom + numer / frac ! c
                if (frac == 0._RKG) frac = FPMIN
                coef = 1._RKG / coef
                delta = frac * coef
                gammaInc = gammaInc * delta
                if(abs(delta - 1._RKG) < EPS) return
            end do
            info = -info ! LCOV_EXCL_LINE

        else

            call setGammaIncAM(gammaInc, x, log_gamma(kappa), nint(kappa, IK), info)

        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGammaIncAM_ENABLED && NXIK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG), parameter :: EPS = epsilon(x) ! * 10 !< The default relative accuracy.
        real(RKG) :: invx, invxsq, summ, frac, delta
        integer(IK) :: decrement, maxit, kappaHalf
        logical(LK) :: isodd

        CHECK_ASSERTION(__LINE__, x < 0, MODULE_NAME//SK_"@setGammaIncAM(): The condition `x < 0` must hold. x = "//getStr(x)) ! fpp
        CHECK_ASSERTION(__LINE__, 0 < kappa, MODULE_NAME//SK_"@setGammaIncAM(): The condition `0 < kappa` must hold. kappa = "//getStr(kappa)) ! fpp
        CHECK_ASSERTION(__LINE__, max(1_IK, kappa - 1_IK) < abs(x), SK_"@setGammaIncAM(): The condition `max(1, kappa - 1) < abs(x)` must hold. kappa, x = "//getStr([real(RKG) :: kappa, x])) ! fpp

        ! Algorithm 2.
        invx = 1._RKG / x
        invxsq = invx**2
        decrement = kappa - 1_IK ! d
        summ = invx * (decrement + x)
        kappaHalf = kappa / 2_IK
        maxit = kappaHalf - 1_IK
        do  info = 1, maxit
            frac = decrement * (decrement - 1_IK) * invxsq
            decrement = decrement - 2_IK
            delta = -frac * (x + decrement)
            summ = summ + delta
            if (delta < summ * EPS) exit
        end do
        if (maxit < info) info = -info
        isodd = logical(kappaHalf * 2_IK < kappa)
        if (info < 0_IK .and. isodd) summ = summ - decrement * frac * invx
        gammaInc = exp(x + logGammaKappa - (kappa - 1._RKG) * log(-x))
        if (isodd) gammaInc = -gammaInc
        gammaInc = -invx * (gammaInc + summ)

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif