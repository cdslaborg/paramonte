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
!>  This include file contains the implementation of procedures in [pm_mathBeta](@ref pm_mathBeta).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%
#if     getLogBeta_ENABLED
        !%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0._RKG < beta, SK_"@getLogBeta(): The condition `0._RKG < beta` must hold. beta = "//getStr(beta)) ! fpp
        CHECK_ASSERTION(__LINE__, 0._RKG < alpha, SK_"@getLogBeta(): The condition `0._RKG < alpha` must hold. alpha = "//getStr(alpha)) ! fpp
        logFuncBeta = log_gamma(alpha) + log_gamma(beta) - log_gamma(alpha + beta)

        !%%%%%%%%%%%%%%%%%
#elif   getBetaInc_ENABLED
        !%%%%%%%%%%%%%%%%%

        real(RKG)   , parameter :: reltol = epsilon(0._RKG)**.7, THRESH = 3000._RKG ! Based on the suggestion of Press et al (1992). This may need fine-tuning for modern hardware.
        logical(LK) :: dengis
        integer(IK) :: info
        if (present(signed)) then
            dengis = signed
        else
            dengis = .false._LK
        end if
        if (alpha < THRESH .or. beta < THRESH) then
            call setBetaInc(betaInc, x, alpha, beta, getLogBeta(alpha, beta), dengis, info)
        else
            info = 0_IK
        end if
        if (info /= 0_IK) call setBetaInc(betaInc, x, alpha, beta, getLogBeta(alpha, beta), reltol, dengis, info)
        !if (betaInc < 0._RKG) betaInc = betaInc + 1._RKG
        if (info /= 0_IK) error stop MODULE_NAME//SK_"@getBetaInc(): Failed to converge in computing the continued fraction representation of the Beta Function." ! LCOV_EXCL_LINE

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setBetaInc_ENABLED && GK21_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG), parameter :: EPS = epsilon(0._RKG), lb = 0._RKG
        integer(IK), parameter :: nintmax = 100
        logical(LK) :: mirrored
        integer(IK) :: sindex(nintmax), neval, nint
        real(RKG) :: abserr, sinfo(4, nintmax), betaMinus1, alphaMinus1
        CHECK_ASSERTION(__LINE__, 0._RKG < beta, SK_"@setBetaInc(): The condition `0. < beta` must hold. beta = "//getStr(beta)) ! fpp
        CHECK_ASSERTION(__LINE__, 0._RKG < alpha, SK_"@setBetaInc(): The condition `0. < alpha` must hold. alpha = "//getStr(alpha)) ! fpp
        mirrored = signed .and. x < 0._RKG
        if (mirrored) x = -x
        if(0._RKG < x .and. x < 1._RKG) then
            if (.not. mirrored) then
                mirrored = (alpha + 1._RKG) / (alpha + beta + 2._RKG) < x
                if (mirrored) x = 1._RKG - x
            endif
            if (mirrored) then ! swap
                abserr = beta
                beta = alpha
                alpha = abserr
            endif
            ! Prepare the integrand parameters.
            betaMinus1 = beta - 1._RKG
            alphaMinus1 = alpha - 1._RKG
            ! integrate.
            info =getQuadErr( getFunc & ! LCOV_EXCL_LINE
                            , lb = lb & ! LCOV_EXCL_LINE
                            , ub = x & ! LCOV_EXCL_LINE
                            , abstol = EPS & ! LCOV_EXCL_LINE
                            , reltol = reltol & ! LCOV_EXCL_LINE
                            , qrule = qrule & ! LCOV_EXCL_LINE
                            , integral = betaInc & ! LCOV_EXCL_LINE
                            , abserr = abserr & ! LCOV_EXCL_LINE
                            , sinfo = sinfo & ! LCOV_EXCL_LINE
                            , sindex = sindex & ! LCOV_EXCL_LINE
                            , neval = neval & ! LCOV_EXCL_LINE
                            , nint = nint & ! LCOV_EXCL_LINE
                            , help = weps & ! LCOV_EXCL_LINE
                            )
            if (mirrored) then
                if (signed) then
                    betaInc = -betaInc
                else
                    betaInc = 1._RKG - betaInc
                end if
            end if
        elseif (0._RKG == x .or. x == 1._RKG) then
            info = 0_IK
            betaInc = x
        else
            info = 1_IK
            CHECK_ASSERTION(__LINE__, info == 0_IK, SK_"@setBetaInc(): The condition `0. <= x .and. x <= 1.` must hold. x = "//getStr(x)) ! fpp
        endif
    contains
        pure function getFunc(xx) result(func)
            real(RKG), intent(in) :: xx
            real(RKG) :: func
            func = exp(log(xx) * alphaMinus1 + log(1._RKG - xx) * betaMinus1 - logFuncBeta)
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setBetaInc_ENABLED && Def_ENABLED && 1
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! This implementation is partially inspired by the proposed approach of Numerical Recipes (1992) by the Great Bill Press et al.
        logical(LK)                 :: mirrored
        integer(IK)                 :: iter, iterTwice
        integer(IK) , parameter     :: MAX_ITER = 10000_IK
        real(RKG)   , parameter     :: EPS = epsilon(0._RKG), FPMIN = tiny(0._RKG) / EPS
        REAL(RKG)                   :: aa, c, d, delta, sumAlphaBeta, alphamMinusOne, alphamPlusOne
        CHECK_ASSERTION(__LINE__, 0._RKG < beta, SK_"@setBetaInc(): The condition `0._RKG < beta` must hold. beta = "//getStr(beta)) ! fpp
        CHECK_ASSERTION(__LINE__, 0._RKG < alpha, SK_"@setBetaInc(): The condition `0._RKG < alpha` must hold. alpha = "//getStr(alpha)) ! fpp
        mirrored = signed .and. x < 0._RKG
        if (mirrored) x = -x
        if(0._RKG < x .and. x < 1._RKG) then
            sumAlphaBeta = alpha + beta
            if (.not. mirrored) then
                mirrored = (alpha + 1._RKG) / (sumAlphaBeta + 2._RKG) < x
                if (mirrored) x = 1._RKG - x
            endif
            if (mirrored) then
                d = beta
                beta = alpha
                alpha = d
            endif
            alphamPlusOne = alpha + 1._RKG
            alphamMinusOne = alpha - 1._RKG
            d = 1._RKG - sumAlphaBeta * x / alphamPlusOne
            if (abs(d) < FPMIN) d = FPMIN
            c = 1._RKG
            d = 1._RKG / d
            betaInc = d
            do iter = 1, MAX_ITER
                iterTwice = 2_IK * iter
                aa = iter * (beta - iter) * x / ((alphamMinusOne + iterTwice) * (alpha + iterTwice))
                d = 1._RKG + aa * d
                if (abs(d) < FPMIN) d = FPMIN
                c = 1._RKG + aa / c
                if (abs(c) < FPMIN) c = FPMIN
                d = 1._RKG / d
                betaInc = betaInc * d * c
                aa = -(alpha + iter) * (sumAlphaBeta + iter) * x / ((alpha + iterTwice) * (alphamPlusOne + iterTwice))
                d = 1._RKG + aa * d
                if (abs(d) < FPMIN) d = FPMIN
                c = 1._RKG + aa / c
                if (abs(c) < FPMIN) c = FPMIN
                d = 1._RKG / d
                delta = d * c
                betaInc = betaInc * delta
                if (EPS < abs(delta - 1._RKG)) cycle
               !betaInc = betaInc * exp(alpha * log(x) + beta * getLog1p(-x) - logFuncBeta) / alpha
                betaInc = betaInc * exp(alpha * log(x) + beta * log(1._RKG - x) - logFuncBeta) / alpha
                !if (mirrored) betaInc = 1._RKG - betaInc ! 1._RKG is the regularization term.
                !if (mirrored) betaInc = -betaInc ! 1._RKG is the regularization term.
                if (mirrored) then
                    if (signed) then
                        betaInc = -betaInc
                    else
                        betaInc = 1._RKG - betaInc
                    end if
                end if
                info = 0_IK
                return
            end do
            info = 1_IK
        elseif (0._RKG == x .or. x == 1._RKG) then
            info = 0_IK
            betaInc = x
        else
            info = 1_IK
            CHECK_ASSERTION(__LINE__, info == 0_IK, SK_"@setBetaInc(): The condition `0. <= x .and. x <= 1.` must hold. x = "//getStr(x)) ! fpp
        endif

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   setBetaInc_ENABLED && 1
        !%%%%%%%%%%%%%%%%%%%%%%

        ! This implementation is based on netlib BETAIN.
        ! As far as evidenced by tests, it requires many cycles for convergence for some precisions other than 64-bit.
        ! Its performance and accuracy against the alternative above is yet to be tested.
        ! KL Majumder, GP Bhattacharjee. Modifications by John Burkardt.
        ! Reference:
        ! KL Majumder, GP Bhattacharjee,
        ! Algorithm AS 63:
        ! The incomplete Beta Integral,
        ! Applied Statistics,
        ! Volume 22, Number 3, 1973, pages 409-411.
        integer(IK) :: ns, iter
        logical(LK) :: mirrored
        real(RKG)   :: ai, sumAlphaBeta, rx, temp, term
        real(RKG)   , parameter :: EPS = epsilon(0._RKG)**(7._RKG/8._RKG)
        integer(IK) , parameter :: MAX_ITER = 10000_IK
        CHECK_ASSERTION(__LINE__, 0._RKG < beta, SK_"@setBetaInc(): The condition `0._RKG < beta` must hold. beta = "//getStr(beta)) ! fpp
        CHECK_ASSERTION(__LINE__, 0._RKG < alpha, SK_"@setBetaInc(): The condition `0._RKG < alpha` must hold. alpha = "//getStr(alpha)) ! fpp
        mirrored = signed .and. x < 0._RKG
        if (mirrored) x = -x
        if (0._RKG < x .and. x < 1._RKG) then

            ! Change tail if necessary.

            sumAlphaBeta = alpha + beta
            if (.not. mirrored) then
                !mirrored = (alpha + 1._RKG) / (sumAlphaBeta + 2._RKG) < x ! this is the original netlib condition.
                mirrored = alpha < sumAlphaBeta * x ! This is the implementation by press et al.
                if (mirrored) x = 1._RKG - x
            endif
            if (mirrored) then
                temp = beta
                beta = alpha
                alpha = temp
            endif

            ai = 1._RKG
            term = 1._RKG
            betaInc = 1._RKG
            ns = int(beta + (1._RKG - x) * sumAlphaBeta, IK)
            ! Use the Soper reduction formula.
            rx = x / (1._RKG - x)
            temp = beta - ai
            if (ns == 0_IK) rx = x
            do iter = 1, MAX_ITER
                term = term * temp * rx / (alpha + ai)
                betaInc = betaInc + term
                temp = abs(term)
                if (temp <= EPS .and. temp <= EPS * betaInc) then
                   !betaInc = betaInc * exp(alpha * log(x) + (beta - 1._RKG) * getLog1p(-x) - logFuncBeta) / alpha
                    betaInc = betaInc * exp(alpha * log(x) + (beta - 1._RKG) * log(1._RKG - x) - logFuncBeta) / alpha
                    if (mirrored) then
                        if (signed) then
                            betaInc = -betaInc
                        else
                            betaInc = 1._RKG - betaInc
                        end if
                    end if
                    info = 0_IK
                    return
                end if
                ai = ai + 1._RKG
                ns = ns - 1
                if (ns < 0_IK) then
                    temp = sumAlphaBeta
                    sumAlphaBeta = sumAlphaBeta + 1._RKG
                else
                    temp = beta - ai
                    if (ns == 0_IK) rx = x
                end if
            end do
            info = 1_IK
        elseif (0._RKG == x .or. x == 1._RKG) then
            info = 0_IK
            betaInc = x
        else
            info = 1_IK
            CHECK_ASSERTION(__LINE__, info == 0_IK, SK_"@setBetaInc(): The condition `0. <= x .and. x <= 1.` must hold. x = "//getStr(x)) ! fpp
        endif

        !%%%%%%%%%%%%%%%%%
#elif   getBetaInv_ENABLED
        !%%%%%%%%%%%%%%%%%

        integer(IK) :: info
        logical(LK) :: dengis
        if (present(signed)) then
            dengis = signed
        else
            dengis = .false._LK
        end if
        call setBetaInv(betaInv, betaInc, alpha, beta, getLogBeta(alpha, beta), dengis, info)
        if (info /= 0_IK) error stop MODULE_NAME//SK_"@getBetaInv(): Failed to converge in computing the Regularized Inverse Incomplete Beta Function." ! LCOV_EXCL_LINE

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   setBetaInv_ENABLED && 1
        !%%%%%%%%%%%%%%%%%%%%%%

        ! This implementation follows the approach of:
        ! GW Cran, KJ Martin, GE Thomas,
        ! Remark AS R19 and Algorithm AS 109:
        ! A Remark on Algorithms AS 63: The Incomplete Beta Integral
        ! and AS 64: Inverse of the Incomplete Beta Integeral,
        ! Applied Statistics,
        ! Volume 26, Number 1, 1977, pages 111-114.
        logical(LK)             :: mirrored
        real(RKG)   , parameter :: underflow = tiny(0._RKG)
        real(RKG)   , parameter :: experr = log10(underflow)
        real(RKG)   , parameter :: ONE_THIRD = 1._RKG / 3._RKG
        real(RKG)   , parameter :: ONE_SIXTH = 1._RKG / 6._RKG
        real(RKG)   , parameter :: FIVE_SIXTH = 5._RKG / 6._RKG
        real(RKG)   , parameter :: ONE_MEPS = 1._RKG - sqrt(epsilon(0._RKG))
        real(RKG)               :: tolerance, adj, g, h, prev, r, s, sq, t, tx, w, betaIncOld, betaIncNew
        CHECK_ASSERTION(__LINE__, 0._RKG < alpha, SK_"@setBetaInv(): The condition `0. < alpha` must hold. alpha = "//getStr(alpha)) ! fpp
        CHECK_ASSERTION(__LINE__, 0._RKG < beta, SK_"@setBetaInv(): The condition `0. < beta` must hold. alpha = "//getStr(beta)) ! fpp
        mirrored = signed .and. betaInc < 0._RKG
        if (mirrored) betaInc = -betaInc
        if (0._RKG < betaInc .and. betaInc < 1._RKG) then

            ! Change tail if necessary.

            if (mirrored) then
                r = alpha
                alpha = beta
                beta = r
            else
                mirrored = 0.5_RKG < betaInc
                if (mirrored) then
                    betaInc = 1._RKG - betaInc
                    r = alpha
                    alpha = beta
                    beta = r
                end if
            end if

            ! Compute the initial approximation.

            r = sqrt(-log(betaInc**2))
            betaIncNew = r - (2.30753_RKG + 0.27061_RKG * r) / (1._RKG + (0.99229_RKG + 0.04481_RKG * r) * r)
            if (1._RKG < alpha .and. 1._RKG < beta) then
                r = (betaIncNew * betaIncNew - 3._RKG) * ONE_SIXTH
                s = 1._RKG / (alpha + alpha - 1._RKG)
                t = 1._RKG / (beta + beta - 1._RKG)
                h = 2._RKG / (s + t)
                w = betaIncNew * sqrt(h + r) / h - (t - s) * (r + FIVE_SIXTH - 2._RKG / (3._RKG * h))
                betaInv = alpha / (alpha + beta * exp(w + w))
            else
                r = beta + beta
                t = 1._RKG / (9._RKG * beta)
                t = r * (1._RKG - t + betaIncNew * sqrt(t))**3
                if (t <= 0._RKG) then
                    betaInv = 1._RKG - exp((log((1._RKG - betaInc) * beta) + logFuncBeta) / beta)
                else
                    t = (4._RKG * alpha + r - 2._RKG) / t
                    if (t <= 1._RKG) then
                        betaInv = exp((log(betaInc * alpha) + logFuncBeta) / alpha)
                    else
                        betaInv = 1._RKG - 2._RKG / (t + 1._RKG)
                    end if
                end if
            end if

            ! Solve for X via modified Newton-Raphson method, using the function `setBetaInc()`.

            t = 1._RKG - beta
            r = 1._RKG - alpha
            betaIncOld = 0._RKG
            prev = 1._RKG
            sq = 1._RKG
            if (betaInv < 0.0001_RKG) betaInv = 0.0001_RKG
            if (0.9999_RKG < betaInv) betaInv = 0.9999_RKG
            tolerance = 10._RKG**int(max(-5._RKG / alpha**2 - 1._RKG / betaInc**0.2 - 13._RKG, experr), IK)

            ! Begin iteration.

            loopFindRoot: do ! 10
                call setBetaInc(betaIncNew, betaInv, alpha, beta, logFuncBeta, signed, info)
                if (info /= 0_IK) return ! LCOV_EXCL_LINE
                if (betaIncNew < 0._RKG) then
                    !if (abs(betaIncNew) < epsilon(0._RKG)) error stop "Underflow detected."
                    betaIncNew = betaIncNew - betaInc
                    !if (abs(betaIncNew) < epsilon(0._RKG)) error stop "Underflow detected."
                    betaIncNew = betaIncNew + 1._RKG
                else
                    betaIncNew = betaIncNew - betaInc
                end if
               !betaIncNew = betaIncNew * exp(logFuncBeta + r * log(betaInv) + t * getLog1p(-betaInv))
                betaIncNew = betaIncNew * exp(logFuncBeta + r * log(betaInv) + t * log(1._RKG - betaInv))
                if (betaIncNew * betaIncOld <= 0._RKG) prev = max(sq, underflow)
                g = 1._RKG
                loopRefine: do ! 20
                    adj = g * betaIncNew
                    sq = adj * adj
                    if (sq < prev) then
                        tx = betaInv - adj
                        if (tx < 0._RKG .or. 1._RKG < tx) then
                            g = g * ONE_THIRD
                            cycle loopRefine
                        end if
                    else ! 30
                        g = g * ONE_THIRD
                        cycle loopRefine
                    end if
                    ! Check whether current estimate is acceptable. The change "VALUE = TX" was suggested by Ivan Ukhov.
                    if (prev <= tolerance .or. betaIncNew * betaIncNew <= tolerance) then ! 40
                        betaInv = tx
                        exit loopFindRoot
                    end if
                    if (tx /= 0._RKG .and. tx /= 1._RKG) exit loopRefine
                    g = g * ONE_THIRD
                end do loopRefine
                if (betaInv == tx) exit loopFindRoot
                betaInv = tx
                betaIncOld = betaIncNew
            end do loopFindRoot
            if (mirrored) then
                if (signed) then
                    betaInv = -betaInv
                else
                    betaInv = 1._RKG - betaInv
                end if
            end if

        elseif (0._RKG == betaInc .or. betaInc == 1._RKG) then

            info = 0_IK
            betaInv = betaInc

        else

            info = 1_IK
            CHECK_ASSERTION(__LINE__, info == 0_IK, SK_"@setBetaInc(): The condition `0. <= betaInc .and. betaInc <= 1.` must hold. betaInc = "//getStr(betaInc)) ! fpp

        end if

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   setBetaInv_ENABLED && 0
        !%%%%%%%%%%%%%%%%%%%%%%

        ! This implementation follows the proposed approach of Numerical Recipes by Press et al. (2007) using Halley approximation.
        ! No sign corrections are implemented in this algorithm. As such the output is prone to numerical roundup to `1`.
        integer(IK) :: iter
        real(RKG)   , parameter :: EPS = sqrt(epsilon(0._RKG))
        real(RKG)   :: pp, t, u, betaIncNew, al, h, w, alphaMinus1, betaMinus1
        if (signed .and. betaInc < 0._RKG) betaInc = 1 + betaInc
        if (0._RKG < betaInc .and. betaInc < 1._RKG) then
            betaMinus1 = beta - 1._RKG
            alphaMinus1 = alpha - 1._RKG
            if (alpha < 1._RKG .or. beta < 1._RKG) then
                u = exp(beta * log(beta / (alpha + beta))) / beta
                t = exp(alpha * log(alpha / (alpha + beta))) / alpha
                w = t + u
                if (betaInc < t / w) then
                    betaInv = (alpha * w * betaInc)**(1._RKG / alpha)
                else
                    betaInv = 1._RKG - (beta * w * (1._RKG - betaInc))**(1._RKG / beta)
                end if
            else
                if (betaInc < 0.5_RKG) then
                    pp = betaInc
                else
                    pp = 1._RKG - betaInc
                end if
                t = sqrt(-2._RKG * log(pp))
                betaInv = (2.30753_RKG + t * 0.27061_RKG) / (1._RKG + t * (0.99229_RKG + t * 0.04481_RKG)) - t
                if (betaInc < .5_RKG) betaInv = -betaInv
                al = (betaInv**2 - 3._RKG) / 6._RKG
                h = 2._RKG / (1._RKG / (2._RKG * alpha - 1._RKG) + 1._RKG / (2._RKG * beta - 1._RKG))
                w = (betaInv * sqrt(al + h) / h) - (1._RKG / (2._RKG * beta - 1._RKG) - 1._RKG / (2._RKG * alpha - 1._RKG)) * (al + 5._RKG / 6._RKG - 2._RKG / (3._RKG * h))
                betaInv = alpha / (alpha + beta * exp(2._RKG * w))
            end if
            loopFindRoot: do iter = 1, 10
                if (betaInv == 0._RKG .or. betaInv == 1._RKG) exit loopFindRoot
                call setBetaInc(betaIncNew, betaInv, alpha, beta, logFuncBeta, signed, info)
                if (info /= 0_IK) return
                betaIncNew = betaIncNew - betaInc
                t = exp(alphaMinus1 * log(betaInv) + betaMinus1 * log(1._RKG - betaInv) - logFuncBeta)
                u = betaIncNew / t
                t = u / (1._RKG - .5_RKG * min(1._RKG, u * (alphaMinus1 / betaInv - betaMinus1 / (1._RKG - betaInv))))
                betaInv = betaInv - t
                if (betaInv <= 0._RKG) betaInv = .5_RKG * (betaInv + t)
                if (betaInv >= 1._RKG) betaInv = .5_RKG * (betaInv + t + 1._RKG)
                if (abs(t) < EPS * betaInv .and. 1_IK < iter) exit loopFindRoot
            end do loopFindRoot
            info = 0_IK
        elseif (0._RKG == betaInc .or. betaInc == 1._RKG) then
            info = 0_IK
            betaInv = betaInc
        else
            info = 1_IK
            CHECK_ASSERTION(__LINE__, .false., SK_"@setBetaInv(): The condition `0. <= betaInc .and. betaInc <= 1.` must hold. betaInc = "//getStr(betaInc)) ! fpp
        end if

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   setBetaInv_ENABLED && 0
        !%%%%%%%%%%%%%%%%%%%%%%

        !   Naive root-finding implementation.
        !   This implementation is exploratory and should never be used.
        integer(IK) :: neval
        logical(LK) :: mirrored
        real(RKG)   :: lb, ub, lf, uf, abstol, h, w, r, s, t, betaIncNew
        real(RKG)   , parameter :: ONE_SIXTH = 1._RKG / 6._RKG
        !real(RKG)   , parameter :: abstol = epsilon(0._RKG)
        CHECK_ASSERTION(__LINE__, 0._RKG < alpha, SK_"@setBetaInv(): The condition `0. < alpha` must hold. alpha = "//getStr(alpha)) ! fpp
        CHECK_ASSERTION(__LINE__, 0._RKG < beta, SK_"@setBetaInv(): The condition `0. < beta` must hold. alpha = "//getStr(beta)) ! fpp
        info = 0_IK
        mirrored = signed .and. betaInc < 0._RKG
        if (mirrored) betaInc = -betaInc
        if (0._RKG < betaInc .and. betaInc < 1._RKG) then
            if (.not. mirrored) then
                mirrored = 0.5_RKG < betaInc
                if (mirrored) betaInc = 1._RKG - betaInc
            end if
            if (mirrored) then
                r = alpha
                alpha = beta
                beta = r
            end if
            ! Compute the initial approximation.
            r = sqrt(-log(betaInc**2))
            betaIncNew = r - (2.30753_RKG + 0.27061_RKG * r) / (1._RKG + (0.99229_RKG + 0.04481_RKG * r) * r)
            if (1._RKG < alpha .and. 1._RKG < beta) then
                r = (betaIncNew * betaIncNew - 3._RKG) * ONE_SIXTH
                s = 1._RKG / (alpha + alpha - 1._RKG)
                t = 1._RKG / (beta + beta - 1._RKG)
                h = 2._RKG / (s + t)
                w = betaIncNew * sqrt(h + r) / h - (t - s) * (r + 5._RKG / 6._RKG - 2._RKG / (3._RKG * h))
                betaInv = alpha / (alpha + beta * exp(w + w))
            else
                r = beta + beta
                t = 1._RKG / (9._RKG * beta)
                t = r * (1._RKG - t + betaIncNew * sqrt(t))**3
                if (t <= 0._RKG) then
                    betaInv = 1._RKG - exp((log((1._RKG - betaInc) * beta) + logFuncBeta) / beta)
                else
                    t = (4._RKG * alpha + r - 2._RKG) / t
                    if (t <= 1._RKG) then
                        betaInv = exp((log(betaInc * alpha) + logFuncBeta) / alpha)
                    else
                        betaInv = 1._RKG - 2._RKG / (t + 1._RKG)
                    end if
                end if
            end if
            abstol = epsilon(0._RKG)
            if (betaInv < 0._RKG) betaInv = sqrt(abstol)
            if (1._RKG < betaInv) betaInv = 1._RKG - sqrt(abstol)
            lb = 0._RKG
            ub = 1._RKG
            lf = betaInc
            uf = betaInc - 1._RKG
            if (abstol < ub - lb) then
                call setRoot(newton, getFunc, betaInv, getFuncDiff, lb, ub, lf, uf, abstol, neval)
                if (neval < 0_IK) call setRoot(bisection, getFunc, betaInv, lb, ub, lf, uf, abstol, neval)
                if (neval < 0_IK) info = 1_IK
            else
                if (abs(lf) < abs(uf)) then
                    betaInv = lb
                else
                    betaInv = ub
                end if
            end if
            if (mirrored) then
                if (signed) then
                    betaInv = -betaInv
                else
                    betaInv = 1._RKG - betaInv
                end if
            end if

        elseif (0._RKG == betaInc .or. betaInc == 1._RKG) then
            info = 0_IK
            betaInv = betaInc
        else
            info = 1_IK
            CHECK_ASSERTION(__LINE__, info == 0_IK, SK_"@setBetaInc(): The condition `0. <= betaInc .and. betaInc <= 1.` must hold. betaInc = "//getStr(betaInc)) ! fpp
        end if
    contains
        PURE function getFunc(x) result(func)
            real(RKG)   , intent(in)    :: x ! betaInv guess.
            real(RKG)                   :: func
            integer(IK)                 :: info
            call setBetaInc(func, x, alpha, beta, logFuncBeta, signed, info)
            if (info /= 0_IK) error stop ! LCOV_EXCL_LINE
            func = betaInc - func
        end function
        PURE function getFuncDiff(x, order) result(betaPDF)
            integer(IK) , intent(in)    :: order
            real(RKG)   , intent(in)    :: x
            real(RKG)   , parameter     :: SQRT_HUGE = sqrt(huge(x))
            real(RKG)                   :: betaPDF
            if(order == 0_IK) then
                betaPDF = getFunc(x)
            else
                if (x /= 0._RKG .and. x /= 1._RKG) then
                    call setBetaLogPDF(betaPDF, x, alpha, beta)
                    betaPDF = exp(betaPDF)
                elseif ((x == 0._RKG .and. alpha < 1._RKG) .or. (x == 1._RKG .and. beta < 1._RKG)) then
                    betaPDF = SQRT_HUGE
                else
                    betaPDF = 0._RKG
                end if
            end if
        end function

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif