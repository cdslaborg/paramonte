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

        CHECK_ASSERTION(__LINE__, 0._RKC < beta, SK_"@getLogBeta(): The condition `0._RKC < beta` must hold. beta = "//getStr(beta)) ! fpp
        CHECK_ASSERTION(__LINE__, 0._RKC < alpha, SK_"@getLogBeta(): The condition `0._RKC < alpha` must hold. alpha = "//getStr(alpha)) ! fpp
        logFuncBeta = log_gamma(alpha) + log_gamma(beta) - log_gamma(alpha + beta)

        !%%%%%%%%%%%%%%%%%
#elif   getBetaInc_ENABLED
        !%%%%%%%%%%%%%%%%%

        real(RKC)   , parameter :: reltol = epsilon(0._RKC)**.7, THRESH = 3000._RKC ! Based on the suggestion of Press et al (1992). This may need fine-tuning for modern hardware.
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
        !if (betaInc < 0._RKC) betaInc = betaInc + 1._RKC
        if (info /= 0_IK) error stop MODULE_NAME//SK_"@getBetaInc(): Failed to converge in computing the continued fraction representation of the Beta Function." ! LCOV_EXCL_LINE

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setBetaInc_ENABLED && GK21_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKC), parameter :: EPS = epsilon(0._RKC), lb = 0._RKC
        integer(IK), parameter :: nintmax = 100
        logical(LK) :: mirrored
        integer(IK) :: sindex(nintmax), neval, nint
        real(RKC) :: abserr, sinfo(4, nintmax), betaMinus1, alphaMinus1
        CHECK_ASSERTION(__LINE__, 0._RKC < beta, SK_"@setBetaInc(): The condition `0. < beta` must hold. beta = "//getStr(beta)) ! fpp
        CHECK_ASSERTION(__LINE__, 0._RKC < alpha, SK_"@setBetaInc(): The condition `0. < alpha` must hold. alpha = "//getStr(alpha)) ! fpp
        mirrored = signed .and. x < 0._RKC
        if (mirrored) x = -x
        if(0._RKC < x .and. x < 1._RKC) then
            if (.not. mirrored) then
                mirrored = (alpha + 1._RKC) / (alpha + beta + 2._RKC) < x
                if (mirrored) x = 1._RKC - x
            endif
            if (mirrored) then ! swap
                abserr = beta
                beta = alpha
                alpha = abserr
            endif
            ! Prepare the integrand parameters.
            betaMinus1 = beta - 1._RKC
            alphaMinus1 = alpha - 1._RKC
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
                    betaInc = 1._RKC - betaInc
                end if
            end if
        elseif (0._RKC == x .or. x == 1._RKC) then
            info = 0_IK
            betaInc = x
        else
            info = 1_IK
            CHECK_ASSERTION(__LINE__, info == 0_IK, SK_"@setBetaInc(): The condition `0. <= x .and. x <= 1.` must hold. x = "//getStr(x)) ! fpp
        endif
    contains
        pure function getFunc(xx) result(func)
            real(RKC), intent(in) :: xx
            real(RKC) :: func
            func = exp(log(xx) * alphaMinus1 + log(1._RKC - xx) * betaMinus1 - logFuncBeta)
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setBetaInc_ENABLED && Def_ENABLED && 1
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! This implementation is partially inspired by the proposed approach of Numerical Recipes (1992) by the Great Bill Press et al.
        logical(LK)                 :: mirrored
        integer(IK)                 :: iter, iterTwice
        integer(IK) , parameter     :: MAX_ITER = 10000_IK
        real(RKC)   , parameter     :: EPS = epsilon(0._RKC), FPMIN = tiny(0._RKC) / EPS
        REAL(RKC)                   :: aa, c, d, delta, sumAlphaBeta, alphamMinusOne, alphamPlusOne
        CHECK_ASSERTION(__LINE__, 0._RKC < beta, SK_"@setBetaInc(): The condition `0._RKC < beta` must hold. beta = "//getStr(beta)) ! fpp
        CHECK_ASSERTION(__LINE__, 0._RKC < alpha, SK_"@setBetaInc(): The condition `0._RKC < alpha` must hold. alpha = "//getStr(alpha)) ! fpp
        mirrored = signed .and. x < 0._RKC
        if (mirrored) x = -x
        if(0._RKC < x .and. x < 1._RKC) then
            sumAlphaBeta = alpha + beta
            if (.not. mirrored) then
                mirrored = (alpha + 1._RKC) / (sumAlphaBeta + 2._RKC) < x
                if (mirrored) x = 1._RKC - x
            endif
            if (mirrored) then
                d = beta
                beta = alpha
                alpha = d
            endif
            alphamPlusOne = alpha + 1._RKC
            alphamMinusOne = alpha - 1._RKC
            d = 1._RKC - sumAlphaBeta * x / alphamPlusOne
            if (abs(d) < FPMIN) d = FPMIN
            c = 1._RKC
            d = 1._RKC / d
            betaInc = d
            do iter = 1, MAX_ITER
                iterTwice = 2_IK * iter
                aa = iter * (beta - iter) * x / ((alphamMinusOne + iterTwice) * (alpha + iterTwice))
                d = 1._RKC + aa * d
                if (abs(d) < FPMIN) d = FPMIN
                c = 1._RKC + aa / c
                if (abs(c) < FPMIN) c = FPMIN
                d = 1._RKC / d
                betaInc = betaInc * d * c
                aa = -(alpha + iter) * (sumAlphaBeta + iter) * x / ((alpha + iterTwice) * (alphamPlusOne + iterTwice))
                d = 1._RKC + aa * d
                if (abs(d) < FPMIN) d = FPMIN
                c = 1._RKC + aa / c
                if (abs(c) < FPMIN) c = FPMIN
                d = 1._RKC / d
                delta = d * c
                betaInc = betaInc * delta
                if (EPS < abs(delta - 1._RKC)) cycle
               !betaInc = betaInc * exp(alpha * log(x) + beta * getLog1p(-x) - logFuncBeta) / alpha
                betaInc = betaInc * exp(alpha * log(x) + beta * log(1._RKC - x) - logFuncBeta) / alpha
                !if (mirrored) betaInc = 1._RKC - betaInc ! 1._RKC is the regularization term.
                !if (mirrored) betaInc = -betaInc ! 1._RKC is the regularization term.
                if (mirrored) then
                    if (signed) then
                        betaInc = -betaInc
                    else
                        betaInc = 1._RKC - betaInc
                    end if
                end if
                info = 0_IK
                return
            end do
            info = 1_IK
        elseif (0._RKC == x .or. x == 1._RKC) then
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
        real(RKC)   :: ai, sumAlphaBeta, rx, temp, term
        real(RKC)   , parameter :: EPS = epsilon(0._RKC)**(7._RKC/8._RKC)
        integer(IK) , parameter :: MAX_ITER = 10000_IK
        CHECK_ASSERTION(__LINE__, 0._RKC < beta, SK_"@setBetaInc(): The condition `0._RKC < beta` must hold. beta = "//getStr(beta)) ! fpp
        CHECK_ASSERTION(__LINE__, 0._RKC < alpha, SK_"@setBetaInc(): The condition `0._RKC < alpha` must hold. alpha = "//getStr(alpha)) ! fpp
        mirrored = signed .and. x < 0._RKC
        if (mirrored) x = -x
        if (0._RKC < x .and. x < 1._RKC) then

            ! Change tail if necessary.

            sumAlphaBeta = alpha + beta
            if (.not. mirrored) then
                !mirrored = (alpha + 1._RKC) / (sumAlphaBeta + 2._RKC) < x ! this is the original netlib condition.
                mirrored = alpha < sumAlphaBeta * x ! This is the implementation by press et al.
                if (mirrored) x = 1._RKC - x
            endif
            if (mirrored) then
                temp = beta
                beta = alpha
                alpha = temp
            endif

            ai = 1._RKC
            term = 1._RKC
            betaInc = 1._RKC
            ns = int(beta + (1._RKC - x) * sumAlphaBeta, IK)
            ! Use the Soper reduction formula.
            rx = x / (1._RKC - x)
            temp = beta - ai
            if (ns == 0_IK) rx = x
            do iter = 1, MAX_ITER
                term = term * temp * rx / (alpha + ai)
                betaInc = betaInc + term
                temp = abs(term)
                if (temp <= EPS .and. temp <= EPS * betaInc) then
                   !betaInc = betaInc * exp(alpha * log(x) + (beta - 1._RKC) * getLog1p(-x) - logFuncBeta) / alpha
                    betaInc = betaInc * exp(alpha * log(x) + (beta - 1._RKC) * log(1._RKC - x) - logFuncBeta) / alpha
                    if (mirrored) then
                        if (signed) then
                            betaInc = -betaInc
                        else
                            betaInc = 1._RKC - betaInc
                        end if
                    end if
                    info = 0_IK
                    return
                end if
                ai = ai + 1._RKC
                ns = ns - 1
                if (ns < 0_IK) then
                    temp = sumAlphaBeta
                    sumAlphaBeta = sumAlphaBeta + 1._RKC
                else
                    temp = beta - ai
                    if (ns == 0_IK) rx = x
                end if
            end do
            info = 1_IK
        elseif (0._RKC == x .or. x == 1._RKC) then
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
        real(RKC)   , parameter :: underflow = tiny(0._RKC)
        real(RKC)   , parameter :: experr = log10(underflow)
        real(RKC)   , parameter :: ONE_THIRD = 1._RKC / 3._RKC
        real(RKC)   , parameter :: ONE_SIXTH = 1._RKC / 6._RKC
        real(RKC)   , parameter :: FIVE_SIXTH = 5._RKC / 6._RKC
        real(RKC)   , parameter :: ONE_MEPS = 1._RKC - sqrt(epsilon(0._RKC))
        real(RKC)               :: tolerance, adj, g, h, prev, r, s, sq, t, tx, w, betaIncOld, betaIncNew
        CHECK_ASSERTION(__LINE__, 0._RKC < alpha, SK_"@setBetaInv(): The condition `0. < alpha` must hold. alpha = "//getStr(alpha)) ! fpp
        CHECK_ASSERTION(__LINE__, 0._RKC < beta, SK_"@setBetaInv(): The condition `0. < beta` must hold. alpha = "//getStr(beta)) ! fpp
        mirrored = signed .and. betaInc < 0._RKC
        if (mirrored) betaInc = -betaInc
        if (0._RKC < betaInc .and. betaInc < 1._RKC) then

            ! Change tail if necessary.

            if (mirrored) then
                r = alpha
                alpha = beta
                beta = r
            else
                mirrored = 0.5_RKC < betaInc
                if (mirrored) then
                    betaInc = 1._RKC - betaInc
                    r = alpha
                    alpha = beta
                    beta = r
                end if
            end if

            ! Compute the initial approximation.

            r = sqrt(-log(betaInc**2))
            betaIncNew = r - (2.30753_RKC + 0.27061_RKC * r) / (1._RKC + (0.99229_RKC + 0.04481_RKC * r) * r)
            if (1._RKC < alpha .and. 1._RKC < beta) then
                r = (betaIncNew * betaIncNew - 3._RKC) * ONE_SIXTH
                s = 1._RKC / (alpha + alpha - 1._RKC)
                t = 1._RKC / (beta + beta - 1._RKC)
                h = 2._RKC / (s + t)
                w = betaIncNew * sqrt(h + r) / h - (t - s) * (r + FIVE_SIXTH - 2._RKC / (3._RKC * h))
                betaInv = alpha / (alpha + beta * exp(w + w))
            else
                r = beta + beta
                t = 1._RKC / (9._RKC * beta)
                t = r * (1._RKC - t + betaIncNew * sqrt(t))**3
                if (t <= 0._RKC) then
                    betaInv = 1._RKC - exp((log((1._RKC - betaInc) * beta) + logFuncBeta) / beta)
                else
                    t = (4._RKC * alpha + r - 2._RKC) / t
                    if (t <= 1._RKC) then
                        betaInv = exp((log(betaInc * alpha) + logFuncBeta) / alpha)
                    else
                        betaInv = 1._RKC - 2._RKC / (t + 1._RKC)
                    end if
                end if
            end if

            ! Solve for X via modified Newton-Raphson method, using the function `setBetaInc()`.

            t = 1._RKC - beta
            r = 1._RKC - alpha
            betaIncOld = 0._RKC
            prev = 1._RKC
            sq = 1._RKC
            if (betaInv < 0.0001_RKC) betaInv = 0.0001_RKC
            if (0.9999_RKC < betaInv) betaInv = 0.9999_RKC
            tolerance = 10._RKC**int(max(-5._RKC / alpha**2 - 1._RKC / betaInc**0.2 - 13._RKC, experr), IK)

            ! Begin iteration.

            loopFindRoot: do ! 10
                call setBetaInc(betaIncNew, betaInv, alpha, beta, logFuncBeta, signed, info)
                if (info /= 0_IK) return ! LCOV_EXCL_LINE
                if (betaIncNew < 0._RKC) then
                    !if (abs(betaIncNew) < epsilon(0._RKC)) error stop "Underflow detected."
                    betaIncNew = betaIncNew - betaInc
                    !if (abs(betaIncNew) < epsilon(0._RKC)) error stop "Underflow detected."
                    betaIncNew = betaIncNew + 1._RKC
                else
                    betaIncNew = betaIncNew - betaInc
                end if
               !betaIncNew = betaIncNew * exp(logFuncBeta + r * log(betaInv) + t * getLog1p(-betaInv))
                betaIncNew = betaIncNew * exp(logFuncBeta + r * log(betaInv) + t * log(1._RKC - betaInv))
                if (betaIncNew * betaIncOld <= 0._RKC) prev = max(sq, underflow)
                g = 1._RKC
                loopRefine: do ! 20
                    adj = g * betaIncNew
                    sq = adj * adj
                    if (sq < prev) then
                        tx = betaInv - adj
                        if (tx < 0._RKC .or. 1._RKC < tx) then
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
                    if (tx /= 0._RKC .and. tx /= 1._RKC) exit loopRefine
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
                    betaInv = 1._RKC - betaInv
                end if
            end if

        elseif (0._RKC == betaInc .or. betaInc == 1._RKC) then

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
        real(RKC)   , parameter :: EPS = sqrt(epsilon(0._RKC))
        real(RKC)   :: pp, t, u, betaIncNew, al, h, w, alphaMinus1, betaMinus1
        if (signed .and. betaInc < 0._RKC) betaInc = 1 + betaInc
        if (0._RKC < betaInc .and. betaInc < 1._RKC) then
            betaMinus1 = beta - 1._RKC
            alphaMinus1 = alpha - 1._RKC
            if (alpha < 1._RKC .or. beta < 1._RKC) then
                u = exp(beta * log(beta / (alpha + beta))) / beta
                t = exp(alpha * log(alpha / (alpha + beta))) / alpha
                w = t + u
                if (betaInc < t / w) then
                    betaInv = (alpha * w * betaInc)**(1._RKC / alpha)
                else
                    betaInv = 1._RKC - (beta * w * (1._RKC - betaInc))**(1._RKC / beta)
                end if
            else
                if (betaInc < 0.5_RKC) then
                    pp = betaInc
                else
                    pp = 1._RKC - betaInc
                end if
                t = sqrt(-2._RKC * log(pp))
                betaInv = (2.30753_RKC + t * 0.27061_RKC) / (1._RKC + t * (0.99229_RKC + t * 0.04481_RKC)) - t
                if (betaInc < .5_RKC) betaInv = -betaInv
                al = (betaInv**2 - 3._RKC) / 6._RKC
                h = 2._RKC / (1._RKC / (2._RKC * alpha - 1._RKC) + 1._RKC / (2._RKC * beta - 1._RKC))
                w = (betaInv * sqrt(al + h) / h) - (1._RKC / (2._RKC * beta - 1._RKC) - 1._RKC / (2._RKC * alpha - 1._RKC)) * (al + 5._RKC / 6._RKC - 2._RKC / (3._RKC * h))
                betaInv = alpha / (alpha + beta * exp(2._RKC * w))
            end if
            loopFindRoot: do iter = 1, 10
                if (betaInv == 0._RKC .or. betaInv == 1._RKC) exit loopFindRoot
                call setBetaInc(betaIncNew, betaInv, alpha, beta, logFuncBeta, signed, info)
                if (info /= 0_IK) return
                betaIncNew = betaIncNew - betaInc
                t = exp(alphaMinus1 * log(betaInv) + betaMinus1 * log(1._RKC - betaInv) - logFuncBeta)
                u = betaIncNew / t
                t = u / (1._RKC - .5_RKC * min(1._RKC, u * (alphaMinus1 / betaInv - betaMinus1 / (1._RKC - betaInv))))
                betaInv = betaInv - t
                if (betaInv <= 0._RKC) betaInv = .5_RKC * (betaInv + t)
                if (betaInv >= 1._RKC) betaInv = .5_RKC * (betaInv + t + 1._RKC)
                if (abs(t) < EPS * betaInv .and. 1_IK < iter) exit loopFindRoot
            end do loopFindRoot
            info = 0_IK
        elseif (0._RKC == betaInc .or. betaInc == 1._RKC) then
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
        real(RKC)   :: lb, ub, lf, uf, abstol, h, w, r, s, t, betaIncNew
        real(RKC)   , parameter :: ONE_SIXTH = 1._RKC / 6._RKC
        !real(RKC)   , parameter :: abstol = epsilon(0._RKC)
        CHECK_ASSERTION(__LINE__, 0._RKC < alpha, SK_"@setBetaInv(): The condition `0. < alpha` must hold. alpha = "//getStr(alpha)) ! fpp
        CHECK_ASSERTION(__LINE__, 0._RKC < beta, SK_"@setBetaInv(): The condition `0. < beta` must hold. alpha = "//getStr(beta)) ! fpp
        info = 0_IK
        mirrored = signed .and. betaInc < 0._RKC
        if (mirrored) betaInc = -betaInc
        if (0._RKC < betaInc .and. betaInc < 1._RKC) then
            if (.not. mirrored) then
                mirrored = 0.5_RKC < betaInc
                if (mirrored) betaInc = 1._RKC - betaInc
            end if
            if (mirrored) then
                r = alpha
                alpha = beta
                beta = r
            end if
            ! Compute the initial approximation.
            r = sqrt(-log(betaInc**2))
            betaIncNew = r - (2.30753_RKC + 0.27061_RKC * r) / (1._RKC + (0.99229_RKC + 0.04481_RKC * r) * r)
            if (1._RKC < alpha .and. 1._RKC < beta) then
                r = (betaIncNew * betaIncNew - 3._RKC) * ONE_SIXTH
                s = 1._RKC / (alpha + alpha - 1._RKC)
                t = 1._RKC / (beta + beta - 1._RKC)
                h = 2._RKC / (s + t)
                w = betaIncNew * sqrt(h + r) / h - (t - s) * (r + 5._RKC / 6._RKC - 2._RKC / (3._RKC * h))
                betaInv = alpha / (alpha + beta * exp(w + w))
            else
                r = beta + beta
                t = 1._RKC / (9._RKC * beta)
                t = r * (1._RKC - t + betaIncNew * sqrt(t))**3
                if (t <= 0._RKC) then
                    betaInv = 1._RKC - exp((log((1._RKC - betaInc) * beta) + logFuncBeta) / beta)
                else
                    t = (4._RKC * alpha + r - 2._RKC) / t
                    if (t <= 1._RKC) then
                        betaInv = exp((log(betaInc * alpha) + logFuncBeta) / alpha)
                    else
                        betaInv = 1._RKC - 2._RKC / (t + 1._RKC)
                    end if
                end if
            end if
            abstol = epsilon(0._RKC)
            if (betaInv < 0._RKC) betaInv = sqrt(abstol)
            if (1._RKC < betaInv) betaInv = 1._RKC - sqrt(abstol)
            lb = 0._RKC
            ub = 1._RKC
            lf = betaInc
            uf = betaInc - 1._RKC
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
                    betaInv = 1._RKC - betaInv
                end if
            end if

        elseif (0._RKC == betaInc .or. betaInc == 1._RKC) then
            info = 0_IK
            betaInv = betaInc
        else
            info = 1_IK
            CHECK_ASSERTION(__LINE__, info == 0_IK, SK_"@setBetaInc(): The condition `0. <= betaInc .and. betaInc <= 1.` must hold. betaInc = "//getStr(betaInc)) ! fpp
        end if
    contains
        PURE function getFunc(x) result(func)
            real(RKC)   , intent(in)    :: x ! betaInv guess.
            real(RKC)                   :: func
            integer(IK)                 :: info
            call setBetaInc(func, x, alpha, beta, logFuncBeta, signed, info)
            if (info /= 0_IK) error stop ! LCOV_EXCL_LINE
            func = betaInc - func
        end function
        PURE function getFuncDiff(x, order) result(betaPDF)
            integer(IK) , intent(in)    :: order
            real(RKC)   , intent(in)    :: x
            real(RKC)   , parameter     :: SQRT_HUGE = sqrt(huge(x))
            real(RKC)                   :: betaPDF
            if(order == 0_IK) then
                betaPDF = getFunc(x)
            else
                if (x /= 0._RKC .and. x /= 1._RKC) then
                    call setBetaLogPDF(betaPDF, x, alpha, beta)
                    betaPDF = exp(betaPDF)
                elseif ((x == 0._RKC .and. alpha < 1._RKC) .or. (x == 1._RKC .and. beta < 1._RKC)) then
                    betaPDF = SQRT_HUGE
                else
                    betaPDF = 0._RKC
                end if
            end if
        end function

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif