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
!>  This include file contains the implementation of procedures in [pm_mathRoot](@ref pm_mathRoot).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CHECK_LF(LF) \
CHECK_ASSERTION(__LINE__, abs(LF - lf) < abstol, SK_"@setRoot(): The condition `abs(getFunc(lb), lf) < abstol` must hold. lb, getFunc(lb), lf, abstol = "//getStr([lb, LF, lf, abstol]))
#define CHECK_UF(UF) \
CHECK_ASSERTION(__LINE__, abs(UF - uf) < abstol, SK_"@setRoot(): The condition `abs(getFunc(ub), uf) < abstol` must hold. ub, getFunc(ub), uf, abstol = "//getStr([ub, UF, uf, abstol]))
#define CHECK_BRACKET \
CHECK_ASSERTION(__LINE__, sign(lf, 1._RKG) /= sign(uf, 1._RKG), \
SK_"@setRoot(): The condition `sign(lf, 1.) /= sign(uf, 1.)` must hold. lf, uf = "//getStr([lf, uf])) ! fpp
#define CHECK_LB_LT_UB \
CHECK_ASSERTION(__LINE__, lb < ub, SK_"@setRoot(): The condition `lb < ub` must hold. lb, ub = "//getStr([lb, ub])) ! fpp
#define CHECK_ZERO_LT_ABSTOL \
CHECK_ASSERTION(__LINE__, 0._RKG < abstol, SK_"@setRoot(): The condition `0. < abstol` must hold. abstol = "//getStr(abstol)) ! fpp
#define CHECK_ABSTOL_LT_LB_UB_DIFF \
CHECK_ASSERTION(__LINE__, abstol < ub - lb, \
SK_"@setRoot(): The condition `abstol < ub - lb` must hold. abstol, ub, lb = "//getStr([abstol, ub, lb])) ! fpp
        real(RKG), parameter :: EPS10 = 10 * epsilon(0._RKG)
#if     Fixed_ENABLED
#define CHECK_ZERO_LT_NITER
        integer(IK), parameter :: NITER = ceiling(50 * precision(0._RKG) / log10(2._RKG))
#elif   Niter_ENABLED
#define CHECK_ZERO_LT_NITER \
CHECK_ASSERTION(__LINE__, 0_IK < niter, SK_"@setRoot(): The condition `0 < niter` must hold. niter = "//getStr(niter)) ! fpp
#elif   !getRoot_ENABLED
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getRoot_ENABLED && Def_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        root = getRoot(brent, getFunc, lb, ub, abstol, neval, niter)

        !%%%%%%%%%%%%%%
#elif   getRoot_ENABLED
        !%%%%%%%%%%%%%%

        real(RKG)   :: abstol_def, lf, uf
        integer(IK) :: neval_def
        if (present(abstol)) then
            abstol_def = abstol
        else
            abstol_def = epsilon(0._RKG)**.8 * (abs(lb) + abs(ub))
        end if
#if     Newton_ENABLED || Halley_ENABLED || Schroder_ENABLED
        lf = getFunc(lb, 0_IK)
        uf = getFunc(ub, 0_IK)
        if (present(init)) then
            root = init
        else
            root = 0.5_RKG * (lb + ub)
        end if
#elif   False_ENABLED || Bisection_ENABLED || Secant_ENABLED || Brent_ENABLED || Ridders_ENABLED || TOMS748_ENABLED
        lf = getFunc(lb)
        uf = getFunc(ub)
#else
#error  "Unrecognized interface."
#endif
        if (present(niter)) then
            call setRoot(method, getFunc, root, lb, ub, lf, uf, abstol_def, neval_def, niter)
        else
            call setRoot(method, getFunc, root, lb, ub, lf, uf, abstol_def, neval_def)
        end if
        if (present(neval)) then
            neval = neval_def + 2_IK
        elseif (neval_def < 0_IK) then
            error stop MODULE_NAME//SK_"@getRoot(): Failed to converge in search for the root of the specified function." ! LCOV_EXCL_LINE
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setRoot_ENABLED && False_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: diff, interval, func
        if (abs(lf) <= EPS10) then
            neval = 0_IK
            root = lb
        elseif (abs(uf) <= EPS10) then
            neval = 0_IK
            root = ub
        elseif (abs(ub - lb) < abstol) then
            root = .5_RKG * (ub + lb)
            neval = 0_IK
        else
            CHECK_BRACKET
            CHECK_LB_LT_UB
            CHECK_ZERO_LT_ABSTOL
            CHECK_LF(getFunc(lb))
            CHECK_UF(getFunc(ub))
            CHECK_ABSTOL_LT_LB_UB_DIFF
            if (0._RKG <= lf) then
                diff = lb
                lb = ub
                ub = diff
                diff = lf
                lf = uf
                uf = diff
            end if
            interval = ub - lb
            do neval = 1_IK, niter
                root = lb + interval * lf / (lf - uf)
                func = getFunc(root)
                if (func < 0._RKG) then
                    diff = lb - root
                    lb = root
                    lf = func
                else
                    diff = ub - root
                    ub = root
                    uf = func
                end if
                interval = ub - lb
                if (abs(diff) < abstol .or. func == 0._RKG) return
            end do
            ! Error occurred.
            root = .5_RKG * (lb + ub) ! for the sake of defining `root` on output (important for tests).
            neval = -neval
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setRoot_ENABLED && Bisection_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: delta, center
        if (abs(lf) <= EPS10) then
            neval = 0_IK
            root = lb
        elseif (abs(uf) <= EPS10) then
            neval = 0_IK
            root = ub
        elseif (abs(ub - lb) < abstol) then
            root = .5_RKG * (ub + lb)
            neval = 0_IK
        else
            CHECK_BRACKET
            CHECK_LB_LT_UB
            CHECK_ZERO_LT_NITER
            CHECK_ZERO_LT_ABSTOL
            CHECK_LF(getFunc(lb))
            CHECK_UF(getFunc(ub))
            CHECK_ABSTOL_LT_LB_UB_DIFF
            if (lf < 0._RKG) then
                root = lb
                delta = ub - lb
            else
                root = ub
                delta = lb - ub
            end if
            do neval = 3_IK, NITER + 2_IK
                delta = delta * 0.5_RKG
                center = root + delta
                uf = getFunc(center)
                if (uf <= 0._RKG) root = center
                if (abs(delta) < abstol .or. uf == 0._RKG) return
            end do
            ! Error occurred.
            root = .5_RKG * (lb + ub) ! for the sake of defining `root` on output (important for tests).
            neval = -neval
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setRoot_ENABLED && Secant_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: x, delta
        if (abs(lf) <= EPS10) then
            neval = 0_IK
            root = lb
        elseif (abs(uf) <= EPS10) then
            neval = 0_IK
            root = ub
        elseif (abs(ub - lb) < abstol) then
            root = .5_RKG * (ub + lb)
            neval = 0_IK
        else
            CHECK_BRACKET
            CHECK_LB_LT_UB
            CHECK_ZERO_LT_NITER
            CHECK_ZERO_LT_ABSTOL
            CHECK_LF(getFunc(lb))
            CHECK_UF(getFunc(ub))
            CHECK_ABSTOL_LT_LB_UB_DIFF
            if (abs(lf) < abs(uf)) then
                root = lb
                x = ub
                delta = lf
                lf = uf
                uf = delta
            else
                x = lb
                root = ub
            end if
            do neval = 1_IK, niter
                delta = (x - root) * uf / (uf - lf)
                x = root
                lf = uf
                root = root + delta
                uf = getFunc(root)
                if (abs(delta) < abstol .or. uf == 0._RKG) return
            end do
            ! Error occurred.
            root = .5_RKG * (lb + ub) ! for the sake of defining `root` on output (important for tests).
            neval = -neval
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setRoot_ENABLED && Brent_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG)   :: halfAbsTol
        real(RKG)   :: middle, interval, funMid, tol1, halfint, rnew, p, q, r, ratio
        real(RKG)   , parameter :: EPS = epsilon(lb)
        if (abs(lf) <= EPS10) then
            neval = 0_IK
            root = lb
        elseif (abs(uf) <= EPS10) then
            neval = 0_IK
            root = ub
        elseif (abs(ub - lb) < abstol) then
            root = .5_RKG * (ub + lb)
            neval = 0_IK
        else
            CHECK_BRACKET
            CHECK_LB_LT_UB
            CHECK_ZERO_LT_NITER
            CHECK_ZERO_LT_ABSTOL
            CHECK_LF(getFunc(lb))
            CHECK_UF(getFunc(ub))
            CHECK_ABSTOL_LT_LB_UB_DIFF
            halfAbsTol = 0.5_RKG * abstol
            root = ub
            middle = ub
            funMid = uf
            do neval = 0, niter
                ! adjust the search interval, if needed.
                if ((funMid < 0._RKG .and. uf < 0._RKG) .or. (0._RKG < uf .and. 0._RKG < funMid)) then ! same sign
                    middle = lb
                    funMid = lf
                    interval = root - lb
                    rnew = interval
                end if
                if (abs(funMid) < abs(uf)) then
                    lb = root
                    root = middle
                    middle = lb
                    lf = uf
                    uf = funMid
                    funMid = lf
                end if
                halfint = .5_RKG * (middle - root)
                tol1 = 2._RKG * EPS * abs(root) + halfAbsTol
                if (abs(halfint) < tol1 .or. uf == 0._RKG) return
                ! See if a bisection is needed.
                if ((abs(rnew) >= tol1) .and. (abs(lf) > abs(uf))) then
                    ratio = uf / lf ! 50
                    if (lb /= middle) then ! Try inverse quadratic interpolation.
                        q = lf / funMid
                        r = uf / funMid
                        p = ratio * (2._RKG * halfint * q * (q - r) - (root - lb) * (r - 1._RKG))
                        q = (q - 1._RKG) * (r - 1._RKG) * (ratio - 1._RKG)
                    else ! Try linear interpolation.
                        p = 2._RKG * halfint * ratio
                        q = 1._RKG - ratio
                    end if
                    if (0._RKG < p) then ! Check inbounds.
                        q = -q
                    else
                        p = -p
                    end if
                    if (((2._RKG * p) >= (3._RKG * halfint * q - abs(tol1 * q))) .or. (p >= abs(0.5_RKG * rnew * q))) then ! Interpolation failed, use bisection.
                        interval = halfint
                        rnew = interval
                    else ! Accept interpolation.
                        rnew = interval
                        interval = p / q
                    end if
                else ! Bounds decreasing too slowly, use bisection.
                    interval = halfint
                    rnew = interval
                end if
                ! Move last best guess to a.
                lb = root ! 110
                lf = uf
                ! Evaluate new trial root.
                if (tol1 < abs(interval)) then
                    root = root + interval
                elseif (0._RKG < halfint) then
                    root = root + tol1
                else
                    root = root - tol1
                end if
                uf = getFunc(root)
            end do
            ! Error occurred.
            root = .5_RKG * (lb + root) ! for the sake of defining `root` on output (important for tests).
            neval = -neval
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setRoot_ENABLED && Ridders_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: halfint, rold, fold, func, multiplicity
        neval = 0_IK
        if (abs(lf) <= EPS10) then
            root = lb
        elseif (abs(uf) <= EPS10) then
            root = ub
        elseif (abs(ub - lb) < abstol) then
            root = .5_RKG * (ub + lb)
        else
            CHECK_BRACKET
            CHECK_LB_LT_UB
            CHECK_ZERO_LT_NITER
            CHECK_ZERO_LT_ABSTOL
            CHECK_LF(getFunc(lb))
            CHECK_UF(getFunc(ub))
            CHECK_ABSTOL_LT_LB_UB_DIFF
            !rnew = huge(rold)**.66
            do
                halfint = 0.5_RKG * (ub - lb)
                root = lb + halfint
                func = getFunc(root)
                neval = neval + 1_IK
                if (niter < neval) exit
                multiplicity = sqrt(func**2 - lf * uf)
                if (multiplicity == 0._RKG) return ! LCOV_EXCL_LINE ! This should rarely occur.
                rold = root
                fold = func
                root = rold + halfint * sign(fold, lf - uf) / multiplicity
                if (abs(root - rold) < abstol) return ! This should never occur in the first iteration.
                func = getFunc(root)
                neval = neval + 1_IK
                if (niter < neval) exit
                if (func == 0._RKG) return ! LCOV_EXCL_LINE
                if (sign(fold, func) /= fold) then
                    lb = rold
                    lf = fold
                    uf = func
                    ub = root
                elseif (sign(lf, func) /= lf) then
                    ub = root
                    uf = func
                elseif (sign(uf, func) /= uf) then
                    lb = root
                    lf = func
                else
                    error stop MODULE_NAME//SK_"@setRootRidders(): Internal library error detected. This condition should never occur."
                end if
                if (abs(ub - lb) < abstol) exit
            end do
            if (niter < neval) then
                neval = -neval ! Error occurred.
            else
                root = .5_RKG * (lb + ub)
            end if
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setRoot_ENABLED && TOMS748_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) , parameter :: EPS5 = 5 * epsilon(0._RKG), SMALL = tiny(0._RKG) / epsilon(0._RKG) !* 32
        real(RKG) :: rmid, rupp, fupp, rold, fold, rnew, fnew, interval
        neval = 0_IK
        if (abs(lf) <= EPS10) then
            root = lb
        elseif (abs(uf) <= EPS10) then
            root = ub
        elseif (abs(ub - lb) < abstol) then
            root = .5_RKG * (ub + lb)
        else
            CHECK_BRACKET
            CHECK_LB_LT_UB
            CHECK_ZERO_LT_NITER
            CHECK_ZERO_LT_ABSTOL
            CHECK_LF(getFunc(lb))
            CHECK_UF(getFunc(ub))
            CHECK_ABSTOL_LT_LB_UB_DIFF
            !fnew = 1.e5_RKG
            !rnew = 1.e5_RKG
            !fold = 1.e5_RKG
            blockConvergence: block
                ! On the first step we take lb secant step.
                rmid = lb - lf * (ub - lb) / (uf - lf)
                if (rmid <= lb + abs(lb) * EPS5 .or. ub - abs(ub) * EPS5 <= rmid) rmid = .5_RKG * (lb + ub)
                call setBracket(getFunc, lb, ub, rmid, lf, uf, rold, fold)
                neval = neval + 1_IK
                ! Take a quadratic interpolation on the second step.
                call setPolIntQuad(rmid, lb, ub, rold, lf, uf, fold, nstep = 2_IK)
                rnew = rold
                fnew = fold
                neval = neval + 1_IK
                call setBracket(getFunc, lb, ub, rmid, lf, uf, rold, fold)
                loopConvergence: do
                    interval = ub - lb
                    if (niter < neval .or. lf == 0._RKG .or. abs(ub - lb) < abstol) exit blockConvergence
                    ! Starting with the third step taken use either quadratic or cubic interpolation.
                    ! Cubic interpolation requires that all four function values lf, uf, fold, and fnew are distinct.
                    ! Should that not be the case, take a quadratic step instead.
                    if ((abs(lf - uf) < SMALL) .or. & ! LCOV_EXCL_LINE
                        (abs(lf - fold) < SMALL) .or. & ! LCOV_EXCL_LINE
                        (abs(lf - fnew) < SMALL) .or. & ! LCOV_EXCL_LINE
                        (abs(uf - fold) < SMALL) .or. & ! LCOV_EXCL_LINE
                        (abs(uf - fnew) < SMALL) .or. & ! LCOV_EXCL_LINE
                        (abs(fold - fnew) < SMALL)) then
                        call setPolIntQuad(rmid, lb, ub, rold, lf, uf, fold, nstep = 2_IK)
                    else
                        call setPolIntCubic(rmid, lb, ub, rold, rnew, lf, uf, fold, fnew)
                    end if
                    ! Re-bracket and check for termination.
                    rnew = rold
                    fnew = fold
                    neval = neval + 1_IK
                    call setBracket(getFunc, lb, ub, rmid, lf, uf, rold, fold)
                    if (niter < neval .or. lf == 0._RKG .or. abs(ub - lb) < abstol) exit blockConvergence
                    ! Next interpolated step.
                    if ((abs(lf - uf) < SMALL) .or. & ! LCOV_EXCL_LINE
                        (abs(lf - fold) < SMALL) .or. & ! LCOV_EXCL_LINE
                        (abs(lf - fnew) < SMALL) .or. & ! LCOV_EXCL_LINE
                        (abs(uf - fold) < SMALL) .or. & ! LCOV_EXCL_LINE
                        (abs(uf - fnew) < SMALL) .or. & ! LCOV_EXCL_LINE
                        (abs(fold - fnew) < SMALL)) then
                        call setPolIntQuad(rmid, lb, ub, rold, lf, uf, fold, nstep = 3_IK)
                    else
                        call setPolIntCubic(rmid, lb, ub, rold, rnew, lf, uf, fold, fnew)
                    end if
                    ! Bracket, check termination condition, and update rnew.
                    neval = neval + 1_IK
                    call setBracket(getFunc, lb, ub, rmid, lf, uf, rold, fold)
                    if (niter < neval .or. lf == 0._RKG .or. abs(ub - lb) < abstol) exit blockConvergence
                    ! Take a double-length secant step.
                    if (abs(lf) < abs(uf)) then
                        rupp = lb
                        fupp = lf
                    else
                        rupp = ub
                        fupp = uf
                    end if
                    rmid = rupp - 2 * fupp * (ub - lb) / (uf - lf)
                    if (abs(rmid - rupp) > .5_RKG * (ub - lb)) rmid = lb + .5_RKG * (ub - lb)
                    ! Bracket and check termination condition.
                    rnew = rold
                    fnew = fold
                    neval = neval + 1_IK
                    call setBracket(getFunc, lb, ub, rmid, lf, uf, rold, fold)
                    if (niter < neval .or. lf == 0._RKG .or. abs(ub - lb) < abstol) exit blockConvergence
                    if (2 * (ub - lb) < interval) cycle
                    ! Take an additional bisection step due to slow convergence.
                    rnew = rold
                    fnew = fold
                    neval = neval + 1_IK
                    call setBracket(getFunc, lb, ub, lb + .5_RKG * (ub - lb), lf, uf, rold, fold)
                    if (niter < neval .or. lf == 0._RKG .or. abs(ub - lb) < abstol) exit blockConvergence
                end do loopConvergence
            end block blockConvergence
            if (niter < neval) neval = -neval
            if (abs(lf) < EPS10) ub = lb
            if (abs(uf) < EPS10) lb = ub
            root = .5_RKG * (ub + lb)
        end if

    contains

        subroutine setBracket(getFunc, lb, ub, rmid, lf, uf, rold, fold)
            !   Given a point `rmid` inside the existing enclosing interval `[lb, ub]`:
            !   1.  set lb = rmid if `getFunc(rmid) == 0`, otherwise,
            !   2.  find the new enclosing interval, either `[lb, rmid]` or `[rmid, ub]`,
            !       and set `rold` and `fold` to the point that has just been removed from the interval.
            !       In other words `rold` is the third best guess to the root.
            real(RKG) , parameter :: EPS2 = 2 * epsilon(0._RKG)
            real(RKG), intent(inout) :: lb, ub, lf, uf
            real(RKG), intent(out) :: rold, fold
            procedure(real(RKG)) :: getFunc
            real(RKG), value :: rmid
            real(RKG) :: fmid
            ! Adjust the location of `rmid` accordingly if the `[lb, ub]` is small, or `rmid` is too close to one interval end.
            if (ub - lb < 2 * EPS2 * lb) then
                rmid = lb + .5_RKG * (ub - lb)
            elseif (rmid <= lb + abs(lb) * EPS2) then
                rmid = lb + abs(lb) * EPS2
            elseif (ub - abs(ub) * EPS2 <= rmid) then
                rmid = ub - abs(ub) * EPS2
            end if
            fmid = getFunc(rmid)
            if (fmid == 0._RKG) then
                rold = 0._RKG
                fold = 0._RKG
                lf = 0._RKG
                lb = rmid
                return
            end if
            ! Update the interval.
            if ((lf < 0._RKG .and. 0._RKG < fmid) .or. (fmid < 0._RKG .and. 0._RKG < lf)) then
                rold = ub
                fold = uf
                ub = rmid
                uf = fmid
            else
                rold = lb
                fold = lf
                lb = rmid
                lf = fmid
            end if
        end subroutine

        ! Perform quadratic interpolation to determine the next point,
        ! by taking `niter - neval` Newton steps to find the location of the quadratic polynomial.
        ! `rold`, the third best approximation to the root, after `lb` and `ub`, must lie outside of the interval `[lb, ub]`.
        ! Fall back to a secant step should the result be out of range.
        ! Start by obtaining the coefficients of the quadratic polynomial.
        ! Compute division without undue over/under-flow.
        pure subroutine setPolIntQuad(rmid, lb, ub, rold, lf, uf, fold, nstep)
            real(RKG), parameter :: HUGE_RKG = huge(0._RKG)
            real(RKG), intent(in) :: lb, ub, rold, lf, uf, fold
            integer(IK), intent(in) :: nstep
            real(RKG), intent(out) :: rmid
            integer(IK) :: istep
            real(RKG) :: ratio1, ratio2, numer, denom, divres
#define SET_DIV(DIVRES,NUMER,DENOM,OFLOW) \
if (abs(DENOM) < 1._RKG) then; if (abs((DENOM) * HUGE_RKG) <= abs(NUMER)) then; \
DIVRES = OFLOW; else; DIVRES = (NUMER) / (DENOM); end if; else; DIVRES = (NUMER) / (DENOM); end if;
            SET_DIV(ratio1,uf - lf,ub - lb,HUGE_RKG)
            SET_DIV(ratio2,fold - uf,rold - ub,HUGE_RKG)
            SET_DIV(ratio2,ratio2 - ratio1,rold - lb,0._RKG)
            if (ratio2 == 0._RKG) then ! Failed to determine coefficients. Try a secant step.
                call setPolIntSecant(rmid, lb, ub, lf, uf)
                return
            end if
            ! Determine the starting point of the Newton steps.
            if ((ratio2 < 0._RKG .and. lf < 0._RKG) .or. (0._RKG < ratio2 .and. 0._RKG < lf)) then
                rmid = lb
            else
                rmid = ub
            end if
            ! Take the Newton steps.
            do istep = 1, nstep
                numer = lf + (ratio1 + ratio2 * (rmid - ub)) * (rmid - lb)
                denom = ratio1 + ratio2 * (2 * rmid - lb - ub)
                SET_DIV(divres,numer,denom,1 + rmid - lb)
                rmid = rmid - divres
            end do
            if (rmid <= lb .or. ub <= rmid) call setPolIntSecant(rmid, lb, ub, lf, uf) ! Failed. Try a secant step.
#undef      SET_DIV
        end subroutine

       ! Perform standard secant interpolation of `[lb, ub]` given function evaluations `(lf, uf)`.
       ! Perform a bisection if secant interpolation is very close to either `lb` or `ub`.
       ! This interpolation is needed when at least one other form of interpolation has already failed,
       ! implying that the function is unlikely to be smooth with a root near `lb` or `ub`.
        pure subroutine setPolIntSecant(rmid, lb, ub, lf, uf)
            real(RKG), intent(in) :: lb, ub, lf, uf
            real(RKG), intent(out) :: rmid
           rmid = lb - lf * (ub - lb) / (uf - lf)
           if (rmid <= lb + abs(lb) * EPS5 .or. ub - abs(ub) * EPS5 <= rmid) rmid = .5_RKG * (lb + ub)
        end subroutine

        ! Use inverse cubic interpolation of getFunc(x) at points `[lb, ub, rold, rnew]` to obtain an approximate root of getFunc(x).
        ! `rold` and `rnew` lie out of `[lb, ub]` and are the third and forth best approximations to the root that we have found so far.
        ! Fall back to quadratic interpolation in case of an erroneous result out of `[lb, ub]`.
        pure subroutine setPolIntCubic(rmid, lb, ub, rold, rnew, lf, uf, fold, fnew)
            real(RKG), intent(in) :: lb, ub, rold, rnew, lf, uf, fold, fnew
            real(RKG), intent(out) :: rmid
            real(RKG) :: q11, q21, q31, d21, d31, q22, q32, d32, q33
            q11 = (rold - rnew) * fold / (fnew - fold)
            q21 = (ub - rold) * uf / (fold - uf)
            q31 = (lb - ub) * lf / (uf - lf)
            d21 = (ub - rold) * fold / (fold - uf)
            d31 = (lb - ub) * uf / (uf - lf)

            q22 = (d21 - q11) * uf / (fnew - uf)
            q32 = (d31 - q21) * lf / (fold - lf)
            d32 = (d31 - q21) * fold / (fold - lf)
            q33 = (d32 - q22) * lf / (fnew - lf)
            rmid = q31 + q32 + q33 + lb
            if (rmid <= lb .or. ub <= rmid) call setPolIntQuad(rmid, lb, ub, rold, lf, uf, fold, nstep = 3_IK) ! Out of bounds step. Us quadratic interpolation.
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setRoot_ENABLED && Newton_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: iter
        real(RKG) :: func, grad, interval, func2grad, temp
        neval = 0_IK
        if (abs(lf) <= EPS10) then
            root = lb
        elseif (abs(uf) <= EPS10) then
            root = ub
        elseif (abs(ub - lb) < abstol) then
            root = .5_RKG * (ub + lb)
        else
            CHECK_BRACKET
            CHECK_LB_LT_UB
            CHECK_ZERO_LT_NITER
            CHECK_ZERO_LT_ABSTOL
            CHECK_LF(getFunc(lb, 0_IK))
            CHECK_UF(getFunc(ub, 0_IK))
            CHECK_ABSTOL_LT_LB_UB_DIFF
            if (root <= lb .or. ub <= root) root = .5_RKG * (lb + ub)
            if (0._RKG <= lf) then
                root = lb
                lb = ub
                ub = root
            end if
            interval = abs(ub - lb)
            func2grad = interval
            neval = neval + 1_IK
            func = getFunc(root, 0_IK)
            do iter = 1, niter
                neval = neval + 1_IK
                grad = getFunc(root, 1_IK)
                if (0._RKG < ((root - ub) * grad - func) * ((root - lb) * grad - func) .or. abs(interval * grad) < abs(2._RKG * func)) then
                    interval = func2grad
                    func2grad = .5_RKG * (ub - lb)
                    root = lb + func2grad
                    if (lb == root) return
                else
                    interval = func2grad
                    func2grad = func / grad
                    temp = root
                    root = root - func2grad
                    if (temp == root) return
                end if
                if (abs(func2grad) < abstol) return
                neval = neval + 1_IK
                func = getFunc(root, 0_IK)
                if (func < 0._RKG) then
                    lb = root
                else
                    ub = root
                end if
            end do
            ! Error occurred.
            root = .5_RKG * (lb + ub) ! for the sake of defining `root` on output (important for tests).
            neval = -neval
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setRoot_ENABLED && (Halley_ENABLED || Schroder_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        logical(LK) :: out_of_bounds_sentry, varset
        real(RKG) :: func, grad, hess, rold, fold, fmax, fmin, delta, delta1, delta2, temp, diff, ratio
        neval = 0_IK
        if (abs(lf) <= EPS10) then
            root = lb
        elseif (abs(uf) <= EPS10) then
            root = ub
        elseif (abs(ub - lb) < abstol) then
            root = .5_RKG * (ub + lb)
        else
            CHECK_BRACKET
            CHECK_LB_LT_UB
            CHECK_ZERO_LT_NITER
            CHECK_ZERO_LT_ABSTOL
            CHECK_LF(getFunc(lb, 0_IK))
            CHECK_UF(getFunc(ub, 0_IK))
            CHECK_ABSTOL_LT_LB_UB_DIFF
            !print *, "lb, root, ub", lb, root, ub
            if (root <= lb .or. ub <= root) root = .5_RKG * (lb + ub)
            out_of_bounds_sentry = .false._LK
            delta = huge(delta)
            fmax = 0._RKG
            fmin = 0._RKG
            rold = root
            func = 0._RKG
            delta1 = delta
            delta2 = delta
            neval = 0_IK
            do
                fold = func
                delta2 = delta1
                delta1 = delta
                neval = neval + 3_IK
                if (niter < neval) exit
                func = getFunc(root, 0_IK)
                grad = getFunc(root, 1_IK)
                hess = getFunc(root, 2_IK)
                if (func == 0._RKG) return
                if (grad == 0._RKG) then ! Handle zero derivative and hessian.
                    if (fold == 0._RKG) then
                        ! First iteration: pretend that we had a previous one at either `lb` or `ub`.
                        if (root == lb) then
                            rold = ub
                        else
                            rold = lb
                        end if
                        neval = neval + 1_IK
                        if (niter < neval) exit
                        fold = getFunc(rold, 0_IK)
                        delta = rold - root
                    end if
                    if ((fold < 0._RKG .and. 0._RKG < func) .or. (func < 0._RKG .and. 0._RKG < fold)) then
                        ! Crossed over. Move in opposite direction to last step.
                        if (delta < 0._RKG) then
                            delta = .5_RKG * (root - lb)
                        else
                            delta = .5_RKG * (root - ub)
                        end if
                    else ! Move in same direction as last step.
                        if (delta < 0._RKG) then
                            delta = .5_RKG * (root - ub)
                        else
                            delta = .5_RKG * (root - lb)
                        end if
                    end if
                elseif (hess /= 0._RKG) then
#if                 Halley_ENABLED || Schroder_ENABLED
                    ratio = 2 * func
                    temp = 2 * grad - func * (hess / grad)
                    varset = abs(temp) < 1._RKG
                    if (varset) then
                        varset = abs(temp) * huge(temp) <= abs(ratio)
                        if (varset) delta = func / grad ! Possible overflow, use Newton step.
                    end if
                    if (.not. varset) delta = ratio / temp ! Use Halley step.
#elif               Schroder_ENABLED
                    ! The Schroder stepper appears unstable for certain class of test problems where the Halley method does.
                    ! For now, Schroder is switched to Halley, but the origins of this instability must be explored.
                    ratio = func / grad
                    varset = root /= 0._RKG
                    if (varset) then
                        varset = abs(ratio / root) < 0.1_RKG
                        if (varset) then
                            delta = ratio + (hess / (2 * grad)) * ratio**2
                            ! Fall back to Newton iteration if hess contribution is larger than grad.
                            varset = 0._RKG <= delta * ratio
                        end if
                    end if
                    ! Fall back to Newton iteration.
                    if (.not. varset) delta = ratio
#else
#error              "Unrecognized interface."
#endif
                    if (delta * grad / func < 0._RKG) then
                        ! The Newton and Halley steps disagree about which way we should move,
                        ! likely due to cancelation error in the calculation of the Halley step,
                        ! or else the derivatives are so small that their values are basically trash.
                        ! Move in the direction indicated by a Newton step, but by no more than twice the
                        ! current rold value, otherwise the search can jump out of bounds.
                        delta = func / grad
                        if (2 * abs(rold) < abs(delta)) delta = sign(2._RKG, delta) * abs(rold)
                    end if
                else
                    delta = func / grad
                end if
                temp = abs(delta / delta2)
                if (0.8_RKG < temp .and. temp < 2._RKG) then
                    ! The last two steps did not converge.
                    if (0._RKG < delta) then
                        delta = .5_RKG * (root - lb)
                    else
                        delta = .5_RKG * (root - ub)
                    end if
                    if (root /= 0._RKG .and. root < abs(delta)) delta = sign(.9_RKG, delta) * abs(root) ! protect against huge jumps.
                    ! reset delta2 so that this branch will *not* be taken on the next iteration.
                    delta1 = delta * 3._RKG
                    delta2 = delta1
                end if
                rold = root
                root = root - delta
                ! Check for out of bounds step.
                if (root < lb) then
                    varset = abs(lb) < 1._RKG .and. 1._RKG < abs(root)
                    if (varset) then
                        varset = huge(root) / abs(root) < abs(lb)
                        if (varset) diff = 1000._RKG
                    end if
                    if (.not. varset) then
                        varset = abs(lb) < 1._RKG
                        if (varset) then
                            varset = huge(lb) * abs(lb) < abs(root)
                            if (varset) then
                                if (lb < 0._RKG .neqv. root < 0._RKG) then
                                    diff = -huge(diff)
                                else
                                    diff = +huge(diff)
                                end if
                            end if
                        end if
                        if (.not. varset) diff = root / lb
                    end if
                    if (abs(diff) < 1._RKG) diff = 1 / diff
                    if (0._RKG < diff .and. diff < 3._RKG .and. .not. out_of_bounds_sentry) then
                        ! Only a small out of bounds step, lets assume that the root is probably approximately at `lb`.
                        out_of_bounds_sentry = .true. ! only take this branch once!
                        delta = 0.99_RKG * (rold - lb)
                        root = rold - delta
                    else
                        root = .5_RKG * (lb + ub)
                        if (abs(lb - ub) < 2 * spacing(root)) return
                        call setBracketTowardMin(delta, getFunc, rold, func, lb, ub, neval, niter)
                        if (niter < neval) exit
                        root = rold - delta
                        rold = lb
                        cycle
                    end if
                elseif (ub < root) then
                    varset = abs(ub) < 1._RKG .and. 1._RKG < abs(root)
                    if (varset) then
                        varset = huge(root) / abs(root) < abs(ub)
                        if (varset) diff = 1000._RKG
                    end if
                    if (.not. varset) diff = root / ub
                    if (abs(diff) < 1._RKG) diff = 1._RKG / diff
                    if (0._RKG < diff .and. diff < 3._RKG .and. .not. out_of_bounds_sentry) then
                        ! Only a small out of bounds step, assume that the root is approximately at `lb`.
                        out_of_bounds_sentry = .true. ! Take this branch only once.
                        delta = 0.99_RKG * (rold - ub)
                        root = rold - delta
                    else
                        root = .5_RKG * (lb + ub)
                        if (abs(lb - ub) < 2 * spacing(root)) return
                        call setBracketTowardMax(delta, getFunc, rold, func, lb, ub, neval, niter)
                        if (niter < neval) exit
                        root = rold - delta
                        rold = lb
                        cycle
                    end if
                end if
                ! update brackets:
                if (0._RKG < delta) then
                    ub = rold
                    fmax = func
                else
                    lb = rold
                    fmin = func
                end if
                ! Sanity check that we bracket the rold:
                if (0._RKG < fmax * fmin) exit ! no root, possibly a minimum around root.
                if (abs(delta) < abs(root * abstol)) return
            end do
            ! Error occurred.
            root = .5_RKG * (lb + ub) ! for the sake of defining `root` on output (important for tests).
            neval = -neval
        end if

    contains

        subroutine setBracketTowardMax(delta, getFunc, rold, func, lb, ub, neval, niter)
            ! Move rold towards ub until we bracket the root, updating `lb` and `ub` on the fly.
            real(RKG), intent(inout) :: rold, lb, ub
            integer(IK), intent(inout) :: neval
            integer(IK), intent(in) :: niter
            real(RKG), intent(out) :: delta
            procedure(real(RKG)) :: getFunc
            real(RKG), intent(in) :: func
            real(RKG) :: guess0, funcurrent
            integer(IK) :: multiplier
            multiplier = 2_IK
            funcurrent = func
            guess0 = rold
            if (abs(lb) < abs(ub)) then
                loopSearch1: do
                    if (funcurrent < 0._RKG .neqv. func < 0._RKG) exit loopSearch1
                    lb = rold
                    rold = rold * multiplier
                    if (ub < rold) then
                        rold = ub
                        funcurrent = -funcurrent ! There must be a change of sign.
                        exit loopSearch1
                    end if
                    neval = neval + 1_IK
                    if (niter < neval) return
                    funcurrent = getFunc(rold, 0_IK)
                    multiplier = multiplier * 2_IK
                end do loopSearch1
            else ! If `lb` and `ub` are negative, divide to head towards `ub`.
                loopSearch2: do
                    if (funcurrent < 0._RKG .neqv. func < 0._RKG) exit loopSearch2
                    lb = rold
                    rold = rold / multiplier
                    if (ub < rold) then
                        rold = ub
                        funcurrent = -funcurrent ! There must be a change of sign.
                        exit loopSearch2
                    end if
                    neval = neval + 1_IK
                    if (niter < neval) return
                    funcurrent = getFunc(rold, 0_IK)
                    multiplier = multiplier * 2_IK
                end do loopSearch2
            end if
            if (neval < niter) then
                ub = rold
                if (16_IK < multiplier) then
                    call setBracketTowardMin(delta, getFunc, rold, funcurrent, lb, ub, neval, niter)
                    delta = delta + (guess0 - rold)
                    return
                end if
            end if
            delta = guess0 - .5_RKG * (lb + ub)
        end subroutine

        subroutine setBracketTowardMin(delta, getFunc, rold, func, lb, ub, neval, niter)
            ! Move rold towards lb until we bracket the root, updating `lb` and `ub` on the fly.
            real(RKG), intent(inout) :: rold, lb, ub
            integer(IK), intent(inout) :: neval
            integer(IK), intent(in) :: niter
            real(RKG), intent(out) :: delta
            procedure(real(RKG)) :: getFunc
            real(RKG), intent(in) :: func
            real(RKG) :: guess0, funcurrent
            integer(IK) :: multiplier
            multiplier = 2_IK
            funcurrent = func
            guess0 = rold
            if (abs(lb) < abs(ub)) then
                loopSearch1: do
                    if (funcurrent < 0._RKG .neqv. func < 0._RKG) exit loopSearch1
                    ub = rold
                    rold = rold / multiplier
                    if (rold < lb) then
                        rold = lb
                        funcurrent = -funcurrent ! There must be a change of sign.
                        exit loopSearch1
                    end if
                    neval = neval + 1_IK
                    if (niter < neval) return
                    funcurrent = getFunc(rold, 0_IK)
                    multiplier = multiplier * 2_IK
                end do loopSearch1
            else ! If `lb` and `ub` are negative, multiply to head towards `lb`.
                loopSearch2: do
                    if (funcurrent < 0._RKG .neqv. func < 0._RKG) exit loopSearch2
                    ub = rold
                    rold = rold * multiplier
                    if (rold < lb) then
                        rold = lb
                        funcurrent = -funcurrent ! There must be a change of sign.
                        exit loopSearch2
                    end if
                    neval = neval + 1_IK
                    if (niter < neval) return
                    funcurrent = getFunc(rold, 0_IK)
                    multiplier = multiplier * 2_IK
                end do loopSearch2
            end if
            if (neval < niter) then
                lb = rold
                if (16_IK < multiplier) then
                    call setBracketTowardMax(delta, getFunc, rold, funcurrent, lb, ub, neval, niter)
                    delta = delta + (guess0 - rold)
                    return
                end if
            end if
            delta = guess0 - .5_RKG * (lb + ub)
        end subroutine

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  CHECK_ABSTOL_LT_LB_UB_DIFF
#undef  CHECK_ZERO_LT_ABSTOL
#undef  CHECK_ZERO_LT_NITER
#undef  CHECK_LB_LT_UB
#undef  CHECK_BRACKET
#undef  CHECK_LF
#undef  CHECK_UF