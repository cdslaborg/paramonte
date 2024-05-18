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
!>  This include file contains procedure implementations of [pm_optimization](@ref pm_optimization).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SWAP(X,Y,TEMP) TEMP = X; X = Y; Y = TEMP

        !%%%%%%%%%%%%%%%%%%%
#if     isBracketMax_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        bracketed = xlow <= xmax .and. xmax <= xupp .and. flow <= fmax .and. fmax >= fupp

        !%%%%%%%%%%%%%%%%%%%
#elif   isBracketMin_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        bracketed = xlow <= xmin .and. xmin <= xupp .and. flow >= fmin .and. fmin <= fupp

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setBracketMax_ENABLED || setBracketMin_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     setBracketMax_ENABLED
#define COMPARES_WITH >
#define xmin xmax
#define fmin fmax
#elif   setBracketMin_ENABLED
#define COMPARES_WITH <
#else
#error  "Unrecognized interface."
#endif
        real(RKC)               :: lowf, uppf, fnew, xnew, ulim, q, r
        real(RKC)   , parameter :: GLIMIT = 100._RKC, GOLDEN = real(GOLDEN_RATIO, RKC)
        real(RKC)   , parameter :: SMALL = sqrt(epsilon(1._RKC))**3 ! Protects against trying to achieve fractional accuracy for a minimum whose value is exactly zero.
        CHECK_ASSERTION(__LINE__, xlow <= xupp, SK_"@setBracketMin(): The condition `xlow < xupp` must hold. xlow, xupp = "//getStr([xlow, xupp]))
        ! Swap `xlow` and `xupp` if needed, to ensure downhill direction from `xlow` to xupp`.
        lowf = getFunc(xlow)
        uppf = getFunc(xupp)
        if (lowf COMPARES_WITH uppf) then
            ulim = xlow
            xlow = xupp
            xupp = ulim
            ulim = lowf
            lowf = uppf
            uppf = ulim
        end if
        xmin = xupp + GOLDEN * (xupp - xlow)
        fmin = getFunc(xmin)
        loopRefine: do niter = 1, niter
            if (uppf COMPARES_WITH fmin) exit loopRefine
            r = (xupp - xlow) * (uppf - fmin)
            q = (xupp - xmin) * (uppf - lowf)
            xnew = xupp - ((xupp - xmin) * q - (xupp - xlow) * r) / (2 * sign(max(abs(q - r), SMALL), q - r))
            ulim = xupp + GLIMIT * (xmin - xupp)
            if (0._RKC < (xupp - xnew) * (xnew - xmin)) then
                fnew = getFunc(xnew)
                if (fnew COMPARES_WITH fmin) then
                    xlow = xupp
                    lowf = uppf
                    xupp = xnew
                    uppf = fnew
                    exit loopRefine
                else if (uppf COMPARES_WITH fnew) then
                    xmin = xnew
                    fmin = fnew
                    exit loopRefine
                end if
                xnew = xmin + GOLDEN * (xmin - xupp)
                fnew = getFunc(xnew)
            else if (0._RKC < (xmin - xnew) * (xnew - ulim)) then
                fnew = getFunc(xnew)
                if (fnew COMPARES_WITH fmin) then
                    xupp = xmin
                    xmin = xnew
                    xnew = xmin + GOLDEN * (xmin - xupp)
                    ! shift.
                    uppf = fmin
                    fmin = fnew
                    fnew = getFunc(xnew)
                end if
            else if ((xnew - ulim) * (ulim - xmin) < 0._RKC) then
                xnew = xmin + GOLDEN * (xmin - xupp)
                fnew = getFunc(xnew)
            else
                xnew = ulim
                fnew = getFunc(xnew)
            end if
            ! shift.
            xlow = xupp
            xupp = xmin
            xmin = xnew
            ! shift.
            lowf = uppf
            uppf = fmin
            fmin = fnew
        end do loopRefine
        ! final swap.
        SWAP(xmin,xupp,ulim)
        SWAP(fmin,uppf,ulim)
        if (xupp < xlow) then
            SWAP(xlow,xupp,ulim)
        end if
        if (present(flow)) flow = lowf
        if (present(fupp)) fupp = uppf
#undef  xmin
#undef  fmin

        !%%%%%%%%%%%%%%%%%%
#elif   getMinBrent_ENABLED
        !%%%%%%%%%%%%%%%%%%

        integer(IK) :: retin, niter_init
        real(RKC) :: lowx, uppx, minf, lot
        integer(IK), parameter :: MAXITER = int(100 * precision(xmin) / 53._RKC, IK)
        if (present(xlow)) then
            lowx = xlow
        else
            lowx = 0.1_RKC
        end if
        if (present(xupp)) then
            uppx = xupp
        else
            uppx = 0.9_RKC
        end if
        if (present(tol)) then
            lot = tol
        else
            lot = sqrt(epsilon(xmin))
        end if
        if (present(niter)) then
            niter_init = niter
        else
            niter_init = MAXITER
        end if
        if (present(xlow) .and. present(xupp)) then
            xmin = .5_RKC * (xupp - xlow)
            minf = getFunc(xmin)
        else
            retin = niter_init
            call setBracketMin(getFunc, retin, xmin, lowx, uppx, minf)
            if (niter_init < retin) then
                if (present(niter)) then
                    niter = retin
                    return
                else
                    error stop "@getMinBrent(): Bracketing failed. The specified bracket `[xlow, xupp]` may be too wide given the specified `niter` or the input function is concave without minimum."
                end if
            end if
        end if
        retin = niter_init
        call setMinBrent(getFunc, xmin, lowx, uppx, minf, lot, retin)
        if (niter_init < retin) then
            if (present(niter)) then
                niter = retin
                return
            else
                error stop "@getMinBrent(): The Brent minimizer failed to converge. Specifying a larger value for `niter` might resolve the error."
            end if
        end if
        if (present(fmin)) fmin = minf

        !%%%%%%%%%%%%%%%%%%
#elif   setMinBrent_ENABLED
        !%%%%%%%%%%%%%%%%%%

        real(RKC)   , parameter :: SMALL = sqrt(epsilon(1._RKC))**3 ! Protects against trying to achieve fractional accuracy for a minimum whose value is exactly zero.
        real(RKC)   , parameter :: GOLEN_MEAN = .5_RKC * (3._RKC - sqrt(5._RKC)) ! Golden Section switch criterion
        real(RKC)               :: xz, xw, xold, fold, fv, fw, p, q, r, d, etemp, relerr, abstol, tolby2, center

        CHECK_ASSERTION(__LINE__, 0 < tol, SK_"@setMinBrent(): The condition `0 < tol` must hold. tol = "//getStr(tol))
        CHECK_ASSERTION(__LINE__, 0 < niter, SK_"@setMinBrent(): The condition `0 < niter` must hold. niter = "//getStr(niter))
        CHECK_ASSERTION(__LINE__, xlow <= xmin .and. xmin <= xupp, SK_"@setMinBrent(): The condition `xlow < xmin .and. xmin < xupp` must hold. xlow, xmin, xupp = "//getStr([xlow, xmin, xupp]))
        CHECK_ASSERTION(__LINE__, getFunc(xlow) >= fmin .or. fmin <= getFunc(xupp), SK_"@setMinBrent(): The condition `getFunc(xlow) >= fmin .or. fmin =< getFunc(xupp)` must hold. xlow, xmin, xupp, flow, fmin, fupp = "//getStr([xlow, xmin, xupp, getFunc(xlow), fmin, getFunc(xupp)]))

        xw = xmin
        xmin = xmin
        xold = xmin
        relerr = 0._RKC
        !fmin = getFunc(xmin)
        fv = fmin
        fw = fmin
        do niter = 1, niter
            center = .5_RKC * (xlow + xupp)
            abstol = tol * abs(xmin) + SMALL
            tolby2 = 2._RKC * abstol
            if (abs(xmin - center) <= (tolby2 - 0.5_RKC * (xupp - xlow))) return
            if (abstol < abs(relerr)) then
                r = (xmin - xw) * (fmin - fv)
                q = (xmin - xold) * (fmin - fw)
                p = (xmin - xold) * q - (xmin - xw) * r
                q = 2.0_RKC * (q - r)
                if (0._RKC < q) p =  - p
                q = abs(q)
                etemp = relerr
                relerr = d
                if (abs(0.5_RKC * q * etemp) <= abs(p) .or. p <= q * (xlow - xmin) .or. q * (xupp - xmin) <= p) then
                    if (xmin < center) then
                        relerr = xupp - xmin
                    else
                        relerr = xlow - xmin
                    end if
                    d = GOLEN_MEAN * relerr
                else
                    d = p / q
                    xz = xmin + d
                    if (xz - xlow < tolby2 .or. xupp - xz < tolby2) d = sign(abstol, center - xmin)
                end if
            else
                if (xmin < center) then
                    relerr = xupp - xmin
                else
                    relerr = xlow - xmin
                end if
                d = GOLEN_MEAN * relerr
            end if
            if (abs(d) < abstol) then
                xz = xmin + sign(abstol, d)
            else
                xz = xmin + d
            end if
            fold = getFunc(xz)
            if (fmin < fold) then
                if (xz < xmin) then
                    xlow = xz
                else
                    xupp = xz
                end if
                if (fold <= fw .or. xw == xmin) then
                    xold = xw
                    fv = fw
                    xw = xz
                    fw = fold
                else if (fold <= fv .or. xold == xmin .or. xold == xw) then
                    xold = xz
                    fv = fold
                end if
            else
                if (xz < xmin) then
                    xupp = xmin
                else
                    xlow = xmin
                end if
                ! shift.
                xold = xw
                xw = xmin
                xmin = xz
                ! shift.
                fv = fw
                fw = fmin
                fmin = fold
            end if
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   isFailedMinPowell_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: retin
        real(RKC) :: dirset(size(xmin, 1, IK), size(xmin, 1, IK)), lot, minf
        integer(IK) , parameter :: MAXITER = int(100 * precision(xmin) / 53._RKC, IK)
        call setMatInit(dirset, uppLowDia, vupp = 0._RKC, vlow = 0._RKC, vdia = 1._RKC, nrow = size(xmin, 1, IK), ncol = size(xmin, 1, IK), roff = 0_IK, coff = 0_IK)
        if (present(tol)) then
            lot = tol
        else
            lot = sqrt(epsilon(xmin))
        end if
        if (present(niter)) then
            retin = niter
        else
            retin = MAXITER
        end if
        if (present(fmin)) then
            minf = fmin
        else
            minf = getFunc(xmin)
        end if
        call setMinPowell(getFunc, xmin, minf, dirset, lot, retin)
        if (present(fmin)) fmin = minf
        failed = MAXITER < retin

        !%%%%%%%%%%%%%%%%%%%
#elif   setMinPowell_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        integer(IK)             :: ndim, idim, ibig, retin
        real(RKC)               :: diff, fold, fnew, t, xold(size(xmin, 1, IK)), xnew(size(xmin, 1, IK)), xdir(size(xmin, 1, IK)), xtmp(size(xmin, 1, IK))
        real(RKC)   , parameter :: SMALL = sqrt(epsilon(1._RKC))**3 ! Protects against trying to achieve fractional accuracy for a minimum whose value is exactly zero.
        ndim = size(xmin, 1, IK)
        CHECK_ASSERTION(__LINE__, 0 < tol, SK_"@setMinPowell(): The condition `0 < tol` must hold. tol = "//getStr(tol))
        CHECK_ASSERTION(__LINE__, 0 < niter, SK_"@setMinPowell(): The condition `0 < niter` must hold. niter = "//getStr(niter))
        CHECK_ASSERTION(__LINE__, all(ndim == shape(dirset, IK)), SK_"@setMinPowell(): The condition `all(size(xmin) == shape(dirset))` must hold. size(xmin), shape(dirset) = "//getStr([size(xmin, 1, IK), shape(dirset, IK)]))
        CHECK_ASSERTION(__LINE__, abs(fmin - getFunc(xmin)) < 100 * epsilon(0._RKC), SK_"@setMinPowell(): The condition `fmin - getFunc(xmin)` must hold. fmin, getFunc(xmin) = "//getStr([fmin, getFunc(xmin)]))
        retin = niter
        xold = xmin
        do
            fold = fmin
            ibig = 0_IK
            diff = 0._RKC
            do idim = 1, ndim
                xdir = dirset(1 : ndim, idim)
                fnew = fmin
                call setMinDir() ! minimize along the specified direction.
                if (retin < niter) return ! error occurred.
                if (diff < fnew - fmin) then
                    diff = fnew - fmin
                    ibig = idim
                end if
            end do
            if (2 * (fold - fmin) <= tol * (abs(fold) + abs(fmin)) + SMALL) return
            !if (niter < iter) return
            xnew = 2 * xmin - xold
            xdir = xmin - xold
            xold = xmin
            fnew = getFunc(xnew)
            CHECK_ASSERTION(__LINE__, ibig /= 0_IK, SK_"@setMinPowell(): Internal error occurred: ibig = 0")
            if (fold < fnew) then
                t = 2 * (fold - 2 * fmin + fnew) * (fold - fmin - diff)**2 - diff * (fold - fnew)**2
                if (t < 0._RKC) then
                    call setMinDir() ! minimize along the specified direction.
                    if (retin < niter) return ! error occurred.
                    dirset(1 : ndim, ibig) = dirset(1 : ndim, ndim)
                    dirset(1 : ndim, ndim) = xdir
                end if
            end if
        end do

    contains

        function getFuncDir(x) result(func)
            real(RKC), intent(in)    :: x
            real(RKC)                :: func
            xtmp = xmin + x * xdir
            func = getFunc(xtmp)
        end function

        subroutine setMinDir()
            integer(IK) :: jdim
            real(RKC) :: minx, xlow, xupp
            ! \todo: There must be a better initial bracketing guess at each iteration.
            xlow = 0._RKC
            xupp = 1._RKC
            niter = 1000
            call setBracketMin(getFuncDir, niter, minx, xlow, xupp, fmin)
            if (1000 < niter) then
                niter = retin + 1
            else
                call setMinBrent(getFuncDir, minx, xlow, xupp, fmin, tol, niter)
                do concurrent(jdim = 1 : ndim)
                    xdir(jdim) = xdir(jdim) * minx
                    xmin(jdim) = xmin(jdim) + xdir(jdim)
                end do
            end if
        end subroutine

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  COMPARES_WITH
#undef  SWAP