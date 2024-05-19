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
!>  This include file contains the implementations of procedures of [pm_quadPack](@ref pm_quadPack).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the help argument.
#if     QAGD_ENABLED
#define HELP_ARG
#elif   QAGP_ENABLED && (IF_ENABLED || FI_ENABLED || II_ENABLED)
#define HELP_ARG, breakTrans
#elif   QAGS_ENABLED || QAGP_ENABLED || QAWC_ENABLED || QAWFS_ENABLED || QAWFC_ENABLED
#define HELP_ARG, help
#elif   !(wcauchy_typer_ENABLED || wsin_typer_ENABLED || wcos_typer_ENABLED || getQuadGK_ENABLED || setNodeWeightGK_ENABLED || setErrSorted_ENABLED || setSeqLimEps_ENABLED || setChebExpan_ENABLED)
#error  "Unrecognized interface."
#endif
        ! Define the GK quadrature arguments.
#if     GKXX_ENABLED
#define QRULE_ARG nodeK, weightK, weightG
#elif   GK15_ENABLED || GK21_ENABLED || GK31_ENABLED || GK41_ENABLED || GK51_ENABLED || GK61_ENABLED
#define QRULE_ARG qrule
#elif   !(wcauchy_typer_ENABLED || wsin_typer_ENABLED || wcos_typer_ENABLED || setNodeWeightGK_ENABLED || setErrSorted_ENABLED || setSeqLimEps_ENABLED || setChebExpan_ENABLED || isFailedQuad_ENABLED)
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%
#if     wcauchy_typer_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        wcauchy%cs = real(cs, RKG)

        !%%%%%%%%%%%%%%%%%
#elif   wsin_typer_ENABLED
        !%%%%%%%%%%%%%%%%%

        wsin%omega = real(omega, RKG)

        !%%%%%%%%%%%%%%%%%
#elif   wcos_typer_ENABLED
        !%%%%%%%%%%%%%%%%%

        wcos%omega = real(omega, RKG)

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   setNodeWeightGK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        logical(LK)             :: orderIsEven
        integer(IK)             :: i, k, kk, l, ll, halfOrder, halfOrderPlusOne
        real(RKG)               :: ak, an, bb, c, coefNode, coefWeight, d, x1, xx, s, y
        real(RKG)   , parameter :: TOL = epsilon(0._RKG) * 10._RKG
#if     Fixed_ENABLED
        real(RKG)               :: B((size(nodeK, kind = IK) / 2_IK) + 1_IK), Tau(size(nodeK, kind = IK) / 2_IK)
        integer(IK)             :: order
        order = size(nodeK, kind = IK) - 1_IK
#elif   Alloc_ENABLED
        real(RKG)               :: B(((order + 1_IK) / 2_IK) + 1_IK), Tau((order + 1_IK) / 2_IK)
        allocate(nodeK(order + 1_IK), weightK(order + 1_IK), weightG((order + 1_IK) / 2_IK))
#else
#error  "Unrecognized intefrace."
#endif
        halfOrder = size(nodeK, kind = IK) / 2_IK
        halfOrderPlusOne = halfOrder + 1_IK
        orderIsEven = (2_IK * halfOrder == order)

        CHECK_ASSERTION(__LINE__, 0_IK <= order, SK_"@setNodeWeightGK(): The condition `0 <= order` must hold. order = "//getStr(order))
        CHECK_ASSERTION(__LINE__, size(weightK, kind = IK) == order + 1_IK, SK_"@setNodeWeightGK(): The condition `size(weightK) == order + 1` must hold. size(weightK), order = "//getStr([size(weightK, kind = IK), order]))
        CHECK_ASSERTION(__LINE__, size(weightG, kind = IK) == halfOrder, SK_"@setNodeWeightGK(): The condition `size(weightG) == halfOrder` must hold. size(weightG), halfOrder = "//getStr([size(weightG, kind = IK), halfOrder]))

        d = 2._RKG
        an = 0._RKG
        do k = 1_IK, order
          an = an + 1._RKG
          d = d * an / (an + 0.5_RKG)
        end do

        ! Compute the Chebyshev coefficients of the orthogonal polynomial.

        Tau(1) = (an + 2._RKG) / (an + an + 3._RKG)
        B(halfOrder) = Tau(1) - 1._RKG
        ak = an
        do l = 1_IK, halfOrder - 1_IK
          ak = ak + 2._RKG
          Tau(l + 1_IK) = ((ak - 1._RKG) * ak - an * (an + 1._RKG)) * (ak + 2._RKG) * Tau(l) / (ak * ((ak + 3._RKG) * (ak + 2._RKG) - an * (an + 1._RKG)))
          B(halfOrder - l) = Tau(l + 1_IK)
          do ll = 1_IK, l
            B(halfOrder - l) = B(halfOrder - l) + Tau(ll) * B(halfOrder - l + ll)
          end do
        end do
        B(halfOrderPlusOne) = 1._RKG

        ! Calculation of approximate values for the abscissas.

        bb = sin(1.570796_RKG / (an + an + 1._RKG))
        x1 = sqrt(1._RKG - bb * bb)
        s = 2._RKG * bb * x1
        c = sqrt(1._RKG - s * s)
        coefNode = 1._RKG - (1._RKG - 1._RKG / an) / (8._RKG * an * an)
        xx = coefNode * x1

        ! Coefficient needed for weights.

        ! COEF2 = 2^(2*order+1) * order! * order! / (2n+1)!

        coefWeight = 2._RKG / real(2_IK * order + 1_IK, RKG)
        do i = 1_IK, order
          coefWeight = coefWeight * 4._RKG * real(i, RKG) / real(order + i, RKG)
        end do

        ! Calculation of the K-th abscissa (a Kronrod abscissa) and the corresponding weight.

        !do k = 1_IK, order, 2_IK
        do k = 1_IK, (order + 1_IK) / 2_IK

            kk = 2_IK * k - 1_IK
            call setNodeWeightK(xx, weightK(kk))
            !weightG(k) = 0._RKG
            nodeK(kk) = xx
            y = x1
            x1 = y * c - bb * s
            bb = y * s + bb * c
            if (kk == order) then
              xx = 0._RKG
            else
              xx = coefNode * x1
            end if

            ! Calculation of the K+1 abscissa (a Gaussian abscissa) and the corresponding weights.

            call setNodeWeightGK(xx, weightK(kk + 1_IK), weightG(k))

            nodeK(kk + 1_IK) = xx
            y = x1
            x1 = y * c - bb * s
            bb = y * s + bb * c
            xx = coefNode * x1

        end do

        ! Compute the Kronrod abscissa at the origin if order is even.

        if (orderIsEven) then
            xx = 0._RKG
            call setNodeWeightK(xx, weightK(order + 1_IK))
            !weightG(order + 1_IK) = 0._RKG
            nodeK(order + 1_IK) = xx
        end if

    contains

        pure subroutine setNodeWeightK(node, weight)

            real(RKG)   , intent(inout) :: node
            real(RKG)   , intent(out)   :: weight
            real(RKG)   :: ai, b0, ub1, ub2, d0, d1, d2, delta, dif, getFunc, fd, yy
            logical(LK) :: convergenceOccurred
            integer(IK) :: ii, iter, kk

            convergenceOccurred = node == 0._RKG

            ! Iteratively compute a Kronrod node.

            do iter = 1_IK, 50_IK

                ub1 = 0._RKG
                ub2 = B(halfOrderPlusOne)
                yy = 4._RKG * node * node - 2._RKG
                d1 = 0._RKG

                if (orderIsEven) then
                    ai = halfOrder + halfOrderPlusOne
                    d2 = ai * B(halfOrderPlusOne)
                    dif = 2._RKG
                else
                    ai = real(halfOrderPlusOne, RKG)
                    d2 = 0._RKG
                    dif = 1._RKG
                end if

                do kk = 1_IK, halfOrder
                    ai = ai - dif
                    ii = halfOrderPlusOne - kk
                    b0 = ub1
                    ub1 = ub2
                    d0 = d1
                    d1 = d2
                    ub2 = yy * ub1 - b0 + B(ii)
                    if (.not. orderIsEven) ii = ii + 1_IK
                    d2 = yy * d1 - d0 + ai * B(ii)
                end do

                if (orderIsEven) then
                    getFunc = node * (ub2 - ub1)
                    fd = d2 + d1
                else
                    getFunc = 0.5_RKG * (ub2 - b0)
                    fd = 4._RKG * node * d2
                end if

                ! Newton correction.

                delta = getFunc / fd
                node = node - delta
                if (convergenceOccurred) then
                    d0 = 1._RKG
                    d1 = node
                    ai = 0._RKG
                    do kk = 2_IK, order
                        ai = ai + 1._RKG
                        d2 = ((ai + ai + 1._RKG) * node * d1 - ai * d0) / (ai + 1._RKG)
                        d0 = d1
                        d1 = d2
                    end do
                    weight = coefWeight / (fd * d2)
                    return
                end if
                convergenceOccurred = logical(abs(delta) <= TOL, LK)

            end do

            error stop "Iterative computation of the Kronrod node and weight failed."

        end subroutine

        pure subroutine setNodeWeightGK(node, weightKronrod, weightGauss)

            real(RKG)   , intent(inout) :: node
            real(RKG)   , intent(out)   :: weightKronrod, weightGauss
            real(RKG)   :: ai, delta, p0, p1, p2, pd0, pd1, pd2, yy
            logical(LK) :: convergenceOccurred
            integer(IK) :: ii, iter, kk
            convergenceOccurred = logical(node == 0._RKG, LK)

            ! Iteratively compute a Kronrod node.

            do iter = 1_IK, 50_IK

                p1 = node
                p0 = 1._RKG
                pd0 = 0._RKG
                pd1 = 1._RKG

                if (order <= 1) then ! Initialize P2 and PD2 to avoid problems with `delta`.
                    if (epsilon(node) < abs(node)) then
                        p2 = 0.5_RKG * (3._RKG * node * node - 1._RKG)
                        pd2 = 3._RKG * node
                    else
                        p2 = 3._RKG * node
                        pd2 = 3._RKG
                    end if
                end if

                ai = 0._RKG
                do kk = 2_IK, order
                    ai = ai + 1._RKG
                    p2 = ((ai + ai + 1._RKG) * node * p1 - ai * p0) / (ai + 1._RKG)
                    pd2 = ((ai + ai + 1._RKG) * (p1 + node * pd1) - ai * pd0) / (ai + 1._RKG)
                    p0 = p1
                    p1 = p2
                    pd0 = pd1
                    pd1 = pd2
                end do

                ! Newton correction.

                delta = p2 / pd2
                node = node - delta
                if (convergenceOccurred) then
                    weightGauss = 2._RKG / (order * pd2 * p0)
                    yy = 4._RKG * node * node - 2._RKG
                    p2 = B(halfOrderPlusOne)
                    p1 = 0._RKG
                    do kk = 1_IK, halfOrder
                        ii = halfOrderPlusOne - kk
                        p0 = p1
                        p1 = p2
                        p2 = yy * p1 - p0 + B(ii)
                    end do
                    if (orderIsEven) then
                        weightKronrod = weightGauss + coefWeight / (pd2 * node * (p2 - p1))
                    else
                        weightKronrod = weightGauss + 2._RKG * coefWeight / (pd2 * (p2 - p0))
                    end if
                    return
                end if
                convergenceOccurred = logical(abs(delta) <= TOL, LK)

            end do

            error stop "Iterative computation of Gauss-Kronrod nodes and weights failed."

        end subroutine

        !%%%%%%%%%%%%%%%%%%%
#elif   setErrSorted_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        ! set maxErrLoc and maxErrVal.
#define SET_ERROR_VALUE maxErrLoc = sindex(nrmax); maxErrVal = sinfoErr(maxErrLoc);

        real(RKG)   :: errmin
        integer(IK) :: i, ibeg, ido, isucc, j, jbnd, jupbn, k
        integer(IK) :: nint

        nint = size(sindex, kind = IK)
        CHECK_ASSERTION(__LINE__, nintmax >= nint, SK_"@setErrSorted(): The condition `nintmax >= nint` must hold. nintmax, nint = "//getStr([nintmax, nint]))
        CHECK_ASSERTION(__LINE__, size(sinfoErr, kind = IK) == nint, SK_"@setErrSorted(): The condition `size(sinfoErr, kind = IK) == nint` must hold. size(sinfoErr, kind = IK), nint = "//getStr([size(sinfoErr, kind = IK), nint]))

        ! Check if the list contains more than two error estimates.

        if (nint > 2_IK) then

            ! \remark
            ! This part of the routine is only executed if, due to a difficult integrand, subdivision increased the error estimate.
            ! The insert procedure should start after the nrmax-th largest error estimate, in the normal case.

            maxErrVal = sinfoErr(maxErrLoc)
            if (nrmax /= 1_IK) then
                ido = nrmax - 1_IK
                do i = 1_IK, ido
                    isucc = sindex(nrmax - 1_IK)
                    if (maxErrVal <= sinfoErr(isucc)) exit
                    sindex(nrmax) = isucc
                    nrmax = nrmax - 1_IK
                end do
            end if

            ! Compute the number of elements in the list to be maintained in descending order.
            ! This number depends on the number of subdivisions still allowed.

            jupbn = nint
            if (nint > (nintmax / 2_IK + 2_IK)) jupbn = nintmax + 3_IK - nint
            errmin = sinfoErr(nint)

            ! Insert maxErrVal by traversing the list top-down, starting comparison from the element sinfoErr(sindex(nrmax + 1_IK)).

            jbnd = jupbn - 1_IK
            ibeg = nrmax + 1_IK
            if (ibeg <= jbnd) then
                do i = ibeg, jbnd
                    isucc = sindex(i)
                    if (maxErrVal >= sinfoErr(isucc)) then
                        ! insert errmin by traversing the list bottom-up.
                        sindex(i - 1) = maxErrLoc
                        k = jbnd
                        do j = i, jbnd
                            isucc = sindex(k)
                            if (errmin < sinfoErr(isucc)) then
                                sindex(k + 1) = nint
                                SET_ERROR_VALUE ! fpp
                                return
                            end if
                            sindex(k + 1) = isucc
                            k = k - 1_IK
                        end do
                        sindex(i) = nint
                        SET_ERROR_VALUE ! fpp
                        return
                    end if
                    sindex(i - 1) = isucc
                end do
            end if
            sindex(jbnd) = maxErrLoc
            sindex(jupbn) = nint
        else
            sindex(1) = 1
            sindex(2) = 2
        end if
        SET_ERROR_VALUE ! fpp
#undef  SET_ERROR_VALUE

        !%%%%%%%%%%%%%%%%%%%
#elif   setSeqLimEps_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        real(RKG)   , parameter :: HUGE_RKG = huge(0._RKG)
        real(RKG)   , parameter :: EPS_RKG = epsilon(0._RKG)
        real(RKG)               :: delta1, delta2, delta3, epsinf, err1, err2, err3, e0, e1, e1abs, e2, e3, res, ss, tol1, tol2, tol3
        real(RKG)               :: error ! `error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)`
        integer(IK)             :: newelm ! number of elements to be computed in the new diagonal
        integer(IK)             :: i, ib, ib2, indx, k1, k2, k3, num
        ! `seqlim` is the element in the new diagonal with least value of error.
        ical = ical + 1_IK
        abserr = HUGE_RKG
        seqlim = EpsTable(inew)
        blockNewCondition: if (inew >= 3_IK) then
            EpsTable(inew + 2) = EpsTable(inew)
            newelm = (inew - 1_IK) / 2_IK
            EpsTable(inew) = HUGE_RKG
            num = inew
            k1 = inew
            loopNewElement: do i = 1_IK, newelm
                k2 = k1 - 1_IK
                k3 = k1 - 2_IK
                res = EpsTable(k1 + 2_IK)
                e0 = EpsTable(k3)
                e1 = EpsTable(k2)
                e2 = res
                e1abs = abs(e1)
                delta2 = e2 - e1
                err2 = abs(delta2)
                tol2 = max(abs(e2), e1abs) * EPS_RKG
                delta3 = e1 - e0
                err3 = abs(delta3)
                tol3 = max(e1abs, abs(e0)) * EPS_RKG
                if (err2 > tol2 .or. err3 > tol3) then
                    e3 = EpsTable(k1)
                    EpsTable(k1) = e1
                    delta1 = e1 - e3
                    err1 = abs(delta1)
                    tol1 = max(e1abs, abs(e3)) * EPS_RKG
                    ! If two elements are very close to each other, omit a part of the table by adjusting the value of `inew`.
                    if (err1 > tol1 .and. err2 > tol2 .and. err3 > tol3) then
                        ss = 1._RKG / delta1 + 1._RKG / delta2 - 1._RKG / delta3
                        epsinf = abs(ss * e1)
                        ! Test to detect irregular behavior in the table and eventually omit a part of the table adjusting the value of `inew`.
                        if (epsinf > 1.e-4_RKG) then
                            res = e1 + 1._RKG / ss
                            EpsTable(k1) = res
                            k1 = k1 - 2_IK
                            error = err2 + abs(res - e2) + err3
                            if (error <= abserr) then
                                abserr = error
                                seqlim = res
                            end if
                            cycle loopNewElement
                        end if
                    end if
                    inew = i + i - 1_IK
                    exit loopNewElement
                else ! e0, e1 and e2 are equal to within machine accuracy. Convergence is assumed.
                    ! abserr = abs(e1 - e0) + abs(e2 - e1)
                    ! seqlim = e2
                    seqlim = res
                    abserr = err2 + err3
                    exit blockNewCondition
                end if
            end do loopNewElement
            ! shift the table.
            !if (inew > MAXLEN_EPSTAB) error stop "inew > MAXLEN_EPSTAB cannot happen."
            if (inew == MAXLEN_EPSTAB) inew = 2_IK * (MAXLEN_EPSTAB / 2_IK) - 1_IK
            ib = 1_IK; if ((num / 2_IK) * 2_IK == num) ib = 2_IK
            do i = 1_IK, newelm + 1_IK
                ib2 = ib + 2_IK
                EpsTable(ib) = EpsTable(ib2)
                ib = ib2
            end do
            if (num /= inew) then
                indx = num - inew + 1_IK
                do i = 1_IK, inew
                    EpsTable(i) = EpsTable(indx)
                    indx = indx + 1_IK
                end do
            end if
            if (ical < 4_IK) then
                seqlims(ical) = seqlim
                abserr = HUGE_RKG
            else
                abserr = sum(abs(seqlim - seqlims)) ! abs(seqlim - seqlims(3)) + abs(seqlim - seqlims(2)) + abs(seqlim - seqlims(1))
                seqlims(1) = seqlims(2)
                seqlims(2) = seqlims(3)
                seqlims(3) = seqlim
            end if
        end if blockNewCondition
        abserr = max(abserr, 5._RKG * EPS_RKG * abs(seqlim)) ! \warning This is set incorrectly in QuadPack of John Burkardt.

        !%%%%%%%%%%%%%%%%%%%
#elif   isFailedQuad_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        real(RKG), parameter :: vinf = huge(0._RKG) / 1000._RKG ! virtual infinity.
        integer(IK) , parameter :: MAX_NUM_INTERVAL = 2000_IK
        integer(IK) :: err, nint_def, neval_def, sindex(MAX_NUM_INTERVAL)
        real(RKG) :: abserr_def, abstol_def, reltol_def
        real(RKG) :: sinfo(4, MAX_NUM_INTERVAL)
        abstol_def = 0._RKG; if (present(abstol)) abstol_def = abstol
        reltol_def = epsilon(0._RKG)**(2./3.); if (present(reltol)) reltol_def = reltol
        if (-vinf < lb .and. ub < vinf)  then ! finite
            err = getQuadErr(getFunc, lb, ub, abstol_def, reltol_def, GK21 HELP_ARG, integral, abserr_def, sinfo, sindex, neval_def, nint_def)
        elseif (lb <= -vinf .and. vinf <= ub) then ! both infinite
            err = getQuadErr(getFunc, ninf, pinf, abstol_def, reltol_def, GK21 HELP_ARG, integral, abserr_def, sinfo, sindex, neval_def, nint_def)
        elseif (vinf <= ub)  then ! positive semi-infinite
            err = getQuadErr(getFunc, lb, pinf, abstol_def, reltol_def, GK21 HELP_ARG, integral, abserr_def, sinfo, sindex, neval_def, nint_def)
        elseif (lb <= -vinf)  then ! negative semi-infinite
            err = getQuadErr(getFunc, ninf, ub, abstol_def, reltol_def, GK21 HELP_ARG, integral, abserr_def, sinfo, sindex, neval_def, nint_def)
        else ! possibly NAN
            failed = .true._LK
            if (present(msg)) msg = SK_"The specified integration range is invalid. lb, ub = "//getStr([lb, ub])
            return
        end if
        failed = logical(err /= 0_IK, kind = LK)
        if (present(abserr)) abserr = abserr_def
        if (present(neval)) neval = neval_def
        if (present(nint)) nint = nint_def
        if (failed) then
            if (present(msg)) then
                if (err == 1_IK) then
                    msg = SK_"The maximum allowed number of subdivisions was reached without convergence."
                elseif (err == 2_IK) then
                    msg = SK_"Roundoff error was detected, which prevents the requested tolerance from being achieved."
                elseif (err == 3_IK) then
                    msg = SK_"Roundoff error was detected, which prevents the requested tolerance from being achieved."
                elseif (err == 4_IK) then
                    msg = SK_"An extremely bad integrand behavior occurs at some points of the integration interval."
                elseif (err == 5_IK) then
                    msg = SK_"The algorithm did not converge because the integral is likely divergent or slowly convergent."
                else
                    error stop SK_"Internal error occurred. Invalid error code detected: "//getStr(err)
                end if
            end if
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getQuadErr_ENABLED && (QAGD_ENABLED || QAGS_ENABLED || QAGP_ENABLED || QAWC_ENABLED) && (FI_ENABLED || IF_ENABLED || II_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Setup the transformation coefficient.
#if     QAWC_ENABLED
        real(RKG) :: cstrans
#if     II_ENABLED
#define TRANS_COEF invx * (-cstrans)
        cstrans = 1._RKG / (1._RKG + real(help%cs, RKG))
#elif   IF_ENABLED
#define TRANS_COEF invx * cstrans
        cstrans = 1._RKG / (1._RKG - real(help%cs, RKG) + ub)
#elif   FI_ENABLED
#define TRANS_COEF invx * (-cstrans)
        cstrans = 1._RKG / (1._RKG + real(help%cs, RKG) - lb)
#endif
#elif   QAGD_ENABLED || QAGS_ENABLED || QAGP_ENABLED
#define TRANS_COEF invx**2
#else
#error  "Unrecognized interface."
#endif
        ! Transform QAGP break points.
#if     QAGP_ENABLED && II_ENABLED
        real(RKG) :: breakTrans(size(help, kind = IK))
        ! Double infinite range requires special care to remove duplicate breaks.
        integer(IK) :: locLastNegBreak, lenBreak
        CHECK_ASSERTION(__LINE__, all(help /= -1._RKG), SK_"@getQuadErr(): The condition `all(help /= -1.)` must hold. help = "//getStr(help))
        lenBreak = size(help, kind = IK)
        locLastNegBreak = maxloc(help, mask = help < -1._RKG, dim = 1)
        call setMerged(breakTrans, 1._RKG / (1._RKG - help(1:locLastNegBreak)), 1._RKG / (1._RKG + help(lenBreak : locLastNegBreak + 1 : -1)))
        !print *, help(1:locLastNegBreak); print *, help(locLastNegBreak+1:); print *, breakTrans
#elif   QAGP_ENABLED && IF_ENABLED
        real(RKG) :: breakTrans(size(help, kind = IK))
        breakTrans = 1._RKG / (1._RKG - help + ub)
#elif   QAGP_ENABLED && FI_ENABLED
        real(RKG) :: breakTrans(size(help, kind = IK))
        breakTrans = 1._RKG / (1._RKG + help(size(help) : 1 : -1) - lb)
#elif   !(QAGD_ENABLED || QAGS_ENABLED || QAWC_ENABLED || QAWFS_ENABLED || QAWFC_ENABLED)
#error  "Unrecognized interface."
#endif
        err = getQuadErr(getFuncTrans, 0._RKG, 1._RKG, abstol, reltol, QRULE_ARG HELP_ARG, integral, abserr, sinfo, sindex, neval, nint)
    contains
        function getFuncTrans(x) result(func)
            implicit none
            real(RKG), intent(in)   :: x
            real(RKG)               :: func
            real(RKG)               :: invx
            real(RKG)               :: xMinuxOne2x
            invx = 1._RKG / x
            xMinuxOne2x = (x - 1._RKG) * invx
#if         II_ENABLED
            func = (getFunc(xMinuxOne2x) + getFunc(-xMinuxOne2x)) * TRANS_COEF
#elif       IF_ENABLED
            func = getFunc(ub + xMinuxOne2x) * TRANS_COEF
#elif       FI_ENABLED
            func = getFunc(lb - xMinuxOne2x) * TRANS_COEF
#else
#error      "Unrecognized interface."
#endif
        end function

        !%%%%%%%%%%%%%%%%
#elif   getQuadGK_ENABLED
        !%%%%%%%%%%%%%%%%

#if     GKXX_ENABLED
#define SET_PARAMETER(x) x
#elif   GK15_ENABLED || GK21_ENABLED || GK31_ENABLED || GK41_ENABLED || GK51_ENABLED || GK61_ENABLED
#define SET_PARAMETER(x) parameter(x)
#else
#error  "Unrecognized interface."
#endif
        ! Define the nodes for different rules.
#if     GK15_ENABLED
       !real(RKG)   , parameter :: nodeG(*) = real(nodeG7(size(nodeG7):1:-1), RKG)
        real(RKG)   , parameter :: nodeK(*) = real(nodeK15(size(nodeK15):1:-1), RKG)
        real(RKG)   , parameter :: weightG(*) = real(weightG7(size(weightG7):1:-1), RKG)
        real(RKG)   , parameter :: weightK(*) = real(weightK15(size(weightK15):1:-1), RKG)
#elif   GK21_ENABLED
       !real(RKG)   , parameter :: nodeG(*) = real(nodeG10(size(nodeG10):1:-1), RKG)
        real(RKG)   , parameter :: nodeK(*) = real(nodeK21(size(nodeK21):1:-1), RKG)
        real(RKG)   , parameter :: weightG(*) = real(weightG10(size(weightG10):1:-1), RKG)
        real(RKG)   , parameter :: weightK(*) = real(weightK21(size(weightK21):1:-1), RKG)
#elif   GK31_ENABLED
       !real(RKG)   , parameter :: nodeG(*) = real(nodeG15(size(nodeG15):1:-1), RKG)
        real(RKG)   , parameter :: nodeK(*) = real(nodeK31(size(nodeK31):1:-1), RKG)
        real(RKG)   , parameter :: weightG(*) = real(weightG15(size(weightG15):1:-1), RKG)
        real(RKG)   , parameter :: weightK(*) = real(weightK31(size(weightK31):1:-1), RKG)
#elif   GK41_ENABLED
       !real(RKG)   , parameter :: nodeG(*) = real(nodeG20(size(nodeG20):1:-1), RKG)
        real(RKG)   , parameter :: nodeK(*) = real(nodeK41(size(nodeK41):1:-1), RKG)
        real(RKG)   , parameter :: weightG(*) = real(weightG20(size(weightG20):1:-1), RKG)
        real(RKG)   , parameter :: weightK(*) = real(weightK41(size(weightK41):1:-1), RKG)
#elif   GK51_ENABLED
       !real(RKG)   , parameter :: nodeG(*) = real(nodeG25(size(nodeG25):1:-1), RKG)
        real(RKG)   , parameter :: nodeK(*) = real(nodeK51(size(nodeK51):1:-1), RKG)
        real(RKG)   , parameter :: weightG(*) = real(weightG25(size(weightG25):1:-1), RKG)
        real(RKG)   , parameter :: weightK(*) = real(weightK51(size(weightK51):1:-1), RKG)
#elif   GK61_ENABLED
       !real(RKG)   , parameter :: nodeG(*) = real(nodeG30(size(nodeG30):1:-1), RKG)
        real(RKG)   , parameter :: nodeK(*) = real(nodeK61(size(nodeK61):1:-1), RKG)
        real(RKG)   , parameter :: weightG(*) = real(weightG30(size(weightG30):1:-1), RKG)
        real(RKG)   , parameter :: weightK(*) = real(weightK61(size(weightK61):1:-1), RKG)
#elif   !GKXX_ENABLED
#error  "Unrecognized interface."
#endif
        real(RKG)   , parameter :: TINY_RKG = tiny(0._RKG)
        real(RKG)   , parameter :: EPS_RKG = epsilon(0._RKG)
        integer(IK)             :: i, ishared, sizeNodeG
        real(RKG)               :: sumFunc
        real(RKG)               :: midFunc      ! function value at the center of the interval.
        real(RKG)               :: avgFunc      ! approximation to the mean value of `getFunc` over `(lb, ub)`, i.e. to `quadGK / (ub - lb)`
        real(RKG)               :: quadG        ! integral via the 7-point gauss formula
        real(RKG)               :: quadK        ! integral via the 15-point kronrod formula
        real(RKG)               :: abscissa
        real(RKG)               :: midInterval
        real(RKG)               :: halfInterval
        real(RKG)               :: absHalfInterval, NegFunc(size(nodeK, kind = IK) - 1_IK), PosFunc(size(nodeK, kind = IK) - 1_IK)
        integer(IK)             :: offset

        ! integration limits.
#if     FF_ENABLED
#define LBT lb
#define UBT ub
#else
        real(RKG) :: invx, xMinuxOne2x
#define LBT 0._RKG
#define UBT 1._RKG
#endif
        ! integration function.
#if     FF_ENABLED
#define EVAL_EXPR(y, x) y = getFunc(x)
#elif   FI_ENABLED
#define EVAL_EXPR(y, x) invx = 1._RKG / (x); xMinuxOne2x = (x - 1._RKG) * invx; y = getFunc(lb - xMinuxOne2x) * invx**2
#elif   IF_ENABLED
#define EVAL_EXPR(y, x) invx = 1._RKG / (x); xMinuxOne2x = (x - 1._RKG) * invx; y = getFunc(ub + xMinuxOne2x) * invx**2
#elif   II_ENABLED
#define EVAL_EXPR(y, x) invx = 1._RKG / (x); xMinuxOne2x = (x - 1._RKG) * invx; y = (getFunc(xMinuxOne2x) + getFunc(-xMinuxOne2x)) * invx**2
#else
#error  "Unrecognized interface."
#endif
        SET_PARAMETER(sizeNodeG = int(0.5 * real(size(nodeK, kind = IK)), IK))
        SET_PARAMETER(offset = merge(1_IK, 0_IK, mod(size(nodeK, kind = IK) - 1_IK, 2_IK) == 1_IK)) ! is 1 when the Gauss abscissas includes the origin 0.
#if     GKXX_ENABLED
        CHECK_ASSERTION(__LINE__, size(weightK) == size(nodeK), SK_"@getQuadGK(): The condition `size(weightK) == size(nodeK)` must hold. size(weightK), size(nodeK) = "//getStr([size(weightK), size(nodeK)]))
        CHECK_ASSERTION(__LINE__, size(weightG) == size(nodeK) / 2, SK_"@getQuadGK(): The condition `size(weightG) == size(nodeK)/2` must hold. size(weightG), size(nodeK)/2 = "//getStr([size(weightG), size(nodeK)/2]))
#else
        CHECK_ASSERTION(__LINE__, precision(nodeK) <= 100, SK_"@getQuadGK(): The condition `precision(nodeK) <= 100` must hold. precision(nodeK) = "//getStr(precision(nodeK)))
#endif
        midInterval     = 0.5_RKG * (LBT + UBT)
        halfInterval    = 0.5_RKG * (UBT - LBT)
        absHalfInterval = abs(halfInterval)

        ! Initialize the summation depending on whether zero is in the abscissas of the Gauss summation.

        EVAL_EXPR(midFunc, midInterval) ! fpp
        if (offset == 1_IK) then
            quadG = midFunc * weightG(size(weightG))
        else
            quadG = 0._RKG
        end if
        quadK = midFunc * weightK(size(weightK))
        intAbsFunc = abs(quadK)

        ! Gauss summation.

        do i = 1_IK, sizeNodeG - offset
            ishared = 2_IK * i
            abscissa = halfInterval * nodeK(ishared)
            EVAL_EXPR(NegFunc(ishared), midInterval - abscissa) ! fpp
            EVAL_EXPR(PosFunc(ishared), midInterval + abscissa) ! fpp
            sumFunc = NegFunc(ishared) + PosFunc(ishared)
            quadG = quadG + sumFunc * weightG(i)
            quadK = quadK + sumFunc * weightK(ishared)
            intAbsFunc = intAbsFunc + weightK(ishared) * (abs(NegFunc(ishared)) + abs(PosFunc(ishared)))
        end do

        ! Kronrod summation.

        do i = 1_IK, 2_IK * sizeNodeG - 1_IK, 2_IK
            abscissa = halfInterval * nodeK(i)
            EVAL_EXPR(NegFunc(i), midInterval - abscissa) ! fpp
            EVAL_EXPR(PosFunc(i), midInterval + abscissa) ! fpp
            sumFunc = NegFunc(i) + PosFunc(i)
            quadK = quadK + weightK(i) * sumFunc
            intAbsFunc = intAbsFunc + weightK(i) * (abs(NegFunc(i)) + abs(PosFunc(i)))
        end do

        ! Estimate the error.

        avgFunc = quadK * 0.5_RKG
        smoothness = weightK(1) * abs(midFunc - avgFunc)
        do i = 1_IK, size(nodeK, kind = IK) - 1_IK
            smoothness = smoothness + weightK(i) * (abs(NegFunc(i) - avgFunc) + abs(PosFunc(i) - avgFunc))
        end do
        quadGK = quadK * halfInterval
        intAbsFunc = intAbsFunc*absHalfInterval
        smoothness = smoothness*absHalfInterval
        abserr = abs((quadK - quadG) * halfInterval)
        !write(*,*) "======================================================================================"
        !print *, "nodeK", nodeK
        !print *, "weightK", weightK
        !print *, "weightG", weightG
        !print *, "size(nodeK), size(weightK), size(weightG)", size(nodeK), size(weightK), size(weightG)
        !print *, "abserr, quadK, quadG", abserr, quadK, quadG, halfInterval
        !write(*,*) "======================================================================================"
        !write(*,*)
        if (smoothness /= 0._RKG .and. abserr /= 0._RKG) then
            abserr = 200._RKG * abserr / smoothness ! `abserr` represents relative error here.
            if (abserr < 1._RKG) then
                abserr = smoothness * sqrt(abserr)**3
            else
                abserr = smoothness
            end if
            !abserr = smoothness * min(1._RKG, sqrt(200._RKG * abserr / smoothness)**3) ! the sqrt and exponentiation can be avoided here.
        end if
        if (intAbsFunc > TINY_RKG / (50._RKG * EPS_RKG)) abserr = max((EPS_RKG * 50._RKG) * intAbsFunc, abserr)
#undef  SET_PARAMETER
#undef  EVAL_EXPR
#undef  LBT
#undef  UBT

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getQuadErr_ENABLED && (QAGD_ENABLED || QAGS_ENABLED || QAGP_ENABLED || QAWC_ENABLED) && FF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     QAWC_ENABLED
        real(RKG) :: cs
#endif
        real(RKG)   , parameter :: HUGE_RKG = huge(0._RKG)
        real(RKG)   , parameter :: TINY_RKG = tiny(0._RKG)
        real(RKG)   , parameter :: EPS_RKG = epsilon(0._RKG)
        integer(IK)             :: iroff1, iroff2
        integer(IK)             :: nintmax, nrmax, maxErrLoc    ! pointer to the interval with largest error estimate.
        real(RKG)               :: result1, lb1, ub1, abserr1   ! variables for the left  subinterval.
        real(RKG)               :: result2, lb2, ub2, abserr2   ! variables for the right subinterval.
        real(RKG)               :: result                       ! sum of the integrals over the subintervals
        real(RKG)               :: result12                     ! `result1 + result2`
        real(RKG)               :: abserr12                     ! `abserr1 + abserr2`
        real(RKG)               :: maxErrVal                    ! `elist(maxErrLoc)`.
        real(RKG)               :: errsum                       ! sum of the errors over the subintervals.
        real(RKG)               :: errbnd                       ! requested accuracy `max(abstol, reltol * abs(integral))`.
        real(RKG)               :: intAbsFunc1, intAbsFunc2
        real(RKG)               :: smoothness1, smoothness2
#if     QAGD_ENABLED || QAWC_ENABLED
#define FINALIZE_INTEGRATION integral = sum(sinfo(3,1:nint)); abserr = errsum; return
        real(RKG)   , parameter :: FACTOR = 1._RKG
#elif   QAGS_ENABLED || QAGP_ENABLED
#define FINALIZE_INTEGRATION integral = sum(sinfo(3,1:nint)); abserr = errsum; if (err > 2_IK) err = err - 1_IK; return
        real(RKG)   , parameter :: FACTOR = 2._RKG
        integer(IK)             :: ksgn
        real(RKG)               :: EpsTable(MAXLEN_EPSTAB + 2)
        real(RKG)               :: absErrEps, seqlim, seqlims(3)
        ! error on the interval currently subdivided;
        ! sum of the errors over the intervals larger than the smallest interval considered up to now.
        real(RKG)               :: erlast, erlarg, correc, ertest
        integer(IK)             :: lenEpsTable, nextrap, iroffnew, k, ktmin, jupbnd, ierro ! roundoff error flag
        logical(LK)             :: extrapolating, extrapAllowed
#if     QAGS_ENABLED
#define CYCLE_NEEDED abs(sinfo(2,maxErrLoc) - sinfo(1,maxErrLoc)) > small
#define EPS_TABLE_START_OFFSET 1_IK
        real(RKG)               :: small
#elif   QAGP_ENABLED
#define CYCLE_NEEDED level(maxErrLoc) < maximumLevel
#define EPS_TABLE_START_OFFSET 0_IK
        integer(IK)             :: i, j, jlow, ind1, ind2, lenBreak, lenBreakPlusOne
        integer(IK)             :: currentLevel, maximumLevel, level(size(sindex, kind = IK))
        logical(LK)             :: refinementNeeded(size(help, kind = IK) + 1_IK)
#else
#error  "Unrecognized interface."
#endif
#else
#error  "Unrecognized interface."
#endif
#if     GK15_ENABLED
        integer(IK) , parameter :: NEVAL0 = 15_IK
#elif   GK21_ENABLED
        integer(IK) , parameter :: NEVAL0 = 21_IK
#elif   GK31_ENABLED
        integer(IK) , parameter :: NEVAL0 = 31_IK
#elif   GK41_ENABLED
        integer(IK) , parameter :: NEVAL0 = 41_IK
#elif   GK51_ENABLED
        integer(IK) , parameter :: NEVAL0 = 51_IK
#elif   GK61_ENABLED
        integer(IK) , parameter :: NEVAL0 = 61_IK
#elif   GKXX_ENABLED
        integer(IK)             :: NEVAL0
        NEVAL0 = 2_IK * size(nodeK, kind = IK) - 1_IK
#else
#error  "Unrecognized interface."
#endif
        nintmax = size(sindex, kind = IK)
#if     QAWC_ENABLED
        cs = real(help%cs, RKG)
        CHECK_ASSERTION(__LINE__, lb < cs, SK_"@getQuadErr(): The condition `lb < cs` must hold. lb, cs = "//getStr([lb, cs]))
        CHECK_ASSERTION(__LINE__, cs < ub, SK_"@getQuadErr(): The condition `cs < ub` must hold. cs, ub = "//getStr([cs, ub]))
#endif
        !   \bug Intel ifort 2021.6
        !   The intel ifort yields an internal compiler error which appears to be related to the following macro placements.
        !   As such they are commented out and implemented explicitly as all-in-one checking block.
        !   Update Dec 2022: This bug is likely related to the Intel bug with the maximum number
        !   of `use` statements in subdmoules in [pm_araryRebill](@ref pm_araryRebill).
        CHECK_ASSERTION(__LINE__, lb < ub, SK_"@getQuadErr(): The condition `lb < ub` must hold. lb, ub = "//getStr([lb, ub]))
        CHECK_ASSERTION(__LINE__, abstol >= 0._RKG, SK_"@getQuadErr(): The condition `abstol >= 0.` must hold. abstol = "//getStr(abstol))
        CHECK_ASSERTION(__LINE__, reltol >= 0._RKG, SK_"@getQuadErr(): The condition `reltol >= 0.` must hold. reltol = "//getStr(reltol))
        CHECK_ASSERTION(__LINE__, all(shape(sinfo, kind = IK) == [4_IK, nintmax]), SK_"@getQuadErr(): The condition `all(shape(sinfo) == [4, size(sindex)])` must hold. shape(sinfo), size(sindex) = "//getStr([shape(sinfo, kind = IK), size(sindex, kind = IK)]))
        CHECK_ASSERTION(__LINE__, size(sinfo, 2, IK) == nintmax, SK_"@getQuadErr(): The condition `size(sinfo,2) == size(sindex)` must hold. size(sinfo,2), size(sindex) = "//getStr([size(sinfo, 2, IK), nintmax]))
        CHECK_ASSERTION(__LINE__, nintmax >= 1_IK, SK_"@getQuadErr(): The condition `size(sindex, kind = IK) >= 1` must hold. size(sinfo,2) = "//getStr(size(sindex, kind = IK)))
#if     QAGP_ENABLED
        neval = 0_IK
        lenBreak = size(help, kind = IK)
        lenBreakPlusOne = lenBreak + 1_IK
        CHECK_ASSERTION(__LINE__, size(help, kind = IK) > 0_IK, SK_"@getQuadErr(): The condition `size(help) > 0` must hold. size(help) = "//getStr(size(help, kind = IK)))
        CHECK_ASSERTION(__LINE__, size(help, kind = IK) < nintmax, SK_"@getQuadErr(): The condition `size(help) < size(sindex)` must hold. size(help), size(sindex) = "//getStr([size(help, kind = IK), nintmax]))
        CHECK_ASSERTION(__LINE__, all(help > lb), SK_"@getQuadErr(): The condition `all(lb < help)` must hold for the input arguments `help` and `lb`. lb, help = "//getStr([lb, help]))
        CHECK_ASSERTION(__LINE__, all(help < ub), SK_"@getQuadErr(): The condition `all(help < ub)` must hold for the input arguments `help` and `ub`. help, ub = "//getStr([help, ub]))
        CHECK_ASSERTION(__LINE__, size(sinfo, 2, IK) == nintmax, SK_"@getQuadErr(): The condition `size(sinfo,2) == size(sindex)` must hold. size(sinfo,2), size(sindex) = "//getStr([size(sinfo, 2, IK), nintmax]))
        CHECK_ASSERTION(__LINE__, isAscending(help), SK_"@getQuadErr(): The input argument `help` must be in ascending order. help = "//getStr(help))
        ! Define the first set of interval edges.
        sinfo(1,1) = lb
        do i = 1_IK, lenBreak
            sinfo(2, i) = help(i)
            sinfo(1, i + 1) = sinfo(2, i)
        end do
        sinfo(2,i) = ub
        ! Compute the first approximation to the integral.
        abserr = 0._RKG
        integral = 0._RKG
        intAbsFunc1 = 0._RKG
        do i = 1_IK, lenBreakPlusOne ! loop over the first set of intervals.
            result1 = getQuadGK(getFunc, sinfo(1,i), sinfo(2,i), QRULE_ARG, abserr1, intAbsFunc2, smoothness1) ! Compute the first approximation to the integral.
            neval = neval + NEVAL0
            refinementNeeded(i) = abserr1 /= 0._RKG .and. abserr1 == smoothness1
            intAbsFunc1 = intAbsFunc1 + intAbsFunc2
            integral = integral + result1
            abserr = abserr + abserr1
            sinfo(3,i) = result1
            sinfo(4,i) = abserr1
            level(i) = 0_IK
            sindex(i) = i
        end do
        ! Set the overall error estimate.
        errsum = 0._RKG
        do i = 1_IK, lenBreakPlusOne
            if (refinementNeeded(i)) sinfo(4,i) = abserr
            errsum = errsum + sinfo(4,i)
        end do
        nint = lenBreakPlusOne
        ! Sort the error estimate indices in descending order of error estimates.
        do i = 1_IK, lenBreak
            jlow = i + 1_IK
            ind1 = sindex(i)
            do j = jlow, lenBreakPlusOne
                ind2 = sindex(j)
                if (sinfo(4,ind1) <= sinfo(4,ind2)) then
                    ind1 = ind2
                    k = j
                end if
            end do
            if (ind1 /= sindex(i)) then
                sindex(k) = sindex(i)
                sindex(i) = ind1
            end if
        end do
#elif   QAGD_ENABLED || QAGS_ENABLED || QAWC_ENABLED
#if     QAWC_ENABLED
        neval = 0_IK
        call setQuadQAWC(integral, lb, ub, abserr, intAbsFunc1, smoothness1, neval)
#elif   QAGD_ENABLED || QAGS_ENABLED
        integral = getQuadGK(getFunc, lb, ub, QRULE_ARG, abserr, intAbsFunc1, smoothness1)
        neval = NEVAL0
#endif
        sinfo(1,1) = lb
        sinfo(2,1) = ub
        sinfo(3,1) = integral
        sinfo(4,1) = abserr
        sindex(1) = 1_IK
        nint = 1_IK
#else
#error  "Unrecognized interface."
#endif
        ! Test the accuracy of the zeroth order integration.
        err = 0_IK
        errbnd = max(abstol, reltol * abs(integral))
        if (nint == nintmax) then; err = 1_IK; return; end if
        if (errbnd < abserr .and. abserr <= FACTOR * 50._RKG * EPS_RKG * intAbsFunc1) then; err = 2_IK; return; end if
#if     QAGP_ENABLED
        if (abserr <= errbnd) return
#elif   QAGD_ENABLED || QAGS_ENABLED || QAWC_ENABLED
        if (abserr == 0._RKG) return
        if (abserr <= errbnd .and. abserr /= smoothness1) return
#else
#error  "Unrecognized interface."
#endif
!#if    QAWC_ENABLED
!       if (abserr <= min(0.01_RKG * abs(integral), errbnd)) return
!#endif
        ! Refine the integration.
#if     QAGP_ENABLED
        erlarg = errsum
        ertest = errbnd
        maximumLevel = 1_IK
#elif   QAGD_ENABLED || QAGS_ENABLED || QAWC_ENABLED
        errsum = abserr
#else
#error  "Unrecognized interface."
#endif
        result = integral
        maxErrLoc = sindex(1)
        maxErrVal = sinfo(4,maxErrLoc)
        iroff1 = 0_IK
        iroff2 = 0_IK
        nrmax = 1_IK
#if     QAGS_ENABLED || QAGP_ENABLED
        ierro = 0_IK
        ktmin = 0_IK
        nextrap = 0_IK
        iroffnew = 0_IK
        abserr = HUGE_RKG
        EpsTable(1) = integral
        lenEpsTable = 1_IK + EPS_TABLE_START_OFFSET ! remark: The second element of `EpsTable` for QAGS will be set later on.
        extrapAllowed = .true._LK
        extrapolating = .false._LK
        ! \warning \bug The quadpack.F90 version of John Burkardt is different from the original quadpack, where 50 is replaced with 0.5 in QAGP.
        ksgn = -1_IK; if (abs(integral) >= (1._RKG - 50._RKG * EPS_RKG) * intAbsFunc1) ksgn = 1_IK
#elif   !(QAGD_ENABLED || QAWC_ENABLED)
#error  "Unrecognized interface."
#endif
        !print *, "abserr, errsum", abserr, errsum
        ! Adaptively integrate the integrand.
        loopAdaptation: do nint = nint + 1_IK, nintmax
            ! Bisect the subinterval that has the largest error estimate.
            lb1 = sinfo(1,maxErrLoc)
            ub2 = sinfo(2,maxErrLoc)
            ub1 = 0.5_RKG * (lb1 + ub2)
#if         QAWC_ENABLED
            if (lb1 < cs .and. cs <= ub1) then
                ub1 = 0.5_RKG * (cs + ub2)
            elseif (ub1 < cs .and. cs <  ub2) then
                ub1 = 0.5_RKG * (lb1 + cs)
            end if
#endif
            lb2 = ub1
#if         QAGS_ENABLED || QAGP_ENABLED
            erlast = maxErrVal
#if         QAGP_ENABLED
            currentLevel = level(maxErrLoc) + 1_IK
#endif
#elif       !(QAGD_ENABLED || QAWC_ENABLED)
#error      "Unrecognized interface."
#endif
#if         QAWC_ENABLED
            call setQuadQAWC(result1, lb1, ub1, abserr1, intAbsFunc2, smoothness1, neval)
            call setQuadQAWC(result2, lb2, ub2, abserr2, intAbsFunc2, smoothness2, neval)
            !print *, abserr1, abserr2, lb1 < cs .and. cs <= ub1
#elif       QAGD_ENABLED || QAGS_ENABLED || QAGP_ENABLED
            result1 = getQuadGK(getFunc, lb1, ub1, QRULE_ARG, abserr1, intAbsFunc2, smoothness1)
            result2 = getQuadGK(getFunc, lb2, ub2, QRULE_ARG, abserr2, intAbsFunc2, smoothness2)
            neval = neval + 2_IK * NEVAL0
#else
#error      "Unrecognized interface."
#endif
            ! Improve the previous approximations to integral and error and verify the accuracy.
            result12 = result1 + result2
            abserr12 = abserr1 + abserr2
            errsum = errsum + abserr12 - maxErrVal
            result = result + result12 - sinfo(3,maxErrLoc)
            if (smoothness1 /= abserr1 .and. smoothness2 /= abserr2) then
                if (nint > 10_IK .and. abserr12 > maxErrVal) iroff2 = iroff2 + 1_IK
                if (abs(sinfo(3,maxErrLoc) - result12) <= 1.e-05_RKG * abs(result12) .and. abserr12 >= 0.99_RKG * maxErrVal) then
#if                 QAGD_ENABLED || QAWC_ENABLED
                    iroff1 = iroff1 + 1_IK
#elif               QAGS_ENABLED || QAGP_ENABLED
                    if (extrapolating) then
                        iroffnew = iroffnew + 1_IK
                    else
                        iroff1 = iroff1 + 1_IK
                    end if
#else
#error              "Unrecognized interface."
#endif
                end if
            end if
#if         QAGP_ENABLED
            level(maxErrLoc) = currentLevel
            level(nint) = currentLevel
#endif
            ! Append the newly-created intervals to the list.
            CHECK_ASSERTION(__LINE__, nint /= maxErrLoc, SK_"@getQuadErr(): The condition `nint /= maxErrLoc` happened. This in an internal error. Please report this issue to the developers.")
            if (abserr2 > abserr1) then ! \todo is this order necessary?
                sinfo(1,maxErrLoc) = lb2
                sinfo(1,nint) = lb1
                sinfo(2,nint) = ub1
                sinfo(3,maxErrLoc) = result2
                sinfo(3,nint) = result1
                sinfo(4,maxErrLoc) = abserr2
                sinfo(4,nint) = abserr1
            else
                sinfo(1,nint) = lb2
                sinfo(2,maxErrLoc) = ub1
                sinfo(2,nint) = ub2
                sinfo(3,maxErrLoc) = result1
                sinfo(3,nint) = result2
                sinfo(4,maxErrLoc) = abserr1
                sinfo(4,nint) = abserr2
            end if
            ! Ensure the descending ordering of the list of error estimates and select the subinterval with the largest error estimate (to be bisected next).
            call setErrSorted(nintmax, sinfo(4,1:nint), sindex(1:nint), nrmax = nrmax, maxErrLoc = maxErrLoc, maxErrVal = maxErrVal)
            ! Check for convergence/divergence.
            errbnd = max(abstol, reltol * abs(result))
            if (errsum <= errbnd) then ! convergence occurred.
                FINALIZE_INTEGRATION
            end if
#if         QAGD_ENABLED || QAWC_ENABLED
            if (iroff1 >= 6_IK .or. iroff2 >= 20_IK) err = 2_IK ! test for roundoff error and eventually set error flag.
            if (nint == nintmax) err = 1_IK ! set error flag in the case that the number of subintervals equals nintmax.
            ! \warning \bug QuadPack version of John Burkardt, replaces 1000 in the original code with 10000 in the condition below.
            if (max(abs(lb1), abs(ub2)) <= (1._RKG + NEVAL0 * 100._RKG * EPS_RKG) * (abs(lb2) + 1000._RKG * TINY_RKG)) err = 3_IK !>  \warning \todo The impact of NEVAL0 in this condition may need further investigation for high-order > 61 rules.
#elif       QAGS_ENABLED || QAGP_ENABLED
            if (iroff1 + iroffnew >= 10_IK .or. iroff2 >= 20) err = 2_IK ! roundoff error.
            if (iroffnew >= 5_IK) ierro = 3_IK
            !if (iroffnew >= 5_IK) error stop __LINE__
            if (nint == nintmax) err = 1_IK
            ! \warning \bug QuadPack version of John Burkardt, replaces 100 in the original code with 1000 in the condition below.
            if (max(abs(lb1), abs(ub2)) <= (1._RKG + 100._RKG * EPS_RKG) * (abs(lb2) + 1000._RKG * TINY_RKG)) err = 4_IK ! bad integrand behavior.
#else
#error      "Unrecognized interface."
#endif
            if (err /= 0_IK) exit loopAdaptation
#if         QAGS_ENABLED
            ! Extrapolate, if allowed.
            if (nint == 2_IK) then
                small = abs(ub - lb) * 0.375_RKG
                erlarg = errsum
                ertest = errbnd
                EpsTable(2) = result
                cycle loopAdaptation
            end if
#endif
#if         QAGS_ENABLED || QAGP_ENABLED
            if (extrapAllowed) then
                erlarg = erlarg - erlast
                if  ( & ! LCOV_EXCL_LINE
#if                 QAGS_ENABLED
                    small < abs(ub1 - lb1) & ! LCOV_EXCL_LINE
#elif               QAGP_ENABLED
                    currentLevel < maximumLevel & ! LCOV_EXCL_LINE
#else
#error              "Unrecognized intefrace."
#endif
                ) erlarg = erlarg + abserr12
                if (.not. extrapolating) then ! test whether the interval to be bisected next is the smallest interval.
                    if (CYCLE_NEEDED) cycle loopAdaptation
                    extrapolating = .true._LK
                    nrmax = 2_IK
                end if
                if (ierro /= 3_IK .and. erlarg > ertest) then
                    ! The smallest interval has the largest error. Before bisecting decrease the
                    ! sum of the errors over the larger intervals (erlarg) and perform extrapolation.
                    jupbnd = nint
                    if (nint > 2_IK + nintmax / 2_IK) jupbnd = nintmax + 3_IK - nint
                    do k = nrmax, jupbnd
                        maxErrLoc = sindex(nrmax)
                        maxErrVal = sinfo(4,maxErrLoc)
                        if (CYCLE_NEEDED) cycle loopAdaptation
                        nrmax = nrmax + 1_IK
                    end do
                end if
                ! perform extrapolation.
                lenEpsTable = lenEpsTable + 1_IK
                EpsTable(lenEpsTable) = result
#if             QAGP_ENABLED
                if (lenEpsTable > 2_IK) then
#elif           !QAGS_ENABLED
#error          "Unrecognized intefrace."
#endif
                    !print *, "lenEpsTable, nextrap = ", lenEpsTable, nextrap
                    call setSeqLimEps(lenEpsTable, nextrap, EpsTable, seqlims, seqlim, absErrEps)
                    ktmin = ktmin + 1_IK
                    !print *, "end0", "ktmin", ktmin
                    !print *, "end0", "errsum", errsum
                    !print *, "end0", "abserr", abserr
                    if (ktmin > 5_IK .and. abserr < 1.e-03_RKG * errsum) err = 5_IK
                    !print *, "end1", MAXLEN_EPSTAB, lenEpsTable
                    if (absErrEps < abserr) then
                        ktmin = 0_IK
                        abserr = absErrEps
                        integral = seqlim
                        correc = erlarg
                        ertest = max(abstol, reltol * abs(seqlim))
                        if (abserr <= ertest) exit loopAdaptation
                    end if
                    ! prepare bisection of the smallest interval.
                    if (lenEpsTable == 1_IK) extrapAllowed = .false._LK
                    if (err == 5_IK) exit loopAdaptation
#if             QAGP_ENABLED
                end if
                maximumLevel = maximumLevel + 1_IK
#elif           QAGS_ENABLED
                small = small * 0.5_RKG
#else
#error          "Unrecognized intefrace."
#endif
                maxErrLoc = sindex(1)
                maxErrVal = sinfo(4,maxErrLoc)
                extrapolating = .false._LK
                erlarg = errsum
                nrmax = 1_IK
            end if
#endif

        end do loopAdaptation
        CHECK_ASSERTION(__LINE__, nint <= nintmax, \
        SK_"@getQuadErr(): Internal library error: The condition `nint <= nintmax` must hold. nint, nintmax = ["//\
        getStr([nint, nintmax])//SK_"]. Please report this error to the ParaMonte library developers.")

#if     QAGS_ENABLED || QAGP_ENABLED
        if (abserr /= HUGE_RKG) then
            if (err + ierro /= 0_IK) then
                if (ierro == 3_IK) abserr = abserr + correc
                if (err == 0_IK) err = 3_IK
                !if (err == 3) error stop __LINE__
                if (integral == 0._RKG .or. result == 0._RKG) then
                    if (abserr > errsum) then
                        FINALIZE_INTEGRATION
                    end if
                    if (result == 0._RKG) then
                        if (err > 2_IK) err = err - 1_IK
                        return
                    end if
                elseif (abserr / abs(integral) > errsum / abs(result)) then
                    FINALIZE_INTEGRATION
                end if
            end if
            if (ksgn /= -1_IK .or. max(abs(integral), abs(result)) > intAbsFunc1 * 0.01_RKG) then ! test on divergence.
                if (0.01_RKG > integral / result .or. integral / result > 100._RKG .or. errsum > abs(result)) err = 6_IK
            end if
            if (err > 2_IK) err = err - 1_IK
            return
        end if
#elif   !(QAGD_ENABLED || QAWC_ENABLED)
#error  "Unrecognized interface."
#endif
        FINALIZE_INTEGRATION

    contains

#if     QAWC_ENABLED
        !   Approximate the integral of the input function using,
        !   <ul>
        !       <li>    a generalized clenshaw-curtis method if `cs` lies within ten percent of the integration limits, otherwise,
        !       <li>    a user-specified Gauss-Kronrod quadrature rule.
        !   </ul>
        !   This routine uses `getFunc`, `getFuncWC`, `NEVAL0`, and `cs` objects from the parent routine.
        !   However, it does not have any side effects (does not change any global variables).
        subroutine setQuadQAWC(quadQAWC, lb, ub, abserr, intAbsFunc, smoothness, neval)
            real(RKG)           , intent(in)    :: lb, ub           ! the integration lower/upper bound.
            real(RKG)           , intent(out)   :: abserr           ! the abslue value of the integral error estimate, which equals or exceeds `abs(i - quadQAWC)`.
            real(RKG)           , intent(out)   :: intAbsFunc       ! integral of the absolute of the function.
            real(RKG)           , intent(out)   :: smoothness       ! An integrand smoothness measure.
            integer(IK)         , intent(inout) :: neval            ! number of function evaluations.
            real(RKG)           , intent(out)   :: quadQAWC         ! approximation to the integral.
            integer(IK)                         :: i, j
            real(RKG)           , parameter     :: NODE(11) = [(cos(i * acos(-1._RKG) / 24._RKG), i = 1_IK, 11_IK)]
            real(RKG)                           :: ak22, amom0, amom1, amom2, u, cstrans
            real(RKG)                           :: cheb12(13)           ! Chebyshev series expansion of degree 12.
            real(RKG)                           :: cheb24(25)           ! Chebyshev series expansion of degree 24.
            real(RKG)                           :: func(25)             ! Function values at `cos(k*pi/24)`, `k = 0, ..., 24`.
            real(RKG)                           :: result12, result24   ! approximations to the integral Via Clenshaw-Curtis of orders 13 and 25.
            real(RKG)                           :: midInterval, halfInterval
            cstrans = (2._RKG * cs - ub - lb) / (ub - lb) ! Compute the position of cs.
            !print *, cstrans
            if (abs(cstrans) < 1.1_RKG) then ! Use the generalized Clenshaw-Curtis method.
                smoothness = huge(smoothness)
                halfInterval = 0.5_RKG * (ub - lb)
                midInterval = 0.5_RKG * (ub + lb)
                func(1) = 0.5_RKG * getFunc(halfInterval + midInterval)
                func(13) = getFunc(midInterval)
                func(25) = 0.5_RKG * getFunc(midInterval - halfInterval)
                do i = 2_IK, 12_IK
                    u = halfInterval * NODE(i - 1_IK)
                    func(i) = getFunc(u + midInterval)
                    func(26_IK - i) = getFunc(midInterval - u)
                end do
                call setChebExpan(func, cheb12, cheb24) ! compute the Chebyshev series expansion.
                ! the modified Chebyshev moments are computed by forward recursion, using amom0 and amom1 as starting values.
                amom0 = log(abs((1._RKG - cstrans) / (1._RKG + cstrans)))
                amom1 = 2._RKG + cstrans * amom0
                result12 = cheb12(1) * amom0 + cheb12(2) * amom1
                result24 = cheb24(1) * amom0 + cheb24(2) * amom1
                intAbsFunc = abs(cheb24(1))
                j = 1_IK
                do i = 3_IK, 13_IK
                    amom2 = 2._RKG * cstrans * amom1 - amom0
                    ak22 = real((i - 2_IK) * (i - 2_IK), RKG)
                    if (i == (i / 2_IK) * 2_IK) amom2 = amom2 - 4._RKG / (ak22 - 1._RKG)
                    result12 = result12 + cheb12(i) * amom2
                    result24 = result24 + cheb24(i) * amom2
                    amom0 = amom1
                    amom1 = amom2
                    j = j + 2_IK
                    intAbsFunc = intAbsFunc + abs(cheb24(j - 1)) + abs(cheb24(j))
                end do
                do i = 14_IK, 25_IK
                    amom2 = 2._RKG * cstrans * amom1 - amom0
                    ak22 = real((i - 2_IK) * (i - 2_IK), RKG)
                    if (i == (i / 2_IK) * 2_IK) amom2 = amom2 - 4._RKG / (ak22 - 1._RKG)
                    result24 = result24 + cheb24(i) * amom2
                    amom0 = amom1
                    amom1 = amom2
                end do
                abserr = abs(result24 - result12)
                quadQAWC = result24
                neval = neval + 25_IK
            else ! Use the default Gauss-Kronrod quadrature.
                quadQAWC = getQuadGK(getFuncWC, lb, ub, QRULE_ARG, abserr, intAbsFunc, smoothness)
                neval = neval + NEVAL0
            end if
        end subroutine

        ! Cauchy-Weighted function.

        function getFuncWC(x) result(funcWC)
            real(RKG)   , intent(in)    :: x
            real(RKG)                   :: funcWC
            funcWC = getFunc(x) / (x - cs)
        end function

#elif   QAWFS_ENABLED || QAWFC_ENABLED

        subroutine getQuadWF(quadWF, lb, ub, currentLevel, ksave, abserr, intAbsFunc, smoothness, chebmom, momcom, neval)
            real(RKG)   , intent(in)                    :: lb, ub           !<  The current integration limits.
            integer(IK) , intent(in)                    :: ksave            !<  The key which is one when the moments for the current interval have been computed.
            integer(IK) , intent(in)                    :: currentLevel     !<  The length of interval `(lb,ub)` is equal to the length of the original integration interval divided by `2**currentLevel`.
                                                                            !!  We assume that the routine is used in an adaptive integration process, otherwise set `currentLevel = 0`.
                                                                            !!  `currentLevel` must be zero at the first call.
            real(RKG)   , intent(out)                   :: abserr           !<  estimate of the modulus of the absolute error, which should equal or exceed `abs(i-quadWF)`.
            real(RKG)   , intent(out)                   :: intAbsFunc       !<  Approximation to the integral of the absolute of the function `abs(f(x))`.
            real(RKG)   , intent(out)                   :: smoothness       !<  Approximation to the integral of `abs(f - i / (ub - lb))`.
            real(RKG)   , intent(inout) , contiguous    :: chebmom(:, 25)   !<  Array of shape `(maxLenChebMom,25)` containing the modified Chebyshev moments for the first `momcom` interval lengths.
            integer(IK) , intent(inout)                 :: momcom           !<  For each interval length we need to compute the Chebyshev moments.
                                                                            !!  `momcom` counts the number of intervals for which these moments have already been computed.
                                                                            !!  If `currentLevel < momcom` or `ksave = 1`, the Chebyshev moments for the interval `(lb, ub)` have already been computed.
                                                                            !!  Otherwise we compute them and we increase `momcom`.
            integer(IK) , intent(out)                   :: neval            !<  Number of integrand evaluations.
            real(RKG)                                   :: quadWF           !<  approximation to the integral i

            real(RKG)                   , parameter     :: NODE(11) = [(cos(i * acos(-1._RKG) / 24._RKG), i = 1_IK, 11_IK)]
            integer(IK)                                 :: i, iers, isym, j, k, m, noequ, noeq1
            integer(IK)                                 :: maxLenChebMom        !<  The upper bound on the number of Chebyshev moments which can be stored,
                                                                                !<  i.e. for the intervals of lengths `abs(bb-aa)*2**(-l)`, `l = 0,1,2, ..., maxLenChebMom-2`.
            real(RKG)                                   :: ac, an, ansq, invansq, as, asap, ass, conc, cons, cospar, d(25), d1(25), d2(25), estc, ests, parint, par2, par22, p2, p3, p4, sinpar, Vector(28)
            real(RKG)                                   :: halfInterval
            real(RKG)                                   :: midInterval
            real(RKG)                                   :: cheb12(13)       !<  Chebyshev series expansion of degree 12.
            real(RKG)                                   :: cheb24(25)       !<  Chebyshev series expansion of degree 24.
            real(RKG)                                   :: func(25)         !<  Function values at `(ub - lb) * 0.5 * cos(k * pi / 12) + (ub + lb) * 0.5`, `k = 0, ..., 24`.
            real(RKG)                                   :: result12cos      !<  approximation to the integral of `cos(0.5 * (ub - lb) * omega * NODE) * getFunc(0.5 * (ub - lb) * NODE + 0.5 * (ub + lb))` over `(-1, +1)`, using the Chebyshev series expansion of degree 12.
            real(RKG)                                   :: result24cos      !<  approximation to the integral of `cos(0.5 * (ub - lb) * omega * NODE) * getFunc(0.5 * (ub - lb) * NODE + 0.5 * (ub + lb))` over `(-1, +1)`, using the Chebyshev series expansion of degree 12.
            real(RKG)                                   :: result12sin      !<  approximation to the integral of `sin(0.5 * (ub - lb) * omega * NODE) * getFunc(0.5 * (ub - lb) * NODE + 0.5 * (ub + lb))` over `(-1, +1)`, using the Chebyshev series expansion of degree 12.
            real(RKG)                                   :: result24sin      !<  approximation to the integral of `sin(0.5 * (ub - lb) * omega * NODE) * getFunc(0.5 * (ub - lb) * NODE + 0.5 * (ub + lb))` over `(-1, +1)`, using the Chebyshev series expansion of degree 12.

            midInterval = 0.5_RKG * (ub + lb)
            halfInterval = 0.5_RKG * (ub - lb)
            parint = Weight%omega * halfInterval
            maxLenChebMom = size(chebmom, 1_IK, IK)

            ! Compute the integral using the 15-point gauss-kronrod formula if the value of the parameter in the integrand is small.

            if (abs(parint) > 2._RKG) then

                ! Compute the integral using the generalized clenshaw-curtis method.

                conc = halfInterval * cos(midInterval * Weight%omega)
                cons = halfInterval * sin(midInterval * Weight%omega)
                smoothness = huge(smoothness)
                neval = 25_IK

                ! Check whether the Chebyshev moments for this interval have already been computed.

                if (momcom <= currentLevel .and. ksave /= 1_IK) then ! Compute a new set of Chebyshev moments.

                    m = momcom + 1_IK
                    par2 = parint*parint
                    par22 = par2 + 2._RKG
                    sinpar = sin(parint)
                    cospar = cos(parint)

                    ! compute the Chebyshev moments with respect to cosine.

                    Vector(1) = 2._RKG*sinpar / parint
                    Vector(2) = (8._RKG*cospar + (par2 + par2 - 8._RKG) * sinpar / parint) / par2
                    Vector(3) = (32._RKG*(par2 - 12._RKG)*cospar + (2._RKG * ((par2 - 80._RKG) * par2 + 192._RKG) * sinpar) / parint) / (par2*par2)
                    ac = 8._RKG * cospar
                    as = 24._RKG * parint * sinpar
                    if (abs(parint) > 24._RKG) then ! Compute the Chebyshev moments by means of forward recursion.
                        an = 4._RKG
                        do i = 4_IK, 13_IK
                            ansq = an * an
                            Vector(i) = ((ansq - 4._RKG) * (2._RKG * (par22 - ansq - ansq) * Vector(i - 1) - ac) + as - par2 * (an + 1._RKG) * (an + 2._RKG) * Vector(i - 2)) / (par2 * (an - 1._RKG) * (an - 2._RKG))
                            an = an + 2._RKG
                        end do
                    else ! Compute the Chebyshev moments as the solutions of a boundary value problem with 1 initial value (Vector(3)) and 1 end value (computed using an asymptotic formula).
                        noequ = 25_IK
                        noeq1 = noequ - 1_IK
                        an = 6._RKG
                        do k = 1_IK, noeq1
                            ansq = an * an
                            d(k) = -2._RKG * (ansq - 4._RKG) * (par22 - ansq - ansq)
                            d2(k) = (an - 1._RKG) * (an - 2._RKG) * par2
                            d1(k + 1) = (an + 3._RKG) * (an + 4._RKG) * par2
                            Vector(k + 3) = as - (ansq - 4._RKG) * ac
                            an = an + 2._RKG
                        end do
                        ansq = an * an
                        invansq = 1._RKG / ansq
                        d(noequ) = -2._RKG * (ansq - 4._RKG) * (par22 - ansq - ansq)
                        Vector(noequ + 3) = as - (ansq - 4._RKG) * ac
                        Vector(4) = Vector(4) - 56._RKG * par2 * Vector(3)
                        ass = parint * sinpar
                        asap = (((((210._RKG * par2 - 1._RKG) * cospar - (105._RKG * par2 - 63._RKG) * ass) * invansq - (1._RKG - 15._RKG * par2) * cospar + 15._RKG * ass) * invansq - cospar + 3._RKG * ass) * invansq - cospar) * invansq
                        Vector(noequ + 3) = Vector(noequ + 3) - 2._RKG * asap * par2 * (an - 1._RKG) * (an - 2._RKG)
                        ! Solve the tridiagonal system by means of Gaussian elimination with partial pivoting.
                        call dgtsl(noequ, d1, d, d2, Vector(4), iers)
                    end if
                    do j = 1_IK, 13_IK
                        chebmom(m, 2_IK * j - 1_IK) = Vector(j)
                    end do

                    ! Compute the Chebyshev moments with respect to sine.

                    Vector(1) = 2._RKG * (sinpar - parint * cospar) / par2
                    Vector(2) = (18._RKG - 48._RKG / par2) * sinpar / par2 + (-2._RKG + 48._RKG / par2) * cospar / parint
                    ac = -24._RKG * parint * cospar
                    as = -8._RKG * sinpar
                    if (abs(parint) > 24._RKG) then ! Compute the Chebyshev moments by means of forward recursion.
                        an = 3._RKG
                        do i = 3_IK, 12_IK
                            ansq = an * an
                            Vector(i) = ((ansq - 4._RKG) * (2._RKG * (par22 - ansq - ansq) * Vector(i - 1) + as) + ac - par2 * (an + 1._RKG) * (an + 2._RKG) * Vector(i - 2)) / (par2 * (an - 1._RKG) * (an - 2._RKG))
                            an = an + 2._RKG
                        end do
                    else ! Compute the Chebyshev moments as the solutions of a boundary value problem with 1 initial value (Vector(2)) and 1 end value (computed using an asymptotic formula).
                        an = 5._RKG
                        do k = 1_IK, noeq1
                            ansq = an * an
                            d(k) = -2._RKG * (ansq - 4._RKG) * (par22 - ansq - ansq)
                            d2(k) = (an - 1._RKG) * (an - 2._RKG) * par2
                            d1(k + 1) = (an + 3._RKG) * (an + 4._RKG) * par2
                            Vector(k + 2) = ac + (ansq - 4._RKG) * as
                            an = an + 2._RKG
                        end do
                        ansq = an * an
                        d(noequ) = -2._RKG * (ansq - 4._RKG) * (par22 - ansq - ansq)
                        Vector(noequ + 2) = ac + (ansq - 4._RKG) * as
                        Vector(3) = Vector(3) - 42._RKG * par2 * Vector(2)
                        ass = parint * cospar
                        invansq = 1._RKG / ansq
                        asap = (((((105._RKG * par2 - 63._RKG) * ass + (210._RKG * par2 - 1._RKG) * sinpar) * invansq + (15._RKG * par2 - 1._RKG) * sinpar - 15._RKG * ass) * invansq - 3._RKG * ass - sinpar) * invansq - sinpar) * invansq
                        Vector(noequ + 2) = Vector(noequ + 2) - 2._RKG * asap * par2 * (an - 1._RKG) * (an - 2._RKG)
                        ! Solve the tridiagonal system by means of Gaussian elimination with partial pivoting.
                        call dgtsl(noequ, d1, d, d2, Vector(3), iers)
                    end if
                    do j = 1_IK, 12_IK
                        chebmom(m, 2 * j) = Vector(j)
                    end do
                end if
                if (currentLevel < momcom) m = currentLevel + 1_IK
                if (momcom < maxLenChebMom - 1_IK .and. momcom <= currentLevel) momcom = momcom + 1_IK

                ! Compute the coefficients of the Chebyshev expansions of degrees 12 and 24 of the function `f(x)`.

                func(1) = 0.5_RKG * f(midInterval + halfInterval)
                func(13) = f(midInterval)
                func(25) = 0.5_RKG * f(midInterval - halfInterval)
                do i = 2_IK, 12_IK
                    isym = 26_IK - i
                    func(i) = f(halfInterval * NODE(i - 1) + midInterval)
                    func(isym) = f(midInterval - halfInterval * NODE(i - 1))
                end do
                call dqcheb(x, func, cheb12, cheb24)

                ! Compute the integral and error estimates.

                result12cos = cheb12(13) * chebmom(m, 13)
                result12sin = 0._RKG
                k = 11_IK
                do j = 1_IK, 6_IK
                    result12cos = result12cos + cheb12(k) * chebmom(m, k)
                    result12sin = result12sin + cheb12(k + 1) * chebmom(m, k + 1)
                    k = k - 2_IK
                end do
                result24cos = cheb24(25) * chebmom(m, 25)
                result24sin = 0._RKG
                intAbsFunc = abs(cheb24(25))
                k = 23_IK
                do j = 1_IK, 12_IK
                    result24cos = result24cos + cheb24(k) * chebmom(m, k)
                    result24sin = result24sin + cheb24(k + 1) * chebmom(m, k + 1)
                    intAbsFunc = intAbsFunc + abs(cheb24(k)) + abs(cheb24(k + 1))
                    k = k - 2_IK
                end do
                estc = abs(result24cos - result12cos)
                ests = abs(result24sin - result12sin)
                intAbsFunc = intAbsFunc * abs(halfInterval)
                if (Integr == 2_IK) then
                    quadWF = conc * result24sin + cons * result24cos
                    abserr = abs(conc * ests) + abs(cons * estc)
                else
                    quadWF = conc * result24cos - cons * result24sin
                    abserr = abs(conc * estc) + abs(cons * ests)
                end if
            else
                call dqk15w(f, dqwgtf, Weight%omega, p2, p3, p4, Integr, lb, ub, quadWF, abserr, intAbsFunc, smoothness)
                neval = 15_IK
            end if
        end function

        function getFuncWeightedSin(x) result(funcWeightedFour)
            real(RKG)   , intent(in)    :: x
            real(RKG)                   :: funcWeightedFour
            funcWeightedFour = sin(Weight%omega * x) * getFunc(x)
        end function

        function getFuncWeightedSin(x) result(funcWeightedFour)
            real(RKG)   , intent(in)    :: x
            real(RKG)                   :: funcWeightedFour
            funcWeightedFour = sin(Weight%omega * x) * getFunc(x)
        end function

#endif

        !%%%%%%%%%%%%%%%%%%%
#elif   setChebExpan_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        integer(IK)             :: i, j
        real(RKG)               :: alam, alam1, alam2, part1, part2, part3, Vector(12)
        real(RKG)   , parameter :: NODE(11) = [(cos(i * acos(-1._RKG) / 24._RKG), i = 1, size(NODE))]

        do i = 1_IK, 12_IK
            j = 26_IK - i
            Vector(i) = func(i) - func(j)
            func(i) = func(i) + func(j)
        end do
        alam1 = Vector(1) - Vector(9)
        alam2 = NODE(6) * (Vector(3) - Vector(7) - Vector(11))
        cheb12(4) = alam1 + alam2
        cheb12(10) = alam1 - alam2
        alam1 = Vector(2) - Vector(8) - Vector(10)
        alam2 = Vector(4) - Vector(6) - Vector(12)
        alam = NODE(3) * alam1 + NODE(9) * alam2
        cheb24(4) = cheb12(4) + alam
        cheb24(22) = cheb12(4) - alam
        alam = NODE(9) * alam1 - NODE(3) * alam2
        cheb24(10) = cheb12(10) + alam
        cheb24(16) = cheb12(10) - alam
        part1 = NODE(4) * Vector(5)
        part2 = NODE(8) * Vector(9)
        part3 = NODE(6) * Vector(7)
        alam1 = Vector(1) + part1 + part2
        alam2 = NODE(2) * Vector(3) + part3 + NODE(10) * Vector(11)
        cheb12(2) = alam1 + alam2
        cheb12(12) = alam1 - alam2
        alam = NODE(1) * Vector(2) + NODE(3) * Vector(4) + NODE(5) * Vector(6) + NODE(7) * Vector(8) + NODE(9) * Vector(10) + NODE(11) * Vector(12)
        cheb24(2) = cheb12(2) + alam
        cheb24(24) = cheb12(2) - alam
        alam = NODE(11) * Vector(2) - NODE(9) * Vector(4) + NODE(7) * Vector(6) - NODE(5) * Vector(8) + NODE(3) * Vector(10) - NODE(1) * Vector(12)
        cheb24(12) = cheb12(12) + alam
        cheb24(14) = cheb12(12) - alam
        alam1 = Vector(1) - part1 + part2
        alam2 = NODE(10) * Vector(3) - part3 + NODE(2) * Vector(11)
        cheb12(6) = alam1 + alam2
        cheb12(8) = alam1 - alam2
        alam = NODE(5) * Vector(2) - NODE(9) * Vector(4) - NODE(1) * Vector(6) - NODE(11) * Vector(8) + NODE(3) * Vector(10) + NODE(7) * Vector(12)
        cheb24(6) = cheb12(6) + alam
        cheb24(20) = cheb12(6) - alam
        alam = NODE(7) * Vector(2) - NODE(3) * Vector(4) - NODE(11) * Vector(6) + NODE(1) * Vector(8) - NODE(9) * Vector(10) - NODE(5) * Vector(12)
        cheb24(8) = cheb12(8) + alam
        cheb24(18) = cheb12(8) - alam
        do i = 1_IK, 6_IK
            j = 14_IK - i
            Vector(i) = func(i) - func(j)
            func(i) = func(i) + func(j)
        end do
        alam1 = Vector(1) + NODE(8) * Vector(5)
        alam2 = NODE(4) * Vector(3)
        cheb12(3) = alam1 + alam2
        cheb12(11) = alam1 - alam2
        cheb12(7) = Vector(1) - Vector(5)
        alam = NODE(2) * Vector(2) + NODE(6) * Vector(4) + NODE(10) * Vector(6)
        cheb24(3) = cheb12(3) + alam
        cheb24(23) = cheb12(3) - alam
        alam = NODE(6) * (Vector(2) - Vector(4) - Vector(6))
        cheb24(7) = cheb12(7) + alam
        cheb24(19) = cheb12(7) - alam
        alam = NODE(10) * Vector(2) - NODE(6) * Vector(4) + NODE(2) * Vector(6)
        cheb24(11) = cheb12(11) + alam
        cheb24(15) = cheb12(11) - alam
        do i = 1_IK, 3_IK
            j = 8_IK - i
            Vector(i) = func(i) - func(j)
            func(i) = func(i) + func(j)
        end do
        cheb12(5) = Vector(1) + NODE(8) * Vector(3)
        cheb12(9) = func(1) - NODE(8) * func(3)
        alam = NODE(4) * Vector(2)
        cheb24(5) = cheb12(5) + alam
        cheb24(21) = cheb12(5) - alam
        alam = NODE(8) * func(2) - func(4)
        cheb24(9) = cheb12(9) + alam
        cheb24(17) = cheb12(9) - alam
        cheb12(1) = func(1) + func(3)
        alam = func(2) + func(4)
        cheb24(1) = cheb12(1) + alam
        cheb24(25) = cheb12(1) - alam
        cheb12(13) = Vector(1) - Vector(3)
        cheb24(13) = cheb12(13)
        alam = 1._RKG / 6._RKG
        do i = 2_IK, 12_IK
            cheb12(i) = cheb12(i) * alam
        end do
        alam = 0.5_RKG * alam
        cheb12(1) = cheb12(1) * alam
        cheb12(13) = cheb12(13) * alam
        do i = 2_IK, 24_IK
            cheb24(i) = cheb24(i) * alam
        end do
        cheb24(1) = 0.5_RKG * alam * cheb24(1)
        cheb24(25) = 0.5_RKG * alam * cheb24(25)

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  EPS_TABLE_START_OFFSET
#undef  FINALIZE_INTEGRATION
#undef  CYCLE_NEEDED
#undef  TRANS_COEF
#undef  QRULE_ARG
#undef  HELP_ARG