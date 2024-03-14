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
!>  This include file contains procedure implementation of [pm_ziggurat](@ref pm_ziggurat).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, April 25, 2015, 2:21 PM, National Institute for Fusion Studies, The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%
#if     getZig_ENABLED
        !%%%%%%%%%%%%%

        real(RKC)   :: ub, area
        real(RKC)   , parameter :: SQRT_HUGE = sqrt(huge(0._RKC))
        real(RKC)   , parameter :: lb = 10 * epsilon(0._RKC)
        CHECK_ASSERTION(__LINE__, 1_IK < nlay, SK_"@getZig(): The condition `1 <= nlay` must hold. shape(nlay) = "//getStr(nlay))
        zig(1, nlay) = 0._RKC ! The highest rectangle point (origin).
        zig(2, nlay) = getFunc(zig(1, nlay)) ! maxFunc (corresponding to the highest rectangle point)
        ! Define the upper bound of the search limit.
        ub = lb
        CHECK_ASSERTION(__LINE__, 0 <= getZ(lb), SK_"@getZig(): Internal library error occurred. The condition `0 <= getZ(lb)` must hold. lb, getZ(lb) = "//getStr([lb, getZ(lb)]))
        do
            ub = 2._RKC * ub
            if (0._RKC <= getZ(ub)) cycle ! \todo should we compare the sign change from zlb to zub here instead?
            CHECK_ASSERTION(__LINE__, ub < SQRT_HUGE, SK_"@getZig(): The specified `getFunc()` does not appear to be monotonically decreasing. The condition `ub < SQRT_HUGE` must hold. ub, SQRT_HUGE = "//getStr([ub, SQRT_HUGE]))
            exit
        end do
        zig(1, 1) = getRoot(getZ, lb, ub, abstol = abstol) ! r: first/lowest/largest rectangle-corner point.
        zig(1, 0) = area / zig(2, 1)
        zig(2, 0) = 0._RKC

    contains

        ! This function returns the difference between the first computed area (from the lowest layer) and the last computed area from the highest layer.
        impure function getZ(r) result(z)
            real(RKC)   , intent(in)    :: r
            real(RKC)                   :: z
            real(RKC)                   :: func
            integer(IK)                 :: ilay
            zig(1, 1) = r
            zig(2, 1) = getFunc(r)
            area = getZigArea(r) ! v: partition area. Store it in the first (always zero) element to save space.
            do ilay = 2, nlay - 1
                func = zig(2, ilay - 1) + area / zig(1, ilay - 1)
                if (zig(2, nlay) < func) then
                    z = SQRT_HUGE
                    abserr = z
                    return
                end if
                zig(1, ilay) = getFuncInv(func)
                zig(2, ilay) = getFunc(zig(1, ilay))
            end do
            ilay = nlay - 1_IK
            z = area - zig(1, ilay) * (zig(2, nlay) - zig(2, ilay))
            abserr = abs(z)
        end function

!        !%%%%%%%%%%%%%%%%%%
!#elif   getZig_ENABLED && 0
!        !%%%%%%%%%%%%%%%%%%
!
!        real(RKC)   :: ub, area
!        real(RKC)   , parameter :: SQRT_HUGE = sqrt(huge(0._RKC))
!        real(RKC)   , parameter :: lb = 10 * epsilon(0._RKC)
!        integer(IK) :: nlayTwice
!        nlayTwice = 2 * nlay
!        CHECK_ASSERTION(__LINE__, 1_IK < nlay, SK_"@getZig(): The condition `1 <= nlay` must hold. shape(nlay) = "//getStr(nlay))
!        zig(nlay) = 0._RKC ! The highest rectangle point (origin).
!        zig(nlayTwice) = getFunc(zig(nlay)) ! maxFunc (corresponding to the highest rectangle point)
!        ! Define the upper bound of the search limit.
!        ub = lb
!        if (getZ(lb) < 0._RKC) error stop "@getZig(): Internal error occurred. Please report to the ParaMonte library developers."
!        do
!            ub = 2._RKC * ub
!            if (getZ(ub) < 0._RKC) exit
!        end do
!        zig(1) = getRoot(getZ, lb, ub) ! r: first/lowest/largest rectangle-corner point.
!        zig(0) = area / zig(1 + nlay)
!        !zig(nlayTwice + 1) = getRoot(getZ, lb, ub) ! r: first/lowest/largest rectangle-corner point.
!        !zig(2 : nlay) = zig(2 : nlay) / zig(1 : nlay - 1)
!        !zig(1) = zig(1) * zig(1 + nlay) / area ! x1 / x0 = x1 / (area / f(r)) = x1 * f(r) / area.
!        !!zig(1 : nlay) = zig(1 : nlay) * 2._RKC**31
!
!    contains
!
!        impure function getZ(r) result(z)
!            real(RKC)   , intent(in)    :: r
!            real(RKC)                   :: z
!            real(RKC)                   :: func
!            integer(IK)                 :: ilay
!            zig(1) = r
!            zig(1 + nlay) = getFunc(r)
!            area = getZigArea(r) ! v: partition area. Store it in the first (always zero) element to save space.
!            do ilay = 2_IK, nlay - 1_IK
!                func = area / zig(ilay - 1) + zig(ilay - 1 + nlay)
!                if (zig(nlayTwice) < func) then
!                    z = SQRT_HUGE
!                    abserr = z
!                    return
!                end if
!                zig(ilay) = getFuncInv(func)
!                zig(ilay + nlay) = getFunc(zig(ilay))
!            end do
!            ilay = nlay - 1_IK
!            z = area + zig(ilay) * (zig(ilay + nlay) - 1._RKC)
!            abserr = z
!        end function
!
!        !%%%%%%%%%%%%%%%%
!#elif   getZigCRD_ENABLED
!        !%%%%%%%%%%%%%%%%
!
!        integer(IK) :: lenZig ! ilay, nlay,
!        lenZig = size(zig, 1, IK)
!        CHECK_ASSERTION(__LINE__, mod(lenZig, 2_IK) == 1_IK, SK_"@getZigCRD(): The condition `mod(size(zig), 2) == 1` must hold. shape(size(zig)) = "//getStr(lenZig))
!        crd(1,:) = zig(1 : lenZig / 2)
!        crd(2,:) = zig(1 + lenZig / 2 : lenZig - 1)
!        !nlay = lenZig / 2_IK
!        !crd(1, 1) = zig(lenZig)
!        !crd(2, 1) = zig(nlay + 1)
!        !do ilay = 2_IK, nlay
!        !    crd(1, ilay) = crd(1, ilay - 1) * zig(ilay) ! * 0.5_RKC**31
!        !    crd(2, ilay) = zig(nlay + ilay)
!        !end do

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unlayognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
