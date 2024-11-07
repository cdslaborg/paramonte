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
!>  This include file contains procedure implementation of the generic interfaces of [pm_polation](@ref pm_polation).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CHECK_CRDX_LEN(PROC) \
CHECK_ASSERTION(__LINE__, size(crdx, 1, IK) == size(func, 1, IK), PROC//SK_": The condition `size(crdx) == size(func)` must hold. size(crdx), size(func) = "//getStr([size(crdx, 1, IK), size(func, 1, IK)]))
#define CHECK_QBOUNDED(PROC) \
CHECK_ASSERTION(__LINE__, minval(crdx, 1) <= queryx .and. queryx <= maxval(crdx, 1), PROC//SK_": The condition `minval(crdx) <= queryx .and. queryx <= maxval(crdx)` must hold. minval(crdx), queryx, maxval(crdx) = "//getStr([minval(crdx), queryx, maxval(crdx)]))
        ! Define the procedure name.
#if     getExtrap_ENABLED || setExtrap_ENABLED
#define SET_POLATION setExtrap
#define POLATION extrap
#define POLSTR SK_"extrap"
        character(*, SK), parameter :: PROC_NAME = "@setExtrap()"
#elif   getInterp_ENABLED || setInterp_ENABLED
#define SET_POLATION setInterp
#define POLATION interp
#define POLSTR SK_"interp"
        character(*, SK), parameter :: PROC_NAME = "@setInterp()"
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     (getExtrap_ENABLED || getInterp_ENABLED) && ND1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call SET_POLATION(method, crdx, func, queryx, POLATION)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (setExtrap_ENABLED || setInterp_ENABLED) && ND1_ENABLED && QD1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ique
        CHECK_ASSERTION(__LINE__, size(queryx, 1, IK) == size(POLATION, 1, IK), PROC_NAME//SK_": The condition `size(queryx) == size("//POLSTR//SK_")` must hold. shape(queryx), shape("//POLSTR//SK_") = "//getStr([shape(queryx, IK), shape(POLATION, IK)]))
#if     MNPLE_ENABLED
        CHECK_ASSERTION(__LINE__, size(queryx, 1, IK) == size(relerr, 1, IK), PROC_NAME//SK_": The condition `size(queryx) == size(relerr)` must hold. shape(queryx), shape(relerr) = "//getStr([shape(queryx, IK), shape(relerr, IK)]))
#endif
        do ique = 1, size(queryx, 1, IK)
#if         MNPLD_ENABLED || PWLN_ENABLED || MEAN_ENABLED || NEAR_ENABLED || NEXT_ENABLED || PREV_ENABLED
            call SET_POLATION(method, crdx, func, queryx(ique), POLATION(ique))
#elif       MNPLE_ENABLED
            call SET_POLATION(method, crdx, func, queryx(ique), POLATION(ique), relerr(ique))
#else
#error      "Unrecognized interface."
#endif
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (setExtrap_ENABLED || setInterp_ENABLED) && ND1_ENABLED && (MNPLD_ENABLED || MNPLE_ENABLED) && QD0_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     MNPLD_ENABLED
        real(RKG) :: relerr
#elif   !MNPLE_ENABLED
#error  "Unrecognized interface."
#endif
        real(RKG) :: diff, difX, correctionC(size(crdx, kind = IK)), correctionD(size(crdx, kind = IK))
        integer(IK) :: i, level, ns
        integer(IK) :: lenx
        lenx = size(crdx, kind = IK)
        CHECK_CRDX_LEN(PROC_NAME)
#if     setInterp_ENABLED
        CHECK_QBOUNDED(PROC_NAME)
#endif
        ns = 1_IK
        ! Find the nearest neighbor.
        diff = abs(queryx - crdx(1))
        do i = 1_IK, lenx
            difX = abs(queryx - crdx(i))
            if (difX < diff) then
                ns = i
                diff = difX
            endif
            correctionC(i) = func(i)
            correctionD(i) = func(i)
        end do
        POLATION = func(ns)
        ns = ns - 1_IK
        do level = 1_IK, lenx - 1_IK
            do i = 1_IK, lenx - level
                diff = (correctionC(i + 1) - correctionD(i)) / (crdx(i) - crdx(i + level))
                correctionD(i) = (crdx(i + level) - queryx) * diff
                correctionC(i) = (crdx(i) - queryx) * diff
            end do
            if (2 * ns < lenx - level)then
                relerr = correctionC(ns + 1)
            else
                relerr = correctionD(ns)
                ns = ns - 1_IK
            endif
            POLATION = POLATION + relerr
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (setExtrap_ENABLED || setInterp_ENABLED) && ND1_ENABLED && PWLN_ENABLED && QD0_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ibin
        real(RKG) :: const, slope
        CHECK_CRDX_LEN(PROC_NAME)
        CHECK_ASSERTION(__LINE__, isAscendingAll(crdx), PROC_NAME//SK_": The condition `isAscendingAll(crdx)` must hold. crdx = "//getStr(crdx))
#if     setInterp_ENABLED
        CHECK_QBOUNDED(PROC_NAME)
        ibin = getBin(crdx, queryx)
#elif   setExtrap_ENABLED
        ibin = size(crdx, 1, IK)
        if (crdx(ibin) < queryx) then
            ibin = ibin - 1_IK
        elseif (queryx < crdx(1)) then
            ibin = 1_IK
        else
            ibin = getBin(crdx, queryx)
        end if
#else
#error  "Unrecognized interface."
#endif
        if (crdx(ibin) /= queryx) then
            const = 1._RKG / (crdx(ibin) - crdx(ibin + 1))
            slope = -const * func(ibin + 1)
            const = +const * func(ibin)
            POLATION = const * (queryx - crdx(ibin + 1)) + slope * (queryx - crdx(ibin))
        else
            POLATION = func(ibin)
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (setExtrap_ENABLED || setInterp_ENABLED) && ND1_ENABLED && (MEAN_ENABLED || NEAR_ENABLED || NEXT_ENABLED || PREV_ENABLED) && QD0_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ibin
        CHECK_CRDX_LEN(PROC_NAME)
        CHECK_ASSERTION(__LINE__, isAscending(crdx), PROC_NAME//SK_": The condition `isAscending(crdx)` must hold. crdx = "//getStr(crdx))
#if     setInterp_ENABLED
        CHECK_QBOUNDED(PROC_NAME)
#elif   setExtrap_ENABLED
        if (crdx(size(crdx, 1, IK)) < queryx) then
            POLATION = func(size(crdx, 1, IK))
            return
        elseif (queryx < crdx(1)) then
            POLATION = func(1)
            return
        end if
#else
#error  "Unrecognized interface."
#endif
        ibin = getBin(crdx, queryx)
        if (crdx(ibin) /= queryx) then
            ! ibin cannot be the last element of `crdx`, unless `queryx` is out of the range covered by `crdx`.
#if         MEAN_ENABLED
            POLATION = .5_RKG * (func(ibin) + func(ibin + 1))
#elif       NEAR_ENABLED
            if (queryx - crdx(ibin) < crdx(ibin + 1) - queryx) then
                POLATION = func(ibin)
            else
                POLATION = func(ibin + 1)
            end if
#elif       NEXT_ENABLED
            POLATION = func(ibin + 1)
#elif       PREV_ENABLED
            POLATION = func(ibin)
#else
#error      "Unrecognized interface."
#endif
        else
            POLATION = func(ibin)
        end if
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  SET_POLATION
#undef  POLATION
#undef  POLSTR