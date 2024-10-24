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
!>  This file contains procedure implementations of [pm_ellipsoid](@ref pm_ellipsoid).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Wednesday 00:57 AM, September 22, 2021, Dallas TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     (getLogVolUnitBall_ENABLED || setLogVolUnitBall_ENABLED) && Gamm_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: ndimHalf
        real(RKG), parameter :: LOG_PI = log(acos(-1._RKG))
        CHECK_ASSERTION(__LINE__, real(0, RKG) <= ndim, SK_"@setLogVolUnitBall(): The condition `0 <= ndim` must hold. ndim = "//getStr(ndim))
        ndimHalf = 0.5_RKG * ndim
        if (0._RKG < ndim) then
            logVolUnitBall = ndimHalf * LOG_PI - log_gamma(1._RKG + ndimHalf)
        else
            logVolUnitBall = 0._RKG
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (getLogVolUnitBall_ENABLED || setLogVolUnitBall_ENABLED) && Iter_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setVolUnitBall(logVolUnitBall, ndim)
        logVolUnitBall = log(logVolUnitBall)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (getVolUnitBall_ENABLED || setVolUnitBall_ENABLED) && Iter_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim, hdim, ndimint
        real(RKG), parameter :: PI = acos(-1._RKG)
        real(RKG), parameter :: FOUR_PI = PI * 4._RKG
        ndimint = ndim
        CHECK_ASSERTION(__LINE__, 0_IK <= ndimint, SK_"@setVolUnitBall(): The condition `0 <= ndim` must hold. ndim = "//getStr(ndim))
        if (0_IK < ndimint) then
            hdim = ndimint / 2_IK
            if (ndimint == 2_IK * hdim) then
                ! ndim is even
                ! ndim = 2 * hdim
                ! volUnitBall = PI^(hdim) / Factorial(hdim)
                volUnitBall = PI
                do idim = 2_IK, ndimint / 2_IK
                    volUnitBall = volUnitBall * PI / idim
                end do
            else
                ! ndim is an odd integer
                ! ndim = 2 * hdim - 1
                ! gamma(ndim / 2 + 1) = gamma(hdim + 1 / 2);
                ! gamma(hdim + 1 / 2) = sqrt(PI) * (2 * hdim)! / (4**hdim * hdim!)
                !hdim = hdim + 1_IK ! (ndim - 1) / 2
                volUnitBall = 2._RKG! / (hdim + 1_IK) ! This is to avoid an extra unnecessary division of `volUnitBall` by `PI`.
                do idim = 1_IK, hdim !hdim + 2_IK, 2_IK * hdim
                    ! volUnitBall
                    ! = PI**(hdim - 1 / 2) / gamma(hdim + 1 / 2)
                    ! = PI**(hdim + 1) * 4**hdim * hdim! / (2 * hdim)!
                    volUnitBall = volUnitBall * FOUR_PI / (hdim + idim)
                end do
                volUnitBall = volUnitBall / ndim
            end if
        else
            volUnitBall = 1._RKG
        end if

        !%%%%%%%%%%%%%%%%%%%
#elif   getLogVolEll_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, size(gramian, 1, IK) == size(gramian, 2, IK), SK_"@setLogVolEllipsoid(): The condition `size(gramian, 1) == size(gramian, 2)` must hold. shape(gramian) = "//getStr(shape(gramian, IK)))
        logVolEll = getLogVolUnitBall(real(size(gramian, 1, IK), RKG)) + getMatDetSqrtLog(gramian)

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getCountMemberEll_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the center argument.
#if     Org_ENABLED
#define CENTER
#elif   Cen_ENABLED
#define CENTER center,
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: ndim, ipnt
        countMemberEll = 0_IK
        ndim = size(point, 1, IK)
        do ipnt = 1_IK, size(point, 2, IK)
            if  ( & ! LCOV_EXCL_LINE
#if             Sph_ENABLED
                isMemberEll(radius, CENTER point(1:ndim, ipnt)) & ! LCOV_EXCL_LINE
#elif           Ell_ENABLED
                isMemberEll(invGram, CENTER point(1:ndim, ipnt)) & ! LCOV_EXCL_LINE
#else
#error          "Unrecognized interface."
#endif
                ) countMemberEll = countMemberEll + 1_IK
        end do

        !%%%%%%%%%%%%%%%%%%
#elif   isMemberEll_ENABLED
        !%%%%%%%%%%%%%%%%%%

        integer(IK) :: ndim
#if     Cen_ENABLED
        real(RKG)   :: normedPoint(size(point, 1, IK))
        CHECK_ASSERTION(__LINE__, size(point, 1, IK) == size(center, 1, IK), SK_"@isMemberEll(): The condition `size(point) == size(center)` must hold. size(point), size(center) = "//getStr([size(point, 1, IK), size(center, 1, IK)]))
#elif   !Org_ENABLED
#error  "Unrecognized interface."
#endif
#if     Sph_ENABLED
        CHECK_ASSERTION(__LINE__, 0._RKG < radius, SK_"@isMemberEll(): The condition `0. < radius` must hold. radius = "//getStr(radius))
#elif   Ell_ENABLED
        CHECK_ASSERTION(__LINE__, all(size(point, 1, IK) == shape(invGram, IK)), SK_"@isMemberEll(): The condition `all(size(point, 1) == shape(invGram))` must hold. size(point, 1), shape(invGram) = "//getStr([size(point, 1, IK), shape(invGram, IK)]))
#endif
        ndim = size(point, 1, IK)
#if     Org_ENABLED
#define NORMEDPOINT point(1 : ndim)
#elif   Cen_ENABLED
#define NORMEDPOINT normedPoint
        normedPoint(1 : ndim) = point(1 : ndim) - center
#else
#error  "Unrecognized interface."
#endif
#if     Sph_ENABLED
        isMember = logical(dot_product(NORMEDPOINT, NORMEDPOINT) <= radius, LK)
#elif   Ell_ENABLED
        isMember = logical(dot_product(NORMEDPOINT, matmul(invGram, NORMEDPOINT)) <= 1._RKG, LK)
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

#undef  NORMEDPOINT
#undef  CENTER