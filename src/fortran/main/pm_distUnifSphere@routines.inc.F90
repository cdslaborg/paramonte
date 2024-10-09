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
!>  This include file contains the implementation of procedures in [pm_distUnifSphere](@ref pm_distUnifSphere).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Set the uniform RNG
#if     RNGD_ENABLED
#define RNG
#elif   RNGF_ENABLED || RNGX_ENABLED
#define RNG rng,
#elif   getUSR_ENABLED || setUSR_ENABLED
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getUnifSphereLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     D0_ENABLED
        CHECK_ASSERTION(__LINE__, 0_IK < ndim, SK_"@getUnifSphereLogPDF(): The condition `0 < ndim` must hold. ndim = "//getStr(ndim))
        logPDF = (1 - ndim) * logRadius - getLogVolUnitBall(real(ndim, TKG)) - log(real(ndim, TKG))
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   D1_ENABLED && setUSR_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim, ndim
        real(TKG)   :: normfac
#if     AC_ENABLED
#define UNIFSPHERERAND unifSphereRand
        real(TKG)   :: unifSphereRand(size(rand, 1, IK))
#elif   DC_ENABLED
#define UNIFSPHERERAND rand
#else
#error  "Unrecognized interface."
#endif
        ! Perform runtime bound checks.
#if     AM_ENABLED
        CHECK_ASSERTION(__LINE__, size(rand, 1, IK) == size(mean, 1, IK), SK_"@setUnifSphereRand(): The condition `size(rand, 1) == size(mean)` must hold. size(rand, 1), size(mean) = "//getStr([size(rand, 1, IK), size(mean, 1, IK)]))
#elif   AC_ENABLED
        CHECK_ASSERTION(__LINE__, all(size(rand, 1, IK) == shape(chol, IK)), SK_"@setUnifSphereRand(): The condition `all(size(rand, 1) == shape(chol))` must hold. size(rand, 1), shape(chol) = "//getStr([size(rand, 1, IK), shape(chol, IK)]))
#elif   !(DM_ENABLED || DC_ENABLED)
#error  "Unrecognized interface."
#endif
        ndim = size(rand, 1, IK)
        do
            normfac = 0._TKG
            do idim = 1_IK, ndim
                call setNormRand(RNG UNIFSPHERERAND(idim))
                normfac = normfac + UNIFSPHERERAND(idim)**2
            end do
            ! Ensure the vector is not origin. Highly unlikely but possible.
            if (0._TKG < normfac) exit
        end do
        normfac = 1._TKG / sqrt(normfac)
#if     DM_ENABLED && !AC_ENABLED
        UNIFSPHERERAND = UNIFSPHERERAND * normfac ! a uniform random point on the surface of nd-sphere.
#elif   AM_ENABLED && !AC_ENABLED
        UNIFSPHERERAND = UNIFSPHERERAND * normfac + mean ! a uniform random point on the surface of nd-sphere with given mean.
#elif   AC_ENABLED
        ! Define the indexing rules.
#if     XLD_ENABLED
#define GET_INDEX(I,J)I,J
#elif   UXD_ENABLED
#define GET_INDEX(I,J)J,I
#elif   D1_ENABLED && setUSR_ENABLED
#error  "Unrecognized interface."
#endif
        ! Define the default mean.
#if     DM_ENABLED
#define MEAN_PLUS(I)
#elif   AM_ENABLED
#define MEAN_PLUS(I) mean(I) +
#else
#error  "Unrecognized interface."
#endif
        ! Separate the first to allow the possibility of adding an optional `mean`.
        rand(1 : ndim) = MEAN_PLUS(1 : ndim) chol(GET_INDEX(1 : ndim, 1)) * unifSphereRand(1) * normfac
        do idim = 2_IK, ndim
            rand(idim : ndim) = rand(idim : ndim) + chol(GET_INDEX(idim : ndim, idim)) * unifSphereRand(idim) * normfac
        end do
#else
#error  "Unrecognized interface."
#endif
#undef  UNIFSPHERERAND
#undef  MEAN_PLUS
#undef  GET_INDEX
#undef  MEAN

        !%%%%%%%%%%%%%
#elif   getUSR_ENABLED
        !%%%%%%%%%%%%%

#if     AM_ENABLED && DC_ENABLED
        call setUnifSphereRand(RNG rand, mean)
#elif   DM_ENABLED && AC_ENABLED
        call setUnifSphereRand(RNG rand, chol, subset)
#elif   AM_ENABLED && AC_ENABLED
        call setUnifSphereRand(RNG rand, mean, chol, subset)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   D2_ENABLED && setUSR_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ipnt, ndim
        ndim = size(rand, 1, IK)
        do ipnt = 1_IK, size(rand, 2, IK)
#if         DM_ENABLED && DC_ENABLED
            call setUnifSphereRand(RNG rand(1:ndim, ipnt))
#elif       AM_ENABLED && DC_ENABLED
            call setUnifSphereRand(RNG rand(1:ndim, ipnt), mean)
#elif       DM_ENABLED && AC_ENABLED
            call setUnifSphereRand(RNG rand(1:ndim, ipnt), chol, subset)
#elif       AM_ENABLED && AC_ENABLED
            call setUnifSphereRand(RNG rand(1:ndim, ipnt), mean, chol, subset)
#else
#error  "Unrecognized interface."
#endif
        end do

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  RNG