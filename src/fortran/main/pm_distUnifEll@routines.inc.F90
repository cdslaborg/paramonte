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
!>  This include file contains the implementation of procedures in [pm_distUnifEll](@ref pm_distUnifEll).
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
#elif   getMUR_ENABLED || setMUR_ENABLED
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%
#if     getUnifEllLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

#if     D0_ENABLED
        CHECK_ASSERTION(__LINE__, 0_IK < ndim, SK_"@getUnifEllLogPDF(): The condition `0 < ndim` must hold. ndim = "//getStr(ndim))
        logPDF = -ndim * logChoDia - getLogVolUnitBall(real(ndim, RKG))
#elif   D1_ENABLED
        logPDF = -sum(logChoDia) - getLogVolUnitBall(real(size(logChoDia, 1, IK), RKG))
#elif   D2_ENABLED && AIP_ENABLED
        integer(IK) :: info
        real(RKG) :: chol(size(gramian, 1, IK), size(gramian, 2, IK))
        CHECK_ASSERTION(__LINE__, size(gramian, 1, IK) == size(gramian, 2, IK), SK_"@getUnifEllLogPDF(): The condition `size(gramian, 1) == size(gramian, 2)` must hold. shape(gramian) = "//getStr(shape(gramian,IK)))
        call setMatCopy(chol, rdpack, gramian, rdpack, uppDia, doff = 0_IK)
        call setMatDetSqrtLog(chol, uppDia, logPDF, info, chol, transHerm)
        if (info /= 0_IK) error stop "The Cholesky factorization of Gramian failed. Graming is not positive definite."
        logPDF = -(logPDF + getLogVolUnitBall(real(size(gramian, 1, IK), RKG)))
        !logPDF = -getLogVolUnitBall(real(size(gramian, 1, IK), RKG)) - getMatMulTraceLog(chol)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   D1_ENABLED && setMUR_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim, ndim
        real(RKG)   :: unifrnd, ndimInv
        real(RKG)   :: sumSqUnifBallRand
#if     AC_ENABLED
#define UNIFBALLRAND unifBallRand
        real(RKG)   :: unifBallRand(size(rand, 1, IK))
#elif   DC_ENABLED
#define UNIFBALLRAND rand
#else
#error  "Unrecognized interface."
#endif
        ! Perform runtime bound checks.
#if     AM_ENABLED
        CHECK_ASSERTION(__LINE__, size(rand, 1, IK) == size(mean, 1, IK), SK_"@setUnifEllRand(): The condition `size(rand, 1) == size(mean, 1)` must hold. size(rand, 1), size(mean) = "//getStr([size(rand, 1, IK), size(mean, 1, IK)]))
#endif
#if     AC_ENABLED
        CHECK_ASSERTION(__LINE__, all(size(rand, 1, IK) == shape(chol, IK)), SK_"@setUnifEllRand(): The condition `all(size(rand, 1) == shape(chol))` must hold. size(rand, 1), shape(chol) = "//getStr([size(rand, 1, IK), shape(chol, IK)]))
        ! Define the indexing rules.
#if     XLD_ENABLED
#define GET_INDEX(I,J)I,J
#elif   UXD_ENABLED
#define GET_INDEX(I,J)J,I
#elif   D1_ENABLED && setMUR_ENABLED
#error  "Unrecognized interface."
#endif
#endif
        ! Define the default mean.
#if     DM_ENABLED
#define MEAN_PLUS(I)
#elif   AM_ENABLED
#define MEAN_PLUS(I) mean(I) +
#else
#error  "Unrecognized interface."
#endif
        ndim = size(rand, kind = IK)
        ndimInv = 1._RKG / real(ndim, RKG)
        do
            sumSqUnifBallRand = 0._RKG
            do idim = 1_IK, ndim
                call setNormRand(RNG UNIFBALLRAND(idim))
                sumSqUnifBallRand = sumSqUnifBallRand + UNIFBALLRAND(idim)**2
            end do
            ! Ensure the vector is not origin. Highly unlikely but possible.
            if (0._RKG < sumSqUnifBallRand) exit
        end do
#if     RNGD_ENABLED || RNGF_ENABLED
        call random_number(unifrnd)
#elif   RNGX_ENABLED
        call setUnifRand(rng, unifrnd)
#endif
        unifrnd = unifrnd**ndimInv / sqrt(sumSqUnifBallRand)
#if     AC_ENABLED
        ! Separate the first to allow the possibility of adding an optional `mean`.
        rand(1 : ndim) = MEAN_PLUS(1 : ndim) chol(GET_INDEX(1 : ndim, 1)) * unifBallRand(1) * unifrnd
        do idim = 2_IK, ndim
            rand(idim : ndim) = rand(idim : ndim) + chol(GET_INDEX(idim : ndim, idim)) * unifBallRand(idim) * unifrnd
        end do
#else
        UNIFBALLRAND = MEAN_PLUS(1 : ndim) UNIFBALLRAND * unifrnd ! a uniform random point from inside of an ndim-sphere.
#endif
#undef  UNIFBALLRAND
#undef  MEAN_PLUS
#undef  GET_INDEX
#undef  MEAN

        !%%%%%%%%%%%%%
#elif   getMUR_ENABLED
        !%%%%%%%%%%%%%

#if     AM_ENABLED && DC_ENABLED
        call setUnifEllRand(RNG rand, mean)
#elif   DM_ENABLED && AC_ENABLED
        call setUnifEllRand(RNG rand, chol, subset)
#elif   AM_ENABLED && AC_ENABLED
        call setUnifEllRand(RNG rand, mean, chol, subset)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   D2_ENABLED && setMUR_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ipnt, ndim
        ndim = size(rand, 1, IK)
        do ipnt = 1_IK, size(rand, 2, IK)
#if         DM_ENABLED && DC_ENABLED
            call setUnifEllRand(RNG rand(1:ndim, ipnt))
#elif       AM_ENABLED && DC_ENABLED
            call setUnifEllRand(RNG rand(1:ndim, ipnt), mean)
#elif       DM_ENABLED && AC_ENABLED
            call setUnifEllRand(RNG rand(1:ndim, ipnt), chol, subset)
#elif       AM_ENABLED && AC_ENABLED
            call setUnifEllRand(RNG rand(1:ndim, ipnt), mean, chol, subset)
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