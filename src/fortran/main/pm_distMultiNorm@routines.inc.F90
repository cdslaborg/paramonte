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
!>  This include file contains procedure implementation of [pm_distMultiNorm](@ref pm_distMultiNorm).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Set the uniform RNG
#if     RNGD_ENABLED
#define RNG
#elif   RNGF_ENABLED || RNGX_ENABLED
#define RNG rng,
#elif   getMNR_ENABLED || setMNR_ENABLED
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getMultiNormLogPDFNF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG), parameter :: LOG_INVERSE_SQRT_TWO_PI = -0.5_RKG * log(2 * acos(-1._RKG))
#if     I_ENABLED
        integer(IK) :: ndim
        real(RKG) :: logSqrtDetInvCov
        logSqrtDetInvCov = getMatDetSqrtLog(invCov, uppDia)
        ndim = size(invCov, 1, IK)
#elif   IF_ENABLED
        integer(IK) :: ndim
        real(RKG)   :: logSqrtDetInvCov, chol(size(invCov, 1, IK), size(invCov, 1, IK))
        CHECK_ASSERTION(__LINE__, size(invCov, 1, IK) == size(invCov, 2, IK), SK_"@getMultiNormLogPDFNF(): The condition `size(invCov, 1) == size(invCov, 2)` must hold. shape(invCov) = "//getStr(shape(invCov, IK)))
        call setMatDetSqrtLog(invCov, uppDia, logSqrtDetInvCov, info, chol, nothing)
        ndim = size(invCov, 1, IK)
        if (info /= 0_IK) return
#elif   ND_ENABLED
        CHECK_ASSERTION(__LINE__, 0_IK < ndim, SK_"@getMultiNormLogPDFNF(): The input `ndim` must be positive. ndim = "//getStr(ndim))
#else
#error  "Unrecognized interface."
#endif
        logPDFNF = ndim * LOG_INVERSE_SQRT_TWO_PI + logSqrtDetInvCov

        !%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getMultiNormLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Set the normalization factor for different interfaces.
#if     DDD_ENABLED || MDD_ENABLED
#define LOGNORMFAC getMultiNormLogPDFNF(size(X,1,IK), 0._RKG)
#elif   DID_ENABLED || MID_ENABLED
#define LOGNORMFAC getMultiNormLogPDFNF(invCov)
#elif   DDN_ENABLED || MDN_ENABLED || DIN_ENABLED || MIN_ENABLED
#define LOGNORMFAC logPDFNF
#else
#error  "Unrecognized interface."
#endif
        ! Set the Mahalanobis distance.
#if     DDD_ENABLED || DDN_ENABLED
#define MAHAL_SQ sum(X**2)
#elif   DID_ENABLED || DIN_ENABLED
#define MAHAL_SQ getDisMahalSq(X, invCov)
        CHECK_ASSERTION(__LINE__, all(shape(invCov, IK) == size(X, 1, IK)), SK_"@getMultiNormLogPDF(): The condition `all(shape(invCov) == size(X))` must hold. shape(invCov), size(X) = "//getStr([shape(invCov, IK), size(X, 1, IK)]))
#elif   MID_ENABLED || MIN_ENABLED
#define MAHAL_SQ getDisMahalSq(X, invCov, mean)
        CHECK_ASSERTION(__LINE__, all(shape(invCov, IK) == size(X, 1, IK)), SK_"@getMultiNormLogPDF(): The condition `all(shape(invCov) == size(X))` must hold. shape(invCov), size(X) = "//getStr([shape(invCov, IK), size(X, 1, IK)]))
        CHECK_ASSERTION(__LINE__, size(mean, 1, IK) == size(X, 1, IK), SK_"@getMultiNormLogPDF(): The condition `size(mean) == size(X)` must hold. size(mean), size(X) = "//getStr([size(mean, 1, IK), size(X, 1, IK)]))
#elif   MDD_ENABLED || MDN_ENABLED
#define MAHAL_SQ sum((X - mean)**2)
        CHECK_ASSERTION(__LINE__, size(mean, 1, IK) == size(X, 1, IK), SK_"@getMultiNormLogPDF(): The condition `size(mean) == size(X)` must hold. size(mean), size(X) = "//getStr([size(mean, 1, IK), size(X, 1, IK)]))
#else
#error  "Unrecognized interface."
#endif
        logPDF = LOGNORMFAC - 0.5_RKG * MAHAL_SQ

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   D1_ENABLED && (getMNR_ENABLED || setMNR_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ndim
        ndim = size(rand, kind = IK)
#if     AM_ENABLED
        CHECK_ASSERTION(__LINE__, size(rand, 1, IK) == size(mean, 1, IK), SK_"@setMultiNormRand(): The condition `size(rand, 1) == size(mean, 1)` must hold. shape(rand), shape(mean) = "//getStr([shape(rand, IK), shape(mean, IK)]))
#elif   !DM_ENABLED
#error  "Unrecognized interface."
#endif
        ! Generate random MVN.
#if     DC_ENABLED && DM_ENABLED
        call setNormRand(RNG rand)
#elif   DC_ENABLED && AM_ENABLED
        call setNormRand(RNG rand)
        rand = rand + mean
#elif   AC_ENABLED && (AM_ENABLED || DM_ENABLED)
        CHECK_ASSERTION(__LINE__, all(size(rand, 1, IK) == shape(chol, IK)), SK_"@setMultiNormRand(): The condition `all(size(rand, 1) == shape(chol))` must hold. size(rand, 1), shape(chol) = "//getStr([size(rand, 1, IK), shape(chol, IK)]))
        ! Define the indexing rules.
#if     XLD_ENABLED
#define GET_INDEX(I,J)I,J
#elif   UXD_ENABLED
#define GET_INDEX(I,J)J,I
#elif   D1_ENABLED && setMUR_ENABLED
#error  "Unrecognized interface."
#endif
        ! Define the default mean.
#if     DM_ENABLED
#define XPLUS(X)
#elif   AM_ENABLED
#define XPLUS(X)X +
#else
#error  "Unrecognized interface."
#endif
        block
            integer(IK) :: idim
            real(RKG) :: normrnd
            ! Separate the first to allow the possibility of adding an optional `mean`.
            call setNormRand(RNG normrnd)
            rand(1 : ndim) = XPLUS(mean) chol(GET_INDEX(1 : ndim, 1)) * normrnd
            do idim = 2_IK, ndim
                call setNormRand(RNG normrnd)
                rand(idim : ndim) = rand(idim : ndim) + chol(GET_INDEX(idim : ndim, idim)) * normrnd
            end do
        end block
#undef  GET_INDEX
#undef  XPLUS
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   D2_ENABLED && getMNR_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     AM_ENABLED && DC_ENABLED
        call setMultiNormRand(RNG rand, mean)
#elif   DM_ENABLED && AC_ENABLED
        call setMultiNormRand(RNG rand, chol, subset)
#elif   AM_ENABLED && AC_ENABLED
        call setMultiNormRand(RNG rand, mean, chol, subset)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   D2_ENABLED && setMNR_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ipnt, ndim
        ndim = size(rand, 1, IK)
        do ipnt = 1_IK, size(rand, 2, IK)
#if         DM_ENABLED && DC_ENABLED
            call setMultiNormRand(RNG rand(1:ndim, ipnt))
#elif       AM_ENABLED && DC_ENABLED
            call setMultiNormRand(RNG rand(1:ndim, ipnt), mean)
#elif       DM_ENABLED && AC_ENABLED
            call setMultiNormRand(RNG rand(1:ndim, ipnt), chol, subset)
#elif       AM_ENABLED && AC_ENABLED
            call setMultiNormRand(RNG rand(1:ndim, ipnt), mean, chol, subset)
#else
#error  "Unrecognized interface."
#endif
        end do

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  LOGNORMFAC
#undef  MAHAL_SQ
#undef  RNG