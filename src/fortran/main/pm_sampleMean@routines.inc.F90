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
!>  This file contains the implementation details of the routines under the generic interface [pm_sampleMean](@ref pm_sampleMean).
!>
!>  \author
!>  \FatemehBagheri, Saturday 4:40 PM, August 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define weight rules.
#if     WNO_ENABLED
#define INCREMENT(X,Y)
#define GET_WEIGHTED(X,W)X
#elif   WTI_ENABLED || WTR_ENABLED
#define GET_WEIGHTED(X,W)X * W
#define INCREMENT(X,Y)X = X + Y
#elif   !(getMeanMerged_ENABLED || setMeanMerged_ENABLED)
#error  "Unrecognized interface."
#endif
        ! Define weight type and kind.
#if     WTI_ENABLED
#define TYPE_OF_WEIGHT integer(IK)
#elif   WTR_ENABLED
#define TYPE_OF_WEIGHT real(TKG)
#elif   !(WNO_ENABLED || getMeanMerged_ENABLED || setMeanMerged_ENABLED)
#error  "Unrecognized interface."
#endif
        ! Define the runtime checks.
#define CHECK_LEN_MEAN(DIM) \
CHECK_ASSERTION(__LINE__, size(sample, DIM, IK) == size(mean, 1, IK), \
SK_"@setMean(): The condition `size(sample, dim, 1) == size(mean, 1)` must hold. dim, shape(sample), size(mean) = "\
//getStr([DIM, shape(sample, IK), size(mean, 1, IK)]))
#define CHECK_LEN_WEI(DIM) \
CHECK_ASSERTION(__LINE__, size(sample, DIM, IK) == size(weight, 1, IK), \
SK_"@setMean(): The condition `size(sample, dim, 1) == size(weight, 1)` must hold. dim, shape(sample), size(weight) = "\
//getStr([DIM, shape(sample, IK), size(weight, 1, IK)]))
#define CHECK_VAL_WEI \
CHECK_ASSERTION(__LINE__, all(0 <= weight), \
SK_"@setMean(): The condition `all(0. <= weight)` must hold. weight = "//getStr(weight))
#define CHECK_SUM_WEI \
CHECK_ASSERTION(__LINE__, 0 < sum(weight), \
SK_"@setMean(): The condition `0 < sum(weight)` must hold. weight = "//getStr(weight))
#define CHECK_SHAPE_SAMPLE \
CHECK_ASSERTION(__LINE__, all(0_IK < shape(sample, IK)), \
SK_"@setMean(): The condition `all(0 < shape(sample))` must hold. shape(sample) = "//getStr(shape(sample, IK)))
#if     setMean_ENABLED && DIM_ENABLED
#define CHECK_VAL_DIM \
CHECK_ASSERTION(__LINE__, 1 <= dim .and. dim <= rank(sample), \
SK_"@setMean(): The condition `1 <= dim .and. dim <= rank(sample)` must hold. dim, rank(sample) = "//getStr([integer(IK) :: dim, rank(sample)]))
#elif   setMean_ENABLED && ALL_ENABLED
#define CHECK_VAL_DIM
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getMeanMerged_ENABLED && New_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setMeanMerged(meanMerged, meanB, meanA, fracA)

        !%%%%%%%%%%%%%%%%%%%%
#elif   setMeanMerged_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        ! Define the output value for setMeanMerged.
#if     Old_ENABLED
#define TARGET_AVG meanB
#define CHECK_LEN_TARGET_MEAN
#elif   New_ENABLED
#define CHECK_LEN_TARGET_MEAN \
CHECK_ASSERTION(__LINE__, size(meanMerged, 1, IK) == size(meanB, 1, IK), SK_"@setMeanMerged(): The condition `size(meanMerged) == size(meanB)` must hold. size(meanMerged), size(meanB) = "//getStr([size(meanMerged, 1, IK), size(meanB, 1, IK)]))
#define TARGET_AVG meanMerged
#else
#error  "Unrecognized interface."
#endif
#if     D1_ENABLED
        real(TKG) :: fracB
        integer(IK) :: idim
#endif
        CHECK_ASSERTION(__LINE__, 0._TKG < fracA .and. fracA < 1._TKG, SK_"@setMeanMerged(): The condition `0 < fracA .and. fracA < 1` must hold. fracA = "//getStr(fracA))
#if     D0_ENABLED
        TARGET_AVG = fracA * meanA + (1._TKG - fracA) * meanB
#elif   D1_ENABLED
        CHECK_LEN_TARGET_MEAN
        CHECK_ASSERTION(__LINE__, size(meanA, 1, IK) == size(meanB, 1, IK), SK_"@setMeanMerged(): The condition `size(meanA) == size(meanB)` must hold. size(meanA), size(meanB) = "//getStr([size(meanA, 1, IK), size(meanB, 1, IK)]))
        fracB = 1._TKG - fracA
        do concurrent(idim = 1 : size(meanB, 1, IK))
            TARGET_AVG(idim) = fracA * meanA(idim) + fracB * meanB(idim)
        end do
#else
#error  "Unrecognized interface."
#endif
#undef  TARGET_AVG

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getMean_ENABLED && XY_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     ALL_ENABLED && WNO_ENABLED
        call setMean(mean, x, y)
#elif   ALL_ENABLED && (WTI_ENABLED || WTR_ENABLED)
        TYPE_OF_WEIGHT :: weisum
        call setMean(mean, x, y, weight, weisum)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getMean_ENABLED && (D1_ENABLED || D2_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     ALL_ENABLED && WNO_ENABLED
        call setMean(mean, sample)
#elif   DIM_ENABLED && WNO_ENABLED
        call setMean(mean, sample, dim)
#elif   ALL_ENABLED && (WTI_ENABLED || WTR_ENABLED)
        TYPE_OF_WEIGHT :: weisum
        call setMean(mean, sample, weight, weisum)
#elif   DIM_ENABLED && (WTI_ENABLED || WTR_ENABLED)
        TYPE_OF_WEIGHT :: weisum
        call setMean(mean, sample, dim, weight, weisum)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMean_ENABLED && XY_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: isam
#if     WNO_ENABLED
        real(TKG) :: weisum
        weisum = size(x, 1, IK)
#elif   WTI_ENABLED || WTR_ENABLED
        CHECK_ASSERTION(__LINE__, size(x, 1, IK) == size(weight, 1, IK), SK_"@setMean(): The condition `size(x) == size(weight)` must hold. size(x), size(weight) = "//getStr([size(x, 1, IK), size(weight, 1, IK)]))
        CHECK_SUM_WEI
        CHECK_VAL_WEI
        weisum = 0
#else
#error  "Unrecognized interface."
#endif
        mean = 0._TKG
        CHECK_ASSERTION(__LINE__, size(x, 1, IK) == size(y, 1, IK), SK_"@setMean(): The condition `size(x) == size(y)` must hold. size(x), size(y) = "//getStr([size(x, 1, IK), size(y, 1, IK)]))
        do isam = 1_IK, size(x, 1, IK)
            INCREMENT(weisum,weight(isam))
            mean(1) = mean(1) + GET_WEIGHTED(x(isam),weight(isam))
            mean(2) = mean(2) + GET_WEIGHTED(y(isam),weight(isam))
        end do
        mean = mean / weisum

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMean_ENABLED && D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: isam
#if     WNO_ENABLED
        real(TKG) :: weisum
        weisum = size(sample, 1, IK)
#elif   WTI_ENABLED || WTR_ENABLED
        CHECK_LEN_WEI(1_IK)
        CHECK_SUM_WEI
        CHECK_VAL_WEI
        weisum = 0
#else
#error  "Unrecognized interface."
#endif
        mean = 0._TKG
        CHECK_VAL_DIM
        CHECK_SHAPE_SAMPLE
        do isam = 1_IK, size(sample)
            INCREMENT(weisum,weight(isam))
            mean = mean + GET_WEIGHTED(sample(isam),weight(isam))
        end do
        mean = mean / weisum

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMean_ENABLED && D2_ENABLED && ALL_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim, jdim
        integer(IK) :: mdim, ndim
#if     WNO_ENABLED
        real(TKG) :: weisum
        weisum = size(sample, kind = IK)
#elif   WTI_ENABLED || WTR_ENABLED
        integer(IK) :: iwei
        CHECK_ASSERTION(__LINE__, size(sample, kind = IK) == size(weight, 1, IK), SK_"@setMean(): The condition `size(sample) == size(weight)` must hold. shape(sample), size(weight) = "//getStr([shape(sample, IK), shape(weight, IK)]))
        CHECK_SUM_WEI
        CHECK_VAL_WEI
        weisum = 0
        iwei = 0
#else
#error  "Unrecognized interface."
#endif
        mdim = size(sample, 1, IK)
        ndim = size(sample, 2, IK)
        mean = 0._TKG
        do jdim = 1, ndim
            do idim = 1, mdim
                INCREMENT(iwei,1_IK)
                INCREMENT(weisum,weight(iwei))
                mean = mean + GET_WEIGHTED(sample(idim, jdim),weight(iwei))
            end do
        end do
        mean = mean / weisum

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMean_ENABLED && D2_ENABLED && WNO_ENABLED && DIM_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKG) :: normfac
        integer(IK) :: ndim, nsam, idim, isam
        CHECK_VAL_DIM
        CHECK_SHAPE_SAMPLE
        CHECK_LEN_MEAN(3 - dim)
        nsam = size(sample, dim, IK)
        ndim = size(sample, 3 - dim, IK)
        normfac = 1._TKG / nsam
        if (dim == 2_IK) then
            ! attributes are along the columns.
            mean = 0._TKG
            do isam = 1, nsam
                mean = mean + sample(1 : ndim, isam)
            end do
            mean = mean * normfac
        else
            ! attributes are along the rows.
            do concurrent(idim = 1 : ndim)
                mean(idim) = sum(sample(1 : nsam, idim)) * normfac
            end do
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMean_ENABLED && D2_ENABLED && DIM_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKG) :: normfac
        integer(IK) :: ndim, nsam, idim, isam
        ndim = size(sample, 3 - dim, IK)
        nsam = size(sample, dim, IK)
        CHECK_LEN_MEAN(3 - dim)
        CHECK_SHAPE_SAMPLE
        CHECK_VAL_DIM
#if     WNO_ENABLED
        normfac = 1._TKG / nsam
#elif   WTI_ENABLED || WTR_ENABLED
        CHECK_LEN_WEI(dim)
        CHECK_SUM_WEI
        CHECK_VAL_WEI
        weisum = 0
#else
#error  "Unrecognized interface."
#endif
        if (dim == 2_IK) then ! attributes are along the columns.
            mean = 0._TKG
            do isam = 1, nsam
                INCREMENT(weisum,weight(isam))
                mean = mean + GET_WEIGHTED(sample(1 : ndim, isam),weight(isam))
            end do
#if         WTI_ENABLED || WTR_ENABLED
            normfac = 1._TKG / weisum
#endif
            mean = mean * normfac
        else ! attributes are along the rows.
#if         WNO_ENABLED
            do concurrent(idim = 1 : ndim)
                mean(idim) = sum(sample(1 : nsam, idim)) * normfac
            end do
#elif       WTI_ENABLED || WTR_ENABLED
            mean(1) = 0._TKG
            do isam = 1, nsam
                INCREMENT(weisum,weight(isam))
                mean(1) = mean(1) + GET_WEIGHTED(sample(isam, 1),weight(isam))
            end do
            normfac = 1._TKG / weisum
            mean(1) = mean(1) * normfac
            do concurrent(idim = 2 : ndim)
                ! warning: dot_product() complex-conjugates its first argument.
                mean(idim) = dot_product(weight, sample(1 : nsam, idim)) * normfac
            end do
        end if
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef CHECK_LEN_TARGET_MEAN
#undef CHECK_SHAPE_SAMPLE
#undef CHECK_LEN_MEAN
#undef TYPE_OF_WEIGHT
#undef CHECK_LEN_WEI
#undef CHECK_SUM_WEI
#undef CHECK_VAL_WEI
#undef CHECK_VAL_DIM
#undef GET_WEIGHTED
#undef INCREMENT