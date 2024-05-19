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
!>  This file contains the implementation details of the routines of [pm_sampleVar](@ref pm_sampleVar).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Saturday 4:40 PM, August 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the looping ranges for the corresponding matrix subsets.
#if     (setCov_ENABLED || setCovMean_ENABLED || setCovMerged_ENABLED || setCovMeanMerged_ENABLED || setCovUpdated_ENABLED || setCovMeanUpdated_ENABLED) && (XLD_ENABLED || XLD_XLD_ENABLED || XLD_UXD_ENABLED)
        ! Start from the beginning of the lower-triangle.
#define OFF_RANGE(I,J,OFFSET)I + OFFSET, J
#define ROW_RANGE(I,J,K)J + 1_IK, K
#define COL_RANGE(I,J)I, J
#define FIRST 1
#elif   (setCov_ENABLED || setCovMean_ENABLED || setCovMerged_ENABLED || setCovMeanMerged_ENABLED || setCovUpdated_ENABLED || setCovMeanUpdated_ENABLED) && (UXD_ENABLED || UXD_UXD_ENABLED || UXD_XLD_ENABLED)
        ! Start from the end of the upper-triangle.
#define OFF_RANGE(I,J,OFFSET)J - OFFSET, I, -1_IK
#define ROW_RANGE(I,J,K)J - 1_IK, I, -1_IK
#define COL_RANGE(I,J)J, I, -1_IK
#define FIRST ndim
#elif   setCovMerged_ENABLED || (setCov_ENABLED && !(XY_ENABLED || CorStd_ENABLED))
#error  "Unrecognized interface."
#endif
        ! Set the conjugation rule.
#if     CK_ENABLED
#define GET_RE(X)X%re
#define SET_CONJG(X)X = conjg(X)
#define GET_CONJG(X)conjg(X)
#define TYPE_OF_SAMPLE complex(TKG)
#define GET_ABSQ(X)(real(X)**2 + aimag(X)**2)
#define GET_PROD(X,Y)(X%re * Y%re + X%im * Y%im)
#elif   RK_ENABLED
#define GET_RE(X)X
#define SET_CONJG(X)
#define GET_CONJG(X)X
#define TYPE_OF_SAMPLE real(TKG)
#define GET_PROD(X,Y)X * Y
#define GET_ABSQ(X)X**2
#else
#error  "Unrecognized interface."
#endif
        ! Set the shifting rule.
#if     setCov_ENABLED && Org_ENABLED
#define GET_SHIFTED(X,Y)X
#elif   setCov_ENABLED && Avg_ENABLED
#define GET_SHIFTED(X,Y)(X - Y)
#elif   setCov_ENABLED && !CorStd_ENABLED
#error  "Unrecognized interface."
#endif
        ! Set the weighting rule.
#if     (setCov_ENABLED || setCovMean_ENABLED) && WNO_ENABLED
#define GET_WEIGHTED(X,Y)X
#elif   (setCov_ENABLED || setCovMean_ENABLED) && (WTI_ENABLED || WTR_ENABLED)
#define GET_WEIGHTED(X,Y)X * Y
#elif   setCovMean_ENABLED || (setCov_ENABLED && !CorStd_ENABLED)
#error  "Unrecognized interface."
#endif
        ! Define weight type and kind and ZEROW.
#if     WTI_ENABLED
#define TYPE_OF_WEIGHT integer(IK)
#define ZEROW 0_IK
#elif   WTR_ENABLED || WNO_ENABLED
#define TYPE_OF_WEIGHT real(TKG)
#define ZEROW 0._TKG
#elif   (getCov_ENABLED || setCov_ENABLED) && !(WNO_ENABLED || CorStd_ENABLED)
#error  "Unrecognized interface."
#endif
        ! Define runtime sanity checks.
#define CHECK_VAL_NSAM(PROC,DIM) \
CHECK_ASSERTION(__LINE__, 1 < size(sample, DIM, IK), \
PROC//SK_": The condition `1 < size(sample, dim)` must hold. dim, shape(sample) = "//getStr([DIM, shape(sample, IK)]))
#define CHECK_VAL_NDIM(PROC,DIM) \
CHECK_ASSERTION(__LINE__, 0 < size(sample, DIM, IK), \
PROC//SK_": The condition `0 < size(sample, dim)` must hold. dim, shape(sample) = "//getStr([DIM, shape(sample, IK)]))
#define CHECK_VAL_DIM(PROC) \
CHECK_ASSERTION(__LINE__, 1 <= dim .and. dim <= rank(sample), \
PROC//SK_": The condition `1 <= dim .and. dim <= rank(sample)` must hold. dim, rank(sample) = "\
//getStr([integer(IK) :: dim, rank(sample)]))
#define CHECK_SUM_WEI(PROC) \
CHECK_ASSERTION(__LINE__, ZEROW < sum(weight), \
PROC//SK_": The condition `0 < sum(weight)` must hold. weight = "//getStr(weight))
#define CHECK_VAL_WEI(PROC) \
CHECK_ASSERTION(__LINE__, all(0._TKG <= weight), \
PROC//SK_": The condition `all(0. <= weight)` must hold. weight = "//getStr(weight))
#define CHECK_SHAPE_COV(PROC) \
CHECK_ASSERTION(__LINE__, all(size(sample, 3 - dim, IK) == shape(cov, IK)), \
PROC//SK_": The condition `all(size(sample, 3 - dim) == shape(cov))` must hold. dim, size(sample, 3 - dim), shape(cov) = "\
//getStr([dim, size(sample, 3 - dim, IK), shape(cov, IK)]))
#define CHECK_VAL_MEANG(PROC,dim) \
CHECK_ASSERTION(__LINE__, all([minval(sample, dim) <= meang .and. meang <= maxval(sample, dim)]), \
PROC//SK_": The condition `all([minval(sample, dim) <= meang .and. meang <= maxval(sample, dim)])` must hold. dim, minval(sample, dim), meang, maxval(sample, dim) = "//\
getStr(dim)//SK_"; "//getStr(minval(sample, dim))//SK_"; "//getStr(meang)//SK_"; "//getStr(maxval(sample, dim)))
#define CHECK_LEN_MEANG(PROC) \
CHECK_ASSERTION(__LINE__, size(sample, 3 - dim, IK) == size(meang, 1, IK), \
PROC//SK_": The condition `size(sample, 3 - dim) == size(meang)` must hold. dim, size(sample, 3 - dim), size(meang, 1) = "\
//getStr([dim, size(sample, 3 - dim, IK), size(meang, 1, IK)]))
#define CHECK_LEN_MEAN(PROC) \
CHECK_ASSERTION(__LINE__, size(sample, 3 - dim, IK) == size(mean, 1, IK), \
PROC//SK_": The condition `size(sample, 3 - dim) == size(mean)` must hold. dim, size(sample, 3 - dim), size(mean, 1) = "\
//getStr([dim, size(sample, 3 - dim, IK), size(mean, 1, IK)]))
#define CHECK_LEN_WEI(PROC,DIM) \
CHECK_ASSERTION(__LINE__, size(sample, DIM, IK) == size(weight, 1, IK), \
PROC//SK_": The condition `size(sample, dim) == size(weight)` must hold. dim, size(sample, dim), size(weight, 1) = "\
//getStr([DIM, size(sample, DIM, IK), size(weight, 1, IK)]))
#define CHECK_WEISUM(PROC) \
CHECK_ASSERTION(__LINE__, abs(weisum - sum(weight)) < 1000 * epsilon(0._TKG), \
PROC//SK_": The condition `0 < sum(weight)` must hold. weisum, sum(weight) = "//getStr([weisum, sum(weight)]))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getCovMerged_ENABLED && New_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim
        call setCovMerged(cov, covB, covA, meanDiff, fracA, uppDia)
        !do concurrent(idim = 2 : size(cov, 1, IK))
        do idim = 2, size(cov, 1, IK)
            cov(idim, 1 : idim - 1) = GET_CONJG(cov(1 : idim - 1, idim))
        end do

        !%%%%%%%%%%%%%%%%%%%%
#elif   setCovUpdated_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim, jdim, ndim
        real(TKG) :: fracB, fracAB
        fracB = 1._TKG - fracA
        ! Define the output value for setCovMerged.
        CHECK_ASSERTION(__LINE__, 0._TKG < fracA .and. fracA < 1._TKG, SK_"@setCovMerged(): The condition `0 < fracA .and. fracA < 1` must hold. fracA = "//getStr(fracA))
        CHECK_ASSERTION(__LINE__, all(size(meanDiff, 1, IK) == shape(covA, IK)), SK_"@setCovMerged(): The condition `all(size(meanDiff) == shape(covA))` must hold. size(meanDiff), shape(covA) = "//getStr([size(meanDiff, 1, IK), shape(covA, IK)]))
        fracAB = fracA * fracB
        ndim = size(covA, 1, IK)
        do jdim = COL_RANGE(1_IK,ndim)
            covA(jdim, jdim) = fracA * GET_RE(covA(jdim, jdim)) + fracAB * GET_ABSQ(meanDiff(jdim))
            do idim = ROW_RANGE(1_IK,jdim,ndim)
                covA(idim, jdim) = fracA * covA(idim, jdim) + fracAB * (meanDiff(idim) * GET_CONJG(meanDiff(jdim)))
            end do
        end do

        !%%%%%%%%%%%%%%%%%%%
#elif   setCovMerged_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim, jdim, ndim
        real(TKG) :: fracB, fracAB
        fracB = 1._TKG - fracA
        ! Define the output value for setCovMerged.
#if     Old_ENABLED
#define cov covB
#elif   New_ENABLED
        CHECK_ASSERTION(__LINE__, all(shape(cov, IK) == shape(covA, IK)), SK_"@setCovMerged(): The condition `all(shape(cov) == shape(covA))` must hold. shape(cov), shape(covA) = "//getStr([shape(cov, IK), shape(covA, IK)]))
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, 0._TKG < fracA .and. fracA < 1._TKG, SK_"@setCovMerged(): The condition `0 < fracA .and. fracA < 1` must hold. fracA = "//getStr(fracA))
        CHECK_ASSERTION(__LINE__, all(shape(covB, IK) == shape(covA, IK)), SK_"@setCovMerged(): The condition `all(shape(covB, IK) == shape(covA, IK))` must hold. shape(covB), shape(covA) = "//getStr([shape(covB, IK), shape(covA, IK)]))
        CHECK_ASSERTION(__LINE__, all(size(meanDiff, 1, IK) == shape(covA, IK)), SK_"@setCovMerged(): The condition `all(size(meanDiff) == shape(covA))` must hold. size(meanDiff), shape(covA) = "//getStr([size(meanDiff, 1, IK), shape(covA, IK)]))
        fracAB = fracA * fracB
        ndim = size(cov, 1, IK)
        do jdim = COL_RANGE(1_IK,ndim)
            cov(jdim, jdim) = fracB * GET_RE(covB(jdim, jdim)) + fracA * GET_RE(covA(jdim, jdim)) + fracAB * GET_ABSQ(meanDiff(jdim))
            do idim = ROW_RANGE(1_IK,jdim,ndim)
                cov(idim, jdim) = fracB * covB(idim, jdim) + fracA * covA(idim, jdim) + fracAB * (meanDiff(idim) * GET_CONJG(meanDiff(jdim)))
            end do
        end do
#undef  cov

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCovMeanUpdated_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_OF_SAMPLE :: idiff, temp
        integer(IK) :: idim, jdim, ndim
        real(TKG) :: fracB, fracAB
        fracB = 1._TKG - fracA
        ! Define the output value for setCovMerged.
        CHECK_ASSERTION(__LINE__, 0._TKG < fracA .and. fracA < 1._TKG, SK_"@setCovMerged(): The condition `0 < fracA .and. fracA < 1` must hold. fracA = "//getStr(fracA))
        CHECK_ASSERTION(__LINE__, all(size(meanA, 1, IK) == shape(covA, IK)), SK_"@setCovMerged(): The condition `all(size(meanA) == shape(covA))` must hold. size(meanA), shape(covA) = "//getStr([size(meanA, 1, IK), shape(covA, IK)]))
        CHECK_ASSERTION(__LINE__, all(size(meanB, 1, IK) == shape(covA, IK)), SK_"@setCovMerged(): The condition `all(size(meanB) == shape(covA))` must hold. size(meanB), shape(covA) = "//getStr([size(meanB, 1, IK), shape(covA, IK)]))
        fracAB = fracA * fracB
        ndim = size(covA, 1, IK)
        do jdim = COL_RANGE(1_IK,ndim)
            temp = meanA(jdim) - meanB(jdim)
            covA(jdim, jdim) = fracA * GET_RE(covA(jdim, jdim)) + fracAB * GET_ABSQ(temp)
            SET_CONJG(temp)
            do idim = ROW_RANGE(1_IK,jdim,ndim)
                idiff = meanA(idim) - meanB(idim)
                covA(idim, jdim) = fracA * covA(idim, jdim) + fracAB * idiff * temp
            end do
        end do
        !do concurrent(idim = 1 : size(covA, 1, IK))
        do idim = 1, size(covA, 1, IK)
            meanA(idim) = fracA * meanA(idim) + fracB * meanB(idim)
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCovMeanMerged_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_OF_SAMPLE :: idiff, temp
        integer(IK) :: idim, jdim, ndim
        real(TKG) :: fracB, fracAB
        fracB = 1._TKG - fracA
        ! Define the output value for setCovMerged.
#if     Old_ENABLED
#define mean meanB
#define cov covB
#elif   New_ENABLED
        CHECK_ASSERTION(__LINE__, all(shape(cov, IK) == shape(covA, IK)), SK_"@setCovMerged(): The condition `all(shape(cov) == shape(covA))` must hold. shape(cov), shape(covA) = "//getStr([shape(cov, IK), shape(covA, IK)]))
        CHECK_ASSERTION(__LINE__, all(size(mean, 1, IK) == shape(covA, IK)), SK_"@setCovMerged(): The condition `all(size(mean) == shape(covA))` must hold. size(mean), shape(covA) = "//getStr([size(mean, 1, IK), shape(covA, IK)]))
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, 0._TKG < fracA .and. fracA < 1._TKG, SK_"@setCovMerged(): The condition `0 < fracA .and. fracA < 1` must hold. fracA = "//getStr(fracA))
        CHECK_ASSERTION(__LINE__, all(size(meanA, 1, IK) == shape(covA, IK)), SK_"@setCovMerged(): The condition `all(size(meanA) == shape(covA))` must hold. size(meanA), shape(covA) = "//getStr([size(meanA, 1, IK), shape(covA, IK)]))
        CHECK_ASSERTION(__LINE__, all(size(meanB, 1, IK) == shape(covA, IK)), SK_"@setCovMerged(): The condition `all(size(meanB) == shape(covA))` must hold. size(meanB), shape(covA) = "//getStr([size(meanB, 1, IK), shape(covA, IK)]))
        CHECK_ASSERTION(__LINE__, all(shape(covB, IK) == shape(covA, IK)), SK_"@setCovMerged(): The condition `all(shape(covB, IK) == shape(covA, IK))` must hold. shape(covB), shape(covA) = "//getStr([shape(covB, IK), shape(covA, IK)]))
        fracAB = fracA * fracB
        ndim = size(cov, 1, IK)
        do jdim = COL_RANGE(1_IK,ndim)
            temp = meanA(jdim) - meanB(jdim)
            cov(jdim, jdim) = fracB * GET_RE(covB(jdim, jdim)) + fracA * GET_RE(covA(jdim, jdim)) + fracAB * GET_ABSQ(temp)
            SET_CONJG(temp)
            do idim = ROW_RANGE(1_IK,jdim,ndim)
                idiff = meanA(idim) - meanB(idim)
                cov(idim, jdim) = fracB * covB(idim, jdim) + fracA * covA(idim, jdim) + fracAB * idiff * temp
            end do
        end do
        !do concurrent(idim = 1 : size(covA, 1, IK))
        do idim = 1, size(covA, 1, IK)
            mean(idim) = fracA * meanA(idim) + fracB * meanB(idim)
        end do
#undef  mean
#undef  cov

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getCov_ENABLED && CorStd_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim
        call setCov(cov, uppDia, cor, subsetr, std)
        !do concurrent(idim = 2 : size(cov, 1, IK))
        do idim = 2, size(cov, 1, IK)
            cov(idim, 1 : idim - 1) = GET_CONJG(cov(1 : idim - 1, idim))
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getCov_ENABLED && CorStd_ENABLED && 0
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     ULD_ULD_ENABLED
        integer(IK) :: idim
        type(uppDia_type), parameter :: subset = uppDia_type(), subsetr = uppDia_type()
#elif   !(UXD_UXD_ENABLED || UXD_XLD_ENABLED || XLD_UXD_ENABLED || XLD_XLD_ENABLED)
#error  "Unrecognized interface."
#endif
        call setCov(cov, subset, cor, subsetr, std)
#if     ULD_ULD_ENABLED
        !do concurrent(idim = 2 : size(cov, 1, IK))
        do idim = 2, size(cov, 1, IK)
            cov(idim, 1 : idim - 1) = GET_CONJG(cov(1 : idim - 1, idim))
        end do
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCov_ENABLED && CorStd_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim, jdim, ndim
        ndim = size(cov, 1, IK)
        CHECK_ASSERTION(__LINE__, ndim == size(cov, 2, IK), SK_"@setCov(): The condition `size(cov, 1) == size(cov, 2)` must hold. shape(cov) = "//getStr(shape(cov, IK)))
        CHECK_ASSERTION(__LINE__, ndim == size(std, 1, IK), SK_"@setCov(): The condition `size(cov, 1) == size(std)` must hold. size(cov, 1), size(std) = "//getStr([ndim, size(std, 1, IK)]))
        CHECK_ASSERTION(__LINE__, all(ndim == shape(cor, IK)), SK_"@setCov(): The condition `all(size(cov, 1) == shape(cor))` must hold. size(cov, 1), shape(cor) = "//getStr([ndim, shape(cor, IK)]))
        CHECK_ASSERTION(__LINE__, all(0._TKG < std), SK_"@setCov(): The condition `all(0. < std)` must hold. std = "//getStr(std))
        if (ndim == 0_IK) return
        do jdim = COL_RANGE(1,ndim)
            cov(jdim, jdim) = std(jdim)**2
            do idim = ROW_RANGE(1,jdim,ndim)
#if             UXD_UXD_ENABLED || XLD_XLD_ENABLED
                cov(idim, jdim) = std(idim) * std(jdim) * cor(idim, jdim)
#elif           UXD_XLD_ENABLED || XLD_UXD_ENABLED
                cov(idim, jdim) = std(idim) * std(jdim) * GET_CONJG(cor(jdim, idim))
#else
#error          "Unrecognized interface."
#endif
            end do
        end do
!#if     UXD_UXD_ENABLED
!        do jdim = 1_IK, ndim
!            do idim = 1_IK, jdim - 1_IK
!                cov(idim, jdim) = cor(idim, jdim) * std(idim) * std(jdim)
!            end do
!            cov(jdim, jdim) = std(jdim)**2
!        end do
!#elif   UXD_XLD_ENABLED
!        cov(1, 1) = std(1)**2
!        do jdim = 2_IK, ndim
!            do idim = 1_IK, jdim - 1_IK
!                cov(idim, jdim) = GET_CONJG(cor(jdim, idim)) * std(idim) * std(jdim)
!            end do
!            cov(jdim, jdim) = std(jdim)**2
!        end do
!#elif   XLD_UXD_ENABLED
!        cov(1, 1) = std(1)**2
!        do jdim = 2_IK, ndim
!            do idim = 1_IK, jdim - 1_IK
!                cov(jdim, idim) = GET_CONJG(cor(idim, jdim)) * std(idim) * std(jdim)
!            end do
!            cov(jdim, jdim) = std(jdim)**2
!        end do
!#elif   XLD_XLD_ENABLED
!        cov(1, 1) = std(1)**2
!        do jdim = 2_IK, ndim
!            do idim = 1_IK, jdim - 1_IK
!                cov(jdim, idim) = cor(jdim, idim) * std(idim) * std(jdim)
!            end do
!            cov(jdim, jdim) = std(jdim)**2
!        end do
!#else
!#error  "Unrecognized interface."
!#endif

        !%%%%%%%%%%%%%
#elif   getCov_ENABLED
        !%%%%%%%%%%%%%

        type(uppDia_type), parameter :: subset = uppDia_type()
        integer(IK) :: ndim, idim, nsam
        real(TKG) :: normfac
#if     WNO_ENABLED
#define WEIGHT_ARGS
#elif   WTI_ENABLED || WTR_ENABLED
#define WEIGHT_ARGS, weight, weisum
        TYPE_OF_WEIGHT :: weisum
#else
#error  "Unrecognized interface."
#endif
#if     XY_ENABLED
        TYPE_OF_SAMPLE :: mean(2)
        call setCovMean(cov, mean, x, y WEIGHT_ARGS, [x(1), y(1)])
        nsam = size(x, 1, IK)
#elif   ULD_ENABLED
        TYPE_OF_SAMPLE, dimension(size(sample, 3 - dim, IK)) :: mean, meang
        nsam = size(sample, dim, IK)
        if (dim == 1_IK) then
            meang = sample(1,:)
        else
            meang = sample(:,1)
        end if
        call setCovMean(cov, subset, mean, sample, dim WEIGHT_ARGS, meang)
#else
#error  "Unrecognized interface."
#endif
        ! Symmetrize.
        ndim = size(cov, 1, IK)
        !do concurrent(idim = 1 : ndim)
        do idim = 1, ndim
            cov(idim, 1 : idim - 1) = GET_CONJG(cov(1 : idim - 1, idim))
        end do
        ! Correct if requested.
        if (present(correction)) then
            CHECK_ASSERTION(__LINE__, same_type_as(correction, fweight) .or. same_type_as(correction, rweight), SK_"@getCov(): The condition `same_type_as(correction, fweight) .or. same_type_as(correction, rweight)` must hold.")
#if         WNO_ENABLED
            normfac = getVarCorrection(real(nsam, TKG))
#elif       (WTI_ENABLED || WTR_ENABLED)
            if (same_type_as(correction, fweight)) then
                normfac = getVarCorrection(real(weisum, TKG))
            elseif (same_type_as(correction, rweight)) then
                normfac = getVarCorrection(real(weisum, TKG), real(sum(weight**2), TKG))
            end if
#else
#error      "Unrecognized interface."
#endif
            cov = cov * normfac
        end if
#undef  WEIGHT_ARGS
#undef  GET_SAM

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCov_ENABLED && XY_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: isam, nsam
#if     WTI_ENABLED || WTR_ENABLED
        TYPE_OF_SAMPLE :: temp
#endif
        TYPE_OF_SAMPLE :: cxy
        real(TKG) :: cxx, cyy
        real(TKG) :: normFac
#if     Avg_ENABLED
        TYPE_OF_SAMPLE :: difx, dify
        CHECK_ASSERTION(__LINE__, size(mean, 1, IK) == 2_IK, SK_"@setCov(): The condition `size(mean) == 2` must hold. size(mean) = "//getStr(size(mean, 1, IK)))
#elif   !Org_ENABLED
#error  "Unrecognized interface."
#endif
        nsam = size(x, 1, IK)
        CHECK_ASSERTION(__LINE__, size(x, 1, IK) == size(y, 1, IK), SK_"@setCov(): The condition `size(x) == size(y)` must hold. size(x), size(y) = "//getStr([size(x, 1, IK), size(y, 1, IK)]))
        CHECK_ASSERTION(__LINE__, all(shape(cov, IK) == 2_IK), SK_"@setCov(): The condition `all(shape(cov) == 2)` must hold. shape(cov) = "//getStr(shape(cov, IK)))
        CHECK_ASSERTION(__LINE__, 1_IK < nsam, SK_"@setCov(): The condition `1 < size(x)` must hold. size(x) = "//getStr(nsam))
        CHECK_ASSERTION(__LINE__, any(x(1) /= x), SK_"@setCov(): The condition `any(x(1) /= x)` must hold. x = "//getStr(x))
        CHECK_ASSERTION(__LINE__, any(y(1) /= y), SK_"@setCov(): The condition `any(y(1) /= y)` must hold. y = "//getStr(y))
#if     WNO_ENABLED
        normFac = 1._TKG / real(nsam, TKG)
#elif   WTI_ENABLED || WTR_ENABLED
        CHECK_ASSERTION(__LINE__, size(x, 1, IK) == size(weight, 1, IK), SK_"@setCov(): The condition `size(x) == size(weight)` must hold. size(x), size(weight) = "//getStr([size(x, 1, IK), size(weight, 1, IK)]))
        CHECK_ASSERTION(__LINE__, all(0 <= weight), SK_"@setCov(): The condition `all(0 <= weight)` must hold. pack(weight, weight < 0) = "//getStr(pack(weight, weight < 0)))
        normFac = 1._TKG / real(weisum, TKG)
#else
#error  "Unrecognized interface."
#endif
        cxx = 0._TKG
        cxy = 0._TKG
        cyy = 0._TKG
        do isam = 1, nsam
#if         Avg_ENABLED && (WTI_ENABLED || WTR_ENABLED)
            difx = x(isam) - mean(1)
            dify = y(isam) - mean(2)
            temp = dify * weight(isam)
            cyy = cyy + GET_PROD(dify,temp)
            cxy = cxy + difx * GET_CONJG(temp)
            cxx = cxx + GET_ABSQ(difx) * weight(isam)
#elif       Org_ENABLED && (WTI_ENABLED || WTR_ENABLED)
            temp = y(isam) * weight(isam)
            cyy = cyy + GET_PROD(y(isam),temp)
            cxy = cxy + x(isam) * GET_CONJG(temp)
            cxx = cxx + GET_ABSQ(x(isam)) * weight(isam)
#elif       Avg_ENABLED && WNO_ENABLED
            difx = x(isam) - mean(1)
            dify = y(isam) - mean(2)
            cyy = cyy + GET_ABSQ(dify)
            cxy = cxy + difx * GET_CONJG(dify)
            cxx = cxx + GET_ABSQ(difx)
#elif       Org_ENABLED && WNO_ENABLED
            cyy = cyy + GET_ABSQ(y(isam))
            cxy = cxy + x(isam) * GET_CONJG(y(isam))
            cxx = cxx + GET_ABSQ(x(isam))
#else
#error      "Unrecognized interface."
#endif
        end do
        cov(1,1) = normFac * cxx
        cov(2,2) = normFac * cyy
        cov(1,2) = normFac * cxy
        cov(2,1) = GET_CONJG(cov(1,2))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCov_ENABLED && (UXD_ENABLED || XLD_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKG) :: normFac
        TYPE_OF_SAMPLE :: temp
        integer(IK) :: idim, jdim, isam, ndim, nsam
#if     Avg_ENABLED
#if     WTI_ENABLED || WTR_ENABLED
        TYPE_OF_SAMPLE :: diff
#elif   !WNO_ENABLED
#error  "Unrecognized interface."
#endif
        CHECK_LEN_MEAN(SK_"@setCov()")
#elif   !Org_ENABLED
#error  "Unrecognized interface."
#endif
        ndim = size(sample, 3 - dim, IK)
        nsam = size(sample, dim, IK)
        CHECK_VAL_DIM(SK_"@setCov()")
        CHECK_SHAPE_COV(SK_"@setCov()")
        CHECK_VAL_NSAM(SK_"@setCov()",dim)
        CHECK_VAL_NDIM(SK_"@setCov()",3 - dim)
#if     WNO_ENABLED
        normFac = 1._TKG / real(nsam, TKG)
#elif   WTI_ENABLED || WTR_ENABLED
        CHECK_WEISUM(SK_"@setCov()")
        CHECK_SUM_WEI(SK_"@setCov()")
        CHECK_VAL_WEI(SK_"@setCov()")
        CHECK_LEN_WEI(SK_"@setCov()",dim)
        normFac = 1._TKG / real(weisum, TKG)
#else
#error  "Unrecognized interface."
#endif
        ! Compute the cov.
        if (dim == 2_IK) then
            do jdim = COL_RANGE(1_IK,ndim)
                CHECK_ASSERTION(__LINE__, any(sample(jdim,1) /= sample(jdim,:)), SK_"@setCov(): The condition `any(sample(jdim,1) /= sample(jdim,:))` must hold. jdim = "//getStr(jdim)//SK_", sample(jdim,:) = "//getStr(reshape(sample(jdim,:),[size(sample,2),1])))
                cov(jdim, jdim) = 0._TKG
                do idim = ROW_RANGE(1_IK,jdim,ndim)
                    cov(idim, jdim) = 0._TKG
                end do
            end do
            do isam = 1, nsam
                do jdim = COL_RANGE(1_IK,ndim)
#if                 Avg_ENABLED && (WTI_ENABLED || WTR_ENABLED)
                    diff = sample(jdim, isam) - mean(jdim)
                    temp = diff * weight(isam)
                    GET_RE(cov(jdim, jdim)) = GET_RE(cov(jdim, jdim)) + GET_PROD(diff,temp)
#elif               Org_ENABLED && (WTI_ENABLED || WTR_ENABLED)
                    temp = sample(jdim, isam) * weight(isam)
                    GET_RE(cov(jdim, jdim)) = GET_RE(cov(jdim, jdim)) + GET_PROD(sample(jdim, isam),temp)
#elif               Avg_ENABLED && WNO_ENABLED
                    temp = sample(jdim, isam) - mean(jdim)
                    GET_RE(cov(jdim, jdim)) = GET_RE(cov(jdim, jdim)) + GET_PROD(temp,temp)
#elif               Org_ENABLED && WNO_ENABLED
                    temp = sample(jdim, isam)
                    GET_RE(cov(jdim, jdim)) = GET_RE(cov(jdim, jdim)) + GET_PROD(sample(jdim, isam),temp)
#else
#error              "Unrecognized interface."
#endif
                    SET_CONJG(temp)
                    do idim = ROW_RANGE(1_IK,jdim,ndim)
                        cov(idim, jdim) = cov(idim, jdim) + GET_SHIFTED(sample(idim, isam),mean(idim)) * temp
                    end do
                end do
            end do
            do jdim = COL_RANGE(1_IK,ndim)
                GET_RE(cov(jdim, jdim)) = GET_RE(cov(jdim, jdim)) * normFac
                do idim = ROW_RANGE(1_IK,jdim,ndim)
                    cov(idim, jdim) = cov(idim, jdim) * normFac
                end do
            end do
        else ! dim = 1_IK
            do jdim = COL_RANGE(1_IK,ndim)
                CHECK_ASSERTION(__LINE__, any(sample(1,jdim) /= sample(:,jdim)), SK_"@setCov(): The condition `any(sample(1,jdim) /= sample(:,jdim))` must hold. jdim = "//getStr(jdim)//SK_", sample(:,jdim) = "//getStr(sample(:,jdim)))
                cov(jdim, jdim) = 0._TKG
                do isam = 1, nsam
#if                 Avg_ENABLED
                    temp = sample(isam, jdim) - mean(jdim)
                    GET_RE(cov(jdim, jdim)) = GET_RE(cov(jdim, jdim)) + GET_WEIGHTED(GET_ABSQ(temp),weight(isam))
#elif               Org_ENABLED
                    GET_RE(cov(jdim, jdim)) = GET_RE(cov(jdim, jdim)) + GET_WEIGHTED(GET_ABSQ(sample(isam, jdim)),weight(isam))
#endif
                end do
                cov(jdim, jdim) = cov(jdim, jdim) * normFac
                do idim = ROW_RANGE(1_IK,jdim,ndim)
                    cov(idim, jdim) = 0._TKG
                    do isam = 1, nsam
                        cov(idim, jdim) = cov(idim, jdim) + GET_WEIGHTED(GET_SHIFTED(sample(isam, idim),mean(idim)) * GET_CONJG(GET_SHIFTED(sample(isam, jdim),mean(jdim))),weight(isam))
                    end do
                    cov(idim, jdim) = cov(idim, jdim) * normFac
                end do
            end do
        end if ! dim

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCovMean_ENABLED && XY_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: isam
        TYPE_OF_SAMPLE :: cxy, difx, dify
        real(TKG) :: cxx, cyy
        real(TKG) :: normFac
#if     WNO_ENABLED
        real(TKG) :: weisum
        weisum = size(x, 1, IK)
#elif   WTI_ENABLED || WTR_ENABLED
        TYPE_OF_SAMPLE :: temp
        CHECK_ASSERTION(__LINE__, size(x, 1, IK) == size(weight, 1, IK), SK_"@setCovMean(): The condition `size(x) == size(weight)` must hold. size(x), size(weight) = "//getStr([size(x, 1, IK), size(weight, 1, IK)]))
        CHECK_ASSERTION(__LINE__, all(0 <= weight), SK_"@setCovMean(): The condition `all(0 <= weight)` must hold. pack(weight, weight < 0) = "//getStr(pack(weight, weight < 0)))
        weisum = ZEROW
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, size(mean, 1, IK) == 2_IK, SK_"@setCovMean(): The condition `size(mean) == 2` must hold. size(mean) = "//getStr(size(mean, 1, IK)))
        CHECK_ASSERTION(__LINE__, size(meang, 1, IK) == 2_IK, SK_"@setCovMean(): The condition `size(meang) == 2` must hold. size(meang) = "//getStr(size(meang, 1, IK)))
        CHECK_ASSERTION(__LINE__, all(shape(cov, IK) == 2_IK), SK_"@setCovMean(): The condition `all(shape(cov) == 2)` must hold. shape(cov) = "//getStr(shape(cov, IK)))
        CHECK_ASSERTION(__LINE__, size(x, 1, IK) == size(y, 1, IK), SK_"@setCovMean(): The condition `size(x) == size(y)` must hold. size(x), size(y) = "//getStr([size(x, 1, IK), size(y, 1, IK)]))
        CHECK_ASSERTION(__LINE__, all([minval(x) <= meang(1) .and. meang(2) <= maxval(y)]), SK_"@setCovMean(): The condition `all([minval(x) <= meang(1) .and. meang(2) <= maxval(y)])` must hold. minval(x), meang, maxval(y) = "//getStr([minval(x), meang, maxval(y)]))
        CHECK_ASSERTION(__LINE__, 1_IK < size(x, 1, IK), SK_"@setCovMean(): The condition `1 < size(x)` must hold. size(x) = "//getStr(size(x, 1, IK)))
        CHECK_ASSERTION(__LINE__, any(x(1) /= x), SK_"@setCovMean(): The condition `any(x(1) /= x)` must hold. x = "//getStr(x))
        CHECK_ASSERTION(__LINE__, any(y(1) /= y), SK_"@setCovMean(): The condition `any(y(1) /= y)` must hold. y = "//getStr(y))
        cxx = 0._TKG
        cyy = 0._TKG
        cxy = 0._TKG
        mean = 0._TKG
        do isam = 1, size(x, 1, IK)
            difx = x(isam) - meang(1)
            dify = y(isam) - meang(2)
#if         WTI_ENABLED || WTR_ENABLED
            weisum = weisum + weight(isam)
            temp = difx * weight(isam)
            mean(1) = mean(1) + temp
            cxx = cxx + GET_PROD(difx,temp)
            temp = dify * weight(isam)
            mean(2) = mean(2) + temp
            cyy = cyy + GET_PROD(dify,temp)
            cxy = cxy + difx * GET_CONJG(temp)
#elif       WNO_ENABLED
            mean(1) = mean(1) + difx
            cxx = cxx + GET_ABSQ(difx)
            mean(2) = mean(2) + dify
            cyy = cyy + GET_ABSQ(dify)
            cxy = cxy + difx * GET_CONJG(dify)
#else
#error      "Unrecognized interface."
#endif
        end do
        normFac = 1._TKG / weisum
        mean = mean * normFac
        cov(1,1) = (cxx - GET_ABSQ(mean(1)) * weisum) * normFac
        cov(2,2) = (cyy - GET_ABSQ(mean(2)) * weisum) * normFac
        cov(1,2) = (cxy - mean(1) * GET_CONJG(mean(2)) * weisum) * normFac
        cov(2,1) = GET_CONJG(cov(1,2))
        mean = mean + meang

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCovMean_ENABLED && (UXD_ENABLED || XLD_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKG) :: normFac
        TYPE_OF_SAMPLE :: diff
        integer(IK) :: idim, jdim, isam, ndim, nsam
#if     WNO_ENABLED
#define SET(X)
#define temp diff
#define INCREMENT(X,Y)
        real(TKG) :: weisum
        weisum = size(sample, dim, IK)
#elif   WTI_ENABLED || WTR_ENABLED
#define INCREMENT(X,Y)X = X + Y
#define SET(X)X
        TYPE_OF_SAMPLE :: temp
        CHECK_LEN_WEI(SK_"@setCovMean()",dim)
        CHECK_VAL_WEI(SK_"@setCovMean()")
        weisum = ZEROW
#else
#error  "Unrecognized interface."
#endif
        nsam = size(sample, dim, IK)
        ndim = size(sample, 3 - dim, IK)
        CHECK_LEN_MEAN(SK_"@setCovMean()")
        CHECK_VAL_NSAM(SK_"@setCovMean()",dim)
        CHECK_VAL_MEANG(SK_"@setCovMean()",dim)
        CHECK_LEN_MEANG(SK_"@setCovMean()")
        CHECK_SHAPE_COV(SK_"@setCovMean()")
        CHECK_VAL_DIM(SK_"@setCovMean()")
        CHECK_VAL_NSAM(SK_"@setCovMean()",dim)
        CHECK_VAL_NDIM(SK_"@setCovMean()",3 - dim)
        ! Compute the cov.
        if (dim == 2_IK) then
            do jdim = COL_RANGE(1_IK,ndim)
                CHECK_ASSERTION(__LINE__, any(sample(jdim,1) /= sample(jdim,:)), SK_"@setCovMean(): The condition `any(sample(jdim,1) /= sample(jdim,:))` must hold. jdim = "//getStr(jdim)//SK_", sample(jdim,:) = "//getStr(sample(jdim,:)))
                mean(jdim) = 0._TKG
                cov(jdim, jdim) = 0._TKG
                do idim = ROW_RANGE(1_IK,jdim,ndim)
                    cov(idim, jdim) = 0._TKG
                end do
            end do
            do isam = 1, nsam
                INCREMENT(weisum,weight(isam))
                do jdim = COL_RANGE(1_IK,ndim)
                    diff = sample(jdim, isam) - meang(jdim)
                    SET(temp = diff * weight(isam))
                    mean(jdim) = mean(jdim) + temp
                    GET_RE(cov(jdim, jdim)) = GET_RE(cov(jdim, jdim)) + GET_PROD(diff,temp)
                    SET_CONJG(temp)
                    do idim = ROW_RANGE(1_IK,jdim,ndim)
                        cov(idim, jdim) = cov(idim, jdim) + (sample(idim, isam) - meang(idim)) * temp
                    end do
                end do
            end do
            normFac = 1._TKG / weisum
            mean = mean * normFac
            do jdim = COL_RANGE(1_IK,ndim)
                GET_RE(cov(jdim, jdim)) = (GET_RE(cov(jdim, jdim)) - GET_ABSQ(mean(jdim)) * weisum) * normFac
                temp = GET_CONJG(mean(jdim))
                do idim = ROW_RANGE(1_IK,jdim,ndim)
                    cov(idim, jdim) = (cov(idim, jdim) - mean(idim) * temp * weisum) * normFac
                end do
            end do
            mean = mean + meang
        else ! dim == 1
            ! Compute `weisum` and the first element of `mean` in the first round.
            CHECK_ASSERTION(__LINE__, any(sample(1,FIRST) /= sample(:,FIRST)), SK_"@setCovMean(): The condition `any(sample(1,idim) /= sample(:,idim))` must hold. idim = "//getStr(FIRST)//SK_", sample(:,idim) = "//getStr(sample(:,FIRST)))
            cov(FIRST, FIRST) = 0._TKG
            mean(FIRST) = 0._TKG
            do isam = 1, nsam
                INCREMENT(weisum,weight(isam))
                diff = sample(isam, FIRST) - meang(FIRST)
                SET(temp = diff * weight(isam))
                GET_RE(cov(FIRST, FIRST)) = GET_RE(cov(FIRST, FIRST)) + GET_PROD(diff,temp)
                mean(FIRST) = mean(FIRST) + temp
            end do
            normFac = 1._TKG / weisum
            mean(FIRST) = mean(FIRST) * normFac
            GET_RE(cov(FIRST, FIRST)) = (GET_RE(cov(FIRST, FIRST)) - GET_ABSQ(mean(FIRST)) * weisum) * normFac
            ! Compute the rest of the `mean` elements in the second round.
            do idim = OFF_RANGE(1_IK,ndim,1_IK)
                CHECK_ASSERTION(__LINE__, any(sample(1,idim) /= sample(:,idim)), SK_"@setCovMean(): The condition `any(sample(1,idim) /= sample(:,idim))` must hold. idim = "//getStr(idim)//SK_", sample(:,idim) = "//getStr(sample(:,idim)))
                cov(idim, FIRST) = 0._TKG
                mean(idim) = 0._TKG
                do isam = 1, nsam
                    diff = sample(isam, idim) - meang(idim)
                    SET(temp = diff * weight(isam))
                    cov(idim, FIRST) = cov(idim, FIRST) + temp * GET_CONJG((sample(isam, FIRST) - meang(FIRST)))
                    mean(idim) = mean(idim) + temp
                end do
                mean(idim) = mean(idim) * normFac
                cov(idim, FIRST) = (cov(idim, FIRST) - mean(idim) * GET_CONJG(mean(FIRST)) * weisum) * normFac
                mean(idim) = mean(idim) + meang(idim)
            end do
            mean(FIRST) = mean(FIRST) + meang(FIRST) ! This normalization must be done right here and not any sooner.
            ! Now use the computed mean to calculate the rest of the covariance matrix using the normal algorithm.
            do jdim = OFF_RANGE(1_IK,ndim,1_IK)
                cov(jdim, jdim) = 0._TKG
                do isam = 1, nsam
                    diff = sample(isam, jdim) - mean(jdim)
                    GET_RE(cov(jdim, jdim)) = GET_RE(cov(jdim, jdim)) + GET_WEIGHTED(GET_ABSQ(diff),weight(isam))
                end do
                GET_RE(cov(jdim, jdim)) = GET_RE(cov(jdim, jdim)) * normFac
                do idim = ROW_RANGE(1_IK,jdim,ndim)
                    cov(idim, jdim) = 0._TKG
                    do isam = 1, nsam
                        cov(idim, jdim) = cov(idim, jdim) + GET_WEIGHTED((sample(isam, idim) - mean(idim)) * GET_CONJG((sample(isam, jdim) - mean(jdim))),weight(isam))
                    end do
                    cov(idim, jdim) = cov(idim, jdim) * normFac
                end do
            end do
        end if
#undef  INCREMENT
#undef  temp
#undef  SET

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCov_ENABLED && Wei_ENABLED && Prob_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#error  "Unrecognized interface. This is a relic of the past from the era of frequentist ignorance."
!       ! Define the looping ranges for the corresponding matrix subsets.
!if     (setCov_ENABLED || setCovMerged_ENABLED) && UXD_ENABLED
!define COL_RANGE 2_IK, ndim
!define ROW_RANGE 1_IK, jdim
!define FIRST 1_IK
!elif   (setCov_ENABLED || setCovMerged_ENABLED) && XLD_ENABLED
!define COL_RANGE ndim - 1_IK, 1_IK, -1_IK
!define ROW_RANGE jdim, ndim
!define FIRST ndim
!elif   setCovMerged_ENABLED || (setCov_ENABLED && !(XY_ENABLED || CorStd_ENABLED))
!error  "Unrecognized interface."
!endif
!        integer(IK) :: zeroWeightCount
!        integer(IK) :: idim, jdim, isam,ndim,nsam
!        real(TKG) :: normFac
!#if     Def_ENABLED
!        TYPE_OF_WEIGHT :: weisum
!        real(TKG), allocatable :: mean(:)
!#elif   Avg_ENABLED
!        CHECK_LEN_MEAN(SK_"@setCov()")
!#else
!#error  "Unrecognized interface."
!#endif
!        CHECK_VAL_DIM(SK_"@setCov()")
!        CHECK_SHAPE_COV(SK_"@setCov()")
!        CHECK_SUM_WEI(SK_"@setCov()")
!        CHECK_VAL_WEI(SK_"@setCov()")
!        CHECK_LEN_WEI(SK_"@setCov()",dim)
!        CHECK_VAL_NSAM(SK_"@setCov()",dim)
!        CHECK_VAL_NDIM(SK_"@setCov()"),3 - dim)
!        nsam = size(sample, dim, IK)
!        ndim = size(sample, 3 - dim, IK)
!        zeroWeightCount = 0_IK
!        ! Compute the cov.
!        if (dim == 2_IK) then
!            do jdim = 1, ndim; do idim = ROW_RANGE; cov(idim, jdim) = 0._TKG; end do; end do;
!#if         Def_ENABLED
!            if (shifted) then
!                ! probability-weight correction and shifted.
!                weisum = ZEROW
!                do isam = 1, nsam
!                    weisum = weisum + weight(isam)
!                    if (weight(isam) == 0._TKG) zeroWeightCount = zeroWeightCount + 1
!                    do jdim = 1, ndim
!                        do idim = ROW_RANGE
!                            cov(idim, jdim) = cov(idim, jdim) + sample(idim, isam) * sample(jdim, isam) * weight(isam)**2
!                        end do
!                    end do
!                end do
!                normFac = real(nsam - zeroWeightCount, TKG) / (real(nsam - zeroWeightCount - 1, TKG) * weisum**2)
!                do jdim = 1, ndim; do idim = ROW_RANGE; cov(idim, jdim) = cov(idim, jdim) * normFac; end do; end do;
!                return
!            end if
!            allocate(mean(ndim))
!            call setMean(mean, sample, dim, weight, weisum)
!#endif
!            ! probability-weight correction and not shifted.
!            do isam = 1, nsam
!                if (weight(isam) == 0._TKG) zeroWeightCount = zeroWeightCount + 1
!                do jdim = 1, ndim
!                    do idim = ROW_RANGE
!                        cov(idim, jdim) = cov(idim, jdim) + (sample(idim, isam) - mean(idim)) * (sample(jdim, isam) - mean(jdim)) * weight(isam)**2
!                    end do
!                end do
!            end do
!            normFac = real(nsam - zeroWeightCount, TKG) / (real(nsam - zeroWeightCount - 1, TKG) * weisum**2)
!            do jdim = 1, ndim; do idim = ROW_RANGE; cov(idim, jdim) = cov(idim, jdim) * normFac; end do; end do;
!        else ! dim = 1_IK
!#if         Def_ENABLED
!            if (shifted) then
!                ! probability-weight correction and shifted.
!                weisum = ZEROW
!                cov(FIRST, FIRST) = 0._TKG
!                do isam = 1, nsam
!                    weisum = weisum + weight(isam)
!                    cov(FIRST, FIRST) = cov(FIRST, FIRST) + (sample(isam, FIRST) * weight(isam))**2
!                    if (weight(isam) == 0._TKG) zeroWeightCount = zeroWeightCount + 1
!                end do
!                normFac = real(nsam - zeroWeightCount, TKG) / (real(nsam - zeroWeightCount - 1, TKG) * weisum**2)
!                cov(FIRST, FIRST) = cov(FIRST, FIRST) * normFac
!                do jdim = COL_RANGE
!                    do idim = ROW_RANGE
!                        cov(idim, jdim) = 0._TKG
!                        do isam = 1, nsam
!                            cov(idim, jdim) = cov(idim, jdim) + sample(isam, idim) * sample(isam, jdim) * weight(isam)**2
!                        end do
!                        cov(idim, jdim) = cov(idim, jdim) * normFac
!                    end do
!                end do
!                return
!            end if
!            allocate(mean(ndim))
!            call setMean(mean, sample, dim, weight, weisum)
!#endif
!            ! probability-weight correction and not shifted.
!            cov(FIRST, FIRST) = 0._TKG
!            do isam = 1, nsam
!                if (weight(isam) == 0._TKG) zeroWeightCount = zeroWeightCount + 1
!                cov(FIRST, FIRST) = cov(FIRST, FIRST) + ((sample(isam, FIRST) - mean(FIRST)) * weight(isam))**2
!            end do
!            normFac = real(nsam - zeroWeightCount, TKG) / (real(nsam - zeroWeightCount - 1, TKG) * weisum**2)
!            cov(FIRST, FIRST) = cov(FIRST, FIRST) * normFac
!            do jdim = COL_RANGE
!                do idim = ROW_RANGE
!                    cov(idim, jdim) = 0._TKG
!                    do isam = 1, nsam
!                        cov(idim, jdim) = cov(idim, jdim) + (sample(isam, idim) - mean(idim) * (sample(isam, jdim) - mean(jdim) * weight(isam)**2
!                    end do
!                    cov(idim, jdim) = cov(idim, jdim) * normFac
!                end do
!            end do
!        end if ! dim
!#if     Def_ENABLED
!        deallocate(mean) ! Bypass gfortran bug for automatic deallocation of local heap arrays.
!#endif
!undef  COL_RANGE

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  TYPE_OF_WEIGHT
#undef  TYPE_OF_SAMPLE
#undef  CHECK_SHAPE_COV
#undef  CHECK_LEN_MEANG
#undef  CHECK_VAL_MEANG
#undef  CHECK_LEN_MEAN
#undef  CHECK_VAL_NSAM
#undef  CHECK_VAL_NDIM
#undef  CHECK_LEN_WEI
#undef  CHECK_VAL_DIM
#undef  CHECK_SUM_WEI
#undef  CHECK_VAL_WEI
#undef  CHECK_WEISUM
#undef  GET_WEIGHTED
#undef  GET_SHIFTED
#undef  GET_CONJG
#undef  SET_CONJG
#undef  OFF_RANGE
#undef  COL_RANGE
#undef  ROW_RANGE
#undef  GET_ABSQ
#undef  GET_PROD
#undef  GET_RE
#undef  FIRST
#undef  ZEROW