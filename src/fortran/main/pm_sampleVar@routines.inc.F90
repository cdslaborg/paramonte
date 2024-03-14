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
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Saturday 4:40 PM, August 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Set the conjugation rule.
#if     CK_ENABLED
#define TYPE_OF_SAMPLE complex(TKC)
#define GET_ABSQ(X)(real(X)**2 + aimag(X)**2)
#define GET_PROD(X,Y)(X%re * Y%re + X%im * Y%im)
#elif   RK_ENABLED
#define TYPE_OF_SAMPLE real(TKC)
#define GET_PROD(X,Y)X * Y
#define GET_ABSQ(X)X**2
#else
#error  "Unrecognized interface."
#endif
        ! Set the shifting rule.
#if     setVar_ENABLED && Org_ENABLED
#define GET_SHIFTED(X,Y)X
#elif   setVar_ENABLED && Avg_ENABLED
#define GET_SHIFTED(X,Y)(X - Y)
#elif   setVar_ENABLED
#error  "Unrecognized interface."
#endif
        ! Set the weighting rule.
#if     setVar_ENABLED && WNO_ENABLED
#define GET_WEIGHTED(X,W)X
#elif   setVar_ENABLED && (WTI_ENABLED || WTR_ENABLED)
#define GET_WEIGHTED(X,W)X * W
#elif   setVar_ENABLED
#error  "Unrecognized interface."
#endif
        ! Define weight type and kind.
#if     WTI_ENABLED
#define TYPE_OF_WEIGHT integer(IK)
#elif   WTR_ENABLED
#define TYPE_OF_WEIGHT real(TKC)
#elif   !WNO_ENABLED && (getVar_ENABLED || setVarMean_ENABLED)
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
CHECK_ASSERTION(__LINE__, 0 < sum(weight), \
PROC//SK_": The condition `0 < sum(weight)` must hold. weight = "//getStr(weight))
#define CHECK_VAL_WEI(PROC) \
CHECK_ASSERTION(__LINE__, all(0._TKC <= weight), \
PROC//SK_": The condition `all(0. <= weight)` must hold. weight = "//getStr(weight))
#define CHECK_LEN_VAR(PROC) \
CHECK_ASSERTION(__LINE__, size(sample, 3 - dim, IK) == size(var, 1, IK), \
PROC//SK_": The condition `size(sample, 3 - dim) == size(var)` must hold. dim, size(sample, 3 - dim), size(var, 1) = "\
//getStr([dim, size(sample, 3 - dim, IK), size(var, 1, IK)]))
#define CHECK_VAL_MEANG(PROC,DIM) \
CHECK_ASSERTION(__LINE__, all([minval(sample, DIM) <= meang .and. meang <= maxval(sample, DIM)]), \
PROC//SK_": The condition `all([minval(sample, dim) <= meang .and. meang <= maxval(sample, dim)])` must hold. dim, minval(sample, dim), meang, maxval(sample, dim) = "//\
getStr(DIM)//SK_"; "//getStr(minval(sample, DIM))//SK_"; "//getStr(meang)//SK_"; "//getStr(maxval(sample, DIM)))
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
! For very large sample, the following test is very crude and sensitive to the default vs. high-precision summation methods.
#define CHECK_WEISUM(PROC) \
CHECK_ASSERTION(__LINE__, abs(weisum - sum(weight)) < 1000 * epsilon(0._TKC), \
PROC//SK_": The condition `0 < sum(weight)` must hold. weisum, sum(weight) = "//getStr([weisum, sum(weight)]))

        !%%%%%%%%%%%%%%%%%%%%%%%
#if     getVarCorrection_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

#if     Freq_ENABLED
        CHECK_ASSERTION(__LINE__, 1._TKC < weisum, SK_"@setVarCorrection(): The condition `0 < weisum` must hold. weisum = "//getStr(weisum))
        correction = weisum / (weisum - 1._TKC)
#elif   Reli_ENABLED
        CHECK_ASSERTION(__LINE__, 0._TKC < weisum, SK_"@setVarCorrection(): The condition `0 < weisum` must hold. weisum = "//getStr(weisum))
        CHECK_ASSERTION(__LINE__, 0._TKC < weisqs, SK_"@setVarCorrection(): The condition `0 < weisqs` must hold. weisqs = "//getStr(weisqs))
        correction = weisum**2 / (weisum**2 - weisqs)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getVarMerged_ENABLED && New_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setVarMerged(varMerged, varB, varA, meanDiff, fracA)

        !%%%%%%%%%%%%%%%%%%%
#elif   setVarMerged_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        ! Define the output value for setVarMerged.
#if     Old_ENABLED
#define CHECK_LEN_TARGET_VAR
#define TARGET_VAR varB
#elif   New_ENABLED
#define TARGET_VAR varMerged
#define CHECK_LEN_TARGET_VAR \
CHECK_ASSERTION(__LINE__, size(varMerged, 1, IK) == size(varB, 1, IK), SK_"@setVarMerged(): The condition `size(varMerged) == size(varB)` must hold. size(varMerged), size(varB) = "//getStr([size(varMerged, 1, IK), size(varB, 1, IK)]))
#else
#error  "Unrecognized interface."
#endif
#if     D1_ENABLED
        integer(IK) :: idim
        real(TKC) :: fracAB
#endif
        real(TKC) :: fracB
        fracB = 1._TKC - fracA
        CHECK_ASSERTION(__LINE__, 0._TKC < fracA .and. fracA < 1._TKC, SK_"@setVarMerged(): The condition `0 < fracA .and. fracA < 1` must hold. fracA = "//getStr(fracA))
#if     D0_ENABLED
        TARGET_VAR = fracA * varA + fracB * varB + fracA * fracB * GET_ABSQ(meanDiff)
#elif   D1_ENABLED
        CHECK_LEN_TARGET_VAR
        CHECK_ASSERTION(__LINE__, size(varB, 1, IK) == size(varA, 1, IK), SK_"@setVarMerged(): The condition `size(varB) == size(varA)` must hold. size(varB), size(varA) = "//getStr([size(varB, 1, IK), size(varA, 1, IK)]))
        CHECK_ASSERTION(__LINE__, size(meanDiff, 1, IK) == size(varA, 1, IK), SK_"@setVarMerged(): The condition `size(meanDiff) == size(varA)` must hold. size(meanDiff), size(varA) = "//getStr([size(meanDiff, 1, IK), size(varA, 1, IK)]))
        fracAB = fracA * fracB
        do concurrent(idim = 1 : size(varA, 1, IK))
            TARGET_VAR(idim) = fracB * varB(idim) + fracA * varA(idim) + fracAB * GET_ABSQ(meanDiff(idim))
        end do
#else
#error  "Unrecognized interface."
#endif
#undef  CHECK_LEN_TARGET_VAR
#undef  TARGET_VAR

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   setVarMeanMerged_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the output value for setVarMerged.
#if     Old_ENABLED
#define CHECK_LEN_TARGET_AVG
#define CHECK_LEN_TARGET_VAR
#define TARGET_AVG meanB
#define TARGET_VAR varB
#elif   New_ENABLED
#define CHECK_LEN_TARGET_AVG \
CHECK_ASSERTION(__LINE__, size(meanMerged, 1, IK) == size(varB, 1, IK), SK_"@setVarMeanMerged(): The condition `size(meanMerged) == size(varB)` must hold. size(meanMerged), size(varB) = "//getStr([size(meanMerged, 1, IK), size(varB, 1, IK)]))
#define CHECK_LEN_TARGET_VAR \
CHECK_ASSERTION(__LINE__, size(varMerged, 1, IK) == size(varB, 1, IK), SK_"@setVarMeanMerged(): The condition `size(varMerged) == size(varB)` must hold. size(varMerged), size(varB) = "//getStr([size(varMerged, 1, IK), size(varB, 1, IK)]))
#define TARGET_AVG meanMerged
#define TARGET_VAR varMerged
#else
#error  "Unrecognized interface."
#endif
        TYPE_OF_SAMPLE :: temp
#if     D1_ENABLED
        integer(IK) :: idim
        real(TKC) :: fracAB
#endif
        real(TKC) :: fracB
        fracB = 1._TKC - fracA
        CHECK_ASSERTION(__LINE__, 0._TKC < fracA .and. fracA < 1._TKC, SK_"@setVarMeanMerged(): The condition `0 < fracA .and. fracA < 1` must hold. fracA = "//getStr(fracA))
#if     D0_ENABLED
        temp = meanA - meanB
        TARGET_VAR = fracA * varA + fracB * varB + fracA * fracB * GET_ABSQ(temp)
#elif   D1_ENABLED
        CHECK_LEN_TARGET_VAR
        CHECK_ASSERTION(__LINE__, size(varB, 1, IK) == size(varA, 1, IK), SK_"@setVarMeanMerged(): The condition `size(varB) == size(varA)` must hold. size(varB), size(varA) = "//getStr([size(varB, 1, IK), size(varA, 1, IK)]))
        CHECK_ASSERTION(__LINE__, size(meanB, 1, IK) == size(varA, 1, IK), SK_"@setVarMeanMerged(): The condition `size(meanB) == size(varA)` must hold. size(meanB), size(varA) = "//getStr([size(meanB, 1, IK), size(varA, 1, IK)]))
        CHECK_ASSERTION(__LINE__, size(meanA, 1, IK) == size(varA, 1, IK), SK_"@setVarMeanMerged(): The condition `size(meanA) == size(varA)` must hold. size(meanA), size(varA) = "//getStr([size(meanA, 1, IK), size(varA, 1, IK)]))
        fracAB = fracA * fracB
        do concurrent(idim = 1 : size(varA, 1, IK))
            temp = meanA(idim) - meanB(idim)
            TARGET_VAR(idim) = fracB * varB(idim) + fracA * varA(idim) + fracAB * GET_ABSQ(temp)
        end do
#else
#error  "Unrecognized interface."
#endif
        TARGET_AVG = fracA * meanA + fracB * meanB
#undef  CHECK_LEN_TARGET_AVG
#undef  CHECK_LEN_TARGET_VAR
#undef  TARGET_AVG
#undef  TARGET_VAR

        !%%%%%%%%%%%%%
#elif   getVar_ENABLED
        !%%%%%%%%%%%%%

        ! Define the default dimension.
#if     ALL_ENABLED
#define DIM_ARG
#elif   DIM_ENABLED
#define DIM_ARG, dim
#else
#error  "Unrecognized interface."
#endif
        real(TKC) :: normfac
        ! Define the weight sum.
#if     WNO_ENABLED
        real(TKC) :: weisum
        weisum = real(product(shape(sample, IK)), TKC)
        call setVar(var, getMean(sample DIM_ARG), sample DIM_ARG)
#elif   WTI_ENABLED || WTR_ENABLED
        TYPE_OF_WEIGHT :: weisum
        weisum = sum(weight)
        call setVar(var, getMean(sample DIM_ARG, weight), sample DIM_ARG, weight, weisum)
#else
#error  "Unrecognized interface."
#endif
        if (present(correction)) then
            CHECK_ASSERTION(__LINE__, same_type_as(correction, fweight) .or. same_type_as(correction, rweight), SK_"@getVar(): The condition `same_type_as(correction, fweight) .or. same_type_as(correction, rweight)` must hold.")
#if         WNO_ENABLED
            normfac = getVarCorrection(real(weisum, TKC))
#elif       (WTI_ENABLED || WTR_ENABLED)
            if (same_type_as(correction, fweight)) then
                normfac = getVarCorrection(real(weisum, TKC))
            elseif (same_type_as(correction, rweight)) then
                normfac = getVarCorrection(real(weisum, TKC), real(sum(weight**2), TKC))
            end if
#else
#error      "Unrecognized interface."
#endif
            var = var * normfac
        end if
#undef  DIM_ARG

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setVar_ENABLED && D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_OF_SAMPLE :: temp
        integer(IK) :: isam
#if     WNO_ENABLED
        real(TKC) :: weisum
        weisum = real(size(sample, 1, IK), TKC)
#elif   WTI_ENABLED || WTR_ENABLED
        CHECK_WEISUM(SK_"@setVar()")
        CHECK_SUM_WEI(SK_"@setVar()")
        CHECK_VAL_WEI(SK_"@setVar()")
        CHECK_LEN_WEI(SK_"@setVar()",1_IK)
#else
#error  "Unrecognized interface."
#endif
        CHECK_VAL_NSAM(SK_"@setVar()",1_IK)
        var = 0._TKC
        do isam = 1, size(sample, 1, IK)
            temp = GET_SHIFTED(sample(isam),mean)
            var = var + GET_WEIGHTED(GET_ABSQ(temp),weight(isam))
        end do
        var = var / weisum

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setVar_ENABLED && D2_ENABLED && ALL_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_OF_SAMPLE :: temp
        integer(IK) :: idim, jdim
        integer(IK) :: mdim, ndim
#if     WNO_ENABLED
        real(TKC) :: weisum
        weisum = size(sample, kind = IK)
#elif   WTI_ENABLED || WTR_ENABLED
        integer(IK) :: iwei
        iwei = 0_IK
        CHECK_ASSERTION(__LINE__, size(sample, kind = IK) == size(weight, 1, IK), SK_"@setVar(): The condition `size(sample) == size(weight)` must hold. shape(sample), size(weight) = "//getStr([shape(sample, IK), size(weight, 1, IK)]))
        CHECK_SUM_WEI(SK_"@setVar()")
        CHECK_VAL_WEI(SK_"@setVar()")
        CHECK_WEISUM(SK_"@setVar()")
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, 1_IK < size(sample, kind = IK), SK_"@setVar(): The condition `1 < size(sample)` must hold. shape(sample) = "//getStr(shape(sample, IK)))
        mdim = size(sample, 1, IK)
        ndim = size(sample, 2, IK)
        var = 0._TKC
        do jdim = 1, ndim
            do idim = 1, mdim
#if             WTI_ENABLED || WTR_ENABLED
                iwei = iwei + 1_IK
#endif
                temp = GET_SHIFTED(sample(idim, jdim),mean)
                var = var + GET_WEIGHTED(GET_ABSQ(temp),weight(iwei))
            end do
        end do
        var = var / weisum

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setVar_ENABLED && D2_ENABLED && DIM_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKC) :: normFac
        TYPE_OF_SAMPLE :: temp
        integer(IK) :: idim, isam, ndim, nsam
        ndim = size(sample, 3 - dim, IK)
        nsam = size(sample, dim, IK)
#if     WNO_ENABLED
        normFac = 1._TKC / real(nsam, TKC)
#elif   WTI_ENABLED || WTR_ENABLED
        CHECK_WEISUM(SK_"@setVar()")
        CHECK_SUM_WEI(SK_"@setVar()")
        CHECK_VAL_WEI(SK_"@setVar()")
        CHECK_LEN_WEI(SK_"@setVar()",dim)
        normFac = 1._TKC / real(weisum, TKC)
#else
#error  "Unrecognized interface."
#endif
        CHECK_VAL_DIM(SK_"@setVar()")
        CHECK_LEN_VAR(SK_"@setVar()")
        CHECK_VAL_NSAM(SK_"@setVar()",dim)
        CHECK_VAL_NDIM(SK_"@setVar()",3 - dim)
#if     Avg_ENABLED
        CHECK_LEN_MEAN(SK_"@setVar()")
#endif
        if (dim == 2_IK) then
            var = 0._TKC
            do isam = 1, nsam
                do concurrent(idim = 1 : ndim)
                    temp = GET_SHIFTED(sample(idim, isam),mean(idim))
                    var(idim) = var(idim) + GET_WEIGHTED(GET_ABSQ(temp),weight(isam))
                end do
            end do
            var = var * normFac
        else ! dim = 1
            do idim = 1, ndim
                var(idim) = 0._TKC
                do isam = 1, nsam
                    temp = GET_SHIFTED(sample(isam, idim),mean(idim))
                    var(idim) = var(idim) + GET_WEIGHTED(GET_ABSQ(temp),weight(isam))
                end do
                var(idim) = var(idim) * normFac
            end do
        end if ! dim

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setVarMean_ENABLED && D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: isam, nsam
        TYPE_OF_SAMPLE :: temp1
#if     WNO_ENABLED
        real(TKC) :: weisum
        weisum = size(sample, 1, IK)
#elif   WTI_ENABLED || WTR_ENABLED
        TYPE_OF_SAMPLE :: temp2
        CHECK_LEN_WEI(SK_"@setVarMean()",1_IK)
        CHECK_SUM_WEI(SK_"@setVarMean()")
        CHECK_VAL_WEI(SK_"@setVarMean()")
        weisum = 0
#else
#error  "Unrecognized interface."
#endif
        CHECK_VAL_MEANG(SK_"@setVarMean()",1_IK)
        CHECK_VAL_NSAM(SK_"@setVarMean()",1_IK)
        nsam = size(sample, 1, IK)
        mean = 0._TKC
        var = 0._TKC
        do isam = 1, nsam
            temp1 = sample(isam) - meang
#if         WNO_ENABLED
            var = var + GET_ABSQ(temp1)
            mean = mean + temp1
#elif       WTI_ENABLED || WTR_ENABLED
            temp2 = temp1 * weight(isam)
            weisum = weisum + weight(isam)
            var = var + GET_PROD(temp1,temp2)
            mean = mean + temp2
#endif
        end do
        mean = mean / weisum
        var = (var - GET_ABSQ(mean) * weisum) / weisum
        mean = mean + meang

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setVarMean_ENABLED && D2_ENABLED && ALL_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim, jdim
        integer(IK) :: mdim, ndim
        TYPE_OF_SAMPLE :: temp1
#if     WNO_ENABLED
        real(TKC) :: weisum
        weisum = size(sample, kind = IK)
#elif   WTI_ENABLED || WTR_ENABLED
        integer(IK) :: iwei
        TYPE_OF_SAMPLE :: temp2
        CHECK_ASSERTION(__LINE__, size(sample, kind = IK) == size(weight, 1, IK), SK_"@setVarMean(): The condition `size(sample) == size(weight)` must hold. shape(sample), size(weight) = "//getStr([shape(sample, IK), size(weight, 1, IK)]))
        CHECK_SUM_WEI(SK_"@setVarMean()")
        CHECK_VAL_WEI(SK_"@setVarMean()")
        iwei = 0_IK
        weisum = 0
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, 1_IK < size(sample, kind = IK), SK_"@setVarMean(): The condition `1 < size(sample)` must hold. shape(sample) = "//getStr(shape(sample, IK)))
        CHECK_ASSERTION(__LINE__, minval(sample) <= meang .and. meang <= maxval(sample), SK_"@setVarMean(): The condition `minval(sample) <= meang .and. meang <= minval(sample)` must hold. minval(sample), meang, maxval(sample) = "//getStr([minval(sample), meang, maxval(sample)]))
        mdim = size(sample, 1, IK)
        ndim = size(sample, 2, IK)
        mean = 0._TKC
        var = 0._TKC
        do jdim = 1, ndim
            do idim = 1, mdim
                temp1 = sample(idim, jdim) - meang
#if             WNO_ENABLED
                var = var + GET_ABSQ(temp1)
                mean = mean + temp1
#elif           WTI_ENABLED || WTR_ENABLED
                iwei = iwei + 1_IK
                temp2 = temp1 * weight(iwei)
                weisum = weisum + weight(iwei)
                var = var + GET_PROD(temp1,temp2)
                mean = mean + temp2
#endif
            end do
        end do
        mean = mean / weisum
        var = (var - GET_ABSQ(mean) * weisum) / weisum
        mean = mean + meang

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setVarMean_ENABLED && D2_ENABLED && DIM_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKC) :: normFac
        TYPE_OF_SAMPLE :: temp1
#if     WTI_ENABLED || WTR_ENABLED
        TYPE_OF_SAMPLE :: temp2
#endif
        integer(IK) :: idim, isam, ndim, nsam
        ndim = size(sample, 3 - dim, IK)
        nsam = size(sample, dim, IK)
        CHECK_VAL_DIM(SK_"@setVarMean()")
        CHECK_LEN_VAR(SK_"@setVarMean()")
        CHECK_VAL_NSAM(SK_"@setVarMean()",dim)
        CHECK_VAL_NDIM(SK_"@setVarMean()",3 - dim)
        CHECK_VAL_MEANG(SK_"@setVarMean()",dim)
        CHECK_LEN_MEAN(SK_"@setVarMean()")
#if     WNO_ENABLED
        normFac = 1._TKC / real(nsam, TKC)
        if (dim == 2_IK) then
            var = 0._TKC
            mean = 0._TKC
            do isam = 1, nsam
                do idim = 1, ndim
                    temp1 = sample(idim, isam) - meang(idim)
                    var(idim) = var(idim) + GET_ABSQ(temp1)
                    mean(idim) = mean(idim) + temp1
                end do
            end do
            do concurrent(idim = 1 : ndim)
                mean(idim) = mean(idim) * normFac
                var(idim) = (var(idim) - GET_ABSQ(mean(idim)) * nsam) * normFac
                mean(idim) = mean(idim) + meang(idim)
            end do
        else
            do concurrent(idim = 1 : ndim)
                var(idim) = 0._TKC
                mean(idim) = 0._TKC
                do isam = 1, nsam
                    temp1 = sample(isam, idim) - meang(idim)
                    var(idim) = var(idim) + GET_ABSQ(temp1)
                    mean(idim) = mean(idim) + temp1
                end do
                mean(idim) = mean(idim) * normFac
                var(idim) = (var(idim) - GET_ABSQ(mean(idim)) * nsam) * normFac
                mean(idim) = mean(idim) + meang(idim)
            end do
        end if
#elif   WTI_ENABLED || WTR_ENABLED
        CHECK_LEN_WEI(SK_"@setVarMean()",dim)
        CHECK_SUM_WEI(SK_"@setVarMean()")
        CHECK_VAL_WEI(SK_"@setVarMean()")
        weisum = 0
        if (dim == 2_IK) then
            var = 0._TKC
            mean = 0._TKC
            do isam = 1, nsam
                weisum = weisum + weight(isam)
                do concurrent(idim = 1 : ndim)
                    temp1 = sample(idim, isam) - meang(idim)
                    temp2 = temp1 * weight(isam)
                    var(idim) = var(idim) + GET_PROD(temp1,temp2)
                    mean(idim) = mean(idim) + temp2
                end do
            end do
            normFac = 1._TKC / weisum
            do concurrent(idim = 1 : ndim)
                mean(idim) = mean(idim) * normFac
                var(idim) = (var(idim) - GET_ABSQ(mean(idim)) * weisum) * normFac
                mean(idim) = mean(idim) + meang(idim)
            end do
        else
            idim = 1_IK
            var(idim) = 0._TKC
            mean(idim) = 0._TKC
            do isam = 1, nsam
                weisum = weisum + weight(isam)
                temp1 = sample(isam, idim) - meang(idim)
                temp2 = temp1 * weight(isam)
                var(idim) = var(idim) + GET_PROD(temp1,temp2)
                mean(idim) = mean(idim) + temp2
            end do
            normFac = 1._TKC / real(weisum, TKC)
            mean(idim) = mean(idim) * normFac
            var(idim) = (var(idim) - GET_ABSQ(mean(idim)) * weisum) * normFac
            mean(idim) = mean(idim) + meang(idim)
            do idim = 2, ndim
                var(idim) = 0._TKC
                mean(idim) = 0._TKC
                do isam = 1, nsam
                    temp1 = sample(isam, idim) - meang(idim)
                    temp2 = temp1 * weight(isam)
                    var(idim) = var(idim) + GET_PROD(temp1,temp2)
                    mean(idim) = mean(idim) + temp2
                end do
                mean(idim) = mean(idim) * normFac
                var(idim) = (var(idim) - GET_ABSQ(mean(idim)) * weisum) * normFac
                mean(idim) = mean(idim) + meang(idim)
            end do
        end if
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  CHECK_VAL_MEANG
#undef  CHECK_LEN_MEANG
#undef  CHECK_LEN_MEAN
#undef  TYPE_OF_SAMPLE
#undef  TYPE_OF_WEIGHT
#undef  CHECK_VAL_NSAM
#undef  CHECK_VAL_NDIM
#undef  CHECK_LEN_WEI
#undef  CHECK_VAL_DIM
#undef  CHECK_LEN_VAR
#undef  CHECK_SUM_WEI
#undef  CHECK_VAL_WEI
#undef  CHECK_WEISUM
#undef  GET_WEIGHTED
#undef  GET_SHIFTED
#undef  GET_ABSQ
#undef  GET_PROD
#undef  FIRST
