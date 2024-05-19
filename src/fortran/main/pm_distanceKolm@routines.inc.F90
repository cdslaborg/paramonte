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
!>  This include file contains the implementation of procedures in [pm_distanceKolm](@ref pm_distanceKolm).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     WII_ENABLED || WID_ENABLED
#define TYPE_OF_WEIGHT integer(IK)
#elif   WRR_ENABLED || WRD_ENABLED
#define TYPE_OF_WEIGHT real(TKG)
#elif   !WDD_ENABLED
#error  "Unrecognized interface."
#endif
#if     SXD_ENABLED || SXA_ENABLED
#define CDF_ARG
#elif   SCD_ENABLED || SCA_ENABLED
#define CDF_ARG, getCDF
#elif   !(SSD_ENABLED || SSA_ENABLED)
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getDisKolm_ENABLED && (SXD_ENABLED || SCD_ENABLED) && WDD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKG) :: copySample1(size(sample1, 1, IK))
        copySample1 = sample1
        call setDisKolm(disKolm, copySample1 CDF_ARG)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setDisKolm_ENABLED && (SXD_ENABLED || SCD_ENABLED) && WDD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setSorted(sample1)
        call setDisKolm(disKolm, sample1 CDF_ARG, ascending)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (getDisKolm_ENABLED || setDisKolm_ENABLED) && (SXA_ENABLED || SCA_ENABLED) && WDD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: isam1, nsam1
        real(TKG) :: invWeiSum1, cdfn, cdfp, cdfr ! next, previous, reference
        !check_assertion(__LINE__, 0 < size(sample1, 1, IK), SK_"@setDisKolm(): The condition `0 < size(sample1)` must hold. size(sample1) = "//getStr(size(sample1, 1, IK)))
        CHECK_ASSERTION(__LINE__, isAscending(sample1), SK_"@setDisKolm(): The condition `isAscending(sample1)` must hold. sample1 = "//getStr(sample1))
#if     SXA_ENABLED
        CHECK_ASSERTION(__LINE__, all(0 <= sample1) .and. all(sample1 <= 1), SK_"@setDisKolm(): The condition `all(0 <= sample1) .and. all(sample1 <= 1)` must hold. sample1 = "//getStr(sample1))
#endif
        nsam1 = size(sample1, 1, IK)
        if (0 < nsam1) invWeiSum1 = 1._TKG / nsam1
        disKolm = 0._TKG
        cdfp = 0._TKG
        do isam1 = 1, nsam1
#if         SXA_ENABLED
            cdfr = sample1(isam1)
#elif       SCA_ENABLED
            cdfr = getCDF(sample1(isam1))
            CHECK_ASSERTION(__LINE__, 0 <= cdfr .and. cdfr <= 1, SK_"@setDisKolm(): The condition `0 <= getCDF(sample1(isam1)) <= 1` must hold. sample1(isam1), getCDF(sample1(isam1)) = "//getStr([sample1(isam1), cdfr]))
#else
#error      "Unrecognized interface."
#endif
            cdfn = isam1 * invWeiSum1
            disKolm = max(disKolm, abs(cdfr - cdfp), abs(cdfr - cdfn))
            cdfp = cdfn
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getDisKolm_ENABLED && (SXD_ENABLED || SCD_ENABLED) && (WID_ENABLED | WRD_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_OF_WEIGHT :: copyWeight1(size(weight1, 1, IK))
        real(TKG) :: copySample1(size(sample1, 1, IK))
        copySample1 = sample1
        copyWeight1 = weight1
        call setDisKolm(disKolm, copySample1, copyWeight1, weisum1 CDF_ARG)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setDisKolm_ENABLED && (SXD_ENABLED || SCD_ENABLED) && (WID_ENABLED | WRD_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: index1(size(sample1, 1, IK))
        call setSorted(sample1, index1)
        sample1 = sample1(index1)
        weight1 = weight1(index1)
        call setDisKolm(disKolm, sample1, weight1, weisum1 CDF_ARG, ascending)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (getDisKolm_ENABLED || setDisKolm_ENABLED) && (SXA_ENABLED || SCA_ENABLED) && (WID_ENABLED | WRD_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: isam1, nsam1
        real(TKG) :: invWeiSum1, cdfn, cdfp, cdfr ! next, previous, reference
        CHECK_ASSERTION(__LINE__, all(0 <= weight1), SK_"@setDisKolm(): The condition `all(0 <= weight1)` must hold. weight1 = "//getStr(sample1))
        CHECK_ASSERTION(__LINE__, isAscending(sample1), SK_"@setDisKolm(): The condition `isAscending(sample1)` must hold. sample1 = "//getStr(sample1))
#if     SXA_ENABLED
        CHECK_ASSERTION(__LINE__, all(0 <= sample1) .and. all(sample1 <= 1), SK_"@setDisKolm(): The condition `all(0 <= sample1) .and. all(sample1 <= 1)` must hold. sample1 = "//getStr(sample1))
#endif
        CHECK_ASSERTION(__LINE__, abs(weisum1 - sum(weight1)) < 100 * epsilon(0._TKG), SK_"@setDisKolm(): The condition `weisum1 == sum(weight1)` must hold. weisum1, sum(weight1) = "//getStr([weisum1, sum(weight1)]))
        CHECK_ASSERTION(__LINE__, size(sample1, 1, IK) == size(weight1, 1, IK), SK_"@setDisKolm(): The condition `size(sample1) == size(weight1)` must hold. size(sample1), size(weight1) = "//getStr([size(sample1, 1, IK), size(weight1, 1, IK)]))
        nsam1 = size(sample1, 1, IK)
        if (0 < weisum1) invWeiSum1 = 1._TKG / weisum1
        disKolm = 0._TKG
        cdfn = 0._TKG
        cdfp = 0._TKG
        do isam1 = 1, nsam1
#if         SXA_ENABLED
            cdfr = sample1(isam1)
#elif       SCA_ENABLED
            cdfr = getCDF(sample1(isam1))
            CHECK_ASSERTION(__LINE__, 0 <= cdfr .and. cdfr <= 1, SK_"@setDisKolm(): The condition `0 <= getCDF(sample1(isam1)) <= 1` must hold. sample1(isam1), getCDF(sample1(isam1)) = "//getStr([sample1(isam1), cdfr]))
#else
#error      "Unrecognized interface."
#endif
            cdfn = cdfn + weight1(isam1) * invWeiSum1
            disKolm = max(disKolm, abs(cdfr - cdfp), abs(cdfr - cdfn))
            cdfp = cdfn
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getDisKolm_ENABLED && SSD_ENABLED && WDD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKG) :: copySample1(size(sample1, 1, IK))
        real(TKG) :: copySample2(size(sample2, 1, IK))
        copySample1 = sample1
        copySample2 = sample2
        call setDisKolm(disKolm, copySample1, copySample2)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setDisKolm_ENABLED && SSD_ENABLED && WDD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setSorted(sample1)
        call setSorted(sample2)
        call setDisKolm(disKolm, sample1, sample2, ascending)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (getDisKolm_ENABLED || setDisKolm_ENABLED) && SSA_ENABLED && WDD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKG) :: cdf1, cdf2
        real(TKG) :: ell1, ell2
        integer(IK) :: isam1, isam2
        integer(IK) :: nsam1, nsam2
        real(TKG) :: invWeiSum1, invWeiSum2
        CHECK_ASSERTION(__LINE__, isAscending(sample1), SK_"@setDisKolm(): The condition `isAscending(sample1)` must hold. sample1 = "//getStr(sample1))
        CHECK_ASSERTION(__LINE__, isAscending(sample2), SK_"@setDisKolm(): The condition `isAscending(sample2)` must hold. sample2 = "//getStr(sample2))
        nsam1 = size(sample1, 1, IK)
        nsam2 = size(sample2, 1, IK)
        if (0_IK < nsam1) invWeiSum1 = 1._TKG / nsam1
        if (0_IK < nsam2) invWeiSum2 = 1._TKG / nsam2
        isam1 = 1_IK
        isam2 = 1_IK
        cdf1 = 0._TKG
        cdf2 = 0._TKG
        disKolm = 0._TKG
        do
            if (nsam1 < isam1 .or. nsam2 < isam2) exit
            ell1 = sample1(isam1)
            ell2 = sample2(isam2)
            if (ell1 <= ell2) then
              cdf1 = isam1 * invWeiSum1
              isam1 = isam1 + 1_IK
            end if
            if (ell2 <= ell1) then
              cdf2 = isam2 * invWeiSum2
              isam2 = isam2 + 1_IK
            end if
            disKolm = max(disKolm, abs(cdf2 - cdf1))
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getDisKolm_ENABLED && SSD_ENABLED && (WID_ENABLED | WRD_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_OF_WEIGHT :: copyWeight1(size(weight1, 1, IK))
        real(TKG) :: copySample1(size(sample1, 1, IK))
        real(TKG) :: copySample2(size(sample2, 1, IK))
        copySample1 = sample1
        copySample2 = sample2
        copyWeight1 = weight1
        call setDisKolm(disKolm, copySample1, copyWeight1, weisum1, copySample2)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setDisKolm_ENABLED && SSD_ENABLED && (WID_ENABLED | WRD_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: index1(size(sample1, 1, IK))
        call setSorted(sample1, index1)
        sample1 = sample1(index1)
        weight1 = weight1(index1)
        call setSorted(sample2)
        call setDisKolm(disKolm, sample1, weight1, weisum1, sample2, ascending)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (getDisKolm_ENABLED || setDisKolm_ENABLED) && SSA_ENABLED && (WID_ENABLED | WRD_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKG) :: cdf1, cdf2
        real(TKG) :: ell1, ell2
        integer(IK) :: isam1, isam2
        integer(IK) :: nsam1, nsam2
        real(TKG) :: invWeiSum1, invWeiSum2
        CHECK_ASSERTION(__LINE__, all(0 <= weight1), SK_"@setDisKolm(): The condition `all(0 <= weight1)` must hold. weight1 = "//getStr(sample1))
        CHECK_ASSERTION(__LINE__, isAscending(sample1), SK_"@setDisKolm(): The condition `isAscending(sample1)` must hold. sample1 = "//getStr(sample1))
        CHECK_ASSERTION(__LINE__, isAscending(sample2), SK_"@setDisKolm(): The condition `isAscending(sample2)` must hold. sample2 = "//getStr(sample2))
        CHECK_ASSERTION(__LINE__, abs(weisum1 - sum(weight1)) < 100 * epsilon(0._TKG), SK_"@setDisKolm(): The condition `weisum1 == sum(weight1)` must hold. weisum1, sum(weight1) = "//getStr([weisum1, sum(weight1)]))
        CHECK_ASSERTION(__LINE__, size(sample1, 1, IK) == size(weight1, 1, IK), SK_"@setDisKolm(): The condition `size(sample1) == size(weight1)` must hold. size(sample1), size(weight1) = "//getStr([size(sample1, 1, IK), size(weight1, 1, IK)]))
        nsam1 = size(sample1, 1, IK)
        nsam2 = size(sample2, 1, IK)
        if (0 < weisum1) invWeiSum1 = 1._TKG / weisum1
        if (0 < nsam2) invWeiSum2 = 1._TKG / nsam2
        disKolm = 0._TKG
        cdf1 = 0._TKG
        cdf2 = 0._TKG
        isam1 = 1_IK
        isam2 = 1_IK
        do
            if (nsam1 < isam1 .or. nsam2 < isam2) exit
            ell1 = sample1(isam1)
            ell2 = sample2(isam2)
            if (ell1 <= ell2) then
              cdf1 = cdf1 + weight1(isam1) * invWeiSum1
              isam1 = isam1 + 1_IK
            end if
            if (ell2 <= ell1) then
              cdf2 = isam2 * invWeiSum2
              isam2 = isam2 + 1_IK
            end if
            disKolm = max(disKolm, abs(cdf2 - cdf1))
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getDisKolm_ENABLED && SSD_ENABLED && (WII_ENABLED | WRR_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_OF_WEIGHT :: copyWeight1(size(weight1, 1, IK))
        TYPE_OF_WEIGHT :: copyWeight2(size(weight2, 1, IK))
        real(TKG) :: copySample1(size(sample1, 1, IK))
        real(TKG) :: copySample2(size(sample2, 1, IK))
        copySample1 = sample1
        copySample2 = sample2
        copyWeight1 = weight1
        copyWeight2 = weight2
        call setDisKolm(disKolm, copySample1, copyWeight1, weisum1, copySample2, copyWeight2, weisum2)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setDisKolm_ENABLED && SSD_ENABLED && (WII_ENABLED | WRR_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: index1(size(sample1, 1, IK))
        integer(IK) :: index2(size(sample2, 1, IK))
        call setSorted(sample1, index1)
        call setSorted(sample2, index2)
        sample1 = sample1(index1)
        sample2 = sample2(index2)
        weight1 = weight1(index1)
        weight2 = weight2(index2)
        call setDisKolm(disKolm, sample1, weight1, weisum1, sample2, weight2, weisum2, ascending)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (getDisKolm_ENABLED || setDisKolm_ENABLED) && SSA_ENABLED && (WII_ENABLED | WRR_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKG) :: cdf1, cdf2
        real(TKG) :: ell1, ell2
        integer(IK) :: isam1, isam2
        integer(IK) :: nsam1, nsam2
        real(TKG) :: invWeiSum1, invWeiSum2
        CHECK_ASSERTION(__LINE__, all(0 <= weight1), SK_"@setDisKolm(): The condition `all(0 <= weight1)` must hold. weight1 = "//getStr(sample1))
        CHECK_ASSERTION(__LINE__, all(0 <= weight2), SK_"@setDisKolm(): The condition `all(0 <= weight2)` must hold. weight2 = "//getStr(sample2))
        CHECK_ASSERTION(__LINE__, isAscending(sample1), SK_"@setDisKolm(): The condition `isAscending(sample1)` must hold. sample1 = "//getStr(sample1))
        CHECK_ASSERTION(__LINE__, isAscending(sample2), SK_"@setDisKolm(): The condition `isAscending(sample2)` must hold. sample2 = "//getStr(sample2))
        CHECK_ASSERTION(__LINE__, abs(weisum1 - sum(weight1)) < 100 * epsilon(0._TKG), SK_"@setDisKolm(): The condition `weisum1 == sum(weight1)` must hold. weisum1, sum(weight1) = "//getStr([weisum1, sum(weight1)]))
        CHECK_ASSERTION(__LINE__, abs(weisum2 - sum(weight2)) < 100 * epsilon(0._TKG), SK_"@setDisKolm(): The condition `weisum2 == sum(weight2)` must hold. weisum2, sum(weight2) = "//getStr([weisum2, sum(weight2)]))
        CHECK_ASSERTION(__LINE__, size(sample1, 1, IK) == size(weight1, 1, IK), SK_"@setDisKolm(): The condition `size(sample1) == size(weight1)` must hold. size(sample1), size(weight1) = "//getStr([size(sample1, 1, IK), size(weight1, 1, IK)]))
        CHECK_ASSERTION(__LINE__, size(sample2, 1, IK) == size(weight2, 1, IK), SK_"@setDisKolm(): The condition `size(sample2) == size(weight2)` must hold. size(sample2), size(weight2) = "//getStr([size(sample2, 1, IK), size(weight2, 1, IK)]))
        nsam1 = size(sample1, 1, IK)
        nsam2 = size(sample2, 1, IK)
        if (0 < weisum1) invWeiSum1 = 1._TKG / weisum1
        if (0 < weisum2) invWeiSum2 = 1._TKG / weisum2
        disKolm = 0._TKG
        cdf1 = 0._TKG
        cdf2 = 0._TKG
        isam1 = 1_IK
        isam2 = 1_IK
        do
            if (nsam1 < isam1 .or. nsam2 < isam2) exit
            ell1 = sample1(isam1)
            ell2 = sample2(isam2)
            if (ell1 <= ell2) then
              cdf1 = cdf1 + weight1(isam1) * invWeiSum1
              isam1 = isam1 + 1_IK
            end if
            if (ell2 <= ell1) then
              cdf2 = cdf2 + weight2(isam2) * invWeiSum2
              isam2 = isam2 + 1_IK
            end if
            disKolm = max(disKolm, abs(cdf2 - cdf1))
        end do
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  TYPE_OF_WEIGHT
#undef  CDF_ARG