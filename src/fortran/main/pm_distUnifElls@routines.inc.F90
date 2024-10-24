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
!>  This include file contains the implementation of procedures in [pm_distUnifElls](@ref pm_distUnifElls).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%
#if     getUnifEllsLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: nell, isim, def_nsim, membership
        real(TKG) :: rand(size(mean, 1, IK)), mahalSq(size(mean, 2, IK)), invmul, invmulsum
#if     AC_ENABLED
        integer(IK) :: iell
        real(TKG) :: maxLogVol, cumPropVol(size(mean, 2, IK))
        CHECK_ASSERTION(__LINE__, size(chol, 1, IK) <= size(chol, 2, IK), SK_"@getUnifEllsLogPDF(): The condition `size(chol, 1) <= size(chol, 2)` must hold. shape(chol) = "//getStr(shape(chol, IK)))
        CHECK_ASSERTION(__LINE__, size(chol, 1, IK) == size(mean, 1, IK) .and. size(chol, 3, IK) == size(mean, 2, IK), SK_"@getUnifEllsLogPDF(): The condition `size(chol, 1) == size(mean, 1) .and. size(chol, 3) == size(mean, 2)` must hold. shape(chol), shape(mean) = "//getStr([shape(chol, IK), shape(mean, IK)]))
        CHECK_ASSERTION(__LINE__, all(shape(invGram, IK) == [size(mean, 1, IK), size(mean, 1, IK), size(mean, 2, IK)]), SK_"@getUnifEllsLogPDF(): The condition `all(shape(invGram) == [size(mean, 1), size(mean, 1), size(mean, 2)])` must hold. shape(invGram), shape(mean) = "//getStr([shape(invGram, IK), shape(mean, IK)]))
        maxLogVol = -huge(maxLogVol)
        do  iell = 1, size(mean, 2, IK)
            cumPropVol(iell) = getMatMulTraceLog(chol(:, 1 : size(mean, 1, IK), iell))
            if (maxLogVol < cumPropVol(iell)) maxLogVol = cumPropVol(iell)
        end do
        call setCumPropExp(cumPropVol, maxArray = maxLogVol, control = sequence)
        logPDF = log(sum(exp(cumPropVol - maxLogVol))) + maxLogVol
#elif   DC_ENABLED
        logPDF = log(real(size(mean, 2, IK), TKG))
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, 0 < product(shape(mean, IK)), SK_"@getUnifEllsLogPDF(): The condition `0 < product(shape(mean))` must hold. shape(mean) = "//getStr(shape(mean, IK)))
        nell = size(mean, 2, IK)
        if (present(nsim)) then
            CHECK_ASSERTION(__LINE__, 0 < nsim, SK_"@getUnifEllsLogPDF(): The condition `0 < nsim` must hold. nsim = "//getStr(nsim))
            def_nsim = nsim
        else
            def_nsim = 10000_IK
        end if
        invmulsum = 0._TKG
        do  isim = 1, def_nsim
#if         AC_ENABLED
            call setUnifEllsRand(rng, rand, mahalSq, invmul, membership, mean, chol, subset, invGram, cumPropVol)
#elif       DC_ENABLED
            call setUnifEllsRand(rng, rand, mahalSq, invmul, membership, mean)
#else
#error      "Unrecognized interface."
#endif
            invmulsum = invmulsum + invmul
        end do
        logPDF = logPDF + log(invmulsum / def_nsim)
        if (present(normed)) then
            if (normed) logPDF = logPDF + getLogVolUnitBall(real(size(rand, 1, IK), TKG))
        end if

!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!#elif   getUnifEllsLogPDF_ENABLED && NELL_ENABLED
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!        integer(IK) :: iell
!        check_assertion(__LINE__, 0 < nell, SK_"@getUnifEllsLogPDF(): The condition `0 < nell` must hold. ndim, nell = "//getStr([ndim, nell]))
!        check_assertion(__LINE__, all(invmul <= 1), SK_"@getUnifEllsLogPDF(): The condition `all(invmul <= 1)` must hold. pack(invmul, mask = invmul > 1) = "//getStr(pack(invmul, mask = invmul > 1)))
!        logPDF = -nell * getLogVolUnitBall(real(ndim, TKG)) - log(sum(invmul) / size(invmul, 1, IK))
!
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!#elif   getUnifEllsLogPDF_ENABLED && CHOL_ENABLED
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!        integer(IK) :: iell
!        real(TKG) :: logMulTrace(size(chol, 3, IK)), logMulTraceMax
!        check_assertion(__LINE__, size(chol, 1, IK) <= size(chol, 2, IK), SK_"@getUnifEllsLogPDF(): The condition `size(chol, 1) == size(chol, 2)` must hold. shape(chol) = "//getStr(shape(chol, IK)))
!        check_assertion(__LINE__, all(invmul <= 1), SK_"@getUnifEllsLogPDF(): The condition `all(invmul <= 1)` must hold. pack(invmul, mask = invmul > 1) = "//getStr(pack(invmul, mask = invmul > 1)))
!        logMulTraceMax = -huge(logMulTraceMax)
!        do  iell = 1, size(chol, 3, IK)
!            logMulTrace(iell) = getMatMulTraceLog(chol(:, :, iell))
!            logMulTraceMax = max(logMulTraceMax, logMulTrace(iell))
!        end do
!        logPDF = -getLogSumExp(logMulTrace, logMulTraceMax) * getLogVolUnitBall(real(size(chol, 1, IK), TKG)) - log(sum(invmul) / size(invmul, 1, IK))
!
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMMUR_ENABLED && D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKG) :: unifrnd
        integer(IK) :: iell, nell, ndim
#if     AC_ENABLED
        CHECK_ASSERTION(__LINE__, size(chol, 1, IK) <= size(chol, 2, IK), SK_"@setUnifEllsRand(): The condition `size(chol, 1) <= size(chol, 2)` must hold. shape(chol) = "//getStr(shape(chol, IK)))
        CHECK_ASSERTION(__LINE__, size(chol, 1, IK) == size(mean, 1, IK) .and. size(chol, 3, IK) == size(mean, 2, IK), SK_"@setUnifEllsRand(): The condition `size(chol, 1) == size(mean, 1) .and. size(chol, 3) == size(mean, 2)` must hold. shape(chol), shape(mean) = "//getStr([shape(chol, IK), shape(mean, IK)]))
        CHECK_ASSERTION(__LINE__, all(shape(invGram, IK) == [size(mean, 1, IK), size(mean, 1, IK), size(mean, 2, IK)]), SK_"@setUnifEllsRand(): The condition `all(shape(chol) == [size(mean, 1), size(mean, 1), size(mean, 2)])` must hold. shape(chol), shape(mean) = "//getStr([shape(chol, IK), shape(mean, IK)]))
        CHECK_ASSERTION(__LINE__, size(mean, 2, IK) == size(cumPropVol, 1, IK), SK_"@setUnifEllsRand(): The condition `size(mean, 2) == size(cumPropVol)` must hold. size(mean, 2), size(cumPropVol) = "//getStr([size(mean, 2, IK) == size(cumPropVol, 1, IK)]))
#elif   !DC_ENABLED
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, size(rand, 1, IK) == size(mean, 1, IK), SK_"@setUnifEllsRand(): The condition `size(rand, 1) == size(mean, 1)` must hold. size(rand, 1), size(mean, 1) = "//getStr([size(rand, 1, IK), size(mean, 1, IK)]))
        CHECK_ASSERTION(__LINE__, 0 < product(shape(mean, IK)), SK_"@setUnifEllsRand(): The condition `0 < product(shape(mean))` must hold. shape(mean) = "//getStr(shape(mean, IK)))
        ndim = size(rand, 1, IK)
        nell = size(mean, 2, IK)
        !!!!
        !!!! Choose an ellipsoid to draw random vector from.
        !!!!
        loopOverSample: do
            invmul = 0._TKG
#if         AC_ENABLED
            call setUnifRand(rng, unifrnd)
            do  membership = 1, nell
                if (unifrnd < cumPropVol(membership)) exit
            end do
            call setUnifEllRand(rng, rand(1 : ndim), mean(1 : ndim, membership), chol(:, 1 : ndim, membership), subset)
            do  iell = 1, nell
                call setDisMahalSq(mahalSq(iell), rand(1 : ndim), invCov = invGram(:, 1 : ndim, iell), center = mean(1 : ndim, iell))
                if (mahalSq(iell) <= 1._TKG) invmul = invmul + 1._TKG
            end do
#elif       DC_ENABLED
            call setUnifRand(rng, membership, 1_IK, nell)
            call setUnifEllRand(rng, rand(1 : ndim), mean(1 : ndim, membership))
            do  iell = 1, nell
                mahalSq(iell) = sum((rand(1 : ndim) - mean(1 : ndim, iell))**2)
                if (mahalSq(iell) <= 1._TKG) invmul = invmul + 1._TKG
            end do
#endif
            if (1._TKG < invmul) then
                !!!!
                !!!! Accept the current proposal only if it is probabilistically feasible.
                !!!!
                invmul = 1._TKG / invmul
                call setUnifRand(rng, unifrnd)
                if (invmul < unifrnd) cycle loopOverSample
            end if
            exit loopOverSample
        end do loopOverSample

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMMUR_ENABLED && D2_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: nell, isam, nsam, ndim
#if     AC_ENABLED
        real(TKG) :: maxLogVol, cumPropVol(size(mean, 2, IK))
        integer(IK) :: iell
        ndim = size(rand, 1, IK)
        maxLogVol = -huge(maxLogVol)
        do  iell = 1, size(mean, 2, IK)
            cumPropVol(iell) = getMatMulTraceLog(chol(:, 1 : ndim, iell))
            if (maxLogVol < cumPropVol(iell)) maxLogVol = cumPropVol(iell)
        end do
        call setCumPropExp(cumPropVol, maxArray = maxLogVol, control = sequence)
#elif   DC_ENABLED
        ndim = size(rand, 1, IK)
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, 0 < product(shape(mean, IK)), SK_"@setUnifEllsRand(): The condition `0 < product(shape(mean))` must hold. shape(mean) = "//getStr(shape(mean, IK)))
        CHECK_ASSERTION(__LINE__, size(rand, 1, IK) == size(mean, 1, IK), SK_"@setUnifEllsRand(): The condition `size(rand, 1) == size(mean, 1)` must hold. size(rand, 1), size(mean, 1) = "//getStr([size(rand, 1, IK), size(mean, 1, IK)]))
        CHECK_ASSERTION(__LINE__, all(shape(mahalSq) == [size(mean, 2), size(rand, 2)]), SK_"@setUnifEllsRand(): The condition `all(shape(mahalSq) == [size(mean, 2), size(rand, 2)])` must hold. shape(mahalSq), shape(mean), shape(rand) = "//getStr([shape(mahalSq, IK), shape(mean, IK), shape(rand, IK)]))
        nsam = size(rand, 2, IK)
        nell = size(mean, 2, IK)
        do  isam = 1, nsam
#if         AC_ENABLED
            call setUnifEllsRand(rng, rand(1 : ndim, isam), mahalSq(1 : nell, isam), invmul(isam), membership(isam), mean, chol, subset, invGram, cumPropVol)
#elif       DC_ENABLED
            call setUnifEllsRand(rng, rand(1 : ndim, isam), mahalSq(1 : nell, isam), invmul(isam), membership(isam), mean)
#else
#error      "Unrecognized interface."
#endif
        end do
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif