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

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getUnifEllsLogPDF_ENABLED && NELL_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: iell
        CHECK_ASSERTION(__LINE__, 0 < nell, SK_"@getUnifEllsLogPDF(): The condition `0 < nell` must hold. ndim, nell = "//getStr([ndim, nell]))
        CHECK_ASSERTION(__LINE__, all(invmul <= 1), SK_"@getUnifEllsLogPDF(): The condition `all(invmul <= 1)` must hold. pack(invmul, mask = invmul > 1) = "//getStr(pack(invmul, mask = invmul > 1)))
        logPDF = -nell * getLogVolUnitBall(real(ndim, TKG)) - log(sum(invmul) / size(invmul, 1, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getUnifEllsLogPDF_ENABLED && CHOL_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: iell
        real(TKG) :: logMulTrace(size(chol, 3, IK)), logMulTraceMax
        CHECK_ASSERTION(__LINE__, size(chol, 1, IK) == size(chol, 2, IK), SK_"@getUnifEllsLogPDF(): The condition `size(chol, 1) == size(chol, 2)` must hold. shape(chol) = "//getStr(shape(chol, IK)))
        CHECK_ASSERTION(__LINE__, all(invmul <= 1), SK_"@getUnifEllsLogPDF(): The condition `all(invmul <= 1)` must hold. pack(invmul, mask = invmul > 1) = "//getStr(pack(invmul, mask = invmul > 1)))
        logMulTraceMax = -huge(logMulTraceMax)
        do  iell = 1, size(chol, 3, IK)
            logMulTrace(iell) = getMatMulTraceLog(chol(:, :, iell))
            logMulTraceMax = max(logMulTraceMax, logMulTrace(iell))
        end do
        logPDF = -getLogSumExp(logMulTrace, logMulTraceMax) * getLogVolUnitBall(real(size(chol, 1, IK), TKG)) - log(sum(invmul) / size(invmul, 1, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   D2_ENABLED && setMMUR_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKG) :: unifrnd
        integer(IK) :: selection, iell, nell, isam, nsam, ndim
#if     AC_ENABLED
        real(TKG) :: maxLogVol, cumPropVol(0 : size(mean, 2, IK))
        CHECK_ASSERTION(__LINE__, all(shape(chol, IK) == shape(invGram, IK)), SK_"@setUnifEllsRand(): The condition `all(shape(chol) == shape(invGram))` must hold. shape(chol), shape(invGram) = "//getStr([shape(chol, IK), shape(invGram, IK)]))
        CHECK_ASSERTION(__LINE__, all(shape(chol, IK) == [size(mean, 1, IK), size(mean, 1, IK), size(mean, 2, IK)]), SK_"@setUnifEllsRand(): The condition `all(shape(chol) == [size(mean, 1), size(mean, 1), size(mean, 2)])` must hold. shape(chol), shape(mean) = "//getStr([shape(chol, IK), shape(mean, IK)]))
        maxLogVol = -huge(maxLogVol)
        cumPropVol(0) = 0._TKG
        do  iell = 1, size(mean, 2, IK)
            cumPropVol(iell) = getMatMulTraceLog(chol(:, :, iell))
            if (maxLogVol < cumPropVol(iell)) maxLogVol = cumPropVol(iell)
        end do
        call setCumPropExp(cumPropVol(1 : size(mean, 2, IK)), maxArray = maxLogVol, control = sequence)
        !#if RNGF_ENABLED
        !block
        !use pm_io, only: disp
        !call disp%show("chol")
        !call disp%show( chol )
        !call disp%show("cumPropVol")
        !call disp%show( cumPropVol )
        !!error stop
        !end block
        !#endif
#elif   !DC_ENABLED
#error  "Unrecognized interface."
#endif
        ndim = size(rand, 1, IK)
        nsam = size(rand, 2, IK)
        nell = size(mean, 2, IK)
        CHECK_ASSERTION(__LINE__, 0 < product(shape(mean, IK)), SK_"@setUnifEllsRand(): The condition `0 < product(shape(mean))` must hold. shape(mean) = "//getStr(shape(mean, IK)))
        CHECK_ASSERTION(__LINE__, size(rand, 1, IK) == size(mean, 1, IK), SK_"@setUnifEllsRand(): The condition `size(rand, 1) == size(mean, 1)` must hold. size(rand, 1), size(mean, 1) = "//getStr([size(rand, 1, IK), size(mean, 1, IK)]))
        CHECK_ASSERTION(__LINE__, all(shape(mahalSq) == [size(mean, 2), size(rand, 2)]), SK_"@setUnifEllsRand(): The condition `all(shape(mahalSq) == [size(mean, 2), size(rand, 2)])` must hold. shape(mahalSq), shape(mean), shape(rand) = "//getStr([shape(mahalSq, IK), shape(mean, IK), shape(rand, IK)]))
        isam = 1_IK
        loopOverSample: do
            !!!!
            !!!! Choose an ellipsoid to draw random vector from.
            !!!!
            invmul(isam) = 0._TKG
#if         AC_ENABLED
            call setUnifRand(rng, unifrnd)
            do  selection = 1, nell
                if (unifrnd < cumPropVol(selection)) exit
            end do
            call setUnifEllRand(rng, rand(:, isam), mean(:, selection), chol(:, :, selection), subset)
            do  iell = 1, nell
                call setDisMahalSq(mahalSq(iell, isam), rand(1 : ndim, isam), invCov = invGram(1 : ndim, 1 : ndim, iell), center = mean(1 : ndim, iell))
                if (mahalSq(iell, isam) <= 1._TKG) invmul(isam) = invmul(isam) + 1._TKG
            end do
#elif       DC_ENABLED
            call setUnifRand(rng, selection, 1_IK, nell)
            call setUnifEllRand(rng, rand(1 : ndim, isam), mean(1 : ndim, selection))
            do  iell = 1, nell
                mahalSq(iell, isam) = sum((rand(1 : ndim, isam) - mean(1 : ndim, iell))**2)
                if (mahalSq(iell, isam) <= 1._TKG) invmul(isam) = invmul(isam) + 1._TKG
            end do
#endif
            if (1._TKG < invmul(isam)) then
                !!!!
                !!!! Accept the current proposal only if it is probabilisitically feasiable.
                !!!!
                invmul(isam) = 1._TKG / invmul(isam)
                call setUnifRand(rng, unifrnd)
                if (invmul(isam) < unifrnd) cycle loopOverSample
            end if
            isam = isam + 1_IK
            if (nsam < isam) exit loopOverSample
        end do loopOverSample
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif