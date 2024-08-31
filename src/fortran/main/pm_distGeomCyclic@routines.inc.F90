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
!>  This include file contains the implementation of procedures in [pm_distGeomCyclic](@ref pm_distGeomCyclic).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getGeomCyclicLogPMF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0._RKG < probSuccess, SK_"@getGeomCyclicLogPMF(): The condition `0 <= probSuccess` must hold. probSuccess = "//getStr(probSuccess))
        CHECK_ASSERTION(__LINE__, probSuccess <= 1._RKG, SK_"@getGeomCyclicLogPMF(): The condition `probSuccess <= 1` must hold. probSuccess = "//getStr(probSuccess))
        call setGeomCyclicLogPMF(logPMF, stepSuccess, log(probSuccess), period)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGeomCyclicLogPMF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG), parameter :: logProbSuccessMin = log(epsilon(0._RKG))
        real(RKG), parameter :: LOGNEGINF = -log(huge(logPMF))
        real(RKG) :: logDenominator
#if     Def_ENABLED
        real(RKG) :: logProbFailure
#endif
        ! Validate the input.
        CHECK_ASSERTION(__LINE__, all(0_IK < [stepSuccess]), SK_"@setGeomCyclicLogPMF(): The condition `all(0 < [stepSuccess])` must hold. stepSuccess = "//getStr(stepSuccess))
        CHECK_ASSERTION(__LINE__, all([stepSuccess] <= period), SK_"@setGeomCyclicLogPMF(): The condition `all([stepSuccess] <= period)` must hold. stepSuccess, period = "//getStr([real(RKG) :: stepSuccess, period]))
        CHECK_ASSERTION(__LINE__, logProbSuccess <= 0._RKG, SK_"@setGeomCyclicLogPMF(): The condition `logProbSuccess <= 0.` must hold. logProbSuccess = "//getStr(logProbSuccess))
        ! Compute the PMF.
        if (logProbSuccess < 0._RKG) then ! imperfect probability of success.
            if (logProbSuccessMin < logProbSuccess) then
#if             Log_ENABLED
                CHECK_ASSERTION(__LINE__, abs(1._RKG - exp(logProbSuccess) - exp(logProbFailure)) < epsilon(0._RKG) * 100, SK_"@setGeomCyclicLogPMF(): The condition `exp(logProbFailure) + exp(logProbSuccess) == 1.` must hold. logProbFailure, logProbSuccess = "//getStr([logProbFailure, logProbSuccess]))
#elif           Def_ENABLED
                logProbFailure = log(get1mexp(logProbSuccess))
#else
#error          "Unrecognized interface."
#endif
                logDenominator = log(get1mexp(period * logProbFailure))
                logPMF = logProbSuccess + (stepSuccess - 1_IK) * logProbFailure - logDenominator
            else
                logPMF = LOGNEGINF
            end if
        else ! 100% probability of success.
#if         D0_ENABLED
            if (1_IK < stepSuccess) then ! stepSuccess can be larger than 1.
                logPMF = LOGNEGINF
            else
                logPMF = 0._RKG
            end if
#elif       D1_ENABLED
            block
                integer(IK) :: istep
                do concurrent(istep = 1_IK : size(stepSuccess, 1, IK))
                    if (1_IK < stepSuccess(istep)) then
                        logPMF(istep) = -log(huge(logPMF))
                    else
                        logPMF(istep) = 0._RKG
                    end if
                end do
            end block
#endif
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getGeomCyclicLogCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0._RKG < probSuccess, SK_"@setGeomCyclicLogCDF(): The condition `0 <= probSuccess` must hold. probSuccess = "//getStr(probSuccess))
        CHECK_ASSERTION(__LINE__, probSuccess <= 1._RKG, SK_"@setGeomCyclicLogCDF(): The condition `probSuccess <= 1` must hold. probSuccess = "//getStr(probSuccess))
        call setGeomCyclicLogCDF(logCDF, stepSuccess, log(probSuccess), period)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGeomCyclicLogCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: logDenominator
#if     Def_ENABLED
        real(RKG) :: logProbFailure
#endif
        ! Validate the input.
        CHECK_ASSERTION(__LINE__, all(0_IK < [stepSuccess]), SK_"@setGeomCyclicLogPMF(): The condition `all(0 < [stepSuccess])` must hold. stepSuccess = "//getStr(stepSuccess))
        CHECK_ASSERTION(__LINE__, all([stepSuccess] <= period), SK_"@setGeomCyclicLogPMF(): The condition `all([stepSuccess] <= period)` must hold. stepSuccess, period = "//getStr([real(RKG) :: stepSuccess, period]))
        CHECK_ASSERTION(__LINE__, logProbSuccess <= 0._RKG, SK_"@setGeomCyclicLogPMF(): The condition `logProbSuccess <= 0.` must hold. logProbSuccess = "//getStr(logProbSuccess))
        ! Compute the CDF.
        if (logProbSuccess < 0._RKG) then ! imperfect probability of success.
#if         Log_ENABLED
            CHECK_ASSERTION(__LINE__, abs(1._RKG - exp(logProbSuccess) - exp(logProbFailure)) < epsilon(0._RKG) * 100, SK_"@setGeomCyclicLogPMF(): The condition `exp(logProbFailure) + exp(logProbSuccess) == 1.` must hold. logProbFailure, logProbSuccess = "//getStr([logProbFailure, logProbSuccess]))
#elif       Def_ENABLED
            logProbFailure = log(get1mexp(logProbSuccess))
#else
#error      "Unrecognized interface."
#endif
            logDenominator = log(get1mexp(logProbFailure * period))
            logCDF = log(get1mexp(logProbFailure * stepSuccess)) - logDenominator
        else ! 100% probability of success.
            logCDF = 0._RKG
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getGeomCyclicRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0._RKG < probSuccess, SK_"@getGeomCyclicRand(): The condition `0 <= probSuccess` must hold. probSuccess = "//getStr(probSuccess))
        if (probSuccess < 1._RKG) then
            call setGeomCyclicRand(rand, log(1._RKG - probSuccess), period)
        else
            rand = 1_IK
            CHECK_ASSERTION(__LINE__, 0_IK < period, SK_"@getGeomCyclicRand(): The condition `0 < period` must hold. probSuccess = "//getStr(period))
            CHECK_ASSERTION(__LINE__, probSuccess <= 1._RKG, SK_"@getGeomCyclicRand(): The condition `probSuccess <= 1` must hold. probSuccess = "//getStr(probSuccess))
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGeomCyclicRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        if (1_IK < period) then
#if         RNGD_ENABLED
            call setGeomRand(rand, logProbFailure)
#elif       RNGF_ENABLED || RNGX_ENABLED
            call setGeomRand(rng, rand, logProbFailure)
#else
#error      "Unrecognized interface."
#endif
#if         D0_ENABLED
            rand = mod(rand, period)
            if (rand == 0_IK) rand = period
#elif       D1_ENABLED
            block
                integer(IK) :: i
                do concurrent(i = 1 : size(rand, 1, IK))
                    rand(i) = mod(rand(i), period)
                    if (rand(i) == 0_IK) rand(i) = period
                end do
            end block
#endif
        else
            rand = 1_IK
            CHECK_ASSERTION(__LINE__, 0_IK < period, SK_"@setGeomCyclicRand(): The condition `0 < period` must hold. probSuccess = "//getStr(period))
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   isFailedGeomCyclicFit_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: niter
        real(RKG) :: sumFreq
        real(RKG) :: logProbFailure
        real(RKG) :: logProbSuccess
        real(RKG) :: diff(size(stepSuccess, 1, IK))
        real(RKG) :: logFreqSuccess(size(freqSuccess, 1, IK))
        real(RKG), parameter :: small = sqrt(epsilon(0._RKG))
        integer(IK), parameter :: niter_def = int(1000 * precision(0._RKG) / 53., IK)
        CHECK_ASSERTION(__LINE__, all(0_IK < freqSuccess), SK_"@isFailedGeomCyclicFit(): The condition `all(0 < freqSuccess)` must hold. freqSuccess = "//getStr(freqSuccess))
        CHECK_ASSERTION(__LINE__, all(0_IK < stepSuccess), SK_"@isFailedGeomCyclicFit(): The condition `all(0 < stepSuccess)` must hold. stepSuccess = "//getStr(stepSuccess))
        CHECK_ASSERTION(__LINE__, 1_IK < size(stepSuccess, 1, IK), SK_"@isFailedGeomCyclicFit(): The condition `1 < size(stepSuccess)` must hold. size(stepSuccess) = "//getStr(size(stepSuccess, 1, IK)))
        CHECK_ASSERTION(__LINE__, maxval(stepSuccess) <= period, SK_"@isFailedGeomCyclicFit(): The condition `maxval(stepSuccess) <= period` must hold. maxval(stepSuccess), period = "//getStr([real(RKG) :: maxval(stepSuccess), period]))
        CHECK_ASSERTION(__LINE__, size(stepSuccess, 1, IK) <= period, SK_"@isFailedGeomCyclicFit(): The condition `size(stepSuccess) <= period` must hold. size(stepSuccess), period = "//getStr([size(stepSuccess, 1, IK), period]))
        CHECK_ASSERTION(__LINE__, size(stepSuccess, 1, IK) == size(freqSuccess, 1, IK), SK_"@isFailedGeomCyclicFit(): The condition `size(stepSuccess) == size(freqSuccess)` must hold. size(stepSuccess), size(freqSuccess) = "//getStr([size(stepSuccess, 1, IK), size(freqSuccess, 1, IK)]))

        ! First find the max likelihood.

        niter = niter_def
        sumFreq = sum(freqSuccess)
        !block
        !    real(RKG) :: xmin(1)
        !    xmin = .5_RKG
        !    failed = isFailedMinPowell(getNegLogLike, xmin)
        !    probSuccess = xmin(1)
        !end block
        probSuccess = getMinBrent(getNegLogLike, xlow = small, xupp = 1._RKG - small, niter = niter)
        failed = niter_def < niter
        if (failed) return

        ! Find the normalization factor.

        logProbFailure = log(1._RKG - probSuccess)
        logProbSuccess = log(probSuccess)

        niter = niter_def
        call setGeomCyclicLogPMF(diff, stepSuccess, logProbSuccess, logProbFailure, period)
        logFreqSuccess = log(real(freqSuccess, RKG))
        diff = logFreqSuccess - diff
        normFac = exp(getMinBrent(getSumDistSq, niter = niter))
        failed = niter_def < niter

    contains

        PURE function getNegLogLike(probSuccess) result(negLogLike)
            real(RKG), intent(in) :: probSuccess!(1)
            real(RKG) :: negLogLike
            if (0._RKG < probSuccess .and. probSuccess < 1._RKG) then
                !negLogLike = -sum(getGeomCyclicLogPMF(stepSuccess, probSuccess, period))
                negLogLike = (1 - probSuccess)**period
                if (negLogLike < 1) then
                    negLogLike = sumFreq * (log(1 - negLogLike) - log(probSuccess)) + sum(freqSuccess * (1 - stepSuccess)) * log(1 - probSuccess)
                else
                    negLogLike = sqrt(huge(0._RKG))
                end if
            else
                error stop "This is a miracle! Please inform the library developers at: https://github.com/cdslaborg/paramonte"
                !negLogLike = sqrt(huge(0._RKG))
            end if
        end function

        PURE function getSumDistSq(logNormFac) result(sumDistSq)
            real(RKG), intent(in) :: logNormFac
            real(RKG) :: sumDistSq
            sumDistSq = sum((diff - logNormFac)**2)
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   isFailedGeomCyclicFit_ENABLED && Def_ENABLED && 0
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, counter
        real(RKG) :: logFreqSuccess(size(freqSuccess, 1, IK))
        integer(IK) :: stepSuccess(size(freqSuccess, 1, IK))
        CHECK_ASSERTION(__LINE__, size(stepSuccess, 1, IK) == size(freqSuccess, 1, IK), SK_"@isFailedGeomCyclicFit(): The condition `size(stepSuccess) == size(freqSuccess)` must hold. size(stepSuccess), size(freqSuccess) = "//getStr([size(stepSuccess, 1, IK), size(freqSuccess, 1, IK)]))
        counter = 0_IK
        do i = 1, size(stepSuccess, 1, IK)
            if (0 < freqSuccess(i)) then
                counter = counter + 1_IK
                logFreqSuccess(counter) = log(real(freqSuccess(counter), RKG))
                stepSuccess(counter) = stepSuccess(i)
            else
                CHECK_ASSERTION(__LINE__, 0 == freqSuccess(i), SK_"@isFailedGeomCyclicFit(): The condition `all(0 <= freqSuccess)` must hold. i, freqSuccess(i) = "//getStr([i, freqSuccess(i)]))
            end if
        end do
        failed = isFailedGeomCyclicFit(stepSuccess(1 : counter), logFreqSuccess(1 : counter), period, probSuccess, normFac)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   isFailedGeomCyclicFit_ENABLED && Per_ENABLED && 0
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: successProbFisherTransLogNormFac(2) ! vector of the two parameters to fit.
        real(RKG), parameter :: SUCCESS_PROB_INIT_GUESS = 0.23_RKG
        real(RKG), parameter :: FISHER_TRANS_SUCCESS_PROB_INIT_GUESS = atanh(2 * (SUCCESS_PROB_INIT_GUESS - 0.5_RKG))
        CHECK_ASSERTION(__LINE__, all(0_IK < stepSuccess), SK_"@isFailedGeomCyclicFit(): The condition `all(0 < stepSuccess)` must hold. stepSuccess = "//getStr(stepSuccess))
        CHECK_ASSERTION(__LINE__, all(0._RKG <= logFreqSuccess), SK_"@isFailedGeomCyclicFit(): The condition `all(0 < stepSuccess)` must hold. stepSuccess = "//getStr(logFreqSuccess))
        CHECK_ASSERTION(__LINE__, 1_IK < size(stepSuccess, 1, IK), SK_"@isFailedGeomCyclicFit(): The condition `1 < size(stepSuccess)` must hold. size(stepSuccess) = "//getStr(size(stepSuccess, 1, IK)))
        CHECK_ASSERTION(__LINE__, maxval(stepSuccess) <= period, SK_"@isFailedGeomCyclicFit(): The condition `maxval(stepSuccess) <= period` must hold. maxval(stepSuccess), period = "//getStr([real(RKG) :: maxval(stepSuccess), period]))
        CHECK_ASSERTION(__LINE__, size(stepSuccess, 1, IK) <= period, SK_"@isFailedGeomCyclicFit(): The condition `size(stepSuccess) <= period` must hold. size(stepSuccess), period = "//getStr([size(stepSuccess, 1, IK), period]))
        CHECK_ASSERTION(__LINE__, size(stepSuccess, 1, IK) == size(logFreqSuccess, 1, IK), SK_"@isFailedGeomCyclicFit(): The condition `size(stepSuccess) == size(logFreqSuccess)` must hold. size(stepSuccess), size(logFreqSuccess) = "//getStr([size(stepSuccess, 1, IK), size(logFreqSuccess, 1, IK)]))

        ! Do Fisher transformation to ensure stability.
        successProbFisherTransLogNormFac(2) = 0._RKG
        successProbFisherTransLogNormFac(1) = FISHER_TRANS_SUCCESS_PROB_INIT_GUESS
        failed = isFailedMinPowell(getSumDistSq, successProbFisherTransLogNormFac)
        if (.not. failed) then
            probSuccess = 0.5_RKG * tanh(successProbFisherTransLogNormFac(1)) + 0.5_RKG ! reverse Fisher-transform
            normFac = exp(successProbFisherTransLogNormFac(2))
        end if

    contains

        PURE function getSumDistSq(successProbFisherTransLogNormFac) result(sumDistSq)
            real(RKG), intent(in) :: successProbFisherTransLogNormFac(2)
            real(RKG) :: logPMF(size(stepSuccess, 1, IK)), logProbSuccess, sumDistSq
            logProbSuccess = log(0.5_RKG * tanh(successProbFisherTransLogNormFac(1)) + 0.5_RKG)
            call setGeomCyclicLogPMF(logPMF, stepSuccess, logProbSuccess, period)
            sumDistSq = sum((logFreqSuccess - logPMF - successProbFisherTransLogNormFac(2))**2)
        end function getSumDistSq
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif