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
!>  This file contains the implementation details of the 2D routines under the generic interface [pm_sampleACT](@ref pm_sampleACT).
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Saturday 4:40 PM, August 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getACT_ENABLED && D1_ENABLED && DEF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        act = getACT(seq, batchMeans)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getACT_ENABLED && D1_ENABLED && BMM_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: smin, smax, step, isize, lenseq
        smin = method%smin
        step = method%step
        smax = method%smax
        lenseq = size(seq, 1, IK)
        CHECK_ASSERTION(__LINE__, 0 < lenseq, SK_"@getACT(): The condition `0 < size(seq)` must hold. size(seq) = "//getStr(lenseq))
        CHECK_ASSERTION(__LINE__, 1_IK < smin, SK_"@getACT(): The condition `1 < method%smin` must hold. method%smin = "//getStr(smin))
        CHECK_ASSERTION(__LINE__, smax == huge(0_IK) .or. (smax <= smin .and. smax <= lenseq / 2_IK), SK_"@getACT(): The condition `method%smin <= method%smax <= size(seq, 1) / 2` must hold. method%smax, size(seq) = "//getStr([smax, lenseq]))
        CHECK_ASSERTION(__LINE__, 0_IK < step, SK_"@getACT(): The condition `0 < method%step` must hold. method%step = "//getStr(step))
        if (1_IK < lenseq) then
            smax = min(smax, lenseq / 2_IK)
            act = -huge(act)
            do isize = smin, smax, step
#if             ONE_ENABLED
                act = max(act, getACT(seq, batchMeans_type(isize)))
#elif           WTI_ENABLED
                act = max(act, getACT(seq, weight, batchMeans_type(isize)))
#else
#error          "Unrecognized interface."
#endif
            end do
        else
            act = lenseq
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getACT_ENABLED && D1_ENABLED && BMD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: batchSize, lenseq, weisum
        real(TKC), allocatable :: batchMean(:)
        real(TKC) :: batchSizeInv, avgSeq, sumDiffSq, varBatchMean, avgBatchMean
        integer(IK) :: nbatch, ibatch, seqLenEffective
#if     WTI_ENABLED
        real(TKC) :: diffSq
        integer(IK) :: cumSumWei(size(seq, 1, IK)), currentBatchEndLoc, iseq, iseqVerbose
        lenseq = size(seq, 1, IK)
        call setCumSum(cumSumWei, weight)
        CHECK_ASSERTION(__LINE__, size(weight, 1, IK) == lenseq, SK_"@getACT(): The condition `size(weight) == size(seq)` must hold. size(weight), size(seq) = "//getStr([size(weight, 1, IK), lenseq]))
        weisum = cumSumWei(lenseq)
#elif   ONE_ENABLED
        integer(IK) :: iseqBeg, iseqEnd
        lenseq = size(seq, 1, IK)
        weisum = lenseq
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, 0 < lenseq, SK_"@getACT(): The condition `0 < size(seq)` must hold. size(seq) = "//getStr(lenseq))
        if (lenseq < 2_IK) then
            act = lenseq
            return
        end if

        ! Compute batch size and count.

        act = 1._TKC
        batchSize = method%size
        if (batchSize == 0_IK) batchSize = int(real(weisum)**real(2./3.), IK)
        batchSizeInv = 1._TKC / real(batchSize, TKC)
        nbatch = weisum / batchSize
        CHECK_ASSERTION(__LINE__, 0_IK < batchSize, SK_"@getACT(): The condition `0 < method%size` must hold. method%size = "//getStr(batchSize))
        if (nbatch < 2_IK) then
            !check_assertion(__LINE__, method%size == 0_IK, SK_"@getACT(): The condition `method%size * 2 < sum(weight)` must hold (there must be more than one batches to compute ACT). method%size, sum(weight) = "//getStr([batchSize, weisum]))
            act = weisum
            return
        end if
        seqLenEffective = nbatch * batchSize
        ! xxx: here goes another GFortran 7.3 bug: `batchMean` is assumed already allocated, despite the first appearance here.
        if (allocated(batchMean)) deallocate(batchMean)
        allocate(batchMean(nbatch))

        ! \todo
        ! iterate from the end to the beginning of the chain to ignore initial points instead of the last points.
        ! This would be beneficial for MCMC samples.

        ! Compute the Batch means vector and mean of the whole (weighted) sequence.

        sumDiffSq = 0._TKC
        avgSeq = 0._TKC
#if     WTI_ENABLED
        iseq = 1
        ibatch = 1
        iseqVerbose = 0
        currentBatchEndLoc = batchSize
        batchMean(ibatch) = 0._TKC
        loopOverWeight: do
            iseqVerbose = iseqVerbose + 1
            if (cumSumWei(iseq) < iseqVerbose) iseq = iseq + 1
            if (currentBatchEndLoc < iseqVerbose) then ! we are done with the current batch
                avgSeq = avgSeq + batchMean(ibatch)
                if (seqLenEffective < iseqVerbose) exit loopOverWeight  ! condition equivalent to currentBatchEndLoc == seqLenEffective.
                currentBatchEndLoc = currentBatchEndLoc + batchSize
                ibatch = ibatch + 1
                batchMean(ibatch) = 0._TKC
            end if
            batchMean(ibatch) = batchMean(ibatch) + seq(iseq)
        end do loopOverWeight
        batchMean = batchMean * batchSizeInv
        avgSeq = avgSeq / real(seqLenEffective, TKC) ! whole sequence mean.
        ! compute whole sequence variance.
        iseq = 1
        iseqVerbose = 0
        diffSq = (seq(iseq) - avgSeq)**2
        loopComputeVariance: do
            iseqVerbose = iseqVerbose + 1
            if (seqLenEffective < iseqVerbose) exit loopComputeVariance
            if (cumSumWei(iseq) < iseqVerbose) then
                iseq = iseq + 1 ! by definition, iseq never become > lenseq, otherwise it leads to disastrous errors
                diffSq = (seq(iseq) - avgSeq)**2
            end if
            sumDiffSq = sumDiffSq + diffSq ! whole sequence variance.
        end do loopComputeVariance
#elif   ONE_ENABLED
        iseqBeg = 1
        iseqEnd = 0
        do ibatch = 1, nbatch
            iseqEnd = iseqEnd + batchSize
            batchMean(ibatch) = sum(seq(iseqBeg : iseqEnd))
            avgSeq = avgSeq + batchMean(ibatch)
            batchMean(ibatch) = batchMean(ibatch) * batchSizeInv
            iseqBeg = iseqEnd + 1_IK
        end do
        avgSeq = avgSeq / real(seqLenEffective, TKC) ! whole sequence mean.
        sumDiffSq = sum((seq(1 : seqLenEffective) - avgSeq)**2) ! whole sequence variance.

#else
#error  "Unrecognized interface."
#endif
        avgBatchMean = sum(batchMean) / real(nbatch, TKC) ! batch means
        varBatchMean = sum((batchMean - avgBatchMean)**2) / real(nbatch - 1, TKC) ! batch variances
        !act = varBatchMean * seqLenEffective * (seqLenEffective - 1) / sumDiffSq
        !act = batchSize * varBatchMean * (seqLenEffective - 1) / sumDiffSq
        act = varBatchMean * seqLenEffective * (seqLenEffective - 1) / sumDiffSq / real(nbatch, TKC)
        !print *, avgSeq, sumDiffSq
        !print *, sum(batchMean) / real(nbatch, TKC), varBatchMean
        !print *, nbatch, batchSize
        !print *, act

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getACT_ENABLED && D1_ENABLED && (CSD_ENABLED || CSM_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKC) :: mean, sample(size(seq, 1, IK))
        integer(IK) :: lenseq, weisum
        lenseq = size(seq, 1, IK)
        CHECK_ASSERTION(__LINE__, 1 < lenseq, SK_"@getACT(): The condition `0 < size(seq)` must hold. size(seq) = "//getStr(lenseq))
        if (lenseq < 2_IK) then
            act = lenseq
            return
        end if
#if     WTI_ENABLED
        CHECK_ASSERTION(__LINE__, size(weight, 1, IK) == size(seq, 1, IK), SK_"@getACT(): The condition `size(weight) == size(seq)` must hold. size(weight), size(seq) = "//getStr([size(weight, 1, IK), size(seq, 1, IK)]))
        weisum = sum(weight)
        sample = getVerbose(seq, weight, weisum)
        mean = sum(sample * weight) / real(weisum, TKC)
        sample = sample - mean
#elif   ONE_ENABLED
        weisum = size(seq, 1, IK)
        mean = sum(seq) / real(weisum, TKC)
        sample = seq - mean
#else
#error  "Unrecognized interface."
#endif
#if     CSD_ENABLED
        block
            integer(IK) :: i
            real(TKC) :: signif
            signif = real(method%signif, TKC)
            CHECK_ASSERTION(__LINE__, 0 <= method%signif, SK_"@getACT(): The condition `0 <= method%signif` must hold. method%signif = "//getStr(signif))
            sample = getACF(sample, norm = stdscale)
            ! For autocorrelation, under the assumption of a completely random series, the ACF standard error reduces to `sqrt(1 / ndata)`.
            act = 0._TKC
            signif = signif * sqrt(1._TKC / weisum)
            do i = 1, size(sample, 1, IK)
                if (sample(i) < signif) exit
                act = act + sample(i)
            end do
            act = 2 * act - 1
        end block
#elif   CSM_ENABLED
        sample = getACF(sample, norm = stdscale)
        act = 2 * maxval(getCumSum(sample)) - 1
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  TYPE_OF_SEQ
#undef  GET_ABSQ
#undef  GET_RE