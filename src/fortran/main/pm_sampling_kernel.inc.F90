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

#if         CAF_ENABLED || MPI_ENABLED || OMP_ENABLED
#define     SET_PARALLEL(X)X
#else
#define     SET_PARALLEL(X)
#endif
#if         CFI_ENABLED
#define     SHAPE_ARG, ndim
#else
#define     SHAPE_ARG
#endif
            !   \warning
            !   The MPI library real kinds do not fully match all available Fortran kinds.
            !   As such, instead of using the MPI intrinsic explicit real kinds
            !   we pass real data as `mpi_byte` data type.
            !   This method, however, has a caveat when MPI communicates data between systems with different endianness.
            !   When data is transferred as bytes, the endianness is not automatically taken into account.
            !   As such, the transferred values may not be meaningful.
            !   For now, we accept this risk as the odds of mixed endianness on the existing contemporary supercomputers is low.
            !   This may, however, need a more robust solution in the future.
            !   One possible solution might be to build an appropriate real data type using the MPI intrinsic `mpi_type_create_struct()`.
            integer             , parameter     :: NBRK = c_sizeof(0._RKC) ! Number of bytes in the current real kind.
#if         OMP_ENABLED && (MATLAB_ENABLED || PYTHON_ENABLED || R_ENABLED)
            real(RKC)                           :: mold
#elif       OMP_ENABLED
            integer(IK)         , parameter     :: CACHELINE = 128_IK * 8_IK
            integer(IK)         , parameter     :: REAL_SIZE = storage_size(0._RKC)
            integer(IK)         , parameter     :: JUMP = 1_IK + CACHELINE / REAL_SIZE
            real(RKC)                           :: logFuncState(JUMP, spec%image%count)
#endif
            character(*, SK)    , parameter     :: PROCEDURE_NAME = MODULE_NAME//SK_"@getErrKernelRun()"
            integer(IK)         , parameter     :: CHAIN_RESTART_OFFSET = 2_IK
            real(RKC)                           :: sumAccrAccRejDel ! sum(accr since start) for accepted/rejected/delayed proposals: used to figure out the average acceptance ratio for the entire chain.
            real(RKC)                           :: sumAccrAccRej    ! sum(accr since start) for accepted/rejected proposals: used to figure out the average acceptance ratio for the entire chain.
            integer             , allocatable   :: pos(:)
            integer                             :: ndimp1
            integer(IK)                         :: ndim
            integer(IK)                         :: numFunCallAcceptedLastAdaptation                     ! number of function calls accepted at Last proposal adaptation occurrence
            integer(IK)                         :: counterPAP                                           ! counter for proposalAdaptationPeriod
            integer(IK)                         :: counterPAC                                           ! counter for proposalAdaptationCount
            integer(IK)                         :: counterDRS                                           ! counter for Delayed Rejection Stages
            integer(IK)                         :: lastStateWeight                                      ! This is used for passing the most recent verbose chain segment to the adaptive updater of the sampler
            integer(IK)                         :: currentStateWeight                                   ! counter for SampleWeight, used only in restart mode
            integer(IK)                         :: numFunCallAcceptedPlusOne                            ! counter for SampleWeight, used only in restart mode
            real(RKC)                           :: meanAccRateSinceStart                                ! used for restart file read. xxx \todo: this could be merged with stat%cfc%meanAcceptanceRate.
            real(RKC)                           :: logFuncDiff                                          ! The difference between the log of the old and the new states. Used to avoid underflow.
            integer(IK)                         :: dumint
            logical(LK)                         :: proposalAdaptationGreedyEnabled
            integer(IK)                         :: proposalAdaptationLen
            logical(LK)                         :: proposalAdaptationSampleUsed
            integer(IK)                         :: acceptedRejectedDelayedUnusedRestartMode
           !real(RKC)                           :: proposalAdaptationDummy
            real(RKC)           , allocatable   :: proposalAdaptation(:)
            integer                             :: imageStartID, imageEndID, imageID ! This must be default integer, at least for MPI.
#if         OMP_ENABLED
#define     CO_LOGFUNCSTATE(I,J,K)co_logFuncState(I,J,K)
#define     CO_ACCR(I,J)co_accr(I,J)
#define     GET_ELL(X,I)X(I)
#define     SET_CONOMP(X)X
            integer                             :: co_proposalFound_samplerUpdateOccurred(2)            ! merging these scalars would reduce the MPI communication overhead cost: co_proposalFound, co_samplerUpdateOccurred, co_counterDRS, 0 means false, 1 means true
            real(RKC)                           :: maxLogFuncRejectedProposal(spec%image%count)         ! used in delayed rejection sampling.
            real(RKC)                           :: unifrnd(spec%image%count)                            ! used for random number generation.
            real(RKC)           , allocatable   :: co_logFuncState(:,:,:)                               ! (0 : ndim, 1 : njob, -1 : proposalDelayedRejectionCount), -1 is the current accepted state.
            real(RKC)           , allocatable   :: co_accr(:,:)                                         ! (1 : njob, -1 : proposalDelayedRejectionCount)
#else
#define     CO_LOGFUNCSTATE(I,J,K)co_logFuncState(I,K)
#define     CO_ACCR(I,J)co_accr(J)
#define     GET_ELL(X,I)X
#define     SET_CONOMP(X)
            real(RKC)                           :: maxLogFuncRejectedProposal                           ! used in delayed rejection sampling.
            real(RKC)                           :: unifrnd                                              ! used for random number generation.
#if         CAF_ENABLED
            integer             , save          :: co_proposalFound_samplerUpdateOccurred(2)[*]         ! merging these scalars would reduce the MPI communication overhead cost: co_proposalFound, co_samplerUpdateOccurred, co_counterDRS, 0 means false, 1 means true.
            real(RKC)           , allocatable   :: co_logFuncState(:,:)[:]                              ! (0 : ndim, -1 : proposalDelayedRejectionCount), -1 is the current accepted state.
            real(RKC)           , allocatable   :: co_accr(:)[:]                                        ! (1 : njob, -1 : proposalDelayedRejectionCount)
#else
            integer                             :: co_proposalFound_samplerUpdateOccurred(2)            ! merging these scalars would reduce the MPI communication overhead cost: co_proposalFound, co_samplerUpdateOccurred, co_counterDRS, 0 means false, 1 means true
            real(RKC)           , allocatable   :: co_logFuncState(:,:)                                 ! (0 : ndim, -1 : proposalDelayedRejectionCount), -1 is the current accepted state.
            real(RKC)           , allocatable   :: co_accr(:)                                           ! (1 : njob, -1 : proposalDelayedRejectionCount)
#endif
#if         CAF_ENABLED || MPI_ENABLED
#if         MPI_ENABLED
            integer                             :: proposalDelayedRejectionCountPlusTwo
            integer                             :: proposalFoundSinglChainModeReduced                   ! the reduced value by summing proposalFoundSinglChainModeReduced over all images.
            real(RKC)           , allocatable   :: accrMat(:,:)                                         ! matrix of size (-1:spec%proposalDelayedRejectionCount%val,1:spec%image%count).
            integer                             :: ierrMPI
#endif
            integer                             :: proposalFoundSinglChainMode                          ! used in singleChain delayed rejection. zero if the proposal is not accepted. 1 if the proposal is accepted.
#endif
#endif
#if         !(CAF_ENABLED || MPI_ENABLED || OMP_ENABLED)
            parameter(imageID = 1, imageStartID = 1, imageEndID = 1) ! needed even in the case of serial run to assign a proper value to stat%cfc%processID
#endif
            stat%avgCommPerFunCall = 0._RKD ! Until the reporting time, this is in reality, sumCommTimePerFunCall. This is meaningful only in singleChain parallelism.
            stat%numFunCallAcceptedRejectedDelayedUnused = spec%image%count ! used only in singleChain parallelism, and relevant only on the first image.
            ndim = spec%ndim%val
            ndimp1 = int(ndim + 1)
            err%occurred = .false._LK
            err%msg = repeat(SK_" ", 255)
            acceptedRejectedDelayedUnusedRestartMode = 0_IK ! used to compute more accurate timings in the restart mode.
            stat%avgTimePerFunCall = 0._RKD
            if (allocated(co_accr)) deallocate(co_accr)
            if (allocated(co_logFuncState)) deallocate(co_logFuncState)
#if         CAF_ENABLED
            allocate(co_logFuncState(0 : ndim, -1 : spec%proposalDelayedRejectionCount%val)[*])
            allocate(co_accr(-1 : spec%proposalDelayedRejectionCount%val)[*]) ! the negative element will contain counterDRS.
#else
            allocate(CO_LOGFUNCSTATE(0 : ndim, 1 : spec%image%count, -1 : spec%proposalDelayedRejectionCount%val))
            allocate(CO_ACCR(1 : spec%image%count, -1 : spec%proposalDelayedRejectionCount%val)) ! the negative element will contain counterDRS.
#endif
            CO_ACCR(:, 0) = 1._RKC ! initial acceptance rate for the first zeroth DR stage.
            CO_ACCR(:, -1) = 0._RKC ! the real-valued counterDRS, indicating the initial delayed rejection stage at which the first point is sampled.
            CO_ACCR(:, 1 : spec%proposalDelayedRejectionCount%val) = 0._RKC ! indicates the very first proposal acceptance on image 1.
#if         MPI_ENABLED
            if (allocated(accrMat)) deallocate(accrMat)
            allocate(accrMat(-1 : spec%proposalDelayedRejectionCount%val, 1 : spec%image%count)) ! the negative element will contain counterDRS.
            accrMat = 0._RKC ! -huge(1._RKC) ! debug
            accrMat(0, 1 : spec%image%count) = 1._RKC ! initial acceptance rate for the first zeroth DR stage.
            proposalDelayedRejectionCountPlusTwo = (spec%proposalDelayedRejectionCount%val + 2) * NBRK
            !proposalDelayedRejectionCountPlusTwo = spec%proposalDelayedRejectionCount%val + 2
#endif
            if (proposal%delRejEnabled) then
                stat%numFunCallAcceptedRejectedDelayed = 0_IK ! Markov Chain counter
                sumAccrAccRejDel = 0._RKC ! sum of acceptance rate
            end if
            ! proposalAdaptation = 0._RKC ! needed for the first output.
            proposalAdaptationSampleUsed = .true._LK ! needed to set up lastStateWeight and numFunCallAcceptedLastAdaptation for the first accepted proposal.
            sumAccrAccRej = 0._RKC ! sum of acceptance rate.
            counterPAC = 0_IK ! counter for pproposalAdaptationCount.
            counterPAP = 0_IK ! counter for proposalAdaptationPeriod.
            stat%numFunCallAccepted = 0_IK ! Markov Chain acceptance counter.
            stat%numFunCallAcceptedRejected = 0_IK ! Markov Chain counter
            numFunCallAcceptedLastAdaptation = 0_IK
            lastStateWeight = -huge(lastStateWeight)
            meanAccRateSinceStart = 1._RKC ! needed for the first restart output in fresh run.

            proposalAdaptationLen = 100_IK
            blockNewDryRun: if (spec%run%is%new) then
                RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(isFailedChainResize(stat%cfc, ndim, spec%outputChainSize%val, err%msg)),SK_"Insufficient memory allocation space for the output chain file contents.")
            else blockNewDryRun
                ! Load the existing Chain file into stat%cfc components.
                err = getErrChainRead(stat%cfc, spec%chainFile%file, spec%outputChainFileFormat%val, pre = spec%outputChainSize%val, pos = pos)
                !err = getErrChainWrite(stat%cfc, spec%chainFile%file//".temp", spec%outputChainFileFormat%val, spec%chainFile%format%rows, 1_IK, stat%cfc%nsam, spec%proposalAdaptationPeriod%val)
                RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(err%occurred),err%msg)
                if (err%iswarned .and. spec%image%is%first) call spec%disp%warn%show(PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SK_": "//err%msg)
                if (0 < stat%cfc%nsam) proposalAdaptationLen = maxval(stat%cfc%sampleWeight(1 : stat%cfc%nsam))
                spec%run%is%new = stat%cfc%nsam <= CHAIN_RESTART_OFFSET
                spec%run%is%dry = .not. spec%run%is%new
                call spec%image%sync()
                ! set up the chain file.
                ! This chain file rewrite is necessary because the last two chains rows must be deleted.
                blockLeaderChainSetup: if (spec%image%is%leader) then
                    block
                        !integer :: ibeg, iend
                        integer(IK) :: irow, iweight, nrow
                        nrow = merge(CHAIN_RESTART_OFFSET, stat%cfc%nsam, spec%run%is%dry)
                        if (spec%outputChainFileFormat%isCompact) then
                            do irow = 0, nrow - 1
                                backspace(spec%chainFile%unit, iostat = err%stat, iomsg = err%msg)
                                RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(err%stat /= 0),err%msg)
                            end do
                        elseif (spec%outputChainFileFormat%isVerbose) then
                            do irow = 0, nrow - 1
                                do iweight = 1, stat%cfc%sampleWeight(stat%cfc%nsam - irow)
                                    backspace(spec%chainFile%unit, iostat = err%stat, iomsg = err%msg)
                                    RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(err%stat /= 0),err%msg)
                                end do
                            end do
                        elseif (spec%outputChainFileFormat%isBinary .and. 0_IK < stat%cfc%nsam) then
                            !print *, pos(2:) - pos(:size(pos) - 1)
                            read(spec%chainFile%unit, pos = pos(stat%cfc%nsam - nrow + 1), iostat = err%stat, iomsg = err%msg)
                            RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(err%stat /= 0),err%msg)
                            ! Infer the record length, then jump two records backward.
                            !inquire(spec%chainFile%unit, pos = iend, iostat = err%stat, iomsg = err%msg)
                            !RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(err%occurred),err%msg)
                            !ibeg = iend
                            !irow = 0_IK
                            !do
                            !    ibeg = ibeg - 1_IK
                            !    RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(ibeg < 1_IK),SK_"Failed to infer the record length in the output binary chain file """//spec%chainFile%file//SK_"""")
                            !    read(spec%chainFile%unit, pos = ibeg, iostat = err%stat ) stat%cfc%processID            (stat%cfc%nsam) &
                            !                                                            , stat%cfc%delayedRejectionStage(stat%cfc%nsam) &
                            !                                                            , stat%cfc%meanAcceptanceRate   (stat%cfc%nsam) &
                            !                                                            , stat%cfc%proposalAdaptation   (stat%cfc%nsam) &
                            !                                                            , stat%cfc%burninLocation       (stat%cfc%nsam) &
                            !                                                            , stat%cfc%sampleWeight         (stat%cfc%nsam) &
                            !                                                            , stat%cfc%sampleLogFunc        (stat%cfc%nsam) &
                            !                                                            , stat%cfc%sampleState          (:, stat%cfc%nsam)
                            !    if (err%stat == 0) then
                            !        irow = irow + 1_IK
                            !        !print *, irow
                            !        if (irow <= nrow) cycle
                            !        read(spec%chainFile%unit, pos = ibeg, iostat = err%stat, iomsg = err%msg)
                            !        RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(err%stat /= 0),err%msg)
                            !        exit
                            !    end if
                            !end do
                            !!
                            !!if (0 < nrow) then
                            !!    inquire(spec%chainFile%unit, pos = ibeg, iostat = err%stat, iomsg = err%msg)
                            !!    RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(err%stat /= 0),err%msg)
                            !!    write(spec%chainFile%unit,iostat = err%stat, iomsg = err%msg) stat%cfc%processID            (stat%cfc%nsam) &
                            !!                                                                , stat%cfc%delayedRejectionStage(stat%cfc%nsam) &
                            !!                                                                , stat%cfc%meanAcceptanceRate   (stat%cfc%nsam) &
                            !!                                                                , stat%cfc%proposalAdaptation   (stat%cfc%nsam) &
                            !!                                                                , stat%cfc%burninLocation       (stat%cfc%nsam) &
                            !!                                                                , stat%cfc%sampleWeight         (stat%cfc%nsam) &
                            !!                                                                , stat%cfc%sampleLogFunc        (stat%cfc%nsam) &
                            !!                                                                , stat%cfc%sampleState          (:, stat%cfc%nsam)
                            !!    RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(err%stat /= 0),err%msg)
                            !!    inquire(spec%chainFile%unit, pos = iend, iostat = err%stat, iomsg = err%msg)
                            !!    RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(err%stat /= 0),err%msg)
                            !!    write(spec%chainFile%unit, pos = iend - (nrow + 1) * (iend - ibeg), iostat = err%stat, iomsg = err%msg)
                            !!    RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(err%stat /= 0),err%msg)
                            !!end if
                        end if
                        !! Create a temporary copy of the chain file, just for the sake of not losing the simulation results if the write action fails at any stage.
                        !character(:, SK), allocatable :: contents, tempfile
                        !call setContentsFrom(spec%chainFile%file, contents, iostat = iostat, iomsg = err%msg)
                        !RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(iostat /= 0_IK),err%msg)
                        !tempfile = spec%chainFile%file//SK_".restart"
                        !call setContentsTo(tempfile, contents, iostat = iostat, iomsg = err%msg)
                        !RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(iostat /= 0_IK),err%msg)
                        !! reopen the chain file to write contents minus the last two rows and resume the simulation.
                        !err = getErrChainWrite(stat%cfc, spec%chainFile%file, spec%outputChainFileFormat%val, spec%chainFile%format%rows, 1_IK, stat%cfc%nsam - CHAIN_RESTART_OFFSET, spec%proposalAdaptationPeriod%val)
                        !if (err%iswarned .and. spec%image%is%leader) call spec%disp%warn%show(PROCEDURE_NAME//getFine(__LINE__)//SK_": "//err%msg)
                        !RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(err%occurred),err%msg)
                        !! remove the temporary copy of the chain file
                        !RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(isFailedRemove(tempfile, errmsg = err%msg)),err%msg)
                        !deallocate(contents, tempfile)
                    end block
                end if blockLeaderChainSetup
                !RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(err%occurred),err%msg)
            end if blockNewDryRun

            stat%timer = timer_type()
            !stat%timer%start = stat%timer%time()

            if (spec%image%is%leader) call setResized(proposalAdaptation, proposalAdaptationLen)
            if (spec%run%is%new) then ! this must be done separately from the above blockNewDryRun.
                stat%cfc%burninLocation(1) = 1_IK
                CO_LOGFUNCSTATE(1 : ndim, 1, 0) = spec%proposalStart%val ! proposal state.
#if             OMP_ENABLED && (MATLAB_ENABLED || PYTHON_ENABLED || R_ENABLED)
                co_logFuncState(1 : ndim, 2 : spec%image%count, 0) = spread(co_logFuncState(1 : ndim, 1, 0), 2, spec%image%count - 1) ! This defines all input states to the user-defined subroutine.
                block
                    real(RKD) :: avgTimePerFunCall, avgCommPerFunCall
                    !avgTimePerFunCall = epsilon(0._RKD); avgCommPerFunCall = epsilon(0._RKD)
                    mold = getLogFunc(co_logFuncState(:, 1 : spec%image%count, 0), avgTimePerFunCall, avgCommPerFunCall)
                    stat%avgTimePerFunCall = stat%avgTimePerFunCall + max(avgTimePerFunCall, epsilon(0._RKD))
                    stat%avgCommPerFunCall = stat%avgCommPerFunCall + max(avgCommPerFunCall, epsilon(0._RKD))
                end block
#else
                stat%timer%clock = stat%timer%time()
                CO_LOGFUNCSTATE(0, 1, 0) = getLogFunc(CO_LOGFUNCSTATE(1 : ndim, 1, 0)) ! proposal logFunc
                stat%avgTimePerFunCall = stat%avgTimePerFunCall + stat%timer%time() - stat%timer%clock
                SET_OMP(co_logFuncState(:, 2 : spec%image%count, 0) = spread(co_logFuncState(:, 1, 0), 2, spec%image%count - 1))
#endif
            else
                CO_LOGFUNCSTATE(1 : ndim, 1, 0) = stat%cfc%sampleState(1 : ndim, 1) ! proposal state
                CO_LOGFUNCSTATE(0, 1, 0) = stat%cfc%sampleLogFunc(1) ! proposal logFunc
                SET_CONOMP(co_logFuncState(:, 2 : spec%image%count, 0) = spread(co_logFuncState(:, 1, 0), 2, spec%image%count - 1))
            end if
            CO_LOGFUNCSTATE(0 : ndim, :, -1) = CO_LOGFUNCSTATE(0 : ndim, :, 0) ! set current logFunc and State equal to the first proposal.
            stat%chain%mode%val = -huge(stat%chain%mode%val)
            stat%chain%mode%loc = 0_IK
            SET_PARALLEL(if (spec%parallelism%is%singleChain) then)
                SET_PARALLEL(imageEndID = spec%image%count)
                SET_PARALLEL(imageStartID = 1)
            SET_PARALLEL(else) ! is%multiChain ! This never happens in the current OpenMP implementation.
#if             OMP_ENABLED
                error stop "Internal error occurred: multiChain parallelism is unsupported in threaded or concurrent simulations."
#endif
                SET_PARALLEL(imageStartID = spec%image%id)
                SET_PARALLEL(imageEndID = spec%image%id)
            SET_PARALLEL(end if)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start of loopMarkovChain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call initProgressReport(spec)

            loopMarkovChain: do

                co_proposalFound_samplerUpdateOccurred(2) = 0 ! at each iteration assume no samplerUpdateOccurred, unless it occurs.
                SET_CAFMPI(blockLeaderImage: if (spec%image%is%leader) then)

                    co_proposalFound_samplerUpdateOccurred(1) = 0 ! co_proposalFound = .false._LK
                    proposalAdaptationGreedyEnabled = counterPAC < spec%proposalAdaptationCountGreedy%val

                    SET_PARALLEL(imageID = imageStartID)
                    SET_PARALLEL(loopOverImages: do)

                        SET_PARALLEL(if (imageEndID < imageID) exit loopOverImages)
                        SET_CAF(stat%timer%clock = stat%timer%time())
                        SET_CAF(if (imageID /= spec%image%id) co_accr(-1 : spec%proposalDelayedRejectionCount%val) = co_accr(-1 : spec%proposalDelayedRejectionCount%val)[imageID]) ! happens only in is%singleChain holds.
                        SET_CAF(stat%avgCommPerFunCall = stat%avgCommPerFunCall + stat%timer%time() - stat%timer%clock)
                        SET_MPI(if (imageID /= spec%image%id) co_accr(-1 : spec%proposalDelayedRejectionCount%val) = accrMat(-1 : spec%proposalDelayedRejectionCount%val, imageID)) ! happens only in is%singleChain holds.
                        SET_CONOMP(if (imageID <= spec%image%count) co_accr(1, -1 : spec%proposalDelayedRejectionCount%val) = co_accr(imageID, -1 : spec%proposalDelayedRejectionCount%val))
                        ! This workaround is essential for efficient Coarray/MPI communications, but redundant for serial, openmp, or concurrent applications.
                        ! But who cares using a few more CPU cycles for serial mode in exchange for saving thousands in parallel?
                        counterDRS = nint(CO_ACCR(1, -1), IK)
                        if (-1_IK < counterDRS) co_proposalFound_samplerUpdateOccurred(1) = 1 ! co_proposalFound = .true._LK

                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% blockProposalAccepted %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        ! On the very first iteration, this block is (and must be) executed for imageID == 1
                        ! since it is for the first (starting) point, which is assumed to have been accepted
                        ! as the first point by the first coarray imageID.

                        blockProposalAccepted: if (co_proposalFound_samplerUpdateOccurred(1) == 1) then ! co_proposalAccepted = .true._LK

                            currentStateWeight = 0_IK

                            ! communicate the accepted logFunc and State from the winning image to the leader/all images: co_logFuncState.
#if                         CAF_ENABLED
                            if (imageID /= spec%image%id) then ! Avoid remote connection for something that is available locally.
                                stat%timer%clock = stat%timer%time()
                                co_logFuncState(0 : ndim, -1) = co_logFuncState(0 : ndim, -1)[imageID]
                                stat%avgCommPerFunCall = stat%avgCommPerFunCall + stat%timer%time() - stat%timer%clock
                            end if
#elif                       MPI_ENABLED
                            if (spec%parallelism%is%singleChain) then
                                stat%timer%clock = stat%timer%time()
                                ! broadcast winning image to all processes: buffer, count, datatype, root (broadcasting rank), comm, ierr
                                call mpi_bcast(imageID, 1, mpi_integer, 0, mpi_comm_world, ierrMPI)
                                ! broadcast co_logFuncState from the winning image to all others: buffer, count, datatype, root (broadcasting rank), comm, ierr
                                call mpi_bcast(co_logFuncState(0 : ndim, -1), ndimp1 * NBRK, mpi_byte, imageID - 1, mpi_comm_world, ierrMPI)
                                !call mpi_bcast(co_logFuncState(0 : ndim, -1), ndimp1, mpi_double_precision, imageID - 1, mpi_comm_world, ierrMPI)
                                stat%avgCommPerFunCall = stat%avgCommPerFunCall + stat%timer%time() - stat%timer%clock
                            end if
#endif
                            ! Note: after every adaptive update of the sampler, counterPAP is reset to 0.
                            if (counterPAP == 0_IK .and. proposalAdaptationSampleUsed) then
                                numFunCallAcceptedLastAdaptation = numFunCallAcceptedLastAdaptation + 1_IK
                                lastStateWeight = 0_IK
                            end if

                            blockFreshDryRun: if (spec%run%is%new) then
                                call writeCFC(spec, stat, proposalAdaptation)
                                stat%numFunCallAccepted = stat%numFunCallAccepted + 1_IK
                                stat%cfc%processID(stat%numFunCallAccepted) = imageID
                                stat%cfc%delayedRejectionStage(stat%numFunCallAccepted) = counterDRS
                                stat%cfc%proposalAdaptation(stat%numFunCallAccepted) = 0._RKC
                                stat%cfc%sampleWeight(stat%numFunCallAccepted) = 0_IK
                                stat%cfc%sampleLogFunc(stat%numFunCallAccepted) = CO_LOGFUNCSTATE(0, imageID, -1)
                                stat%cfc%sampleState(1 : ndim, stat%numFunCallAccepted) = CO_LOGFUNCSTATE(1 : ndim, imageID, -1)
                                ! Find the burnin point.
                                stat%cfc%burninLocation(stat%numFunCallAccepted) = getBurninLoc(stat%numFunCallAccepted, stat%chain%mode%val, stat%cfc%sampleLogFunc(1 : stat%numFunCallAccepted))
                            else blockFreshDryRun ! in restart mode: determine the correct value of co_proposalFound_samplerUpdateOccurred(1).
                                numFunCallAcceptedPlusOne = stat%numFunCallAccepted + 1_IK
                                if (numFunCallAcceptedPlusOne == stat%cfc%nsam) then
                                    spec%run%is%new = .true._LK
                                    call writeCFC(spec, stat, proposalAdaptation)
                                    spec%run%is%dry = .not. spec%run%is%new
                                    stat%cfc%sampleWeight(numFunCallAcceptedPlusOne) = 0_IK
                                    sumAccrAccRej = stat%cfc%meanAcceptanceRate(stat%numFunCallAccepted) * real(stat%numFunCallAcceptedRejected, RKC)
                                    if (proposal%delRejEnabled) sumAccrAccRejDel = stat%cfc%meanAcceptanceRate(stat%numFunCallAccepted) * real(stat%numFunCallAcceptedRejectedDelayed, RKC)
                                end if
                                stat%numFunCallAccepted = numFunCallAcceptedPlusOne
                                numFunCallAcceptedPlusOne = stat%numFunCallAccepted + 1_IK
                            end if blockFreshDryRun
                            if (stat%chain%mode%val < stat%cfc%sampleLogFunc(stat%numFunCallAccepted)) then
                                stat%chain%mode%val = stat%cfc%sampleLogFunc(stat%numFunCallAccepted)
                                stat%chain%mode%loc = stat%numFunCallAccepted
                            end if
                            sumAccrAccRej = sumAccrAccRej + CO_ACCR(imageID, counterDRS)

                        else blockProposalAccepted

                            counterDRS = spec%proposalDelayedRejectionCount%val
                            sumAccrAccRej = sumAccrAccRej + CO_ACCR(1, counterDRS)

                        end if blockProposalAccepted

                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% blockProposalAccepted %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        counterPAP = counterPAP + 1_IK
                        stat%progress%counterPRP = stat%progress%counterPRP + 1_IK
                        ! UPDATE: The following comment is irrelevant and history as of Nov 2023.
                        ! In fact, it is now required to update `stat%numFunCallAcceptedRejected` before calling `reportProgress()`
                        ! because `numFunCallAcceptedRejected` on the sampler status report display will be always off by one.
                        ! OLD COMMENTS: It is critical for this if block to occur before updating `stat%numFunCallAcceptedRejected` otherwise,
                        ! `numFunCallAcceptedRejectedLastReport` will be updated to `stat%numFunCallAcceptedRejected`
                        ! which can lead to "Floating-point exception - erroneous arithmetic operation" in the computation
                        ! of `inverseoutputReportPeriod` when `blockLastSample` happens to be activated.
                        stat%numFunCallAcceptedRejected = stat%numFunCallAcceptedRejected + 1_IK
                        if (stat%progress%counterPRP == spec%outputReportPeriod%val) call reportProgress(spec, stat)
                        currentStateWeight = currentStateWeight + 1_IK
                        if (proposal%delRejEnabled) then
                            sumAccrAccRejDel = sumAccrAccRejDel + sum(CO_ACCR(1, 0 : counterDRS))
                            stat%numFunCallAcceptedRejectedDelayed = stat%numFunCallAcceptedRejectedDelayed + counterDRS + 1_IK
                        end if
                        if (spec%run%is%new) then ! these are used for adaptive proposal updating, so they have to be set on every accepted or rejected iteration (excluding delayed rejections)
                            stat%cfc%meanAcceptanceRate(stat%numFunCallAccepted) = sumAccrAccRej / real(stat%numFunCallAcceptedRejected,kind=RKC)
                            stat%cfc%sampleWeight(stat%numFunCallAccepted) = stat%cfc%sampleWeight(stat%numFunCallAccepted) + 1_IK
                            if (proposalAdaptationLen < stat%cfc%sampleWeight(stat%numFunCallAccepted)) then
                                proposalAdaptationLen = 2_IK * proposalAdaptationLen
                                call setResized(proposalAdaptation, proposalAdaptationLen)
                            end if
                        else
                            SET_CAFMPI(acceptedRejectedDelayedUnusedRestartMode = stat%numFunCallAcceptedRejectedDelayedUnused)
                            SET_CONOMP(acceptedRejectedDelayedUnusedRestartMode = stat%numFunCallAcceptedRejectedDelayedUnused)
                            SET_SERIAL(acceptedRejectedDelayedUnusedRestartMode = stat%numFunCallAcceptedRejectedDelayed)
                        end if

                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% last output write %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        ! in paradise mode, it is imperative to finish the simulation before any further redundant sampler updates occurs.
                        ! This is the reason why blockLastSample appears before blockSamplerAdaptation.

                        blockLastSample: if (stat%numFunCallAccepted == spec%outputChainSize%val) then !co_missionAccomplished = .true._LK
                            ! on 3 images Windows, substituting co_missionAccomplished with the following leads to 10% less communication overhead for 1D Gaussian example
                            ! on 3 images Linux  , substituting co_missionAccomplished with the following leads to 16% less communication overhead for 1D Gaussian example
                            ! on 5 images Linux  , substituting co_missionAccomplished with the following leads to 11% less communication overhead for 1D Gaussian example
                            co_proposalFound_samplerUpdateOccurred(1) = -1 ! equivalent to co_missionAccomplished = .true._LK
                            if (spec%run%is%new) then
                                call writeCFC(spec, stat, proposalAdaptation)
                                flush(spec%chainFile%unit)
                            end if
                            call reportProgress(spec, stat, timeLeft = 0._RKD) ! timeElapsed = stat%timer%time() - stat%timer%start, timeLeft = 0._RKC
                            SET_PARALLEL(exit loopOverImages)
                        end if blockLastSample

                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% last output write %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% proposalAdaptation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        blockSamplerAdaptation: if (counterPAC < spec%proposalAdaptationCount%val .and. counterPAP == spec%proposalAdaptationPeriod%val) then
                            co_proposalFound_samplerUpdateOccurred(2) = 1 ! istart = numFunCallAcceptedLastAdaptation ! = max( numFunCallAcceptedLastAdaptation , stat%cfc%burninLocation(stat%numFunCallAccepted) ) ! this is experimental
                            ! the order in the following two MUST be preserved because occasionally stat%numFunCallAccepted = numFunCallAcceptedLastAdaptation
                            dumint = stat%cfc%sampleWeight(stat%numFunCallAccepted) ! needed for the restart mode, not needed in the fresh run
                            if (stat%numFunCallAccepted == numFunCallAcceptedLastAdaptation) then ! no new point has been accepted since last time
                                stat%cfc%sampleWeight(numFunCallAcceptedLastAdaptation) = currentStateWeight - lastStateWeight
                                CHECK_ASSERTION(__LINE__, mod(stat%cfc%sampleWeight(numFunCallAcceptedLastAdaptation), spec%proposalAdaptationPeriod%val) == 0_IK, PROCEDURE_NAME//SK_": Internal error occurred: "//getStr([spec%proposalAdaptationPeriod%val, stat%cfc%sampleWeight(numFunCallAcceptedLastAdaptation), currentStateWeight, lastStateWeight]))
                            else
                                stat%cfc%sampleWeight(numFunCallAcceptedLastAdaptation) = stat%cfc%sampleWeight(numFunCallAcceptedLastAdaptation) - lastStateWeight
                                stat%cfc%sampleWeight(stat%numFunCallAccepted) = currentStateWeight ! needed for the restart mode, not needed in the fresh run
                            end if
                            meanAccRateSinceStart = stat%cfc%meanAcceptanceRate(stat%numFunCallAccepted) ! used only in fresh runs, but is not worth putting it in a conditional fresh run block.

                            call setProposalAdapted ( spec, proposal &
                                                    , sampleState = stat%cfc%sampleState(:, numFunCallAcceptedLastAdaptation : stat%numFunCallAccepted) &
                                                    , sampleWeight = stat%cfc%sampleWeight(numFunCallAcceptedLastAdaptation : stat%numFunCallAccepted) &
                                                    , proposalAdaptationGreedyEnabled = proposalAdaptationGreedyEnabled &
                                                    , meanAccRateSinceStart = meanAccRateSinceStart &
                                                    , proposalAdaptationSampleUsed = proposalAdaptationSampleUsed &
                                                    , proposalAdaptation = proposalAdaptation(dumint) &
                                                    , err = err &
                                                    )
                            RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(err%occurred),err%msg)
                            if (spec%run%is%dry) sumAccrAccRej = meanAccRateSinceStart * stat%numFunCallAcceptedRejected
                            stat%cfc%sampleWeight(stat%numFunCallAccepted) = dumint   ! needed for the restart mode, not needed in the fresh run, but is not worth fencing it.
                            if (stat%numFunCallAccepted == numFunCallAcceptedLastAdaptation) then
                                !proposalAdaptation = proposalAdaptation + proposalAdaptationDummy ! this is the worst-case upper-bound
                                stat%cfc%proposalAdaptation(stat%numFunCallAccepted) = min(1._RKC, stat%cfc%proposalAdaptation(stat%numFunCallAccepted) + proposalAdaptation(dumint)) ! this is the worst-case upper-bound
                            else
                                !proposalAdaptation = proposalAdaptationDummy
                                stat%cfc%proposalAdaptation(stat%numFunCallAccepted) = proposalAdaptation(dumint)
                                stat%cfc%sampleWeight(numFunCallAcceptedLastAdaptation) = stat%cfc%sampleWeight(numFunCallAcceptedLastAdaptation) + lastStateWeight
                            end if
                            if (proposalAdaptationSampleUsed) then
                                lastStateWeight = currentStateWeight ! stat%cfc%sampleWeight(stat%numFunCallAccepted) ! informative, do not remove
                                numFunCallAcceptedLastAdaptation = stat%numFunCallAccepted
                            end if
                            counterPAP = 0_IK
                            counterPAC = counterPAC + 1_IK
                            !if (counterPAC==spec%proposalAdaptationCount%val) proposalAdaptation = 0._RKC
                        else blockSamplerAdaptation
                            proposalAdaptation(stat%cfc%sampleWeight(stat%numFunCallAccepted)) = 0._RKC
                        end if blockSamplerAdaptation

                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% proposal proposalAdaptation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        SET_PARALLEL(if (co_proposalFound_samplerUpdateOccurred(1) == 1) exit loopOverImages)
                        SET_PARALLEL(imageID = imageID + 1)

                    SET_PARALLEL(end do loopOverImages)

                    SET_CAF(if (spec%parallelism%is%singleChain) then)
                        SET_CAF(if (co_proposalFound_samplerUpdateOccurred(2) == 1) call bcastProposalAdaptation(spec, proposal))
                        SET_CAF(stat%timer%clock = stat%timer%time())
                        SET_CAF(sync images(*))
                        SET_CAF(stat%avgCommPerFunCall = stat%avgCommPerFunCall + stat%timer%time() - stat%timer%clock)
                    SET_CAF(end if)
                    SET_MPI(if (spec%parallelism%is%singleChain .and. co_proposalFound_samplerUpdateOccurred(1) == 0) then)
                        ! broadcast rank #0 to all processes indicating unsuccessful sampling.
                        SET_MPI(imageID = 0)
                        ! buffer; count; datatype; root: broadcasting rank; comm; ierr
                        SET_MPI(call mpi_bcast(imageID, 1, mpi_integer, 0, mpi_comm_world, ierrMPI))
                    SET_MPI(end if)

                SET_CAFMPI(else blockLeaderImage) ! ATTN: This block should be executed only when singleChain parallelism is requested

#if                 CAF_ENABLED
                    sync images(1)
                    ! get the accepted proposal from the first image
                    stat%timer%clock = stat%timer%time()
                    co_proposalFound_samplerUpdateOccurred(1:2) = co_proposalFound_samplerUpdateOccurred(1:2)[1]
                    if (co_proposalFound_samplerUpdateOccurred(1) == 1) co_logFuncState(0 : ndim, -1) = co_logFuncState(0 : ndim, -1)[1]
                    if (co_proposalFound_samplerUpdateOccurred(2) == 1) call bcastProposalAdaptation(spec, proposal)
                    stat%avgCommPerFunCall = stat%avgCommPerFunCall + stat%timer%time() - stat%timer%clock
                    if (spec%run%is%dry) then
                        if (co_proposalFound_samplerUpdateOccurred(1) == 1) then
                            stat%numFunCallAccepted = stat%numFunCallAccepted + 1_IK
                            numFunCallAcceptedPlusOne = stat%numFunCallAccepted + 1_IK
                            currentStateWeight = 1_IK
                            if (stat%numFunCallAccepted == stat%cfc%nsam) then
                                spec%run%is%dry = .false._LK
                                spec%run%is%new = .true._LK
                            end if
                        else
                            currentStateWeight = currentStateWeight + spec%image%count
                        end if
                    else
                        SET_DEBUG(stat%numFunCallAccepted = stat%numFunCallAccepted + co_proposalFound_samplerUpdateOccurred(1))
                    end if
#elif               MPI_ENABLED
                    ! fetch the winning rank from the main process
                    stat%timer%clock = stat%timer%time()
                    ! buffer, count, datatype, root: broadcasting rank, comm, ierr
                    call mpi_bcast(imageID, 1, mpi_integer, 0, mpi_comm_world, ierrMPI)
                    stat%avgCommPerFunCall = stat%avgCommPerFunCall + stat%timer%time() - stat%timer%clock
                    if (0 < imageID) then ! co_proposalFound = .true._LK, sampling successful
                        ! broadcast co_logFuncState from the winning image to all others
                        stat%timer%clock = stat%timer%time()
                        ! buffer, count, datatype, root: broadcasting rank, comm, ierr
                        !call mpi_bcast(co_logFuncState(0 : ndim, -1), ndimp1, mpi_double_precision, imageID - 1, mpi_comm_world, ierrMPI)
                        call mpi_bcast(co_logFuncState(0 : ndim, -1), ndimp1 * NBRK, mpi_byte, imageID - 1, mpi_comm_world, ierrMPI)
                        stat%avgCommPerFunCall = stat%avgCommPerFunCall + stat%timer%time() - stat%timer%clock
                    end if
                    if (spec%run%is%dry) then
                        if (0 < imageID) then ! equivalent to co_proposalFound_samplerUpdateOccurred(1) == 1_IK
                            stat%numFunCallAccepted = stat%numFunCallAccepted + 1_IK
                            numFunCallAcceptedPlusOne = stat%numFunCallAccepted + 1_IK
                            currentStateWeight = 1_IK
                            if (stat%numFunCallAccepted == stat%cfc%nsam) then
                                spec%run%is%new = .true._LK
                                spec%run%is%dry = .false._LK
                            end if
                        else
                            currentStateWeight = currentStateWeight + spec%image%count
                        end if
                    else
                        SET_DEBUG(stat%numFunCallAccepted = stat%numFunCallAccepted + co_proposalFound_samplerUpdateOccurred(1))
                    end if
#endif
                SET_CAFMPI(end if blockLeaderImage)

                RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(err%occurred),err%msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% begin Common block between all images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if             MPI_ENABLED
                ! broadcast samplerUpdateOccurred from the root process to all others, and broadcast proposal adaptations if needed.
                if (spec%parallelism%is%singleChain) then
                    stat%timer%clock = stat%timer%time()
                    ! buffer: XXX: first element of co_proposalFound_samplerUpdateOccurred not needed except to end the simulation.
                    ! This could perhaps be enhanced in the further to only pass one elemtn.
                    ! buffer, count, datatype, root: broadcasting rank, comm, ierr
                    call mpi_bcast(co_proposalFound_samplerUpdateOccurred, 2, mpi_integer, 0, mpi_comm_world, ierrMPI)
                    if (co_proposalFound_samplerUpdateOccurred(2) == 1) call bcastProposalAdaptation(spec, proposal)
                    stat%avgCommPerFunCall = stat%avgCommPerFunCall + stat%timer%time() - stat%timer%clock
                end if
#endif
                if (co_proposalFound_samplerUpdateOccurred(1) == -1) exit loopMarkovChain   ! we are done: co_missionAccomplished = .true._LK
                CO_ACCR(:, -1) = -1._RKC ! counterDRS at which new proposal is accepted. This null initialization is essential for all serial and parallel modes.
                maxLogFuncRejectedProposal = -huge(maxLogFuncRejectedProposal) * .1_RKC

                ! This CAFMPI_SINGLCHAIN_ENABLED preprocessing avoids unnecessary communication in serial or multichain-parallel modes.

#if             CAF_ENABLED || MPI_ENABLED
                blockSingleChainParallelism: if (spec%parallelism%is%singleChain) then
                    proposalFoundSinglChainMode = 0
                    ! This CAF syncing is necessary since the co_logFuncState has to be fetched from the
                    ! first image by all other images before it is updated again below by the first image.
                    SET_CAF(if (spec%image%is%leader) then)
                        SET_CAF(stat%timer%clock = stat%timer%time())
                        SET_CAF(sync images(*))
                        SET_CAF(stat%avgCommPerFunCall = stat%avgCommPerFunCall + stat%timer%time() - stat%timer%clock)
                    SET_CAF(else)
                        SET_CAF(sync images(1))
                    SET_CAF(end if)
#define             CAFMPI_SINGLCHAIN_ENABLED 1
#include            "pm_sampling_kernel.inc.inc.F90"
#undef              CAFMPI_SINGLCHAIN_ENABLED
                    ! gather all accr on the leader image.
                    SET_MPI(stat%timer%clock = stat%timer%time())
                    ! send buffer, send count, send datatype, recv buffer, recv count, recv datatype, root rank, comm, mpi error flag
                    !SET_MPI(call mpi_gather(co_accr, proposalDelayedRejectionCountPlusTwo, mpi_double_precision, accrMat(:,:), proposalDelayedRejectionCountPlusTwo, mpi_double_precision, 0, mpi_comm_world, ierrMPI))
                    SET_MPI(call mpi_gather(co_accr, proposalDelayedRejectionCountPlusTwo, mpi_byte, accrMat(:,:), proposalDelayedRejectionCountPlusTwo, mpi_byte, 0, mpi_comm_world, ierrMPI))
                    SET_MPI(stat%avgCommPerFunCall = stat%avgCommPerFunCall + stat%timer%time() - stat%timer%clock)
                    SET_CAF(if (spec%image%is%leader) then)
                        SET_CAF(stat%timer%clock = stat%timer%time())
                        SET_CAF(sync images(*))
                        SET_CAF(stat%avgCommPerFunCall = stat%avgCommPerFunCall + stat%timer%time() - stat%timer%clock)
                    SET_CAF(else)
                        SET_CAF(sync images(1))
                    SET_CAF(end if)
                else blockSingleChainParallelism ! multichain.
                    block
#include            "pm_sampling_kernel.inc.inc.F90"
                    end block
                end if blockSingleChainParallelism
#elif           OMP_ENABLED
#include        "pm_sampling_kernel.inc.inc.F90"
#else
                ! serial.
#include        "pm_sampling_kernel.inc.inc.F90"
#endif
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end Common block between all images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end do loopMarkovChain

#if         (MPI_ENABLED || CAF_ENABLED) && (CODECOV_ENABLED || SAMPLER_TEST_ENABLED)
            block; use pm_err, only: isFailedImage; call isFailedImage(err%occurred); end block
#endif
            if (err%occurred) return

            stat%timer%clock = stat%timer%time() - stat%timer%start ! This one last timing is crucial for postprocessing report.

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%* end of loopMarkovChain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            stat%cfc%nsam = stat%numFunCallAccepted
            stat%cfc%sumw = stat%numFunCallAcceptedRejected

            if (.not. proposal%delRejEnabled) then
                stat%numFunCallAcceptedRejectedDelayed = stat%numFunCallAcceptedRejected
                sumAccrAccRejDel = sumAccrAccRej
            end if

            if (spec%image%is%leader) then
                block
                    integer(IK) :: i
                    if (0.999999_RKC < spec%proposalAdaptationBurnin%val) then
                        stat%burninLocDRAM%compact = 1_IK
                        stat%burninLocDRAM%verbose = 1_IK
                    else
                        stat%burninLocDRAM%compact = stat%numFunCallAccepted
                        stat%burninLocDRAM%verbose = stat%numFunCallAcceptedRejected - stat%cfc%sampleWeight(stat%numFunCallAccepted) + 1_IK
                        loopAdaptationBurnin: do i = stat%numFunCallAccepted - 1, 1, -1
                            if (spec%proposalAdaptationBurnin%val < stat%cfc%proposalAdaptation(i)) exit loopAdaptationBurnin
                            stat%burninLocDRAM%compact = stat%burninLocDRAM%compact - 1_IK
                            stat%burninLocDRAM%verbose = stat%burninLocDRAM%verbose - stat%cfc%sampleWeight(i)
                        end do loopAdaptationBurnin
                    end if
                end block
                stat%burninLocMCMC%compact = stat%cfc%burninLocation(stat%numFunCallAccepted)
                stat%burninLocMCMC%verbose = sum(stat%cfc%sampleWeight(1 : stat%burninLocMCMC%compact - 1)) + 1_IK
                stat%chain%mode%crd = stat%cfc%sampleState(1 : ndim, stat%chain%mode%loc)
                !stat%chain%mode%Loc%verbose = sum(stat%cfc%sampleWeight(1:stat%chain%mode%loc-1)) + 1_IK
                close(spec%restartFile%unit)
                close(spec%chainFile%unit)
            endif

            ! LCOV_EXCL_START
            ! multichain parallelism should never happen for serial, concurrent, or openmp applications.
            SET_PARALLEL(if (spec%parallelism%is%multiChain) then)
                stat%avgCommPerFunCall = 0._RKC
                stat%numFunCallAcceptedRejectedDelayedUnused = stat%numFunCallAcceptedRejectedDelayed
                dumint = stat%numFunCallAcceptedRejectedDelayedUnused - acceptedRejectedDelayedUnusedRestartMode ! this is needed to avoid division-by-zero undefined behavior
                if (dumint /= 0_IK) stat%avgTimePerFunCall = stat%avgTimePerFunCall / dumint
            SET_PARALLEL(elseif(spec%image%is%first) then)
                SET_PARALLEL(stat%avgCommPerFunCall = stat%avgCommPerFunCall / stat%numFunCallAcceptedRejectedDelayed)
                SET_PARALLEL(dumint = stat%numFunCallAcceptedRejectedDelayedUnused - acceptedRejectedDelayedUnusedRestartMode) ! this is needed to avoid division-by-zero undefined behavior
                SET_PARALLEL(if (dumint /= 0_IK) stat%avgTimePerFunCall = (stat%avgTimePerFunCall / dumint) * spec%image%count)
            SET_PARALLEL(end if)
            ! LCOV_EXCL_STOP
#undef      CO_LOGFUNCSTATE
#undef      SET_PARALLEL
#undef      SET_CONOMP
#undef      SHAPE_ARG
#undef      GET_ELL
#undef      CO_ACCR