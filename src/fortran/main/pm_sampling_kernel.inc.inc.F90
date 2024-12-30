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

                    loopNextMove: do counterDRS = 0, spec%proposalDelayedRejectionCount%val
#if                     OMP_ENABLED || CAFMPI_SINGLCHAIN_ENABLED
                        ! stat%numFunCallAcceptedRejectedDelayedUnused is relevant only on the first image, despite being updated by all images
                        stat%numFunCallAcceptedRejectedDelayedUnused = stat%numFunCallAcceptedRejectedDelayedUnused + spec%image%count
#endif
#if                     OMP_ENABLED
                        do imageID = 1, spec%image%count
                            call setProposalStateNew(spec, proposal, counterDRS, co_logFuncState(1 : ndim, imageID, counterDRS - 1), co_logFuncState(1 : ndim, imageID, counterDRS), spec%rng, err)
                            if(err%occurred) then; err%msg = getFine(__FILE__, __LINE__)//PROCEDURE_NAME//SK_": "//err%msg; return; end if
                        end do
#else
                        call setProposalStateNew(spec, proposal, counterDRS, co_logFuncState(1 : ndim, counterDRS - 1), co_logFuncState(1 : ndim, counterDRS), spec%rng, err)
                        if(err%occurred) then; err%msg = getFine(__FILE__, __LINE__)//PROCEDURE_NAME//SK_": "//err%msg; return; end if
#endif
!!#if                    (MATLAB_ENABLED || PYTHON_ENABLED || R_ENABLED) && !CAF_ENABLED && !MPI_ENABLED
!!                       if(err%occurred) then; err%msg = PROCEDURE_NAME//SK_": "//err%msg; return; end if
!!#elif                  CODECOV_ENABLED || SAMPLER_TEST_ENABLED
!#if                     CODECOV_ENABLED || SAMPLER_TEST_ENABLED
!                        ! This block is exclusively used to test the deterministic restart functionality of ParaXXXX samplers.
!                        ! This block must not be activated under any other circumstances.
!                        ! This block must be executed by all images.
!                        !spec%testSamplingCounter = spec%testSamplingCounter + 1_IK
!                        !err%occurred = err%occurred .or. spec%testSamplingCountTarget < spec%testSamplingCounter
!                        !if (err%occurred) then
!                        !    if (spec%testSamplingCounter > spec%testSamplingCountTarget) err%msg = "The simulation was interrupted at the requested sampling count for restart testing purposes."
!                        !    SET_CAFMPI(call isFailedImage(err%occurred))
!                        !    return
!                        !end if
!#endif
                        ! The following random function call is only needed in fresh runs to evaluate the acceptance of a proposal.
                        ! However, it is taken out of the subsequent loops to achieve 100% deterministic reproducibility when
                        ! a simulation is restarted.
                        call setUnifRand(spec%rng, unifrnd)
                        if (spec%run%is%new .or. numFunCallAcceptedPlusOne == stat%cfc%nsam) then
#if                         OMP_ENABLED && (MATLAB_ENABLED || PYTHON_ENABLED || R_ENABLED)
                            block
                                real(RKD) :: avgTimePerFunCall, avgCommPerFunCall
                                !avgTimePerFunCall = epsilon(0._RKD); avgCommPerFunCall = epsilon(0._RKD)
                                !block
                                !use pm_io
                                !call disp%skip(count = 2)
                                !call disp%show("before")
                                !call disp%show("co_logFuncState(:, 1 : spec%image%count, counterDRS)")
                                !call disp%show( co_logFuncState(:, 1 : spec%image%count, counterDRS) )
                                !call disp%skip(count = 2)
                                !end block
                                !stat%timer%clock = stat%timer%time()
                                mold = getLogFunc(co_logFuncState(:, 1 : spec%image%count, counterDRS), avgTimePerFunCall, avgCommPerFunCall)
                                !avgTimePerFunCall = max(stat%timer%resol * 10, stat%timer%time() - stat%timer%clock - avgCommPerFunCall)
                                !block
                                !use pm_io
                                !call disp%skip(count = 2)
                                !call disp%show("after")
                                !call disp%show("co_logFuncState(:, 1 : spec%image%count, counterDRS)")
                                !call disp%show( co_logFuncState(:, 1 : spec%image%count, counterDRS) )
                                !call disp%skip(count = 2)
                                !end block
                                stat%avgTimePerFunCall = stat%avgTimePerFunCall + avgTimePerFunCall
                                stat%avgCommPerFunCall = stat%avgCommPerFunCall + avgCommPerFunCall
                            end block
#elif                       OMP_ENABLED
                            block
                                use pm_timer, only: getTime => getTimeOMP
                                real(RKD) :: start, avgTimePerFunCall
                                avgTimePerFunCall = 0._RKD
                                stat%timer%clock = stat%timer%time()
                                ! preliminary tests indicate that false sharing is not a concern here.
                                ! nevertheless, the current implementation attempts to avoid it by assuming a 128-byte cacheline.
                               !!$omp parallel do reduction(+ : avgTimePerFunCall) num_threads(spec%image%count) default(none) shared(co_logFuncState, counterDRS, ndim, spec) private(imageID, start, logFuncState)
                                !$omp parallel do reduction(+ : avgTimePerFunCall) num_threads(spec%image%count) default(none) shared(co_logFuncState, counterDRS, ndim, spec, logFuncState) private(imageID, start)
                                do imageID = 1, spec%image%count
                                    start = getTime()
                                    logFuncState(1, imageID) = getLogFunc(co_logFuncState(1 : ndim, imageID, counterDRS))
                                    !avgTimePerFunCall(imageID * JUMP) = getTime() - start
                                    !co_logFuncState(0, imageID, counterDRS) = getLogFunc(co_logFuncState(1 : ndim, imageID, counterDRS))
                                    avgTimePerFunCall = avgTimePerFunCall + (getTime() - start)
                                end do
                                !$omp end parallel do
                                co_logFuncState(0, :, counterDRS) = logFuncState(1, :)
                                avgTimePerFunCall = avgTimePerFunCall / spec%image%count
                                stat%avgTimePerFunCall = stat%avgTimePerFunCall + avgTimePerFunCall
                                stat%avgCommPerFunCall = stat%avgCommPerFunCall + stat%timer%time() - stat%timer%clock - avgTimePerFunCall
                            end block
#else
                            stat%timer%clock = stat%timer%time()
                            co_logFuncState(0, counterDRS) = getLogFunc(co_logFuncState(1 : ndim, counterDRS))
                            stat%avgTimePerFunCall = stat%avgTimePerFunCall + stat%timer%time() - stat%timer%clock
#endif
                            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            ! accept or reject the proposed state
                            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if                         OMP_ENABLED
                            do imageID = 1, spec%image%count
#endif
                                if (CO_LOGFUNCSTATE(0, imageID, -1) <= CO_LOGFUNCSTATE(0, imageID, counterDRS)) then ! accept the proposed state.
                                    CO_ACCR(imageID, counterDRS) = 1._RKG
                                    CO_ACCR(imageID, -1) = real(counterDRS, RKG)
                                    CO_LOGFUNCSTATE(0 : ndim, imageID, -1) = CO_LOGFUNCSTATE(0 : ndim, imageID, counterDRS)
#if                                 OMP_ENABLED || !CAFMPI_SINGLCHAIN_ENABLED
                                    exit loopNextMove
#endif
                                elseif (CO_LOGFUNCSTATE(0, imageID, counterDRS) < GET_ELL(maxLogFuncRejectedProposal, imageID)) then ! reject the proposed state. This step should be reachable only when proposalDelayedRejectionCount > 0
                                    CO_ACCR(imageID, counterDRS) = 0._RKG ! proposal rejected XXX is this correct? could `co_accr` be absolute zero?
                                else ! accept with probability co_accr.
                                    if (counterDRS == 0_IK) then ! This should be equivalent to GET_ELL(maxLogFuncRejectedProposal, imageID) == NEGBIG_RK
                                        logFuncDiff = CO_LOGFUNCSTATE(0, imageID, counterDRS) - CO_LOGFUNCSTATE(0, imageID, -1)
                                        if (logFuncDiff < log(tiny(0._RKG))) then ! xxx should the condition for LOGHUGE_RK be also added?
                                            CO_ACCR(imageID, counterDRS) = 0._RKG
                                        else
                                            CO_ACCR(imageID, counterDRS) = exp(logFuncDiff)
                                        end if
                                    else    ! ensure no arithmetic overflow/underflow. ATTN: GET_ELL(maxLogFuncRejectedProposal, imageID) < co_logFuncState(0, counterDRS) < co_logFuncState(0, -1)
                                        CO_ACCR(imageID, counterDRS) & ! LCOV_EXCL_LINE
                                        = getLogSubExp(smaller = GET_ELL(maxLogFuncRejectedProposal, imageID), larger = CO_LOGFUNCSTATE(0, imageID, counterDRS)) & ! LCOV_EXCL_LINE
                                        - getLogSubExp(smaller = GET_ELL(maxLogFuncRejectedProposal, imageID), larger = CO_LOGFUNCSTATE(0, imageID, -1))
                                        if (CO_ACCR(imageID, counterDRS) < log(tiny(0._RKG))) then
                                            CO_ACCR(imageID, counterDRS) = 0._RKG
                                        else
                                            CO_ACCR(imageID, counterDRS) = exp(CO_ACCR(imageID, counterDRS))
                                        end if
                                    end if
                                    if (GET_ELL(unifrnd, imageID) < CO_ACCR(imageID, counterDRS)) then ! accept the proposed state
                                        CO_LOGFUNCSTATE(0 : ndim, imageID, -1) = CO_LOGFUNCSTATE(0 : ndim, imageID, counterDRS)
                                        CO_ACCR(imageID, -1) = real(counterDRS, RKG)
#if                                     OMP_ENABLED || !CAFMPI_SINGLCHAIN_ENABLED
                                        exit loopNextMove
#endif
                                    end if
                                end if
                                GET_ELL(maxLogFuncRejectedProposal, imageID) = max(GET_ELL(maxLogFuncRejectedProposal, imageID), CO_LOGFUNCSTATE(0, imageID, counterDRS))
#if                         OMP_ENABLED
                            end do
#endif
                        else ! if dryrun
#if                         CAFMPI_SINGLCHAIN_ENABLED
                            if (spec%image%id == stat%cfc%processID(numFunCallAcceptedPlusOne) .and. &
                                currentStateWeight + spec%image%id - 1_IK == stat%cfc%sampleWeight(stat%numFunCallAccepted) .and. &
                                counterDRS == stat%cfc%delayedRejectionStage(numFunCallAcceptedPlusOne)) then
                                co_accr(-1) = real(counterDRS, RKG)
                                co_logFuncState(0, counterDRS) = stat%cfc%sampleLogFunc(numFunCallAcceptedPlusOne)
                                co_logFuncState(1 : ndim, counterDRS) = stat%cfc%sampleState(1 : ndim,numFunCallAcceptedPlusOne)
                                co_logFuncState(0 : ndim, -1) = co_logFuncState(0 : ndim, counterDRS)
                            end if
#elif                       OMP_ENABLED
                            imageID = stat%cfc%processID(numFunCallAcceptedPlusOne)
                            if (currentStateWeight + imageID - 1_IK == stat%cfc%sampleWeight(stat%numFunCallAccepted) .and. counterDRS == stat%cfc%delayedRejectionStage(numFunCallAcceptedPlusOne)) then
                                co_accr(imageID, -1) = real(counterDRS, RKG)
                                co_logFuncState(0, imageID, counterDRS) = stat%cfc%sampleLogFunc(numFunCallAcceptedPlusOne)
                                co_logFuncState(1 : ndim, imageID, counterDRS) = stat%cfc%sampleState(1 : ndim,numFunCallAcceptedPlusOne)
                                co_logFuncState(0 : ndim, imageID, -1) = co_logFuncState(0 : ndim, imageID, counterDRS)
                                exit loopNextMove
                            end if
#else
                            if (currentStateWeight == stat%cfc%sampleWeight(stat%numFunCallAccepted) .and. counterDRS == stat%cfc%delayedRejectionStage(numFunCallAcceptedPlusOne)) then
                                co_accr(-1) = real(counterDRS, RKG)
                                co_logFuncState(0, -1) = stat%cfc%sampleLogFunc(numFunCallAcceptedPlusOne)
                                co_logFuncState(1 : ndim, -1) = stat%cfc%sampleState(1 : ndim, numFunCallAcceptedPlusOne)
                                exit loopNextMove
                            end if
#endif
                        end if
#if                     CAFMPI_SINGLCHAIN_ENABLED
                        ! broadcast the sampling status from the all images to all images.
                        ! This is needed in singleChain parallelism to avoid unnecessary delayed rejection if a proposal is already found.
                        ! Note that this is not necessary when the delayed rejection is deactivated.
                        if (0 < spec%proposalDelayedRejectionCount%val) then
                            if (-0.5_RKG < co_accr(-1)) proposalFoundSinglChainMode = 1 ! -0.5 instead of -1. avoids possible real roundoff errors.
                            stat%timer%clock = stat%timer%time()
                            SET_CAF(call co_sum(proposalFoundSinglChainMode))
                            ! send buffer, recv buffer, buffer size, datatype, mpi reduction operation, comm, ierr
                            SET_MPI(call mpi_allreduce(proposalFoundSinglChainMode, proposalFoundSinglChainModeReduced, 1, mpi_integer, mpi_sum, mpi_comm_world, ierrMPI))
                            SET_MPI(proposalFoundSinglChainMode = proposalFoundSinglChainModeReduced)
                            stat%avgCommPerFunCall = stat%avgCommPerFunCall + stat%timer%time() - stat%timer%clock
                            if (0_IK < proposalFoundSinglChainMode) exit loopNextMove
                        end if
#endif
                    end do loopNextMove