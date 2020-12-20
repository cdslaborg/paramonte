!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!   MIT License
!!!!
!!!!   ParaMonte: plain powerful parallel Monte Carlo library.
!!!!
!!!!   Copyright (C) 2012-present, The Computational Data Science Lab
!!!!
!!!!   This file is part of the ParaMonte library.
!!!!
!!!!   Permission is hereby granted, free of charge, to any person obtaining a
!!!!   copy of this software and associated documentation files (the "Software"),
!!!!   to deal in the Software without restriction, including without limitation
!!!!   the rights to use, copy, modify, merge, publish, distribute, sublicense,
!!!!   and/or sell copies of the Software, and to permit persons to whom the
!!!!   Software is furnished to do so, subject to the following conditions:
!!!!
!!!!   The above copyright notice and this permission notice shall be
!!!!   included in all copies or substantial portions of the Software.
!!!!
!!!!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!!!!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!!!!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!!!!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
!!!!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
!!!!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
!!!!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!!!!
!!!!   ACKNOWLEDGMENT
!!!!
!!!!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!!!!   As per the ParaMonte library license agreement terms, if you use any parts of
!!!!   this library for any purposes, kindly acknowledge the use of ParaMonte in your
!!!!   work (education/research/industry/development/...) by citing the ParaMonte
!!!!   library as described on this page:
!!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined SINGLCHAIN_PARALLELISM
#define LOOP_NEXT_MOVE loopNextMoveSinglChain
#else
#define LOOP_NEXT_MOVE loopNextMove
#endif

#if defined SINGLCHAIN_PARALLELISM
                proposalFoundSinglChainMode = 0_IK
#if defined CAF_ENABLED
                ! This syncing is necessary since the co_LogFuncState has to be fetched from the 
                ! first image by all other images before it is updated again below by the first image.
                if (self%Image%isLeader) then
                    call self%Timer%toc()
                    sync images(*)
                    call self%Timer%toc(); self%Stats%avgCommTimePerFunCall = self%Stats%avgCommTimePerFunCall + self%Timer%Time%delta
                else
                    sync images(1)
                end if
#endif
#endif

                LOOP_NEXT_MOVE : do counterDRS = 0, self%SpecDRAM%DelayedRejectionCount%val

#if defined SINGLCHAIN_PARALLELISM
                    ! self%Stats%NumFunCall%acceptedRejectedDelayedUnused is relevant only on the first image, despite being updated by all images
                    self%Stats%NumFunCall%acceptedRejectedDelayedUnused = self%Stats%NumFunCall%acceptedRejectedDelayedUnused + self%Image%count
#endif

                    co_LogFuncState(1:nd,counterDRS) = self%Proposal%getNew ( nd            = nd & ! LCOV_EXCL_LINE
                                                                            , counterDRS    = counterDRS & ! LCOV_EXCL_LINE
                                                                            , StateOld      = co_LogFuncState(1:nd,counterDRS-1) & ! LCOV_EXCL_LINE
                                                                            )
#if (defined MATLAB_ENABLED || defined PYTHON_ENABLED || defined R_ENABLED) && !defined CAF_ENABLED && !defined MPI_ENABLED
                    if(ProposalErr%occurred) then; self%Err%occurred = .true.; self%Err%msg = ProposalErr%msg; return; end if
#elif defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
                    ! This block is exclusively used to test the deterministic restart functionality of ParaDXXX samplers. 
                    ! This block must not be activated under any other circumstances.
                    ! This block must be executed by all images.
                    self%testSamplingCounter = self%testSamplingCounter + 1_IK
                    self%Err%occurred = self%Err%occurred .or. self%testSamplingCounter > self%testSamplingCountTarget
                    if (self%Err%occurred) then
                        if (self%testSamplingCounter > self%testSamplingCountTarget) self%Err%msg = "The simulation was interrupted at the requested sampling count for restart testing purposes."
#if defined MPI_ENABLED || defined CAF_ENABLED
                        block; use Err_mod, only: bcastErr; call bcastErr(self%Err); end block
#endif
                        return
                    end if
#endif

                    ! The following random function call is only needed fresh runs to evaluate the acceptance of a proposal.
                    ! However, it is taken out of the subsequent loops to achieve 100% deterministic reproducibility when
                    ! a simulation is restarted.

                    call random_number(uniformRnd)

                    if (self%isFreshRun .or. numFunCallAcceptedPlusOne==self%Chain%Count%compact) then

                        call self%Timer%toc()
                        co_LogFuncState(0,counterDRS) = getLogFunc(nd,co_LogFuncState(1:nd,counterDRS))
                        call self%Timer%toc(); self%Stats%avgTimePerFunCalInSec = self%Stats%avgTimePerFunCalInSec + self%Timer%Time%delta

                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        ! accept or reject the proposed state
                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        if ( co_LogFuncState(0,counterDRS) >= co_LogFuncState(0,-1) ) then ! accept the proposed state

                            co_AccRate(counterDRS) = 1._RK
                            co_AccRate(-1) = real(counterDRS,kind=RK)
                            co_LogFuncState(0:nd,-1) = co_LogFuncState(0:nd,counterDRS)
#if !defined SINGLCHAIN_PARALLELISM
                            exit LOOP_NEXT_MOVE
#endif

                            elseif ( co_LogFuncState(0,counterDRS) < maxLogFuncRejectedProposal ) then  ! reject the proposed state. This step should be reachable only when delayedRejectionCount > 0

                                co_AccRate(counterDRS) = 0._RK  ! proposal rejected XXX is this correct? could co_AccRatebe absolute zero?

                            else    ! accept with probability co_AccRate

                            if ( counterDRS == 0_IK ) then ! This should be equivalent to maxLogFuncRejectedProposal == NEGINF_RK
                                logFuncDiff = co_LogFuncState(0,counterDRS) - co_LogFuncState(0,-1)
                                if (logFuncDiff < NEGLOGINF_RK) then ! xxx should the condition for LOGHUGE_RK be also added?
                                    co_AccRate(counterDRS) = 0._RK
                                else
                                    co_AccRate(counterDRS) = exp( logFuncDiff )
                                end if
                            else    ! ensure no arithmetic overflow/underflow. ATT: co_LogFuncState(0,-1) > co_LogFuncState(0,counterDRS) > maxLogFuncRejectedProposal
                                co_AccRate(counterDRS)  = getLogSubExp( co_LogFuncState(0,counterDRS)   , maxLogFuncRejectedProposal ) & ! LCOV_EXCL_LINE
                                                        - getLogSubExp( co_LogFuncState(0,-1)           , maxLogFuncRejectedProposal )
                                if (co_AccRate(counterDRS) < NEGLOGINF_RK) then
                                    co_AccRate(counterDRS) = 0._RK
                                else
                                    co_AccRate(counterDRS) = exp(co_AccRate(counterDRS))
                                end if
                            end if

                            if (uniformRnd<co_AccRate(counterDRS)) then ! accept the proposed state
                                co_AccRate(-1) = real(counterDRS,kind=RK)
                                co_LogFuncState(0:nd,-1) = co_LogFuncState(0:nd,counterDRS)
#if !defined SINGLCHAIN_PARALLELISM
                                exit LOOP_NEXT_MOVE
#endif
                            end if

                        end if

                        maxLogFuncRejectedProposal = max( maxLogFuncRejectedProposal, co_LogFuncState(0,counterDRS) )

                    else ! if dryrun

#if defined SINGLCHAIN_PARALLELISM
                        if (self%Image%id==self%Chain%ProcessID(numFunCallAcceptedPlusOne) .and. &
                            currentStateWeight+self%Image%id-1_IK==self%Chain%Weight(self%Stats%NumFunCall%accepted) .and. &
                            counterDRS==self%Chain%DelRejStage(numFunCallAcceptedPlusOne)) then
                            co_AccRate(-1) = real(counterDRS,kind=RK)
                            co_LogFuncState(   0,counterDRS) = self%Chain%LogFunc   (numFunCallAcceptedPlusOne)
                            co_LogFuncState(1:nd,counterDRS) = self%Chain%State(1:nd,numFunCallAcceptedPlusOne)
                            co_LogFuncState(0:nd,-1) = co_LogFuncState(0:nd,counterDRS)
                        end if
#else
                        if (currentStateWeight==self%Chain%Weight(self%Stats%NumFunCall%accepted) .and. counterDRS==self%Chain%DelRejStage(numFunCallAcceptedPlusOne)) then
                            co_AccRate(-1) = real(counterDRS,kind=RK)
                            co_LogFuncState(   0,-1) = self%Chain%LogFunc   (numFunCallAcceptedPlusOne)
                            co_LogFuncState(1:nd,-1) = self%Chain%State(1:nd,numFunCallAcceptedPlusOne)
                            exit LOOP_NEXT_MOVE
                        end if
#endif
                    end if

#if defined SINGLCHAIN_PARALLELISM
                    ! broadcast the sampling status from the all images to all images.
                    ! This is needed in singlChain parallelism to avoid unnecessary delayed rejection if a proposal is already found. 
                    ! Note that this is not necessary when the delayed rejection is deactivated.
                    if (delayedRejectionRequested) then
                        if (co_AccRate(-1) > -0.5_RK) proposalFoundSinglChainMode = 1_IK ! -0.5 instead of -1. avoids possible real roundoff errors.
                        call self%Timer%toc()
#if defined CAF_ENABLED
                        call co_sum(proposalFoundSinglChainMode)
#elif defined MPI_ENABLED
                        call mpi_allreduce  ( proposalFoundSinglChainMode           & ! send buffer
                                            , proposalFoundSinglChainModeReduced    & ! recv buffer
                                            , 1                                     & ! buffer size
                                            , mpi_integer                           & ! datatype
                                            , mpi_sum                               & ! mpi reduction operation
                                            , mpi_comm_world                        & ! comm
                                            , ierrMPI                               & ! ierr
                                            )
                        proposalFoundSinglChainMode = proposalFoundSinglChainModeReduced
#endif
                        call self%Timer%toc(); self%Stats%avgCommTimePerFunCall = self%Stats%avgCommTimePerFunCall + self%Timer%Time%delta
                        if (proposalFoundSinglChainMode>0_IK) exit LOOP_NEXT_MOVE
                    end if
#endif

                end do LOOP_NEXT_MOVE

#if defined SINGLCHAIN_PARALLELISM && defined MPI_ENABLED
                ! This SINGLCHAIN_PARALLELISM preprocessing avoids unnecessary communication in serial or multichain-parallel modes.
                ! gather all accr on the leader image. 
                call self%Timer%toc()
                call mpi_gather ( co_AccRate(:)                 & ! send buffer
                                , delayedRejectionCountPlusTwo  & ! send count
                                , mpi_double_precision          & ! send datatype
                                , AccRateMatrix(:,:)            & ! recv buffer
                                , delayedRejectionCountPlusTwo  & ! recv count
                                , mpi_double_precision          & ! recv datatype
                                , 0                             & ! root rank
                                , mpi_comm_world                & ! comm
                                , ierrMPI                       & ! mpi error flag
                                )
                call self%Timer%toc(); self%Stats%avgCommTimePerFunCall = self%Stats%avgCommTimePerFunCall + self%Timer%Time%delta
#elif defined SINGLCHAIN_PARALLELISM && defined CAF_ENABLED
                if (self%Image%isLeader) then
                    call self%Timer%toc()
                    sync images(*)
                    call self%Timer%toc(); self%Stats%avgCommTimePerFunCall = self%Stats%avgCommTimePerFunCall + self%Timer%Time%delta
                else
                    sync images(1)
                end if
#endif

#undef LOOP_NEXT_MOVE
