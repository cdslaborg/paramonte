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
!!!!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined SINGLCHAIN_PARALLELISM
#define LOOP_NEXT_MOVE loopNextMoveSinglChain
#else
#define LOOP_NEXT_MOVE loopNextMove
#endif

                LOOP_NEXT_MOVE : do counterDRS = 0, self%SpecDRAM%DelayedRejectionCount%val

#if defined SINGLCHAIN_PARALLELISM
                    ! self%Stats%NumFunCall%acceptedRejectedDelayedUnused is relevant only on the first image, despite being updated by all images
                    self%Stats%NumFunCall%acceptedRejectedDelayedUnused = self%Stats%NumFunCall%acceptedRejectedDelayedUnused + self%Image%count
#endif

                    co_LogFuncState(1:nd,counterDRS) = self%Proposal%getNew ( nd            = nd &
                                                                            , counterDRS    = counterDRS &
                                                                            , StateOld      = co_LogFuncState(1:nd,counterDRS-1) &
                                                                            )
#if (defined MATLAB_ENABLED || defined PYTHON_ENABLED) && !defined CAF_ENABLED && !defined MPI_ENABLED
                    if(ProposalErr%occurred) then; self%Err%occurred = .true.; return; end if
#endif

                    call random_number(uniformRnd) ! only for the purpose of restart mode reproducibility

#if defined SINGLCHAIN_PARALLELISM && defined CAF_ENABLED
                    ! this is necessary to avoid racing condition on co_LogFuncState and co_proposalFoundSinglChainMode
                    if (self%Image%isMaster) then
                        call self%Timer%toc()
                        sync images(*)
                        call self%Timer%toc(); self%Stats%avgCommTimePerFunCall = self%Stats%avgCommTimePerFunCall + self%Timer%Time%delta
                    else
                        sync images(1)
                    end if
#endif

                    if (self%isFreshRun .or. numFunCallAcceptedPlusOne==self%Chain%Count%compact) then

#if defined SINGLCHAIN_PARALLELISM
                        if (co_AccRate(-1)<-0.5_RK) then ! accept or reject the proposed state, only if no acceptance has occurred yet
#endif

                            call self%Timer%toc()
                            co_LogFuncState(0,counterDRS) = getLogFunc(nd,co_LogFuncState(1:nd,counterDRS))
                            call self%Timer%toc(); self%Stats%avgTimePerFunCalInSec = self%Stats%avgTimePerFunCalInSec + self%Timer%Time%delta

                            ! accept or reject the proposed state

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
                                    co_AccRate(counterDRS) = exp( co_LogFuncState(0,counterDRS) - co_LogFuncState(0,-1) )
                                else    ! ensure no arithmetic overflow/underflow. ATT: co_LogFuncState(0,-1) > co_LogFuncState(0,counterDRS) > maxLogFuncRejectedProposal
                                    co_AccRate(counterDRS) = exp( getLogSubExp( co_LogFuncState(0,counterDRS)   , maxLogFuncRejectedProposal ) &
                                                                - getLogSubExp( co_LogFuncState(0,-1)           , maxLogFuncRejectedProposal ) )
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

#if defined SINGLCHAIN_PARALLELISM
                        end if
#endif

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
                    co_proposalFoundSinglChainMode = co_AccRate(-1) > -1._RK
                    if (delayedRejectionRequested) then ! broadcast the sampling status from the first image to all others
#if defined CAF_ENABLED
                        if (self%Image%isMaster) then ! this is necessary to avoid racing condition on the value of co_proposalFoundSinglChainMode
                            call self%Timer%toc()
                            sync images(*)
                            call self%Timer%toc(); self%Stats%avgCommTimePerFunCall = self%Stats%avgCommTimePerFunCall + self%Timer%Time%delta
                        else
                            sync images(1)
                        end if
                        co_proposalFoundSinglChainMode = co_proposalFoundSinglChainMode[1]
#elif defined MPI_ENABLED
                        ! broadcast winning image to all processes
                        call mpi_bcast  ( co_proposalFoundSinglChainMode    &   ! buffer
                                        , 1                                 &   ! count
                                        , mpi_logical                       &   ! datatype
                                        , 0                                 &   ! root: broadcasting rank
                                        , mpi_comm_world                    &   ! comm
                                        , ierrMPI                           &   ! ierr
                                        )
#endif
                    end if
                    if (co_proposalFoundSinglChainMode) exit LOOP_NEXT_MOVE
#endif

                end do LOOP_NEXT_MOVE

#if defined SINGLCHAIN_PARALLELISM && defined MPI_ENABLED
                ! gather all accr on the master image. Avoid unnecessary communication when proposal is on the master image
                if ( noDelayedRejectionRequested .or. (delayedRejectionRequested .and. .not.co_proposalFoundSinglChainMode) ) then 
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
                end if
#endif

#undef LOOP_NEXT_MOVE

