!***********************************************************************************************************************************
!***********************************************************************************************************************************
!
!   ParaMonte: plain powerful parallel Monte Carlo library.
!
!   Copyright (C) 2012-present, The Computational Data Science Lab
!
!   This file is part of the ParaMonte library.
!
!   ParaMonte is free software: you can redistribute it and/or modify it
!   under the terms of the GNU Lesser General Public License as published
!   by the Free Software Foundation, version 3 of the License.
!
!   ParaMonte is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the ParaMonte library. If not, see,
!
!       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
!
!   ACKNOWLEDGMENT
!
!   As per the ParaMonte library license agreement terms,
!   if you use any parts of this library for any purposes,
!   we ask you to acknowledge the use of the ParaMonte library
!   in your work (education/research/industry/development/...)
!   by citing the ParaMonte library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!***********************************************************************************************************************************
!***********************************************************************************************************************************

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
#if MATLAB_ENABLED && !defined CAF_ENABLED && !defined MPI_ENABLED
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

#if defined ParaDISE
                !if (counterDRS/=0_IK) then
                !    write(*,*) "counterDRS/=0_IK: ", counterDRS
                !    stop
                !end if
                if (maxNumFunCallAcceptedRejected==self%Stats%NumFunCall%acceptedRejected) then

                    maxNumFunCallAcceptedRejectedNew = ceiling(maxNumFunCallAcceptedRejected * (1._RK + getRemainingSimulationFraction()))
!write(*,*) maxNumFunCallAcceptedRejectedNew, 2 * maxNumFunCallAcceptedRejected
!block
!integer :: i
!real(RK), allocatable :: Dummy(:,:)
!write(*,*) 
!write(*,*) size(LogFuncStateDISE)
!write(*,*) lbound(LogFuncStateDISE,1), ubound(LogFuncStateDISE,1)
!write(*,*) lbound(LogFuncStateDISE,2), ubound(LogFuncStateDISE,2)
!if (allocated(Dummy)) deallocate(Dummy); allocate(Dummy, source = LogFuncStateDISE)
                    allocate(LogFuncStateDISE_placeHolder(0:nd,1:maxNumFunCallAcceptedRejectedNew))!, source = 9999999._RK)
                    LogFuncStateDISE_placeHolder(0:nd,1:maxNumFunCallAcceptedRejected) = LogFuncStateDISE(0:nd,1:maxNumFunCallAcceptedRejected)
                    call move_alloc(from=LogFuncStateDISE_placeHolder,to=LogFuncStateDISE)
!do i = 1, maxNumFunCallAcceptedRejected
!    if (any(Dummy(0:nd,i)/=LogFuncStateDISE(0:nd,i))) then
!        write(*,"(*(g0,:,' '))") "Dummy/=LogFuncStateDISE", Dummy(0:nd,i:i)
!        write(*,"(*(g0,:,' '))") "Dummy/=LogFuncStateDISE", LogFuncStateDISE(0:nd,i:i)
!        stop
!    end if
!end do
!write(*,*) size(LogFuncStateDISE)
!write(*,*) lbound(LogFuncStateDISE,1), ubound(LogFuncStateDISE,1)
!write(*,*) lbound(LogFuncStateDISE,2), ubound(LogFuncStateDISE,2)
!write(*,*) 
!write(*,*) size(Dummy), size(LogFuncStateDISE)
!write(*,*) 
!end block

!block
!integer :: i
!real(RK), allocatable :: Dummy(:)
!write(*,*) 
!write(*,*) size(LogProbTrans)
!write(*,*) lbound(LogProbTrans,1), ubound(LogProbTrans,1)
!if (allocated(Dummy)) deallocate(Dummy); Dummy = LogProbTrans(1:maxNumFunCallAcceptedRejected)
                    allocate(LogProbTrans_placeHolder(1:maxNumFunCallAcceptedRejectedNew))!, source = 9999999._RK)
                    LogProbTrans_placeHolder(1:maxNumFunCallAcceptedRejected) = LogProbTrans!(1:maxNumFunCallAcceptedRejected)
                    call move_alloc(from=LogProbTrans_placeHolder,to=LogProbTrans)
!do i = 1, maxNumFunCallAcceptedRejected
!write(*,"(*(g0,:,' '))") "i", i, Dummy(i), LogProbTrans(i)
!    if (Dummy(i)/=LogProbTrans(i)) then
!        write(*,"(*(g0,:,' '))") "Dummy/=LogProbTrans", Dummy(i:i)
!        write(*,"(*(g0,:,' '))") "Dummy/=LogProbTrans", LogProbTrans(i:i)
!        stop
!    end if
!end do
!write(*,*) size(LogProbTrans)
!write(*,*) lbound(LogProbTrans,1), ubound(LogProbTrans,1)
!write(*,*) 
!end block

!block
!integer :: i
!real(RK), allocatable :: Dummy(:)
!if (allocated(Dummy)) deallocate(Dummy); Dummy = LogProbActionAdaptive(1:maxNumFunCallAcceptedRejected)
                    allocate(LogProbActionAdaptive_placeHolder(1:maxNumFunCallAcceptedRejectedNew))!, source = 9999999._RK)
                    LogProbActionAdaptive_placeHolder(1:maxNumFunCallAcceptedRejected) = LogProbActionAdaptive!(1:maxNumFunCallAcceptedRejected)
                    call move_alloc(from=LogProbActionAdaptive_placeHolder,to=LogProbActionAdaptive)
!do i = 1, maxNumFunCallAcceptedRejected
!    if ((Dummy(i)/=LogProbActionAdaptive(i))) then
!        write(*,"(*(g0,:,' '))") "Dummy/=LogProbActionAdaptive", Dummy(i:i)
!        write(*,"(*(g0,:,' '))") "Dummy/=LogProbActionAdaptive", LogProbActionAdaptive(i:i)
!        stop
!    end if
!end do
!end block

!block
!integer :: i
!integer(IK), allocatable :: Dummy(:)
!if (allocated(Dummy)) deallocate(Dummy); Dummy = SampleWeightDISE(1:maxNumFunCallAcceptedRejected)
                    allocate(SampleWeightDISE_placeHolder(1:maxNumFunCallAcceptedRejectedNew), source = 0_IK) ! initialization to 0 is essential
                    SampleWeightDISE_placeHolder(1:maxNumFunCallAcceptedRejected) = SampleWeightDISE!(1:maxNumFunCallAcceptedRejected)
                    call move_alloc(from=SampleWeightDISE_placeHolder,to=SampleWeightDISE)
!do i = 1, maxNumFunCallAcceptedRejected
!    if ((Dummy(i)/=SampleWeightDISE(i))) then
!        write(*,"(*(g0,:,' '))") "Dummy/=SampleWeightDISE", Dummy(i:i)
!        write(*,"(*(g0,:,' '))") "Dummy/=SampleWeightDISE", SampleWeightDISE(i:i)
!        stop
!    end if
!end do
!end block

!block
!integer :: i
!integer(IK), allocatable :: Dummy(:,:)
!if (allocated(Dummy)) deallocate(Dummy); Dummy = TransitionIndexFromTo(1:2,1:maxNumFunCallAcceptedRejected)
                    allocate(TransitionIndexFromTo_placeHolder(1:2,1:maxNumFunCallAcceptedRejectedNew), source = 0_IK)
                    TransitionIndexFromTo_placeHolder(1:2,1:maxNumFunCallAcceptedRejected) = TransitionIndexFromTo!(1:2,1:maxNumFunCallAcceptedRejected)
                    call move_alloc(from=TransitionIndexFromTo_placeHolder,to=TransitionIndexFromTo)
!do i = 1, maxNumFunCallAcceptedRejected
!    if (any(Dummy(:,i)/=TransitionIndexFromTo(:,i))) then
!        write(*,"(*(g0,:,' '))") "Dummy/=TransitionIndexFromTo", Dummy(:,i:i)
!        write(*,"(*(g0,:,' '))") "Dummy/=TransitionIndexFromTo", TransitionIndexFromTo(:,i:i)
!        stop
!    end if
!end do
!end block

                    maxNumFunCallAcceptedRejected = maxNumFunCallAcceptedRejectedNew

                end if

                numFunCallAcceptedRejectedPlusOne = self%Stats%NumFunCall%acceptedRejected + 1_IK

                LogFuncStateDISE(0:nd,numFunCallAcceptedRejectedPlusOne) = co_LogFuncState(0:nd,0) !counterDRS)

                ! for non-symmetric proposal, this will have to computed separately under acceptance and rejection scenarios
                LogProbTrans(self%Stats%NumFunCall%acceptedRejected)  = self%Proposal%getLogProb( nd = nd &
                                                                                                , counterDRS = 0_IK & !counterDRS &
                                                                                                , StateOld = LogFuncStateDISE(1:nd,numFunCallAcceptedRejectedPlusOne) &
                                                                                                , StateNew = LogFuncStateDISE(1:nd,lastAcceptedStateIndex) &
                                                                                                )

                if (co_AccRate(-1) > -1._RK) then
                    TransitionIndexFromTo(1:2,self%Stats%NumFunCall%acceptedRejected) = [ numFunCallAcceptedRejectedPlusOne, lastAcceptedStateIndex ] 
                    logProbActionAdaptive(self%Stats%NumFunCall%acceptedRejected) = min( 0._RK, LogFuncStateDISE(0,lastAcceptedStateIndex) - LogFuncStateDISE(0,numFunCallAcceptedRejectedPlusOne) )
                    SampleWeightDISE(numFunCallAcceptedRejectedPlusOne) = 1_IK
                    lastAcceptedStateIndex = numFunCallAcceptedRejectedPlusOne
                else
                    TransitionIndexFromTo(1:2,self%Stats%NumFunCall%acceptedRejected) = [ lastAcceptedStateIndex, numFunCallAcceptedRejectedPlusOne ] 
                    SampleWeightDISE(lastAcceptedStateIndex) = SampleWeightDISE(lastAcceptedStateIndex) + 1_IK
                    logProbActionAdaptive(self%Stats%NumFunCall%acceptedRejected) = exp(LogFuncStateDISE(0,numFunCallAcceptedRejectedPlusOne) - LogFuncStateDISE(0,lastAcceptedStateIndex))
                    if (logProbActionAdaptive(self%Stats%NumFunCall%acceptedRejected)>1._RK) then
                        write(*,"(A)") "logProbActionAdaptive(self%Stats%NumFunCall%acceptedRejected)>0._RK"
                        write(*,*) logProbActionAdaptive(self%Stats%NumFunCall%acceptedRejected)
                    end if
                    if (logProbActionAdaptive(self%Stats%NumFunCall%acceptedRejected)/=co_AccRate(0)) then
                        write(*,"(A)") "logProbActionAdaptive(self%Stats%NumFunCall%acceptedRejected)/=co_AccRate(0)"
                        write(*,*) logProbActionAdaptive(self%Stats%NumFunCall%acceptedRejected), co_AccRate(0)
                    end if
                    logProbActionAdaptive(self%Stats%NumFunCall%acceptedRejected) = log(1._RK-logProbActionAdaptive(self%Stats%NumFunCall%acceptedRejected))
                    !logProbActionAdaptive(self%Stats%NumFunCall%acceptedRejected) = log(1._RK-co_AccRate(0)) !counterDRS))
                end if
                LogProbTrans(self%Stats%NumFunCall%acceptedRejected) = LogProbTrans(self%Stats%NumFunCall%acceptedRejected) + logProbActionAdaptive(self%Stats%NumFunCall%acceptedRejected)
#endif

        