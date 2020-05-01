!**********************************************************************************************************************************
!**********************************************************************************************************************************
!
!  ParaMonte: plain powerful parallel Monte Carlo library.
!
!  Copyright (C) 2012-present, The Computational Data Science Lab
!
!  This file is part of the ParaMonte library. 
!
!  ParaMonte is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, version 3 of the License.
!
!  ParaMonte is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************

submodule (ParaDRAM_mod) Kernel_smod

    use, intrinsic :: iso_fortran_env, only: output_unit
    !use Constants_mod, only: IK, RK ! gfortran 9.3 compile crashes with this line
    use Err_mod, only: Err_type
#if defined MPI_ENABLED
    use mpi
#endif

    implicit none

    character(*), parameter :: SUBMODULE_NAME = MODULE_NAME // "@Kernel_smod"

    type            :: SumAccRateSinceStart_type
        real(RK)    :: acceptedRejected
        real(RK)    :: acceptedRejectedDelayed
    end type SumAccRateSinceStart_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    module subroutine runKernel ( PD            &
                                , getLogFunc    &
                                )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: runKernel
#endif
        use Constants_mod, only: RK, IK, NEGINF_RK, NLC
        use String_mod, only: num2str
        use Math_mod, only: getCumSum, getLogSubExp
        use ParaMonteLogFunc_mod, only: getLogFunc_proc

        implicit none

        character(*), parameter             :: PROCEDURE_NAME = SUBMODULE_NAME//"@runKernel()"
        integer(IK) , parameter             :: CHAIN_RESTART_OFFSET = 2_IK

        class(ParaDRAM_type), intent(inout) :: PD
        procedure(getLogFunc_proc)          :: getLogFunc

#if defined CAF_ENABLED
        real(RK)    , allocatable           :: co_AccRate(:)[:]
        real(RK)    , allocatable           :: co_LogFuncState(:,:)[:]                          ! (0:nd,-1:delayedRejectionCount), -1 corresponds to the current accepted state
        integer(IK) , save                  :: co_proposalFound_samplerUpdateOccurred(2)[*]     ! merging these scalars would reduce the MPI communication overhead cost: co_proposalFound, co_samplerUpdateOccurred, co_counterDRS, 0 means false, 1 means true
#else
        real(RK)    , allocatable           :: co_LogFuncState(:,:)
        real(RK)    , allocatable           :: co_AccRate(:)
        integer(IK) , save                  :: co_proposalFound_samplerUpdateOccurred(2)        ! merging these scalars would reduce the MPI communication overhead cost: co_proposalFound, co_samplerUpdateOccurred, co_counterDRS, 0 means false, 1 means true
#endif
        type(SumAccRateSinceStart_type)     :: SumAccRateSinceStart                             ! used to figure out the average acceptance ratio for the entire chain.
        integer(IK)                         :: numFunCallAcceptedLastAdaptation                 ! number of function calls accepted at Last proposal adaptation occurrence
        integer(IK)                         :: counterAUP                                       ! counter for adaptiveUpdatePeriod
        integer(IK)                         :: counterAUC                                       ! counter for adaptiveUpdateCount
        integer(IK)                         :: counterPRP                                       ! counter for progressReportPeriod
        integer(IK)                         :: counterDRS                                       ! counter for Delayed Rejection Stages
        integer(IK)                         :: lastState                                        ! dummy temporary argument to hold the value of PD%Stats%NumFunCall%accepted - 1_IK
        integer(IK)                         :: lastStateWeight                                  ! This is used for passing the most recent verbose chain segment to the adaptive updater of the sampler
        integer(IK)                         :: currentStateWeight                               ! counter for SampleWeight, used only in in restart mode
        integer(IK)                         :: numFunCallAcceptedPlusOne                        ! counter for SampleWeight, used only in in restart mode
        integer(IK)                         :: numFunCallAcceptedRejectedLastReport             ! used for progress-report: must be initialized to zero upon entry to the procedure
        real(RK)                            :: timeElapsedUntilLastReportInSeconds              ! used for progress-report: must be initialized to zero upon entry to the procedure
        real(RK)                            :: inverseProgressReportPeriod                      ! used for progress-report: inverse of progressReportPeriod
        real(RK)                            :: sumAccRateLastReport                             ! used for progress-report: must be initialized to zero upon entry to the procedure
        real(RK)                            :: uniformRnd                                       ! used for random number generation.
        real(RK)                            :: meanAccRateSinceStart                            ! used for restart file read
        real(RK)                            :: adaptationMeasure
        real(RK)                            :: adaptationMeasureDummy
        real(RK)                            :: maxLogFuncRejectedProposal
        logical                             :: samplerUpdateIsGreedy
        logical                             :: samplerUpdateSucceeded
        logical                             :: delayedRejectionRequested
        logical                             :: noDelayedRejectionRequested
        integer(IK)                         :: acceptedRejectedDelayedUnusedRestartMode
        integer(IK)                         :: j, imageID, dummy
        integer(IK)                         :: nd
#if defined CAF_ENABLED || defined MPI_ENABLED
        integer(IK)                         :: imageStartID, imageEndID
#if defined CAF_ENABLED
        logical , save                      :: co_proposalFoundSinglChainMode[*]                ! used in the delayed rejection section
#elif defined MPI_ENABLED
        logical                             :: co_proposalFoundSinglChainMode                   ! used in the delayed rejection section
        real(RK), allocatable               :: AccRateMatrix(:,:)                               ! matrix of size (-1:PD%SpecDRAM%DelayedRejectionCount%val,1:PD%Image%count)
        integer(IK)                         :: ndPlusOne
        integer(IK)                         :: ierrMPI
        integer(IK)                         :: delayedRejectionCountPlusTwo
#endif
        PD%Stats%avgCommTimePerFunCall = 0._RK                                                  ! Until the reporting time, this is in reality, sumCommTimePerFunCall. This is meaningful only in singlChain parallelism
        PD%Stats%NumFunCall%acceptedRejectedDelayedUnused = PD%Image%count                      ! used only in singlChain parallelism, and relevant only on the first image
#endif
        acceptedRejectedDelayedUnusedRestartMode = 0_IK                                         ! used to compute more accurate timings in the restart mode
        PD%Stats%avgTimePerFunCalInSec = 0._RK
        numFunCallAcceptedRejectedLastReport = 0_IK
        timeElapsedUntilLastReportInSeconds = 0._RK
        inverseProgressReportPeriod = 1._RK / real(PD%SpecBase%ProgressReportPeriod%val,kind=RK) ! this remains a constant except for the last the last report of the simulation
        sumAccRateLastReport = 0._RK
        nd = PD%nd%val
        
        if (allocated(co_AccRate)) deallocate(co_AccRate)
        if (allocated(co_LogFuncState)) deallocate(co_LogFuncState)
#if defined CAF_ENABLED
        allocate(co_LogFuncState(0:nd,-1:PD%SpecDRAM%DelayedRejectionCount%val)[*])
        allocate(co_AccRate(-1:PD%SpecDRAM%DelayedRejectionCount%val)[*])                       ! the negative element will contain counterDRS
#else
        allocate(co_LogFuncState(0:nd,-1:PD%SpecDRAM%DelayedRejectionCount%val))
        allocate(co_AccRate(-1:PD%SpecDRAM%DelayedRejectionCount%val))
#endif
        co_AccRate(-1)  = 0._RK                                                                 ! the real-value counterDRS, indicating the initial delayed rejection stage at which the first point is sampled
        co_AccRate(0)   = 1._RK                                                                 ! initial acceptance rate for the first zeroth DR stage.
        co_AccRate(1:PD%SpecDRAM%DelayedRejectionCount%val) = 0._RK                             ! indicates the very first proposal acceptance on image 1

#if defined MPI_ENABLED
        if (allocated(AccRateMatrix)) deallocate(AccRateMatrix)
        allocate(AccRateMatrix(-1:PD%SpecDRAM%DelayedRejectionCount%val,1:PD%Image%count))      ! the negative element will contain counterDRS
        AccRateMatrix = 0._RK                                                                   ! -huge(1._RK)  ! debug
        AccRateMatrix(0,1:PD%Image%count) = 1._RK                                               ! initial acceptance rate for the first zeroth DR stage.
        ndPlusOne = nd + 1_IK
        delayedRejectionCountPlusTwo = PD%SpecDRAM%DelayedRejectionCount%val + 2_IK
#endif

        delayedRejectionRequested                       = PD%SpecDRAM%DelayedRejectionCount%val > 0_IK
        noDelayedRejectionRequested                     = .not. delayedRejectionRequested
        if (delayedRejectionRequested) then
        PD%Stats%NumFunCall%acceptedRejectedDelayed     = 0_IK                                  ! Markov Chain counter
        SumAccRateSinceStart%acceptedRejectedDelayed    = 0._RK                                 ! sum of acceptance rate
        end if

        adaptationMeasure                               = 0._RK                                 ! needed for the first output
        SumAccRateSinceStart%acceptedRejected           = 0._RK                                 ! sum of acceptance rate
        PD%Stats%NumFunCall%acceptedRejected            = 0_IK                                  ! Markov Chain counter
        counterAUC                                      = 0_IK                                  ! counter for padaptiveUpdateCount.
        counterPRP                                      = 0_IK                                  ! counter for progressReportPeriod.
        counterAUP                                      = 0_IK                                  ! counter for adaptiveUpdatePeriod. 
        PD%Stats%NumFunCall%accepted                    = 0_IK                                  ! Markov Chain acceptance counter.
        samplerUpdateSucceeded                          = .true.                                ! needed to set up lastStateWeight and numFunCallAcceptedLastAdaptation for the first accepted proposal
        numFunCallAcceptedLastAdaptation                = 0_IK
        lastStateWeight                                 = -huge(lastStateWeight)

        call PD%Timer%tic()

        blockDryRunSetup: if (PD%isFreshRun) then

            allocate(PD%Chain%ProcessID     (   PD%SpecMCMC%ChainSize%val))
            allocate(PD%Chain%DelRejStage   (   PD%SpecMCMC%ChainSize%val))
            allocate(PD%Chain%MeanAccRate   (   PD%SpecMCMC%ChainSize%val))
            allocate(PD%Chain%Adaptation    (   PD%SpecMCMC%ChainSize%val))
            allocate(PD%Chain%BurninLoc     (   PD%SpecMCMC%ChainSize%val))
            allocate(PD%Chain%Weight        (   PD%SpecMCMC%ChainSize%val))
            allocate(PD%Chain%State         (nd,PD%SpecMCMC%ChainSize%val))
            allocate(PD%Chain%LogFunc       (   PD%SpecMCMC%ChainSize%val))

        else blockDryRunSetup

            ! load the existing Chain file into PD%Chain components

            call PD%Chain%get   ( chainFilePath = PD%ChainFile%Path%original        &
                                , chainFileForm = PD%SpecBase%ChainFileFormat%val   &
                                , Err = PD%Err                                      &
                                , resizedChainSize = PD%SpecMCMC%ChainSize%val      &
                                , lenHeader = PD%Chain%lenHeader                    &
                                , ndim = nd                                         &
                                , delimiter = PD%SpecBase%OutputDelimiter%val       &
                                )
            if (PD%Err%occurred) then
                PD%Err%msg = PROCEDURE_NAME//PD%Err%msg
                call PD%abort( Err = PD%Err, prefix = PD%brand, newline = NLC, outputUnit = PD%LogFile%unit )
            end if

            if (PD%Chain%Count%compact<=CHAIN_RESTART_OFFSET) then
                PD%isFreshRun = .true.
                PD%isDryRun = .not. PD%isFreshRun
            end if

#if defined CAF_ENABLED
            sync all
#elif defined MPI_ENABLED
            call mpi_barrier(mpi_comm_world,ierrMPI)
#endif

            if (PD%Image%isMaster) then

                ! set up the chain file

                block

                    use System_mod, only: RandomFileName_type, copyFile, removeFile
                    type(RandomFileName_type) :: RFN

                    ! create a copy of the chain file, just for the sake of not losing the simulation results

                    RFN = RandomFileName_type(dir = "", key = PD%ChainFile%Path%original//".temporary_restart_copy", ext="") 
                    call copyFile(pathOld=PD%ChainFile%Path%original,pathNew=RFN%path,isWindows=PD%OS%isWindows,Err=PD%err)
                    if (PD%Err%occurred) then
                        PD%Err%msg = PROCEDURE_NAME//PD%Err%msg
                        call PD%abort( Err = PD%Err, prefix = PD%brand, newline = NLC, outputUnit = PD%LogFile%unit )
                    end if

                    ! reopen the chain file to resume the simulation

                    open( newunit = PD%ChainFile%unit               &
                        , file = PD%ChainFile%Path%original         &
                        , form = PD%ChainFile%Form%value            &
                        , status = PD%ChainFile%status              &
                        , iostat = PD%ChainFile%Err%stat            &
                        , position = PD%ChainFile%Position%value    )
                    PD%Err = PD%ChainFile%getOpenErr(PD%ChainFile%Err%stat)
                    if (PD%Err%occurred) then
                        PD%Err%msg = PROCEDURE_NAME//": Error occurred while opening "//PD%name//" "//PD%ChainFile%suffix//" file='"//PD%ChainFile%Path%original//"'."
                        call PD%abort( Err = PD%Err, prefix = PD%brand, newline = NLC, outputUnit = PD%LogFile%unit )
                    end if

                    ! rewrite the chain file

                    call PD%Chain%writeChainFile( ndim = nd &
                                                , compactStartIndex = 1_IK &
                                                , compactEndIndex = PD%Chain%Count%compact - CHAIN_RESTART_OFFSET &
                                                , chainFileUnit = PD%ChainFile%unit &
                                                , chainFileForm = PD%SpecBase%ChainFileFormat%val &
                                                , chainFileFormat = PD%ChainFile%format &
                                                )

                    ! remove the temporary copy of the chain file

                    call removeFile(path=RFN%path,isWindows=PD%OS%isWindows,Err=PD%Err)
                    if (PD%Err%occurred) then
                        PD%Err%msg = PROCEDURE_NAME//PD%Err%msg
                        call PD%abort( Err = PD%Err, prefix = PD%brand, newline = NLC, outputUnit = PD%LogFile%unit )
                    end if

                end block

            end if

        end if blockDryRunSetup

        if (PD%isFreshRun) then ! this must be done separately from the above blockDryRunSetup
            
            PD%Chain%BurninLoc(1)  = 1_IK

            co_LogFuncState(1:nd,0) = PD%SpecMCMC%StartPointVec%Val    ! proposal state
            call PD%Timer%toc()
            co_LogFuncState(0,0) = getLogFunc( nd, co_LogFuncState(1:nd,0) )    ! proposal logFunc
            call PD%Timer%toc(); PD%Stats%avgTimePerFunCalInSec = PD%Stats%avgTimePerFunCalInSec + PD%Timer%Time%delta

        else

            co_LogFuncState(1:nd,0)     = PD%Chain%State(1:nd,1)    ! proposal logFunc
            co_LogFuncState(0,0)        = PD%Chain%LogFunc(1)       ! proposal logFunc

        end if

        co_LogFuncState(0:nd,-1) = co_LogFuncState(0:nd,0)  ! set current logFunc and State equal to the first proposal
        PD%Stats%LogFuncMode%val = -huge(PD%Stats%LogFuncMode%val)
        PD%Stats%LogFuncMode%Loc%compact = 0_IK

        if (PD%Image%isFirst) then
            write(output_unit,"(25X,*(A25,3X))" ) "Accepted/Total Func. Call" &
                                                , "Dynamic/Overall Acc. Rate" &
                                                , "Elapsed/Remained Time [s]"
            write(output_unit,"(25X,*(A25,3X))" ) "=========================" &
                                                , "=========================" &
                                                , "========================="
            !call execute_command_line(" ")
            flush(output_unit)
        end if

#if defined CAF_ENABLED || defined MPI_ENABLED
        if (PD%SpecBase%ParallelizationModel%isSinglChain) then
            imageStartID = 1_IK
            imageEndID = PD%Image%count
        else ! isMultiChain
            imageStartID = PD%Image%id
            imageEndID = PD%Image%id
        end if
#else
        imageID = 1_IK ! needed even in the case of serial run to assign a proper value to PD%Chain%ProcessID
#endif

        !******************************************** start of loopMarkovChain *****************************************************
        !***************************************************************************************************************************
        !***************************************************************************************************************************

        loopMarkovChain: do

            co_proposalFound_samplerUpdateOccurred(2) = 0_IK ! at each iteration assume no samplerUpdateOccurred, unless it occurs

#if defined CAF_ENABLED || defined MPI_ENABLED
            blockMasterImage: if (PD%Image%isMaster) then
#endif

                co_proposalFound_samplerUpdateOccurred(1) = 0_IK ! co_proposalFound = .false.
                samplerUpdateIsGreedy = counterAUC < PD%SpecDRAM%GreedyAdaptationCount%val

#if defined CAF_ENABLED || defined MPI_ENABLED
                loopOverImages: do imageID = imageStartID, imageEndID
#if defined CAF_ENABLED
                    call PD%Timer%toc()
                    if (imageID/=PD%Image%id) co_AccRate(-1:PD%SpecDRAM%DelayedRejectionCount%val) = co_AccRate(-1:PD%SpecDRAM%DelayedRejectionCount%val)[imageID] ! happens only in isSinglChain=TRUE
                    call PD%Timer%toc(); PD%Stats%avgCommTimePerFunCall = PD%Stats%avgCommTimePerFunCall + PD%Timer%Time%delta
#elif defined MPI_ENABLED
                    if (imageID/=PD%Image%id) co_AccRate(-1:PD%SpecDRAM%DelayedRejectionCount%val) = AccRateMatrix(-1:PD%SpecDRAM%DelayedRejectionCount%val,imageID) ! happens only in isSinglChain=TRUE
#endif
#endif

                    counterDRS = nint(co_AccRate(-1),kind=IK)
                    if (counterDRS > -1_IK) co_proposalFound_samplerUpdateOccurred(1) = 1_IK ! co_proposalFound = .true.

                    !**************************************** blockProposalAccepted ************************************************
                    !***************************************************************************************************************

                    ! On the very first iteration, this block is (and must be) executed for imageID==1,
                    ! since it is for the first (starting) point, which is assumed to have been accepted
                    ! as the first point by the first coarray imageID.

                    blockProposalAccepted: if ( co_proposalFound_samplerUpdateOccurred(1) == 1_IK ) then ! co_proposalAccepted = .true.

                        lastState = PD%Stats%NumFunCall%accepted
                        PD%Stats%NumFunCall%accepted = PD%Stats%NumFunCall%accepted + 1_IK
                        currentStateWeight = 0_IK

                        ! communicate the accepted logFunc and State from the winning image to master/all images: co_LogFuncState

#if defined MPI_ENABLED
                        if (PD%SpecBase%ParallelizationModel%isSinglChain) then
                            call PD%Timer%toc()
                            ! broadcast winning image to all processes
                            call mpi_bcast  ( imageID           &   ! buffer
                                            , 1                 &   ! count
                                            , mpi_integer       &   ! datatype
                                            , 0                 &   ! root: broadcasting rank
                                            , mpi_comm_world    &   ! comm
                                            , ierrMPI           &   ! ierr
                                            )
                            ! broadcast co_LogFuncState from the winning image to all others
                            call mpi_bcast  ( co_LogFuncState(0:nd,-1)          &   ! buffer
                                            , ndPlusOne                         &   ! count
                                            , mpi_double_precision              &   ! datatype
                                            , imageID - 1_IK                    &   ! root: broadcasting rank
                                            , mpi_comm_world                    &   ! comm
                                            , ierrMPI                           &   ! ierr
                                            )
                            call PD%Timer%toc(); PD%Stats%avgCommTimePerFunCall = PD%Stats%avgCommTimePerFunCall + PD%Timer%Time%delta
                        end if
#elif defined CAF_ENABLED
                        if (imageID/=PD%Image%id) then ! Avoid remote connection for something that is available locally.
                            call PD%Timer%toc()
                            co_LogFuncState(0:nd,-1) = co_LogFuncState(0:nd,-1)[imageID]
                            call PD%Timer%toc(); PD%Stats%avgCommTimePerFunCall = PD%Stats%avgCommTimePerFunCall + PD%Timer%Time%delta
                        end if
#endif

                        ! Note: after every adaptive update of the sampler, counterAUP is reset to 0.
                        if (counterAUP==0_IK .and. samplerUpdateSucceeded) then
                            numFunCallAcceptedLastAdaptation = numFunCallAcceptedLastAdaptation + 1_IK
                            lastStateWeight = 0_IK
                        end if

                        blockFreshDryRun: if (PD%isFreshRun) then

                            PD%Chain%ProcessID(PD%Stats%NumFunCall%accepted)    = imageID
                            PD%Chain%DelRejStage(PD%Stats%NumFunCall%accepted)  = counterDRS
                            PD%Chain%Adaptation(PD%Stats%NumFunCall%accepted)   = adaptationMeasure
                            PD%Chain%Weight(PD%Stats%NumFunCall%accepted)       = 0_IK
                            PD%Chain%LogFunc(PD%Stats%NumFunCall%accepted)      = co_LogFuncState(0,-1)
                            PD%Chain%State(1:nd,PD%Stats%NumFunCall%accepted)   = co_LogFuncState(1:nd,-1)

                            ! find the burnin point

                            PD%Chain%BurninLoc(PD%Stats%NumFunCall%accepted) = getBurninLoc ( lenLogFunc    = PD%Stats%NumFunCall%accepted                      &
                                                                                            , refLogFunc    = PD%Stats%LogFuncMode%val                          &
                                                                                            , LogFunc       = PD%Chain%LogFunc(1:PD%Stats%NumFunCall%accepted)  &
                                                                                            )

                        else blockFreshDryRun ! in restart mode: determine the correct value of co_proposalFound_samplerUpdateOccurred(1)

                            if (PD%Stats%NumFunCall%accepted==PD%Chain%Count%compact) then
                                PD%isFreshRun = .true.
                                PD%isDryRun = .not. PD%isFreshRun
                                PD%Chain%Weight(PD%Stats%NumFunCall%accepted) = 0_IK
                                SumAccRateSinceStart%acceptedRejected = PD%Chain%MeanAccRate(PD%Stats%NumFunCall%accepted-1) * real(PD%Stats%NumFunCall%acceptedRejected,kind=RK)
                                if (delayedRejectionRequested) SumAccRateSinceStart%acceptedRejectedDelayed = PD%Chain%MeanAccRate(PD%Stats%NumFunCall%accepted-1) * real(PD%Stats%NumFunCall%acceptedRejectedDelayed,kind=RK)
                            end if
                            numFunCallAcceptedPlusOne = PD%Stats%NumFunCall%accepted + 1_IK

                        end if blockFreshDryRun

                        if (PD%Stats%LogFuncMode%val < PD%Chain%LogFunc(PD%Stats%NumFunCall%accepted)) then
                            PD%Stats%LogFuncMode%val = PD%Chain%LogFunc(PD%Stats%NumFunCall%accepted)
                            PD%Stats%LogFuncMode%Loc%compact = PD%Stats%NumFunCall%accepted
                        end if

                        SumAccRateSinceStart%acceptedRejected = SumAccRateSinceStart%acceptedRejected + co_AccRate(counterDRS)

                    else blockProposalAccepted

                        counterDRS = PD%SpecDRAM%DelayedRejectionCount%val
                        SumAccRateSinceStart%acceptedRejected = SumAccRateSinceStart%acceptedRejected + co_AccRate(counterDRS)

                    end if blockProposalAccepted

                    !***************************************************************************************************************
                    !**************************************** blockProposalAccepted ************************************************

                    counterAUP = counterAUP + 1_IK
                    counterPRP = counterPRP + 1_IK
                    currentStateWeight = currentStateWeight + 1_IK
                    PD%Stats%NumFunCall%acceptedRejected = PD%Stats%NumFunCall%acceptedRejected + 1_IK

                    if (delayedRejectionRequested) then
                        SumAccRateSinceStart%acceptedRejectedDelayed = SumAccRateSinceStart%acceptedRejectedDelayed + sum(co_AccRate(0:counterDRS))
                        PD%Stats%NumFunCall%acceptedRejectedDelayed = PD%Stats%NumFunCall%acceptedRejectedDelayed + counterDRS + 1_IK
                    end if

                    if (PD%isFreshRun) then ! these are used for adaptive proposal updating, so they have to be set on every accepted or rejected iteration (excluding delayed rejections)
                        PD%Chain%MeanAccRate(PD%Stats%NumFunCall%accepted)  = SumAccRateSinceStart%acceptedRejected / real(PD%Stats%NumFunCall%acceptedRejected,kind=RK)
                        PD%Chain%Weight(PD%Stats%NumFunCall%accepted)       = PD%Chain%Weight(PD%Stats%NumFunCall%accepted) + 1_IK
                    else
#if defined CAF_ENABLED || defined MPI_ENABLED
                        acceptedRejectedDelayedUnusedRestartMode = PD%Stats%NumFunCall%acceptedRejectedDelayedUnused
#else
                        acceptedRejectedDelayedUnusedRestartMode = PD%Stats%NumFunCall%acceptedRejectedDelayed
#endif
                    end if

                    if (counterPRP == PD%SpecBase%ProgressReportPeriod%val) then
                        counterPRP = 0_IK
                        call reportProgress()
                    end if

                    !******************************************** Proposal Adaptation **********************************************
                    !***************************************************************************************************************

                    blockSamplerAdaptation: if ( counterAUC < PD%SpecDRAM%AdaptiveUpdateCount%val .and. counterAUP == PD%SpecDRAM%AdaptiveUpdatePeriod%val ) then

                        co_proposalFound_samplerUpdateOccurred(2) = 1_IK ! istart = numFunCallAcceptedLastAdaptation ! = max( numFunCallAcceptedLastAdaptation , PD%Chain%BurninLoc(PD%Stats%NumFunCall%accepted) ) ! this is experimental

                        ! the order in the following two MUST be preserved as occasionally PD%Stats%NumFunCall%accepted = numFunCallAcceptedLastAdaptation

                        dummy = PD%Chain%Weight(PD%Stats%NumFunCall%accepted) ! needed for the restart mode, not needed in the fresh run
                        if (PD%Stats%NumFunCall%accepted==numFunCallAcceptedLastAdaptation) then    ! no new point has been accepted since last time
                            PD%Chain%Weight(numFunCallAcceptedLastAdaptation) = currentStateWeight - lastStateWeight
#if defined DBG_ENABLED && !defined CAF_ENABLED && !defined MPI_ENABLED
                            if (mod(PD%Chain%Weight(numFunCallAcceptedLastAdaptation),PD%SpecDRAM%AdaptiveUpdatePeriod%val)/=0) then
                                write(output_unit,"(*(g0,:,' '))"   ) PROCEDURE_NAME//": Internal error occurred: ", PD%SpecDRAM%AdaptiveUpdatePeriod%val &
                                                                    , PD%Chain%Weight(numFunCallAcceptedLastAdaptation), currentStateWeight, lastStateWeight
                                !call execute_command_line(" ")
                                flush(output_unit)
                                error stop
                            end if
#endif
                        else
                            PD%Chain%Weight(numFunCallAcceptedLastAdaptation) = PD%Chain%Weight(numFunCallAcceptedLastAdaptation) - lastStateWeight
                            PD%Chain%Weight(PD%Stats%NumFunCall%accepted) = currentStateWeight ! needed for the restart mode, not needed in the fresh run
                        end if

                        if (PD%isFreshRun) then
                            meanAccRateSinceStart = PD%Chain%MeanAccRate(PD%Stats%NumFunCall%accepted)
                            if (PD%SpecBase%RestartFileFormat%isBinary) then
                                write(PD%RestartFile%unit) meanAccRateSinceStart
                            else
                                write(PD%RestartFile%unit,PD%RestartFile%format) "meanAccRateSinceStart", meanAccRateSinceStart
                                call PD%Proposal%writeRestartFile()
                            end if
                            flush(PD%RestartFile%unit)
                        else
                            if (PD%SpecBase%RestartFileFormat%isBinary) then
                                read(PD%RestartFile%unit) meanAccRateSinceStart
                            else
                                !read(PD%RestartFile%unit,PD%RestartFile%format)
                                !read(PD%RestartFile%unit,PD%RestartFile%format) meanAccRateSinceStart
                                read(PD%RestartFile%unit,*)
                                read(PD%RestartFile%unit,*) meanAccRateSinceStart
                                call PD%Proposal%readRestartFile()
                            end if
                            SumAccRateSinceStart%acceptedRejected = meanAccRateSinceStart * real(PD%Stats%NumFunCall%acceptedRejected,kind=RK)
                        end if

                        call PD%Proposal%doAdaptation   ( nd                        = nd                                                                                    &
                                                        , chainSize                 = PD%Stats%NumFunCall%accepted - numFunCallAcceptedLastAdaptation + 1_IK                &
                                                        , Chain                     = PD%Chain%State(1:nd,numFunCallAcceptedLastAdaptation:PD%Stats%NumFunCall%accepted)    &
                                                        , ChainWeight               = PD%Chain%Weight(numFunCallAcceptedLastAdaptation:PD%Stats%NumFunCall%accepted)        &
                                                        , samplerUpdateIsGreedy     = samplerUpdateIsGreedy                                                                 &
                                                        , meanAccRateSinceStart     = meanAccRateSinceStart                                                                 &
                                                        , samplerUpdateSucceeded    = samplerUpdateSucceeded                                                                &
                                                        , adaptationMeasure         = adaptationMeasureDummy                                                                &
                                                        )

                        PD%Chain%Weight(PD%Stats%NumFunCall%accepted) = dummy   ! needed for the restart mode, not needed in the fresh run
                        if (PD%Stats%NumFunCall%accepted==numFunCallAcceptedLastAdaptation) then
                            adaptationMeasure = adaptationMeasure + adaptationMeasureDummy ! this is the worst-case upper-bound
                        else
                            adaptationMeasure = adaptationMeasureDummy
                            PD%Chain%Weight(numFunCallAcceptedLastAdaptation) = PD%Chain%Weight(numFunCallAcceptedLastAdaptation) + lastStateWeight
                        end if
                        if (samplerUpdateSucceeded) then
                            lastStateWeight = currentStateWeight ! PD%Chain%Weight(PD%Stats%NumFunCall%accepted) ! informative, do not remove
                            numFunCallAcceptedLastAdaptation = PD%Stats%NumFunCall%accepted
                        end if

                        counterAUP = 0_IK
                        counterAUC = counterAUC + 1_IK
                        if (counterAUC==PD%SpecDRAM%AdaptiveUpdateCount%val) adaptationMeasure = 0._RK

                    !else blockSamplerAdaptation

                        !adaptationMeasure = 0._RK

                    end if blockSamplerAdaptation

                    !***************************************************************************************************************
                    !******************************************** Proposal Adaptation **********************************************

#if defined CAF_ENABLED || defined MPI_ENABLED
                    if (co_proposalFound_samplerUpdateOccurred(1)==1_IK) exit loopOverImages

                end do loopOverImages
#endif

#if defined MPI_ENABLED
                if (PD%SpecBase%ParallelizationModel%isSinglChain .and. co_proposalFound_samplerUpdateOccurred(1)==0_IK) then
                    imageID = 0_IK  ! broadcast rank # 0 to all processes, indicating unsuccessful sampling
                    call mpi_bcast  ( imageID           &   ! buffer
                                    , 1                 &   ! count
                                    , mpi_integer       &   ! datatype
                                    , 0                 &   ! root: broadcasting rank
                                    , mpi_comm_world    &   ! comm
                                    , ierrMPI           &   ! ierr
                                    )
                end if
#endif

                !************************************************* output write ****************************************************
                !*******************************************************************************************************************

                ! if new point has been sampled, write the previous sampled point to output file

                blockOutputWrite: if (PD%isFreshRun .and. co_proposalFound_samplerUpdateOccurred(1)==1_IK .and. lastState>0_IK) then
                    if (PD%SpecBase%ChainFileFormat%isCompact) then
                        write(PD%ChainFile%unit,PD%ChainFile%format ) PD%Chain%ProcessID(lastState)     &
                                                                    , PD%Chain%DelRejStage(lastState)   &
                                                                    , PD%Chain%MeanAccRate(lastState)   &
                                                                    , PD%Chain%Adaptation(lastState)    &
                                                                    , PD%Chain%BurninLoc(lastState)     &
                                                                    , PD%Chain%Weight(lastState)        &
                                                                    , PD%Chain%LogFunc(lastState)       &
                                                                    , PD%Chain%State(1:nd,lastState)
                    elseif (PD%SpecBase%ChainFileFormat%isBinary) then
                        write(PD%ChainFile%unit                     ) PD%Chain%ProcessID(lastState)     &
                                                                    , PD%Chain%DelRejStage(lastState)   &
                                                                    , PD%Chain%MeanAccRate(lastState)   &
                                                                    , PD%Chain%Adaptation(lastState)    &
                                                                    , PD%Chain%BurninLoc(lastState)     &
                                                                    , PD%Chain%Weight(lastState)        &
                                                                    , PD%Chain%LogFunc(lastState)       &
                                                                    , PD%Chain%State(1:nd,lastState)
                    elseif (PD%SpecBase%ChainFileFormat%isVerbose) then
                        do j = 1, PD%Chain%Weight(lastState)
                        write(PD%ChainFile%unit,PD%ChainFile%format ) PD%Chain%ProcessID(lastState)     &
                                                                    , PD%Chain%DelRejStage(lastState)   &
                                                                    , PD%Chain%MeanAccRate(lastState)   &
                                                                    , PD%Chain%Adaptation(lastState)    &
                                                                    , PD%Chain%BurninLoc(lastState)     &
                                                                    , 1_IK                              &
                                                                    , PD%Chain%LogFunc(lastState)       &
                                                                    , PD%Chain%State(1:nd,lastState)
                        end do
                    end if
                end if blockOutputWrite

                !*******************************************************************************************************************
                !************************************************* output write ****************************************************

                !********************************************* last output write ***************************************************
                !*******************************************************************************************************************

                blockLastSample: if (PD%Stats%NumFunCall%accepted==PD%SpecMCMC%ChainSize%val) then !co_missionAccomplished = .true.

                    ! on 3 images Windows, sustituting co_missionAccomplished with the following leads to 10% less communication overhead for 1D Gaussian example
                    ! on 3 images Linux  , sustituting co_missionAccomplished with the following leads to 16% less communication overhead for 1D Gaussian example
                    ! on 5 images Linux  , sustituting co_missionAccomplished with the following leads to 11% less communication overhead for 1D Gaussian example

                    co_proposalFound_samplerUpdateOccurred(1) = -1_IK  ! equivalent to co_missionAccomplished = .true.
                    inverseProgressReportPeriod = 1._RK / (PD%Stats%NumFunCall%acceptedRejected-numFunCallAcceptedRejectedLastReport)

                    blockLastOutputWrite: if (PD%isFreshRun) then

                        ! write the last sampled point to the output files

                        lastState = PD%Stats%NumFunCall%accepted
                        if (PD%SpecBase%ChainFileFormat%isCompact .or. PD%SpecBase%ChainFileFormat%isVerbose) then
                            write(PD%ChainFile%unit,PD%ChainFile%format ) PD%Chain%ProcessID(lastState)     &
                                                                        , PD%Chain%DelRejStage(lastState)   &
                                                                        , PD%Chain%MeanAccRate(lastState)   &
                                                                        , PD%Chain%Adaptation(lastState)    &
                                                                        , PD%Chain%BurninLoc(lastState)     &
                                                                        , 1_IK                              &
                                                                        , PD%Chain%LogFunc(lastState)       &
                                                                        , PD%Chain%State(1:nd,lastState)
                        elseif (PD%SpecBase%ChainFileFormat%isBinary) then
                            write(PD%ChainFile%unit                     ) PD%Chain%ProcessID(lastState)     &
                                                                        , PD%Chain%DelRejStage(lastState)   &
                                                                        , PD%Chain%MeanAccRate(lastState)   &
                                                                        , PD%Chain%Adaptation(lastState)    &
                                                                        , PD%Chain%BurninLoc(lastState)     &
                                                                        , 1_IK                              &
                                                                        , PD%Chain%LogFunc(lastState)       &
                                                                        , PD%Chain%State(1:nd,lastState)
                        end if
                        flush(PD%ChainFile%unit)

                    end if blockLastOutputWrite

                    call reportProgress()

                end if blockLastSample

                !*******************************************************************************************************************
                !********************************************* last output write ***************************************************

#if defined CAF_ENABLED

                if (PD%SpecBase%ParallelizationModel%isSinglChain) then
                    call PD%Timer%toc()
                    sync images(*)
                    call PD%Timer%toc(); PD%Stats%avgCommTimePerFunCall = PD%Stats%avgCommTimePerFunCall + PD%Timer%Time%delta
                end if

            else blockMasterImage   ! ATTN: This block should be executed only when singlChain parallelizationModel is requested

                sync images(1)

                ! get the accepted proposal from the first image

                call PD%Timer%toc()
                co_proposalFound_samplerUpdateOccurred(1:2) = co_proposalFound_samplerUpdateOccurred(1:2)[1]
                if (co_proposalFound_samplerUpdateOccurred(1)==1_IK) co_LogFuncState(0:nd,-1) = co_LogFuncState(0:nd,-1)[1]
                if (co_proposalFound_samplerUpdateOccurred(2)==1_IK) call PD%Proposal%getAdaptation()
                call PD%Timer%toc(); PD%Stats%avgCommTimePerFunCall = PD%Stats%avgCommTimePerFunCall + PD%Timer%Time%delta

                if (PD%isDryRun) then
                    if (co_proposalFound_samplerUpdateOccurred(1)==1_IK) then
                        PD%Stats%NumFunCall%accepted = PD%Stats%NumFunCall%accepted + 1_IK
                        numFunCallAcceptedPlusOne = PD%Stats%NumFunCall%accepted + 1_IK
                        currentStateWeight = 1_IK
                        if (PD%Stats%NumFunCall%accepted==PD%Chain%Count%compact) then
                            PD%isFreshRun = .true.
                            PD%isDryRun = .false.
                        end if
                    else
                        currentStateWeight = currentStateWeight + PD%Image%count
                    end if
#if defined DBG_ENABLED
                elseif (co_proposalFound_samplerUpdateOccurred(1)==1_IK) then
                    PD%Stats%NumFunCall%accepted = PD%Stats%NumFunCall%accepted + 1_IK
#endif
                end if

            end if blockMasterImage

#elif defined MPI_ENABLED

            else blockMasterImage   ! This block should be executed only when singlChain parallelizationModel is requested

                ! fetch winning rank to all processes

                call PD%Timer%toc()
                call mpi_bcast  ( imageID           &   ! buffer
                                , 1                 &   ! count
                                , mpi_integer       &   ! datatype
                                , 0                 &   ! root: broadcasting rank
                                , mpi_comm_world    &   ! comm
                                , ierrMPI           &   ! ierr
                                )
                call PD%Timer%toc(); PD%Stats%avgCommTimePerFunCall = PD%Stats%avgCommTimePerFunCall + PD%Timer%Time%delta

                if (imageID>0_IK) then ! co_ProposalFound = .true., sampling successful

                    ! broadcast co_LogFuncState from the winning image to all others

                    call PD%Timer%toc()
                    call mpi_bcast  ( co_LogFuncState(0:nd,-1)  &   ! buffer
                                    , ndPlusOne                 &   ! count
                                    , mpi_double_precision      &   ! datatype
                                    , imageID - 1               &   ! root: broadcasting rank
                                    , mpi_comm_world            &   ! comm
                                    , ierrMPI                   &   ! ierr
                                    )
                    call PD%Timer%toc(); PD%Stats%avgCommTimePerFunCall = PD%Stats%avgCommTimePerFunCall + PD%Timer%Time%delta

                end if

                if (PD%isDryRun) then
                    if (imageID>0_IK) then ! equivalent to co_proposalFound_samplerUpdateOccurred(1)==1_IK
                        PD%Stats%NumFunCall%accepted = PD%Stats%NumFunCall%accepted + 1_IK
                        numFunCallAcceptedPlusOne = PD%Stats%NumFunCall%accepted + 1_IK
                        currentStateWeight = 1_IK
                        if (PD%Stats%NumFunCall%accepted==PD%Chain%Count%compact) then
                            PD%isFreshRun = .true.
                            PD%isDryRun = .false.
                        end if
                    else
                        currentStateWeight = currentStateWeight + PD%Image%count
                    end if
#if defined DBG_ENABLED
                elseif (co_proposalFound_samplerUpdateOccurred(1)==1_IK) then
                    PD%Stats%NumFunCall%accepted = PD%Stats%NumFunCall%accepted + 1_IK
#endif
                end if

            end if blockMasterImage

            !***********************************************************************************************************************
            !**************************************** begin Common block between all images ****************************************

            ! broadcast samplerUpdateOccurred from the root process to all others, and broadcast proposal adaptations if needed

            if (PD%SpecBase%ParallelizationModel%isSinglChain) then
                call PD%Timer%toc()
                call mpi_bcast  ( co_proposalFound_samplerUpdateOccurred  &   ! buffer: XXX: first element not needed except to end the simulation.
                                , 2                     &   ! count
                                , mpi_integer           &   ! datatype
                                , 0                     &   ! root: broadcasting rank
                                , mpi_comm_world        &   ! comm
                                , ierrMPI               &   ! ierr
                                )
                if (co_proposalFound_samplerUpdateOccurred(2)==1_IK) call PD%Proposal%getAdaptation()
                call PD%Timer%toc(); PD%Stats%avgCommTimePerFunCall = PD%Stats%avgCommTimePerFunCall + PD%Timer%Time%delta
            end if

#endif

            if (co_proposalFound_samplerUpdateOccurred(1) == -1_IK) exit loopMarkovChain   ! we are done: co_missionAccomplished = .true.

            co_AccRate(-1) = -1._RK ! counterDRS at which new proposal is accepted
            maxLogFuncRejectedProposal = NEGINF_RK

#if defined CAF_ENABLED || defined MPI_ENABLED

            blockSingleChainParallelism: if (PD%SpecBase%ParallelizationModel%isSinglChain) then

                !*******************************************************************************************************************

                loopDelayedRejectionSingleChain: do counterDRS = 0, PD%SpecDRAM%DelayedRejectionCount%val

                    ! PD%Stats%NumFunCall%acceptedRejectedDelayedUnused is relevant only on the first image, despite being updated by all images
                    PD%Stats%NumFunCall%acceptedRejectedDelayedUnused = PD%Stats%NumFunCall%acceptedRejectedDelayedUnused + PD%Image%count

                    co_LogFuncState(1:nd,counterDRS) = PD%Proposal%getNew   ( nd            = nd                                    &
                                                                            , counterDRS    = counterDRS                            &
                                                                            , StateOld      = co_LogFuncState(1:nd,counterDRS-1)    &
                                                                            )

                    call random_number(uniformRnd) ! only for the purpose of restart mode reproducibility

#if defined CAF_ENABLED
                    ! this is necessary to avoid racing condition on co_LogFuncState and co_proposalFoundSinglChainMode
                    if (PD%Image%isMaster) then
                        call PD%Timer%toc()
                        sync images(*)
                        call PD%Timer%toc(); PD%Stats%avgCommTimePerFunCall = PD%Stats%avgCommTimePerFunCall + PD%Timer%Time%delta
                    else
                        sync images(1)
                    end if
#endif

                    if (PD%isFreshRun .or. numFunCallAcceptedPlusOne==PD%Chain%Count%compact) then

                        ! accept or reject the proposed state, only if no acceptance has occurred yet
                        if (co_AccRate(-1)<-0.5_RK) then

                            call PD%Timer%toc()
                            co_LogFuncState(0,counterDRS) = getLogFunc(nd,co_LogFuncState(1:nd,counterDRS))
                            call PD%Timer%toc(); PD%Stats%avgTimePerFunCalInSec = PD%Stats%avgTimePerFunCalInSec + PD%Timer%Time%delta

                            if ( co_LogFuncState(0,counterDRS) >= co_LogFuncState(0,-1) ) then ! accept the proposed state

                                co_AccRate(counterDRS) = 1._RK
                                co_AccRate(-1) = real(counterDRS,kind=RK)
                                co_LogFuncState(0:nd,-1) = co_LogFuncState(0:nd,counterDRS)

                            elseif ( co_LogFuncState(0,counterDRS) < maxLogFuncRejectedProposal ) then  ! reject the proposed state. This step should be reachable only when delayedRejectionCount > 0

                                co_AccRate(counterDRS) = 0._RK  ! proposal rejected

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
                                end if

                            end if

                            maxLogFuncRejectedProposal = max( maxLogFuncRejectedProposal, co_LogFuncState(0,counterDRS) )

                        end if


                    else

                        if (PD%Image%id==PD%Chain%ProcessID(numFunCallAcceptedPlusOne) .and. &
                            currentStateWeight+PD%Image%id-1_IK==PD%Chain%Weight(PD%Stats%NumFunCall%accepted) .and. &
                            counterDRS==PD%Chain%DelRejStage(numFunCallAcceptedPlusOne)) then

                            co_AccRate(-1) = real(counterDRS,kind=RK)
                            co_LogFuncState(   0,counterDRS) = PD%Chain%LogFunc   (numFunCallAcceptedPlusOne)
                            co_LogFuncState(1:nd,counterDRS) = PD%Chain%State(1:nd,numFunCallAcceptedPlusOne)
                            co_LogFuncState(0:nd,-1) = co_LogFuncState(0:nd,counterDRS)

                        end if

                    end if

                    co_proposalFoundSinglChainMode = co_AccRate(-1) > -1._RK

                    if (delayedRejectionRequested) then ! broadcast the sampling status from the first image to all others
#if defined CAF_ENABLED
                        if (PD%Image%isMaster) then ! this is necessary to avoid racing condition on the value of co_proposalFoundSinglChainMode
                            call PD%Timer%toc()
                            sync images(*)
                            call PD%Timer%toc(); PD%Stats%avgCommTimePerFunCall = PD%Stats%avgCommTimePerFunCall + PD%Timer%Time%delta
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

                    if (co_proposalFoundSinglChainMode) exit loopDelayedRejectionSingleChain

                end do loopDelayedRejectionSingleChain

#if defined MPI_ENABLED
                ! gather all accr on the master image. Avoid unnecessary communication when proposal is on the master image
                if ( noDelayedRejectionRequested .or. (delayedRejectionRequested .and. .not.co_proposalFoundSinglChainMode) ) then 
                    call PD%Timer%toc()
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
                    call PD%Timer%toc(); PD%Stats%avgCommTimePerFunCall = PD%Stats%avgCommTimePerFunCall + PD%Timer%Time%delta
                endif
#endif

                !*******************************************************************************************************************

            else blockSingleChainParallelism

                !*******************************************************************************************************************
#endif

                loopDelayedRejection: do counterDRS = 0, PD%SpecDRAM%DelayedRejectionCount%val

                    co_LogFuncState(1:nd,counterDRS) = PD%Proposal%getNew   ( nd            = nd                                    &
                                                                            , counterDRS    = counterDRS                            &
                                                                            , StateOld      = co_LogFuncState(1:nd,counterDRS-1)    &
                                                                            )

                    call random_number(uniformRnd) ! only for the purpose of restart mode reproducibility

                    if (PD%isFreshRun .or. numFunCallAcceptedPlusOne==PD%Chain%Count%compact) then

                        call PD%Timer%toc()
                        co_LogFuncState(0,counterDRS) = getLogFunc(nd,co_LogFuncState(1:nd,counterDRS))
                        call PD%Timer%toc(); PD%Stats%avgTimePerFunCalInSec = PD%Stats%avgTimePerFunCalInSec + PD%Timer%Time%delta

                        ! accept or reject the proposed state

                        if ( co_LogFuncState(0,counterDRS) >= co_LogFuncState(0,-1) ) then ! accept the proposed state

                            co_AccRate(counterDRS) = 1._RK
                            co_AccRate(-1) = real(counterDRS,kind=RK)
                            co_LogFuncState(0:nd,-1) = co_LogFuncState(0:nd,counterDRS)
                            exit loopDelayedRejection

                        elseif ( co_LogFuncState(0,counterDRS) < maxLogFuncRejectedProposal ) then  ! reject the proposed state. This step should be reachable only when delayedRejectionCount > 0

                            co_AccRate(counterDRS) = 0._RK  ! proposal rejected

                        else    ! accept with probability co_AccRate

                            if ( counterDRS == 0_IK ) then ! This should be equivalent to maxLogFuncRejectedProposal == NEGINF_RK
                                co_AccRate(counterDRS) = exp( co_LogFuncState(0,counterDRS) - co_LogFuncState(0,-1) )
                            else    ! ensure no arithmetic overflow/underflow. ATTN: co_LogFuncState(0,-1) > co_LogFuncState(0,counterDRS) > maxLogFuncRejectedProposal
                                co_AccRate(counterDRS) = exp( getLogSubExp( co_LogFuncState(0,counterDRS)   , maxLogFuncRejectedProposal ) &
                                                            - getLogSubExp( co_LogFuncState(0,-1)           , maxLogFuncRejectedProposal ) )
                            end if

                            if (uniformRnd<co_AccRate(counterDRS)) then ! accept the proposed state
                                co_AccRate(-1) = real(counterDRS,kind=RK)
                                co_LogFuncState(0:nd,-1) = co_LogFuncState(0:nd,counterDRS)
                                exit loopDelayedRejection
                            end if

                        end if

                        maxLogFuncRejectedProposal = max( maxLogFuncRejectedProposal, co_LogFuncState(0,counterDRS) )

                    else

                        if (currentStateWeight==PD%Chain%Weight(PD%Stats%NumFunCall%accepted) .and. counterDRS==PD%Chain%DelRejStage(numFunCallAcceptedPlusOne)) then
                            co_AccRate(-1) = real(counterDRS,kind=RK)
                            co_LogFuncState(   0,-1) = PD%Chain%LogFunc   (numFunCallAcceptedPlusOne)
                            co_LogFuncState(1:nd,-1) = PD%Chain%State(1:nd,numFunCallAcceptedPlusOne)
                            exit loopDelayedRejection
                        end if

                    end if

                end do loopDelayedRejection

                !*******************************************************************************************************************

#if defined CAF_ENABLED || defined MPI_ENABLED
            end if blockSingleChainParallelism
#endif

            cycle loopMarkovChain

            !***************************************** end Common block between all images *****************************************
            !***********************************************************************************************************************

        end do loopMarkovChain

        call PD%Timer%toc()

        !***************************************************************************************************************************
        !***************************************************************************************************************************
        !********************************************** end of loopMarkovChain *****************************************************

        PD%chain%Count%resized = PD%Stats%NumFunCall%accepted
        PD%chain%Count%compact = PD%Stats%NumFunCall%accepted
        PD%chain%Count%verbose = PD%Stats%NumFunCall%acceptedRejected

        if (noDelayedRejectionRequested) then
            PD%Stats%NumFunCall%acceptedRejectedDelayed = PD%Stats%NumFunCall%acceptedRejected
            SumAccRateSinceStart%acceptedRejectedDelayed = SumAccRateSinceStart%acceptedRejected
        end if
            

        if (PD%Image%isMaster) then
            block
                integer(IK) :: i
                if (PD%SpecDRAM%BurninAdaptationMeasure%val>0.999999_RK) then
                    PD%Stats%AdaptationBurninLoc%compact = 1_IK
                    PD%Stats%AdaptationBurninLoc%verbose = 1_IK
                else
                    PD%Stats%AdaptationBurninLoc%compact = PD%Stats%NumFunCall%accepted
                    PD%Stats%AdaptationBurninLoc%verbose = PD%Stats%NumFunCall%acceptedRejected - PD%Chain%Weight(PD%Stats%NumFunCall%accepted) + 1_IK
                    loopAdaptationBurnin: do i = PD%Stats%NumFunCall%accepted-1, 1, -1
                        if (PD%Chain%Adaptation(i)>PD%SpecDRAM%BurninAdaptationMeasure%val) exit loopAdaptationBurnin
                        PD%Stats%AdaptationBurninLoc%compact = PD%Stats%AdaptationBurninLoc%compact - 1_IK
                        PD%Stats%AdaptationBurninLoc%verbose = PD%Stats%AdaptationBurninLoc%verbose - PD%Chain%Weight(i)
                    end do loopAdaptationBurnin
                end if
            end block
            PD%Stats%BurninLoc%compact          = PD%Chain%BurninLoc(PD%Stats%NumFunCall%accepted)
            PD%Stats%BurninLoc%verbose          = sum(PD%Chain%Weight(1:PD%Stats%BurninLoc%compact-1)) + 1_IK
            PD%Stats%LogFuncMode%Crd            = PD%Chain%State(1:nd,PD%Stats%LogFuncMode%Loc%compact)
            PD%Stats%LogFuncMode%Loc%verbose    = sum(PD%Chain%Weight(1:PD%Stats%LogFuncMode%Loc%compact-1)) + 1_IK
            close(PD%RestartFile%unit)
            close(PD%ChainFile%unit)
        endif

        if (PD%Image%isFirst) then
            write(output_unit,*)
            !call execute_command_line(" ")
            flush(output_unit)
        end if
#if defined CAF_ENABLED || defined MPI_ENABLED
        if (PD%SpecBase%ParallelizationModel%isMultiChain) then
#endif
            PD%Stats%avgCommTimePerFunCall = 0._RK
            PD%Stats%NumFunCall%acceptedRejectedDelayedUnused = PD%Stats%NumFunCall%acceptedRejectedDelayed
            PD%Stats%avgTimePerFunCalInSec =  PD%Stats%avgTimePerFunCalInSec / (PD%Stats%NumFunCall%acceptedRejectedDelayedUnused-acceptedRejectedDelayedUnusedRestartMode)
#if defined CAF_ENABLED || defined MPI_ENABLED
        elseif(PD%Image%isFirst) then
            PD%Stats%avgCommTimePerFunCall =  PD%Stats%avgCommTimePerFunCall / PD%Stats%NumFunCall%acceptedRejectedDelayed
            PD%Stats%avgTimePerFunCalInSec = (PD%Stats%avgTimePerFunCalInSec / (PD%Stats%NumFunCall%acceptedRejectedDelayedUnused-acceptedRejectedDelayedUnusedRestartMode)) * PD%Image%count
            return
        end if
#endif

        !***************************************************************************************************************************
        !***************************************************************************************************************************

    contains

        !***************************************************************************************************************************
        !***************************************************************************************************************************

        ! Objects that change state in this subroutine are: Timer, timeElapsedUntilLastReportInSeconds, sumAccRateLastReport
        subroutine reportProgress()

            use Constants_mod, only: CARRIAGE_RETURN
            implicit none

            real(RK)                            :: meanAccRateSinceStart
            real(RK)                            :: meanAccRateSinceLastReport
            real(RK)                            :: estimatedTimeToFinishInSeconds
            real(RK)                            :: timeElapsedSinceLastReportInSeconds
            integer(IK)                         :: numFunCallAccepted_dummy ! used in restart mode


            if (PD%isFreshRun) then

                call PD%Timer%toc()
                timeElapsedSinceLastReportInSeconds = PD%Timer%Time%total - timeElapsedUntilLastReportInSeconds
                timeElapsedUntilLastReportInSeconds = PD%Timer%Time%total
                meanAccRateSinceStart = SumAccRateSinceStart%acceptedRejected / real(PD%Stats%NumFunCall%acceptedRejected,kind=RK)
                meanAccRateSinceLastReport = (SumAccRateSinceStart%acceptedRejected-sumAccRateLastReport) * inverseProgressReportPeriod
                estimatedTimeToFinishInSeconds = real(PD%SpecMCMC%ChainSize%val-PD%Stats%NumFunCall%accepted,kind=RK) * timeElapsedUntilLastReportInSeconds / real(PD%Stats%NumFunCall%accepted,kind=RK)

                write( PD%TimeFile%unit, PD%TimeFile%format ) PD%Stats%NumFunCall%acceptedRejected  &
                                                            , PD%Stats%NumFunCall%accepted          &
                                                            , meanAccRateSinceStart                 &
                                                            , meanAccRateSinceLastReport            &
                                                            , timeElapsedSinceLastReportInSeconds   &
                                                            , timeElapsedUntilLastReportInSeconds   &
                                                            , estimatedTimeToFinishInSeconds
                flush(PD%TimeFile%unit)

            else

                block
                    use String_mod, only: String_type
                    type(String_type) :: Record
                    allocate( character(600) :: Record%value )
                    read(PD%TimeFile%unit, "(A)" ) Record%value
                    Record%Parts = Record%SplitStr(trim(adjustl(Record%value)),PD%SpecBase%OutputDelimiter%val,Record%nPart)
                    read(Record%Parts(1)%record,*) numFunCallAcceptedRejectedLastReport
                    read(Record%Parts(2)%record,*) numFunCallAccepted_dummy
                    read(Record%Parts(3)%record,*) meanAccRateSinceStart
                    read(Record%Parts(4)%record,*) meanAccRateSinceLastReport
                    read(Record%Parts(5)%record,*) timeElapsedSinceLastReportInSeconds
                    read(Record%Parts(6)%record,*) timeElapsedUntilLastReportInSeconds
                    read(Record%Parts(7)%record,*) estimatedTimeToFinishInSeconds
                end block
                SumAccRateSinceStart%acceptedRejected = meanAccRateSinceStart * numFunCallAcceptedRejectedLastReport

            end if

            ! report progress in the standard output
            if (PD%Image%isFirst) then
                write(output_unit   , "(A,25X,*(A25,3X))",advance="no") CARRIAGE_RETURN &
                                    , num2str(PD%Stats%NumFunCall%accepted)//" / "//num2str(PD%Stats%NumFunCall%acceptedRejected,formatIn="(1I10)") &
                                    , num2str(meanAccRateSinceLastReport,"(1F11.3)")//" / "//num2str(SumAccRateSinceStart%acceptedRejected / real(PD%Stats%NumFunCall%acceptedRejected,kind=RK),"(1F10.4)") &
                                    , num2str(timeElapsedUntilLastReportInSeconds,"(1F10.4)")//" / "//num2str(estimatedTimeToFinishInSeconds,"(1F11.3)")
                !call execute_command_line(" ")
                flush(output_unit)
            end if

            numFunCallAcceptedRejectedLastReport = PD%Stats%NumFunCall%acceptedRejected
            sumAccRateLastReport = SumAccRateSinceStart%acceptedRejected

        end subroutine reportProgress

        !***************************************************************************************************************************
        !***************************************************************************************************************************

    end subroutine runKernel

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    pure function getBurninLoc(lenLogFunc,refLogFunc,LogFunc) result(burninLoc)
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK), intent(in) :: lenLogFunc
        real(RK)   , intent(in) :: refLogFunc,LogFunc(lenLogFunc)
        real(RK)                :: negLogIncidenceProb
        integer(IK)             :: burninLoc
        negLogIncidenceProb = log( real(lenLogFunc,kind=RK) )
        burninLoc = 0_IK
        do
            burninLoc = burninLoc + 1_IK
            if ( burninLoc<lenLogFunc .and. refLogFunc-LogFunc(burninLoc)>negLogIncidenceProb ) cycle
            !burninLoc = burninLoc - 1
            exit
        end do
    end function getBurninLoc

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end submodule Kernel_smod