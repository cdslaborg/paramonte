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

#if defined PARADRAM

#define SAMPLER_TYPE ParaDRAM_type
#define SAMPLER_PROPOSAL_ABSTRACT_MOD ParaDRAMProposalAbstract_mod

#elif defined PARADISE

#define SAMPLER_TYPE ParaDISE_type
#define SAMPLER_PROPOSAL_ABSTRACT_MOD ParaDISEProposalAbstract_mod

#else

#error "Unrecognized sampler in ParaDRAM_mod.inc.f90"

#endif

    use, intrinsic :: iso_fortran_env, only: output_unit
    !use Constants_mod, only: IK, RK ! gfortran 9.3 compile crashes with this line
#if (defined MATLAB_ENABLED || defined PYTHON_ENABLED) && !defined CAF_ENABLED && !defined MPI_ENABLED
    use SAMPLER_PROPOSAL_ABSTRACT_MOD, only: ProposalErr
#endif
#if defined MPI_ENABLED
    use mpi
#endif

    implicit none

    character(*), parameter :: SUBMODULE_NAME = MODULE_NAME // "@Kernel_smod"

    type            :: SumAccRateSinceStart_type
        real(RK)    :: acceptedRejected
        real(RK)    :: acceptedRejectedDelayed
    end type SumAccRateSinceStart_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine runKernel ( self          &
                                , getLogFunc    &
                                )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: runKernel
#endif
        use ParaMonteLogFunc_mod, only: getLogFunc_proc
        use Decoration_mod, only: write
        use Constants_mod, only: RK, IK, NEGINF_RK, NLC
        use String_mod, only: num2str
        use Math_mod, only: getLogSubExp

        implicit none

        character(*), parameter             :: PROCEDURE_NAME = SUBMODULE_NAME//"@runKernel()"
        integer(IK) , parameter             :: CHAIN_RESTART_OFFSET = 2_IK

        class(SAMPLER_TYPE), intent(inout)  :: self
        procedure(getLogFunc_proc)          :: getLogFunc

#if defined CAF_ENABLED
        real(RK)    , allocatable           :: co_AccRate(:)[:]
        real(RK)    , allocatable           :: co_LogFuncState(:,:)[:]                              ! (0:nd,-1:delayedRejectionCount), -1 corresponds to the current accepted state
        integer(IK) , save                  :: co_proposalFound_samplerUpdateOccurred(2)[*]         ! merging these scalars would reduce the MPI communication overhead cost: co_proposalFound, co_samplerUpdateOccurred, co_counterDRS, 0 means false, 1 means true
#else
        real(RK)    , allocatable           :: co_LogFuncState(:,:)
        real(RK)    , allocatable           :: co_AccRate(:)
        integer(IK) , save                  :: co_proposalFound_samplerUpdateOccurred(2)            ! merging these scalars would reduce the MPI communication overhead cost: co_proposalFound, co_samplerUpdateOccurred, co_counterDRS, 0 means false, 1 means true
#endif
        type(SumAccRateSinceStart_type)     :: SumAccRateSinceStart                                 ! used to figure out the average acceptance ratio for the entire chain.
        integer(IK)                         :: numFunCallAcceptedLastAdaptation                     ! number of function calls accepted at Last proposal adaptation occurrence
        integer(IK)                         :: counterAUP                                           ! counter for adaptiveUpdatePeriod
        integer(IK)                         :: counterAUC                                           ! counter for adaptiveUpdateCount
        integer(IK)                         :: counterPRP                                           ! counter for progressReportPeriod
        integer(IK)                         :: counterDRS                                           ! counter for Delayed Rejection Stages
        integer(IK)                         :: lastStateWeight                                      ! This is used for passing the most recent verbose chain segment to the adaptive updater of the sampler
        integer(IK)                         :: currentStateWeight                                   ! counter for SampleWeight, used only in in restart mode
        integer(IK)                         :: numFunCallAcceptedPlusOne                            ! counter for SampleWeight, used only in in restart mode
        integer(IK)                         :: numFunCallAcceptedRejectedLastReport                 ! used for progress-report: must be initialized to zero upon entry to the procedure
        real(RK)                            :: timeElapsedUntilLastReportInSeconds                  ! used for progress-report: must be initialized to zero upon entry to the procedure
        real(RK)                            :: inverseProgressReportPeriod                          ! used for progress-report: inverse of progressReportPeriod
        real(RK)                            :: sumAccRateLastReport                                 ! used for progress-report: must be initialized to zero upon entry to the procedure
        real(RK)                            :: uniformRnd                                           ! used for random number generation.
        real(RK)                            :: meanAccRateSinceStart                                ! used for restart file read
       !real(RK)                            :: adaptationMeasureDummy
        real(RK)                            :: maxLogFuncRejectedProposal
        logical                             :: samplerUpdateIsGreedy
        logical                             :: samplerUpdateSucceeded
        logical                             :: delayedRejectionRequested
        logical                             :: noDelayedRejectionRequested
        real(RK)    , allocatable           :: AdaptationMeasure(:), adaptationMeasurePlaceHolder(:)
        integer(IK)                         :: adaptationMeasureLen
        integer(IK)                         :: acceptedRejectedDelayedUnusedRestartMode
        integer(IK)                         :: imageID, dummy
        integer(IK)                         :: nd
        integer(IK), parameter              :: STDOUT_SEGLEN = 25
        character(4*STDOUT_SEGLEN+2*3+1)    :: txt
#if defined CAF_ENABLED || defined MPI_ENABLED
        integer(IK)                         :: imageStartID, imageEndID
#if defined CAF_ENABLED
        logical     , save                  :: co_proposalFoundSinglChainMode[*]                    ! used in the delayed rejection section
#elif defined MPI_ENABLED
        logical                             :: co_proposalFoundSinglChainMode                       ! used in the delayed rejection section
        real(RK)    , allocatable           :: AccRateMatrix(:,:)                                   ! matrix of size (-1:self%SpecDRAM%DelayedRejectionCount%val,1:self%Image%count)
        integer(IK)                         :: ndPlusOne
        integer(IK)                         :: ierrMPI
        integer(IK)                         :: delayedRejectionCountPlusTwo
#endif
        self%Stats%avgCommTimePerFunCall = 0._RK                                                    ! Until the reporting time, this is in reality, sumCommTimePerFunCall. This is meaningful only in singlChain parallelism
        self%Stats%NumFunCall%acceptedRejectedDelayedUnused = self%Image%count                      ! used only in singlChain parallelism, and relevant only on the first image
#endif
        acceptedRejectedDelayedUnusedRestartMode = 0_IK                                             ! used to compute more accurate timings in the restart mode
        self%Stats%avgTimePerFunCalInSec = 0._RK
        numFunCallAcceptedRejectedLastReport = 0_IK
        timeElapsedUntilLastReportInSeconds = 0._RK
        inverseProgressReportPeriod = 1._RK/real(self%SpecBase%ProgressReportPeriod%val,kind=RK)    ! this remains a constant except for the last the last report of the simulation
        sumAccRateLastReport = 0._RK
        nd = self%nd%val
        
        if (allocated(co_AccRate)) deallocate(co_AccRate)
        if (allocated(co_LogFuncState)) deallocate(co_LogFuncState)
#if defined CAF_ENABLED
        allocate(co_LogFuncState(0:nd,-1:self%SpecDRAM%DelayedRejectionCount%val)[*])
        allocate(co_AccRate(-1:self%SpecDRAM%DelayedRejectionCount%val)[*])                         ! the negative element will contain counterDRS
#else
        allocate(co_LogFuncState(0:nd,-1:self%SpecDRAM%DelayedRejectionCount%val))
        allocate(co_AccRate(-1:self%SpecDRAM%DelayedRejectionCount%val))
#endif
        co_AccRate(-1)  = 0._RK                                                                     ! the real-value counterDRS, indicating the initial delayed rejection stage at which the first point is sampled
        co_AccRate(0)   = 1._RK                                                                     ! initial acceptance rate for the first zeroth DR stage.
        co_AccRate(1:self%SpecDRAM%DelayedRejectionCount%val) = 0._RK                               ! indicates the very first proposal acceptance on image 1

#if defined MPI_ENABLED
        if (allocated(AccRateMatrix)) deallocate(AccRateMatrix)
        allocate(AccRateMatrix(-1:self%SpecDRAM%DelayedRejectionCount%val,1:self%Image%count))      ! the negative element will contain counterDRS
        AccRateMatrix = 0._RK                                                                       ! -huge(1._RK)  ! debug
        AccRateMatrix(0,1:self%Image%count) = 1._RK                                                 ! initial acceptance rate for the first zeroth DR stage.
        ndPlusOne = nd + 1_IK
        delayedRejectionCountPlusTwo = self%SpecDRAM%DelayedRejectionCount%val + 2_IK
#endif

        delayedRejectionRequested                       = self%SpecDRAM%DelayedRejectionCount%val > 0_IK
        noDelayedRejectionRequested                     = .not. delayedRejectionRequested
        if (delayedRejectionRequested) then
        self%Stats%NumFunCall%acceptedRejectedDelayed   = 0_IK                                      ! Markov Chain counter
        SumAccRateSinceStart%acceptedRejectedDelayed    = 0._RK                                     ! sum of acceptance rate
        end if

       !adaptationMeasure                               = 0._RK                                     ! needed for the first output
        SumAccRateSinceStart%acceptedRejected           = 0._RK                                     ! sum of acceptance rate
        self%Stats%NumFunCall%acceptedRejected          = 0_IK                                      ! Markov Chain counter
        counterAUC                                      = 0_IK                                      ! counter for padaptiveUpdateCount.
        counterPRP                                      = 0_IK                                      ! counter for progressReportPeriod.
        counterAUP                                      = 0_IK                                      ! counter for adaptiveUpdatePeriod. 
        self%Stats%NumFunCall%accepted                  = 0_IK                                      ! Markov Chain acceptance counter.
        samplerUpdateSucceeded                          = .true.                                    ! needed to set up lastStateWeight and numFunCallAcceptedLastAdaptation for the first accepted proposal
        numFunCallAcceptedLastAdaptation                = 0_IK
        lastStateWeight                                 = -huge(lastStateWeight)
        meanAccRateSinceStart                           = 1._RK ! needed for the first restart output in fresh run.

        call self%Timer%tic()

        blockDryRunSetup: if (self%isFreshRun) then

            adaptationMeasureLen = 100_IK

            allocate(self%Chain%ProcessID   (   self%SpecMCMC%ChainSize%val))
            allocate(self%Chain%DelRejStage (   self%SpecMCMC%ChainSize%val))
            allocate(self%Chain%MeanAccRate (   self%SpecMCMC%ChainSize%val))
            allocate(self%Chain%Adaptation  (   self%SpecMCMC%ChainSize%val))
            allocate(self%Chain%BurninLoc   (   self%SpecMCMC%ChainSize%val))
            allocate(self%Chain%Weight      (   self%SpecMCMC%ChainSize%val))
            allocate(self%Chain%State       (nd,self%SpecMCMC%ChainSize%val))
            allocate(self%Chain%LogFunc     (   self%SpecMCMC%ChainSize%val))

        else blockDryRunSetup

            ! load the existing Chain file into self%Chain components

            call self%Chain%get ( chainFilePath = self%ChainFile%Path%original      &
                                , chainFileForm = self%SpecBase%ChainFileFormat%val &
                                , Err = self%Err                                    &
                                , targetChainSize = self%SpecMCMC%ChainSize%val     &
                                , lenHeader = self%Chain%lenHeader                  &
                                , ndim = nd                                         &
                                , delimiter = self%SpecBase%OutputDelimiter%val     &
                                )
            if (self%Err%occurred) then
                self%Err%msg = PROCEDURE_NAME//self%Err%msg
                call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                return
            end if

            if (self%Chain%Count%compact<=CHAIN_RESTART_OFFSET) then
                self%isFreshRun = .true.
                self%isDryRun = .not. self%isFreshRun
            end if

#if defined CAF_ENABLED
            sync all
#elif defined MPI_ENABLED
            call mpi_barrier(mpi_comm_world,ierrMPI)
#endif

            blockMasterSetup: if (self%Image%isMaster) then

                ! set up the chain file

                block

                    use System_mod, only: RandomFileName_type, copyFile, removeFile
                    type(RandomFileName_type) :: RFN

                    ! create a copy of the chain file, just for the sake of not losing the simulation results

                    RFN = RandomFileName_type(dir = "", key = self%ChainFile%Path%original//".temporary_restart_copy", ext="") 
                    call copyFile(pathOld=self%ChainFile%Path%original,pathNew=RFN%path,isWindows=self%OS%isWindows,Err=self%err)
                    if (self%Err%occurred) then
                        self%Err%msg = PROCEDURE_NAME//self%Err%msg
                        call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                        return
                    end if

                    ! reopen the chain file to resume the simulation

                    open( newunit = self%ChainFile%unit             &
                        , file = self%ChainFile%Path%original       &
                        , form = self%ChainFile%Form%value          &
                        , status = self%ChainFile%status            &
                        , iostat = self%ChainFile%Err%stat          &
#if defined IFORT_ENABLED && defined OS_IS_WINDOWS
                        , SHARED                                    &
#endif
                        , position = self%ChainFile%Position%value  )
                    self%Err = self%ChainFile%getOpenErr(self%ChainFile%Err%stat)
                    if (self%Err%occurred) then
                        self%Err%msg = PROCEDURE_NAME//": Error occurred while opening "//self%name//" "//self%ChainFile%suffix//" file='"//self%ChainFile%Path%original//"'."
                        call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                        return
                    end if

                    ! rewrite the chain file

                    call self%Chain%writeChainFile  ( ndim = nd &
                                                    , compactStartIndex = 1_IK &
                                                    , compactEndIndex = self%Chain%Count%compact-CHAIN_RESTART_OFFSET &
                                                    , chainFileUnit = self%ChainFile%unit &
                                                    , chainFileForm = self%SpecBase%ChainFileFormat%val &
                                                    , chainFileFormat = self%ChainFile%format &
                                                    , adaptiveUpdatePeriod = self%SpecDRAM%AdaptiveUpdatePeriod%val &
                                                    )

                    ! remove the temporary copy of the chain file

                    call removeFile(path=RFN%path,isWindows=self%OS%isWindows,Err=self%Err)
                    if (self%Err%occurred) then
                        self%Err%msg = PROCEDURE_NAME//self%Err%msg
                        call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                        return
                    end if

                end block

                adaptationMeasureLen = maxval(self%Chain%Weight(1:self%Chain%Count%compact-CHAIN_RESTART_OFFSET))

            end if blockMasterSetup

        end if blockDryRunSetup

        if (self%Image%isMaster) then
            if (allocated(AdaptationMeasure)) deallocate(AdaptationMeasure); allocate(AdaptationMeasure(adaptationMeasureLen))
        end if

        if (self%isFreshRun) then ! this must be done separately from the above blockDryRunSetup

            self%Chain%BurninLoc(1)  = 1_IK

            co_LogFuncState(1:nd,0) = self%SpecMCMC%StartPointVec%Val   ! proposal state
            call self%Timer%toc()
            co_LogFuncState(0,0) = getLogFunc( nd, co_LogFuncState(1:nd,0) )    ! proposal logFunc
            call self%Timer%toc(); self%Stats%avgTimePerFunCalInSec = self%Stats%avgTimePerFunCalInSec + self%Timer%Time%delta

        else

            co_LogFuncState(1:nd,0)     = self%Chain%State(1:nd,1)      ! proposal logFunc
            co_LogFuncState(0,0)        = self%Chain%LogFunc(1)         ! proposal logFunc

        end if

        co_LogFuncState(0:nd,-1) = co_LogFuncState(0:nd,0)  ! set current logFunc and State equal to the first proposal
        self%Stats%LogFuncMode%val = -huge(self%Stats%LogFuncMode%val)
        self%Stats%LogFuncMode%Loc%compact = 0_IK

        if (self%Image%isFirst) then
            txt =   repeat(" ",STDOUT_SEGLEN) &
                //  "Accepted/Total Func. Call   " &
                //  "Dynamic/Overall Acc. Rate   " &
                //  "Elapsed/Remained Time [s] "
            call write(string=txt)
            txt =   repeat(" ",STDOUT_SEGLEN) &
                //  repeat("=",STDOUT_SEGLEN) // "   " &
                //  repeat("=",STDOUT_SEGLEN) // "   " &
                //  repeat("=",STDOUT_SEGLEN) // " "
            call write(string=txt)
            !call execute_command_line(" ")
            flush(output_unit)
        end if

#if defined CAF_ENABLED || defined MPI_ENABLED
        if (self%SpecBase%ParallelizationModel%isSinglChain) then
            imageStartID = 1_IK
            imageEndID = self%Image%count
        else ! isMultiChain
            imageStartID = self%Image%id
            imageEndID = self%Image%id
        end if
#else
        imageID = 1_IK ! needed even in the case of serial run to assign a proper value to self%Chain%ProcessID
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start of loopMarkovChain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        loopMarkovChain: do

            co_proposalFound_samplerUpdateOccurred(2) = 0_IK ! at each iteration assume no samplerUpdateOccurred, unless it occurs

#if defined CAF_ENABLED || defined MPI_ENABLED
            blockMasterImage: if (self%Image%isMaster) then
#endif

                co_proposalFound_samplerUpdateOccurred(1) = 0_IK ! co_proposalFound = .false.
                samplerUpdateIsGreedy = counterAUC < self%SpecDRAM%GreedyAdaptationCount%val

#if defined CAF_ENABLED || defined MPI_ENABLED
                loopOverImages: do imageID = imageStartID, imageEndID
#if defined CAF_ENABLED
                    call self%Timer%toc()
                    if (imageID/=self%Image%id) co_AccRate(-1:self%SpecDRAM%DelayedRejectionCount%val) = co_AccRate(-1:self%SpecDRAM%DelayedRejectionCount%val)[imageID] ! happens only in isSinglChain=TRUE
                    call self%Timer%toc(); self%Stats%avgCommTimePerFunCall = self%Stats%avgCommTimePerFunCall + self%Timer%Time%delta
#elif defined MPI_ENABLED
                    if (imageID/=self%Image%id) co_AccRate(-1:self%SpecDRAM%DelayedRejectionCount%val) = AccRateMatrix(-1:self%SpecDRAM%DelayedRejectionCount%val,imageID) ! happens only in isSinglChain=TRUE
#endif
#endif

                    counterDRS = nint(co_AccRate(-1),kind=IK)
                    if (counterDRS > -1_IK) co_proposalFound_samplerUpdateOccurred(1) = 1_IK ! co_proposalFound = .true.

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% blockProposalAccepted %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    ! On the very first iteration, this block is (and must be) executed for imageID==1,
                    ! since it is for the first (starting) point, which is assumed to have been accepted
                    ! as the first point by the first coarray imageID.

                    blockProposalAccepted: if ( co_proposalFound_samplerUpdateOccurred(1) == 1_IK ) then ! co_proposalAccepted = .true.

                        currentStateWeight = 0_IK

                        ! communicate the accepted logFunc and State from the winning image to master/all images: co_LogFuncState

#if defined MPI_ENABLED
                        if (self%SpecBase%ParallelizationModel%isSinglChain) then
                            call self%Timer%toc()
                            ! broadcast winning image to all processes
                            call mpi_bcast  ( imageID           &   ! buffer
                                            , 1                 &   ! count
                                            , mpi_integer       &   ! datatype
                                            , 0                 &   ! root: broadcasting rank
                                            , mpi_comm_world    &   ! comm
                                            , ierrMPI           &   ! ierr
                                            )
                            ! broadcast co_LogFuncState from the winning image to all others
                            call mpi_bcast  ( co_LogFuncState(0:nd,-1)  &   ! buffer
                                            , ndPlusOne                 &   ! count
                                            , mpi_double_precision      &   ! datatype
                                            , imageID - 1_IK            &   ! root: broadcasting rank
                                            , mpi_comm_world            &   ! comm
                                            , ierrMPI                   &   ! ierr
                                            )
                            call self%Timer%toc(); self%Stats%avgCommTimePerFunCall = self%Stats%avgCommTimePerFunCall + self%Timer%Time%delta
                        end if
#elif defined CAF_ENABLED
                        if (imageID/=self%Image%id) then ! Avoid remote connection for something that is available locally.
                            call self%Timer%toc()
                            co_LogFuncState(0:nd,-1) = co_LogFuncState(0:nd,-1)[imageID]
                            call self%Timer%toc(); self%Stats%avgCommTimePerFunCall = self%Stats%avgCommTimePerFunCall + self%Timer%Time%delta
                        end if
#endif

                        ! Note: after every adaptive update of the sampler, counterAUP is reset to 0.
                        if (counterAUP==0_IK .and. samplerUpdateSucceeded) then
                            numFunCallAcceptedLastAdaptation = numFunCallAcceptedLastAdaptation + 1_IK
                            lastStateWeight = 0_IK
                        end if

                        blockFreshDryRun: if (self%isFreshRun) then

                            call writeOutput()

                            self%Stats%NumFunCall%accepted = self%Stats%NumFunCall%accepted + 1_IK

                            self%Chain%ProcessID(self%Stats%NumFunCall%accepted)    = imageID
                            self%Chain%DelRejStage(self%Stats%NumFunCall%accepted)  = counterDRS
                            self%Chain%Adaptation(self%Stats%NumFunCall%accepted)   = 0._RK ! adaptationMeasure
                            self%Chain%Weight(self%Stats%NumFunCall%accepted)       = 0_IK
                            self%Chain%LogFunc(self%Stats%NumFunCall%accepted)      = co_LogFuncState(0,-1)
                            self%Chain%State(1:nd,self%Stats%NumFunCall%accepted)   = co_LogFuncState(1:nd,-1)

                            ! find the burnin point

                            self%Chain%BurninLoc(self%Stats%NumFunCall%accepted) = getBurninLoc ( lenLogFunc    = self%Stats%NumFunCall%accepted                        &
                                                                                                , refLogFunc    = self%Stats%LogFuncMode%val                            &
                                                                                                , LogFunc       = self%Chain%LogFunc(1:self%Stats%NumFunCall%accepted)  &
                                                                                                )

                        else blockFreshDryRun ! in restart mode: determine the correct value of co_proposalFound_samplerUpdateOccurred(1)

                            numFunCallAcceptedPlusOne = self%Stats%NumFunCall%accepted + 1_IK
                            if (numFunCallAcceptedPlusOne==self%Chain%Count%compact) then
                                self%isFreshRun = .true.
                                call writeOutput()
                                self%isDryRun = .not. self%isFreshRun
                                self%Chain%Weight(numFunCallAcceptedPlusOne) = 0_IK
                                SumAccRateSinceStart%acceptedRejected = self%Chain%MeanAccRate(self%Stats%NumFunCall%accepted) * real(self%Stats%NumFunCall%acceptedRejected,kind=RK)
                                if (delayedRejectionRequested) SumAccRateSinceStart%acceptedRejectedDelayed = self%Chain%MeanAccRate(self%Stats%NumFunCall%accepted) * real(self%Stats%NumFunCall%acceptedRejectedDelayed,kind=RK)
                            end if
                            self%Stats%NumFunCall%accepted = numFunCallAcceptedPlusOne
                            numFunCallAcceptedPlusOne = self%Stats%NumFunCall%accepted + 1_IK

                        end if blockFreshDryRun

                        if (self%Stats%LogFuncMode%val < self%Chain%LogFunc(self%Stats%NumFunCall%accepted)) then
                            self%Stats%LogFuncMode%val = self%Chain%LogFunc(self%Stats%NumFunCall%accepted)
                            self%Stats%LogFuncMode%Loc%compact = self%Stats%NumFunCall%accepted
                        end if

                        SumAccRateSinceStart%acceptedRejected = SumAccRateSinceStart%acceptedRejected + co_AccRate(counterDRS)

                    else blockProposalAccepted

                        counterDRS = self%SpecDRAM%DelayedRejectionCount%val
                        SumAccRateSinceStart%acceptedRejected = SumAccRateSinceStart%acceptedRejected + co_AccRate(counterDRS)

                    end if blockProposalAccepted

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% blockProposalAccepted %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    counterAUP = counterAUP + 1_IK
                    counterPRP = counterPRP + 1_IK
                    currentStateWeight = currentStateWeight + 1_IK
                    self%Stats%NumFunCall%acceptedRejected = self%Stats%NumFunCall%acceptedRejected + 1_IK

                    if (delayedRejectionRequested) then
                        SumAccRateSinceStart%acceptedRejectedDelayed = SumAccRateSinceStart%acceptedRejectedDelayed + sum(co_AccRate(0:counterDRS))
                        self%Stats%NumFunCall%acceptedRejectedDelayed = self%Stats%NumFunCall%acceptedRejectedDelayed + counterDRS + 1_IK
                    end if

                    if (self%isFreshRun) then ! these are used for adaptive proposal updating, so they have to be set on every accepted or rejected iteration (excluding delayed rejections)
                        self%Chain%MeanAccRate(self%Stats%NumFunCall%accepted)  = SumAccRateSinceStart%acceptedRejected / real(self%Stats%NumFunCall%acceptedRejected,kind=RK)
                        self%Chain%Weight(self%Stats%NumFunCall%accepted)       = self%Chain%Weight(self%Stats%NumFunCall%accepted) + 1_IK
                        if (adaptationMeasureLen<self%Chain%Weight(self%Stats%NumFunCall%accepted)) then
                            allocate(adaptationMeasurePlaceHolder(2*adaptationMeasureLen))
                            adaptationMeasurePlaceHolder(1:adaptationMeasureLen) = AdaptationMeasure
                            call move_alloc(from=adaptationMeasurePlaceHolder,to=AdaptationMeasure)
                            adaptationMeasureLen = 2_IK * adaptationMeasureLen
                        end if
                    else
#if defined CAF_ENABLED || defined MPI_ENABLED
                        acceptedRejectedDelayedUnusedRestartMode = self%Stats%NumFunCall%acceptedRejectedDelayedUnused
#else
                        acceptedRejectedDelayedUnusedRestartMode = self%Stats%NumFunCall%acceptedRejectedDelayed
#endif
                    end if

                    if (counterPRP == self%SpecBase%ProgressReportPeriod%val) then
                        counterPRP = 0_IK
                        call reportProgress()
                    end if

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% last output write %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    ! in paradise mode, it is imperative to finish the simulation before any further redundant sampler updates.
                    blockLastSample: if (self%Stats%NumFunCall%accepted==self%SpecMCMC%ChainSize%val) then !co_missionAccomplished = .true.

                        ! on 3 images Windows, substituting co_missionAccomplished with the following leads to 10% less communication overhead for 1D Gaussian example
                        ! on 3 images Linux  , substituting co_missionAccomplished with the following leads to 16% less communication overhead for 1D Gaussian example
                        ! on 5 images Linux  , substituting co_missionAccomplished with the following leads to 11% less communication overhead for 1D Gaussian example

                        co_proposalFound_samplerUpdateOccurred(1) = -1_IK  ! equivalent to co_missionAccomplished = .true.
                        inverseProgressReportPeriod = 1._RK / (self%Stats%NumFunCall%acceptedRejected-numFunCallAcceptedRejectedLastReport)

                        if (self%isFreshRun) then
                            call writeOutput()
                            flush(self%ChainFile%unit)
                        end if

                        call reportProgress()

#if defined CAF_ENABLED || defined MPI_ENABLED
                        exit loopOverImages
#endif
                    end if blockLastSample

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% last output write %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Proposal Adaptation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    blockSamplerAdaptation: if ( counterAUC < self%SpecDRAM%AdaptiveUpdateCount%val .and. counterAUP == self%SpecDRAM%AdaptiveUpdatePeriod%val ) then

                        co_proposalFound_samplerUpdateOccurred(2) = 1_IK ! istart = numFunCallAcceptedLastAdaptation ! = max( numFunCallAcceptedLastAdaptation , self%Chain%BurninLoc(self%Stats%NumFunCall%accepted) ) ! this is experimental

                        ! the order in the following two MUST be preserved as occasionally self%Stats%NumFunCall%accepted = numFunCallAcceptedLastAdaptation

                        dummy = self%Chain%Weight(self%Stats%NumFunCall%accepted) ! needed for the restart mode, not needed in the fresh run
                        if (self%Stats%NumFunCall%accepted==numFunCallAcceptedLastAdaptation) then    ! no new point has been accepted since last time
                            self%Chain%Weight(numFunCallAcceptedLastAdaptation) = currentStateWeight - lastStateWeight
#if defined DBG_ENABLED && !defined CAF_ENABLED && !defined MPI_ENABLED
                            if (mod(self%Chain%Weight(numFunCallAcceptedLastAdaptation),self%SpecDRAM%AdaptiveUpdatePeriod%val)/=0) then
                                write(output_unit,"(*(g0,:,' '))"   ) PROCEDURE_NAME//": Internal error occurred: ", self%SpecDRAM%AdaptiveUpdatePeriod%val &
                                                                    , self%Chain%Weight(numFunCallAcceptedLastAdaptation), currentStateWeight, lastStateWeight
                                !call execute_command_line(" ")
                                flush(output_unit)
                                error stop
                            end if
#endif
                        else
                            self%Chain%Weight(numFunCallAcceptedLastAdaptation) = self%Chain%Weight(numFunCallAcceptedLastAdaptation) - lastStateWeight
                            self%Chain%Weight(self%Stats%NumFunCall%accepted) = currentStateWeight ! needed for the restart mode, not needed in the fresh run
                        end if

                        meanAccRateSinceStart = self%Chain%MeanAccRate(self%Stats%NumFunCall%accepted) ! used only in fresh run, but not worth putting it in a conditional block.
                        call self%Proposal%doAdaptation ( nd                        = nd                                                                                        &
                                                        , chainSize                 = self%Stats%NumFunCall%accepted - numFunCallAcceptedLastAdaptation + 1_IK                  &
                                                        , Chain                     = self%Chain%State(1:nd,numFunCallAcceptedLastAdaptation:self%Stats%NumFunCall%accepted)    &
                                                        , ChainWeight               = self%Chain%Weight(numFunCallAcceptedLastAdaptation:self%Stats%NumFunCall%accepted)        &
                                                        , isFreshRun                = self%isFreshRun                                                                           &
                                                        , samplerUpdateIsGreedy     = samplerUpdateIsGreedy                                                                     &
                                                        , meanAccRateSinceStart     = meanAccRateSinceStart                                                                     &
                                                        , samplerUpdateSucceeded    = samplerUpdateSucceeded                                                                    &
                                                        , adaptationMeasure         = AdaptationMeasure(dummy)                                                                  &
                                                        )
#if (defined MATLAB_ENABLED || defined PYTHON_ENABLED) && !defined CAF_ENABLED && !defined MPI_ENABLED
                        if(ProposalErr%occurred) then; self%Err%occurred = .true.; return; end if
#endif
                        if (self%isDryRun) SumAccRateSinceStart%acceptedRejected = meanAccRateSinceStart * real(self%Stats%NumFunCall%acceptedRejected,kind=RK)

                        self%Chain%Weight(self%Stats%NumFunCall%accepted) = dummy   ! needed for the restart mode, not needed in the fresh run
                        if (self%Stats%NumFunCall%accepted==numFunCallAcceptedLastAdaptation) then
                            !adaptationMeasure = adaptationMeasure + adaptationMeasureDummy ! this is the worst-case upper-bound
                            self%Chain%Adaptation(self%Stats%NumFunCall%accepted) = min(1._RK, self%Chain%Adaptation(self%Stats%NumFunCall%accepted) + AdaptationMeasure(dummy)) ! this is the worst-case upper-bound
                        else
                            !adaptationMeasure = adaptationMeasureDummy
                            self%Chain%Adaptation(self%Stats%NumFunCall%accepted) = AdaptationMeasure(dummy)
                            self%Chain%Weight(numFunCallAcceptedLastAdaptation) = self%Chain%Weight(numFunCallAcceptedLastAdaptation) + lastStateWeight
                        end if
                        if (samplerUpdateSucceeded) then
                            lastStateWeight = currentStateWeight ! self%Chain%Weight(self%Stats%NumFunCall%accepted) ! informative, do not remove
                            numFunCallAcceptedLastAdaptation = self%Stats%NumFunCall%accepted
                        end if

                        counterAUP = 0_IK
                        counterAUC = counterAUC + 1_IK
                        !if (counterAUC==self%SpecDRAM%AdaptiveUpdateCount%val) adaptationMeasure = 0._RK

                    else blockSamplerAdaptation

                        AdaptationMeasure(self%Chain%Weight(self%Stats%NumFunCall%accepted)) = 0._RK

                    end if blockSamplerAdaptation

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Proposal Adaptation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined CAF_ENABLED || defined MPI_ENABLED
                    if (co_proposalFound_samplerUpdateOccurred(1)==1_IK) exit loopOverImages

                end do loopOverImages
#endif

#if defined MPI_ENABLED
                if (self%SpecBase%ParallelizationModel%isSinglChain .and. co_proposalFound_samplerUpdateOccurred(1)==0_IK) then
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

#if defined CAF_ENABLED

                if (self%SpecBase%ParallelizationModel%isSinglChain) then
                    call self%Timer%toc()
                    sync images(*)
                    call self%Timer%toc(); self%Stats%avgCommTimePerFunCall = self%Stats%avgCommTimePerFunCall + self%Timer%Time%delta
                end if

            else blockMasterImage   ! ATTN: This block should be executed only when singlChain parallelizationModel is requested

                sync images(1)

                ! get the accepted proposal from the first image

                call self%Timer%toc()
                co_proposalFound_samplerUpdateOccurred(1:2) = co_proposalFound_samplerUpdateOccurred(1:2)[1]
                if (co_proposalFound_samplerUpdateOccurred(1)==1_IK) co_LogFuncState(0:nd,-1) = co_LogFuncState(0:nd,-1)[1]
                if (co_proposalFound_samplerUpdateOccurred(2)==1_IK) call self%Proposal%bcastAdaptation()
                call self%Timer%toc(); self%Stats%avgCommTimePerFunCall = self%Stats%avgCommTimePerFunCall + self%Timer%Time%delta

                if (self%isDryRun) then
                    if (co_proposalFound_samplerUpdateOccurred(1)==1_IK) then
                        self%Stats%NumFunCall%accepted = self%Stats%NumFunCall%accepted + 1_IK
                        numFunCallAcceptedPlusOne = self%Stats%NumFunCall%accepted + 1_IK
                        currentStateWeight = 1_IK
                        if (self%Stats%NumFunCall%accepted==self%Chain%Count%compact) then
                            self%isFreshRun = .true.
                            self%isDryRun = .false.
                        end if
                    else
                        currentStateWeight = currentStateWeight + self%Image%count
                    end if
#if defined DBG_ENABLED
                elseif (co_proposalFound_samplerUpdateOccurred(1)==1_IK) then
                    self%Stats%NumFunCall%accepted = self%Stats%NumFunCall%accepted + 1_IK
#endif
                end if

            end if blockMasterImage

#elif defined MPI_ENABLED

            else blockMasterImage   ! This block should be executed only when singlChain parallelizationModel is requested

                ! fetch winning rank to all processes

                call self%Timer%toc()
                call mpi_bcast  ( imageID           &   ! buffer
                                , 1                 &   ! count
                                , mpi_integer       &   ! datatype
                                , 0                 &   ! root: broadcasting rank
                                , mpi_comm_world    &   ! comm
                                , ierrMPI           &   ! ierr
                                )
                call self%Timer%toc(); self%Stats%avgCommTimePerFunCall = self%Stats%avgCommTimePerFunCall + self%Timer%Time%delta

                if (imageID>0_IK) then ! co_ProposalFound = .true., sampling successful

                    ! broadcast co_LogFuncState from the winning image to all others

                    call self%Timer%toc()
                    call mpi_bcast  ( co_LogFuncState(0:nd,-1)  &   ! buffer
                                    , ndPlusOne                 &   ! count
                                    , mpi_double_precision      &   ! datatype
                                    , imageID - 1               &   ! root: broadcasting rank
                                    , mpi_comm_world            &   ! comm
                                    , ierrMPI                   &   ! ierr
                                    )
                    call self%Timer%toc(); self%Stats%avgCommTimePerFunCall = self%Stats%avgCommTimePerFunCall + self%Timer%Time%delta

                end if

                if (self%isDryRun) then
                    if (imageID>0_IK) then ! equivalent to co_proposalFound_samplerUpdateOccurred(1)==1_IK
                        self%Stats%NumFunCall%accepted = self%Stats%NumFunCall%accepted + 1_IK
                        numFunCallAcceptedPlusOne = self%Stats%NumFunCall%accepted + 1_IK
                        currentStateWeight = 1_IK
                        if (self%Stats%NumFunCall%accepted==self%Chain%Count%compact) then
                            self%isFreshRun = .true.
                            self%isDryRun = .false.
                        end if
                    else
                        currentStateWeight = currentStateWeight + self%Image%count
                    end if
#if defined DBG_ENABLED
                elseif (co_proposalFound_samplerUpdateOccurred(1)==1_IK) then
                    self%Stats%NumFunCall%accepted = self%Stats%NumFunCall%accepted + 1_IK
#endif
                end if

            end if blockMasterImage

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% begin Common block between all images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ! broadcast samplerUpdateOccurred from the root process to all others, and broadcast proposal adaptations if needed

            if (self%SpecBase%ParallelizationModel%isSinglChain) then
                call self%Timer%toc()
                call mpi_bcast  ( co_proposalFound_samplerUpdateOccurred  &   ! buffer: XXX: first element not needed except to end the simulation.
                                , 2                     &   ! count
                                , mpi_integer           &   ! datatype
                                , 0                     &   ! root: broadcasting rank
                                , mpi_comm_world        &   ! comm
                                , ierrMPI               &   ! ierr
                                )
                if (co_proposalFound_samplerUpdateOccurred(2)==1_IK) call self%Proposal%bcastAdaptation()
                call self%Timer%toc(); self%Stats%avgCommTimePerFunCall = self%Stats%avgCommTimePerFunCall + self%Timer%Time%delta
            end if

#endif

            if (co_proposalFound_samplerUpdateOccurred(1) == -1_IK) exit loopMarkovChain   ! we are done: co_missionAccomplished = .true.

            co_AccRate(-1) = -1._RK ! counterDRS at which new proposal is accepted. this is essential for all serial and parallel modes
            maxLogFuncRejectedProposal = NEGINF_RK

#if defined CAF_ENABLED || defined MPI_ENABLED
            blockSingleChainParallelism: if (self%SpecBase%ParallelizationModel%isSinglChain) then
#define SINGLCHAIN_PARALLELISM
#include "ParaDRAM_mod@Kernel_smod@nextMove.inc.f90"
#undef SINGLCHAIN_PARALLELISM
            else blockSingleChainParallelism
#endif
#include "ParaDRAM_mod@Kernel_smod@nextMove.inc.f90"
#if defined CAF_ENABLED || defined MPI_ENABLED
            end if blockSingleChainParallelism
#endif

            cycle loopMarkovChain

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end Common block between all images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do loopMarkovChain

        call self%Timer%toc()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%* end of loopMarkovChain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        self%chain%Count%target  = self%Stats%NumFunCall%accepted
        self%chain%Count%compact = self%Stats%NumFunCall%accepted
        self%chain%Count%verbose = self%Stats%NumFunCall%acceptedRejected

        if (noDelayedRejectionRequested) then
            self%Stats%NumFunCall%acceptedRejectedDelayed = self%Stats%NumFunCall%acceptedRejected
            SumAccRateSinceStart%acceptedRejectedDelayed = SumAccRateSinceStart%acceptedRejected
        end if
            

        if (self%Image%isMaster) then
            block
                integer(IK) :: i
                if (self%SpecDRAM%BurninAdaptationMeasure%val>0.999999_RK) then
                    self%Stats%AdaptationBurninLoc%compact = 1_IK
                    self%Stats%AdaptationBurninLoc%verbose = 1_IK
                else
                    self%Stats%AdaptationBurninLoc%compact = self%Stats%NumFunCall%accepted
                    self%Stats%AdaptationBurninLoc%verbose = self%Stats%NumFunCall%acceptedRejected - self%Chain%Weight(self%Stats%NumFunCall%accepted) + 1_IK
                    loopAdaptationBurnin: do i = self%Stats%NumFunCall%accepted-1, 1, -1
                        if (self%Chain%Adaptation(i)>self%SpecDRAM%BurninAdaptationMeasure%val) exit loopAdaptationBurnin
                        self%Stats%AdaptationBurninLoc%compact = self%Stats%AdaptationBurninLoc%compact - 1_IK
                        self%Stats%AdaptationBurninLoc%verbose = self%Stats%AdaptationBurninLoc%verbose - self%Chain%Weight(i)
                    end do loopAdaptationBurnin
                end if
            end block
            self%Stats%BurninLoc%compact        = self%Chain%BurninLoc(self%Stats%NumFunCall%accepted)
            self%Stats%BurninLoc%verbose        = sum(self%Chain%Weight(1:self%Stats%BurninLoc%compact-1)) + 1_IK
            self%Stats%LogFuncMode%Crd          = self%Chain%State(1:nd,self%Stats%LogFuncMode%Loc%compact)
            self%Stats%LogFuncMode%Loc%verbose  = sum(self%Chain%Weight(1:self%Stats%LogFuncMode%Loc%compact-1)) + 1_IK
            close(self%RestartFile%unit)
            close(self%ChainFile%unit)
        endif

        if (self%Image%isFirst) then
            call write()
            !call execute_command_line(" ")
            flush(output_unit)
        end if
#if defined CAF_ENABLED || defined MPI_ENABLED
        if (self%SpecBase%ParallelizationModel%isMultiChain) then
#endif
            self%Stats%avgCommTimePerFunCall = 0._RK
            self%Stats%NumFunCall%acceptedRejectedDelayedUnused = self%Stats%NumFunCall%acceptedRejectedDelayed
            self%Stats%avgTimePerFunCalInSec =  self%Stats%avgTimePerFunCalInSec / (self%Stats%NumFunCall%acceptedRejectedDelayedUnused-acceptedRejectedDelayedUnusedRestartMode)
#if defined CAF_ENABLED || defined MPI_ENABLED
        elseif(self%Image%isFirst) then
            self%Stats%avgCommTimePerFunCall =  self%Stats%avgCommTimePerFunCall / self%Stats%NumFunCall%acceptedRejectedDelayed
            self%Stats%avgTimePerFunCalInSec = (self%Stats%avgTimePerFunCalInSec / (self%Stats%NumFunCall%acceptedRejectedDelayedUnused-acceptedRejectedDelayedUnusedRestartMode)) * self%Image%count
            return
        end if
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine writeOutput()

            implicit none
            integer(IK) :: j

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output write %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ! if new point has been sampled, write the previous sampled point to output file

            blockOutputWrite: if (self%Stats%NumFunCall%accepted>0_IK) then
                if (self%SpecBase%ChainFileFormat%isCompact) then
                    write(self%ChainFile%unit,self%ChainFile%format ) self%Chain%ProcessID(self%Stats%NumFunCall%accepted)      &
                                                                    , self%Chain%DelRejStage(self%Stats%NumFunCall%accepted)    &
                                                                    , self%Chain%MeanAccRate(self%Stats%NumFunCall%accepted)    &
                                                                    , self%Chain%Adaptation(self%Stats%NumFunCall%accepted)     &
                                                                    , self%Chain%BurninLoc(self%Stats%NumFunCall%accepted)      &
                                                                    , self%Chain%Weight(self%Stats%NumFunCall%accepted)         &
                                                                    , self%Chain%LogFunc(self%Stats%NumFunCall%accepted)        &
                                                                    , self%Chain%State(1:nd,self%Stats%NumFunCall%accepted)
                elseif (self%SpecBase%ChainFileFormat%isBinary) then
                    write(self%ChainFile%unit                       ) self%Chain%ProcessID(self%Stats%NumFunCall%accepted)      &
                                                                    , self%Chain%DelRejStage(self%Stats%NumFunCall%accepted)    &
                                                                    , self%Chain%MeanAccRate(self%Stats%NumFunCall%accepted)    &
                                                                    , self%Chain%Adaptation(self%Stats%NumFunCall%accepted)     &
                                                                    , self%Chain%BurninLoc(self%Stats%NumFunCall%accepted)      &
                                                                    , self%Chain%Weight(self%Stats%NumFunCall%accepted)         &
                                                                    , self%Chain%LogFunc(self%Stats%NumFunCall%accepted)        &
                                                                    , self%Chain%State(1:nd,self%Stats%NumFunCall%accepted)
                elseif (self%SpecBase%ChainFileFormat%isVerbose) then
                    do j = 1, self%Chain%Weight(self%Stats%NumFunCall%accepted)
                    write(self%ChainFile%unit,self%ChainFile%format ) self%Chain%ProcessID(self%Stats%NumFunCall%accepted)      &
                                                                    , self%Chain%DelRejStage(self%Stats%NumFunCall%accepted)    &
                                                                    , self%Chain%MeanAccRate(self%Stats%NumFunCall%accepted)    &
                                                                    , AdaptationMeasure(j)                                      &
                                                                   !, self%Chain%Adaptation(self%Stats%NumFunCall%accepted)     &
                                                                    , self%Chain%BurninLoc(self%Stats%NumFunCall%accepted)      &
                                                                    , 1_IK                                                      &
                                                                    , self%Chain%LogFunc(self%Stats%NumFunCall%accepted)        &
                                                                    , self%Chain%State(1:nd,self%Stats%NumFunCall%accepted)
                    end do
                end if
            end if blockOutputWrite

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output write %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end subroutine writeOutput

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Objects that change state in this subroutine are: Timer, timeElapsedUntilLastReportInSeconds, sumAccRateLastReport
        subroutine reportProgress()

            use Constants_mod, only: CARRIAGE_RETURN
            implicit none

            real(RK)                            :: meanAccRateSinceStart
            real(RK)                            :: meanAccRateSinceLastReport
            real(RK)                            :: estimatedTimeToFinishInSeconds
            real(RK)                            :: timeElapsedSinceLastReportInSeconds
            integer(IK)                         :: numFunCallAccepted_dummy ! used in restart mode


            if (self%isFreshRun) then

                call self%Timer%toc()
                timeElapsedSinceLastReportInSeconds = self%Timer%Time%total - timeElapsedUntilLastReportInSeconds
                timeElapsedUntilLastReportInSeconds = self%Timer%Time%total
                meanAccRateSinceStart = SumAccRateSinceStart%acceptedRejected / real(self%Stats%NumFunCall%acceptedRejected,kind=RK)
                meanAccRateSinceLastReport = (SumAccRateSinceStart%acceptedRejected-sumAccRateLastReport) * inverseProgressReportPeriod
                estimatedTimeToFinishInSeconds = getRemainingSimulationFraction() * self%Timer%Time%total

                write( self%TimeFile%unit, self%TimeFile%format ) self%Stats%NumFunCall%acceptedRejected    &
                                                                , self%Stats%NumFunCall%accepted            &
                                                                , meanAccRateSinceStart                     &
                                                                , meanAccRateSinceLastReport                &
                                                                , timeElapsedSinceLastReportInSeconds       &
                                                                , self%Timer%Time%total                     & ! timeElapsedUntilLastReportInSeconds
                                                                , estimatedTimeToFinishInSeconds
                flush(self%TimeFile%unit)

            else

                block
                    use String_mod, only: String_type
                    type(String_type) :: Record
                    allocate( character(600) :: Record%value )
                    read(self%TimeFile%unit, "(A)" ) Record%value
                    Record%Parts = Record%SplitStr(trim(adjustl(Record%value)),self%SpecBase%OutputDelimiter%val,Record%nPart)
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
            if (self%Image%isFirst) then
                write( &
#if defined MEXPRINT_ENABLED
                txt, &
#else
                output_unit, advance = "no", &
#endif
                fmt = "(A,25X,*(A25,3X))" ) &
                CARRIAGE_RETURN, &
                num2str(self%Stats%NumFunCall%accepted)//" / "//num2str(self%Stats%NumFunCall%acceptedRejected,formatIn="(1I10)"), &
                num2str(meanAccRateSinceLastReport,"(1F11.4)")//" / "//num2str(SumAccRateSinceStart%acceptedRejected / real(self%Stats%NumFunCall%acceptedRejected,kind=RK),"(1F11.4)"), &
                num2str(timeElapsedUntilLastReportInSeconds,"(1F11.4)")//" / "//num2str(estimatedTimeToFinishInSeconds,"(1F11.4)")
#if defined MEXPRINT_ENABLED
                call write(string=txt,advance=.false.)
#else
                !call execute_command_line(" ")
                flush(output_unit)
#endif
            end if

            numFunCallAcceptedRejectedLastReport = self%Stats%NumFunCall%acceptedRejected
            sumAccRateLastReport = SumAccRateSinceStart%acceptedRejected

        end subroutine reportProgress

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pure function getRemainingSimulationFraction() result(remainingSimulationFraction)
            use Constants_mod, only: IK, RK
            implicit none
            real(RK) :: remainingSimulationFraction
            remainingSimulationFraction = real(self%SpecMCMC%ChainSize%val-self%Stats%NumFunCall%accepted,kind=RK) / self%Stats%NumFunCall%accepted
        end function getRemainingSimulationFraction

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end subroutine runKernel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SAMPLER_TYPE
#undef SAMPLER_PROPOSAL_ABSTRACT_MOD
