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

!> \brief
!> This file implements the body of the `Kernel` routine of the ParaNest sampler.
!>
!> \author Amir Shahmoradi

submodule (ParaNest_mod) Kernel_smod

    use, intrinsic :: iso_fortran_env, only: output_unit
    !use Constants_mod, only: IK, RK ! gfortran 9.3 compile crashes with this line

#if defined MPI_ENABLED
    use mpi
#endif

    implicit none

    character(*), parameter :: SUBMODULE_NAME = MODULE_NAME // "@Kernel_smod"

    type                    :: SumAccRateSinceStart_type
        real(RK)            :: acceptedRejected
        real(RK)            :: acceptedRejectedDelayed
    end type SumAccRateSinceStart_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a method of [ParaDRAM_type](@ref paradram_type) and [ParaDISE_type](@ref paradise_type) classes.
    !> Run the sampler and return the sampled states.
    !>
    !> @param[inout]    self        :   An object of class [ParaDRAM_type](@ref paradram_type) or [ParaDISE_type](@ref paradise_type).
    !> @param[in]       getLogFunc  :   The target objective function that is to be sampled.
    !>
    !> \remark
    !> This procedure requires preprocessing.
    module subroutine runKernel ( self          &
                                , getLogFunc    &
                                )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: runKernel
#endif
        use Constants_mod, only: RK, IK, NEGINF_RK, NLC, LOGTINY_RK, NEGLOGINF_RK
        use ParaMonteLogFunc_mod, only: getLogFunc_proc
        use Math_mod, only: getLogSubExp
        use Decoration_mod, only: write
        use String_mod, only: num2str

        implicit none

        character(*), parameter             :: PROCEDURE_NAME = SUBMODULE_NAME//"@runKernel()"
        integer(IK) , parameter             :: CHAIN_RESTART_OFFSET = 2_IK

        class(ParaXXXX_type), intent(inout) :: self
        procedure(getLogFunc_proc)          :: getLogFunc

#if defined CAF_ENABLED
        real(RK)    , allocatable           :: co_AccRate(:)[:]
        real(RK)    , allocatable           :: co_LogFuncState(:,:)[:]                              ! (0:nd,-1:delayedRejectionCount), -1 corresponds to the current accepted state
        integer(IK) , save                  :: co_proposalFound_samplerUpdateOccurred(2)[*]         ! merging these scalars would reduce the MPI communication overhead cost: co_proposalFound, co_samplerUpdateOccurred, co_counterDRS, 0 means false, 1 means true
#else
        real(RK)    , allocatable           :: co_LogFuncState(:,:)
        real(RK)    , allocatable           :: co_AccRate(:)
        integer(IK)                         :: co_proposalFound_samplerUpdateOccurred(2)            ! merging these scalars would reduce the MPI communication overhead cost: co_proposalFound, co_samplerUpdateOccurred, co_counterDRS, 0 means false, 1 means true
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
        real(RK)                            :: logFuncDiff                                          ! The difference between the log of the old and the new states. Used to avoid underflow.
       !real(RK)                            :: adaptationMeasureDummy
        real(RK)                            :: maxLogFuncRejectedProposal
        logical                             :: samplerUpdateIsGreedy
        logical                             :: samplerUpdateSucceeded
        logical                             :: delayedRejectionRequested
        logical                             :: noDelayedRejectionRequested
        real(RK)    , allocatable           :: AdaptationMeasure(:), adaptationMeasurePlaceHolder(:)
        integer(IK)                         :: adaptationMeasureLen
        integer(IK)                         :: imageID, dumint
        integer(IK)                         :: nd
        integer(IK) , parameter             :: STDOUT_SEGLEN = 25
        character(4*STDOUT_SEGLEN+2*3+1)    :: txt

        integer(IK)                         :: i,minLogFuncIndex
        integer(IK)                         :: LogFuncIndex(np)                               ! index array that gives the ascending order of the elements of LogFunc active sample.
        real(RK)                            :: Point(nd,np)                                   ! Point is the nd-by-np array of the stochastic variable for the active sample.
        real(RK)                            :: LogFunc(np)                                    ! The vector containing LogFunc values for the active sample
        real(RK)                            :: domainSizeNormed                               ! domainSizeNormed is the prior mass
        real(RK)                            :: LogFuncSwap(0:1)
        real(RK)                            :: domainDiffNormed
        real(RK)                            :: logEvidence
        real(RK)                            :: maxLogFunc
        real(RK)                            :: logMaxError,accr,accrSum,accrMean
        real(RK)                            :: domainShrinkageExponent

#if defined CAF_ENABLED || defined MPI_ENABLED
        integer(IK)                         :: imageStartID, imageEndID
        integer(IK)                         :: proposalFoundSinglChainMode                          ! used in singlChain delayed rejection. zero if the proposal is not accepted. 1 if the proposal is accepted.
#if defined MPI_ENABLED
        integer(IK)                         :: proposalFoundSinglChainModeReduced                   ! the reduced value by summing proposalFoundSinglChainModeReduced over all images.
        real(RK)    , allocatable           :: AccRateMatrix(:,:)                                   ! matrix of size (-1:self%SpecDRAM%DelayedRejectionCount%val,1:self%Image%count)
        integer(IK)                         :: ndPlusOne
        integer(IK)                         :: ierrMPI
        integer(IK)                         :: delayedRejectionCountPlusTwo
#endif
        self%Stats%avgCommTimePerFunCall = 0._RK                                                    ! Until the reporting time, this is in reality, sumCommTimePerFunCall. This is meaningful only in singlChain parallelism
        self%Stats%NumFunCall%acceptedRejectedDelayedUnused = self%Image%count                      ! used only in singlChain parallelism, and relevant only on the first image
#endif
        self%Stats%avgTimePerFunCalInSec = 0._RK
        numFunCallAcceptedRejectedLastReport = 0_IK
        timeElapsedUntilLastReportInSeconds = 0._RK
        inverseProgressReportPeriod = 1._RK/real(self%SpecBase%ProgressReportPeriod%val,kind=RK)    ! this remains a constant except for the last the last report of the simulation
        sumAccRateLastReport = 0._RK
        nd = self%nd%val

        if (allocated(co_AccRate)) deallocate(co_AccRate)
        if (allocated(co_LogFuncState)) deallocate(co_LogFuncState)
#if defined CAF_ENABLED
        allocate(co_LogFuncState(0:nd,self%SpecNest%LiveSampleSize%val)[*])
#else
        allocate(co_LogFuncState(0:nd,self%SpecNest%LiveSampleSize%val))
        !allocate(co_AccRate(-1:self%SpecDRAM%DelayedRejectionCount%val))                            ! the negative element will contain counterDRS
#endif
        !co_AccRate(-1)  = 0._RK                                                                     ! the real-value counterDRS, indicating the initial delayed rejection stage at which the first point is sampled
        !co_AccRate(0)   = 1._RK                                                                     ! initial acceptance rate for the first zeroth DR stage.
        !co_AccRate(1:self%SpecDRAM%DelayedRejectionCount%val) = 0._RK                               ! indicates the very first proposal acceptance on image 1

#if defined MPI_ENABLED
        !if (allocated(AccRateMatrix)) deallocate(AccRateMatrix)
        !allocate(AccRateMatrix(-1:self%SpecDRAM%DelayedRejectionCount%val,1:self%Image%count))      ! the negative element will contain counterDRS
        !AccRateMatrix = 0._RK                                                                       ! -huge(1._RK)  ! debug
        !AccRateMatrix(0,1:self%Image%count) = 1._RK                                                 ! initial acceptance rate for the first zeroth DR stage.
        ndPlusOne = nd + 1_IK
#endif

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

        logMaxError     = 0._RK
        accrSum         = 0._RK
        accr            = 1._RK
        accrMean        = 1._RK
        LogFuncSwap(0)  = -huge(LogFuncSwap(0))    ! Required by incrementEvidence()
        minLogFuncIndex = minloc(LogFunc,1)

        domainShrinkageExponent = 1._RK / real(np,kind=RK)

        call self%Timer%tic()

        blockDryFreshRunSetup: if (self%isFreshRun) then

            do concurrent(ip = 1:self%SpecNset%LiveSampleSize%val)
                co_LogFuncState(1:nd,ip) = self%Proposal%initializeActiveSample(nd,np)
            end do

            allocate(self%Chain%ProcessID   (   self%SpecMCMC%ChainSize%val))
            allocate(self%Chain%MeanAccRate (   self%SpecMCMC%ChainSize%val))
            allocate(self%Chain%Weight      (   self%SpecMCMC%ChainSize%val))
            allocate(self%Chain%State       (nd,self%SpecMCMC%ChainSize%val))
            allocate(self%Chain%LogFunc     (   self%SpecMCMC%ChainSize%val))

        else blockDryFreshRunSetup

!            ! load the existing Chain file into self%Chain components
!
!            call self%Chain%get ( chainFilePath = self%ChainFile%Path%original      & ! LCOV_EXCL_LINE
!                                , chainFileForm = self%SpecBase%ChainFileFormat%val & ! LCOV_EXCL_LINE
!                                , Err = self%Err                                    & ! LCOV_EXCL_LINE
!                                , targetChainSize = self%SpecMCMC%ChainSize%val     & ! LCOV_EXCL_LINE
!                                , lenHeader = self%Chain%lenHeader                  & ! LCOV_EXCL_LINE
!                                , ndim = nd                                         & ! LCOV_EXCL_LINE
!                                , delimiter = self%SpecBase%OutputDelimiter%val     & ! LCOV_EXCL_LINE
!                                )
!            if (self%Err%occurred) then
!                ! LCOV_EXCL_START
!                self%Err%msg = PROCEDURE_NAME//self%Err%msg
!                call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
!                return
!                ! LCOV_EXCL_STOP
!            end if
!
!            if (self%Chain%Count%compact<=CHAIN_RESTART_OFFSET) then
!                self%isFreshRun = .true.
!                self%isDryRun = .not. self%isFreshRun
!            end if
!
!            call self%Image%sync()
!
!            blockLeaderSetup: if (self%Image%isLeader) then
!
!                ! set up the chain file
!
!                block
!
!                    use System_mod, only: RandomFileName_type, copyFile, removeFile
!                    type(RandomFileName_type) :: RFN
!
!                    ! create a copy of the chain file, just for the sake of not losing the simulation results
!
!                    RFN = RandomFileName_type(dir = "", key = self%ChainFile%Path%original//".rst", ext="") ! temporary_restart_copy
!                    call copyFile(pathOld=self%ChainFile%Path%original,pathNew=RFN%path,isUnixShell=self%OS%Shell%isUnix,Err=self%err)
!                    if (self%Err%occurred) then
!                        ! LCOV_EXCL_START
!                        self%Err%msg = PROCEDURE_NAME//self%Err%msg
!                        call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
!                        exit blockLeaderSetup
!                        return
!                        ! LCOV_EXCL_STOP
!                    end if
!
!                    ! reopen the chain file to resume the simulation
!
!                    open( newunit = self%ChainFile%unit             & ! LCOV_EXCL_LINE
!                        , file = self%ChainFile%Path%original       & ! LCOV_EXCL_LINE
!                        , form = self%ChainFile%Form%value          & ! LCOV_EXCL_LINE
!                        , status = self%ChainFile%status            & ! LCOV_EXCL_LINE
!                        , iostat = self%ChainFile%Err%stat          & ! LCOV_EXCL_LINE
!#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
!                        , SHARED                                    & ! LCOV_EXCL_LINE
!#endif
!                        , position = self%ChainFile%Position%value  )
!                    self%Err = self%ChainFile%getOpenErr(self%ChainFile%Err%stat)
!                    if (self%Err%occurred) then
!                        ! LCOV_EXCL_START
!                        self%Err%msg = PROCEDURE_NAME//": Error occurred while opening "//self%name//" "//self%ChainFile%suffix//" file='"//self%ChainFile%Path%original//"'."
!                        call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
!                        exit blockLeaderSetup
!                        return
!                        ! LCOV_EXCL_STOP
!                    end if
!
!                    ! rewrite the chain file
!
!                    call self%Chain%writeChainFile  ( ndim = nd & ! LCOV_EXCL_LINE
!                                                    , compactStartIndex = 1_IK & ! LCOV_EXCL_LINE
!                                                    , compactEndIndex = self%Chain%Count%compact-CHAIN_RESTART_OFFSET & ! LCOV_EXCL_LINE
!                                                    , chainFileUnit = self%ChainFile%unit & ! LCOV_EXCL_LINE
!                                                    , chainFileForm = self%SpecBase%ChainFileFormat%val & ! LCOV_EXCL_LINE
!                                                    , chainFileFormat = self%ChainFile%format & ! LCOV_EXCL_LINE
!                                                    , adaptiveUpdatePeriod = self%SpecDRAM%AdaptiveUpdatePeriod%val & ! LCOV_EXCL_LINE
!                                                    )
!
!                    ! remove the temporary copy of the chain file
!
!                    call removeFile(RFN%path, self%Err) ! Passing the optional Err argument, handles exceptions should any occur.
!
!                end block
!
            end if blockLeaderSetup

#if (defined MPI_ENABLED || defined CAF_ENABLED) && (defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED)
            block; use Err_mod, only: bcastErr; call bcastErr(self%Err); end block
#endif
            if (self%Err%occurred) return

        end if blockDryFreshRunSetup

        if (self%isFreshRun) then ! this must be done separately from the above blockDryFreshRunSetup

            co_LogFuncState(1:nd) = self%SpecMCMC%StartPointVec%Val   ! proposal state
            call self%Timer%toc()

            call self%Timer%toc(); self%Stats%avgTimePerFunCalInSec = self%Stats%avgTimePerFunCalInSec + self%Timer%Time%delta

        else

            co_LogFuncState(1:nd,0)     = self%Chain%State(1:nd,1)      ! proposal state
            co_LogFuncState(0,0)        = self%Chain%LogFunc(1)         ! proposal logFunc

        end if

        if (self%Image%isFirst) then
            ! LCOV_EXCL_START
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
            ! LCOV_EXCL_STOP
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

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start of loopNestSampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Start sampling and logEvidence calculation
        loopNestedSampling: do

            if (mod(nCall2LogFuncAccepted,iverbose)==0) write(logFileUnit,*) nCall2LogFuncAccepted,' points sampled. Total function calls: ', nCall2LogFuncTotal+np-1
            if (mod(nCall2LogFuncAccepted,ireport)==0) then
                call writeRestartFile(nd,np,Point,LogFunc)    ! write restart file
            end if


            LogFuncSwap(1)   = LogFunc(minLogFuncIndex)     ! LogFuncSwap(0) is the minimum LogFunc in the old active set.

            call getDomainSize(domainShrinkageExponent,domainSizeNormed,domainDiffNormed)      ! Calculates domainDiffNormed. Swaps domainSizeNormed.
            call incrementEvidence(domainDiffNormed,LogFuncSwap,logEvidence)

            LogFuncSwap(0) = LogFuncSwap(1)

            ! Now write out the sample
            nCall2LogFuncAccepted = nCall2LogFuncAccepted + 1
            accrSum  = accrSum + accr
            accrMean = accrSum / real(nCall2LogFuncAccepted,kind=RK)
            accr = real(nCall2LogFuncAccepted,kind=RK) / real(nCall2LogFuncTotal,kind=RK)
            call passResult( nd, nCall2LogFuncAccepted, nCall2LogFuncTotal, accr, accrMean, domainSizeNormed, domainDiffNormed &
                            , LogFuncSwap(1), logDomainSize+logEvidence, logMaxError, Point(1:nd,minLogFuncIndex) )

            ! check for stopping criterion. Exit loop if integration accuracy achieved.
            maxLogFunc = maxval(LogFunc)
            logMaxError = maxLogFunc + log(domainSizeNormed) - logEvidence
            if (logMaxError < logTolerance) then  ! Write out the remaining part of the logEvidence and exit sampling
                write(*,*) 'The desired accuracy for logEvidence achieved'
                call indexArrayReal(np,LogFunc,LogFuncIndex)   ! get the ascending order of the elements of LogFunc
                domainDiffNormed = 0.5_RK*domainSizeNormed/real(np)
                do i = 2,np    ! Add the final contribution to logEvidence
                    LogFuncSwap(1) = LogFunc(LogFuncIndex(i))
                    call incrementEvidence(domainDiffNormed,LogFuncSwap,logEvidence)
                    LogFuncSwap(0)  = LogFuncSwap(1)
                    nCall2LogFuncAccepted = nCall2LogFuncAccepted + 1
                    nCall2LogFuncTotal = nCall2LogFuncTotal + 1
                    accrSum  = accrSum + accr
                    accrMean = accrSum / real(nCall2LogFuncAccepted,kind=RK)
                    accr = real(nCall2LogFuncAccepted,kind=RK) / real(nCall2LogFuncTotal,kind=RK)
                    call passResult( nd, nCall2LogFuncAccepted, nCall2LogFuncTotal, accr, accrMean, domainSizeNormed, domainDiffNormed, LogFuncSwap(1) &
                                    , logDomainSize+logEvidence, logMaxError, Point(1:nd,LogFuncIndex(i)) )
                end do
                exit loopNestedSampling
            end if

            ! Now draw a new point subject to getLogFunc > LogFunc(minLogFuncIndex)
            ! Updates LogFunc and Point arrays with new point replacing the lowest-likelihood point.
            ! Replaces the lowest likelihood and the associated Parameter set with the new proposal values.
            call getNewPointCARS(nd,np,iverbose,nCall2LogFuncTotal,minLogFuncIndex,Point,LogFunc,getLogFunc)

            cycle loopNestedSampling

        end do loopNestedSampling

    end subroutine runLMC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine getDomainSize(domainShrinkageExponent,domainSizeNormed,domainDiffNormed)
        implicit none
        real(RK), intent(in)    :: domainShrinkageExponent
        real(RK), intent(out)   :: domainSizeNormed
        real(RK), intent(out)   :: domainDiffNormed
        real(RK), save          :: domainSizeNormedCurrent = 1._RK
        real(RK)                :: rnd
        call random_number(rnd)
        domainSizeNormed = domainSizeNormedCurrent * rnd**domainShrinkageExponent
        domainDiffNormed = 0.5_RK * ( domainSizeNormedCurrent - domainSizeNormed )
        domainSizeNormedCurrent = domainSizeNormed
    end subroutine getDomainSize

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! Increments logEvidence using the minimum logFunc and the corresponding stochastic prior volume associated with it.
    subroutine incrementEvidence(domainDiffNormed,LogFuncSwap,logEvidence)
        use Constants_mod only: RK
        implicit none
        real(RK), intent(in)  :: domainDiffNormed,LogFuncSwap(0:1)
        real(RK), intent(out) :: logEvidence
        real(RK), save        :: evidenceNormed = 0._RK
        real(RK)              :: logFuncRatio
        logFuncRatio = exp( LogFuncSwap(0) - LogFuncSwap(1) )        ! This ratio is by definition always < 1.
        evidenceNormed  = ( evidenceNormed + domainDiffNormed ) * logFuncRatio + domainDiffNormed  ! evidence = evidenceNormed * exp(LogFuncSwap(1))
        logEvidence     = log( evidenceNormed ) + LogFuncSwap(1)
        !!!XX write(*,*) 'LogFuncSwap: ', LogFuncSwap; read(*,*)
        !!!XX write(*,*) 'logFuncRatio: ', logFuncRatio; read(*,*)
        !!!XX write(*,*) 'evidenceNormed, logEvidence: ', evidenceNormed, logEvidence; read(*,*)
        !!!XX write(*,*) 'log[ domainDiffNormed*exp(LogFuncSwap(1)) ] : ', log(domainDiffNormed*exp(LogFuncSwap(1))); read(*,*)
    end subroutine incrementEvidence

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule Kernel_smod