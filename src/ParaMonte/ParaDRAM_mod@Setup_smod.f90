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

submodule (ParaDRAM_mod) Setup_smod

    !use Constants_mod, only: IK, RK ! gfortran 9.3 compile crashes with this line
    implicit none

    character(*), parameter :: SUBMODULE_NAME = MODULE_NAME // "@Setup_smod"

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    module subroutine runSampler( PD                                    &
                                , ndim                                  &
                                , getLogFunc                            &
                                , inputFile                             &
                                ! ParaMonte variables
                                , sampleSize                            &
                                , randomSeed                            &
                                , description                           &
                                , outputFileName                        &
                                , outputDelimiter                       &
                                , chainFileFormat                       &
                                , variableNameList                      &
                                , restartFileFormat                     &
                                , outputColumnWidth                     &
                                , outputRealPrecision                   &
                                , silentModeRequested                   &
                                , domainLowerLimitVec                   &
                                , domainUpperLimitVec                   &
                                , parallelizationModel                  &
                                , progressReportPeriod                  &
                                , TargetAcceptanceRate                  &
                                , mpiFinalizeRequested                  &
                                , maxNumDomainCheckToWarn               &
                                , maxNumDomainCheckToStop               &
                                ! ParaMCMC variables
                                , chainSize                             &
                                , startPointVec                         &
                                , sampleRefinementCount                 &
                                , sampleRefinementMethod                &
                                , randomStartPointRequested             &
                                , randomStartPointDomainLowerLimitVec   &
                                , randomStartPointDomainUpperLimitVec   &
                                ! ParaDRAM variables
                                , scaleFactor                           &
                                , proposalModel                         &
                                , ProposalStartCovMat                   &
                                , ProposalStartCorMat                   &
                                , ProposalStartStdVec                   &
                                , adaptiveUpdateCount                   &
                                , adaptiveUpdatePeriod                  &
                                , greedyAdaptationCount                 &
                                , delayedRejectionCount                 &
                                , burninAdaptationMeasure               &
                                , delayedRejectionScaleFactorVec        &
                                ) !  result(PD)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !!DEC$ ATTRIBUTES DLLEXPORT :: ParaDRAM_type
        !DEC$ ATTRIBUTES DLLEXPORT :: runSampler
#endif
        use, intrinsic :: iso_fortran_env, only: output_unit, stat_stopped_image
        use ParaDRAMProposalSymmetric_mod, only: ProposalSymmetric_type
        use ParaMonteLogFunc_mod, only: getLogFunc_proc
        use Statistics_mod, only: getUpperCorMatFromUpperCovMat
        use Statistics_mod, only: getWeiSamCovUppMeanTrans
        use Statistics_mod, only: getQuantile
        use ParaMonte_mod, only: QPROB
        use Constants_mod, only: IK, RK, NLC, PMSM
        use DateTime_mod, only: getNiceDateTime
        use String_mod, only: num2str

        implicit none

        ! self
        class(ParaDRAM_type), intent(inout) :: PD

        ! mandatory variables
        integer(IK) , intent(in)            :: ndim
        procedure(getLogFunc_proc)          :: getLogFunc

        character(*), intent(in), optional  :: inputFile
        
        ! ParaMonte variables
        integer(IK) , intent(in), optional  :: sampleSize
        integer(IK) , intent(in), optional  :: randomSeed
        character(*), intent(in), optional  :: description
        character(*), intent(in), optional  :: outputFileName
        character(*), intent(in), optional  :: outputDelimiter
        character(*), intent(in), optional  :: chainFileFormat
        character(*), intent(in), optional  :: variableNameList(ndim)
        character(*), intent(in), optional  :: restartFileFormat
        integer(IK) , intent(in), optional  :: outputColumnWidth
        integer(IK) , intent(in), optional  :: outputRealPrecision
        real(RK)    , intent(in), optional  :: domainLowerLimitVec(ndim)
        real(RK)    , intent(in), optional  :: domainUpperLimitVec(ndim)
        logical     , intent(in), optional  :: silentModeRequested
        character(*), intent(in), optional  :: parallelizationModel
        integer(IK) , intent(in), optional  :: progressReportPeriod
        real(RK)    , intent(in), optional  :: TargetAcceptanceRate(2)
        logical     , intent(in), optional  :: mpiFinalizeRequested
        integer(IK) , intent(in), optional  :: maxNumDomainCheckToWarn
        integer(IK) , intent(in), optional  :: maxNumDomainCheckToStop

        ! ParaMCMC variables
        integer(IK) , intent(in), optional  :: chainSize
        real(RK)    , intent(in), optional  :: startPointVec(ndim)
        integer(IK) , intent(in), optional  :: sampleRefinementCount
        character(*), intent(in), optional  :: sampleRefinementMethod
        logical     , intent(in), optional  :: randomStartPointRequested
        real(RK)    , intent(in), optional  :: randomStartPointDomainLowerLimitVec(ndim)
        real(RK)    , intent(in), optional  :: randomStartPointDomainUpperLimitVec(ndim)

        ! ParaDRAM variables
        character(*), intent(in), optional  :: scaleFactor
        character(*), intent(in), optional  :: proposalModel
        real(RK)    , intent(in), optional  :: proposalStartCovMat(ndim,ndim)
        real(RK)    , intent(in), optional  :: proposalStartCorMat(ndim,ndim)
        real(RK)    , intent(in), optional  :: proposalStartStdVec(ndim)
        integer(IK) , intent(in), optional  :: adaptiveUpdateCount
        integer(IK) , intent(in), optional  :: adaptiveUpdatePeriod
        integer(IK) , intent(in), optional  :: greedyAdaptationCount
        integer(IK) , intent(in), optional  :: delayedRejectionCount
        real(RK)    , intent(in), optional  :: burninAdaptationMeasure
        real(RK)    , intent(in), optional  :: delayedRejectionScaleFactorVec(:)

        character(*), parameter             :: PROCEDURE_NAME = SUBMODULE_NAME // "@runSampler()"

       !type(ProposalSymmetric_type), target   :: ProposalSymmetric
        integer(IK)                         :: i, iq
        character(:), allocatable           :: formatStr, logFileColWidthStr, msg
        real(RK)                            :: mcmcSamplingEfficiency

        !***************************************************************************************************************************
        ! Initialize SpecBase variables, then check existence of inputFile and open it and return the unit file, if it exists
        !***************************************************************************************************************************

        call PD%setupParaMonte  ( nd = ndim               &
                                , name = PMSM%ParaDRAM    &
                                , date = "May 23 2018"    &
                                , version = "1.0.0"       &
                                , inputFile = inputFile   &
                                )
        if (PD%Err%occurred) return

        !***************************************************************************************************************************
        ! Initialize ParaMCMC variables
        !***************************************************************************************************************************

        call PD%setupParaMCMC()

        !***************************************************************************************************************************
        ! Initialize ParaDRAM variables
        !***************************************************************************************************************************

        PD%SpecDRAM = SpecDRAM_type( PD%nd%val, PD%name ) ! , PD%SpecMCMC%ChainSize%def )

        !***************************************************************************************************************************
        ! read variables from input file if it exists
        !***************************************************************************************************************************

        call PD%getSpecFromInputFile(ndim)
        if (PD%Err%occurred) then
            PD%Err%msg = PROCEDURE_NAME // PD%Err%msg
            call PD%abort( Err = PD%Err, prefix = PD%brand, newline = "\n", outputUnit = PD%LogFile%unit )
            return
        end if

        !***************************************************************************************************************************
        ! read variables from argument list if needed
        !***************************************************************************************************************************

        if (PD%Image%isFirst) call PD%setWarnAboutProcArgHasPriority()
        if (PD%procArgNeeded) then
            call PD%SpecBase%setFromInputArgs   ( Err                                   = PD%Err                                &
                                                , sampleSize                            = sampleSize                            &
                                                , randomSeed                            = randomSeed                            &
                                                , description                           = description                           &
                                                , outputFileName                        = outputFileName                        &
                                                , outputDelimiter                       = outputDelimiter                       &
                                                , chainFileFormat                       = chainFileFormat                       &
                                                , variableNameList                      = variableNameList                      &
                                                , restartFileFormat                     = restartFileFormat                     &
                                                , outputColumnWidth                     = outputColumnWidth                     &
                                                , outputRealPrecision                   = outputRealPrecision                   &
                                                , silentModeRequested                   = silentModeRequested                   &
                                                , domainLowerLimitVec                   = domainLowerLimitVec                   &
                                                , domainUpperLimitVec                   = domainUpperLimitVec                   &
                                                , parallelizationModel                  = parallelizationModel                  &
                                                , progressReportPeriod                  = progressReportPeriod                  &
                                                , TargetAcceptanceRate                  = TargetAcceptanceRate                  &
                                                , mpiFinalizeRequested                  = mpiFinalizeRequested                  &
                                                , maxNumDomainCheckToWarn               = maxNumDomainCheckToWarn               &
                                                , maxNumDomainCheckToStop               = maxNumDomainCheckToStop               &
                                                )
            call PD%SpecMCMC%setFromInputArgs   ( domainLowerLimitVec                   = PD%SpecBase%DomainLowerLimitVec%Val   &
                                                , domainUpperLimitVec                   = PD%SpecBase%DomainUpperLimitVec%Val   &
                                                , chainSize                             = chainSize                             &
                                                , scaleFactor                           = scaleFactor                           &
                                                , startPointVec                         = startPointVec                         &
                                                , proposalModel                         = proposalModel                         &
                                                , proposalStartCovMat                   = proposalStartCovMat                   &
                                                , proposalStartCorMat                   = proposalStartCorMat                   &
                                                , proposalStartStdVec                   = proposalStartStdVec                   &
                                                , sampleRefinementCount                 = sampleRefinementCount                 &
                                                , sampleRefinementMethod                = sampleRefinementMethod                &
                                                , randomStartPointRequested             = randomStartPointRequested             &
                                                , randomStartPointDomainLowerLimitVec   = randomStartPointDomainLowerLimitVec   &
                                                , randomStartPointDomainUpperLimitVec   = randomStartPointDomainUpperLimitVec   &
                                                )
            call PD%SpecDRAM%setFromInputArgs   ( adaptiveUpdateCount                   = adaptiveUpdateCount                   &
                                                , adaptiveUpdatePeriod                  = adaptiveUpdatePeriod                  &
                                                , greedyAdaptationCount                 = greedyAdaptationCount                 &
                                                , delayedRejectionCount                 = delayedRejectionCount                 &
                                                , burninAdaptationMeasure               = burninAdaptationMeasure               &
                                                , delayedRejectionScaleFactorVec        = delayedRejectionScaleFactorVec        &
                                                )
        end if
        if (PD%Err%occurred) then
            PD%Err%msg = PROCEDURE_NAME // PD%Err%msg
            call PD%abort( Err = PD%Err, prefix = PD%brand, newline = "\n", outputUnit = PD%LogFile%unit )
            return
        end if

        !***************************************************************************************************************************
        ! Now depending on the requested parallelism type, determine the process master/slave type and open output files
        !***************************************************************************************************************************

        PD%Image%isMaster = .false.
        if (PD%SpecBase%ParallelizationModel%isSinglChain) then
            if (PD%Image%isFirst) PD%Image%isMaster = .true.
        elseif (PD%SpecBase%ParallelizationModel%isMultiChain) then
            PD%Image%isMaster = .true.
        else
            PD%Err%msg = PROCEDURE_NAME//": Error occurred. Unknown parallelism requested via the input variable parallelizationModel='"//PD%SpecBase%ParallelizationModel%val//"'."
            call PD%abort( Err = PD%Err, prefix = PD%brand, newline = "\n", outputUnit = PD%LogFile%unit )
            return
        end if
        PD%Image%isNotMaster = .not. PD%Image%isMaster

        !***************************************************************************************************************************
        ! setup log file and open output files
        !***************************************************************************************************************************

        call PD%setupOutputFiles()
        if (PD%Err%occurred) then
            PD%Err%msg = PROCEDURE_NAME // PD%Err%msg
            call PD%abort( Err = PD%Err, prefix = PD%brand, newline = NLC, outputUnit = PD%LogFile%unit )
            return
        end if

        !***************************************************************************************************************************
        ! report variable values to the log file
        !***************************************************************************************************************************

        if (PD%isFreshRun) then
            call PD%SpecBase%reportValues(prefix=PD%brand,outputUnit=PD%LogFile%unit,isMasterImage=PD%Image%isMaster)
            call PD%SpecMCMC%reportValues(prefix=PD%brand,outputUnit=PD%LogFile%unit,isMasterImage=PD%Image%isMaster,methodName=PD%name,splashModeRequested=PD%SpecBase%SilentModeRequested%isFalse)
            call PD%SpecDRAM%reportValues(prefix=PD%brand,outputUnit=PD%LogFile%unit,isMasterImage=PD%Image%isMaster,splashModeRequested=PD%SpecBase%SilentModeRequested%isFalse)
        end if

        !***************************************************************************************************************************
        ! perform variable value sanity checks
        !***************************************************************************************************************************

        PD%Err%msg = ""
        PD%Err%occurred = .false.

        call PD%SpecBase%checkForSanity ( Err = PD%Err          &
                                        , methodName = PD%name  &
                                        )

        call PD%SpecMCMC%checkForSanity ( Err = PD%Err                                              &
                                        , methodName = PD%name                                      &
                                        , nd = ndim                                                 &
                                        , domainLowerLimitVec = PD%SpecBase%DomainLowerLimitVec%Val &
                                        , domainUpperLimitVec = PD%SpecBase%DomainUpperLimitVec%Val &
                                        )

        call PD%SpecDRAM%checkForSanity ( Err = PD%Err          &
                                        , methodName = PD%name  &
                                        , nd = ndim             &
                                        )

        if (PD%Image%isMaster) then

            if (PD%Err%occurred) then
                PD%Err%msg = PROCEDURE_NAME // PD%Err%msg
                call PD%abort( Err = PD%Err, prefix = PD%brand, newline = "\n", outputUnit = PD%LogFile%unit )
                return
            end if

        end if

        !***************************************************************************************************************************
        ! setup output files
        !***************************************************************************************************************************

        block
            use ParaDRAMChainFileContents_mod, only: ChainFileContents_type
            PD%Chain = ChainFileContents_type( ndim = ndim, variableNameList = PD%SpecBase%VariableNameList%Val )
        end block
        if (PD%isFreshRun) then
            if (PD%Image%isMaster) then
                write( PD%TimeFile%unit, PD%TimeFile%format ) "NumFuncCallTotal" &
                                                            , "NumFuncCallAccepted" &
                                                            , "MeanAcceptanceRateSinceStart" &
                                                            , "MeanAcceptanceRateSinceLastReport" &
                                                            , "TimeElapsedSinceLastReportInSeconds" &
                                                            , "TimeElapsedSinceStartInSeconds" &
                                                            , "TimeRemainedToFinishInSeconds"
                call PD%Chain%writeHeader   ( ndim              = ndim &
                                            , chainFileUnit     = PD%ChainFile%unit &
                                            , isBinary          = PD%SpecBase%ChainFileFormat%isBinary &
                                            , chainFileFormat   = PD%ChainFile%format &
                                            )
            end if
        else
            if (PD%Image%isMaster) read(PD%TimeFile%unit,*) ! read the header line of the time file, only by master images
            call PD%Chain%getLenHeader(ndim=ndim,isBinary=PD%SpecBase%ChainFileFormat%isBinary,chainFileFormat=PD%ChainFile%format)
        end if

        ! The format setup of setupOutputFiles() uses the generic g0 edit descriptor. Here the format is revised to be more specific.
        ! g0 edit descriptor's format is slightly more arbitrary and compiler-dependent.

        if (PD%SpecBase%OutputColumnWidth%val>0) then
            PD%TimeFile%format  =   "("     // &
                                    "2I"    // PD%SpecBase%OutputColumnWidth%str // ",'" // PD%SpecBase%OutputDelimiter%val // &
                                    "',5(E" // PD%SpecBase%OutputColumnWidth%str // "." // PD%SpecBase%OutputRealPrecision%str // &
                                    ",:,'"  // PD%SpecBase%OutputDelimiter%val // "'))"
            PD%ChainFile%format =   "("     // &
                                    "2(I"   // PD%SpecBase%OutputColumnWidth%str // &
                                    ",'"    // PD%SpecBase%OutputDelimiter%val // "')," // &
                                    "2(E"   // PD%SpecBase%OutputColumnWidth%str // "." // PD%SpecBase%OutputRealPrecision%str // &
                                    ",'"    // PD%SpecBase%OutputDelimiter%val // "')," // &
                                    "2(I"   // PD%SpecBase%OutputColumnWidth%str // &
                                    ",'"    // PD%SpecBase%OutputDelimiter%val // "')" // &
                                    ","     // num2str(3_IK+PD%nd%val) // &
                                    "(E"    // PD%SpecBase%OutputColumnWidth%str // "." // PD%SpecBase%OutputRealPrecision%str // &
                                    ",:,'"  // PD%SpecBase%OutputDelimiter%val // "')" // &
                                    ")"
        end if

        !***************************************************************************************************************************
        ! setup the proposal distribution
        !***************************************************************************************************************************

        if (PD%SpecMCMC%ProposalModel%isNormal .or. PD%SpecMCMC%ProposalModel%isUniform) then
            !ProposalSymmetric = ProposalSymmetric_type(ndim = ndim, PD = PD)
            !PD%Proposal => ProposalSymmetric
            allocate( PD%Proposal, source = ProposalSymmetric_type  ( ndim          = ndim &
                                                                    , SpecBase      = PD%SpecBase &
                                                                    , SpecMCMC      = PD%SpecMCMC &
                                                                    , SpecDRAM      = PD%SpecDRAM &
                                                                    , Image         = PD%Image &
                                                                    , name          = PD%name &
                                                                    , brand         = PD%brand &
                                                                    , LogFile       = PD%LogFile &
                                                                    , RestartFile   = PD%RestartFile &
                                                                    ) )
#if MATLAB_ENABLED && !defined CAF_ENABLED && !defined MPI_ENABLED
            block
                use ParaDRAMProposal_mod, only: ProposalErr
                if(ProposalErr%occurred) return
            end block
#endif
        !elseif (PD%SpecMCMC%ProposalModel%isUniform) then
        !    allocate( PD%Proposal, source = ProposalUniform_type(ndim = ndim, PD = PD) )
        else
            PD%Err%msg = PROCEDURE_NAME // ": Internal error occurred. Unsupported proposal distribution for " // PD%name // "."
            call PD%abort( Err = PD%Err, prefix = PD%brand, newline = "\n", outputUnit = PD%LogFile%unit )
            return
        end if

        !***************************************************************************************************************************
        ! run ParaDRAM kernel
        !***************************************************************************************************************************

        if (PD%isFreshRun .and. PD%Image%isMaster) then
            call PD%Decor%writeDecoratedText( text = "\nStarting " // PD%name // " sampling - " // getNiceDateTime() // "\n" &
                                            , marginTop = 1_IK  &
                                            , marginBot = 1_IK  &
                                            , newline = "\n"    &
                                            , outputUnit = PD%LogFile%unit )
        end if

        call PD%runKernel( getLogFunc = getLogFunc )
#if MATLAB_ENABLED && !defined CAF_ENABLED && !defined MPI_ENABLED
        if(PD%Err%occurred) return
#endif

        if (PD%isFreshRun .and. PD%Image%isMaster) then
            call PD%Decor%writeDecoratedText( text = "\nExiting " // PD%name // " sampling - " // getNiceDateTime() // "\n" &
                                            , marginTop = 1_IK  &
                                            , marginBot = 1_IK  &
                                            , newline = "\n"    &
                                            , outputUnit = PD%LogFile%unit )
        end if

        !***************************************************************************************************************************
        ! start ParaDRAM post-processing
        !***************************************************************************************************************************

        blockMasterPostProcessing: if (PD%Image%isMaster) then

            logFileColWidthStr = num2str( max( PD%SpecBase%OutputRealPrecision%val, PD%SpecBase%VariableNameList%MaxLen%val ) + 9_IK )

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Total number of accepted function calls (unique samples):" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str(PD%Stats%NumFunCall%accepted) )

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Total number of accepted or rejected function calls:" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str(PD%Stats%NumFunCall%acceptedRejected) )

            !if (PD%SpecDRAM%DelayedRejectionCount%val>0_IK) then
            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Total number of accepted or rejected or delayed-rejection (if any requested) function calls:" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( PD%Stats%NumFunCall%acceptedRejectedDelayed ) )
            !end if

#if defined CAF_ENABLED || defined MPI_ENABLED
            !if (PD%SpecBase%ParallelizationModel%isSinglChain) then
            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Total number of accepted or rejected or unused function calls (by all processes, including delayed rejections, if any requested):" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str(PD%Stats%NumFunCall%acceptedRejectedDelayedUnused) )
            !end if
#endif

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Average MCMC acceptance rate:" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str  ( PD%Chain%MeanAccRate(PD%Stats%NumFunCall%accepted) ) )

            mcmcSamplingEfficiency = real(PD%Stats%NumFunCall%accepted,kind=RK) / real(PD%Stats%NumFunCall%acceptedRejected,kind=RK)
            call PD%Decor%write(PD%LogFile%unit,1,0,1, "MCMC sampling efficiency [ = acceptedFunctionCalls / acceptedPlusRejectedFunctionCalls ]:" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str  ( mcmcSamplingEfficiency ) )

            !if (PD%SpecDRAM%DelayedRejectionCount%val>0_IK) then
            call PD%Decor%write(PD%LogFile%unit,1,0,1, "MCMC sampling efficiency (including delayed rejections, if any requested) [ = acceptedFunctionCalls / acceptedPlusRejectedPlusDelayedRejectionFunctionCalls ]:" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str  ( real(PD%Stats%NumFunCall%accepted,kind=RK) / real(PD%Stats%NumFunCall%acceptedRejectedDelayed,kind=RK) ) )
            !end if

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Total runtime in seconds:" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str(PD%Timer%Time%total) )

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Average effective time cost of each accepted function call, in seconds:" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( PD%Timer%Time%total / real(PD%Stats%NumFunCall%accepted,kind=RK) ) )

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Average effective time cost of each accepted or rejected function call, in seconds:" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( PD%Timer%Time%total / real(PD%Stats%NumFunCall%acceptedRejected,kind=RK) ) )

            if (PD%SpecDRAM%DelayedRejectionCount%val>0_IK) then
            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Average effective time cost of each accepted or rejected or delayed-rejection function call, in seconds:" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( PD%Timer%Time%total / real(PD%Stats%NumFunCall%acceptedRejectedDelayed,kind=RK) ) )
            end if

!#if defined CAF_ENABLED || defined MPI_ENABLED
            !if (PD%SpecBase%ParallelizationModel%isSinglChain) then
            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Average effective time cost of each accepted or rejected or unused function call (including delayed-rejections, if any requested), in seconds:" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( PD%Timer%Time%total / real(PD%Stats%NumFunCall%acceptedRejectedDelayedUnused,kind=RK) ) )

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Average time cost of inter-process communications per used (accepted or rejected or delayed-rejection) function call, in seconds:" )
            if (PD%Image%count==1_IK) then
                call PD%Decor%write(PD%LogFile%unit,0,1,1, "UNDEFINED" )
            else
                call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( PD%Stats%avgCommTimePerFunCall ) )
            end if
            !end if
!#endif

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Average pure time cost of each function call, in seconds:" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( PD%Stats%avgTimePerFunCalInSec ) )

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Number of processes (images):" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( PD%Image%count ) )

            if (PD%Image%count==1_IK .or. PD%SpecBase%ParallelizationModel%isMultiChain) then

                call PD%Decor%write(PD%LogFile%unit,1,0,1, "Estimated maximum speedup gained via "//PD%SpecBase%ParallelizationModel%val//" parallelization model compared to serial mode:" )
                call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( PD%Image%count ) )

                call PD%Decor%write(PD%LogFile%unit,1,0,1, "Predicted optimal number of physical computing processes for "//PD%SpecBase%ParallelizationModel%val//" parallelization model, given the current MCMC sampling efficiency:" )
                if (PD%SpecBase%ParallelizationModel%isMultiChain) then
                    call PD%Decor%write(PD%LogFile%unit,0,1,1, "+INFINITY" )
                else
                    call PD%Decor%write(PD%LogFile%unit,0,1,1, "UNDEFINED" )
                end if

                call PD%Decor%write(PD%LogFile%unit,1,0,1, "Predicted optimal maximum speedup gained via "//PD%SpecBase%ParallelizationModel%val//" parallelization model, given the current MCMC sampling efficiency:" )
                if (PD%SpecBase%ParallelizationModel%isMultiChain) then
                    call PD%Decor%write(PD%LogFile%unit,0,1,1, "+INFINITY" )
                else
                    call PD%Decor%write(PD%LogFile%unit,0,1,1, "UNDEFINED" )
                end if

                call PD%Decor%write(PD%LogFile%unit,1,0,1, "Predicted absolute optimal number of physical computing processes for "//PD%SpecBase%ParallelizationModel%val//" parallelization model, under any MCMC sampling efficiency:" )
                if (PD%SpecBase%ParallelizationModel%isMultiChain) then
                    call PD%Decor%write(PD%LogFile%unit,0,1,1, "+INFINITY" )
                else
                    call PD%Decor%write(PD%LogFile%unit,0,1,1, "UNDEFINED" )
                end if

                call PD%Decor%write(PD%LogFile%unit,1,0,1, "Predicted absolute optimal maximum speedup gained via "//PD%SpecBase%ParallelizationModel%val//" parallelization model, under any MCMC sampling efficiency:" )
                if (PD%SpecBase%ParallelizationModel%isMultiChain) then
                    call PD%Decor%write(PD%LogFile%unit,0,1,1, "+INFINITY" )
                else
                    call PD%Decor%write(PD%LogFile%unit,0,1,1, "UNDEFINED" )
                end if

#if defined CAF_ENABLED || defined MPI_ENABLED
                if (PD%Image%count==1_IK) then
                    msg = PD%name// " is being used in parallel mode but with only one processor.\n&
                                    &This is computationally inefficient.\n&
                                    &Consider using the serial version of the code or provide more processes at runtime."
                    call PD%note( prefix     = PD%brand         &
                                , outputUnit = PD%LogFile%unit  &
                                , newline    = "\n"             &
                                , marginTop  = 0_IK             &
                                , marginBot  = 1_IK             &
                                , msg        = msg              )
                    call PD%note( prefix     = PD%brand         &
                                , outputUnit = output_unit      &
                                , newline    = "\n"             &
                                , marginTop  = 3_IK             &
                                , marginBot  = 0_IK             &
                                , msg        = msg              )
                end if

            else

                block

                    use Statistics_mod, only: getGeoPDF
                    integer(IK) :: imageCount, maxSpeedupImageCount, lenGeoPDF
                    real(RK)    :: seqSecTime, parSecTime, comSecTime, serialTime, avgCommTimePerFunCallPerNode
                    real(RK)    :: speedup, currentSpeedup, maxSpeedup
                    real(RK)    :: firstImageWeight
                    real(RK)    , allocatable :: GeoPDF(:)
                    character(:), allocatable :: formatIn

                    GeoPDF = getGeoPDF(successProb=mcmcSamplingEfficiency,minSeqLen=PD%Image%count)
                    lenGeoPDF = size(GeoPDF)

                    ! compute the serial and sequential runtime of the code per function call

                    seqSecTime = epsilon(seqSecTime) ! time cost of the sequential section of the code, which is negligible here
                    serialTime = PD%Stats%avgTimePerFunCalInSec + seqSecTime ! serial runtime of the code per function call

                    ! compute the communication overhead for each additional image beyond master

                    avgCommTimePerFunCallPerNode = PD%Stats%avgCommTimePerFunCall / real(PD%Image%count-1_IK,kind=RK)

                    !***************************************************************************************************************
                    ! compute the optimal parallelism effciency with the current MCMC sampling efficiency
                    !***************************************************************************************************************

                    ! compute fraction of points sampled by the first image

                    maxSpeedup = 1._RK
                    maxSpeedupImageCount = 1_IK
                    loopOptimalImageCount: do imageCount = 2, max(PD%Image%count,lenGeoPDF)
                        firstImageWeight = sum(GeoPDF(1:lenGeoPDF:imageCount))
                        parSecTime = PD%Stats%avgTimePerFunCalInSec * firstImageWeight ! parallel-section runtime of the code per function call
                        comSecTime = (imageCount-1_IK) * avgCommTimePerFunCallPerNode  ! assumption: communication time grows linearly with the number of nodes
                        speedup = serialTime / (seqSecTime+parSecTime+comSecTime)
                        maxSpeedup = max( maxSpeedup , speedup )
                        if (maxSpeedup==speedup) maxSpeedupImageCount = imageCount
                        if (imageCount==PD%Image%count) then
                            currentSpeedup = speedup
                            if (maxSpeedup/=speedup) exit loopOptimalImageCount
                        end if
                    end do loopOptimalImageCount

                    call PD%Decor%write(PD%LogFile%unit,1,0,1, "Estimated maximum speedup gained via "//PD%SpecBase%ParallelizationModel%val//" parallelization model compared to serial mode:" )
                    call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( currentSpeedup ) )

                    call PD%Decor%write(PD%LogFile%unit,1,0,1, "Predicted optimal number of physical computing processes for "//PD%SpecBase%ParallelizationModel%val//" parallelization model, given the current MCMC sampling efficiency:" )
                    call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( maxSpeedupImageCount ) )

                    call PD%Decor%write(PD%LogFile%unit,1,0,1, "Predicted optimal maximum speedup gained via "//PD%SpecBase%ParallelizationModel%val//" parallelization model, given the current MCMC sampling efficiency:" )
                    call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( maxSpeedup ) )

                    if (currentSpeedup<1._RK) then
                        formatIn = "(g0.6)"
                        msg =   "The time cost of calling the user-provided objective function must be at least " // &
                                num2str(1._RK/currentSpeedup,formatIn) // " times more (that is, ~" // &
                                num2str(10**6*PD%Stats%avgTimePerFunCalInSec/currentSpeedup,formatIn) // &
                                " microseconds) to see any performance benefits from " // &
                                PD%SpecBase%ParallelizationModel%val // " parallelization model for this simulation. "
                        if (maxSpeedupImageCount==1_IK) then
                            msg = msg// "Consider simulating in the serial mode in the future (if used within &
                                        &the same computing environment and with the same configuration as used here)."
                        else
                            msg = msg// "Consider simulating on " // num2str(maxSpeedupImageCount) // " processors in the future &
                                        &(if used within the same computing environment and with the same configuration as used here)."
                        end if
                        call PD%note( prefix     = PD%brand         &
                                    , outputUnit = PD%LogFile%unit  &
                                    , newline    = NLC              &
                                    , marginTop  = 0_IK             &
                                    , marginBot  = 1_IK             &
                                    , msg        = msg              )
                        call PD%note( prefix     = PD%brand         &
                                    , outputUnit = output_unit      &
                                    , newline    = NLC              &
                                    , marginTop  = 3_IK             &
                                    , marginBot  = 0_IK             &
                                    , msg        = msg              )
                    end if

                    !***************************************************************************************************************
                    ! compute the absolute optimal parallelism effciency under any MCMC sampling efficiency
                    !***************************************************************************************************************

                    imageCount = 1_IK
                    maxSpeedup = 1._RK
                    maxSpeedupImageCount = 1_IK
                    loopAbsoluteOptimalImageCount: do 
                        imageCount = imageCount + 1_IK
                        parSecTime = PD%Stats%avgTimePerFunCalInSec / imageCount
                        comSecTime = (imageCount-1_IK) * avgCommTimePerFunCallPerNode  ! assumption: communication time grows linearly with the number of nodes
                        speedup = serialTime / (seqSecTime+parSecTime+comSecTime)
                        if (maxSpeedup>speedup) then
                            maxSpeedupImageCount = imageCount - 1_IK
                            exit loopAbsoluteOptimalImageCount
                        end if
                        maxSpeedup = max( maxSpeedup , speedup )
                    end do loopAbsoluteOptimalImageCount

                    call PD%Decor%write(PD%LogFile%unit,1,0,1, "Predicted absolute optimal number of physical computing processes for "//PD%SpecBase%ParallelizationModel%val//" parallelization model, under any MCMC sampling efficiency:" )
                    call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( maxSpeedupImageCount ) )

                    call PD%Decor%write(PD%LogFile%unit,1,0,1, "Predicted absolute optimal maximum speedup gained via "//PD%SpecBase%ParallelizationModel%val//" parallelization model, under any MCMC sampling efficiency:" )
                    call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( maxSpeedup ) )

                    msg = "This simulation will likely NOT benefit from any additional computing processor beyond the predicted absolute optimal number, "//num2str(maxSpeedupImageCount) &
                        //", in the above. This is true for any value of MCMC sampling efficiency."
                    call PD%note( prefix     = PD%brand         &
                                , outputUnit = PD%LogFile%unit  &
                                , newline    = NLC              &
                                , marginTop  = 0_IK             &
                                , marginBot  = 1_IK             &
                                , msg        = msg              )

                end block
#endif
            end if

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Burnin location in the compact chain, based on the occurrence likelihood:" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( PD%Stats%BurninLoc%compact ) )

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Burnin location in the compact chain, based on the value of burninAdaptationMeasure:" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( PD%Stats%AdaptationBurninLoc%compact ) )

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Burnin location in the verbose (Markov) chain, based on the occurrence likelihood:" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( PD%Stats%BurninLoc%verbose ) )

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Burnin location in the verbose (Markov) chain, based on the value of burninAdaptationMeasure:" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( PD%Stats%AdaptationBurninLoc%verbose ) )

            ! reset BurninLoc to the maximum value

            if (PD%Stats%AdaptationBurninLoc%compact>PD%Stats%BurninLoc%compact) then
                PD%Stats%BurninLoc%compact = PD%Stats%AdaptationBurninLoc%compact
                PD%Stats%BurninLoc%verbose = PD%Stats%AdaptationBurninLoc%verbose
            end if

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Location of the first occurrence of maximum-logFunc in the compact chain:" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( PD%Stats%LogFuncMode%Loc%compact ) )

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Location of the first occurrence of maximum-logFunc in the verbose (Markov) chain:" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( PD%Stats%LogFuncMode%Loc%verbose ) )

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Maximum-logFunc value (maximum of the user-specified objective function):" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( PD%Stats%LogFuncMode%val ) )

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Maximum-logFunc coordinates (mode of the user-specified objective function):" )
            PD%SpecBase%VariableNameList%MaxLen%val = max( PD%SpecBase%VariableNameList%MaxLen%val , PD%SpecBase%OutputColumnWidth%val )
            PD%SpecBase%VariableNameList%MaxLen%val = max( PD%SpecBase%VariableNameList%MaxLen%val , PD%SpecBase%OutputRealPrecision%val + 7 )
            PD%SpecBase%VariableNameList%MaxLen%str = num2str(PD%SpecBase%VariableNameList%MaxLen%val+1)
            formatStr = "(*(A"//logFileColWidthStr//"))" !formatStr = "(*(g0,:,', '))"
            write(PD%LogFile%unit,formatStr) (trim(adjustl(PD%SpecBase%VariableNameList%Val(i))),i=1,ndim)
            formatStr = "(*(E"//logFileColWidthStr//"."//PD%SpecBase%OutputRealPrecision%str//"))"
            write(PD%LogFile%unit,formatStr) PD%Stats%LogFuncMode%Crd
            call PD%Decor%write(PD%LogFile%unit,0,1)

            !***********************************************************************************************************************
            ! Compute the MCMC chain's statistical properties
            !***********************************************************************************************************************

            if (PD%Image%isFirst) then
                call PD%note( prefix        = PD%brand          &
                            , outputUnit    = output_unit       &
                            , newline       = NLC               &
                            , marginTop     = 3_IK              &
                            , msg           = "Computing the Markov chain's statistical properties..." )
            end if

            call PD%Decor%writeDecoratedText( text = "\nMarkov chain's statistical properties\n" &
                                            , marginTop = 1_IK  &
                                            , marginBot = 1_IK  &
                                            , newline = "\n"    &
                                            , outputUnit = PD%LogFile%unit )

            PD%Stats%Chain%count = sum(PD%Chain%Weight(PD%Stats%BurninLoc%compact:PD%Chain%count%compact))

            ! compute the covariance and correlation upper-triangle matrices

            if (allocated(PD%Stats%Chain%Mean)) deallocate(PD%Stats%Chain%Mean); allocate(PD%Stats%Chain%Mean(ndim))
            if (allocated(PD%Stats%Chain%CovMat)) deallocate(PD%Stats%Chain%CovMat); allocate(PD%Stats%Chain%CovMat(ndim,ndim))
            call getWeiSamCovUppMeanTrans   ( np = PD%Chain%count%compact - PD%Stats%BurninLoc%compact + 1_IK &
                                            , sumWeight = PD%Stats%Chain%count &
                                            , nd = ndim &
                                            , Point = PD%Chain%State(1:ndim,PD%Stats%BurninLoc%compact:PD%Chain%count%compact) &
                                            , Weight = PD%Chain%Weight(PD%Stats%BurninLoc%compact:PD%Chain%count%compact) &
                                            , CovMatUpper = PD%Stats%Chain%CovMat &
                                            , Mean = PD%Stats%Chain%Mean &
                                            )

            PD%Stats%Chain%CorMat = getUpperCorMatFromUpperCovMat(nd=ndim,CovMatUpper=PD%Stats%Chain%CovMat)

            ! transpose the covariance and correlation matrices

            do i = 1, ndim
                PD%Stats%Chain%CorMat(i,i) = 1._RK
                PD%Stats%Chain%CorMat(i+1:ndim,i) = PD%Stats%Chain%CorMat(i,i+1:ndim)
                PD%Stats%Chain%CovMat(i+1:ndim,i) = PD%Stats%Chain%CovMat(i,i+1:ndim)
            end do

            ! compute the quantiles

            if (allocated(PD%Stats%Chain%Quantile)) deallocate(PD%Stats%Chain%Quantile)
            allocate(PD%Stats%Chain%Quantile(QPROB%count,ndim))
            do i = 1, ndim
                PD%Stats%Chain%Quantile(1:QPROB%count,i) = getQuantile  ( np = PD%Chain%count%compact - PD%Stats%BurninLoc%compact + 1_IK &
                                                                        , nq = QPROB%count &
                                                                        , SortedQuantileProbability = QPROB%Value &
                                                                        , Point = PD%Chain%State(i,PD%Stats%BurninLoc%compact:PD%Chain%count%compact) &
                                                                        , Weight = PD%Chain%Weight(PD%Stats%BurninLoc%compact:PD%Chain%count%compact) &
                                                                        , sumWeight = PD%Stats%Chain%count &
                                                                        )
            end do

            ! report the MCMC chain statistics

            formatStr = "(1A" // logFileColWidthStr // ",*(E" // logFileColWidthStr // "." // PD%SpecBase%OutputRealPrecision%str // "))"

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "length of the Markov Chain excluding burnin:" )
            call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str(PD%Stats%Chain%count) )

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Mean and standard deviation of the Markov chain variables:")
            write(PD%LogFile%unit, "(*(A"//logFileColWidthStr//"))") "", "Mean", "Standard Deviation"
            do i = 1, ndim
                write( PD%LogFile%unit , formatStr ) trim(adjustl(PD%SpecBase%VariableNameList%Val(i))), PD%Stats%Chain%Mean(i), sqrt(PD%Stats%Chain%CovMat(i,i))
            end do
            call PD%Decor%write(PD%LogFile%unit,0,1)

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Covariance matrix of the Markov chain:")
            write(PD%LogFile%unit, "(*(A"//logFileColWidthStr//"))") "", (trim(adjustl(PD%SpecBase%VariableNameList%Val(i))),i=1,ndim)
            do i = 1, ndim
                write( PD%LogFile%unit , formatStr ) trim(adjustl(PD%SpecBase%VariableNameList%Val(i))), PD%Stats%Chain%CovMat(1:ndim,i)
            end do
            call PD%Decor%write(PD%LogFile%unit,0,1)

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Correlation matrix of the Markov chain:")
            write(PD%LogFile%unit, "(*(A"//logFileColWidthStr//"))") "", (trim(adjustl(PD%SpecBase%VariableNameList%Val(i))),i=1,ndim)
            do i = 1, ndim
                write( PD%LogFile%unit , formatStr ) trim(adjustl(PD%SpecBase%VariableNameList%Val(i))), PD%Stats%Chain%CorMat(1:ndim,i)
            end do
            call PD%Decor%write(PD%LogFile%unit,0,1)

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Quantiles of the Markov chain variables:")
            write(PD%LogFile%unit, "(*(A"//logFileColWidthStr//"))") "Quantile", (trim(adjustl(PD%SpecBase%VariableNameList%Val(i))),i=1,ndim)
            do iq = 1, QPROB%count
                write( PD%LogFile%unit , formatStr ) trim(adjustl(QPROB%Name(iq))), (PD%Stats%Chain%Quantile(iq,i),i=1,ndim)
            end do
            call PD%Decor%write(PD%LogFile%unit,0,1)

            !***********************************************************************************************************************
            ! Generate the i.i.d. sample statistics and output file (if requested)
            !***********************************************************************************************************************

            ! report refined sample statistics, and generate output refined sample if requested.

            if (PD%Image%isFirst) then
                call PD%note( prefix        = PD%brand          &
                            , outputUnit    = output_unit       &
                            , newline       = NLC               &
                            , msg           = "Computing the final decorrelated sample size..." )
            end if

            call PD%RefinedChain%get( CFC                       = PD%Chain                                  &
                                    , Err                       = PD%Err                                    &
                                    , burninLoc                 = PD%Stats%BurninLoc%compact                &
                                    , sampleRefinementCount     = PD%SpecMCMC%SampleRefinementCount%val     &
                                    , sampleRefinementMethod    = PD%SpecMCMC%SampleRefinementMethod%val    &
                                    )

            if (PD%Err%occurred) then
                PD%Err%msg = PROCEDURE_NAME // PD%Err%msg
                call PD%abort( Err = PD%Err, prefix = PD%brand, newline = NLC, outputUnit = PD%LogFile%unit )
                return
            end if

            ! compute the maximum integrated autocorrelation times for each variable

            formatStr = "(2I" // logFileColWidthStr // ",*(E" // logFileColWidthStr // "." // PD%SpecBase%OutputRealPrecision%str // "))"

            call PD%Decor%write(PD%LogFile%unit,1,0,1, "Integrated Autocorrelation (IAC) of the Markov chain:")
            write(PD%LogFile%unit, "(*(A"//logFileColWidthStr//"))") "RefinementStage","SampleSize","IAC(SampleLogFunc)",("IAC("//trim(adjustl(PD%SpecBase%VariableNameList%Val(i)))//")",i=1,ndim)
            do i = 0, PD%RefinedChain%numRefinement
                write( PD%LogFile%unit , formatStr ) i, PD%RefinedChain%Count(i)%Verbose, PD%RefinedChain%IAC(0:ndim,i)
            end do
            call PD%Decor%write(PD%LogFile%unit,0,1)

            if (PD%RefinedChain%numRefinement==0_IK) then
                call PD%warn( prefix = PD%brand &
                            , newline = "\n" &
                            , outputUnit = PD%LogFile%unit &
                            , msg = "The user-specified sampleRefinementCount ("// num2str(PD%SpecMCMC%SampleRefinementCount%val) // ") &
                                    &is too small to ensure an accurate computation of the decorrelated i.i.d. effective sample size. &
                                    &No refinement of the Markov chain was performed." &
                            )
            end if

            ! report the final Effective Sample Size (ESS) based on IAC

            blockEffectiveSampleSize: associate( effectiveSampleSize => sum(PD%RefinedChain%Weight(1:PD%RefinedChain%Count(PD%RefinedChain%numRefinement)%compact)) )
                call PD%Decor%write(PD%LogFile%unit,1,0,1, "Estimated Effective (decorrelated) Sample Size (ESS):" )
                call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( effectiveSampleSize ) )
    
                if (PD%SpecDRAM%DelayedRejectionCount%val==0_IK) then
                    call PD%Decor%write(PD%LogFile%unit,1,0,1, "Effective sampling efficiency ( = effectiveSampleSize / acceptedPlusRejectedFunctionCalls ):" )
                    call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( real(effectiveSampleSize,kind=RK) / real(PD%Stats%NumFunCall%acceptedRejected,kind=RK) ) )
                else
                    call PD%Decor%write(PD%LogFile%unit,1,0,1, "Effective sampling efficiency (including delayed rejections) [ = effectiveSampleSize / acceptedPlusRejectedPlusDelayedRejectionFunctionCalls ]:" )
                    call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str( real(effectiveSampleSize,kind=RK) / real(PD%Stats%NumFunCall%acceptedRejectedDelayed,kind=RK) ) )
                end if
            end associate blockEffectiveSampleSize

            ! generate output refined sample if requested

            blockSampleFileGeneration: if (PD%SpecBase%SampleSize%val==0_IK) then

                call PD%note( prefix        = PD%brand          &
                            , outputUnit    = PD%LogFile%unit   &
                            , newline       = "\n"              &
                            , msg           = "Skipping the generation of the decorrelated sample and output file, as requested by the user..." )

            else blockSampleFileGeneration

                !*******************************************************************************************************************
                ! sample file generation report
                !*******************************************************************************************************************

                ! report to the report-file(s)

                call PD%note( prefix        = PD%brand          &
                            , outputUnit    = PD%LogFile%unit   &
                            , newline       = NLC               &
                            , msg           = "Generating the output " // PD%SampleFile%suffix // " file:"//NLC // PD%SampleFile%Path%original )


                if (PD%Image%isFirst) then 

                    ! print the message for the generating the output sample file on the first image

                    call PD%note( prefix     = PD%brand         &
                                , outputUnit = output_unit      &
                                , newline    = NLC              &
                                , marginBot  = 0_IK             &
                                , msg        = "Generating the output " // PD%SampleFile%suffix // " file:" )

                    call PD%note( prefix     = PD%brand         &
                                , outputUnit = output_unit      &
                                , newline    = NLC              &
                                , marginTop  = 0_IK             &
                                , marginBot  = 0_IK             &
                                , msg = PD%SampleFile%Path%original )

                    ! print the message for the generating the output sample file on the rest of the images in order

#if defined CAF_ENABLED || defined MPI_ENABLED
                    if (PD%SpecBase%ParallelizationModel%isMultiChain) then
                        block
                            use String_mod, only: replaceStr, num2str
                            integer(IK) :: imageID
                            do imageID = 2, PD%Image%count
                                call PD%note( prefix = PD%brand         &
                                            , outputUnit = output_unit  &
                                            , newline = NLC             &
                                            , marginTop = 0_IK          &
                                            , marginBot = 0_IK          &
                                            , msg = replaceStr( string = PD%SampleFile%Path%original, search = "process_1", substitute = "process_"//num2str(imageID) ) )
                            end do
                        end block
                    end if
#endif
                    call PD%Decor%write(output_unit,0,1)

                end if

                !*******************************************************************************************************************
                ! begin sample file generation
                !*******************************************************************************************************************

                if (PD%SpecBase%SampleSize%val/=-1) then

                    if (PD%SpecBase%SampleSize%val<0_IK) PD%SpecBase%SampleSize%val = abs(PD%SpecBase%SampleSize%val) * PD%RefinedChain%Count(PD%RefinedChain%numRefinement)%verbose

                    if (PD%SpecBase%SampleSize%val/=0_IK) then
                        if (PD%SpecBase%SampleSize%val<PD%RefinedChain%Count(PD%RefinedChain%numRefinement)%verbose) then
                            call PD%warn    ( prefix = PD%brand &
                                            , newline = "\n" &
                                            , marginTop = 0_IK &
                                            , outputUnit = PD%LogFile%unit &
                                            , msg = "The user-specified sample size ("// num2str(PD%SpecBase%SampleSize%val) // ") &
                                                    &is smaller than the potentially-optimal i.i.d. sample size &
                                                    &(" // num2str(PD%RefinedChain%Count(PD%RefinedChain%numRefinement)%verbose) // "). &
                                                    &The output sample contains i.i.d. samples, however, the sample-size &
                                                    &could have been larger if it had been set to the optimal size. &
                                                    &To get the optimal size in the future runs, set sampleSize = -1, or drop&
                                                    &it from the input list." &
                                            )
                        elseif (PD%SpecBase%SampleSize%val>PD%RefinedChain%Count(PD%RefinedChain%numRefinement)%verbose) then
                            call PD%warn    ( prefix = PD%brand &
                                            , newline = "\n" &
                                            , marginTop = 0_IK &
                                            , outputUnit = PD%LogFile%unit &
                                            , msg = "The user-specified sample size ("// num2str(PD%SpecBase%SampleSize%val) // ") &
                                                    &is larger than the potentially-optimal i.i.d. sample size &
                                                    &(" // num2str(PD%RefinedChain%Count(PD%RefinedChain%numRefinement)%verbose) // "). &
                                                    &The resulting sample likely contains duplicates and is not independently &
                                                    &and identically distributed (i.i.d.).\nTo get the optimal &
                                                    &size in the future runs, set sampleSize = -1, or drop &
                                                    &it from the input list." &
                                            )
                        else
                            call PD%warn    ( prefix = PD%brand &
                                            , newline = "\n" &
                                            , marginTop = 0_IK &
                                            , outputUnit = PD%LogFile%unit &
                                            , msg = "How lucky that could be! &
                                                    &The user-specified sample size (" // num2str(PD%SpecBase%SampleSize%val) // ") &
                                                    &is equal to the potentially-optimal i.i.d. sample size determined by "//PD%name//"." &
                                            )
                        end if
                    end if

                    ! regenerate the refined sample, this time with the user-specified sample size.

                    call PD%RefinedChain%get( CFC                       = PD%Chain &
                                            , Err                       = PD%Err &
                                            , burninLoc                 = PD%Stats%BurninLoc%compact &
                                            , refinedChainSize          = PD%SpecBase%SampleSize%val &
                                            , sampleRefinementCount     = PD%SpecMCMC%SampleRefinementCount%val     &
                                            , sampleRefinementMethod    = PD%SpecMCMC%SampleRefinementMethod%val    &
                                            )
                    if (PD%Err%occurred) then
                        PD%Err%msg = PROCEDURE_NAME // PD%Err%msg
                        call PD%abort( Err = PD%Err, prefix = PD%brand, newline = NLC, outputUnit = PD%LogFile%unit )
                        return
                    end if

                end if

                ! open the output sample file

                PD%SampleFile%unit = 3000001 + PD%Image%id  ! for some unknown reason, if newunit is used, GFortran opens the file as an internal file
                open( unit      = PD%SampleFile%unit           &
                    , file      = PD%SampleFile%Path%original  &
                    , status    = PD%SampleFile%status         &
                    , iostat    = PD%SampleFile%Err%stat       &
                    , position  = PD%SampleFile%Position%value )
                PD%Err = PD%SampleFile%getOpenErr(PD%SampleFile%Err%stat)
                if (PD%Err%occurred) then
                    PD%Err%msg = PROCEDURE_NAME//": Error occurred while opening "//PD%name//" "//PD%SampleFile%suffix//" file='"//PD%SampleFile%Path%original//"'."
                    call PD%abort( Err = PD%Err, prefix = PD%brand, newline = NLC, outputUnit = PD%LogFile%unit )
                    return
                end if

                ! determine the sample file's contents' format

                if (PD%SpecBase%OutputColumnWidth%val>0) then
                    formatStr = "(*(E"//PD%SpecBase%OutputColumnWidth%str//"."//PD%SpecBase%OutputRealPrecision%str//",:,'"//PD%SpecBase%OutputDelimiter%val//"'))"
                else
                    formatStr = PD%SampleFile%format
                end if

                ! write to the output sample file

                call PD%RefinedChain%write(PD%SampleFile%unit,PD%SampleFile%format,formatStr)
                close(PD%SampleFile%unit)

                !***********************************************************************************************************************
                ! Compute the refined sample's statistical properties
                !***********************************************************************************************************************

                if (PD%Image%isFirst) then
                    call PD%note( prefix        = PD%brand          &
                                , outputUnit    = output_unit       &
                                , newline       = NLC               &
                                , marginTop     = 2_IK              &
                                , msg           = "Computing the output sample's statistical properties..." )
                end if

                call PD%Decor%writeDecoratedText( text = "\nOutput sample's statistical properties\n" &
                                                , marginTop = 1_IK  &
                                                , marginBot = 1_IK  &
                                                , newline = "\n"    &
                                                , outputUnit = PD%LogFile%unit )

                PD%Stats%Sample%count = PD%RefinedChain%Count(PD%RefinedChain%numRefinement)%verbose

                ! compute the covariance and correlation upper-triangle matrices

                if (allocated(PD%Stats%Sample%Mean)) deallocate(PD%Stats%Sample%Mean); allocate(PD%Stats%Sample%Mean(ndim))
                if (allocated(PD%Stats%Sample%CovMat)) deallocate(PD%Stats%Sample%CovMat); allocate(PD%Stats%Sample%CovMat(ndim,ndim))
                call getWeiSamCovUppMeanTrans   ( np = PD%RefinedChain%Count(PD%RefinedChain%numRefinement)%compact &
                                                , sumWeight = PD%RefinedChain%Count(PD%RefinedChain%numRefinement)%verbose &
                                                , nd = ndim &
                                                , Point = PD%RefinedChain%LogFuncState(1:ndim,1:PD%RefinedChain%Count(PD%RefinedChain%numRefinement)%compact) &
                                                , Weight = PD%RefinedChain%Weight(1:PD%RefinedChain%Count(PD%RefinedChain%numRefinement)%compact) &
                                                , CovMatUpper = PD%Stats%Sample%CovMat &
                                                , Mean = PD%Stats%Sample%Mean &
                                                )

                PD%Stats%Sample%CorMat = getUpperCorMatFromUpperCovMat(nd=ndim,CovMatUpper=PD%Stats%Sample%CovMat)

                ! transpose the covariance and correlation matrices

                do i = 1, ndim
                    PD%Stats%Sample%CorMat(i,i) = 1._RK
                    PD%Stats%Sample%CorMat(i+1:ndim,i) = PD%Stats%Sample%CorMat(i,i+1:ndim)
                    PD%Stats%Sample%CovMat(i+1:ndim,i) = PD%Stats%Sample%CovMat(i,i+1:ndim)
                end do

                ! compute the quantiles

                if (allocated(PD%Stats%Sample%Quantile)) deallocate(PD%Stats%Sample%Quantile)
                allocate(PD%Stats%Sample%Quantile(QPROB%count,ndim))
                do i = 1, ndim
                    PD%Stats%Sample%Quantile(1:QPROB%count,i) = getQuantile ( np = PD%RefinedChain%Count(PD%RefinedChain%numRefinement)%compact &
                                                                            , nq = QPROB%count &
                                                                            , SortedQuantileProbability = QPROB%Value &
                                                                            , Point = PD%RefinedChain%LogFuncState(i,1:PD%RefinedChain%Count(PD%RefinedChain%numRefinement)%compact) &
                                                                            , Weight = PD%RefinedChain%Weight(1:PD%RefinedChain%Count(PD%RefinedChain%numRefinement)%compact) &
                                                                            , sumWeight = PD%RefinedChain%Count(PD%RefinedChain%numRefinement)%verbose &
                                                                            )
                end do

                ! report the MCMC chain statistics

                formatStr = "(1A" // logFileColWidthStr // ",*(E" // logFileColWidthStr // "." // PD%SpecBase%OutputRealPrecision%str // "))"

                call PD%Decor%write(PD%LogFile%unit,1,0,1, "Final output sample size:" )
                call PD%Decor%write(PD%LogFile%unit,0,1,1, num2str(PD%Stats%Sample%count) )

                call PD%Decor%write(PD%LogFile%unit,1,0,1, "Mean and standard deviation of the output sample:")
                write(PD%LogFile%unit, "(*(A"//logFileColWidthStr//"))") "", "Mean", "Standard Deviation"
                do i = 1, ndim
                    write( PD%LogFile%unit , formatStr ) trim(adjustl(PD%SpecBase%VariableNameList%Val(i))), PD%Stats%Sample%Mean(i), sqrt(PD%Stats%Sample%CovMat(i,i))
                end do
                call PD%Decor%write(PD%LogFile%unit,0,1)

                call PD%Decor%write(PD%LogFile%unit,1,0,1, "Covariance matrix of the output sample:")
                write(PD%LogFile%unit, "(*(A"//logFileColWidthStr//"))") "", (trim(adjustl(PD%SpecBase%VariableNameList%Val(i))),i=1,ndim)
                do i = 1, ndim
                    write( PD%LogFile%unit , formatStr ) trim(adjustl(PD%SpecBase%VariableNameList%Val(i))), PD%Stats%Sample%CovMat(1:ndim,i)
                end do
                call PD%Decor%write(PD%LogFile%unit,0,1)

                call PD%Decor%write(PD%LogFile%unit,1,0,1, "Correlation matrix of the output sample:")
                write(PD%LogFile%unit, "(*(A"//logFileColWidthStr//"))") "", (trim(adjustl(PD%SpecBase%VariableNameList%Val(i))),i=1,ndim)
                do i = 1, ndim
                    write( PD%LogFile%unit , formatStr ) trim(adjustl(PD%SpecBase%VariableNameList%Val(i))), PD%Stats%Sample%CorMat(1:ndim,i)
                end do
                call PD%Decor%write(PD%LogFile%unit,0,1)

                call PD%Decor%write(PD%LogFile%unit,1,0,1, "Quantiles of the output sample's variables:")
                write(PD%LogFile%unit, "(*(A"//logFileColWidthStr//"))") "Quantile", (trim(adjustl(PD%SpecBase%VariableNameList%Val(i))),i=1,ndim)
                do iq = 1, QPROB%count
                    write( PD%LogFile%unit , formatStr ) trim(adjustl(QPROB%Name(iq))), (PD%Stats%Sample%Quantile(iq,i),i=1,ndim)
                end do
                call PD%Decor%write(PD%LogFile%unit,0,1)

                !***********************************************************************************************************************
                ! Begin inter-chain convergence test in multiChain parallelization mode
                !***********************************************************************************************************************

#if defined CAF_ENABLED || defined MPI_ENABLED

                blockInterChainConvergence: if (PD%SpecBase%ParallelizationModel%isMultiChain .and. PD%Image%count>1_IK) then

                    call PD%note( prefix            = PD%brand          &
                                , outputUnit        = PD%LogFile%unit   &
                                , newline           = NLC               &
                                , marginTop         = 1_IK              &
                                , marginBot         = 1_IK              &
                                , msg               = "Computing the inter-chain convergence probabilities..." )
                    if (PD%Image%isFirst) then ! print the message on stdout
                        call PD%note( prefix        = PD%brand          &
                                    , outputUnit    = output_unit       &
                                    , newline       = NLC               &
                                    , marginTop     = 1_IK              &
                                    , marginBot     = 1_IK              &
                                    , msg           = "Computing the inter-chain convergence probabilities..." )
                    end if

#if defined CAF_ENABLED
                    sync all
#elif defined MPI_ENABLED
                    block
                        use mpi
                        integer :: ierrMPI
                        call mpi_barrier(mpi_comm_world,ierrMPI)
                    end block
#endif

                    ! read the sample files generated by other images

                    multiChainConvergenceTest: block

                        use Sort_mod, only: sortAscending
                        use String_mod, only: replaceStr
                        use Statistics_mod, only: doSortedKS2
                        use ParaDRAMRefinedChain_mod, only: readRefinedChain
                        type(RefinedChain_type)     :: RefinedChainThisImage, RefinedChainThatImage
                        integer(IK)                 :: imageID,indexMinProbKS,imageMinProbKS
                        real(RK)                    :: statKS, minProbKS
                        real(RK)    , allocatable   :: ProbKS(:)
                        character(:), allocatable   :: inputSamplePath

                        minProbKS = huge(minProbKS)
                        allocate(ProbKS(0:ndim))

                        ! read the refined chain on the current image

                        RefinedChainThisImage = readRefinedChain( sampleFilePath=PD%SampleFile%Path%original, delimiter=PD%SpecBase%OutputDelimiter%val, ndim=ndim )
                        if (RefinedChainThisImage%Err%occurred) then
                            PD%Err%msg = PROCEDURE_NAME//RefinedChainThisImage%Err%msg
                            call PD%abort( Err = PD%Err, prefix = PD%brand, newline = NLC, outputUnit = PD%LogFile%unit )
                            return
                        end if

                        ! sort the refined chain on the current image

                        do i = 0, ndim
                            call sortAscending  ( np = RefinedChainThisImage%Count(RefinedChainThisImage%numRefinement)%verbose &
                                                , Point = RefinedChainThisImage%LogFuncState(1:RefinedChainThisImage%Count(RefinedChainThisImage%numRefinement)%verbose,i) &
                                                )
                        end do

                        ! compute and report the KS convergence probabilities for all images

                        formatStr = "(1I" // logFileColWidthStr // ",*(E" // logFileColWidthStr // "." // PD%SpecBase%OutputRealPrecision%str // "))"

                        call PD%Decor%write(PD%LogFile%unit,1,0,1, "Pairwise inter-chain Kolmogorov-Smirnov (KS) convergence (similarity) probabilities:" )
                        write(PD%LogFile%unit, "(*(A"//logFileColWidthStr//"))") "ProcessID",("probKS("//trim(adjustl(PD%RefinedChain%ColHeader(i)%record))//")",i=0,ndim)
                        do imageID = 1, PD%Image%count

                            if (imageID/=PD%Image%id) then

                                ! read the refined chain on the other image

                                inputSamplePath = replaceStr( string = PD%SampleFile%Path%original &
                                                            , search = "process_"//num2str(PD%Image%id) &
                                                            , substitute = "process_"//num2str(imageID) )

                                RefinedChainThatImage = readRefinedChain( sampleFilePath=inputSamplePath, delimiter=PD%SpecBase%OutputDelimiter%val, ndim=ndim )
                                if (RefinedChainThatImage%Err%occurred) then
                                    PD%Err%msg = PROCEDURE_NAME//RefinedChainThatImage%Err%msg
                                    call PD%abort( Err = PD%Err, prefix = PD%brand, newline = NLC, outputUnit = PD%LogFile%unit )
                                    return
                                end if

                                do i = 0, ndim

                                    ! sort the refined chain on the other image

                                    call sortAscending  ( np = RefinedChainThatImage%Count(RefinedChainThatImage%numRefinement)%verbose &
                                                        , Point = RefinedChainThatImage%LogFuncState(1:RefinedChainThatImage%Count(RefinedChainThatImage%numRefinement)%verbose,i) &
                                                        )

                                    ! compute the inter-chain KS probability table

                                    call doSortedKS2( np1 = RefinedChainThisImage%Count(RefinedChainThisImage%numRefinement)%verbose &
                                                    , np2 = RefinedChainThatImage%Count(RefinedChainThatImage%numRefinement)%verbose &
                                                    , SortedPoint1 = RefinedChainThisImage%LogFuncState(1:RefinedChainThisImage%Count(RefinedChainThisImage%numRefinement)%verbose,i) &
                                                    , SortedPoint2 = RefinedChainThatImage%LogFuncState(1:RefinedChainThatImage%Count(RefinedChainThatImage%numRefinement)%verbose,i) &
                                                    , statKS = statKS &
                                                    , probKS = ProbKS(i) &
                                                    )

                                    ! determine the minimum KS

                                    if (ProbKS(i)<minProbKS) then
                                        minProbKS = ProbKS(i)
                                        indexMinProbKS = i
                                        imageMinProbKS = imageID
                                    end if

                                end do

                                ! write the inter-chain KS probability table row

                                write( PD%LogFile%unit , formatStr ) imageID, (ProbKS(i), i = 0, ndim)

                            end if

                        end do

                        call PD%note( prefix = PD%brand, outputUnit = PD%LogFile%unit, newline = "\n", msg = "Higher KS probabilities are better, indicating less evidence for a lack of convergence." )

                        ! write the smallest KS probabilities

                        call PD%Decor%write(PD%LogFile%unit,0,0,1, "Smallest KS probability for the inter-chain sampling convergence:" )
                        call PD%Decor%write(PD%LogFile%unit,0,1,1,  num2str(minProbKS)//" for "//PD%RefinedChain%ColHeader(indexMinProbKS)%record//&
                                                                    " on the chains generated by processes "//num2str(PD%Image%id)//&
                                                                    " and "//num2str(imageMinProbKS)//"." )

                        ! report the smallest KS probabilities on stdout

                        if (PD%Image%isFirst) call PD%note  ( prefix     = PD%brand         &
                                                            , outputUnit = output_unit      &
                                                            , newline    = NLC              &
                                                            , marginTop  = 2_IK             &
                                                            , marginBot  = 0_IK             &
                                                            , msg        = "The smallest KS probabilities for the inter-chain sampling convergence:" )

                        do imageID = 1, PD%Image%count
                            if (imageID==PD%Image%id) then
                                call PD%note( prefix     = PD%brand         &
                                            , outputUnit = output_unit      &
                                            , newline    = NLC              &
                                            , marginTop  = 0_IK             &
                                            , marginBot  = 0_IK             &
                                            , msg        = num2str(minProbKS)//" for "//PD%RefinedChain%ColHeader(indexMinProbKS)%record//&
                                                           " on the chains generated by processes "//num2str(PD%Image%id)//&
                                                           " and "//num2str(imageMinProbKS)//"." )
                            end if
#if defined CAF_ENABLED
                            sync all
#elif defined MPI_ENABLED
                            block
                                use mpi
                                integer :: ierrMPI
                                call mpi_barrier(mpi_comm_world,ierrMPI)
                            end block
#endif
                        end do

                    end block multiChainConvergenceTest

                end if blockInterChainConvergence
#endif

                !***********************************************************************************************************************
                ! End inter-chain convergence test in multiChain parallelization mode
                !***********************************************************************************************************************

            end if blockSampleFileGeneration

            !***********************************************************************************************************************
            ! End of generating the i.i.d. sample statistics and output file (if requested)
            !***********************************************************************************************************************

            call PD%Decor%write(PD%LogFile%unit,1,1)

            ! Mission accomplished.

            call PD%note( prefix = PD%brand, outputUnit = PD%LogFile%unit, newline = "\n", msg = "Mission Accomplished." )
            if (output_unit/=PD%LogFile%unit .and. PD%Image%isFirst) then
                call PD%Decor%write(output_unit,1,1)
                call PD%note( prefix = PD%brand, outputUnit = output_unit, newline = "\n", msg = "Mission Accomplished." )
                call PD%Decor%write(output_unit,1,1)
            end if

            close(PD%TimeFile%unit)
            close(PD%LogFile%unit)

        end if blockMasterPostProcessing

        !nullify(PD%Proposal)

#if defined CAF_ENABLED
        sync all
#elif defined MPI_ENABLED
        block
            use mpi
            integer :: ierrMPI
            call mpi_barrier(mpi_comm_world,ierrMPI)
            if (PD%SpecBase%MpiFinalizeRequested%val) then
                call mpi_finalize(ierrMPI)
            end if
        end block
#endif

    end subroutine runSampler

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    module subroutine getSpecFromInputFile(PD,nd)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getSpecFromInputFile
#endif
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str

        ! ParaMonte namelist variables
        use SpecBase_SampleSize_mod                         , only: sampleSize
        use SpecBase_RandomSeed_mod                         , only: randomSeed
        use SpecBase_Description_mod                        , only: description
        use SpecBase_OutputFileName_mod                     , only: outputFileName
        use SpecBase_OutputDelimiter_mod                    , only: outputDelimiter
        use SpecBase_ChainFileFormat_mod                    , only: chainFileFormat
        use SpecBase_VariableNameList_mod                   , only: variableNameList
        use SpecBase_RestartFileFormat_mod                  , only: restartFileFormat
        use SpecBase_OutputColumnWidth_mod                  , only: outputColumnWidth
        use SpecBase_OutputRealPrecision_mod                , only: outputRealPrecision
        use SpecBase_SilentModeRequested_mod                , only: silentModeRequested
        use SpecBase_DomainLowerLimitVec_mod                , only: domainLowerLimitVec
        use SpecBase_DomainUpperLimitVec_mod                , only: domainUpperLimitVec
        use SpecBase_ParallelizationModel_mod               , only: parallelizationModel
        use SpecBase_InputFileHasPriority_mod               , only: inputFileHasPriority
        use SpecBase_ProgressReportPeriod_mod               , only: progressReportPeriod
        use SpecBase_TargetAcceptanceRate_mod               , only: TargetAcceptanceRate
        use SpecBase_MpiFinalizeRequested_mod               , only: mpiFinalizeRequested
        use SpecBase_MaxNumDomainCheckToWarn_mod            , only: maxNumDomainCheckToWarn
        use SpecBase_MaxNumDomainCheckToStop_mod            , only: maxNumDomainCheckToStop
#if defined CFI_ENABLED
        use SpecBase_InterfaceType_mod                      , only: interfaceType
#endif

        ! ParaMCMC namelist variables
        use SpecMCMC_ChainSize_mod                          , only: chainSize
        use SpecMCMC_ScaleFactor_mod                        , only: scaleFactor
        use SpecMCMC_StartPointVec_mod                      , only: startPointVec
        use SpecMCMC_ProposalModel_mod                      , only: proposalModel
        use SpecMCMC_ProposalStartCovMat_mod                , only: ProposalStartCovMat
        use SpecMCMC_ProposalStartCorMat_mod                , only: ProposalStartCorMat
        use SpecMCMC_ProposalStartStdVec_mod                , only: ProposalStartStdVec
        use SpecMCMC_SampleRefinementCount_mod              , only: sampleRefinementCount
        use SpecMCMC_sampleRefinementMethod_mod             , only: sampleRefinementMethod
        use SpecMCMC_RandomStartPointRequested_mod          , only: randomStartPointRequested
        use SpecMCMC_RandomStartPointDomainLowerLimitVec_mod, only: randomStartPointDomainLowerLimitVec
        use SpecMCMC_RandomStartPointDomainUpperLimitVec_mod, only: randomStartPointDomainUpperLimitVec

        ! ParaDRAM namelist variables
        use SpecDRAM_AdaptiveUpdateCount_mod                , only: adaptiveUpdateCount
        use SpecDRAM_AdaptiveUpdatePeriod_mod               , only: adaptiveUpdatePeriod
        use SpecDRAM_greedyAdaptationCount_mod              , only: greedyAdaptationCount
        use SpecDRAM_DelayedRejectionCount_mod              , only: delayedRejectionCount
        use SpecDRAM_BurninAdaptationMeasure_mod            , only: burninAdaptationMeasure
        use SpecDRAM_delayedRejectionScaleFactorVec_mod     , only: delayedRejectionScaleFactorVec

        implicit none
        class(ParaDRAM_type), intent(inout) :: PD
        integer(IK), intent(in)             :: nd
        character(*), parameter             :: PROCEDURE_NAME = SUBMODULE_NAME//"@getSpecFromInputFile()"

        ! ParaMonte variables
        namelist /ParaDRAM/ sampleSize
        namelist /ParaDRAM/ randomSeed
        namelist /ParaDRAM/ description
        namelist /ParaDRAM/ outputFileName
        namelist /ParaDRAM/ outputDelimiter
        namelist /ParaDRAM/ ChainFileFormat
        namelist /ParaDRAM/ variableNameList
        namelist /ParaDRAM/ restartFileFormat
        namelist /ParaDRAM/ outputColumnWidth
        namelist /ParaDRAM/ outputRealPrecision
        namelist /ParaDRAM/ silentModeRequested
        namelist /ParaDRAM/ domainLowerLimitVec
        namelist /ParaDRAM/ domainUpperLimitVec
        namelist /ParaDRAM/ parallelizationModel
        namelist /ParaDRAM/ progressReportPeriod
        namelist /ParaDRAM/ inputFileHasPriority
        namelist /ParaDRAM/ TargetAcceptanceRate
        namelist /ParaDRAM/ mpiFinalizeRequested
        namelist /ParaDRAM/ maxNumDomainCheckToWarn
        namelist /ParaDRAM/ maxNumDomainCheckToStop
#if defined CFI_ENABLED
        namelist /ParaDRAM/ interfaceType
#endif

        ! ParaMCMC variables
        namelist /ParaDRAM/ chainSize
        namelist /ParaDRAM/ scaleFactor
        namelist /ParaDRAM/ startPointVec
        namelist /ParaDRAM/ proposalModel
        namelist /ParaDRAM/ proposalStartCovMat
        namelist /ParaDRAM/ proposalStartCorMat
        namelist /ParaDRAM/ proposalStartStdVec
        namelist /ParaDRAM/ sampleRefinementCount
        namelist /ParaDRAM/ sampleRefinementMethod
        namelist /ParaDRAM/ randomStartPointRequested
        namelist /ParaDRAM/ randomStartPointDomainLowerLimitVec
        namelist /ParaDRAM/ randomStartPointDomainUpperLimitVec

        ! ParaDRAM variables
        namelist /ParaDRAM/ adaptiveUpdateCount
        namelist /ParaDRAM/ adaptiveUpdatePeriod
        namelist /ParaDRAM/ greedyAdaptationCount
        namelist /ParaDRAM/ delayedRejectionCount
        namelist /ParaDRAM/ burninAdaptationMeasure
        namelist /ParaDRAM/ delayedRejectionScaleFactorVec

        ! initialize/nullify all general input options

        call PD%SpecBase%nullifyNameListVar(nd)
        call PD%SpecMCMC%nullifyNameListVar(nd)
        call PD%SpecDRAM%nullifyNameListVar(nd)

        ! read input options if input file is provided

        blockReadInputFile: if (PD%inputFileArgIsPresent) then

            blockInputFileType: if (PD%InputFile%isInternal) then

                ! read input file as an internal file
#if defined DBG_ENABLED
                read(PD%InputFile%Path%original,nml=ParaDRAM)
#else
                read(PD%InputFile%Path%original,nml=ParaDRAM,iostat=PD%InputFile%Err%stat)
                PD%Err = PD%InputFile%getReadErr(PD%InputFile%Err%stat,PD%InputFile%Path%modified)
                if (PD%Err%occurred) then
                    if (is_iostat_end(PD%Err%stat)) then
                        call PD%warnUserAboutMissingNamelist(PD%brand,PD%name,"ParaDRAM",PD%LogFile%unit)
                    else
                        PD%Err%msg = PROCEDURE_NAME // PD%Err%msg
                        return
                    end if
                end if
#endif

            else blockInputFileType ! the input file is external

                ! close input file if it is open

                if (PD%InputFile%isOpen) then
                    close(unit=PD%InputFile%unit,iostat=PD%InputFile%Err%stat)
                    PD%Err = PD%InputFile%getCloseErr(PD%InputFile%Err%stat)
                    if (PD%Err%occurred) then
                        PD%Err%msg =    PROCEDURE_NAME // ": Error occurred while attempting to close the user-provided input file='" // &
                                        PD%InputFile%Path%modified // "', unit=" // num2str(PD%InputFile%unit) // ".\n" // &
                                        PD%Err%msg
                        return
                    end if
                end if

                ! open input file

                open( newunit = PD%InputFile%unit       &
                    , file = PD%InputFile%Path%modified &
                    , status = PD%InputFile%status      &
                    , iostat = PD%InputFile%Err%stat    &
                    )
                PD%Err = PD%InputFile%getOpenErr(PD%InputFile%Err%stat)
                if (PD%Err%occurred) then
                    PD%Err%msg =    PROCEDURE_NAME // ": Error occurred while attempting to open the user-provided input file='" // &
                                    PD%InputFile%Path%modified // "', unit=" // num2str(PD%InputFile%unit) // ".\n" // &
                                    PD%Err%msg
                    return
                end if

                ! read input file

#if defined DBG_ENABLED
                read(PD%InputFile%unit,nml=ParaDRAM)
#else
                read(PD%InputFile%unit,nml=ParaDRAM,iostat=PD%InputFile%Err%stat)
                PD%Err = PD%InputFile%getReadErr(PD%InputFile%Err%stat,PD%InputFile%Path%modified)
                if (PD%Err%occurred) then
                    if (is_iostat_end(PD%Err%stat)) then
                        call PD%warnUserAboutMissingNamelist(PD%brand,PD%name,"ParaDRAM",PD%LogFile%unit)
                    else
                        PD%Err%msg = PROCEDURE_NAME // PD%Err%msg
                        return
                    end if
                end if
#endif
                ! close input file

                close(unit=PD%InputFile%unit,iostat=PD%InputFile%Err%stat)
                PD%Err = PD%InputFile%getCloseErr(PD%InputFile%Err%stat)
                if (PD%Err%occurred) then
                    PD%Err%msg =    PROCEDURE_NAME // ": Error occurred while attempting to close the user-provided input file='" // &
                                    PD%InputFile%Path%modified // "', unit=" // num2str(PD%InputFile%unit) // ".\n" // &
                                    PD%Err%msg
                    return
                end if

            end if blockInputFileType

        end if blockReadInputFile

        ! setup SpecBase variables that have been read form the input file

        call PD%SpecBase%setFromInputFile( Err = PD%Err )
        if (PD%Err%occurred) then
            PD%Err%msg = PROCEDURE_NAME // PD%Err%msg
            return
        end if

        ! setup SpecMCMC variables that have been read form the input file

        call PD%SpecMCMC%setFromInputFile( Err = PD%Err                                                 &
                                         , nd = nd                                                      &
                                         , domainLowerLimitVec = PD%SpecBase%DomainLowerLimitVec%Val    &
                                         , domainUpperLimitVec = PD%SpecBase%DomainUpperLimitVec%Val    )
        if (PD%Err%occurred) then
            PD%Err%msg = PROCEDURE_NAME // PD%Err%msg
            return
        end if

        ! setup SpecDRAM variables that have been read form the input file

        call PD%SpecDRAM%setFromInputFile( PD%Err )
        if (PD%Err%occurred) then
            PD%Err%msg = PROCEDURE_NAME // PD%Err%msg
            return
        end if

  end subroutine getSpecFromInputFile

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end submodule Setup_smod