!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (ParaDRAM_mod) Setup_smod

    !use Constants_mod, only: IK, RK ! gfortran 9.3 compile crashes with this line
    implicit none

    character(*), parameter :: SUBMODULE_NAME = MODULE_NAME // "@Setup_smod"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
                                , targetAcceptanceRate                  &
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
                                , proposalStartStdVec                   &
                                , proposalStartCorMat                   &
                                , proposalStartCovMat                   &
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
        use Decoration_mod, only: GENERIC_OUTPUT_FORMAT
        use Decoration_mod, only: GENERIC_TABBED_FORMAT
        use Decoration_mod, only: getGenericFormat
        use Decoration_mod, only: INDENT
        use Statistics_mod, only: getUpperCorMatFromUpperCovMat
        use Statistics_mod, only: getWeiSamCovUppMeanTrans
        use Statistics_mod, only: getQuantile
        use ParaMonte_mod, only: QPROB
        use Constants_mod, only: IK, RK, NLC, PMSM, UNDEFINED
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
        real(RK)    , intent(in), optional  :: targetAcceptanceRate(2)
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
        real(RK)    , intent(in), optional  :: proposalStartStdVec(ndim)
        real(RK)    , intent(in), optional  :: proposalStartCorMat(ndim,ndim)
        real(RK)    , intent(in), optional  :: proposalStartCovMat(ndim,ndim)
        integer(IK) , intent(in), optional  :: adaptiveUpdateCount
        integer(IK) , intent(in), optional  :: adaptiveUpdatePeriod
        integer(IK) , intent(in), optional  :: greedyAdaptationCount
        integer(IK) , intent(in), optional  :: delayedRejectionCount
        real(RK)    , intent(in), optional  :: burninAdaptationMeasure
        real(RK)    , intent(in), optional  :: delayedRejectionScaleFactorVec(:)

        character(*), parameter             :: PROCEDURE_NAME = SUBMODULE_NAME // "@runSampler()"

       !type(ProposalSymmetric_type), target   :: ProposalSymmetric
        integer(IK)                         :: i, iq
        character(:), allocatable           :: msg, formatStr, formatStrInt, formatStrReal, formatAllReal
        real(RK)                            :: mcmcSamplingEfficiency

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Initialize SpecBase variables, then check existence of inputFile and open it and return the unit file, if it exists
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call PD%setupParaMonte  ( nd = ndim               &
                                , name = PMSM%ParaDRAM    &
                               !, date = "May 23 2018"    &
                               !, version = "1.0.0"       &
                                , inputFile = inputFile   &
                                )
        if (PD%Err%occurred) return

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Initialize ParaMCMC variables
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call PD%setupParaMCMC()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Initialize ParaDRAM variables
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        PD%SpecDRAM = SpecDRAM_type( PD%nd%val, PD%name ) ! , PD%SpecMCMC%ChainSize%def )

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! read variables from input file if it exists
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call PD%getSpecFromInputFile(ndim)
        if (PD%Err%occurred) then
            PD%Err%msg = PROCEDURE_NAME // PD%Err%msg
            call PD%abort( Err = PD%Err, prefix = PD%brand, newline = "\n", outputUnit = PD%LogFile%unit )
            return
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! read variables from argument list if needed
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
                                                , targetAcceptanceRate                  = targetAcceptanceRate                  &
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
                                                , proposalStartStdVec                   = proposalStartStdVec                   &
                                                , proposalStartCorMat                   = proposalStartCorMat                   &
                                                , proposalStartCovMat                   = proposalStartCovMat                   &
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

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Now depending on the requested parallelism type, determine the process master/slave type and open output files
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! setup log file and open output files
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call PD%setupOutputFiles()
        if (PD%Err%occurred) then
            PD%Err%msg = PROCEDURE_NAME // PD%Err%msg
            call PD%abort( Err = PD%Err, prefix = PD%brand, newline = NLC, outputUnit = PD%LogFile%unit )
            return
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! report variable values to the log file
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (PD%isFreshRun) then
            call PD%SpecBase%reportValues(prefix=PD%brand,outputUnit=PD%LogFile%unit,isMasterImage=PD%Image%isMaster)
            call PD%SpecMCMC%reportValues(prefix=PD%brand,outputUnit=PD%LogFile%unit,isMasterImage=PD%Image%isMaster,methodName=PD%name,splashModeRequested=PD%SpecBase%SilentModeRequested%isFalse)
            call PD%SpecDRAM%reportValues(prefix=PD%brand,outputUnit=PD%LogFile%unit,isMasterImage=PD%Image%isMaster,splashModeRequested=PD%SpecBase%SilentModeRequested%isFalse)
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! perform variable value sanity checks
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! setup output files
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        ! g0 edit descriptor format is slightly more arbitrary and compiler-dependent.

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

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! setup the proposal distribution
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! run ParaDRAM kernel
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (PD%isFreshRun .and. PD%Image%isMaster) then
            call PD%Decor%writeDecoratedText( text = "\nStarting " // PD%name // " sampling - " // getNiceDateTime() // "\n" &
                                            , marginTop = 1_IK  &
                                            , marginBot = 1_IK  &
                                            , newline = "\n"    &
                                            , outputUnit = PD%LogFile%unit )
        end if

        call PD%runKernel( getLogFunc = getLogFunc )
        if(PD%Err%occurred) return ! relevant only for MATLAB

        if (PD%isFreshRun .and. PD%Image%isMaster) then
            call PD%Decor%writeDecoratedText( text = "\nExiting " // PD%name // " sampling - " // getNiceDateTime() // "\n" &
                                            , marginTop = 1_IK  &
                                            , marginBot = 1_IK  &
                                            , newline = "\n"    &
                                            , outputUnit = PD%LogFile%unit )
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! start ParaDRAM post-processing
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        blockMasterPostProcessing: if (PD%Image%isMaster) then

            formatStrInt  = "('" // INDENT // "',1A" // PD%LogFile%maxColWidth%str // ",*(I" // PD%LogFile%maxColWidth%str // "))"
            formatStrReal = "('" // INDENT // "',1A" // PD%LogFile%maxColWidth%str // ",*(E" // PD%LogFile%maxColWidth%str // "." // PD%SpecBase%OutputRealPrecision%str // "))"

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.numFuncCall.accepted"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Stats%NumFunCall%accepted
            msg = "This is the total number of accepted function calls (unique samples)."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.numFuncCall.acceptedRejected"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Stats%NumFunCall%acceptedRejected
            msg = "This is the total number of accepted or rejected function calls."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.numFuncCall.acceptedRejectedDelayed"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Stats%NumFunCall%acceptedRejectedDelayed
            msg = "This is the total number of accepted or rejected or delayed-rejection (if any requested) function calls."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.numFuncCall.acceptedRejectedDelayedUnused"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Stats%NumFunCall%acceptedRejectedDelayedUnused
            msg = "This is the total number of accepted or rejected or unused function calls (by all processes, including delayed rejections, if any requested)."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.efficiency.meanAcceptanceRate"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Chain%MeanAccRate(PD%Stats%NumFunCall%accepted)
            msg = "This is the average MCMC acceptance rate."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            mcmcSamplingEfficiency = real(PD%Stats%NumFunCall%accepted,kind=RK) / real(PD%Stats%NumFunCall%acceptedRejected,kind=RK)

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.efficiency.acceptedOverAcceptedRejected"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) mcmcSamplingEfficiency
            msg =   "This is the MCMC sampling efficiency given the accepted and rejected function calls, that is, &
                    &the number of accepted function calls divided by the number of (accepted + rejected) function calls."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.efficiency.acceptedOverAcceptedRejectedDelayed"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) real(PD%Stats%NumFunCall%accepted,kind=RK) / real(PD%Stats%NumFunCall%acceptedRejectedDelayed,kind=RK)
            msg =   "This is the MCMC sampling efficiency given the accepted, rejected, and delayed-rejection (if any requested) function calls, that is, &
                    &the number of accepted function calls divided by the number of (accepted + rejected + delayed-rejection) function calls."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.efficiency.acceptedOverAcceptedRejectedDelayedUnused"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) real(PD%Stats%NumFunCall%accepted,kind=RK) / real(PD%Stats%NumFunCall%acceptedRejectedDelayedUnused,kind=RK)
            msg =   "This is the MCMC sampling efficiency given the accepted, rejected, delayed-rejection (if any requested), and unused function calls, that is, &
                    &the number of accepted function calls divided by the number of (accepted + rejected + delayed-rejection + unused) function calls."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.time.total"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Timer%Time%total
            msg = "This is the total runtime in seconds."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.time.perFuncCallAccepted"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Timer%Time%total / real(PD%Stats%NumFunCall%accepted,kind=RK)
            msg = "This is the average effective time cost of each accepted function call, in seconds."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.time.perFuncCallAcceptedRejected"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Timer%Time%total / real(PD%Stats%NumFunCall%acceptedRejected,kind=RK)
            msg = "This is the average effective time cost of each accepted or rejected function call, in seconds."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.time.perFuncCallAcceptedRejectedDelayed"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Timer%Time%total / real(PD%Stats%NumFunCall%acceptedRejectedDelayed,kind=RK)
            msg = "This is the average effective time cost of each accepted or rejected function call (including delayed-rejections, if any requested), in seconds."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.time.perFuncCallAcceptedRejectedDelayedUnused"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Timer%Time%total / real(PD%Stats%NumFunCall%acceptedRejectedDelayedUnused,kind=RK)
            msg = "This is the average effective time cost of each accepted or rejected or unused function call (including delayed-rejections, if any requested), in seconds."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if (PD%Image%count==1_IK) then
                msg = UNDEFINED
            else
                msg = num2str( PD%Stats%avgCommTimePerFunCall )
            end if

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.time.perInterProcessCommunication"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) msg
            msg = "This is the average time cost of inter-process communications per used (accepted or rejected or delayed-rejection) function call, in seconds."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.time.perFuncCall"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Stats%avgTimePerFunCalInSec
            msg = "This is the average pure time cost of each function call, in seconds."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.current.numProcess"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Image%count
            msg = "This is the number of processes (images) used in this simulation."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if (PD%Image%count==1_IK .or. PD%SpecBase%ParallelizationModel%isMultiChain) then

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.current.speedup"
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Image%count
                msg = "This is the estimated maximum speedup gained via "//PD%SpecBase%ParallelizationModel%val//" parallelization model compared to serial mode."
                call PD%note( prefix     = PD%brand         &
                            , outputUnit = PD%LogFile%unit  &
                            , newline    = NLC              &
                            , marginTop  = 1_IK             &
                            , marginBot  = 1_IK             &
                            , msg        = msg              )

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if (PD%SpecBase%ParallelizationModel%isMultiChain) then
                    msg = "+INFINITY"
                else
                    msg = UNDEFINED
                end if

                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.optimal.current.numProcess"
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) msg
                msg = "This is the predicted optimal number of physical computing processes for "//PD%SpecBase%ParallelizationModel%val//" parallelization model, given the current MCMC sampling efficiency."
                call PD%note( prefix     = PD%brand         &
                            , outputUnit = PD%LogFile%unit  &
                            , newline    = NLC              &
                            , marginTop  = 1_IK             &
                            , marginBot  = 1_IK             &
                            , msg        = msg              )

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if (PD%SpecBase%ParallelizationModel%isMultiChain) then
                    msg = "+INFINITY"
                else
                    msg = UNDEFINED
                end if

                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.optimal.current.speedup"
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) msg
                msg = "This is the predicted optimal maximum speedup gained via "//PD%SpecBase%ParallelizationModel%val//" parallelization model, given the current MCMC sampling efficiency."
                call PD%note( prefix     = PD%brand         &
                            , outputUnit = PD%LogFile%unit  &
                            , newline    = NLC              &
                            , marginTop  = 1_IK             &
                            , marginBot  = 1_IK             &
                            , msg        = msg              )

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if (PD%SpecBase%ParallelizationModel%isMultiChain) then
                    msg = "+INFINITY"
                else
                    msg = UNDEFINED
                end if

                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.optimal.absolute.numProcess"
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) msg
                msg = "This is the predicted absolute optimal number of physical computing processes for "//PD%SpecBase%ParallelizationModel%val//" parallelization model, under any MCMC sampling efficiency."
                call PD%note( prefix     = PD%brand         &
                            , outputUnit = PD%LogFile%unit  &
                            , newline    = NLC              &
                            , marginTop  = 1_IK             &
                            , marginBot  = 1_IK             &
                            , msg        = msg              )

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if (PD%SpecBase%ParallelizationModel%isMultiChain) then
                    msg = "+INFINITY"
                else
                    msg = UNDEFINED
                end if

                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.optimal.absolute.speedup"
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) msg
                msg = "This is the predicted absolute optimal maximum speedup gained via "//PD%SpecBase%ParallelizationModel%val//" parallelization model, under any MCMC sampling efficiency."
#if defined CAF_ENABLED || defined MPI_ENABLED
                if (PD%Image%count==1_IK) then
                    msg = msg//NLC//PD%name//   " is being used in parallel mode but with only one processor.\n&
                                                &This is computationally inefficient.\n&
                                                &Consider using the serial version of the code or provide more processes at runtime."
                    call PD%note( prefix     = PD%brand         &
                                , outputUnit = output_unit      &
                                , newline    = "\n"             &
                                , marginTop  = 3_IK             &
                                , marginBot  = 0_IK             &
                                , msg        = msg              )
                end if
#endif
                call PD%note( prefix     = PD%brand         &
                            , outputUnit = PD%LogFile%unit  &
                            , newline    = NLC              &
                            , marginTop  = 1_IK             &
                            , marginBot  = 1_IK             &
                            , msg        = msg              )

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined CAF_ENABLED || defined MPI_ENABLED
            else

                block

                    use Statistics_mod, only: getgeoPDF
                    logical     :: maxSpeedupFound
                    integer(IK) :: imageCount, maxSpeedupImageCount, lengeoPDF
                    real(RK)    :: seqSecTime, parSecTime, comSecTime, serialTime, avgCommTimePerFunCallPerNode
                    real(RK)    :: currentSpeedup, maxSpeedup
                    real(RK)    :: speedup, firstImageWeight
                    real(RK)    , allocatable :: geoPDF(:), speedupVec(:), tempVec(:)
                    character(:), allocatable :: formatIn, formatScaling

                    formatScaling = "('" // INDENT // "',10(E" // PD%LogFile%maxColWidth%str // "." // PD%SpecBase%OutputRealPrecision%str // "))" ! ,:,','

                    geoPDF = getgeoPDF(successProb=mcmcSamplingEfficiency,minSeqLen=10*PD%Image%count)
                    lengeoPDF = size(geoPDF)

                    ! compute the serial and sequential runtime of the code per function call

                    seqSecTime = epsilon(seqSecTime) ! time cost of the sequential section of the code, which is negligible here
                    serialTime = PD%Stats%avgTimePerFunCalInSec + seqSecTime ! serial runtime of the code per function call

                    ! compute the communication overhead for each additional image beyond master

                    avgCommTimePerFunCallPerNode = PD%Stats%avgCommTimePerFunCall / real(PD%Image%count-1_IK,kind=RK)

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ! compute the optimal parallelism efficiency with the current MCMC sampling efficiency
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    ! compute fraction of points sampled by the first image

                    maxSpeedup = 1._RK
                    maxSpeedupImageCount = 1_IK
                    if (allocated(speedupVec)) deallocate(speedupVec); allocate(speedupVec(lengeoPDF))
                    speedupVec(1) = 1._RK
                    loopOptimalImageCount: do imageCount = 2, lengeoPDF
                        firstImageWeight = sum(geoPDF(1:lengeoPDF:imageCount))
                        parSecTime = PD%Stats%avgTimePerFunCalInSec * firstImageWeight ! parallel-section runtime of the code per function call
                        comSecTime = (imageCount-1_IK) * avgCommTimePerFunCallPerNode  ! assumption: communication time grows linearly with the number of nodes
                        speedupVec(imageCount) = serialTime / (seqSecTime+parSecTime+comSecTime)
                        maxSpeedup = max( maxSpeedup , speedupVec(imageCount) )
                        if (maxSpeedup==speedupVec(imageCount)) maxSpeedupImageCount = imageCount
                        if (imageCount==PD%Image%count) then
                            currentSpeedup = speedupVec(imageCount)
                            !if (maxSpeedup/=speedupVec(imageCount)) exit loopOptimalImageCount
                        end if
                    end do loopOptimalImageCount

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.current.speedup"
                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                    write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) currentSpeedup
                    msg = "This is the estimated maximum speedup gained via "//PD%SpecBase%ParallelizationModel%val//" parallelization model compared to serial mode."
                    call PD%note( prefix     = PD%brand         &
                                , outputUnit = PD%LogFile%unit  &
                                , newline    = NLC              &
                                , marginTop  = 1_IK             &
                                , marginBot  = 1_IK             &
                                , msg        = msg              )

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.optimal.current.numProcess"
                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                    write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) maxSpeedupImageCount
                    msg = "This is the predicted optimal number of physical computing processes for "//PD%SpecBase%ParallelizationModel%val//" parallelization model, given the current MCMC sampling efficiency."
                    call PD%note( prefix     = PD%brand         &
                                , outputUnit = PD%LogFile%unit  &
                                , newline    = NLC              &
                                , marginTop  = 1_IK             &
                                , marginBot  = 1_IK             &
                                , msg        = msg              )

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.optimal.current.speedup"
                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                    write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) maxSpeedup
                    msg = ""
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
                                    , outputUnit = output_unit      &
                                    , newline    = NLC              &
                                    , marginTop  = 3_IK             &
                                    , marginBot  = 0_IK             &
                                    , msg        = msg              )
                        msg = NLC // msg
                    end if
                    msg = "This is the predicted optimal maximum speedup gained via "//PD%SpecBase%ParallelizationModel%val//" parallelization model, given the current MCMC sampling efficiency." // msg
                    call PD%note( prefix     = PD%brand         &
                                , outputUnit = PD%LogFile%unit  &
                                , newline    = NLC              &
                                , marginTop  = 1_IK             &
                                , marginBot  = 1_IK             &
                                , msg        = msg              )

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.optimal.current.scaling.strong.speedup"
                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                    !write(PD%LogFile%unit,formatStrInt) "numProcess", (imageCount, imageCount = 1, lengeoPDF)
                    !write(PD%LogFile%unit,formatStrReal) "speedup" , (speedupVec(imageCount), imageCount = 1, lengeoPDF)
                    do imageCount = 1, lengeoPDF, 10
                        write(PD%LogFile%unit,formatScaling) speedupVec(imageCount:min(imageCount+9_IK,lengeoPDF))
                    end do
                    msg =   "This is the predicted strong-scaling speedup behavior of the "//PD%SpecBase%ParallelizationModel%val//" parallelization model, &
                            &given the current MCMC sampling efficiency, for increasing numbers of processes, starting from a single process."
                    call PD%note( prefix     = PD%brand         &
                                , outputUnit = PD%LogFile%unit  &
                                , newline    = NLC              &
                                , marginTop  = 1_IK             &
                                , marginBot  = 1_IK             &
                                , msg        = msg              )

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ! compute the absolute optimal parallelism efficiency under any MCMC sampling efficiency
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    imageCount = 1_IK
                    maxSpeedup = 1._RK
                    maxSpeedupImageCount = 1_IK
                    speedupVec(1) = 1._RK
                    maxSpeedupFound = .false.
                    lenGeoPDF = size(speedupVec)
                    loopAbsoluteOptimalImageCount: do 
                        imageCount = imageCount + 1_IK
                        if (imageCount>lenGeoPDF) then
                            if (maxSpeedupFound) exit loopAbsoluteOptimalImageCount
                            lenGeoPDF = 2_IK * lenGeoPDF
                            if (allocated(tempVec)) deallocate(tempVec); allocate(tempVec(lenGeoPDF))
                            tempVec(1:lenGeoPDF/2_IK) = speedupVec
                            call move_alloc(tempVec,speedupVec)
                        end if
                        parSecTime = PD%Stats%avgTimePerFunCalInSec / imageCount
                        comSecTime = (imageCount-1_IK) * avgCommTimePerFunCallPerNode  ! assumption: communication time grows linearly with the number of nodes
                        speedupVec(imageCount) = serialTime / (seqSecTime+parSecTime+comSecTime)
                        if (.not.maxSpeedupFound .and. maxSpeedup>speedupVec(imageCount)) then
                            maxSpeedupImageCount = imageCount - 1_IK
                            maxSpeedupFound = .true.
                            !exit loopAbsoluteOptimalImageCount
                        end if
                        maxSpeedup = max( maxSpeedup , speedupVec(imageCount) )
                    end do loopAbsoluteOptimalImageCount

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.optimal.absolute.numProcess"
                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                    write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) maxSpeedupImageCount
                    msg = "This is the predicted absolute optimal number of physical computing processes for "//PD%SpecBase%ParallelizationModel%val//" parallelization model, under any MCMC sampling efficiency."
                    call PD%note( prefix     = PD%brand         &
                                , outputUnit = PD%LogFile%unit  &
                                , newline    = NLC              &
                                , marginTop  = 1_IK             &
                                , marginBot  = 1_IK             &
                                , msg        = msg              )

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.optimal.absolute.speedup"
                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                    write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) maxSpeedup
                    msg =   "This is the predicted absolute optimal maximum speedup gained via "//PD%SpecBase%ParallelizationModel%val//" parallelization model, under any MCMC sampling efficiency. &
                            &This simulation will likely NOT benefit from any additional computing processors beyond the predicted absolute optimal number, " // num2str(maxSpeedupImageCount) // ", &
                            &in the above. This is true for any value of MCMC sampling efficiency. Keep in mind that the predicted absolute optimal number of processors is just an estimate &
                            &whose accuracy depends on many runtime factors, including the topology of the communication network being used, the number of processors per node, &
                            &and the number of tasks to each processor or node."
                    call PD%note( prefix     = PD%brand         &
                                , outputUnit = PD%LogFile%unit  &
                                , newline    = NLC              &
                                , marginTop  = 1_IK             &
                                , marginBot  = 1_IK             &
                                , msg        = msg              )

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.optimal.absolute.scaling.strong.speedup"
                    write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                    !write(PD%LogFile%unit,formatStrInt) "numProcess", (imageCount, imageCount = 1, lengeoPDF)
                    !write(PD%LogFile%unit,formatStrReal)"speedup"   , (speedupVec(imageCount), imageCount = 1, lengeoPDF)
                    do imageCount = 1, lengeoPDF, 10
                        write(PD%LogFile%unit,formatScaling) speedupVec(imageCount:min(imageCount+9_IK,lengeoPDF))
                    end do
                    msg =   "This is the predicted absolute strong-scaling speedup behavior of the "//PD%SpecBase%ParallelizationModel%val//" parallelization model, &
                            &under any MCMC sampling efficiency, for increasing numbers of processes, starting from a single process."
                    call PD%note( prefix     = PD%brand         &
                                , outputUnit = PD%LogFile%unit  &
                                , newline    = NLC              &
                                , marginTop  = 1_IK             &
                                , marginBot  = 1_IK             &
                                , msg        = msg              )

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                end block
#endif
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.compact.burnin.location.likelihoodBased"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Stats%BurninLoc%compact
            msg = "This is the burnin location in the compact chain, based on the occurrence likelihood."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.compact.burnin.location.adaptationBased"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Stats%AdaptationBurninLoc%compact
            msg = "This is the burnin location in the compact chain, based on the value of burninAdaptationMeasure simulation specification."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.burnin.location.likelihoodBased"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Stats%BurninLoc%verbose
            msg = "This is the burnin location in the verbose (Markov) chain, based on the occurrence likelihood."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.burnin.location.adaptationBased"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Stats%AdaptationBurninLoc%verbose
            msg = "This is the burnin location in the verbose (Markov) chain, based on the value of burninAdaptationMeasure simulation specification."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ! reset BurninLoc to the maximum value

            if (PD%Stats%AdaptationBurninLoc%compact>PD%Stats%BurninLoc%compact) then
                PD%Stats%BurninLoc%compact = PD%Stats%AdaptationBurninLoc%compact
                PD%Stats%BurninLoc%verbose = PD%Stats%AdaptationBurninLoc%verbose
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.logFunc.max"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Stats%LogFuncMode%val
            msg = "This is the maximum logFunc value (the maximum of the user-specified objective function)."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.compact.logFunc.max.location"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Stats%LogFuncMode%Loc%compact
            msg = "This is the location of the first occurrence of the maximum logFunc in the compact chain."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.logFunc.max.location"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Stats%LogFuncMode%Loc%verbose
            msg = "This is the location of the first occurrence of the maximum logFunc in the verbose (Markov) chain."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            formatAllReal = "('" // INDENT // "',*(E" // PD%LogFile%maxColWidth%str // "." // PD%SpecBase%OutputRealPrecision%str // "))"

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.logFunc.max.state"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,PD%LogFile%format) (trim(adjustl(PD%SpecBase%VariableNameList%Val(i))), i=1, ndim)
            write(PD%LogFile%unit,formatAllReal) PD%Stats%LogFuncMode%Crd
            msg = "This is the coordinates, within the domain of the user-specified objective function, where the maximum logFunc occurs."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Compute the statistical properties of the MCMC chain
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if (PD%Image%isFirst) then
                call PD%note( prefix        = PD%brand          &
                            , outputUnit    = output_unit       &
                            , newline       = NLC               &
                            , marginTop     = 3_IK              &
                            , msg           = "Computing the statistical properties of the Markov chain..." )
            end if

            call PD%Decor%writeDecoratedText( text = "\nThe statistical properties of the Markov chain\n" &
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

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.length.burninExcluded"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Stats%Chain%count
            msg = "This is the length of the verbose (Markov) Chain excluding burnin."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.avgStd"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,PD%LogFile%format) "variableName", "Mean", "Standard Deviation"
            do i = 1, ndim
                write(PD%LogFile%unit,formatStrReal) trim(adjustl(PD%SpecBase%VariableNameList%Val(i))), PD%Stats%Chain%Mean(i), sqrt(PD%Stats%Chain%CovMat(i,i))
            end do
            msg = "This is the mean and standard deviation of the verbose (Markov) chain variables."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.covmat"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,PD%LogFile%format) "", (trim(adjustl(PD%SpecBase%VariableNameList%Val(i))),i=1,ndim)
            do i = 1, ndim
                write(PD%LogFile%unit,formatStrReal) trim(adjustl(PD%SpecBase%VariableNameList%Val(i))), PD%Stats%Chain%CovMat(1:ndim,i)
            end do
            msg = "This is the covariance matrix of the verbose (Markov) chain."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.cormat"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,PD%LogFile%format) "", (trim(adjustl(PD%SpecBase%VariableNameList%Val(i))),i=1,ndim)
            do i = 1, ndim
                write(PD%LogFile%unit,formatStrReal) trim(adjustl(PD%SpecBase%VariableNameList%Val(i))), PD%Stats%Chain%CorMat(1:ndim,i)
            end do
            msg = "This is the correlation matrix of the verbose (Markov) chain."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.quantile"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,PD%LogFile%format) "Quantile", (trim(adjustl(PD%SpecBase%VariableNameList%Val(i))),i=1,ndim)
            do iq = 1, QPROB%count
                write(PD%LogFile%unit,formatStrReal) trim(adjustl(QPROB%Name(iq))), (PD%Stats%Chain%Quantile(iq,i),i=1,ndim)
            end do
            msg = "This is the quantiles table of the variables of the verbose (Markov) chain."
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Generate the i.i.d. sample statistics and output file (if requested)
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            formatStr = "('" // INDENT // "',2I" // PD%LogFile%maxColWidth%str // ",*(E" // PD%LogFile%maxColWidth%str // "." // PD%SpecBase%OutputRealPrecision%str // "))"

            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.iac"
            write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(PD%LogFile%unit,PD%LogFile%format) "RefinementStage","SampleSize","IAC_SampleLogFunc",("IAC_"//trim(adjustl(PD%SpecBase%VariableNameList%Val(i))),i=1,ndim)
            do i = 0, PD%RefinedChain%numRefinement
                write(PD%LogFile%unit,formatStr) i, PD%RefinedChain%Count(i)%Verbose, PD%RefinedChain%IAC(0:ndim,i)
            end do
            msg = "This is the table of the Integrated Autocorrelation (IAC) of individual variables in the verbose (Markov) chain, at increasing stages of chain refinements."
            if (PD%RefinedChain%numRefinement==0_IK) then
                msg = msg // NLC // "The user-specified sampleRefinementCount ("// num2str(PD%SpecMCMC%SampleRefinementCount%val) // ") &
                                    &is too small to ensure an accurate computation of the decorrelated i.i.d. effective sample size. &
                                    &No refinement of the Markov chain was performed."
            end if
            call PD%note( prefix     = PD%brand         &
                        , outputUnit = PD%LogFile%unit  &
                        , newline    = NLC              &
                        , marginTop  = 1_IK             &
                        , marginBot  = 1_IK             &
                        , msg        = msg              )

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ! report the final Effective Sample Size (ESS) based on IAC

            blockEffectiveSampleSize: associate( effectiveSampleSize => sum(PD%RefinedChain%Weight(1:PD%RefinedChain%Count(PD%RefinedChain%numRefinement)%compact)) )

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.ess"
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) effectiveSampleSize
                msg = "This is the estimated Effective (decorrelated) Sample Size (ESS) of the final refined chain."
                call PD%note( prefix     = PD%brand         &
                            , outputUnit = PD%LogFile%unit  &
                            , newline    = NLC              &
                            , marginTop  = 1_IK             &
                            , marginBot  = 1_IK             &
                            , msg        = msg              )

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.efficiency.essOverAccepted"
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) real(effectiveSampleSize,kind=RK) / real(PD%Stats%NumFunCall%accepted,kind=RK)
                msg =   "This is the effective MCMC sampling efficiency given the accepted function calls, that is, &
                        &the final refined effective sample size (ESS) divided by the number of accepted function calls."
                call PD%note( prefix     = PD%brand         &
                            , outputUnit = PD%LogFile%unit  &
                            , newline    = NLC              &
                            , marginTop  = 1_IK             &
                            , marginBot  = 1_IK             &
                            , msg        = msg              )

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.efficiency.essOverAcceptedRejected"
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) real(effectiveSampleSize,kind=RK) / real(PD%Stats%NumFunCall%acceptedRejected,kind=RK)
                msg =   "This is the effective MCMC sampling efficiency given the accepted and rejected function calls, that is, &
                        &the final refined effective sample size (ESS) divided by the number of (accepted + rejected) function calls."
                call PD%note( prefix     = PD%brand         &
                            , outputUnit = PD%LogFile%unit  &
                            , newline    = NLC              &
                            , marginTop  = 1_IK             &
                            , marginBot  = 1_IK             &
                            , msg        = msg              )

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.efficiency.essOverAcceptedRejectedDelayed"
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) real(effectiveSampleSize,kind=RK) / real(PD%Stats%NumFunCall%acceptedRejectedDelayed,kind=RK)
                msg =   "This is the effective MCMC sampling efficiency given the accepted, rejected, and delayed-rejection (if any requested) function calls, that is, &
                        &the final refined effective sample size (ESS) divided by the number of (accepted + rejected + delayed-rejection) function calls."
                call PD%note( prefix     = PD%brand         &
                            , outputUnit = PD%LogFile%unit  &
                            , newline    = NLC              &
                            , marginTop  = 1_IK             &
                            , marginBot  = 1_IK             &
                            , msg        = msg              )

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.efficiency.essOverAcceptedRejectedDelayedUnused"
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) real(effectiveSampleSize,kind=RK) / real(PD%Stats%NumFunCall%acceptedRejectedDelayedUnused,kind=RK)
                msg =   "This is the effective MCMC sampling efficiency given the accepted, rejected, delayed-rejection (if any requested), and unused function calls, that is, &
                        &the final refined effective sample size (ESS) divided by the number of (accepted + rejected + delayed-rejection + unused) function calls."
                call PD%note( prefix     = PD%brand         &
                            , outputUnit = PD%LogFile%unit  &
                            , newline    = NLC              &
                            , marginTop  = 1_IK             &
                            , marginBot  = 1_IK             &
                            , msg        = msg              )

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

                if (PD%SpecBase%SampleSize%val/=-1_IK) then

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
                open( unit      = PD%SampleFile%unit            &
                    , file      = PD%SampleFile%Path%original   &
                    , status    = PD%SampleFile%status          &
                    , iostat    = PD%SampleFile%Err%stat        &
#if defined IFORT_ENABLED && defined OS_IS_WINDOWS
                    , SHARED                                    &
#endif
                    , position  = PD%SampleFile%Position%value  )
                PD%Err = PD%SampleFile%getOpenErr(PD%SampleFile%Err%stat)
                if (PD%Err%occurred) then
                    PD%Err%msg = PROCEDURE_NAME//": Error occurred while opening the "//PD%name//" "//PD%SampleFile%suffix//" file='"//PD%SampleFile%Path%original//"'. "
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

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! Compute the statistical properties of the refined sample
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if (PD%Image%isFirst) then
                    call PD%note( prefix        = PD%brand          &
                                , outputUnit    = output_unit       &
                                , newline       = NLC               &
                                , marginTop     = 2_IK              &
                                , msg           = "Computing the statistical properties of the final refined sample..." )
                end if

                call PD%Decor%writeDecoratedText( text = "\nThe statistical properties of the final refined sample \n" &
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

                ! report the refined chain statistics

                !formatStr = "(1A" // PD%LogFile%maxColWidth%str // ",*(E" // PD%LogFile%maxColWidth%str // "." // PD%SpecBase%OutputRealPrecision%str // "))"

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.length"
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) PD%Stats%Sample%count
                msg = "This is the final output refined sample size."
                call PD%note( prefix     = PD%brand         &
                            , outputUnit = PD%LogFile%unit  &
                            , newline    = NLC              &
                            , marginTop  = 1_IK             &
                            , marginBot  = 1_IK             &
                            , msg        = msg              )

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.avgStd"
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,PD%LogFile%format) "variableName", "Mean", "Standard Deviation"
                do i = 1, ndim
                    write(PD%LogFile%unit,formatStrReal) trim(adjustl(PD%SpecBase%VariableNameList%Val(i))), PD%Stats%Sample%Mean(i), sqrt(PD%Stats%Sample%CovMat(i,i))
                end do
                msg = "This is the Mean and standard deviation table of the final output refined sample."
                call PD%note( prefix     = PD%brand         &
                            , outputUnit = PD%LogFile%unit  &
                            , newline    = NLC              &
                            , marginTop  = 1_IK             &
                            , marginBot  = 1_IK             &
                            , msg        = msg              )

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.covmat"
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,PD%LogFile%format) "", (trim(adjustl(PD%SpecBase%VariableNameList%Val(i))),i=1,ndim)
                do i = 1, ndim
                    write(PD%LogFile%unit,formatStrReal) trim(adjustl(PD%SpecBase%VariableNameList%Val(i))), PD%Stats%Sample%CovMat(1:ndim,i)
                end do
                msg = "This is the covariance matrix of the final output refined sample."
                call PD%note( prefix     = PD%brand         &
                            , outputUnit = PD%LogFile%unit  &
                            , newline    = NLC              &
                            , marginTop  = 1_IK             &
                            , marginBot  = 1_IK             &
                            , msg        = msg              )

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.cormat"
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,PD%LogFile%format) "", (trim(adjustl(PD%SpecBase%VariableNameList%Val(i))),i=1,ndim)
                do i = 1, ndim
                    write(PD%LogFile%unit,formatStrReal) trim(adjustl(PD%SpecBase%VariableNameList%Val(i))), PD%Stats%Sample%CorMat(1:ndim,i)
                end do
                msg = "This is the correlation matrix of the final output refined sample."
                call PD%note( prefix     = PD%brand         &
                            , outputUnit = PD%LogFile%unit  &
                            , newline    = NLC              &
                            , marginTop  = 1_IK             &
                            , marginBot  = 1_IK             &
                            , msg        = msg              )

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.quantile"
                write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(PD%LogFile%unit,PD%LogFile%format) "Quantile", (trim(adjustl(PD%SpecBase%VariableNameList%Val(i))),i=1,ndim)
                do iq = 1, QPROB%count
                    write(PD%LogFile%unit,formatStrReal) trim(adjustl(QPROB%Name(iq))), (PD%Stats%Sample%Quantile(iq,i),i=1,ndim)
                end do
                msg = "This is the quantiles table of the variables of the final output refined sample."
                call PD%note( prefix     = PD%brand         &
                            , outputUnit = PD%LogFile%unit  &
                            , newline    = NLC              &
                            , marginTop  = 1_IK             &
                            , marginBot  = 1_IK             &
                            , msg        = msg              )

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! Begin inter-chain convergence test in multiChain parallelization mode
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

                        formatStr = "(1I" // PD%LogFile%maxColWidth%str // ",*(E" // PD%LogFile%maxColWidth%str // "." // PD%SpecBase%OutputRealPrecision%str // "))"

                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                        write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.kstest.prob"
                        write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                        write(PD%LogFile%unit,PD%LogFile%format) "ProcessID",("probKS("//trim(adjustl(PD%RefinedChain%ColHeader(i)%record))//")",i=0,ndim)

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

                                write(PD%LogFile%unit, formatStr) imageID, (ProbKS(i), i = 0, ndim)

                            end if

                        end do
                        msg =   "This is the table pairwise inter-chain Kolmogorov-Smirnov (KS) convergence (similarity) probabilities. &
                                &Higher KS probabilities are better, indicating less evidence for a lack of convergence."
                        call PD%note( prefix     = PD%brand         &
                                    , outputUnit = PD%LogFile%unit  &
                                    , newline    = NLC              &
                                    , marginTop  = 1_IK             &
                                    , marginBot  = 1_IK             &
                                    , msg        = msg              )

                        ! write the smallest KS probabilities

                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                        write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.kstest.prob.min"
                        write(PD%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                        write(PD%LogFile%unit,GENERIC_TABBED_FORMAT) minProbKS
                        msg =   "This is the smallest KS-test probability for the inter-chain sampling convergence, which has happened between "//PD%RefinedChain%ColHeader(indexMinProbKS)%record// &
                                " on the chains generated by processes "//num2str(PD%Image%id)//" and "//num2str(imageMinProbKS)//"."
                        call PD%note( prefix     = PD%brand         &
                                    , outputUnit = PD%LogFile%unit  &
                                    , newline    = NLC              &
                                    , marginTop  = 1_IK             &
                                    , marginBot  = 1_IK             &
                                    , msg        = msg              )

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

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! End inter-chain convergence test in multiChain parallelization mode
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end if blockSampleFileGeneration

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! End of generating the i.i.d. sample statistics and output file (if requested)
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        use SpecBase_targetAcceptanceRate_mod               , only: targetAcceptanceRate
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
        use SpecMCMC_proposalStartCovMat_mod                , only: proposalStartCovMat
        use SpecMCMC_proposalStartCorMat_mod                , only: proposalStartCorMat
        use SpecMCMC_proposalStartStdVec_mod                , only: proposalStartStdVec
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
        namelist /ParaDRAM/ targetAcceptanceRate
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
        namelist /ParaDRAM/ proposalStartStdVec
        namelist /ParaDRAM/ proposalStartCorMat
        namelist /ParaDRAM/ proposalStartCovMat
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
#if defined IFORT_ENABLED && defined OS_IS_WINDOWS
                    , SHARED                            &
#endif
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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule Setup_smod