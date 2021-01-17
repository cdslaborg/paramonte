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
!> This file implements the body of the `Setup_smod` submodules of the ParaMonte `ParaDRAM_mod`, `ParaDISE_mod`, `ParaNest_mod` modules.
!>
!> \remark
!> This module requires preprocessing, prior to compilation.
!>
!> \author Amir Shahmoradi

#if defined PARADRAM

#define ParaXXXX ParaDRAM
#define ParaXXXX_type ParaDRAM_type
#define ParaXXXX_ChainFileContents_mod ParaDRAM_ChainFileContents_mod

#elif defined PARADISE

#define ParaXXXX ParaDISE
#define ParaXXXX_type ParaDISE_type
#define ParaXXXX_ChainFileContents_mod ParaDISE_ChainFileContents_mod

#elif defined PARANEST

#define ParaXXXX ParaNest
#define ParaXXXX_type ParaNest_type
#define ParaXXXX_ChainFileContents_mod ParaNest_ChainFileContents_mod

#else
#error "Unrecognized sampler in ParaXXXX_mod@Setup_mod.inc.f90"
#endif

!submodule (ParaDRAM_mod) Setup_smod
!submodule (ParaDISE_mod) Setup_smod
!submodule (ParaNest_mod) Setup_smod

    !use Constants_mod, only: IK, RK ! gfortran 9.3 compile crashes with this line
    implicit none

    character(*), parameter :: SUBMODULE_NAME = MODULE_NAME // "@Setup_smod"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a method of [ParaDRAM_type](@ref paradram_type) and [ParaDISE_type](@ref paradise_type) classes.
    !> Setup the sampler and run the simulation.
    !>
    !> \remark
    !> There are only two mandatory pieces of information that must be provided when calling this subroutine.
    !> The rest of the parameters can be also set from within an external input file, whose path can be given
    !> as an input argument `inputFile` to this subroutine.
    !>
    !> @param[in]   self        :   An object instantiated from a ParaMonte sampler class.
    !> @param[in]   ndim        :   The number of dimensions of the domain of the objective function to be sampled.
    !> @param[in]   getLogFunc  :   The first point.
    !> @param[in]   inputFile   :   The path to the external input file (**optional**).
    !>
    !> \remark
    !> For detailed descriptions of the rest of the input parameters of this subroutine,
    !> see: https://www.cdslab.org/paramonte/notes/usage/paradram/specifications/
    !>
    !> \remark
    !> This procedure requires preprocessing.
    module subroutine runSampler( self                                  &
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
                                , overwriteRequested                    &
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
#if defined PARADRAM || defined PARADISE
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
#elif defined ParaNest
                                ! ParaNest variables
                                , tightness                             &
                                , tolerance                             &
                                , scaleFactor                           &
                                , proposalModel                         &
                                , liveSampleSize                        &
                                , inclusionFraction                     &
                                , adaptiveUpdateCount                   &
                                , adaptiveUpdatePeriod                  &
                                , mahalSqWeightExponent                 &
                                , stabilizationRequested                &
                                , MaxAllowedKmeansFailure               &
                                , maxAllowedMinVolFailure               &
                                , maxKvolumeLoopRecursion               &
#endif
                                )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: runSampler
#endif

#if defined PARADRAM
        use ParaDRAM_ProposalUniform_mod, only: ProposalUniform_type => Proposal_type
        use ParaDRAM_ProposalNormal_mod, only: ProposalNormal_type => Proposal_type
#elif defined PARADISE
        use ParaDISE_ProposalUniform_mod, only: ProposalUniform_type => Proposal_type
        use ParaDISE_ProposalNormal_mod, only: ProposalNormal_type => Proposal_type
#elif defined PARANEST
        use ParaNestProposalRejection_mod, only: ProposalRejEll_type => Proposal_type
#endif
        use ParaMonteLogFunc_mod, only: getLogFunc_proc
        use Decoration_mod, only: GENERIC_OUTPUT_FORMAT
        use Decoration_mod, only: GENERIC_TABBED_FORMAT
        use Decoration_mod, only: getGenericFormat
        use Decoration_mod, only: INDENT
        use Statistics_mod, only: getCorMatUpperFromCovMatUpper
        use Statistics_mod, only: getWeiSamCovUppMeanTrans
        use Statistics_mod, only: getQuantile
        use Statistics_mod, only: getMean
        use ParaMonte_mod, only: QPROB
        use Constants_mod, only: IK, RK, NLC, PMSM, UNDEFINED, POSINF_RK
        use DateTime_mod, only: getNiceDateTime
        use String_mod, only: num2str
        use Matrix_mod, only: getEye

        implicit none

        ! self
        class(ParaXXXX_type), intent(inout) :: self

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
        logical     , intent(in), optional  :: overwriteRequested
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
#if defined PARADRAM || defined PARADISE
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
#elif defined ParaNest
        ! ParaNest variables
        real(RK)    , intent(in), optional  :: tightness
        real(RK)    , intent(in), optional  :: tolerance
        real(RK)    , intent(in), optional  :: scaleFactor
        character(*), intent(in), optional  :: proposalModel
        integer(IK) , intent(in), optional  :: liveSampleSize
        real(RK)    , intent(in), optional  :: inclusionFraction
        integer(IK) , intent(in), optional  :: adaptiveUpdateCount
        integer(IK) , intent(in), optional  :: adaptiveUpdatePeriod
        real(RK)    , intent(in), optional  :: mahalSqWeightExponent
        logical     , intent(in), optional  :: stabilizationRequested
        integer(IK) , intent(in), optional  :: maxAllowedKmeansFailure
        integer(IK) , intent(in), optional  :: maxAllowedMinVolFailure
        integer(IK) , intent(in), optional  :: maxKvolumeLoopRecursion
#endif

        character(*), parameter             :: PROCEDURE_NAME = SUBMODULE_NAME // "@runSampler()"

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Initialize SpecBase variables, then check existence of inputFile and open it and return the unit file, if it exists
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call self%setupParaMonte( nd = ndim             &
                                , name = PMSM%ParaXXXX  &
                                , inputFile = inputFile &
                                )
#if (defined MPI_ENABLED || defined CAF_ENABLED) && (defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED)
        block; use Err_mod, only: bcastErr; call bcastErr(self%Err); end block
#endif
        if (self%Err%occurred) then
            ! LCOV_EXCL_START
            self%Err%msg = PROCEDURE_NAME // self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
            ! LCOV_EXCL_STOP
        end if

#if defined PARADRAM || defined PARADISE
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Initialize ParaMCMC variables
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call self%setupParaMCMC()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Initialize ParaXXXX variables
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        self%SpecDRAM = SpecDRAM_type( self%nd%val, self%name ) ! , self%SpecMCMC%ChainSize%def )
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! read variables from input file if it exists
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call self%getSpecFromInputFile(ndim)
#if (defined MPI_ENABLED || defined CAF_ENABLED) && (defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED)
        block; use Err_mod, only: bcastErr; call bcastErr(self%Err); end block
#endif
        if (self%Err%occurred) then
            ! LCOV_EXCL_START
            self%Err%msg = PROCEDURE_NAME // self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = "\n", outputUnit = self%LogFile%unit )
            return
            ! LCOV_EXCL_STOP
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! read variables from argument list if needed
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call self%setWarnAboutProcArgHasPriority()
        if (self%procArgNeeded) then
            call self%SpecBase%setFromInputArgs ( Err                                   = self%Err                              & ! LCOV_EXCL_LINE
                                                , sampleSize                            = sampleSize                            & ! LCOV_EXCL_LINE
                                                , randomSeed                            = randomSeed                            & ! LCOV_EXCL_LINE
                                                , description                           = description                           & ! LCOV_EXCL_LINE
                                                , outputFileName                        = outputFileName                        & ! LCOV_EXCL_LINE
                                                , outputDelimiter                       = outputDelimiter                       & ! LCOV_EXCL_LINE
                                                , chainFileFormat                       = chainFileFormat                       & ! LCOV_EXCL_LINE
                                                , variableNameList                      = variableNameList                      & ! LCOV_EXCL_LINE
                                                , restartFileFormat                     = restartFileFormat                     & ! LCOV_EXCL_LINE
                                                , outputColumnWidth                     = outputColumnWidth                     & ! LCOV_EXCL_LINE
                                                , overwriteRequested                    = overwriteRequested                    & ! LCOV_EXCL_LINE
                                                , outputRealPrecision                   = outputRealPrecision                   & ! LCOV_EXCL_LINE
                                                , silentModeRequested                   = silentModeRequested                   & ! LCOV_EXCL_LINE
                                                , domainLowerLimitVec                   = domainLowerLimitVec                   & ! LCOV_EXCL_LINE
                                                , domainUpperLimitVec                   = domainUpperLimitVec                   & ! LCOV_EXCL_LINE
                                                , parallelizationModel                  = parallelizationModel                  & ! LCOV_EXCL_LINE
                                                , progressReportPeriod                  = progressReportPeriod                  & ! LCOV_EXCL_LINE
                                                , targetAcceptanceRate                  = targetAcceptanceRate                  & ! LCOV_EXCL_LINE
                                                , mpiFinalizeRequested                  = mpiFinalizeRequested                  & ! LCOV_EXCL_LINE
                                                , maxNumDomainCheckToWarn               = maxNumDomainCheckToWarn               & ! LCOV_EXCL_LINE
                                                , maxNumDomainCheckToStop               = maxNumDomainCheckToStop               & ! LCOV_EXCL_LINE
                                                )
#if defined PARADRAM || defined PARADISE
            call self%SpecMCMC%setFromInputArgs ( chainSize                             = chainSize                             & ! LCOV_EXCL_LINE
                                                , scaleFactor                           = scaleFactor                           & ! LCOV_EXCL_LINE
                                                , startPointVec                         = startPointVec                         & ! LCOV_EXCL_LINE
                                                , proposalModel                         = proposalModel                         & ! LCOV_EXCL_LINE
                                                , proposalStartStdVec                   = proposalStartStdVec                   & ! LCOV_EXCL_LINE
                                                , proposalStartCorMat                   = proposalStartCorMat                   & ! LCOV_EXCL_LINE
                                                , proposalStartCovMat                   = proposalStartCovMat                   & ! LCOV_EXCL_LINE
                                                , sampleRefinementCount                 = sampleRefinementCount                 & ! LCOV_EXCL_LINE
                                                , sampleRefinementMethod                = sampleRefinementMethod                & ! LCOV_EXCL_LINE
                                                , randomStartPointRequested             = randomStartPointRequested             & ! LCOV_EXCL_LINE
                                                , randomStartPointDomainLowerLimitVec   = randomStartPointDomainLowerLimitVec   & ! LCOV_EXCL_LINE
                                                , randomStartPointDomainUpperLimitVec   = randomStartPointDomainUpperLimitVec   & ! LCOV_EXCL_LINE
                                                )
            call self%SpecDRAM%setFromInputArgs ( adaptiveUpdateCount                   = adaptiveUpdateCount                   & ! LCOV_EXCL_LINE
                                                , adaptiveUpdatePeriod                  = adaptiveUpdatePeriod                  & ! LCOV_EXCL_LINE
                                                , greedyAdaptationCount                 = greedyAdaptationCount                 & ! LCOV_EXCL_LINE
                                                , delayedRejectionCount                 = delayedRejectionCount                 & ! LCOV_EXCL_LINE
                                                , burninAdaptationMeasure               = burninAdaptationMeasure               & ! LCOV_EXCL_LINE
                                                , delayedRejectionScaleFactorVec        = delayedRejectionScaleFactorVec        & ! LCOV_EXCL_LINE
                                                )
#elif defined PARANEST
            call self%SpecNest%setFromInputArgs ( tightness                             = tightness                             & ! LCOV_EXCL_LINE
                                                , tolerance                             = tolerance                             & ! LCOV_EXCL_LINE
                                                , scaleFactor                           = scaleFactor                           & ! LCOV_EXCL_LINE
                                                , proposalModel                         = proposalModel                         & ! LCOV_EXCL_LINE
                                                , liveSampleSize                        = liveSampleSize                        & ! LCOV_EXCL_LINE
                                                , inclusionFraction                     = inclusionFraction                     & ! LCOV_EXCL_LINE
                                                , adaptiveUpdateCount                   = adaptiveUpdateCount                   & ! LCOV_EXCL_LINE
                                                , adaptiveUpdatePeriod                  = adaptiveUpdatePeriod                  & ! LCOV_EXCL_LINE
                                                , mahalSqWeightExponent                 = mahalSqWeightExponent                 & ! LCOV_EXCL_LINE
                                                , stabilizationRequested                = stabilizationRequested                & ! LCOV_EXCL_LINE
                                                , maxAllowedKmeansFailure               = maxAllowedKmeansFailure               & ! LCOV_EXCL_LINE
                                                , maxAllowedMinVolFailure               = maxAllowedMinVolFailure               & ! LCOV_EXCL_LINE
                                                , maxKvolumeLoopRecursion               = maxKvolumeLoopRecursion               & ! LCOV_EXCL_LINE
                                                )
#endif
        end if

#if (defined MPI_ENABLED || defined CAF_ENABLED) && (defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED)
        block; use Err_mod, only: bcastErr; call bcastErr(self%Err); end block
#endif
        if (self%Err%occurred) then
            ! LCOV_EXCL_START
            self%Err%msg = PROCEDURE_NAME // self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = "\n", outputUnit = self%LogFile%unit )
            return
            ! LCOV_EXCL_STOP
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Now depending on the requested parallelism type, determine the process leader/rooter type and open output files
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        self%Image%isLeader = self%Image%isFirst .or. self%SpecBase%ParallelizationModel%isMultiChain
        self%Image%isRooter = .not. self%Image%isLeader

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! setup log file and open output files
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call self%setupOutputFiles()
#if (defined MPI_ENABLED || defined CAF_ENABLED) && (defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED)
        block; use Err_mod, only: bcastErr; call bcastErr(self%Err); end block
#endif
        if (self%Err%occurred) then
            ! LCOV_EXCL_START
            self%Err%msg = PROCEDURE_NAME // self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
            ! LCOV_EXCL_STOP
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! report variable values to the log file
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (self%isFreshRun) then
            call self%SpecBase%reportValues(prefix=self%brand,outputUnit=self%LogFile%unit,isLeaderImage=self%Image%isLeader)
#if defined PARADRAM || defined PARADISE
            call self%SpecMCMC%reportValues(prefix=self%brand,outputUnit=self%LogFile%unit,isLeaderImage=self%Image%isLeader,methodName=self%name,splashModeRequested=self%SpecBase%SilentModeRequested%isFalse)
            call self%SpecDRAM%reportValues(prefix=self%brand,outputUnit=self%LogFile%unit,isLeaderImage=self%Image%isLeader,splashModeRequested=self%SpecBase%SilentModeRequested%isFalse)
#elif defined PARANEST
            call self%SpecNest%reportValues(prefix=self%brand,outputUnit=self%LogFile%unit,isLeaderImage=self%Image%isLeader,splashModeRequested=self%SpecBase%SilentModeRequested%isFalse)
#endif
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! perform variable value sanity checks
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        self%Err%msg = ""
        self%Err%occurred = .false.

        call self%SpecBase%checkForSanity(Err = self%Err, methodName = self%name)
#if defined PARADRAM || defined PARADISE
        call self%SpecMCMC%checkForSanity(SpecBase = self%SpecBase, methodName = self%name, Err = self%Err, nd = ndim)
        call self%SpecDRAM%checkForSanity(Err = self%Err, methodName = self%name, nd = ndim)
#elif defined PARANEST
        call self%SpecNest%checkForSanity(Err = self%Err, methodName = self%name, nd = ndim)
#endif

        !if (self%Image%isLeader) then
#if (defined MPI_ENABLED || defined CAF_ENABLED) && (defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED)
        block; use Err_mod, only: bcastErr; call bcastErr(self%Err); end block
#endif
        if (self%Err%occurred) then
            self%Err%msg = PROCEDURE_NAME // self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = "\n", outputUnit = self%LogFile%unit )
            return
        end if
        !end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! setup output files
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! setup chain file

        block
            use ParaXXXX_ChainFileContents_mod, only: ChainFileContents_type
            self%Chain = ChainFileContents_type( ndim = ndim, variableNameList = self%SpecBase%VariableNameList%Val )
        end block

        if (self%isFreshRun) then
            if (self%Image%isLeader) then
                call self%Chain%writeHeader ( ndim              = ndim &
                                            , chainFileUnit     = self%ChainFile%unit &
                                            , isBinary          = self%SpecBase%ChainFileFormat%isBinary &
                                            , chainFileFormat   = self%ChainFile%format &
                                            )
                flush(self%ChainFile%unit)
            end if
        else
            call self%Chain%getLenHeader(ndim=ndim,isBinary=self%SpecBase%ChainFileFormat%isBinary,chainFileFormat=self%ChainFile%format)
        end if

        ! setup progress file

        if (self%Image%isLeader) then
            if (.not.self%isFreshRun) then
                read(self%TimeFile%unit, *, iostat = self%Err%stat) ! read the header line of the time file, only by leader images
            end if
            if (self%isFreshRun .or. is_iostat_end(self%Err%stat)) then
                    write( self%TimeFile%unit, self%TimeFile%format ) "NumFuncCallTotal" &
                                                                    , "NumFuncCallAccepted" &
                                                                    , "MeanAcceptanceRateSinceStart" &
                                                                    , "MeanAcceptanceRateSinceLastReport" &
                                                                    , "TimeElapsedSinceLastReportInSeconds" &
                                                                    , "TimeElapsedSinceStartInSeconds" &
                                                                    , "TimeRemainedToFinishInSeconds"
                    flush(self%TimeFile%unit)
            end if
        end if

        ! The format setup of setupOutputFiles() uses the generic g0 edit descriptor. Here the format is revised to be more specific.
        ! g0 edit descriptor format is slightly more arbitrary and compiler-dependent.

        associate(colWidth => self%SpecBase%OutputColumnWidth%str, precision => self%SpecBase%OutputRealPrecision%str, delim => self%SpecBase%OutputDelimiter%val)
            if (self%SpecBase%OutputColumnWidth%val>0) then
                self%TimeFile%format  = "("//"2I"//colWidth//",'"//delim//"',5(E"//colWidth//"."//precision//"E3"//",:,'"//delim//"'))"
#if defined PARADRAM || defined PARADISE
                self%ChainFile%format = "("// &
                                        "2(I"//colWidth//",'"//delim//"')"// &
                                        ","// &
                                        "2(E"//colWidth//"."//precision//"E3"//",'"//delim//"')"// &
                                        ","// &
                                        "2(I"//colWidth//",'"//delim//"')"// &
                                        ","// &
                                        num2str(1+self%nd%val)//"(E"//colWidth//"."//precision//"E3"//",:,'"//delim//"')"// &
                                        ")"
#elif defined PARANEST
                self%ChainFile%format = "("// &
                                        "1(I"//colWidth//",'"//delim//"')"// &
                                        ","// &
                                        num2str(5+self%nd%val)//"(E"//colWidth//"."//precision//"E3"//",:,'"//delim//"')"// &
                                        ")"
#endif
            end if
        end associate

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! setup the proposal distribution
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (allocated(self%Proposal)) deallocate(self%Proposal)
#if defined PARADRAM || defined PARADISE
        if (self%SpecMCMC%ProposalModel%isNormal) then
            allocate( self%Proposal, source = ProposalNormal_type   ( ndim          = ndim & ! LCOV_EXCL_LINE
                                                                    , SpecBase      = self%SpecBase & ! LCOV_EXCL_LINE
                                                                    , SpecMCMC      = self%SpecMCMC & ! LCOV_EXCL_LINE
                                                                    , SpecDRAM      = self%SpecDRAM & ! LCOV_EXCL_LINE
                                                                    , Image         = self%Image & ! LCOV_EXCL_LINE
                                                                    , name          = self%name & ! LCOV_EXCL_LINE
                                                                    , brand         = self%brand & ! LCOV_EXCL_LINE
                                                                    , LogFile       = self%LogFile & ! LCOV_EXCL_LINE
                                                                    , RestartFile   = self%RestartFile & ! LCOV_EXCL_LINE
                                                                    , isFreshRun    = self%isFreshRun & ! LCOV_EXCL_LINE
                                                                    ) )
        elseif (self%SpecMCMC%ProposalModel%isUniform) then
            allocate( self%Proposal, source = ProposalUniform_type  ( ndim          = ndim & ! LCOV_EXCL_LINE
                                                                    , SpecBase      = self%SpecBase & ! LCOV_EXCL_LINE
                                                                    , SpecMCMC      = self%SpecMCMC & ! LCOV_EXCL_LINE
                                                                    , SpecDRAM      = self%SpecDRAM & ! LCOV_EXCL_LINE
                                                                    , Image         = self%Image & ! LCOV_EXCL_LINE
                                                                    , name          = self%name & ! LCOV_EXCL_LINE
                                                                    , brand         = self%brand & ! LCOV_EXCL_LINE
                                                                    , LogFile       = self%LogFile & ! LCOV_EXCL_LINE
                                                                    , RestartFile   = self%RestartFile & ! LCOV_EXCL_LINE
                                                                    , isFreshRun    = self%isFreshRun & ! LCOV_EXCL_LINE
                                                                    ) )
#elif defined PARANEST
        if (self%SpecMCMC%ProposalModel%isRejEll) then
            allocate( self%Proposal, source = ProposalRejEll_type   ( ndim          = ndim & ! LCOV_EXCL_LINE
                                                                    , SpecBase      = self%SpecBase & ! LCOV_EXCL_LINE
                                                                    , SpecMCMC      = self%SpecMCMC & ! LCOV_EXCL_LINE
                                                                    , SpecDRAM      = self%SpecDRAM & ! LCOV_EXCL_LINE
                                                                    , Image         = self%Image & ! LCOV_EXCL_LINE
                                                                    , name          = self%name & ! LCOV_EXCL_LINE
                                                                    , brand         = self%brand & ! LCOV_EXCL_LINE
                                                                    , LogFile       = self%LogFile & ! LCOV_EXCL_LINE
                                                                    , RestartFile   = self%RestartFile & ! LCOV_EXCL_LINE
                                                                    , isFreshRun    = self%isFreshRun & ! LCOV_EXCL_LINE
                                                                    ) )
#endif
        else
            ! LCOV_EXCL_START
            self%Err%occurred = .true.
            self%Err%msg = PROCEDURE_NAME // ": Internal error occurred. Unsupported proposal distribution for " // self%name // "."
            call self%abort( Err = self%Err, prefix = self%brand, newline = "\n", outputUnit = self%LogFile%unit )
            error stop
            return
            ! LCOV_EXCL_STOP
        end if

#if (defined MPI_ENABLED || defined CAF_ENABLED) && (defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED)
        block; use Err_mod, only: bcastErr; call bcastErr(self%Proposal%Err); end block
#endif
        if (self%Proposal%Err%occurred) then
            self%Err%occurred = .true.
            self%Err%msg = PROCEDURE_NAME // self%Proposal%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = "\n", outputUnit = self%LogFile%unit )
            return
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! run ParaXXXX kernel
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (self%isFreshRun .and. self%Image%isLeader) then
            call self%Decor%writeDecoratedText  ( text = "\nStarting the " // self%name // " sampling - " // getNiceDateTime() // "\n" &
                                                , marginTop = 1_IK  &
                                                , marginBot = 1_IK  &
                                                , newline = "\n"    &
                                                , outputUnit = self%LogFile%unit )
        end if

        call self%runKernel( getLogFunc = getLogFunc )
#if (defined MPI_ENABLED || defined CAF_ENABLED) && (defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED)
        block; use Err_mod, only: bcastErr; call bcastErr(self%Err); end block
#endif
        ! relevant only for MATLAB / Python
        if (self%Err%occurred) then
            self%Err%msg = PROCEDURE_NAME // self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if

        if (self%isFreshRun .and. self%Image%isLeader) then
            call self%Decor%writeDecoratedText  ( text = "\nExiting the " // self%name // " sampling - " // getNiceDateTime() // "\n" &
                                                , marginTop = 1_IK  &
                                                , marginBot = 1_IK  &
                                                , newline = "\n"    &
                                                , outputUnit = self%LogFile%unit )
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! start ParaXXXX post-processing
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call self%postprocess()

#if (defined MPI_ENABLED || defined CAF_ENABLED) && (defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED)
        block; use Err_mod, only: bcastErr; call bcastErr(self%Err); end block
#endif
        if (self%Err%occurred) then
            ! LCOV_EXCL_START
            self%Err%msg = PROCEDURE_NAME // self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = "\n", outputUnit = self%LogFile%unit )
            return
            ! LCOV_EXCL_STOP
        end if

        !nullify(self%Proposal)
#if defined CAF_ENABLED || defined MPI_ENABLED
            if (self%SpecBase%MpiFinalizeRequested%val) call self%Image%finalize()
#endif

    end subroutine runSampler

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!end submodule Setup_smod

#undef ParaXXXX
#undef ParaXXXX_type
#undef ParaXXXX_ChainFileContents_mod
