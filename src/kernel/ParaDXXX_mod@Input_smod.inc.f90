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

!> \brief
!> This file implements the body of the `Input_smod` submodules of the `ParaDRAM_mod` and `ParaDISE_mod` modules.
!>
!> \remark
!> This module requires preprocessing, prior to compilation.
!>
!> @author Amir Shahmoradi

#if defined PARADRAM

#define ParaDXXX ParaDRAM
#define ParaDXXX_type ParaDRAM_type

#elif defined PARADISE

#define ParaDXXX ParaDISE
#define ParaDXXX_type ParaDISE_type

#else
#error "Unrecognized sampler in ParaDXXX_mod@Input_mod.inc.f90"
#endif

    implicit none

    character(*), parameter :: SUBMODULE_NAME = MODULE_NAME // "@Input_smod"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a method of [ParaDRAM_type](@ref paradram_type) and [ParaDISE_type](@ref paradise_type) classes.
    !> Read the input file and assign the simulation specification variables.
    !>
    !> @param[inout]    self    :   An object of class [ParaDRAM_type](@ref paradram_type) or [ParaDISE_type](@ref paradise_type).
    !> @param[in]       nd      :   The number of dimensions of the domain of the objective function.
    !>
    !> \remark
    !> This procedure requires preprocessing.
    module subroutine getSpecFromInputFile(self,nd)
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
        use SpecBase_OverwriteRequested_mod                 , only: overwriteRequested
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
        use SpecBase_SystemInfoFilePath_mod                 , only: systemInfoFilePath
        use SpecBase_InterfaceType_mod                      , only: interfaceType

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
        class(ParaDXXX_type), intent(inout) :: self
        integer(IK), intent(in)             :: nd
        character(*), parameter             :: PROCEDURE_NAME = SUBMODULE_NAME//"@getSpecFromInputFile()"

        ! ParaMonte variables
        namelist /ParaDXXX/ sampleSize
        namelist /ParaDXXX/ randomSeed
        namelist /ParaDXXX/ description
        namelist /ParaDXXX/ outputFileName
        namelist /ParaDXXX/ outputDelimiter
        namelist /ParaDXXX/ ChainFileFormat
        namelist /ParaDXXX/ variableNameList
        namelist /ParaDXXX/ restartFileFormat
        namelist /ParaDXXX/ outputColumnWidth
        namelist /ParaDXXX/ overwriteRequested
        namelist /ParaDXXX/ outputRealPrecision
        namelist /ParaDXXX/ silentModeRequested
        namelist /ParaDXXX/ domainLowerLimitVec
        namelist /ParaDXXX/ domainUpperLimitVec
        namelist /ParaDXXX/ parallelizationModel
        namelist /ParaDXXX/ progressReportPeriod
        namelist /ParaDXXX/ inputFileHasPriority
        namelist /ParaDXXX/ targetAcceptanceRate
        namelist /ParaDXXX/ mpiFinalizeRequested
        namelist /ParaDXXX/ maxNumDomainCheckToWarn
        namelist /ParaDXXX/ maxNumDomainCheckToStop
        namelist /ParaDXXX/ systemInfoFilePath
        namelist /ParaDXXX/ interfaceType

        ! ParaMCMC variables
        namelist /ParaDXXX/ chainSize
        namelist /ParaDXXX/ scaleFactor
        namelist /ParaDXXX/ startPointVec
        namelist /ParaDXXX/ proposalModel
        namelist /ParaDXXX/ proposalStartStdVec
        namelist /ParaDXXX/ proposalStartCorMat
        namelist /ParaDXXX/ proposalStartCovMat
        namelist /ParaDXXX/ sampleRefinementCount
        namelist /ParaDXXX/ sampleRefinementMethod
        namelist /ParaDXXX/ randomStartPointRequested
        namelist /ParaDXXX/ randomStartPointDomainLowerLimitVec
        namelist /ParaDXXX/ randomStartPointDomainUpperLimitVec

        ! ParaDRAM variables
        namelist /ParaDXXX/ adaptiveUpdateCount
        namelist /ParaDXXX/ adaptiveUpdatePeriod
        namelist /ParaDXXX/ greedyAdaptationCount
        namelist /ParaDXXX/ delayedRejectionCount
        namelist /ParaDXXX/ burninAdaptationMeasure
        namelist /ParaDXXX/ delayedRejectionScaleFactorVec

        ! initialize/nullify all general input options

        call self%SpecBase%nullifyNameListVar(nd)
        call self%SpecMCMC%nullifyNameListVar(nd)
        call self%SpecDRAM%nullifyNameListVar(nd)

        ! read input options if input file is provided

        blockReadInputFile: if (self%inputFileArgIsPresent) then

            blockInputFileType: if (self%InputFile%isInternal) then

                ! read input file as an internal file

!#if defined DBG_ENABLED
                read(self%InputFile%Path%original,nml=ParaDXXX)
!#else
!                read(self%InputFile%Path%original,nml=ParaDXXX,iostat=self%InputFile%Err%stat)
!                self%Err = self%InputFile%getReadErr(self%InputFile%Err%stat,self%InputFile%Path%modified)
!                if (self%Err%occurred) then
!                    if (is_iostat_end(self%Err%stat)) then
!                        call self%warnUserAboutMissingNamelist(self%brand,self%name,self%name,self%LogFile%unit)
!                    else
!                        self%Err%msg = PROCEDURE_NAME // self%Err%msg
!                        return
!                    end if
!                end if
!#endif

            else blockInputFileType ! the input file is external

                ! close input file if it is open

                if (self%InputFile%isOpen) then
                    close(unit=self%InputFile%unit,iostat=self%InputFile%Err%stat)
                    self%Err = self%InputFile%getCloseErr(self%InputFile%Err%stat)
                    if (self%Err%occurred) then
                        self%Err%msg =  PROCEDURE_NAME // ": Error occurred while attempting to close the user-provided input file='" // &
                                        self%InputFile%Path%modified // "', unit=" // num2str(self%InputFile%unit) // ".\n" // &
                                        self%Err%msg
                        return
                    end if
                end if

                ! open input file

                open( newunit = self%InputFile%unit       &
                    , file = self%InputFile%Path%modified &
                    , status = self%InputFile%status      &
                    , iostat = self%InputFile%Err%stat    &
#if defined IFORT_ENABLED && defined OS_IS_WINDOWS
                    , SHARED &
#endif
                    )
                self%Err = self%InputFile%getOpenErr(self%InputFile%Err%stat)
                if (self%Err%occurred) then
                    self%Err%msg =  PROCEDURE_NAME // ": Error occurred while attempting to open the user-provided input file='" // &
                                    self%InputFile%Path%modified // "', unit=" // num2str(self%InputFile%unit) // ".\n" // &
                                    self%Err%msg
                    return
                end if

                ! read input file

!#if defined DBG_ENABLED
                read(self%InputFile%unit,nml=ParaDXXX)
!#else
!                read(self%InputFile%unit,nml=ParaDXXX,iostat=self%InputFile%Err%stat)
!                self%Err = self%InputFile%getReadErr(self%InputFile%Err%stat,self%InputFile%Path%modified)
!                if (self%Err%occurred) then
!                    if (is_iostat_end(self%Err%stat)) then
!                        call self%warnUserAboutMissingNamelist(self%brand,self%name,self%name,self%LogFile%unit)
!                    else
!                        self%Err%msg = PROCEDURE_NAME // self%Err%msg
!                        return
!                    end if
!                end if
!#endif

                ! close input file

                close(unit=self%InputFile%unit,iostat=self%InputFile%Err%stat)
                self%Err = self%InputFile%getCloseErr(self%InputFile%Err%stat)
                if (self%Err%occurred) then
                    self%Err%msg =  PROCEDURE_NAME // ": Error occurred while attempting to close the user-provided input file='" // &
                                    self%InputFile%Path%modified // "', unit=" // num2str(self%InputFile%unit) // ".\n" // &
                                    self%Err%msg
                    return
                end if

            end if blockInputFileType

        end if blockReadInputFile

        ! setup SpecBase variables that have been read form the input file

        call self%SpecBase%setFromInputFile( Err = self%Err )
        if (self%Err%occurred) then
            self%Err%msg = PROCEDURE_NAME // self%Err%msg
            return
        end if

        ! setup SpecMCMC variables that have been read form the input file

        call self%SpecMCMC%setFromInputFile ( Err = self%Err &
                                            , nd = nd &
                                            , domainLowerLimitVec = self%SpecBase%DomainLowerLimitVec%Val &
                                            , domainUpperLimitVec = self%SpecBase%DomainUpperLimitVec%Val )
        if (self%Err%occurred) then
            self%Err%msg = PROCEDURE_NAME // self%Err%msg
            return
        end if

        ! setup SpecDRAM variables that have been read form the input file

        call self%SpecDRAM%setFromInputFile( self%Err )
        if (self%Err%occurred) then
            self%Err%msg = PROCEDURE_NAME // self%Err%msg
            return
        end if

    end subroutine getSpecFromInputFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ParaDXXX_type
#undef ParaDXXX

