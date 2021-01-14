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
!> This file implements the body of the `Input_smod` submodules of the `ParaDRAM_mod` and `ParaDISE_mod` modules.
!>
!> \remark
!> This module requires preprocessing, prior to compilation.
!>
!> \author Amir Shahmoradi

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

    character(*), parameter :: SUBMODULE_NAME = MODULE_NAME // "@Input_smod"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if defined PARADRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! This will be used in the read statement
#define ParaXXXX ParaDRAM
! This will be used in the declaration of the parent object
#define ParaXXXX_type ParaDRAM_type
! This will be used in the namelist declaration
#define NAMELIST ParaDRAM
#include "ParaXXXX_mod@Input_smod.nml.inc.f90"
#undef NAMELIST

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif defined PARADISE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! This will be used in the read statement
#define ParaXXXX ParaDISE
! This will be used in the declaration of the parent object
#define ParaXXXX_type ParaDISE_type
! This will be used in the namelist declaration
#define NAMELIST ParaDISE
#include "ParaXXXX_mod@Input_smod.nml.inc.f90"
#undef NAMELIST

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#else
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#error "Unrecognized sampler in ParaXXXX_mod@Input_mod.inc.f90"
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NAMELIST paraxxxx
#include "ParaXXXX_mod@Input_smod.nml.inc.f90"
#undef NAMELIST

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

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
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSpecFromInputFile
#endif
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str

        implicit none
        class(ParaXXXX_type), intent(inout) :: self
        integer(IK), intent(in)             :: nd
        character(*), parameter             :: PROCEDURE_NAME = SUBMODULE_NAME//"@getSpecFromInputFile()"

        ! initialize/nullify all general input options

        call self%SpecBase%nullifyNameListVar(nd)
        call self%SpecMCMC%nullifyNameListVar(nd)
        call self%SpecDRAM%nullifyNameListVar(nd)

        ! read input options if input file is provided

        blockReadInputFile: if (self%inputFileArgIsPresent) then

            blockInputFileType: if (self%InputFile%isInternal) then

                ! read input file as an internal file

                read(self%InputFile%Path%original,nml=ParaXXXX,iostat=self%InputFile%Err%stat)
                self%Err = self%InputFile%getReadErr(self%InputFile%Err%stat,self%InputFile%Path%modified)

                if (self%Err%occurred) then

                    if ( is_iostat_end(self%Err%stat) .or. is_iostat_eor(self%Err%stat) ) then

                        ! search for the paraxxxx namelist group in the file.

                        call self%warnUserAboutMissingNamelist(namelist = self%name)
                        read(self%InputFile%Path%original,nml=paraxxxx,iostat=self%InputFile%Err%stat) ! WARNING: "paraxxxx" is NOT the same as fpp macro name "ParaXXXX". This is a real namelist name.
                        self%Err = self%InputFile%getReadErr(self%InputFile%Err%stat,self%InputFile%Path%modified)

                    end if

                    if (self%Err%occurred) then

                        if (is_iostat_end(self%Err%stat) .or. is_iostat_eor(self%Err%stat)) then

                            call self%warnUserAboutMissingNamelist(namelist = "ParaXXXX")

                        else

                            ! LCOV_EXCL_START
                            read(self%InputFile%Path%original,nml=ParaXXXX) ! let the compiler print diagnostic messages, should any error happen.
                            ! LCOV_EXCL_STOP

                        end if

                    end if

                end if

            else blockInputFileType ! the input file is external

                ! close input file if it is open

                if (self%InputFile%isOpen) then
                    ! LCOV_EXCL_START
                    close(unit=self%InputFile%unit,iostat=self%InputFile%Err%stat)
                    self%Err = self%InputFile%getCloseErr(self%InputFile%Err%stat)
                    if (self%Err%occurred) then
                        self%Err%msg =  PROCEDURE_NAME // ": Error occurred while attempting to close the user-provided input file='" // &
                                        self%InputFile%Path%modified // "', unit=" // num2str(self%InputFile%unit) // ".\n" // &
                                        self%Err%msg
                        return
                    end if
                    ! LCOV_EXCL_STOP
                end if

                ! open input file

                open( newunit = self%InputFile%unit       &
                    , file = self%InputFile%Path%modified &
                    , status = self%InputFile%status      &
                    , iostat = self%InputFile%Err%stat    &
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
                    , SHARED &
#endif
                    )
                self%Err = self%InputFile%getOpenErr(self%InputFile%Err%stat)
                if (self%Err%occurred) then
                    ! LCOV_EXCL_START
                    self%Err%msg =  PROCEDURE_NAME // ": Error occurred while attempting to open the user-provided input file='" // &
                                    self%InputFile%Path%modified // "', unit=" // num2str(self%InputFile%unit) // ".\n" // &
                                    self%Err%msg
                    return
                    ! LCOV_EXCL_STOP
                end if

                ! read input file

                read(self%InputFile%unit,nml=ParaXXXX,iostat=self%InputFile%Err%stat)
                self%Err = self%InputFile%getReadErr(self%InputFile%Err%stat,self%InputFile%Path%modified)

                if (self%Err%occurred) then

                    if ( is_iostat_end(self%Err%stat) .or. is_iostat_eor(self%Err%stat) ) then

                        ! search for the paraxxxx namelist group in the file.

                        rewind(self%InputFile%unit)
                        call self%warnUserAboutMissingNamelist(namelist = self%name)
                        read(self%InputFile%unit, nml=paraxxxx, iostat=self%InputFile%Err%stat) ! WARNING: "paraxxxx" is NOT the same as fpp macro name "ParaXXXX"
                        self%Err = self%InputFile%getReadErr(self%InputFile%Err%stat,self%InputFile%Path%modified)

                    end if

                    if (self%Err%occurred) then

                        if (is_iostat_end(self%Err%stat) .or. is_iostat_eor(self%Err%stat)) then

                            call self%warnUserAboutMissingNamelist(namelist = "ParaXXXX")

                        else ! attempt to read the file one more time, without error handling, so that the compiler prints out the error message.

                            ! LCOV_EXCL_START
                            rewind(self%InputFile%unit)
                            read(self%InputFile%unit, nml=ParaXXXX)
                            ! LCOV_EXCL_STOP

                        end if

                    end if

                end if

                ! close input file

                close(unit=self%InputFile%unit,iostat=self%InputFile%Err%stat)
                ! LCOV_EXCL_START
                self%Err = self%InputFile%getCloseErr(self%InputFile%Err%stat)
                if (self%Err%occurred) then
                    self%Err%msg =  PROCEDURE_NAME // ": Error occurred while attempting to close the user-provided input file='" // &
                                    self%InputFile%Path%modified // "', unit=" // num2str(self%InputFile%unit) // ".\n" // &
                                    self%Err%msg
                    return
                end if
                ! LCOV_EXCL_STOP

            end if blockInputFileType

        end if blockReadInputFile

        ! setup SpecBase variables that have been read form the input file

        call self%SpecBase%setFromInputFile(Err = self%Err)
        ! LCOV_EXCL_START
        if (self%Err%occurred) then
            self%Err%msg = PROCEDURE_NAME // self%Err%msg
            return
        end if
        ! LCOV_EXCL_STOP

        ! setup SpecMCMC variables that have been read form the input file

        call self%SpecMCMC%setFromInputFile(Err = self%Err)
        ! LCOV_EXCL_START
        if (self%Err%occurred) then
            self%Err%msg = PROCEDURE_NAME // self%Err%msg
            return
        end if
        ! LCOV_EXCL_STOP

        ! setup SpecDRAM variables that have been read form the input file

        call self%SpecDRAM%setFromInputFile( self%Err )
        ! LCOV_EXCL_START
        if (self%Err%occurred) then
            self%Err%msg = PROCEDURE_NAME // self%Err%msg
            return
        end if
        ! LCOV_EXCL_STOP

    end subroutine getSpecFromInputFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ParaXXXX_type
#undef ParaXXXX

