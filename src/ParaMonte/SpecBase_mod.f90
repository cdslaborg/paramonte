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

module SpecBase_mod

    ! ParaMonte Spec variable types
    use SpecBase_SampleSize_mod                     , only: SampleSize_type
    use SpecBase_RandomSeed_mod                     , only: RandomSeed_type
    use SpecBase_Description_mod                    , only: Description_type
    use SpecBase_OutputFileName_mod                 , only: OutputFileName_type
    use SpecBase_OutputDelimiter_mod                , only: OutputDelimiter_type
    use SpecBase_ChainFileFormat_mod                , only: ChainFileFormat_type
    use SpecBase_VariableNameList_mod               , only: VariableNameList_type
    use SpecBase_RestartFileFormat_mod              , only: RestartFileFormat_type
    use SpecBase_OutputColumnWidth_mod              , only: OutputColumnWidth_type
    use SpecBase_OverwriteRequested_mod             , only: OverwriteRequested_type
    use SpecBase_OutputRealPrecision_mod            , only: OutputRealPrecision_type
    use SpecBase_SilentModeRequested_mod            , only: SilentModeRequested_type
    use SpecBase_DomainLowerLimitVec_mod            , only: DomainLowerLimitVec_type
    use SpecBase_DomainUpperLimitVec_mod            , only: DomainUpperLimitVec_type
    use SpecBase_ParallelizationModel_mod           , only: ParallelizationModel_type
    use SpecBase_InputFileHasPriority_mod           , only: InputFileHasPriority_type
    use SpecBase_ProgressReportPeriod_mod           , only: ProgressReportPeriod_type
    use SpecBase_TargetAcceptanceRate_mod           , only: TargetAcceptanceRate_type
    use SpecBase_MpiFinalizeRequested_mod           , only: MpiFinalizeRequested_type
    use SpecBase_MaxNumDomainCheckToWarn_mod        , only: MaxNumDomainCheckToWarn_type
    use SpecBase_MaxNumDomainCheckToStop_mod        , only: MaxNumDomainCheckToStop_type
    use SpecBase_InterfaceType_mod                  , only: InterfaceType_type
    use SpecBase_SystemInfoFilePath_mod             , only: SystemInfoFilePath_type

    ! ParaMonte namelist variables

    use SpecBase_SampleSize_mod                     , only: sampleSize
    use SpecBase_RandomSeed_mod                     , only: randomSeed
    use SpecBase_Description_mod                    , only: description
    use SpecBase_OutputFileName_mod                 , only: outputFileName
    use SpecBase_OutputDelimiter_mod                , only: outputDelimiter
    use SpecBase_ChainFileFormat_mod                , only: ChainFileFormat
    use SpecBase_VariableNameList_mod               , only: variableNameList
    use SpecBase_RestartFileFormat_mod              , only: restartFileFormat
    use SpecBase_OutputColumnWidth_mod              , only: outputColumnWidth
    use SpecBase_OverwriteRequested_mod             , only: overwriteRequested
    use SpecBase_OutputRealPrecision_mod            , only: outputRealPrecision
    use SpecBase_SilentModeRequested_mod            , only: silentModeRequested
    use SpecBase_DomainLowerLimitVec_mod            , only: domainLowerLimitVec
    use SpecBase_DomainUpperLimitVec_mod            , only: domainUpperLimitVec
    use SpecBase_ParallelizationModel_mod           , only: ParallelizationModel
    use SpecBase_InputFileHasPriority_mod           , only: inputFileHasPriority
    use SpecBase_ProgressReportPeriod_mod           , only: progressReportPeriod
    use SpecBase_TargetAcceptanceRate_mod           , only: targetAcceptanceRate
    use SpecBase_MpiFinalizeRequested_mod           , only: mpiFinalizeRequested
    use SpecBase_MaxNumDomainCheckToWarn_mod        , only: maxNumDomainCheckToWarn
    use SpecBase_MaxNumDomainCheckToStop_mod        , only: maxNumDomainCheckToStop
    use SpecBase_InterfaceType_mod                  , only: interfaceType
    use SpecBase_SystemInfoFilePath_mod             , only: systemInfoFilePath

    implicit none

    character(*), parameter :: MODULE_NAME = "@SpecBase_mod"

    type                                        :: SpecBase_type
        type(SampleSize_type)                   :: SampleSize
        type(RandomSeed_type)                   :: RandomSeed
        type(Description_type)                  :: Description
        type(OutputFileName_type)               :: OutputFileName
        type(OutputDelimiter_type)              :: OutputDelimiter
        type(ChainFileFormat_type)              :: ChainFileFormat
        type(VariableNameList_type)             :: VariableNameList
        type(RestartFileFormat_type)            :: RestartFileFormat
        type(OutputColumnWidth_type)            :: OutputColumnWidth
        type(OverwriteRequested_type)           :: OverwriteRequested
        type(OutputRealPrecision_type)          :: OutputRealPrecision
        type(SilentModeRequested_type)          :: SilentModeRequested
        type(DomainLowerLimitVec_type)          :: domainLowerLimitVec
        type(DomainUpperLimitVec_type)          :: domainUpperLimitVec
        type(ParallelizationModel_type)         :: ParallelizationModel
        type(InputFileHasPriority_type)         :: InputFileHasPriority
        type(ProgressReportPeriod_type)         :: ProgressReportPeriod
        type(TargetAcceptanceRate_type)         :: TargetAcceptanceRate
        type(MpiFinalizeRequested_type)         :: MpiFinalizeRequested
        type(MaxNumDomainCheckToWarn_type)      :: MaxNumDomainCheckToWarn
        type(MaxNumDomainCheckToStop_type)      :: MaxNumDomainCheckToStop
        type(InterfaceType_type)                :: InterfaceType
        type(SystemInfoFilePath_type)           :: SystemInfoFilePath
    contains
        procedure, pass                         :: nullifyNameListVar
        procedure, pass                         :: setFromInputFile
        procedure, pass                         :: setFromInputArgs
        procedure, pass                         :: checkForSanity
        procedure, pass                         :: reportValues
    end type SpecBase_type

    interface SpecBase_type
        module procedure                        :: constructSpecBase
    end interface SpecBase_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructSpecBase(nd,methodName,imageID,imageCount) result(SpecBase)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructSpecBase
#endif
        use Constants_mod, only: IK
        implicit none
        character(*), intent(in)    :: methodName
        integer(IK), intent(in)     :: nd
        integer(IK), intent(in)     :: imageID, imageCount
        type(SpecBase_type)         :: SpecBase
        SpecBase%SampleSize                         = SampleSize_type(methodName)
        SpecBase%RandomSeed                         = RandomSeed_type(methodName,imageID,imageCount)
        SpecBase%Description                        = Description_type(methodName)
        SpecBase%OutputFileName                     = OutputFileName_type(methodName)
        SpecBase%OutputDelimiter                    = OutputDelimiter_type(methodName)
        SpecBase%ChainFileFormat                    = ChainFileFormat_type(methodName)
        SpecBase%VariableNameList                   = VariableNameList_type(nd,methodName)
        SpecBase%RestartFileFormat                  = RestartFileFormat_type(methodName)
        SpecBase%OutputColumnWidth                  = OutputColumnWidth_type(methodName)
        SpecBase%OverwriteRequested                 = OverwriteRequested_type(methodName)
        SpecBase%OutputRealPrecision                = OutputRealPrecision_type(methodName)
        SpecBase%SilentModeRequested                = SilentModeRequested_type(methodName)
        SpecBase%DomainLowerLimitVec                = DomainLowerLimitVec_type(methodName)
        SpecBase%DomainUpperLimitVec                = DomainUpperLimitVec_type(methodName)
        SpecBase%ParallelizationModel               = ParallelizationModel_type(methodName)
        SpecBase%InputFileHasPriority               = InputFileHasPriority_type(methodName)
        SpecBase%ProgressReportPeriod               = ProgressReportPeriod_type()
        SpecBase%TargetAcceptanceRate               = TargetAcceptanceRate_type(methodName)
        SpecBase%MpiFinalizeRequested               = MpiFinalizeRequested_type(methodName)
        SpecBase%MaxNumDomainCheckToWarn            = MaxNumDomainCheckToWarn_type()
        SpecBase%MaxNumDomainCheckToStop            = MaxNumDomainCheckToStop_type()
        SpecBase%InterfaceType                      = InterfaceType_type()
        SpecBase%SystemInfoFilePath                 = SystemInfoFilePath_type()
    end function constructSpecBase

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar( SpecBase, nd )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        use Constants_mod, only: IK
        implicit none
        class(SpecBase_type), intent(inout) :: SpecBase
        integer(IK), intent(in)             :: nd
        ! nullify SpecBase global variables that have to be read form the input namelist file
        call SpecBase%SampleSize                    %nullifyNameListVar()
        call SpecBase%RandomSeed                    %nullifyNameListVar()
        call SpecBase%Description                   %nullifyNameListVar()
        call SpecBase%OutputFileName                %nullifyNameListVar()
        call SpecBase%OutputDelimiter               %nullifyNameListVar()
        call SpecBase%ChainFileFormat               %nullifyNameListVar()
        call SpecBase%VariableNameList              %nullifyNameListVar(nd)
        call SpecBase%RestartFileFormat             %nullifyNameListVar()
        call SpecBase%OutputColumnWidth             %nullifyNameListVar()
        call SpecBase%OverwriteRequested            %nullifyNameListVar()
        call SpecBase%DomainLowerLimitVec           %nullifyNameListVar(nd)
        call SpecBase%DomainUpperLimitVec           %nullifyNameListVar(nd)
        call SpecBase%OutputRealPrecision           %nullifyNameListVar()
        call SpecBase%SilentModeRequested           %nullifyNameListVar()
        call SpecBase%ProgressReportPeriod          %nullifyNameListVar()
        call SpecBase%ParallelizationModel          %nullifyNameListVar()
        call SpecBase%InputFileHasPriority          %nullifyNameListVar()
        call SpecBase%TargetAcceptanceRate          %nullifyNameListVar()
        call SpecBase%MpiFinalizeRequested          %nullifyNameListVar()
        call SpecBase%MaxNumDomainCheckToWarn       %nullifyNameListVar()
        call SpecBase%MaxNumDomainCheckToStop       %nullifyNameListVar()
        call SpecBase%InterfaceType                 %nullifyNameListVar()
        call SpecBase%SystemInfoFilePath            %nullifyNameListVar()
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setFromInputFile( SpecBase, Err )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setFromInputFile
#endif

        use Constants_mod, only: IK, RK
        use Err_mod, only: Err_type

        implicit none

        class(SpecBase_type), intent(inout) :: SpecBase
        type(Err_type), intent(inout)       :: Err

        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME // "@setFromInputFile()"

        call SpecBase%SampleSize                    %set(sampleSize)
        call SpecBase%RandomSeed                    %set(randomSeed,Err)
        call SpecBase%Description                   %set(description)
        call SpecBase%OutputFileName                %set(outputFileName)
        call SpecBase%ChainFileFormat               %set(chainFileFormat)
        call SpecBase%VariableNameList              %set(VariableNameList)
        call SpecBase%RestartFileFormat             %set(restartFileFormat)
        call SpecBase%DomainLowerLimitVec           %set(domainLowerLimitVec)
        call SpecBase%DomainUpperLimitVec           %set(domainUpperLimitVec)

        ! do not change the order with outputDelimiter
        call SpecBase%OutputColumnWidth             %set(outputColumnWidth)
        call SpecBase%OutputDelimiter               %set(outputDelimiter,SpecBase%OutputColumnWidth%val)

        call SpecBase%OverwriteRequested            %set(overwriteRequested)
        call SpecBase%OutputRealPrecision           %set(outputRealPrecision)
        call SpecBase%SilentModeRequested           %set(silentModeRequested)
        call SpecBase%ProgressReportPeriod          %set(progressReportPeriod)
        call SpecBase%ParallelizationModel          %set(parallelizationModel)
        call SpecBase%InputFileHasPriority          %set(inputFileHasPriority)
        call SpecBase%TargetAcceptanceRate          %set(TargetAcceptanceRate)
        call SpecBase%MpiFinalizeRequested          %set(mpiFinalizeRequested)
        call SpecBase%MaxNumDomainCheckToWarn       %set(maxNumDomainCheckToWarn)
        call SpecBase%MaxNumDomainCheckToStop       %set(maxNumDomainCheckToStop)
        call SpecBase%InterfaceType                 %set(interfaceType)
        call SpecBase%SystemInfoFilePath            %set(systemInfoFilePath)

        if (Err%occurred) Err%msg = PROCEDURE_NAME // Err%msg

    end subroutine setFromInputFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setFromInputArgs ( SpecBase, Err &
                                , sampleSize &
                                , randomSeed &
                                , description &
                                , outputFileName &
                                , outputDelimiter &
                                , chainFileFormat &
                                , variableNameList &
                                , domainLowerLimitVec &
                                , domainUpperLimitVec &
                                , restartFileFormat &
                                , outputColumnWidth &
                                , overwriteRequested &
                                , outputRealPrecision &
                                , silentModeRequested &
                                , parallelizationModel &
                                , progressReportPeriod &
                                , TargetAcceptanceRate &
                                , mpiFinalizeRequested &
                                , maxNumDomainCheckToWarn &
                                , maxNumDomainCheckToStop &
                                )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setFromInputArgs
#endif

        use Constants_mod, only: IK, RK
        use Err_mod, only: Err_type

        implicit none

        class(SpecBase_type), intent(inout) :: SpecBase
        type(Err_type), intent(inout)       :: Err

        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME // "@setFromInputArgs()"

        ! ParaMonte variables
        integer(IK) , intent(in), optional  :: sampleSize
        integer(IK) , intent(in), optional  :: randomSeed
        character(*), intent(in), optional  :: description
        character(*), intent(in), optional  :: outputFileName
        character(*), intent(in), optional  :: outputDelimiter
        character(*), intent(in), optional  :: chainFileFormat
        character(*), intent(in), optional  :: variableNameList(:)
        character(*), intent(in), optional  :: restartFileFormat
        integer(IK) , intent(in), optional  :: outputColumnWidth
        logical     , intent(in), optional  :: overwriteRequested
        integer(IK) , intent(in), optional  :: outputRealPrecision
        logical     , intent(in), optional  :: silentModeRequested
        real(RK)    , intent(in), optional  :: domainLowerLimitVec(:)
        real(RK)    , intent(in), optional  :: domainUpperLimitVec(:)
        character(*), intent(in), optional  :: parallelizationModel
        integer(IK) , intent(in), optional  :: progressReportPeriod
        real(RK)    , intent(in), optional  :: TargetAcceptanceRate(2)
        logical     , intent(in), optional  :: mpiFinalizeRequested
        integer(IK) , intent(in), optional  :: maxNumDomainCheckToWarn
        integer(IK) , intent(in), optional  :: maxNumDomainCheckToStop

        if (present(sampleSize))                    call SpecBase%SampleSize                    %set(sampleSize)
        if (present(randomSeed))                    call SpecBase%RandomSeed                    %set(randomSeed,Err)
        if (present(description))                   call SpecBase%Description                   %set(description)
        if (present(outputFileName))                call SpecBase%OutputFileName                %set(outputFileName)
        if (present(chainFileFormat))               call SpecBase%ChainFileFormat               %set(chainFileFormat)
        if (present(variableNameList))              call SpecBase%VariableNameList              %set(variableNameList)
        if (present(restartFileFormat))             call SpecBase%RestartFileFormat             %set(restartFileFormat)
        if (present(domainLowerLimitVec))           call SpecBase%DomainLowerLimitVec           %set(domainLowerLimitVec)
        if (present(domainUpperLimitVec))           call SpecBase%DomainUpperLimitVec           %set(domainUpperLimitVec)

        ! do not change the order with outputDelimiter
        if (present(outputColumnWidth))             call SpecBase%OutputColumnWidth             %set(outputColumnWidth)
        if (present(outputDelimiter))               call SpecBase%OutputDelimiter               %set(outputDelimiter,SpecBase%OutputColumnWidth%val)

        if (present(overwriteRequested))            call SpecBase%OverwriteRequested            %set(overwriteRequested)
        if (present(parallelizationModel))          call SpecBase%ParallelizationModel          %set(parallelizationModel)
        if (present(outputRealPrecision))           call SpecBase%OutputRealPrecision           %set(outputRealPrecision)
        if (present(silentModeRequested))           call SpecBase%SilentModeRequested           %set(silentModeRequested)
        if (present(progressReportPeriod))          call SpecBase%ProgressReportPeriod          %set(progressReportPeriod)
        if (present(TargetAcceptanceRate))          call SpecBase%TargetAcceptanceRate          %set(TargetAcceptanceRate)
        if (present(mpiFinalizeRequested))          call SpecBase%MpiFinalizeRequested          %set(mpiFinalizeRequested)
        if (present(maxNumDomainCheckToWarn))       call SpecBase%MaxNumDomainCheckToWarn       %set(maxNumDomainCheckToWarn)
        if (present(maxNumDomainCheckToStop))       call SpecBase%MaxNumDomainCheckToStop       %set(maxNumDomainCheckToStop)

        if (Err%occurred) Err%msg = PROCEDURE_NAME // Err%msg

    end subroutine setFromInputArgs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine reportValues(SpecBase,prefix,outputUnit,isMasterImage)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: reportValues
#endif
        use Decoration_mod, only: GENERIC_OUTPUT_FORMAT
        use Decoration_mod, only: GENERIC_TABBED_FORMAT
        use Decoration_mod, only: TAB
        use Constants_mod, only: IK, UNDEFINED
        use Err_mod, only: note, informUser
        implicit none
        class(SpecBase_type), intent(in)    :: SpecBase
        character(*), intent(in)            :: prefix
        integer(IK) , intent(in)            :: outputUnit
        logical     , intent(in)            :: isMasterImage
        integer(IK)                         :: i, ndim
        character(:), allocatable           :: msg

        if (isMasterImage) then

            ndim = size(SpecBase%DomainLowerLimitVec%Val(:))
            msg =   "ndim is a 32-bit positive integer, representing the number of dimensions of the domain of the objective function. &
                    &It is the only simulation specification variable that the user must always provide along with the objective function, &
                    &separately from the rest of the simulation specifications. The variable ndim must be always provided directly to the ParaMonte &
                    &routines, along with the objective function. If specified within an input file, its value will be ignored and not used. &
                    &The variable ndim has no default value as it is the only mandatory piece of information that must be provided by the user."
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "ndim"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) ndim
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = msg )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "description"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            call informUser( outputUnit = outputUnit, newline = "\n", wrapWidth = 125_IK, prefix = TAB//TAB, marginTop = 0_IK, marginBot = 0_IK, msg = SpecBase%Description%val )
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%Description%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "inputFileHasPriority"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%InputFileHasPriority%val
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%InputFileHasPriority%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "silentModeRequested"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%SilentModeRequested%val
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%SilentModeRequested%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "domainLowerLimitVec"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            do i = 1, ndim
                write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%DomainLowerLimitVec%Val(i)
            end do
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%DomainLowerLimitVec%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "domainUpperLimitVec"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            do i = 1, ndim
                write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%DomainUpperLimitVec%Val(i)
            end do
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%DomainUpperLimitVec%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "variableNameList"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            do i = 1, ndim
                write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%VariableNameList%Val(i)
            end do
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%VariableNameList%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "parallelizationModel"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%ParallelizationModel%val
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%ParallelizationModel%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "mpiFinalizeRequested"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%MpiFinalizeRequested%val
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%MpiFinalizeRequested%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "outputFileName"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%OutputFileName%modified
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%OutputFileName%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "overwriteRequested"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%OverwriteRequested%val
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%OverwriteRequested%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "targetAcceptanceRate"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            if ( SpecBase%TargetAcceptanceRate%scalingRequested ) then
                if ( SpecBase%TargetAcceptanceRate%Val(1)==SpecBase%TargetAcceptanceRate%Val(2) ) then
                    write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%TargetAcceptanceRate%Val(1)
                else
                    write(outputUnit,GENERIC_TABBED_FORMAT) "[", SpecBase%TargetAcceptanceRate%Val(1), ",", SpecBase%TargetAcceptanceRate%Val(2), "]"
                end if
            else
                write(outputUnit,GENERIC_TABBED_FORMAT) UNDEFINED
            end if
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%TargetAcceptanceRate%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "sampleSize"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%SampleSize%val
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%SampleSize%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "randomSeed"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            if ( SpecBase%RandomSeed%userSeed == SpecBase%RandomSeed%nullSeed ) then
                write(outputUnit,GENERIC_TABBED_FORMAT) UNDEFINED
            else
                write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%RandomSeed%userSeed
            end if


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) "ThisProcessID", "RandomSeedVectorSize", "RandomSeedVectorValues"
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%RandomSeed%imageID, SpecBase%RandomSeed%sizeSeed, SpecBase%RandomSeed%Seed(:,SpecBase%RandomSeed%imageID)
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) "OtherProcessID", "RandomSeedVectorSize", "RandomSeedVectorValues"
            if (SpecBase%RandomSeed%imageCount==1) then
                write(outputUnit,GENERIC_TABBED_FORMAT) "No other processor exists."
#if defined CAF_ENABLED || defined MPI_ENABLED
            else
                do i = 1, SpecBase%RandomSeed%imageCount
                    if (i/=SpecBase%RandomSeed%imageID) write(outputUnit,GENERIC_TABBED_FORMAT) i, SpecBase%RandomSeed%sizeSeed, SpecBase%RandomSeed%Seed(:,i)
                end do
#endif
            end if
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%RandomSeed%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "outputColumnWidth"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%OutputColumnWidth%val
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%OutputColumnWidth%desc )


            !block
            !use Constants_mod, only: TAB
            !character(:), allocatable :: delimiter
            !delimiter = SpecBase%OutputDelimiter%val
            !if (SpecBase%OutputDelimiter%val==TAB) delimiter = "\t"
            !if (SpecBase%OutputDelimiter%val=="\t") delimiter = "\\t"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "outputDelimiter"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%OutputDelimiter%val
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%OutputDelimiter%desc )
            !end block

            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "outputRealPrecision"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%OutputRealPrecision%val
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%OutputRealPrecision%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "chainFileFormat"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%ChainFileFormat%val
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%ChainFileFormat%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "restartFileFormat"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%RestartFileFormat%val
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%RestartFileFormat%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "progressReportPeriod"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%ProgressReportPeriod%val
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%ProgressReportPeriod%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "maxNumDomainCheckToWarn"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%MaxNumDomainCheckToWarn%val
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%MaxNumDomainCheckToWarn%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "maxNumDomainCheckToStop"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecBase%MaxNumDomainCheckToStop%val
            if (SpecBase%SilentModeRequested%isFalse) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecBase%MaxNumDomainCheckToStop%desc )


        end if


    end subroutine reportValues

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine checkForSanity(SpecBase,Err,methodName)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Err_mod, only: Err_type
        implicit none
        class(SpecBase_type), intent(in)    :: SpecBase
        type(Err_type), intent(inout)       :: Err
        character(*), intent(in)            :: methodName
        call SpecBase%ChainFileFormat           %checkForSanity(Err,methodName)
        call SpecBase%OutputDelimiter           %checkForSanity(Err,methodName)
        call SpecBase%DomainLowerLimitVec       %checkForSanity(Err)
        call SpecBase%DomainUpperLimitVec       %checkForSanity(Err,SpecBase%DomainLowerLimitVec%Val)
        call SpecBase%RestartFileFormat         %checkForSanity(Err,methodName)
        call SpecBase%OutputColumnWidth         %checkForSanity(Err,methodName, SpecBase%OutputRealPrecision%val)
        call SpecBase%OutputRealPrecision       %checkForSanity(Err,methodName)
        call SpecBase%ParallelizationModel      %checkForSanity(Err,methodName)
        call SpecBase%ProgressReportPeriod      %checkForSanity(Err,methodName)
        call SpecBase%TargetAcceptanceRate      %checkForSanity(Err)
        call SpecBase%MaxNumDomainCheckToWarn   %checkForSanity(Err,methodName)
        call SpecBase%MaxNumDomainCheckToStop   %checkForSanity(Err,methodName)
    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecBase_mod