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
!> This module contains the classes and procedures for setting up the attributes of samplers of class [ParaNest_type](@ref paranest_mod::paranest_type).
!> For more information, see the description of this attributes in the body of the corresponding modules.
!> \author Amir Shahmoradi

module SpecNest_mod

    ! ParaNest Spec variable types
    use SpecNest_Tightness_mod              , only: Tightness_type
    use SpecNest_Tolerance_mod              , only: Tolerance_type
    use SpecNest_ScaleFactor_mod            , only: ScaleFactor_type
    use SpecNest_ProposalModel_mod          , only: ProposalModel_type
    use SpecNest_LiveSampleSize_mod         , only: LiveSampleSize_type
    use SpecNest_InclusionFraction_mod      , only: InclusionFraction_type
    use SpecNest_AdaptiveUpdateCount_mod    , only: AdaptiveUpdateCount_type
    use SpecNest_AdaptiveUpdatePeriod_mod   , only: AdaptiveUpdatePeriod_type
    use SpecNest_MahalSqWeightExponent_mod  , only: MahalSqWeightExponent_type
    use SpecNest_StabilizationRequested_mod , only: StabilizationRequested_type
    use SpecNest_MaxAllowedKmeansFailure_mod, only: MaxAllowedKmeansFailure_type
    use SpecNest_MaxAllowedMinVolFailure_mod, only: MaxAllowedMinVolFailure_type
    use SpecNest_MaxKvolumeLoopRecursion_mod, only: MaxKvolumeLoopRecursion_type

    ! ParaNest namelist variables
    use SpecNest_Tightness_mod              , only: tightness
    use SpecNest_Tolerance_mod              , only: tolerance
    use SpecNest_ScaleFactor_mod            , only: scaleFactor
    use SpecNest_ProposalModel_mod          , only: proposalModel
    use SpecNest_LiveSampleSize_mod         , only: liveSampleSize
    use SpecNest_InclusionFraction_mod      , only: inclusionFraction
    use SpecNest_AdaptiveUpdateCount_mod    , only: adaptiveUpdateCount
    use SpecNest_AdaptiveUpdatePeriod_mod   , only: adaptiveUpdatePeriod
    use SpecNest_MahalSqWeightExponent_mod  , only: mahalSqWeightExponent
    use SpecNest_StabilizationRequested_mod , only: stabilizationRequested
    use SpecNest_MaxAllowedKmeansFailure_mod, only: MaxAllowedKmeansFailure
    use SpecNest_MaxAllowedMinVolFailure_mod, only: maxAllowedMinVolFailure
    use SpecNest_MaxKvolumeLoopRecursion_mod, only: maxKvolumeLoopRecursion

    implicit none

    type                                    :: SpecNest_type
        type(Tightness_type)                :: Tightness
        type(Tolerance_type)                :: Tolerance
        type(ScaleFactor_type)              :: ScaleFactor
        type(ProposalModel_type)            :: ProposalModel
        type(LiveSampleSize_type)           :: LiveSampleSize
        type(InclusionFraction_type)        :: InclusionFraction
        type(AdaptiveUpdateCount_type)      :: AdaptiveUpdateCount
        type(AdaptiveUpdatePeriod_type)     :: AdaptiveUpdatePeriod
        type(MahalSqWeightExponent_type)    :: MahalSqWeightExponent
        type(StabilizationRequested_type)   :: StabilizationRequested
        type(MaxAllowedKmeansFailure_type)  :: MaxAllowedKmeansFailure
        type(MaxAllowedMinVolFailure_type)  :: MaxAllowedMinVolFailure
        type(MaxKvolumeLoopRecursion_type)  :: MaxKvolumeLoopRecursion
    contains
        procedure, pass                     :: nullifyNameListVar
        procedure, pass                     :: setFromInputFile
        procedure, pass                     :: setFromInputArgs
        procedure, pass                     :: reportValues
        procedure, pass                     :: checkForSanity
    end type SpecNest_type

    interface SpecNest_type
        module procedure                    :: construct
    end interface SpecNest_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function construct ( nd &
                            , methodName &
                            ) result(SpecNest)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct
#endif
        use Constants_mod, only: IK
        implicit none
        integer(IK), intent(in)     :: nd
        character(*), intent(in)    :: methodName
        type(SpecNest_type)         :: SpecNest
        SpecNest%Tightness                  = Tightness_type(methodName)
        SpecNest%Tolerance                  = Tolerance_type(methodName)
        SpecNest%ScaleFactor                = ScaleFactor_type(methodName)
        SpecNest%ProposalModel              = ProposalModel_type(methodName)
        SpecNest%LiveSampleSize             = LiveSampleSize_type(nd,methodName)
        SpecNest%InclusionFraction          = InclusionFraction_type(methodName)
        SpecNest%AdaptiveUpdateCount        = AdaptiveUpdateCount_type(methodName)
        SpecNest%AdaptiveUpdatePeriod       = AdaptiveUpdatePeriod_type(nd,methodName)
        SpecNest%MahalSqWeightExponent      = MahalSqWeightExponent_type(methodName)
        SpecNest%StabilizationRequested     = StabilizationRequested_type(methodName)
        SpecNest%MaxAllowedKmeansFailure    = MaxAllowedKmeansFailure_type(methodName)
        SpecNest%MaxAllowedMinVolFailure    = MaxAllowedMinVolFailure_type(methodName)
        SpecNest%MaxKvolumeLoopRecursion    = MaxKvolumeLoopRecursion_type(methodName)
    end function construct

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(SpecNest, nd)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        use Constants_mod, only: IK
        implicit none
        class(SpecNest_type), intent(in)    :: SpecNest
        integer(IK), intent(in)             :: nd
        call SpecNest%Tightness                 %nullifyNameListVar()
        call SpecNest%Tolerance                 %nullifyNameListVar()
        call SpecNest%ScaleFactor               %nullifyNameListVar()
        call SpecNest%ProposalModel             %nullifyNameListVar()
        call SpecNest%LiveSampleSize            %nullifyNameListVar()
        call SpecNest%InclusionFraction         %nullifyNameListVar()
        call SpecNest%AdaptiveUpdateCount       %nullifyNameListVar()
        call SpecNest%AdaptiveUpdatePeriod      %nullifyNameListVar()
        call SpecNest%MahalSqWeightExponent     %nullifyNameListVar()
        call SpecNest%StabilizationRequested    %nullifyNameListVar()
        call SpecNest%MaxAllowedKmeansFailure   %nullifyNameListVar()
        call SpecNest%MaxAllowedMinVolFailure   %nullifyNameListVar()
        call SpecNest%MaxKvolumeLoopRecursion   %nullifyNameListVar()
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure subroutine setFromInputFile( SpecNest, Err )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFromInputFile
#endif
        use Constants_mod, only: IK, RK
        use Err_mod, only: Err_type
        implicit none
        class(SpecNest_type), intent(inout) :: SpecNest
        type(Err_type), intent(out)         :: Err
        Err%occurred = .false.
        Err%msg = ""
        call SpecNest%Tightness                 %set(tightness              )
        call SpecNest%Tolerance                 %set(tolerance              )
        call SpecNest%ScaleFactor               %set(scaleFactor            )
        call SpecNest%ProposalModel             %set(proposalModel         )
        call SpecNest%LiveSampleSize            %set(liveSampleSize         )
        call SpecNest%InclusionFraction         %set(inclusionFraction      )
        call SpecNest%AdaptiveUpdateCount       %set(adaptiveUpdateCount    )
        call SpecNest%AdaptiveUpdatePeriod      %set(adaptiveUpdatePeriod   )
        call SpecNest%MahalSqWeightExponent     %set(mahalSqWeightExponent  )
        call SpecNest%StabilizationRequested    %set(stabilizationRequested )
        call SpecNest%MaxAllowedKmeansFailure   %set(maxAllowedKmeansFailure)
        call SpecNest%MaxAllowedMinVolFailure   %set(maxAllowedMinVolFailure)
        call SpecNest%MaxKvolumeLoopRecursion   %set(maxKvolumeLoopRecursion)
    end subroutine setFromInputFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure subroutine setFromInputArgs( SpecNest                  &
                                    , tightness                 &
                                    , tolerance                 &
                                    , scaleFactor               &
                                    , proposalModel             &
                                    , liveSampleSize            &
                                    , inclusionFraction         &
                                    , adaptiveUpdateCount       &
                                    , adaptiveUpdatePeriod      &
                                    , mahalSqWeightExponent     &
                                    , stabilizationRequested    &
                                    , maxAllowedKmeansFailure   &
                                    , maxAllowedMinVolFailure   &
                                    , maxKvolumeLoopRecursion   &
                                    )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFromInputArgs
#endif
        use Constants_mod, only: IK, RK
        implicit none
        class(SpecNest_type), intent(inout) :: SpecNest

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

        if (present(tightness               ))    call SpecNest%Tightness               %set(tightness              )
        if (present(tolerance               ))    call SpecNest%Tolerance               %set(tolerance              )
        if (present(scaleFactor             ))    call SpecNest%ScaleFactor             %set(scaleFactor            )
        if (present(proposalModel          ))    call SpecNest%ProposalModel          %set(proposalModel         )
        if (present(liveSampleSize          ))    call SpecNest%LiveSampleSize          %set(liveSampleSize         )
        if (present(inclusionFraction       ))    call SpecNest%InclusionFraction       %set(inclusionFraction      )
        if (present(adaptiveUpdateCount     ))    call SpecNest%AdaptiveUpdateCount     %set(adaptiveUpdateCount    )
        if (present(adaptiveUpdatePeriod    ))    call SpecNest%AdaptiveUpdatePeriod    %set(adaptiveUpdatePeriod   )
        if (present(mahalSqWeightExponent   ))    call SpecNest%MahalSqWeightExponent   %set(mahalSqWeightExponent  )
        if (present(stabilizationRequested  ))    call SpecNest%StabilizationRequested  %set(stabilizationRequested )
        if (present(maxAllowedKmeansFailure ))    call SpecNest%MaxAllowedKmeansFailure %set(maxAllowedKmeansFailure)
        if (present(maxAllowedMinVolFailure ))    call SpecNest%MaxAllowedMinVolFailure %set(maxAllowedMinVolFailure)
        if (present(maxKvolumeLoopRecursion ))    call SpecNest%MaxKvolumeLoopRecursion %set(maxKvolumeLoopRecursion)

    end subroutine setFromInputArgs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine reportValues ( SpecNest              &
                            , prefix                &
                            , outputUnit            &
                            , isLeaderImage         &
                            , splashModeRequested   &
                            )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: reportValues
#endif
        use Decoration_mod, only: GENERIC_OUTPUT_FORMAT
        use Decoration_mod, only: GENERIC_TABBED_FORMAT
        use Constants_mod, only: IK, RK, UNDEFINED
        use Err_mod, only: note
        implicit none
        class(SpecNest_type), intent(in)    :: SpecNest
        character(*), intent(in)            :: prefix
        integer(IK), intent(in)             :: outputUnit
        logical, intent(in)                 :: isLeaderImage, splashModeRequested
        integer(IK)                         :: i

        if (isLeaderImage) then


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "tolerance"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecNest%Tolerance%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecNest%Tolerance%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "tightness"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecNest%Tightness%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecNest%Tightness%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "scaleFactor"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecNest%ScaleFactor%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecNest%ScaleFactor%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "proposalModel"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecNest%ProposalModel%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecNest%ProposalModel%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "liveSampleSize"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecNest%LiveSampleSize%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecNest%LiveSampleSize%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "inclusionFraction"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecNest%InclusionFraction%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecNest%InclusionFraction%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "adaptiveUpdateCount"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecNest%AdaptiveUpdateCount%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecNest%AdaptiveUpdateCount%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "adaptiveUpdatePeriod"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecNest%AdaptiveUpdatePeriod%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecNest%AdaptiveUpdatePeriod%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "mahalSqWeightExponent"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecNest%MahalSqWeightExponent%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecNest%MahalSqWeightExponent%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "stabilizationRequested"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecNest%StabilizationRequested%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecNest%StabilizationRequested%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "maxAllowedKmeansFailure"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecNest%MaxAllowedKmeansFailure%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecNest%MaxAllowedKmeansFailure%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "maxAllowedMinVolFailure"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecNest%maxAllowedMinVolFailure%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecNest%MaxAllowedMinVolFailure%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "maxKvolumeLoopRecursion"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecNest%MaxKvolumeLoopRecursion%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecNest%MaxKvolumeLoopRecursion%desc )

        end if

    end subroutine reportValues

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure subroutine checkForSanity(SpecNest,Err,methodName,nd)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Constants_mod, only: IK, RK
        use Err_mod, only: Err_type
        implicit none
        class(SpecNest_type), intent(inout) :: SpecNest
        integer(IK), intent(in)             :: nd
        character(*), intent(in)            :: methodName
        type(Err_type), intent(inout)       :: Err
        call SpecNest%Tightness                 %checkForSanity (Err,methodName)
        call SpecNest%Tolerance                 %checkForSanity (Err,methodName)
        call SpecNest%ProposalModel            %checkForSanity (Err,methodName)
        call SpecNest%LiveSampleSize            %checkForSanity (Err,methodName,nd)
        call SpecNest%InclusionFraction         %checkForSanity (Err,methodName)
        call SpecNest%AdaptiveUpdateCount       %checkForSanity (Err,methodName)
        call SpecNest%AdaptiveUpdatePeriod      %checkForSanity (Err,methodName)
        call SpecNest%MahalSqWeightExponent     %checkForSanity (Err,methodName)
        call SpecNest%MaxAllowedKmeansFailure   %checkForSanity (Err,methodName)
        call SpecNest%MaxAllowedMinVolFailure   %checkForSanity (Err,methodName)
        call SpecNest%MaxKvolumeLoopRecursion   %checkForSanity (Err,methodName)
    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecNest_mod ! LCOV_EXCL_LINE