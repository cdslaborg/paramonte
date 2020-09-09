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

module SpecDRAM_mod

    ! ParaDRAM Spec variable types
    use SpecDRAM_AdaptiveUpdateCount_mod                , only: AdaptiveUpdateCount_type
    use SpecDRAM_AdaptiveUpdatePeriod_mod               , only: AdaptiveUpdatePeriod_type
    use SpecDRAM_GreedyAdaptationCount_mod              , only: GreedyAdaptationCount_type
    use SpecDRAM_DelayedRejectionCount_mod              , only: DelayedRejectionCount_type
    use SpecDRAM_BurninAdaptationMeasure_mod            , only: BurninAdaptationMeasure_type
    use SpecDRAM_DelayedRejectionScaleFactorVec_mod     , only: DelayedRejectionScaleFactorVec_type

    ! ParaDRAM namelist variables
    use SpecDRAM_AdaptiveUpdateCount_mod                , only: adaptiveUpdateCount
    use SpecDRAM_AdaptiveUpdatePeriod_mod               , only: adaptiveUpdatePeriod
    use SpecDRAM_GreedyAdaptationCount_mod              , only: greedyAdaptationCount
    use SpecDRAM_DelayedRejectionCount_mod              , only: delayedRejectionCount
    use SpecDRAM_BurninAdaptationMeasure_mod            , only: burninAdaptationMeasure
    use SpecDRAM_DelayedRejectionScaleFactorVec_mod     , only: delayedRejectionScaleFactorVec

    implicit none

    type                                                :: SpecDRAM_type
        type(AdaptiveUpdateCount_type)                  :: AdaptiveUpdateCount
        type(AdaptiveUpdatePeriod_type)                 :: AdaptiveUpdatePeriod
        type(GreedyAdaptationCount_type)                :: GreedyAdaptationCount
        type(DelayedRejectionCount_type)                :: DelayedRejectionCount
        type(BurninAdaptationMeasure_type)              :: BurninAdaptationMeasure
        type(DelayedRejectionScaleFactorVec_type)       :: DelayedRejectionScaleFactorVec
    contains
        procedure, pass                                 :: nullifyNameListVar
        procedure, pass                                 :: setFromInputFile
        procedure, pass                                 :: setFromInputArgs
        procedure, pass                                 :: reportValues
        procedure, pass                                 :: checkForSanity
    end type SpecDRAM_type

    interface SpecDRAM_type
        module procedure                                :: constructSpecDRAM
    end interface SpecDRAM_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructSpecDRAM  ( nd &
                                , methodName &
                                !, chainSizeDef
                                ) result(SpecDRAM)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructSpecDRAM
#endif
        use Constants_mod, only: IK
        implicit none
        integer(IK), intent(in)     :: nd
        character(*), intent(in)    :: methodName
       !integer(IK), intent(in)     :: chainSizeDef
        type(SpecDRAM_type)         :: SpecDRAM
        SpecDRAM%AdaptiveUpdatePeriod                   = AdaptiveUpdatePeriod_type             (nd,methodName)
        ! ATTN: AdaptiveUpdateCount has to be constructed after AdaptiveUpdatePeriod. It depends on it.
        SpecDRAM%AdaptiveUpdateCount                    = AdaptiveUpdateCount_type              (methodName) ! ,chainSizeDef,SpecDRAM%AdaptiveUpdatePeriod%def)
        SpecDRAM%GreedyAdaptationCount                  = GreedyAdaptationCount_type            (methodName)
        SpecDRAM%DelayedRejectionCount                  = DelayedRejectionCount_type            (methodName)
        SpecDRAM%BurninAdaptationMeasure                = BurninAdaptationMeasure_type          (methodName)
        SpecDRAM%DelayedRejectionScaleFactorVec         = DelayedRejectionScaleFactorVec_type   (nd,methodName)
    end function constructSpecDRAM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar( SpecDRAM , nd )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        use Constants_mod, only: IK
        implicit none
        class(SpecDRAM_type), intent(in)    :: SpecDRAM
        integer(IK), intent(in)             :: nd
        ! nullify SpecDRAM variables to be read form the input namelist file
        call SpecDRAM%AdaptiveUpdateCount           %nullifyNameListVar()
        call SpecDRAM%AdaptiveUpdatePeriod          %nullifyNameListVar()
        call SpecDRAM%GreedyAdaptationCount         %nullifyNameListVar()
        call SpecDRAM%DelayedRejectionCount         %nullifyNameListVar()
        call SpecDRAM%BurninAdaptationMeasure       %nullifyNameListVar()
        call SpecDRAM%DelayedRejectionScaleFactorVec%nullifyNameListVar()
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setFromInputFile( SpecDRAM, Err )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setFromInputFile
#endif
        use Constants_mod, only: IK, RK
        use Err_mod, only: Err_type
        implicit none
        class(SpecDRAM_type), intent(inout) :: SpecDRAM
        type(Err_type), intent(out)         :: Err
        Err%occurred = .false.
        Err%msg = ""
        call SpecDRAM%AdaptiveUpdateCount           %set(adaptiveUpdateCount)
        call SpecDRAM%AdaptiveUpdatePeriod          %set(adaptiveUpdatePeriod)
        call SpecDRAM%GreedyAdaptationCount         %set(greedyAdaptationCount)
        call SpecDRAM%DelayedRejectionCount         %set(delayedRejectionCount)
        call SpecDRAM%BurninAdaptationMeasure       %set(burninAdaptationMeasure)
        call SpecDRAM%DelayedRejectionScaleFactorVec%set(delayedRejectionScaleFactorVec,delayedRejectionCount)
    end subroutine setFromInputFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setFromInputArgs ( SpecDRAM &
                                ! input arguments to the specific ParaDRAM routine
                                , adaptiveUpdateCount               &
                                , adaptiveUpdatePeriod              &
                                , greedyAdaptationCount             &
                                , delayedRejectionCount             &
                                , burninAdaptationMeasure           &
                                , delayedRejectionScaleFactorVec    &
                                )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setFromInputArgs
#endif

        use Constants_mod, only: IK, RK
        implicit none
        class(SpecDRAM_type), intent(inout) :: SpecDRAM

        ! ParaDRAM variables
        integer(IK) , intent(in), optional  :: adaptiveUpdateCount
        integer(IK) , intent(in), optional  :: adaptiveUpdatePeriod
        integer(IK) , intent(in), optional  :: greedyAdaptationCount
        integer(IK) , intent(in), optional  :: delayedRejectionCount
        real(RK)    , intent(in), optional  :: burninAdaptationMeasure
        real(RK)    , intent(in), optional  :: delayedRejectionScaleFactorVec(:)

        if (present(adaptiveUpdateCount))               call SpecDRAM%AdaptiveUpdateCount           %set(adaptiveUpdateCount)
        if (present(adaptiveUpdatePeriod))              call SpecDRAM%AdaptiveUpdatePeriod          %set(adaptiveUpdatePeriod)
        if (present(greedyAdaptationCount))             call SpecDRAM%GreedyAdaptationCount         %set(greedyAdaptationCount)
        if (present(delayedRejectionCount))             call SpecDRAM%DelayedRejectionCount         %set(delayedRejectionCount)
        if (present(burninAdaptationMeasure))           call SpecDRAM%BurninAdaptationMeasure       %set(burninAdaptationMeasure)
        if (present(delayedRejectionScaleFactorVec))    call SpecDRAM%delayedRejectionScaleFactorVec%set(delayedRejectionScaleFactorVec,SpecDRAM%DelayedRejectionCount%val)

    end subroutine setFromInputArgs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine reportValues ( SpecDRAM              &
                            , prefix                &
                            , outputUnit            &
                            , isMasterImage         &
                            , splashModeRequested   &
                            )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: reportValues
#endif
        use Decoration_mod, only: GENERIC_OUTPUT_FORMAT
        use Decoration_mod, only: GENERIC_TABBED_FORMAT
        use Constants_mod, only: IK, RK, UNDEFINED
        use Err_mod, only: note
        implicit none
        class(SpecDRAM_type), intent(in)    :: SpecDRAM
        character(*), intent(in)            :: prefix
        integer(IK), intent(in)             :: outputUnit
        logical, intent(in)                 :: isMasterImage, splashModeRequested
        integer(IK)                         :: i

        if (isMasterImage) then

            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "adaptiveUpdatePeriod"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecDRAM%AdaptiveUpdatePeriod%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecDRAM%AdaptiveUpdatePeriod%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "adaptiveUpdateCount"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecDRAM%AdaptiveUpdateCount%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecDRAM%AdaptiveUpdateCount%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "greedyAdaptationCount"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecDRAM%GreedyAdaptationCount%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecDRAM%GreedyAdaptationCount%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "burninAdaptationMeasure"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecDRAM%BurninAdaptationMeasure%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecDRAM%BurninAdaptationMeasure%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "delayedRejectionCount"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecDRAM%DelayedRejectionCount%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecDRAM%DelayedRejectionCount%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "delayedRejectionScaleFactorVec"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            if ( size(SpecDRAM%DelayedRejectionScaleFactorVec%Val) == 0 ) then
                write(outputUnit,GENERIC_TABBED_FORMAT) UNDEFINED
            else
                do i = 1, size(SpecDRAM%DelayedRejectionScaleFactorVec%Val)
                    write(outputUnit,GENERIC_TABBED_FORMAT) SpecDRAM%DelayedRejectionScaleFactorVec%Val(i)
                end do
            end if
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecDRAM%DelayedRejectionScaleFactorVec%desc )

        end if

    end subroutine reportValues

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine checkForSanity(SpecDRAM,Err,methodName,nd)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Constants_mod, only: IK, RK
        use Err_mod, only: Err_type
        implicit none
        class(SpecDRAM_type), intent(inout) :: SpecDRAM
        integer(IK), intent(in)             :: nd
        character(*), intent(in)            :: methodName
        type(Err_type), intent(inout)       :: Err
        call SpecDRAM%AdaptiveUpdateCount%checkForSanity            (Err,methodName)
        call SpecDRAM%AdaptiveUpdatePeriod%checkForSanity           (Err,methodName)
        call SpecDRAM%GreedyAdaptationCount%checkForSanity          (Err,methodName)
        call SpecDRAM%DelayedRejectionCount%checkForSanity          (Err,methodName)
        call SpecDRAM%BurninAdaptationMeasure%checkForSanity        (Err,methodName)
        call SpecDRAM%DelayedRejectionScaleFactorVec%checkForSanity (Err,methodName,SpecDRAM%DelayedRejectionCount%val)
    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecDRAM_mod