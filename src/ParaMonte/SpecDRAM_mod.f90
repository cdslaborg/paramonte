!**********************************************************************************************************************************
!**********************************************************************************************************************************
!
!  ParaMonte: plain powerful parallel Monte Carlo library.
!
!  Copyright (C) 2012-present, The Computational Data Science Lab
!
!  This file is part of ParaMonte library. 
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

module SpecDRAM_mod

    ! ParaDRAM Spec variable types
    use SpecDRAM_ScaleFactor_mod                    , only: ScaleFactor_type
    use SpecDRAM_ProposalModel_mod                  , only: ProposalModel_type
    use SpecDRAM_ProposalStartCovMat_mod            , only: ProposalStartCovMat_type
    use SpecDRAM_ProposalStartCorMat_mod            , only: ProposalStartCorMat_type
    use SpecDRAM_ProposalStartStdVec_mod            , only: ProposalStartStdVec_type
    use SpecDRAM_AdaptiveUpdateCount_mod            , only: AdaptiveUpdateCount_type
    use SpecDRAM_AdaptiveUpdatePeriod_mod           , only: AdaptiveUpdatePeriod_type
    use SpecDRAM_GreedyAdaptationCount_mod          , only: GreedyAdaptationCount_type
    use SpecDRAM_DelayedRejectionCount_mod          , only: DelayedRejectionCount_type
    use SpecDRAM_BurninAdaptationMeasure_mod        , only: BurninAdaptationMeasure_type
    use SpecDRAM_DelayedRejectionScaleFactorVec_mod , only: DelayedRejectionScaleFactorVec_type

    ! ParaDRAM namelist variables
    use SpecDRAM_ScaleFactor_mod                    , only: scaleFactor
    use SpecDRAM_ProposalModel_mod                  , only: proposalModel
    use SpecDRAM_AdaptiveUpdateCount_mod            , only: adaptiveUpdateCount
    use SpecDRAM_AdaptiveUpdatePeriod_mod           , only: adaptiveUpdatePeriod
    use SpecDRAM_GreedyAdaptationCount_mod          , only: greedyAdaptationCount
    use SpecDRAM_DelayedRejectionCount_mod          , only: delayedRejectionCount
    use SpecDRAM_BurninAdaptationMeasure_mod        , only: burninAdaptationMeasure
    use SpecDRAM_ProposalStartCovMat_mod            , only: ProposalStartCovMat
    use SpecDRAM_ProposalStartCorMat_mod            , only: ProposalStartCorMat
    use SpecDRAM_ProposalStartStdVec_mod            , only: ProposalStartStdVec
    use SpecDRAM_DelayedRejectionScaleFactorVec_mod , only: delayedRejectionScaleFactorVec

    implicit none

    type                                            :: SpecDRAM_type
        type(ScaleFactor_type)                      :: ScaleFactor
        type(ProposalModel_type)                    :: ProposalModel
        type(ProposalStartCovMat_type)              :: ProposalStartCovMat
        type(ProposalStartCorMat_type)              :: ProposalStartCorMat
        type(ProposalStartStdVec_type)              :: ProposalStartStdVec
        type(AdaptiveUpdateCount_type)              :: AdaptiveUpdateCount
        type(AdaptiveUpdatePeriod_type)             :: AdaptiveUpdatePeriod
        type(GreedyAdaptationCount_type)            :: GreedyAdaptationCount
        type(DelayedRejectionCount_type)            :: DelayedRejectionCount
        type(BurninAdaptationMeasure_type)          :: BurninAdaptationMeasure
        type(DelayedRejectionScaleFactorVec_type)   :: DelayedRejectionScaleFactorVec
    contains
        procedure, pass                             :: nullifyNameListVar
        procedure, pass                             :: setFromInputFile
        procedure, pass                             :: setFromInputArgs
        procedure, pass                             :: reportValues
        procedure, pass                             :: checkForSanity
    end type SpecDRAM_type

    interface SpecDRAM_type
        module procedure                            :: constructSpecDRAM
    end interface SpecDRAM_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

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
        SpecDRAM%ScaleFactor                    = ScaleFactor_type                      (nd,methodName)
        SpecDRAM%ProposalModel                  = ProposalModel_type                    ()
        SpecDRAM%ProposalStartCovMat            = ProposalStartCovMat_type              (nd,methodName)
        SpecDRAM%ProposalStartCorMat            = ProposalStartCorMat_type              (nd,methodName)
        SpecDRAM%ProposalStartStdVec            = ProposalStartStdVec_type              (nd,methodName)
        SpecDRAM%AdaptiveUpdatePeriod           = AdaptiveUpdatePeriod_type             (nd,methodName)
        ! ATTN: AdaptiveUpdateCount has to be constructed after AdaptiveUpdatePeriod. It depends on it.
        SpecDRAM%AdaptiveUpdateCount            = AdaptiveUpdateCount_type              (methodName) ! ,chainSizeDef,SpecDRAM%AdaptiveUpdatePeriod%def)
        SpecDRAM%GreedyAdaptationCount          = GreedyAdaptationCount_type            (methodName)
        SpecDRAM%DelayedRejectionCount          = DelayedRejectionCount_type            (methodName)
        SpecDRAM%BurninAdaptationMeasure        = BurninAdaptationMeasure_type          (methodName)
        SpecDRAM%DelayedRejectionScaleFactorVec = DelayedRejectionScaleFactorVec_type   (nd,methodName)
    end function constructSpecDRAM

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar( SpecDRAM , nd )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        use Constants_mod, only: IK
        implicit none
        class(SpecDRAM_type), intent(in)    :: SpecDRAM
        integer(IK), intent(in)             :: nd
        ! nullify SpecDRAM variables to be read form the input namelist file
        call SpecDRAM%ScaleFactor%nullifyNameListVar()
        call SpecDRAM%ProposalModel%nullifyNameListVar()
        call SpecDRAM%ProposalStartCovMat%nullifyNameListVar(nd)
        call SpecDRAM%ProposalStartCorMat%nullifyNameListVar(nd)
        call SpecDRAM%ProposalStartStdVec%nullifyNameListVar(nd)
        call SpecDRAM%AdaptiveUpdateCount%nullifyNameListVar()
        call SpecDRAM%AdaptiveUpdatePeriod%nullifyNameListVar()
        call SpecDRAM%GreedyAdaptationCount%nullifyNameListVar()
        call SpecDRAM%DelayedRejectionCount%nullifyNameListVar()
        call SpecDRAM%BurninAdaptationMeasure%nullifyNameListVar()
        call SpecDRAM%DelayedRejectionScaleFactorVec%nullifyNameListVar()
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

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
        call SpecDRAM%ScaleFactor                   %set(scaleFactor)
        call SpecDRAM%ProposalModel                 %set(trim(adjustl(proposalModel)))
        call SpecDRAM%ProposalStartCovMat           %set(ProposalStartCovMat)
        call SpecDRAM%ProposalStartCorMat           %set(ProposalStartCorMat)
        call SpecDRAM%ProposalStartStdVec           %set(ProposalStartStdVec)
        call SpecDRAM%AdaptiveUpdateCount           %set(adaptiveUpdateCount)
        call SpecDRAM%AdaptiveUpdatePeriod          %set(adaptiveUpdatePeriod)
        call SpecDRAM%GreedyAdaptationCount         %set(greedyAdaptationCount)
        call SpecDRAM%DelayedRejectionCount         %set(delayedRejectionCount)
        call SpecDRAM%BurninAdaptationMeasure       %set(burninAdaptationMeasure)
        call SpecDRAM%DelayedRejectionScaleFactorVec%set(delayedRejectionScaleFactorVec,delayedRejectionCount)
    end subroutine setFromInputFile

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setFromInputArgs ( SpecDRAM &
                                ! input arguments to the specific ParaDRAM routine
                                , scaleFactor                       &
                                , proposalModel                     &
                                , proposalStartCovMat               &
                                , proposalStartCorMat               &
                                , proposalStartStdVec               &
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
        character(*), intent(in), optional  :: scaleFactor
        character(*), intent(in), optional  :: proposalModel
        real(RK)    , intent(in), optional  :: proposalStartCovMat(:,:)
        real(RK)    , intent(in), optional  :: proposalStartCorMat(:,:)
        real(RK)    , intent(in), optional  :: proposalStartStdVec(:)
        integer(IK) , intent(in), optional  :: adaptiveUpdateCount
        integer(IK) , intent(in), optional  :: adaptiveUpdatePeriod
        integer(IK) , intent(in), optional  :: greedyAdaptationCount
        integer(IK) , intent(in), optional  :: delayedRejectionCount
        real(RK)    , intent(in), optional  :: burninAdaptationMeasure
        real(RK)    , intent(in), optional  :: delayedRejectionScaleFactorVec(:)

        if (present(scaleFactor))                       call SpecDRAM%ScaleFactor                   %set(scaleFactor)
        if (present(proposalModel))                     call SpecDRAM%ProposalModel                 %set(trim(adjustl(proposalModel)))
        if (present(proposalStartCovMat))               call SpecDRAM%ProposalStartCovMat           %set(proposalStartCovMat)
        if (present(proposalStartCorMat))               call SpecDRAM%ProposalStartCorMat           %set(proposalStartCorMat)
        if (present(proposalStartStdVec))               call SpecDRAM%ProposalStartStdVec           %set(proposalStartStdVec)
        if (present(adaptiveUpdateCount))               call SpecDRAM%AdaptiveUpdateCount           %set(adaptiveUpdateCount)
        if (present(adaptiveUpdatePeriod))              call SpecDRAM%AdaptiveUpdatePeriod          %set(adaptiveUpdatePeriod)
        if (present(greedyAdaptationCount))             call SpecDRAM%GreedyAdaptationCount         %set(greedyAdaptationCount)
        if (present(delayedRejectionCount))             call SpecDRAM%DelayedRejectionCount         %set(delayedRejectionCount)
        if (present(burninAdaptationMeasure))           call SpecDRAM%BurninAdaptationMeasure       %set(burninAdaptationMeasure)
        if (present(delayedRejectionScaleFactorVec))    call SpecDRAM%delayedRejectionScaleFactorVec%set(delayedRejectionScaleFactorVec,SpecDRAM%DelayedRejectionCount%val)

    end subroutine setFromInputArgs

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine reportValues(SpecDRAM,prefix,outputUnit,isMasterImage,methodName,splashModeRequested)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: reportValues
#endif
        use Decoration_mod, only: TAB
        use Constants_mod, only: IK, RK
        use Err_mod, only: note, informUser
        implicit none
        class(SpecDRAM_type), intent(in)    :: SpecDRAM
        character(*), intent(in)            :: prefix, methodName
        integer(IK), intent(in)             :: outputUnit
        logical, intent(in)                 :: isMasterImage, splashModeRequested
        character(:), allocatable           :: formatStr, formatVal
        integer(IK)                         :: ndim, i
        real(RK), allocatable               :: Row(:)

        formatStr = "(*(g0,' '))"
        formatVal = "('" // TAB // TAB // "',*(g0,' '))"


        if (isMasterImage) then

            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "adaptiveUpdatePeriod"
            write(outputUnit,formatStr)
            write(outputUnit,formatVal) SpecDRAM%AdaptiveUpdatePeriod%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecDRAM%AdaptiveUpdatePeriod%desc )


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "adaptiveUpdateCount"
            write(outputUnit,formatStr)
            write(outputUnit,formatVal) SpecDRAM%AdaptiveUpdateCount%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecDRAM%AdaptiveUpdateCount%desc )


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "greedyAdaptationCount"
            write(outputUnit,formatStr)
            write(outputUnit,formatVal) SpecDRAM%GreedyAdaptationCount%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecDRAM%GreedyAdaptationCount%desc )


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "burninAdaptationMeasure"
            write(outputUnit,formatStr)
            write(outputUnit,formatVal) SpecDRAM%BurninAdaptationMeasure%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecDRAM%BurninAdaptationMeasure%desc )


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "delayedRejectionCount"
            write(outputUnit,formatStr)
            write(outputUnit,formatVal) SpecDRAM%DelayedRejectionCount%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecDRAM%DelayedRejectionCount%desc )


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "delayedRejectionScaleFactorVec"
            write(outputUnit,formatStr)
            if ( size(SpecDRAM%DelayedRejectionScaleFactorVec%Val) == 0 ) then
                write(outputUnit,formatVal) "UNDEFINED"
            else
                do i = 1, size(SpecDRAM%DelayedRejectionScaleFactorVec%Val)
                    write(outputUnit,formatVal) SpecDRAM%DelayedRejectionScaleFactorVec%Val(i)
                end do
            end if
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecDRAM%DelayedRejectionScaleFactorVec%desc )


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "scaleFactor"
            write(outputUnit,formatStr)
            write(outputUnit,formatVal) SpecDRAM%ScaleFactor%str
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecDRAM%ScaleFactor%desc )


            !***********************************************************************************************************************
            ! proposal distribution
            !***********************************************************************************************************************

            block
                use Decoration_mod, only: writeDecoratedText
                call writeDecoratedText ( text = "\n" // methodName // " proposal specifications\n" &
                                        , marginTop = 1     &
                                        , marginBot = 1     &
                                        , newline = "\n"    &
                                        , outputUnit = outputUnit )
            end block


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "proposalModel"
            write(outputUnit,formatStr)
            write(outputUnit,formatVal) SpecDRAM%ProposalModel%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecDRAM%ProposalModel%desc )


            ndim = size(SpecDRAM%ProposalStartCovMat%Val(:,1))
            allocate( Row(ndim) )


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "proposalStartCovMat"
            write(outputUnit,formatStr)
            if ( SpecDRAM%ProposalStartCovMat%isPresent ) then
                ! User has provided the Start Covariance Matrix
                do i = 1,ndim
                    Row = SpecDRAM%ProposalStartCovMat%Val(i,:)
                    write(outputUnit,formatVal) Row
                end do
            else
                ! User has not provided the Start Covariance Matrix
                call informUser ( prefix = TAB // TAB       &
                                , outputUnit = outputUnit   &
                                , newline = "\n"            &
                                , marginBot = 0_IK          &
                                , msg = "UNDEFINED. It will be constructed from the Correlation Matrix (ProposalStartCorMat) &
                                        &and the Standard Deviation vector (ProposalStartStdVec)." &
                                )
            end if
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecDRAM%ProposalStartCovMat%desc )


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "proposalStartCorMat"
            write(outputUnit,formatStr)
            do i = 1,ndim
                Row = SpecDRAM%ProposalStartCorMat%Val(i,:)
                write(outputUnit,formatVal) Row
            end do
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecDRAM%ProposalStartCorMat%desc )


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "proposalStartStdVec"
            write(outputUnit,formatStr)
            do i = 1,ndim
                write(outputUnit,formatVal) SpecDRAM%ProposalStartStdVec%Val(i)
            end do
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecDRAM%ProposalStartStdVec%desc )

        end if


    end subroutine reportValues

!***********************************************************************************************************************************
!***********************************************************************************************************************************

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
        call SpecDRAM%ScaleFactor%checkForSanity                    (Err,methodName)
        call SpecDRAM%ProposalModel%checkForSanity                  (Err,methodName)
        call SpecDRAM%AdaptiveUpdateCount%checkForSanity            (Err,methodName)
        call SpecDRAM%AdaptiveUpdatePeriod%checkForSanity           (Err,methodName)
        call SpecDRAM%GreedyAdaptationCount%checkForSanity          (Err,methodName)
        call SpecDRAM%DelayedRejectionCount%checkForSanity          (Err,methodName)
        call SpecDRAM%BurninAdaptationMeasure%checkForSanity        (Err,methodName)
        call SpecDRAM%ProposalStartCovMat%checkForSanity            (Err,methodName,nd)
        call SpecDRAM%ProposalStartCorMat%checkForSanity            (Err,methodName,nd)
        call SpecDRAM%ProposalStartStdVec%checkForSanity            (Err,methodName,nd)
        call SpecDRAM%DelayedRejectionScaleFactorVec%checkForSanity (Err,methodName,SpecDRAM%DelayedRejectionCount%val)
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecDRAM_mod