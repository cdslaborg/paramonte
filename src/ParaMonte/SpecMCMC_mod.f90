!***********************************************************************************************************************************
!***********************************************************************************************************************************
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
!   we ask you to acknowledge the ParaMonte library's usage
!   in your work (education/research/industry/development/...)
!   by citing the ParaMonte library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!***********************************************************************************************************************************
!***********************************************************************************************************************************

module SpecMCMC_mod

    ! ParaMCMC Spec variable types
    use SpecMCMC_ChainSize_mod                          , only: ChainSize_type
    use SpecMCMC_ScaleFactor_mod                        , only: ScaleFactor_type
    use SpecMCMC_StartPointVec_mod                      , only: StartPointVec_type
    use SpecMCMC_ProposalModel_mod                      , only: ProposalModel_type
    use SpecMCMC_ProposalStartCovMat_mod                , only: ProposalStartCovMat_type
    use SpecMCMC_ProposalStartCorMat_mod                , only: ProposalStartCorMat_type
    use SpecMCMC_ProposalStartStdVec_mod                , only: ProposalStartStdVec_type
    use SpecMCMC_SampleRefinementCount_mod              , only: SampleRefinementCount_type
    use SpecMCMC_SampleRefinementMethod_mod             , only: SampleRefinementMethod_type
    use SpecMCMC_RandomStartPointRequested_mod          , only: RandomStartPointRequested_type
    use SpecMCMC_RandomStartPointDomainLowerLimitVec_mod, only: RandomStartPointDomainLowerLimitVec_type
    use SpecMCMC_RandomStartPointDomainUpperLimitVec_mod, only: RandomStartPointDomainUpperLimitVec_type

    ! ParaMCMC namelist variables
    use SpecMCMC_ChainSize_mod                          , only: ChainSize
    use SpecMCMC_ScaleFactor_mod                        , only: scaleFactor
    use SpecMCMC_StartPointVec_mod                      , only: startPointVec
    use SpecMCMC_ProposalModel_mod                      , only: proposalModel
    use SpecMCMC_ProposalStartCovMat_mod                , only: ProposalStartCovMat
    use SpecMCMC_ProposalStartCorMat_mod                , only: ProposalStartCorMat
    use SpecMCMC_ProposalStartStdVec_mod                , only: ProposalStartStdVec
    use SpecMCMC_SampleRefinementCount_mod              , only: SampleRefinementCount
    use SpecMCMC_SampleRefinementMethod_mod             , only: sampleRefinementMethod
    use SpecMCMC_RandomStartPointRequested_mod          , only: randomStartPointRequested
    use SpecMCMC_RandomStartPointDomainLowerLimitVec_mod, only: randomStartPointDomainLowerLimitVec
    use SpecMCMC_RandomStartPointDomainUpperLimitVec_mod, only: randomStartPointDomainUpperLimitVec

    implicit none

    type                                                :: SpecMCMC_type
        type(ChainSize_type)                            :: ChainSize
        type(ScaleFactor_type)                          :: ScaleFactor
        type(StartPointVec_type)                        :: StartPointVec
        type(ProposalModel_type)                        :: ProposalModel
        type(ProposalStartCovMat_type)                  :: ProposalStartCovMat
        type(ProposalStartCorMat_type)                  :: ProposalStartCorMat
        type(ProposalStartStdVec_type)                  :: ProposalStartStdVec
        type(SampleRefinementCount_type)                :: SampleRefinementCount
        type(SampleRefinementMethod_type)               :: SampleRefinementMethod
        type(RandomStartPointRequested_type)            :: randomStartPointRequested
        type(RandomStartPointDomainLowerLimitVec_type)  :: RandomStartPointDomainLowerLimitVec
        type(RandomStartPointDomainUpperLimitVec_type)  :: RandomStartPointDomainUpperLimitVec
    contains
        procedure, pass                                 :: nullifyNameListVar
        procedure, pass                                 :: setFromInputFile
        procedure, pass                                 :: setFromInputArgs
        procedure, pass                                 :: checkForSanity
        procedure, pass                                 :: reportValues
    end type SpecMCMC_type

    interface SpecMCMC_type
        module procedure                                :: constructSpecMCMC
    end interface SpecMCMC_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructSpecMCMC  ( nd &
                                , methodName &
                                ) result(SpecMCMC)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructSpecMCMC
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK) , intent(in)    :: nd
        character(*), intent(in)    :: methodName
        type(SpecMCMC_type)         :: SpecMCMC
        SpecMCMC%ChainSize                              = ChainSize_type                            (methodName)
        SpecMCMC%ScaleFactor                            = ScaleFactor_type                          (nd,methodName)
        SpecMCMC%StartPointVec                          = StartPointVec_type                        ()
        SpecMCMC%ProposalModel                          = ProposalModel_type                        ()
        SpecMCMC%ProposalStartCovMat                    = ProposalStartCovMat_type                  (nd,methodName)
        SpecMCMC%ProposalStartCorMat                    = ProposalStartCorMat_type                  (nd,methodName)
        SpecMCMC%ProposalStartStdVec                    = ProposalStartStdVec_type                  (nd,methodName)
        SpecMCMC%SampleRefinementCount                  = SampleRefinementCount_type                (methodName)
        SpecMCMC%SampleRefinementMethod                 = SampleRefinementMethod_type               (methodName)
        SpecMCMC%RandomStartPointRequested              = RandomStartPointRequested_type            (methodName)
        SpecMCMC%RandomStartPointDomainLowerLimitVec    = RandomStartPointDomainLowerLimitVec_type  (methodName)
        SpecMCMC%RandomStartPointDomainUpperLimitVec    = RandomStartPointDomainUpperLimitVec_type  (methodName)
    end function constructSpecMCMC

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar( SpecMCMC, nd )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        use Constants_mod, only: IK, RK
        implicit none
        class(SpecMCMC_type), intent(in)    :: SpecMCMC
        integer(IK), intent(in)             :: nd
        call SpecMCMC%ChainSize                             %nullifyNameListVar()
        call SpecMCMC%ScaleFactor                           %nullifyNameListVar()
        call SpecMCMC%StartPointVec                         %nullifyNameListVar(nd)
        call SpecMCMC%ProposalModel                         %nullifyNameListVar()
        call SpecMCMC%ProposalStartCovMat                   %nullifyNameListVar(nd)
        call SpecMCMC%ProposalStartCorMat                   %nullifyNameListVar(nd)
        call SpecMCMC%ProposalStartStdVec                   %nullifyNameListVar(nd)
        call SpecMCMC%SampleRefinementCount                 %nullifyNameListVar()
        call SpecMCMC%SampleRefinementMethod                %nullifyNameListVar()
        call SpecMCMC%RandomStartPointRequested             %nullifyNameListVar()
        call SpecMCMC%RandomStartPointDomainLowerLimitVec   %nullifyNameListVar(nd)
        call SpecMCMC%RandomStartPointDomainUpperLimitVec   %nullifyNameListVar(nd)
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setFromInputFile( SpecMCMC, Err, nd, domainLowerLimitVec, domainUpperLimitVec )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setFromInputFile
#endif

        use Constants_mod, only: IK, RK
        use Err_mod, only: Err_type
        implicit none
        class(SpecMCMC_type), intent(inout)     :: SpecMCMC
        integer(IK) , intent(in)                :: nd
        real(RK)    , intent(in)                :: domainLowerLimitVec(nd)
        real(RK)    , intent(in)                :: domainUpperLimitVec(nd)
        type(Err_type), intent(out)             :: Err

        Err%occurred = .false.
        Err%msg = ""

        call SpecMCMC%ChainSize                             %set(chainSize)
        call SpecMCMC%ScaleFactor                           %set(scaleFactor)
        call SpecMCMC%ProposalModel                         %set(trim(adjustl(proposalModel)))
        call SpecMCMC%ProposalStartCovMat                   %set(ProposalStartCovMat)
        call SpecMCMC%ProposalStartCorMat                   %set(ProposalStartCorMat)
        call SpecMCMC%ProposalStartStdVec                   %set(ProposalStartStdVec)
        call SpecMCMC%SampleRefinementCount                 %set(sampleRefinementCount)
        call SpecMCMC%SampleRefinementMethod                %set(sampleRefinementMethod)
        call SpecMCMC%RandomStartPointRequested             %set(randomStartPointRequested)
        call SpecMCMC%RandomStartPointDomainLowerLimitVec   %set(randomStartPointDomainLowerLimitVec, domainLowerLimitVec)
        call SpecMCMC%RandomStartPointDomainUpperLimitVec   %set(randomStartPointDomainUpperLimitVec, domainUpperLimitVec)
        call SpecMCMC%StartPointVec                         %set(startPointVec &
                                                                ,randomStartPointDomainLowerLimitVec    = SpecMCMC%RandomStartPointDomainLowerLimitVec%Val  &
                                                                ,randomStartPointDomainUpperLimitVec    = SpecMCMC%RandomStartPointDomainUpperLimitVec%Val  &
                                                                ,randomStartPointRequested              = SpecMCMC%RandomStartPointRequested%val            )

        deallocate(randomStartPointDomainLowerLimitVec)
        deallocate(randomStartPointDomainUpperLimitVec)
        deallocate(startPointVec)

    end subroutine setFromInputFile

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setFromInputArgs ( SpecMCMC, domainLowerLimitVec, domainUpperLimitVec &
                                ! ParaMCMC variables
                                , chainSize                             &
                                , scaleFactor                           &
                                , startPointVec                         &
                                , proposalModel                         &
                                , proposalStartCovMat                   &
                                , proposalStartCorMat                   &
                                , proposalStartStdVec                   &
                                , sampleRefinementCount                 &
                                , sampleRefinementMethod                &
                                , randomStartPointRequested             &
                                , randomStartPointDomainLowerLimitVec   &
                                , randomStartPointDomainUpperLimitVec   &
                                )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setFromInputArgs
#endif

        use Constants_mod, only: IK, RK
        implicit none
        class(SpecMCMC_type), intent(inout) :: SpecMCMC
        real(RK), intent(in)                :: domainLowerLimitVec(:), domainUpperLimitVec(:)

        ! ParaMCMC variables
        integer(IK) , intent(in), optional  :: chainSize
        character(*), intent(in), optional  :: scaleFactor
        real(RK)    , intent(in), optional  :: startPointVec(:)
        character(*), intent(in), optional  :: proposalModel
        real(RK)    , intent(in), optional  :: proposalStartCovMat(:,:)
        real(RK)    , intent(in), optional  :: proposalStartCorMat(:,:)
        real(RK)    , intent(in), optional  :: proposalStartStdVec(:)
        integer(IK) , intent(in), optional  :: sampleRefinementCount
        character(*), intent(in), optional  :: sampleRefinementMethod
        logical     , intent(in), optional  :: randomStartPointRequested
        real(RK)    , intent(in), optional  :: randomStartPointDomainLowerLimitVec(:)
        real(RK)    , intent(in), optional  :: randomStartPointDomainUpperLimitVec(:)

        if (present(chainSize))                             call SpecMCMC%ChainSize                             %set(chainSize)
        if (present(scaleFactor))                           call SpecMCMC%ScaleFactor                           %set(scaleFactor)
        if (present(proposalModel))                         call SpecMCMC%ProposalModel                         %set(trim(adjustl(proposalModel)))
        if (present(proposalStartCovMat))                   call SpecMCMC%ProposalStartCovMat                   %set(proposalStartCovMat)
        if (present(proposalStartCorMat))                   call SpecMCMC%ProposalStartCorMat                   %set(proposalStartCorMat)
        if (present(proposalStartStdVec))                   call SpecMCMC%ProposalStartStdVec                   %set(proposalStartStdVec)
        if (present(sampleRefinementCount))                 call SpecMCMC%SampleRefinementCount                 %set(sampleRefinementCount)
        if (present(sampleRefinementMethod))                call SpecMCMC%SampleRefinementMethod                %set(sampleRefinementMethod)
        if (present(randomStartPointRequested))             call SpecMCMC%RandomStartPointRequested             %set(randomStartPointRequested)
        if (present(randomStartPointDomainLowerLimitVec))   call SpecMCMC%RandomStartPointDomainLowerLimitVec   %set(randomStartPointDomainLowerLimitVec, domainLowerLimitVec)
        if (present(randomStartPointDomainUpperLimitVec))   call SpecMCMC%RandomStartPointDomainUpperLimitVec   %set(randomStartPointDomainUpperLimitVec, domainUpperLimitVec)
        if (present(startPointVec))                         call SpecMCMC%StartPointVec                         %set(startPointVec &
                                                                                                                ,randomStartPointDomainLowerLimitVec    = SpecMCMC%RandomStartPointDomainLowerLimitVec%Val  &
                                                                                                                ,randomStartPointDomainUpperLimitVec    = SpecMCMC%RandomStartPointDomainUpperLimitVec%Val  &
                                                                                                                ,randomStartPointRequested              = SpecMCMC%RandomStartPointRequested%val            )

    end subroutine setFromInputArgs

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine reportValues ( SpecMCMC              &
                            , prefix                &
                            , methodName            &
                            , outputUnit            &
                            , isMasterImage         &
                            , splashModeRequested   &
                            )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: reportValues
#endif
        use Decoration_mod, only: TAB
        use Constants_mod, only: IK, RK
        use Err_mod, only: note, informUser
        implicit none
        class(SpecMCMC_type), intent(in)    :: SpecMCMC
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
            write(outputUnit,formatStr) "chainSize"
            write(outputUnit,formatStr)
            write(outputUnit,formatVal) SpecMCMC%ChainSize%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%ChainSize%desc )


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "randomStartPointDomainLowerLimitVec"
            write(outputUnit,formatStr)
            do i = 1, size(SpecMCMC%RandomStartPointDomainLowerLimitVec%Val(:))
                write(outputUnit,formatVal) SpecMCMC%RandomStartPointDomainLowerLimitVec%Val(i)
            end do
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%RandomStartPointDomainLowerLimitVec%desc )


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "randomStartPointDomainUpperLimitVec"
            write(outputUnit,formatStr)
            do i = 1, size(SpecMCMC%RandomStartPointDomainUpperLimitVec%Val(:))
                write(outputUnit,formatVal) SpecMCMC%RandomStartPointDomainUpperLimitVec%Val(i)
            end do
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%RandomStartPointDomainUpperLimitVec%desc )


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "startPointVec"
            write(outputUnit,formatStr)
            do i = 1, size(SpecMCMC%StartPointVec%Val(:))
                write(outputUnit,formatVal) SpecMCMC%StartPointVec%Val(i)
            end do
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%StartPointVec%desc )


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "randomStartPointRequested"
            write(outputUnit,formatStr)
            write(outputUnit,formatVal) SpecMCMC%RandomStartPointRequested%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%RandomStartPointRequested%desc )


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "sampleRefinementCount"
            write(outputUnit,formatStr)
            write(outputUnit,formatVal) SpecMCMC%SampleRefinementCount%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%SampleRefinementCount%desc )


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "sampleRefinementMethod"
            write(outputUnit,formatStr)
            write(outputUnit,formatVal) SpecMCMC%SampleRefinementMethod%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%SampleRefinementMethod%desc )


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "scaleFactor"
            write(outputUnit,formatStr)
            write(outputUnit,formatVal) SpecMCMC%ScaleFactor%str
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%ScaleFactor%desc )


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
            write(outputUnit,formatVal) SpecMCMC%ProposalModel%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%ProposalModel%desc )


            ndim = size(SpecMCMC%ProposalStartCovMat%Val(:,1))
            allocate( Row(ndim) )


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "proposalStartCovMat"
            write(outputUnit,formatStr)
            if ( SpecMCMC%ProposalStartCovMat%isPresent ) then
                ! User has provided the Start Covariance Matrix
                do i = 1,ndim
                    Row = SpecMCMC%ProposalStartCovMat%Val(i,:)
                    write(outputUnit,formatVal) Row
                end do
            else
                ! User has not provided the Start Covariance Matrix
                call informUser ( prefix = TAB // TAB       &
                                , outputUnit = outputUnit   &
                                , newline = "\n"            &
                                , marginTop = 0_IK          &
                                , marginBot = 0_IK          &
                                , msg = "UNDEFINED. It will be constructed from the Correlation Matrix (ProposalStartCorMat) &
                                        &and the Standard Deviation vector (ProposalStartStdVec)." &
                                )
            end if
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%ProposalStartCovMat%desc )


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "proposalStartCorMat"
            write(outputUnit,formatStr)
            do i = 1,ndim
                Row = SpecMCMC%ProposalStartCorMat%Val(i,:)
                write(outputUnit,formatVal) Row
            end do
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%ProposalStartCorMat%desc )


            write(outputUnit,formatStr)
            write(outputUnit,formatStr) "proposalStartStdVec"
            write(outputUnit,formatStr)
            do i = 1,ndim
                write(outputUnit,formatVal) SpecMCMC%ProposalStartStdVec%Val(i)
            end do
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%ProposalStartStdVec%desc )

        end if


    end subroutine reportValues

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine checkForSanity(SpecMCMC,Err,methodName,nd,domainLowerLimitVec,domainUpperLimitVec)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Constants_mod, only: IK, RK
        use Err_mod, only: Err_type
        implicit none
        class(SpecMCMC_type), intent(inout) :: SpecMCMC
        integer(IK), intent(in)             :: nd
        real(RK), intent(in)                :: domainLowerLimitVec(:), domainUpperLimitVec(:)
        character(*), intent(in)            :: methodName
        type(Err_type), intent(inout)       :: Err
        call SpecMCMC%ChainSize                             %checkForSanity(Err,methodName,nd)
        call SpecMCMC%ScaleFactor                           %checkForSanity(Err,methodName)
        call SpecMCMC%ProposalModel                         %checkForSanity(Err,methodName)
        call SpecMCMC%ProposalStartCovMat                   %checkForSanity(Err,methodName,nd)
        call SpecMCMC%ProposalStartCorMat                   %checkForSanity(Err,methodName,nd)
        call SpecMCMC%ProposalStartStdVec                   %checkForSanity(Err,methodName,nd)
        call SpecMCMC%SampleRefinementCount                 %checkForSanity(Err,methodName)
        call SpecMCMC%SampleRefinementMethod                %checkForSanity(Err,methodName)
        call SpecMCMC%RandomStartPointDomainLowerLimitVec   %checkForSanity(Err,methodName,domainLowerLimitVec)
        call SpecMCMC%RandomStartPointDomainUpperLimitVec   %checkForSanity(Err,methodName,randomStartPointDomainLowerLimitVec=SpecMCMC%RandomStartPointDomainLowerLimitVec%Val,domainUpperLimitVec=domainUpperLimitVec)
        call SpecMCMC%StartPointVec                         %checkForSanity(Err,methodName,randomStartPointDomainLowerLimitVec=SpecMCMC%RandomStartPointDomainLowerLimitVec%Val,randomStartPointDomainUpperLimitVec=SpecMCMC%RandomStartPointDomainUpperLimitVec%Val)
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecMCMC_mod