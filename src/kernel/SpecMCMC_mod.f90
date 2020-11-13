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

module SpecMCMC_mod

    ! ParaMCMC Spec variable types
    use SpecMCMC_ChainSize_mod                          , only: ChainSize_type
    use SpecMCMC_ScaleFactor_mod                        , only: ScaleFactor_type
    use SpecMCMC_StartPointVec_mod                      , only: StartPointVec_type
    use SpecMCMC_ProposalModel_mod                      , only: ProposalModel_type
    use SpecMCMC_ProposalStartStdVec_mod                , only: ProposalStartStdVec_type
    use SpecMCMC_ProposalStartCorMat_mod                , only: ProposalStartCorMat_type
    use SpecMCMC_ProposalStartCovMat_mod                , only: ProposalStartCovMat_type
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
    use SpecMCMC_ProposalStartStdVec_mod                , only: proposalStartStdVec
    use SpecMCMC_ProposalStartCorMat_mod                , only: proposalStartCorMat
    use SpecMCMC_ProposalStartCovMat_mod                , only: proposalStartCovMat
    use SpecMCMC_SampleRefinementCount_mod              , only: SampleRefinementCount
    use SpecMCMC_SampleRefinementMethod_mod             , only: sampleRefinementMethod
    use SpecMCMC_RandomStartPointRequested_mod          , only: randomStartPointRequested
    use SpecMCMC_RandomStartPointDomainLowerLimitVec_mod, only: randomStartPointDomainLowerLimitVec
    use SpecMCMC_RandomStartPointDomainUpperLimitVec_mod, only: randomStartPointDomainUpperLimitVec

    implicit none

    type                                                :: SpecMCMC_type
        type(ChainSize_type)                            :: ChainSize
        type(ScaleFactor_type)                          :: ScaleFactor
        type(StartPointVec_type)                        :: startPointVec
        type(ProposalModel_type)                        :: ProposalModel
        type(ProposalStartStdVec_type)                  :: proposalStartStdVec
        type(ProposalStartCorMat_type)                  :: proposalStartCorMat
        type(ProposalStartCovMat_type)                  :: proposalStartCovMat
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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        SpecMCMC%startPointVec                          = StartPointVec_type                        ()
        SpecMCMC%ProposalModel                          = ProposalModel_type                        ()
        SpecMCMC%proposalStartStdVec                    = ProposalStartStdVec_type                  (nd,methodName)
        SpecMCMC%proposalStartCorMat                    = ProposalStartCorMat_type                  (nd,methodName)
        SpecMCMC%proposalStartCovMat                    = ProposalStartCovMat_type                  (nd,methodName)
        SpecMCMC%SampleRefinementCount                  = SampleRefinementCount_type                (methodName)
        SpecMCMC%SampleRefinementMethod                 = SampleRefinementMethod_type               (methodName)
        SpecMCMC%RandomStartPointRequested              = RandomStartPointRequested_type            (methodName)
        SpecMCMC%RandomStartPointDomainLowerLimitVec    = RandomStartPointDomainLowerLimitVec_type  (methodName)
        SpecMCMC%RandomStartPointDomainUpperLimitVec    = RandomStartPointDomainUpperLimitVec_type  (methodName)
    end function constructSpecMCMC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        call SpecMCMC%startPointVec                         %nullifyNameListVar(nd)
        call SpecMCMC%ProposalModel                         %nullifyNameListVar()
        call SpecMCMC%proposalStartStdVec                   %nullifyNameListVar(nd)
        call SpecMCMC%proposalStartCorMat                   %nullifyNameListVar(nd)
        call SpecMCMC%proposalStartCovMat                   %nullifyNameListVar(nd)
        call SpecMCMC%SampleRefinementCount                 %nullifyNameListVar()
        call SpecMCMC%SampleRefinementMethod                %nullifyNameListVar()
        call SpecMCMC%RandomStartPointRequested             %nullifyNameListVar()
        call SpecMCMC%RandomStartPointDomainLowerLimitVec   %nullifyNameListVar(nd)
        call SpecMCMC%RandomStartPointDomainUpperLimitVec   %nullifyNameListVar(nd)
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        call SpecMCMC%proposalStartStdVec                   %set(proposalStartStdVec)
        call SpecMCMC%proposalStartCorMat                   %set(proposalStartCorMat)
        call SpecMCMC%proposalStartCovMat                   %set(SpecMCMC%proposalStartStdVec%val, SpecMCMC%proposalStartCorMat%val, proposalStartCovMat)
        call SpecMCMC%SampleRefinementCount                 %set(sampleRefinementCount)
        call SpecMCMC%SampleRefinementMethod                %set(sampleRefinementMethod)
        call SpecMCMC%RandomStartPointRequested             %set(randomStartPointRequested)
        call SpecMCMC%RandomStartPointDomainLowerLimitVec   %set(randomStartPointDomainLowerLimitVec, domainLowerLimitVec)
        call SpecMCMC%RandomStartPointDomainUpperLimitVec   %set(randomStartPointDomainUpperLimitVec, domainUpperLimitVec)
        call SpecMCMC%startPointVec                         %set(startPointVec &
                                                                ,randomStartPointDomainLowerLimitVec    = SpecMCMC%RandomStartPointDomainLowerLimitVec%val  &
                                                                ,randomStartPointDomainUpperLimitVec    = SpecMCMC%RandomStartPointDomainUpperLimitVec%val  &
                                                                ,randomStartPointRequested              = SpecMCMC%RandomStartPointRequested%val            &
                                                                ,domainLowerLimitVec                    = domainLowerLimitVec                               &
                                                                ,domainUpperLimitVec                    = domainUpperLimitVec                               )


        deallocate(randomStartPointDomainLowerLimitVec)
        deallocate(randomStartPointDomainUpperLimitVec)
        deallocate(startPointVec)

    end subroutine setFromInputFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setFromInputArgs ( SpecMCMC, domainLowerLimitVec, domainUpperLimitVec &
                                ! ParaMCMC variables
                                , chainSize                             &
                                , scaleFactor                           &
                                , startPointVec                         &
                                , proposalModel                         &
                                , proposalStartStdVec                   &
                                , proposalStartCorMat                   &
                                , proposalStartCovMat                   &
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
        real(RK)    , intent(in), optional  :: proposalStartStdVec(:)
        real(RK)    , intent(in), optional  :: proposalStartCorMat(:,:)
        real(RK)    , intent(in), optional  :: proposalStartCovMat(:,:)
        integer(IK) , intent(in), optional  :: sampleRefinementCount
        character(*), intent(in), optional  :: sampleRefinementMethod
        logical     , intent(in), optional  :: randomStartPointRequested
        real(RK)    , intent(in), optional  :: randomStartPointDomainLowerLimitVec(:)
        real(RK)    , intent(in), optional  :: randomStartPointDomainUpperLimitVec(:)

        logical                             :: proposalStartStdVecIsPresent
        logical                             :: proposalStartCorMatIsPresent
        logical                             :: proposalStartCovMatIsPresent

        proposalStartStdVecIsPresent = present(proposalStartStdVec)
        proposalStartCorMatIsPresent = present(proposalStartCorMat)
        proposalStartCovMatIsPresent = present(proposalStartCovMat) .or. proposalStartCorMatIsPresent .or. proposalStartStdVecIsPresent

        if (present(chainSize))                             call SpecMCMC%ChainSize                             %set(chainSize)
        if (present(scaleFactor))                           call SpecMCMC%ScaleFactor                           %set(scaleFactor)
        if (present(proposalModel))                         call SpecMCMC%ProposalModel                         %set(trim(adjustl(proposalModel)))
        if (proposalStartStdVecIsPresent)                   call SpecMCMC%proposalStartStdVec                   %set(proposalStartStdVec)
        if (proposalStartCorMatIsPresent)                   call SpecMCMC%proposalStartCorMat                   %set(proposalStartCorMat)
        if (proposalStartCovMatIsPresent)                   call SpecMCMC%proposalStartCovMat                   %set(SpecMCMC%proposalStartStdVec%val, SpecMCMC%proposalStartCorMat%val, proposalStartCovMat)
        if (present(sampleRefinementCount))                 call SpecMCMC%SampleRefinementCount                 %set(sampleRefinementCount)
        if (present(sampleRefinementMethod))                call SpecMCMC%SampleRefinementMethod                %set(sampleRefinementMethod)
        if (present(randomStartPointRequested))             call SpecMCMC%RandomStartPointRequested             %set(randomStartPointRequested)
        if (present(randomStartPointDomainLowerLimitVec))   call SpecMCMC%RandomStartPointDomainLowerLimitVec   %set(randomStartPointDomainLowerLimitVec, domainLowerLimitVec)
        if (present(randomStartPointDomainUpperLimitVec))   call SpecMCMC%RandomStartPointDomainUpperLimitVec   %set(randomStartPointDomainUpperLimitVec, domainUpperLimitVec)
        if (present(startPointVec))                         call SpecMCMC%startPointVec                         %set(startPointVec &
                                                                                                                ,randomStartPointDomainLowerLimitVec    = SpecMCMC%RandomStartPointDomainLowerLimitVec%val  &
                                                                                                                ,randomStartPointDomainUpperLimitVec    = SpecMCMC%RandomStartPointDomainUpperLimitVec%val  &
                                                                                                                ,randomStartPointRequested              = SpecMCMC%RandomStartPointRequested%val            &
                                                                                                                ,domainLowerLimitVec                    = domainLowerLimitVec                               &
                                                                                                                ,domainUpperLimitVec                    = domainUpperLimitVec                               )

    end subroutine setFromInputArgs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        use Decoration_mod, only: GENERIC_OUTPUT_FORMAT
        use Decoration_mod, only: GENERIC_TABBED_FORMAT
        use Decoration_mod, only: TAB
        use Constants_mod, only: IK, RK
        use Err_mod, only: note
        implicit none
        class(SpecMCMC_type), intent(in)    :: SpecMCMC
        character(*), intent(in)            :: prefix, methodName
        integer(IK), intent(in)             :: outputUnit
        logical, intent(in)                 :: isMasterImage, splashModeRequested
        integer(IK)                         :: ndim, i

        if (isMasterImage) then

            ndim = size(SpecMCMC%proposalStartCovMat%val(:,1))

            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "chainSize"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecMCMC%ChainSize%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%ChainSize%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "randomStartPointDomainLowerLimitVec"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            do i = 1, size(SpecMCMC%RandomStartPointDomainLowerLimitVec%val(:))
                write(outputUnit,GENERIC_TABBED_FORMAT) SpecMCMC%RandomStartPointDomainLowerLimitVec%val(i)
            end do
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%RandomStartPointDomainLowerLimitVec%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "randomStartPointDomainUpperLimitVec"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            do i = 1, size(SpecMCMC%RandomStartPointDomainUpperLimitVec%val(:))
                write(outputUnit,GENERIC_TABBED_FORMAT) SpecMCMC%RandomStartPointDomainUpperLimitVec%val(i)
            end do
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%RandomStartPointDomainUpperLimitVec%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "startPointVec"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            do i = 1, size(SpecMCMC%startPointVec%val(:))
                write(outputUnit,GENERIC_TABBED_FORMAT) SpecMCMC%startPointVec%val(i)
            end do
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%startPointVec%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "randomStartPointRequested"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecMCMC%RandomStartPointRequested%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%RandomStartPointRequested%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "sampleRefinementCount"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecMCMC%SampleRefinementCount%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%SampleRefinementCount%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "sampleRefinementMethod"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecMCMC%SampleRefinementMethod%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%SampleRefinementMethod%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "scaleFactor"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecMCMC%ScaleFactor%str
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%ScaleFactor%desc )


            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! proposal distribution
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            !block
            !    use Decoration_mod, only: writeDecoratedText
            !    call writeDecoratedText ( text = "\n" // methodName // " proposal specifications\n" &
            !                            , marginTop = 1     &
            !                            , marginBot = 1     &
            !                            , newline = "\n"    &
            !                            , outputUnit = outputUnit )
            !end block


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "proposalModel"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_TABBED_FORMAT) SpecMCMC%ProposalModel%val
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%ProposalModel%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "proposalStartStdVec"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            do i = 1, ndim
                write(outputUnit,GENERIC_TABBED_FORMAT) SpecMCMC%proposalStartStdVec%val(i)
            end do
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%proposalStartStdVec%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "proposalStartCorMat"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            do i = 1, ndim
                write(outputUnit,GENERIC_TABBED_FORMAT) SpecMCMC%proposalStartCorMat%val(:,i)
            end do
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%proposalStartCorMat%desc )


            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            write(outputUnit,GENERIC_OUTPUT_FORMAT) "proposalStartCovMat"
            write(outputUnit,GENERIC_OUTPUT_FORMAT)
            do i = 1, ndim
                write(outputUnit,GENERIC_TABBED_FORMAT) SpecMCMC%proposalStartCovMat%val(:,i)
            end do
            if (splashModeRequested) call note( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = SpecMCMC%proposalStartCovMat%desc )

        end if


    end subroutine reportValues

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        call SpecMCMC%proposalStartCovMat                   %checkForSanity(Err,methodName,nd)
        call SpecMCMC%proposalStartCorMat                   %checkForSanity(Err,methodName,nd)
        call SpecMCMC%proposalStartStdVec                   %checkForSanity(Err,methodName,nd)
        call SpecMCMC%SampleRefinementCount                 %checkForSanity(Err,methodName)
        call SpecMCMC%SampleRefinementMethod                %checkForSanity(Err,methodName)
        call SpecMCMC%RandomStartPointDomainLowerLimitVec   %checkForSanity(Err,methodName,domainLowerLimitVec)
        call SpecMCMC%RandomStartPointDomainUpperLimitVec   %checkForSanity(Err,methodName,randomStartPointDomainLowerLimitVec=SpecMCMC%RandomStartPointDomainLowerLimitVec%val,domainUpperLimitVec=domainUpperLimitVec)
        call SpecMCMC%startPointVec                         %checkForSanity(Err,methodName &
                                                                           ,domainLowerLimitVec = domainLowerLimitVec &
                                                                           ,domainUpperLimitVec = domainUpperLimitVec )

    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecMCMC_mod