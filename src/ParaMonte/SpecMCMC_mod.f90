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

module SpecMCMC_mod

    ! ParaMCMC Spec variable types
    use SpecMCMC_ChainSize_mod                          , only: ChainSize_type
    use SpecMCMC_StartPointVec_mod                      , only: StartPointVec_type
    use SpecMCMC_SampleRefinementCount_mod              , only: SampleRefinementCount_type
    use SpecMCMC_SampleRefinementMethod_mod             , only: SampleRefinementMethod_type
    use SpecMCMC_RandomStartPointRequested_mod          , only: RandomStartPointRequested_type
    use SpecMCMC_RandomStartPointDomainLowerLimitVec_mod, only: RandomStartPointDomainLowerLimitVec_type
    use SpecMCMC_RandomStartPointDomainUpperLimitVec_mod, only: RandomStartPointDomainUpperLimitVec_type

    ! ParaMCMC namelist variables
    use SpecMCMC_ChainSize_mod                          , only: ChainSize
    use SpecMCMC_StartPointVec_mod                      , only: startPointVec
    use SpecMCMC_SampleRefinementCount_mod              , only: SampleRefinementCount
    use SpecMCMC_SampleRefinementMethod_mod             , only: sampleRefinementMethod
    use SpecMCMC_RandomStartPointRequested_mod          , only: randomStartPointRequested
    use SpecMCMC_RandomStartPointDomainLowerLimitVec_mod, only: randomStartPointDomainLowerLimitVec
    use SpecMCMC_RandomStartPointDomainUpperLimitVec_mod, only: randomStartPointDomainUpperLimitVec

    implicit none

    type                                                :: SpecMCMC_type
        type(ChainSize_type)                            :: ChainSize
        type(StartPointVec_type)                        :: StartPointVec
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

    function constructSpecMCMC(methodName) result(SpecMCMC)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructSpecMCMC
#endif
        use Constants_mod, only: IK,RK
        implicit none
        character(*), intent(in)    :: methodName
        type(SpecMCMC_type)         :: SpecMCMC
        SpecMCMC%ChainSize                              = ChainSize_type(methodName)
        SpecMCMC%StartPointVec                          = StartPointVec_type()
        SpecMCMC%SampleRefinementCount                  = SampleRefinementCount_type(methodName)
        SpecMCMC%SampleRefinementMethod                 = SampleRefinementMethod_type(methodName)
        SpecMCMC%RandomStartPointRequested              = RandomStartPointRequested_type(methodName)
        SpecMCMC%RandomStartPointDomainLowerLimitVec    = RandomStartPointDomainLowerLimitVec_type(methodName)
        SpecMCMC%RandomStartPointDomainUpperLimitVec    = RandomStartPointDomainUpperLimitVec_type(methodName)
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
        call SpecMCMC%StartPointVec                         %nullifyNameListVar(nd)
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
                                , chainSize &
                                , startPointVec &
                                , sampleRefinementCount &
                                , sampleRefinementMethod &
                                , randomStartPointRequested &
                                , randomStartPointDomainLowerLimitVec &
                                , randomStartPointDomainUpperLimitVec &
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
        real(RK)    , intent(in), optional  :: startPointVec(:)
        integer(IK) , intent(in), optional  :: sampleRefinementCount
        character(*), intent(in), optional  :: sampleRefinementMethod
        logical     , intent(in), optional  :: randomStartPointRequested
        real(RK)    , intent(in), optional  :: randomStartPointDomainLowerLimitVec(:)
        real(RK)    , intent(in), optional  :: randomStartPointDomainUpperLimitVec(:)

        if (present(chainSize))                             call SpecMCMC%ChainSize                             %set(chainSize)
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

    subroutine reportValues(SpecMCMC,prefix,outputUnit,isMasterImage,splashModeRequested)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: reportValues
#endif
        use Decoration_mod, only: TAB
        use Constants_mod, only: IK
        use Err_mod, only: note
        implicit none
        class(SpecMCMC_type), intent(in)    :: SpecMCMC
        character(*), intent(in)            :: prefix
        integer(IK), intent(in)             :: outputUnit
        logical, intent(in)                 :: isMasterImage, splashModeRequested
        integer(IK)                         :: i
        character(:), allocatable           :: formatStr, formatVal

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
        class(SpecMCMC_type), intent(in)    :: SpecMCMC
        integer(IK), intent(in)             :: nd
        real(RK), intent(in)                :: domainLowerLimitVec(:), domainUpperLimitVec(:)
        character(*), intent(in)            :: methodName
        type(Err_type), intent(inout)       :: Err
        call SpecMCMC%ChainSize                             %checkForSanity(Err,methodName,nd)
        call SpecMCMC%SampleRefinementCount                 %checkForSanity(Err,methodName)
        call SpecMCMC%SampleRefinementMethod                %checkForSanity(Err,methodName)
        call SpecMCMC%RandomStartPointDomainLowerLimitVec   %checkForSanity(Err,methodName,domainLowerLimitVec)
        call SpecMCMC%RandomStartPointDomainUpperLimitVec   %checkForSanity(Err,methodName,randomStartPointDomainLowerLimitVec=SpecMCMC%RandomStartPointDomainLowerLimitVec%Val,domainUpperLimitVec=domainUpperLimitVec)
        call SpecMCMC%StartPointVec                         %checkForSanity(Err,methodName,randomStartPointDomainLowerLimitVec=SpecMCMC%RandomStartPointDomainLowerLimitVec%Val,randomStartPointDomainUpperLimitVec=SpecMCMC%RandomStartPointDomainUpperLimitVec%Val)
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecMCMC_mod