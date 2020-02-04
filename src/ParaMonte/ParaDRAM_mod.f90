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

module ParaDRAM_mod

    use Constants_mod, only: IK, RK
    use ParaMonte_mod, only: ParaMonteNumFunCall_type
    use ParaMCMC_mod, only: ParaMCMC_type, ParaMCMC_Statistics_type, ParaMCMC_Chain_type
    use SpecDRAM_mod, only: SpecDRAM_type
    use ParaDRAMProposal_mod, only: Proposal_type
    use ParaMonteLogFunc_mod, only: getLogFunc_proc
    use ParaDRAMRefinedChain_mod, only: RefinedChain_type
    use ParaDRAMChainFileContents_mod, only: ChainFileContents_type

    implicit none

    character(*), parameter :: MODULE_NAME = "@ParaDRAM_mod"

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    type, extends(ParaMonteNumFunCall_type) :: ParaDRAM_NumFunCall_type
        integer(IK)                         :: acceptedRejectedDelayed          ! accepted + rejected function calls, including the delayed rejections
        integer(IK)                         :: acceptedRejectedDelayedUnused    ! by all processes, used or unused
    end type ParaDRAM_NumFunCall_type

    type, extends(ParaMCMC_Statistics_type) :: ParaDRAM_Statistics_type
        type(ParaMCMC_Chain_type)           :: AdaptationBurninLoc ! burning loc based on the minimum adaptation measure requested
        type(ParaDRAM_NumFunCall_type)      :: NumFunCall
    end type ParaDRAM_Statistics_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    type, extends(ParaMCMC_type)            :: ParaDRAM_type
        type(SpecDRAM_type)                 :: SpecDRAM
        type(ParaDRAM_Statistics_type)      :: Stats
        type(ChainFileContents_type)        :: Chain
        type(RefinedChain_type)             :: RefinedChain
       !class(Proposal_type), pointer       :: Proposal => null()
        class(Proposal_type), allocatable   :: Proposal
    contains
        procedure, pass, private            :: getSpecFromInputFile
        procedure, pass, public             :: runSampler
        procedure, pass, private            :: runKernel
    end type ParaDRAM_type

    !interface ParaDRAM_type
    !    module procedure :: runParaDRAM
    !end interface ParaDRAM_type

    interface
    module subroutine getSpecFromInputFile( PD, nd )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getSpecFromInputFile
#endif
        use Constants_mod, only: IK
        class(ParaDRAM_type), intent(inout) :: PD
        integer(IK), intent(in)             :: nd
    end subroutine getSpecFromInputFile
    end interface

    interface
    module subroutine runKernel( PD, getLogFunc )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: runKernel
#endif
        use ParaMonteLogFunc_mod, only: getLogFunc_proc
        class(ParaDRAM_type), intent(inout) :: PD
        procedure(getLogFunc_proc)          :: getLogFunc
    end subroutine runKernel
    end interface

    interface
    module subroutine runSampler( PD                                    &
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
                                , proposalStartCovMat                   &
                                , proposalStartCorMat                   &
                                , proposalStartStdVec                   &
                                , adaptiveUpdateCount                   &
                                , adaptiveUpdatePeriod                  &
                                , greedyAdaptationCount                 &
                                , delayedRejectionCount                 &
                                , burninAdaptationMeasure               &
                                , delayedRejectionScaleFactorVec        &
                                ) ! result(PD)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !!DEC$ ATTRIBUTES DLLEXPORT :: ParaDRAM_type
        !DEC$ ATTRIBUTES DLLEXPORT :: runSampler
#endif
        use ParaMonteLogFunc_mod, only: getLogFunc_proc
        use Constants_mod, only: IK, RK

        implicit none

        ! self
        class(ParaDRAM_type), intent(inout) :: PD

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
        integer(IK) , intent(in), optional  :: outputRealPrecision
        logical     , intent(in), optional  :: silentModeRequested
        real(RK)    , intent(in), optional  :: domainLowerLimitVec(ndim)
        real(RK)    , intent(in), optional  :: domainUpperLimitVec(ndim)
        character(*), intent(in), optional  :: parallelizationModel
        integer(IK) , intent(in), optional  :: progressReportPeriod
        real(RK)    , intent(in), optional  :: targetAcceptanceRate
        logical     , intent(in), optional  :: mpiFinalizeRequested
        integer(IK) , intent(in), optional  :: maxNumDomainCheckToWarn
        integer(IK) , intent(in), optional  :: maxNumDomainCheckToStop

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
        real(RK)    , intent(in), optional  :: proposalStartCovMat(ndim,ndim)
        real(RK)    , intent(in), optional  :: proposalStartCorMat(ndim,ndim)
        real(RK)    , intent(in), optional  :: proposalStartStdVec(ndim)
        integer(IK) , intent(in), optional  :: adaptiveUpdateCount
        integer(IK) , intent(in), optional  :: adaptiveUpdatePeriod
        integer(IK) , intent(in), optional  :: greedyAdaptationCount
        integer(IK) , intent(in), optional  :: delayedRejectionCount
        real(RK)    , intent(in), optional  :: burninAdaptationMeasure
        real(RK)    , intent(in), optional  :: delayedRejectionScaleFactorVec(:)

    end subroutine runSampler
    end interface

!***********************************************************************************************************************************
!***********************************************************************************************************************************

!contains
!
!!***********************************************************************************************************************************
!!***********************************************************************************************************************************
!
!#if defined CFI_ENABLED
!
!    ! The C/C++ interface of ParaDRAM
!
!    subroutine runParaDRAM  ( ndim          &
!                            , getLogFunc4C  &
!                            , InputFileVec  &
!                            , lenInputFile  &
!                            ) bind(C, name="runParaDRAM")
!#if defined DLL_ENABLED
!        !DEC$ ATTRIBUTES DLLEXPORT :: runParaDRAM
!#endif
!        use, intrinsic :: iso_c_binding, only: c_char, c_funptr, c_f_procpointer
!        use ParaMonteLogFunc_mod, only: getLogFunc_proc
!        use Constants_mod, only: IK, RK
!
!        implicit none
!
!        integer(IK), intent(in), value                          :: ndim
!        integer(IK), intent(in), value                          :: lenInputFile
!        type(c_funptr), intent(in), value                       :: getLogFunc4C
!        character(len=1,kind=c_char), dimension(*), intent(in)  :: InputFileVec
!        procedure(getLogFunc_proc), pointer                     :: getLogFunc
!        character(:), allocatable                               :: inputFileStr
!        type(ParaDRAM_type)                                     :: PD
!        integer                                                 :: i
!
!        ! reconstruct the input file
!
!        if (lenInputFile==0) then
!            inputFileStr = ""
!        else
!            allocate(character(lenInputFile) :: inputFileStr)
!            do i = 1, lenInputFile
!                inputFileStr(i:i) = InputFileVec(i)
!            end do
!        end if
!
!        ! associate the input C procedure pointer to a Fortran procedure pointer
!
!        call c_f_procpointer(cptr=getLogFunc4C, fptr=getLogFunc)
!
!        ! call runParaDRAM
!
!        call PD%runSampler  ( ndim = ndim               &
!                            , getLogFunc = getLogFunc   &
!                            , inputFile = inputFileStr  &
!                            )
!
!        nullify(getLogFunc)
!
!    end subroutine runParaDRAM
!
!#else
!
!    ! The simplistic Fortran interface of ParaDRAM
!
!    subroutine runParaDRAM  ( ndim          &
!                            , getLogFunc    &
!                            , inputFile     &
!                            )
!#if defined DLL_ENABLED
!        !DEC$ ATTRIBUTES DLLEXPORT :: runParaDRAM
!#endif
!        use ParaMonteLogFunc_mod, only: getLogFunc_proc
!        use Constants_mod, only: IK, RK
!
!        implicit none
!
!        integer(IK) , intent(in)    :: ndim
!        character(*), intent(in)    :: inputFile
!        procedure(getLogFunc_proc)  :: getLogFunc
!
!        type(ParaDRAM_type)         :: PD
!        integer                     :: i
!
!        ! call runParaDRAM
!
!        call PD%runSampler  ( ndim = ndim               &
!                            , getLogFunc = getLogFunc   &
!                            , inputFile = inputFile     &
!                            )
!
!    end subroutine runParaDRAM
!
!#endif
!
!!***********************************************************************************************************************************
!!***********************************************************************************************************************************

end module ParaDRAM_mod

!***********************************************************************************************************************************
!***********************************************************************************************************************************

#if defined CFI_ENABLED

    ! The C/C++ interface of ParaDRAM

    subroutine runParaDRAM  ( ndim          &
                            , getLogFunc4C  &
                            , InputFileVec  &
                            , lenInputFile  &
                            ) bind(C, name="runParaDRAM")
#if defined DLL_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: runParaDRAM
#endif
        use, intrinsic :: iso_c_binding, only: c_char, c_funptr, c_f_procpointer
        use ParaMonteLogFunc_mod, only: getLogFunc_proc
        use Constants_mod, only: IK, RK
        use ParaDRAM_mod, only: ParaDRAM_type

        implicit none

        integer(IK), intent(in), value                          :: ndim
        integer(IK), intent(in), value                          :: lenInputFile
        type(c_funptr), intent(in), value                       :: getLogFunc4C
        character(len=1,kind=c_char), dimension(*), intent(in)  :: InputFileVec
        procedure(getLogFunc_proc), pointer                     :: getLogFunc
        character(:), allocatable                               :: inputFileStr
        type(ParaDRAM_type)                                     :: PD
        integer                                                 :: i

        ! reconstruct the input file

        if (lenInputFile==0) then
            inputFileStr = ""
        else
            allocate(character(lenInputFile) :: inputFileStr)
            do i = 1, lenInputFile
                inputFileStr(i:i) = InputFileVec(i)
            end do
        end if

        ! associate the input C procedure pointer to a Fortran procedure pointer

        call c_f_procpointer(cptr=getLogFunc4C, fptr=getLogFunc)

        ! call runParaDRAM

        call PD%runSampler  ( ndim = ndim               &
                            , getLogFunc = getLogFunc   &
                            , inputFile = inputFileStr  &
                            )

        nullify(getLogFunc)

    end subroutine runParaDRAM

#else

    ! The simplistic Fortran interface of ParaDRAM

    subroutine runParaDRAM  ( ndim          &
                            , getLogFunc    &
                            , inputFile     &
                            )
#if defined DLL_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: runParaDRAM
#endif
        use ParaMonteLogFunc_mod, only: getLogFunc_proc
        use Constants_mod, only: IK, RK
        use ParaDRAM_mod, only: ParaDRAM_type

        implicit none

        integer(IK) , intent(in)    :: ndim
        character(*), intent(in)    :: inputFile
        procedure(getLogFunc_proc)  :: getLogFunc

        type(ParaDRAM_type)         :: PD
        integer                     :: i

        ! call runParaDRAM

        call PD%runSampler  ( ndim = ndim               &
                            , getLogFunc = getLogFunc   &
                            , inputFile = inputFile     &
                            )

    end subroutine runParaDRAM

#endif

!***********************************************************************************************************************************
!***********************************************************************************************************************************
