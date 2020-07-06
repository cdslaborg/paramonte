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
!   we ask you to acknowledge the use of the ParaMonte library
!   in your work (education/research/industry/development/...)
!   by citing the ParaMonte library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!***********************************************************************************************************************************
!***********************************************************************************************************************************

!#define CAT(a,b) a##b
!#define XCAT(a,b) CAT(a,b)

#if defined ParaDISE
!module ParaDISE_mod
#elif defined ParaDRAM
module ParaDRAM_mod
#endif

    use Constants_mod, only: IK, RK, CK
    use ParaMonte_mod, only: ParaMonteNumFunCall_type
    use ParaMCMC_mod, only: ParaMCMC_type, ParaMCMC_Statistics_type, ParaMCMC_Chain_type
    use SpecDRAM_mod, only: SpecDRAM_type
    use ParaMonteLogFunc_mod, only: getLogFunc_proc
    use ParaDRAMProposalTemplate_mod, only: ProposalTemplate_type

    implicit none

    character(*), parameter :: MODULE_NAME = "@ParaDRAM_mod"

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    type, extends(ParaMonteNumFunCall_type) :: ParaDRAM_NumFunCall_type
        integer(IK)                         :: acceptedRejectedDelayed          ! accepted + rejected function calls, including the delayed rejections
        integer(IK)                         :: acceptedRejectedDelayedUnused    ! by all processes, used or unused
    end type ParaDRAM_NumFunCall_type

    type, extends(ParaMCMC_Statistics_type) :: ParaDRAM_Statistics_type
        type(ParaMCMC_Chain_type)           :: AdaptationBurninLoc              ! burning loc based on the minimum adaptation measure requested
        type(ParaDRAM_NumFunCall_type)      :: NumFunCall
    end type ParaDRAM_Statistics_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    type, extends(ParaMCMC_type)                    :: ParaDRAM_type
        type(SpecDRAM_type)                         :: SpecDRAM
        type(ParaDRAM_Statistics_type)              :: Stats
       !type(RefinedChain_type)                     :: RefinedChain
       !class(ProposalTemplate_type), pointer       :: Proposal => null()
        class(ProposalTemplate_type), allocatable   :: Proposal
    contains
        procedure, pass, private                    :: getSpecFromInputFile
        procedure, pass, public                     :: runSampler
        procedure, pass, private                    :: runKernel
    end type ParaDRAM_type

    !interface ParaDRAM_type
    !    module procedure :: runParaDRAM
    !end interface ParaDRAM_type

    interface
    module subroutine getSpecFromInputFile( self, nd )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getSpecFromInputFile
#endif
        use Constants_mod, only: IK
        class(ParaDRAM_type), intent(inout) :: self
        integer(IK), intent(in)             :: nd
    end subroutine getSpecFromInputFile
    end interface

    interface
    module subroutine runKernel( self, getLogFunc )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: runKernel
#endif
        use ParaMonteLogFunc_mod, only: getLogFunc_proc
        class(ParaDRAM_type), intent(inout) :: self
        procedure(getLogFunc_proc)          :: getLogFunc
    end subroutine runKernel
    end interface

    interface
    module subroutine runSampler( self                                  &
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
                                , TargetAcceptanceRate                  &
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
                                , proposalStartStdVec                   &
                                , proposalStartCorMat                   &
                                , proposalStartCovMat                   &
                                , adaptiveUpdateCount                   &
                                , adaptiveUpdatePeriod                  &
                                , greedyAdaptationCount                 &
                                , delayedRejectionCount                 &
                                , burninAdaptationMeasure               &
                                , delayedRejectionScaleFactorVec        &
                                ) ! result(self)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !!DEC$ ATTRIBUTES DLLEXPORT :: ParaDRAM_type
        !DEC$ ATTRIBUTES DLLEXPORT :: runSampler
#endif
        use ParaMonteLogFunc_mod, only: getLogFunc_proc
        use Constants_mod, only: IK, RK

        implicit none

        ! self
        class(ParaDRAM_type), intent(inout) :: self

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
        real(RK)    , intent(in), optional  :: TargetAcceptanceRate(2)
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
        real(RK)    , intent(in), optional  :: proposalStartStdVec(ndim)
        real(RK)    , intent(in), optional  :: proposalStartCorMat(ndim,ndim)
        real(RK)    , intent(in), optional  :: proposalStartCovMat(ndim,ndim)
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

#if defined ParaDISE
!end module ParaDISE_mod
#elif defined ParaDRAM
end module ParaDRAM_mod
#endif
