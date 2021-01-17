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
!> This file contains the body of `ParaDRAM_mod` and `ParaDISE_mod` modules.
!>
!> \remark
!> This module requires preprocessing, prior to compilation.
!>
!> \author Amir Shahmoradi

!#define CAT(a,b) a##b
!#define XCAT(a,b) CAT(a,b)
#if defined PARADRAM

#define ParaXXXX_type ParaDRAM_type
#define ParaXXXX_RefinedChain_mod ParaDRAM_RefinedChain_mod
#define ParaXXXX_ProposalAbstract_mod ParaDRAM_ProposalAbstract_mod
#define ParaXXXX_ChainFileContents_mod ParaDRAM_ChainFileContents_mod

#elif defined PARADISE

#define ParaXXXX_type ParaDISE_type
#define ParaXXXX_RefinedChain_mod ParaDISE_RefinedChain_mod
#define ParaXXXX_ProposalAbstract_mod ParaDISE_ProposalAbstract_mod
#define ParaXXXX_ChainFileContents_mod ParaDISE_ChainFileContents_mod

#elif defined PARANEST

#define ParaXXXX_type ParaNest_type
#define ParaXXXX_RefinedChain_mod ParaNest_RefinedChain_mod
#define ParaXXXX_ProposalAbstract_mod ParaNest_ProposalAbstract_mod
#define ParaXXXX_ChainFileContents_mod ParaNest_ChainFileContents_mod

#else
#error "Unrecognized sampler in ParaXXXX_mod.inc.f90"
#endif

    use Constants_mod, only: IK, RK, CK ! LCOV_EXCL_LINE
    use ParaMonte_mod, only: ParaMonteNumFunCall_type
    use ParaMonteLogFunc_mod, only: getLogFunc_proc
    use ParaXXXX_RefinedChain_mod, only: RefinedChain_type
    use ParaXXXX_ProposalAbstract_mod, only: ProposalAbstract_type
    use ParaXXXX_ChainFileContents_mod, only: ChainFileContents_type
#if defined PARADRAM || defined PARADISE
    use ParaMCMC_mod, only: ParaMCMC_Statistics_type
    use ParaMCMC_mod, only: ParaMCMC_Chain_type
    use ParaMCMC_mod, only: ParaMCMC_type
    use SpecDRAM_mod, only: SpecDRAM_type
#endif

    implicit none

#if defined PARADRAM
    character(*), parameter :: MODULE_NAME = "@ParaDRAM_mod"
#elif defined PARADISE
    character(*), parameter :: MODULE_NAME = "@ParaDISE_mod"
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type, extends(ParaMonteNumFunCall_type) :: NumFunCall_type
        integer(IK)                         :: acceptedRejectedDelayed          !< accepted + rejected function calls, including the delayed rejections
        integer(IK)                         :: acceptedRejectedDelayedUnused    !< by all processes, used or unused
    end type NumFunCall_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined PARADRAM || defined PARADISE
    type, extends(ParaMCMC_Statistics_type) :: Statistics_type
        type(ParaMCMC_Chain_type)           :: AdaptationBurninLoc              !< burning loc based on the minimum adaptation measure requested
        type(NumFunCall_type)               :: NumFunCall
    end type Statistics_type

    !> The `ParaDRAM_type` and `ParaDISE_type` classes.
    type, extends(ParaMCMC_type)                    :: ParaXXXX_type
        type(RefinedChain_type)                     :: RefinedChain             !< An object of class `RefinedChain_type` containing the final refined chain.
        type(SpecDRAM_type)                         :: SpecDRAM                 !< An object of class [SpecDRAM_type](@ref specdram_mod::specdram_type) containing DRAM simulation specs.
        type(Statistics_type)                       :: Stats                    !< An object of class [Statistics_type](@ref statistics_type) containing simulation statistics.
        type(ChainFileContents_type)                :: Chain                    !< An object of class `ChainFileContents_type containing information and methods for chain IO.
        class(ProposalAbstract_type), allocatable   :: Proposal                 !< An object of class `ProposalAbstract_type` representing the proposal distribution of the MCMC sampler.
       !class(ProposalAbstract_type), pointer       :: Proposal => null()
    contains
        procedure, pass, private                    :: getSpecFromInputFile
        procedure, pass, public                     :: runSampler
        procedure, pass, private                    :: runKernel
        procedure, pass, private                    :: postprocess
    end type ParaXXXX_type ! ParaDRAM_type or ParaDISE_type
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
    module subroutine getSpecFromInputFile( self, nd )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSpecFromInputFile
#endif
        use Constants_mod, only: IK
        class(ParaXXXX_type), intent(inout) :: self
        integer(IK), intent(in)             :: nd
    end subroutine getSpecFromInputFile
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
    module subroutine runKernel( self, getLogFunc )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: runKernel
#endif
        use ParaMonteLogFunc_mod, only: getLogFunc_proc
        class(ParaXXXX_type), intent(inout) :: self
        procedure(getLogFunc_proc)          :: getLogFunc
    end subroutine runKernel
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
                                , overwriteRequested                    &
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
#if defined PARADRAM || defined PARADISE
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
#elif defined ParaNest
                                ! ParaNest variables
                                , tightness                             &
                                , tolerance                             &
                                , scaleFactor                           &
                                , proposalModel                         &
                                , liveSampleSize                        &
                                , inclusionFraction                     &
                                , adaptiveUpdateCount                   &
                                , adaptiveUpdatePeriod                  &
                                , mahalSqWeightExponent                 &
                                , stabilizationRequested                &
                                , MaxAllowedKmeansFailure               &
                                , maxAllowedMinVolFailure               &
                                , maxKvolumeLoopRecursion               &
#endif
                                )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: runSampler
#endif
        use ParaMonteLogFunc_mod, only: getLogFunc_proc
        use Constants_mod, only: IK, RK

        implicit none

        ! self
        class(ParaXXXX_type), intent(inout) :: self

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
        logical     , intent(in), optional  :: overwriteRequested
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
#if defined PARADRAM || defined PARADISE
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
#elif defined ParaNest
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
#endif

    end subroutine runSampler
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
    module subroutine postprocess(self)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: postprocess
#endif
        class(ParaXXXX_type), intent(inout) :: self
    end subroutine postprocess
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined PARADISE
#elif defined PARADRAM
#elif defined PARANEST
#endif

#undef ParaXXXX_type
#undef ParaXXXX_RefinedChain_mod
#undef ParaXXXX_ProposalAbstract_mod
#undef ParaXXXX_ChainFileContents_mod

