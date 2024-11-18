!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief
!>  This module contains derived types, procedures, and generic interfaces
!>  for the optional specifications of the ParaMonte library sampler routines.<br>
!>
!>  \details
!>  See the documentation of [getErrSampling](@ref pm_sampling::getErrSampling) for example usage.<br>
!>  You may wonder why there are multiple definitions of the derived types of this module in other sampling
!>  modules of the ParaMonte library (that have been intentionally keep hidden from the end users).<br>
!>  This organizational mess is due to the lack of a proper functional support for parameterized
!>  derived types in GNU Fortran compiler as Nov 2024.<br>
!>
!>  \note
!>  End users should normally rarely need any direct access to the contents of this module.<br>
!>  Instead, any application access should normally happen through
!>  the ParaMonte module [pm_sampling](@ref pm_sampling).<br>
!>
!>  \see
!>  [pm_sampling](@ref pm_sampling)<br>
!>
!>  \test
!>  [test_pm_sampling](@ref test_pm_sampling)<br>
!>
!>  \todo
!>  \pvhigh
!>  As soon as the PDT support in \gfortran is complete, the duplicate
!>  specification derived types across the library must be all merged and cleaned up.<br>
!>  To bypass the lack of a robust PDT support in \gfortran, this module stores all `real` components
!>  with the highest precision available which are internally converted to the appropriate-kind variables
!>  within the samplers.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Monday 00:01 AM, January 1, 2018, Institute for Computational Engineering and Sciences, University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_sampling_spec

    use pm_err, only: err_type
    use pm_str, only: getCharSeq
    use, intrinsic :: iso_c_binding, only: c_char, c_funptr, c_null_char, c_f_procpointer
    use pm_kind, only: SK, IK, LK, RKL, RKD, RKH, RKALL, SKG => SK ! LCOV_EXCL_LINE
    use pm_sysPath, only: PATHLEN => MAX_LEN_FILE_PATH

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_sampling_spec"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a base derived type for constructing objects containing
    !>  the basic **optional** simulation specifications of the ParaMonte library explorers and samplers.<br>
    !>
    !>  \details
    !>  See the documentation of [getErrSampling](@ref pm_sampling::getErrSampling) for example usage.<br>
    !>  Objects of this derived type are not meant to be used directly by the end users.<br>
    !>  Instead use one of the appropriate derived types:<br>
    !>  <ol>
    !>      <li>    [specdram_type](@ref pm_sampling_spec::specdram_type)
    !>      <li>    [specnest_type](@ref pm_sampling_spec::specnest_type)
    !>      <li>    [specdise_type](@ref pm_sampling_spec::specdise_type)
    !>  </ol>
    !>  This derived type has no custom constructor.<br>
    !>  All derived type components can be set after creating an object of the type.<br>
    !>
    !>  \note
    !>  For more information on any of the components of this derived type,
    !>  see the output `_report.txt` files from the corresponding sampler simulations.<br>
    !>  Alternatively, see the ParaMonte library cross-language sampler documentation
    !>  for detailed descriptions of the input simulation specifications at:<br>
    !>  <ol>
    !>      <li>    [ParaDRAM simulation options](\pmdoc_usage_sampling/paradram/specifications/)
    !>  </ol>
    !>
    !>  \note
    !>  This derived type has no custom constructor.<br>
    !>
    !>  \note
    !>  All derived type components can be set after creating an object of the type.<br>
    !>
    !>  \see
    !>  [specdram_type](@ref pm_sampling_spec::specdram_type)<br>
    !>  [specnest_type](@ref pm_sampling_spec::specnest_type)<br>
    !>  [specdise_type](@ref pm_sampling_spec::specdise_type)<br>
    !>  [getErrSampling](@ref pm_sampling::getErrSampling)<br>
    !>  [paradram_type](@ref pm_sampling::paradram_type)<br>
    !>  [paranest_type](@ref pm_sampling::paranest_type)<br>
    !>  [paradise_type](@ref pm_sampling::paradise_type)<br>
    !>  [sampler_type](@ref pm_sampling::sampler_type)<br>
    !>
    !>  \final{specbase_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type                                        :: specbase_type
        character(:,SKG)        , allocatable   :: description                          !<  \specdram{description}
        character(:,SKG)        , allocatable   :: domain                               !<  \specdram{domain}
        character(:,SKG)        , allocatable   :: domainAxisName(:)                    !<  \specdram{domainaxisname}
        real(RKH)               , allocatable   :: domainBallAvg(:)                     !<  \specdram{domainballavg}
        real(RKH)               , allocatable   :: domainBallCor(:,:)                   !<  \specdram{domainballcor}
        real(RKH)               , allocatable   :: domainBallCov(:,:)                   !<  \specdram{domainballcov}
        real(RKH)               , allocatable   :: domainBallStd(:)                     !<  \specdram{domainballstd}
        real(RKH)               , allocatable   :: domainCubeLimitLower(:)              !<  \specdram{domaincubelimitlower}
        real(RKH)               , allocatable   :: domainCubeLimitUpper(:)              !<  \specdram{domaincubelimitupper}
        integer(IK)             , allocatable   :: domainErrCount                       !<  \specdram{domainerrcount}
        integer(IK)             , allocatable   :: domainErrCountMax                    !<  \specdram{domainerrcountmax}
        character(:,SKG)        , allocatable   :: outputChainFileFormat                !<  \specdram{outputchainfileformat}
        integer(IK)             , allocatable   :: outputColumnWidth                    !<  \specdram{outputcolumnwidth}
        character(:,SKG)        , allocatable   :: outputFileName                       !<  \specdram{outputfilename}
        character(:,SKG)        , allocatable   :: outputStatus                         !<  \specdram{outputstatus}
        integer(IK)             , allocatable   :: outputPrecision                      !<  \specdram{outputprecision}
        integer(IK)             , allocatable   :: outputReportPeriod                   !<  \specdram{outputreportperiod}
        character(:,SKG)        , allocatable   :: outputRestartFileFormat              !<  \specdram{outputrestartfileformat}
        integer(IK)             , allocatable   :: outputSampleSize                     !<  \specdram{outputsamplesize}
        character(:,SKG)        , allocatable   :: outputSeparator                      !<  \specdram{outputseparator}
        character(:,SKG)        , allocatable   :: outputSplashMode                     !<  \specdram{outputsplashmode}
        character(:,SKG)        , allocatable   :: parallelism                          !<  \specdram{parallelism}
        logical(LK)             , allocatable   :: parallelismMpiFinalizeEnabled        !<  \specdram{parallelismmpifinalizeenabled}
        integer(IK)             , allocatable   :: parallelismNumThread                 !<  \specdram{parallelismnumthread}
        integer(IK)             , allocatable   :: randomSeed                           !<  \specdram{randomseed}
        character(:,SKG)        , allocatable   :: sysInfoFilePath                      !<  \specdram{sysinfofilepath}
        real(RKH)               , allocatable   :: targetAcceptanceRate(:)              !<  \specdram{targetacceptancerate}
    end type

    !>  \brief
    !>  This is a derived type for constructing objects containing
    !>  the **optional** simulation specifications of the ParaMonte library MCMC explorers and samplers.<br>
    !>
    !>  \details
    !>  See the documentation of [getErrSampling](@ref pm_sampling::getErrSampling) for example usage.<br>
    !>  Objects of this derived type are not meant to be used directly by the end users.<br>
    !>  Instead use one of the appropriate derived types:<br>
    !>  <ol>
    !>      <li>    [specdram_type](@ref pm_sampling_spec::specdram_type)
    !>      <li>    [specnest_type](@ref pm_sampling_spec::specnest_type)
    !>      <li>    [specdise_type](@ref pm_sampling_spec::specdise_type)
    !>  </ol>
    !>  This derived type has no custom constructor.<br>
    !>  All derived type components can be set after creating an object of the type.<br>
    !>
    !>  \note
    !>  For more information on any of the components of this derived type,
    !>  see the output `_report.txt` files from the corresponding sampler simulations.<br>
    !>  Alternatively, see the ParaMonte library cross-language sampler documentation
    !>  for detailed descriptions of the input simulation specifications at:<br>
    !>  <ol>
    !>      <li>    [ParaDRAM simulation options](\pmdoc_usage_sampling/paradram/specifications/)
    !>      <li>    [ParaDISE simulation options](\pmdoc_usage_sampling/paradise/specifications/)
    !>      <li>    [ParaNest simulation options](\pmdoc_usage_sampling/paranest/specifications/)
    !>  </ol>
    !>
    !>  \note
    !>  This derived type has no custom constructor.<br>
    !>
    !>  \note
    !>  All derived type components can be set after creating an object of the type.<br>
    !>
    !>  \see
    !>  [specdram_type](@ref pm_sampling_spec::specdram_type)<br>
    !>  [specnest_type](@ref pm_sampling_spec::specnest_type)<br>
    !>  [specdise_type](@ref pm_sampling_spec::specdise_type)<br>
    !>  [getErrSampling](@ref pm_sampling::getErrSampling)<br>
    !>  [paradram_type](@ref pm_sampling::paradram_type)<br>
    !>  [paranest_type](@ref pm_sampling::paranest_type)<br>
    !>  [paradise_type](@ref pm_sampling::paradise_type)<br>
    !>  [sampler_type](@ref pm_sampling::sampler_type)<br>
    !>
    !>  \final{specmcmc_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(specbase_type)                :: specmcmc_type
        integer(IK)             , allocatable   :: outputChainSize                      !<  \specdram{outputchainsize}
        integer(IK)             , allocatable   :: outputSampleRefinementCount          !<  \specdram{outputsamplerefinementcount}
        character(:,SKG)        , allocatable   :: outputSampleRefinementMethod         !<  \specdram{outputsamplerefinementmethod}
        character(:,SKG)        , allocatable   :: proposal                             !<  \specdram{proposal}
        real(RKH)               , allocatable   :: proposalCor(:,:)                     !<  \specdram{proposalcor}
        real(RKH)               , allocatable   :: proposalCov(:,:)                     !<  \specdram{proposalcov}
        character(:,SKG)        , allocatable   :: proposalScale                        !<  \specdram{proposalscale}
        real(RKH)               , allocatable   :: proposalStart(:)                     !<  \specdram{proposalstart}
        real(RKH)               , allocatable   :: proposalStartDomainCubeLimitLower(:) !<  \specdram{proposalstartdomaincubelimitlower}
        real(RKH)               , allocatable   :: proposalStartDomainCubeLimitUpper(:) !<  \specdram{proposalstartdomaincubelimitupper}
        logical(LK)             , allocatable   :: proposalStartRandomized              !<  \specdram{proposalstartrandomized}
        real(RKH)               , allocatable   :: proposalStd(:)                       !<  \specdram{proposalstd}
    end type

    !>  \brief
    !>  This is a derived type for constructing objects containing the **optional** simulation
    !>  properties of the ParaMonte library Delayed-Rejection Adaptive Metropolis (DRAM) MCMC explorers and samplers.<br>
    !>
    !>  \note
    !>  See the documentation of [getErrSampling](@ref pm_sampling::getErrSampling) for example usage.<br>
    !>  For more information on any of the components of this derived type,
    !>  see the output `_report.txt` files from the corresponding sampler simulations.<br>
    !>  Alternatively, see the ParaMonte library cross-language sampler documentation
    !>  for detailed descriptions of the input simulation specifications at:<br>
    !>  <ol>
    !>      <li>    [ParaDRAM simulation options](\pmdoc_usage_sampling/paradram/specifications/)
    !>      <li>    [ParaDISE simulation options](\pmdoc_usage_sampling/paradise/specifications/) (work in progress)
    !>      <li>    [ParaNest simulation options](\pmdoc_usage_sampling/paranest/specifications/) (work in progress)
    !>  </ol>
    !>
    !>  \note
    !>  This derived type has no custom constructor.<br>
    !>
    !>  \note
    !>  All derived type components can be set after creating an object of the type.<br>
    !>
    !>  \see
    !>  [specdram_type](@ref pm_sampling_spec::specdram_type)<br>
    !>  [specnest_type](@ref pm_sampling_spec::specnest_type)<br>
    !>  [specdise_type](@ref pm_sampling_spec::specdise_type)<br>
    !>  [getErrSampling](@ref pm_sampling::getErrSampling)<br>
    !>  [paradram_type](@ref pm_sampling::paradram_type)<br>
    !>  [paranest_type](@ref pm_sampling::paranest_type)<br>
    !>  [paradise_type](@ref pm_sampling::paradise_type)<br>
    !>  [sampler_type](@ref pm_sampling::sampler_type)<br>
    !>
    !>  \final{specdram_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(specmcmc_type)                :: specdram_type
        real(RKH)               , allocatable   :: proposalAdaptationBurnin             !<  \specdram{proposaladaptationburnin}
        integer(IK)             , allocatable   :: proposalAdaptationCount              !<  \specdram{proposaladaptationcount}
        integer(IK)             , allocatable   :: proposalAdaptationCountGreedy        !<  \specdram{proposaladaptationcountgreedy}
        integer(IK)             , allocatable   :: proposalAdaptationPeriod             !<  \specdram{proposaladaptationperiod}
        integer(IK)             , allocatable   :: proposalDelayedRejectionCount        !<  \specdram{proposaldelayedrejectioncount}
        real(RKH)               , allocatable   :: proposalDelayedRejectionScale(:)     !<  \specdram{proposaldelayedrejectionscale}
    end type

    !>  \brief
    !>  This is a derived type for constructing objects containing the **optional** simulation
    !>  properties of the ParaMonte library DRAM-enhanced MCMC explorers and samplers.<br>
    !>
    !>  \note
    !>  See the documentation of [getErrSampling](@ref pm_sampling::getErrSampling) for example usage.<br>
    !>  For more information on any of the components of this derived type,
    !>  see the output `_report.txt` files from the corresponding sampler simulations.<br>
    !>  Alternatively, see the ParaMonte library cross-language sampler documentation
    !>  for detailed descriptions of the input simulation specifications at:<br>
    !>  <ol>
    !>      <li>    [ParaDRAM simulation options](\pmdoc_usage_sampling/paradram/specifications/)
    !>      <li>    [ParaDISE simulation options](\pmdoc_usage_sampling/paradise/specifications/) (work in progress)
    !>      <li>    [ParaNest simulation options](\pmdoc_usage_sampling/paranest/specifications/) (work in progress)
    !>  </ol>
    !>
    !>  \note
    !>  This derived type has no custom constructor.<br>
    !>
    !>  \note
    !>  All derived type components can be set after creating an object of the type.<br>
    !>
    !>  \see
    !>  [specdram_type](@ref pm_sampling_spec::specdram_type)<br>
    !>  [specnest_type](@ref pm_sampling_spec::specnest_type)<br>
    !>  [specdise_type](@ref pm_sampling_spec::specdise_type)<br>
    !>  [getErrSampling](@ref pm_sampling::getErrSampling)<br>
    !>  [paradram_type](@ref pm_sampling::paradram_type)<br>
    !>  [paranest_type](@ref pm_sampling::paranest_type)<br>
    !>  [paradise_type](@ref pm_sampling::paradise_type)<br>
    !>  [sampler_type](@ref pm_sampling::sampler_type)<br>
    !>
    !>  \final{specdise_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(specdram_type)                :: specdise_type
    end type

    !>  \brief
    !>  This is a derived type for constructing objects containing the **optional** simulation
    !>  properties of the ParaMonte library Nested explorers and samplers.<br>
    !>
    !>  \note
    !>  See the documentation of [getErrSampling](@ref pm_sampling::getErrSampling) for example usage.<br>
    !>  For more information on any of the components of this derived type,
    !>  see the output `_report.txt` files from the corresponding sampler simulations.<br>
    !>  Alternatively, see the ParaMonte library cross-language sampler documentation
    !>  for detailed descriptions of the input simulation specifications at:<br>
    !>  <ol>
    !>      <li>    [ParaDRAM simulation options](\pmdoc_usage_sampling/paradram/specifications/)
    !>      <li>    [ParaDISE simulation options](\pmdoc_usage_sampling/paradise/specifications/) (work in progress)
    !>      <li>    [ParaNest simulation options](\pmdoc_usage_sampling/paranest/specifications/) (work in progress)
    !>  </ol>
    !>
    !>  \note
    !>  This derived type has no custom constructor.<br>
    !>
    !>  \note
    !>  All derived type components can be set after creating an object of the type.<br>
    !>
    !>  \see
    !>  [specdram_type](@ref pm_sampling_spec::specdram_type)<br>
    !>  [specnest_type](@ref pm_sampling_spec::specnest_type)<br>
    !>  [specdise_type](@ref pm_sampling_spec::specdise_type)<br>
    !>  [getErrSampling](@ref pm_sampling::getErrSampling)<br>
    !>  [paradram_type](@ref pm_sampling::paradram_type)<br>
    !>  [paranest_type](@ref pm_sampling::paranest_type)<br>
    !>  [paradise_type](@ref pm_sampling::paradise_type)<br>
    !>  [sampler_type](@ref pm_sampling::sampler_type)<br>
    !>
    !>  \final{paranest_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(specbase_type)                :: paranest_type
        integer(IK)             , allocatable   :: domainPartitionAdaptationCount
        integer(IK)             , allocatable   :: domainPartitionAdaptationPeriod
        logical(LK)             , allocatable   :: domainPartitionBiasCorrectionEnabled
        integer(IK)             , allocatable   :: domainPartitionCountMax
        real(RKH)               , allocatable   :: domainPartitionFactorExpansion
        real(RKH)               , allocatable   :: domainPartitionFactorShrinkage
        integer(IK)             , allocatable   :: domainPartitionKmeansClusterCountMax
        integer(IK)             , allocatable   :: domainPartitionKmeansClusterSizeMin
        logical(LK)             , allocatable   :: domainPartitionKmeansNormalizationEnabled
        integer(IK)             , allocatable   :: domainPartitionKmeansNumFailMax
        integer(IK)             , allocatable   :: domainPartitionKmeansNumRecursionMax
        integer(IK)             , allocatable   :: domainPartitionKmeansNumTry
        integer(IK)             , allocatable   :: domainPartitionKvolumeNumRecursionMax
        real(RKH)               , allocatable   :: domainPartitionKvolumeWeightExponent
        character(:,SKG)        , allocatable   :: domainPartitionMethod
        character(:,SKG)        , allocatable   :: domainPartitionObject
        logical(LK)             , allocatable   :: domainPartitionOptimizationScaleEnabled
        logical(LK)             , allocatable   :: domainPartitionOptimizationShapeEnabled
        logical(LK)             , allocatable   :: domainPartitionOptimizationShapeScaleEnabled
        real(RKH)               , allocatable   :: domainPartitionScaleFactor
        character(:,SKG)        , allocatable   :: domainSampler
        integer(IK)             , allocatable   :: liveSampleSize
        real(RKH)               , allocatable   :: tolerance
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_sampling_spec ! LCOV_EXCL_LINE