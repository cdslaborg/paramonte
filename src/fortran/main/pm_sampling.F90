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
!>  This module contains procedures and generic interfaces for the ParaMonte library sampler routines.
!>
!>  \test
!>  [test_pm_sampling](@ref test_pm_sampling)<br>
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Monday 00:01 AM, January 1, 2018, Institute for Computational Engineering and Sciences, University of Texas Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_sampling

    use pm_err, only: err_type
    use pm_str, only: getCharSeq
    use, intrinsic :: iso_c_binding, only: c_char, c_funptr, c_null_char, c_f_procpointer
    use pm_kind, only: SK, IK, LK, RKL, RKD, RKH, RKALL, SKC => SK ! LCOV_EXCL_LINE
    use pm_sysPath, only: PATHLEN => MAX_LEN_FILE_PATH

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_sampling"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a derived type for constructing objects containing
    !>  the **optional** simulation properties of the ParaMonte library explorers and samplers.<br>
    !>
    !>  \details
    !>  Objects of this derived type are not meant to be used directly by the end users.<br>
    !>  Instead use one of the appropriate derived types:<br>
    !>  <ol>
    !>      <li>    [paradram_type](@ref pm_sampling::paradram_type)
    !>      <li>    [paranest_type](@ref pm_sampling::paranest_type)
    !>      <li>    [paradise_type](@ref pm_sampling::paradise_type)
    !>  </ol>
    !>
    !>  \note
    !>  For more information on any of the components of this derived type,
    !>  see the output `_report.txt` files from the corresponding sampler simulations.<br>
    !>  Alternatively, see the ParaMonte library cross-language sampler documentation
    !>  for detailed descriptions of the input simulation specifications at:<br>
    !>  <ol>
    !>      <li>    [ParaDRAM simulation options](https://www.cdslab.org/paramonte/notes/usage/paradram/specifications/)
    !>      <li>    [ParaDISE simulation options](https://www.cdslab.org/paramonte/notes/usage/paradise/specifications/)
    !>      <li>    [ParaNest simulation options](https://www.cdslab.org/paramonte/notes/usage/paranest/specifications/)
    !>  </ol>
    !>
    !>  \see
    !>  [paradram_type](@ref pm_sampling::paradram_type)<br>
    !>  [paranest_type](@ref pm_sampling::paranest_type)<br>
    !>  [paradise_type](@ref pm_sampling::paradise_type)<br>
    !>  [getErrSampling](@ref pm_sampling::getErrSampling)<br>
    !>
    !>  \finmain{sampler_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type                                        :: sampler_type
        character(:,SKC)        , allocatable   :: description
        character(:,SKC)        , allocatable   :: domain
        character(:,SKC)        , allocatable   :: domainAxisName(:)
        real(RKH)               , allocatable   :: domainBallCenter(:)
        real(RKH)               , allocatable   :: domainBallCorMat(:,:)
        real(RKH)               , allocatable   :: domainBallCovMat(:,:)
        real(RKH)               , allocatable   :: domainBallStdVec(:)
        real(RKH)               , allocatable   :: domainCubeLimitLower(:)
        real(RKH)               , allocatable   :: domainCubeLimitUpper(:)
        integer(IK)             , allocatable   :: domainErrCount
        integer(IK)             , allocatable   :: domainErrCountMax
        character(:,SKC)        , allocatable   :: inputFile                !<  \public The input scalar `character` of default kind \SK which be either,
                                                                            !!  <ol>
                                                                            !!      <li>    The path to an external input `namelist` file containing the simulation settings.<br>
                                                                            !!      <li>    A string containing the simulation settings in `namelist` format.<br>
                                                                            !!              The following rules apply to variable assignments in namelist-style input:<br>
                                                                            !!              <ol>
                                                                            !!                  <li>    Comments must begin with an exclamation mark `!`.
                                                                            !!                  <li>    Comments can appear anywhere on an empty line or, after a variable assignment
                                                                            !!                          (but not in the middle of a variable assignment whether in single or multiple lines).
                                                                            !!                  <li>    All variable assignments are optional and can be commented out. In such cases, appropriate default values will be used.
                                                                            !!                  <li>    All variable assignments must appear within a `namelist` group which starts with either
                                                                            !!                          <ol>
                                                                            !!                              <li>    the (case-INsensitive) keyword `&paradise` to be used for specifying [paradram_type](@ref pm_sampling::paradram_type) simulation specifications.<br>
                                                                            !!                              <li>    the (case-INsensitive) keyword `&paradram` to be used for specifying [paradram_type](@ref pm_sampling::paradram_type) simulation specifications.<br>
                                                                            !!                              <li>    the (case-INsensitive) keyword `&paranest` to be used for specifying [paradram_type](@ref pm_sampling::paradram_type) simulation specifications.<br>
                                                                            !!                          </ol>
                                                                            !!                          and ends with a forward slash `/` on the last line of the namelist group.
                                                                            !!                  <li>    Every time a sampler opens an input file, it will first search for the namelist group with the corresponding sampler name.<br>
                                                                            !!                  <li>    The order of the input variables in the namelist groups is irrelevant and unimportant.
                                                                            !!                  <li>    Variables can be defined multiple times, but **only the last definition** will be considered as input.
                                                                            !!                  <li>    All variable names are case-insensitive.
                                                                            !!                          However, for clarity, this software follows the camelCase code-writing practice for simulation specification names.
                                                                            !!                  <li>    String values must be enclosed with either single or double quotation marks.
                                                                            !!                  <li>    Logical values are case-insensitive and can be either<br>
                                                                            !!                          <ol>
                                                                            !!                              <li>    `.true.`, `true`, or `t` (all case-insensitive) representing the truth, or,
                                                                            !!                              <li>    `.false.`, `false`, or `f` (all case-insensitive) representing the falsehood.
                                                                            !!                          </ol>
                                                                            !!                  <li>    All vectors and arrays specified within the input file start with a lower-bound of `1` as opposed to `0` in the C family of languages.<br>
                                                                            !!                          This follows the convention of the majority of science-oriented programming languages: Fortran, Julia, Mathematica, MATLAB, and R.
                                                                            !!                  <li>    Specifying multiple values is as simple as listing the elements as comma or space -separated values.
                                                                            !!                  <li>    Individual vector or element values can be specified by an explicit index or index range:<br>
                                                                            !!                          \code{.F90}
                                                                            !!                              &paranest
                                                                            !!                                  domainAxisName(1 : 2) = "sampleVariable1", "sampleVariable2"
                                                                            !!                                  domainAxisName(3) = "sampleVariable3"
                                                                            !!                              /
                                                                            !!                          \endcode
                                                                            !!                  <li>    If a sequence of elements of a vector or array are equal,
                                                                            !!                          the assignment can be succinctly expressed via the following convention:<br>
                                                                            !!                          \code{.F90}
                                                                            !!                              &paradram
                                                                            !!                                  proposalStart = 4 * 0. ! assign value `0.` to the first four elements of proposalStart.
                                                                            !!                              /
                                                                            !!                          \endcode
                                                                            !!                  <li>    Matrix assignments follow the [column-major](https://en.wikipedia.org/wiki/Row-_and_column-major_order) storage scheme in memory.<br>
                                                                            !!                          \code{.F90}
                                                                            !!                              &paradram
                                                                            !!                                  proposalCorMat = 1, 0, 0, 0
                                                                            !!                                                 , 0, 1, 0, 0
                                                                            !!                                                 , 0, 0, 1, 0
                                                                            !!                                                 , 0, 0, 0, 1
                                                                            !!                              /
                                                                            !!                          \endcode
                                                                            !!              </ol>
                                                                            !!  </ol>
                                                                            !!  For detailed descriptions of the input simulation settings see,
                                                                            !!  <ol>
                                                                            !!      <li>    [ParaDRAM simulation options](https://www.cdslab.org/paramonte/notes/usage/paradram/specifications/)
                                                                            !!      <li>    [ParaDISE simulation options](https://www.cdslab.org/paramonte/notes/usage/paradise/specifications/)
                                                                            !!      <li>    [ParaNest simulation options](https://www.cdslab.org/paramonte/notes/usage/paranest/specifications/)
                                                                            !!  </ol>
                                                                            !!  (**optional**. If missing, the default simulation settings will be used.)
       !logical(LK)             , allocatable   :: inputFileHasPriority
        character(:,SKC)        , allocatable   :: outputChainFileFormat
        integer(IK)             , allocatable   :: outputColumnWidth
        character(:,SKC)        , allocatable   :: outputFileName
        character(:,SKC)        , allocatable   :: outputStatus
        integer(IK)             , allocatable   :: outputPrecision
        integer(IK)             , allocatable   :: outputReportPeriod
        character(:,SKC)        , allocatable   :: outputRestartFileFormat
        integer(IK)             , allocatable   :: outputSampleSize
        character(:,SKC)        , allocatable   :: outputSeparator
        character(:,SKC)        , allocatable   :: outputSplashMode
        character(:,SKC)        , allocatable   :: parallelism
        logical(LK)             , allocatable   :: parallelismMpiFinalizeEnabled
        integer(IK)             , allocatable   :: parallelismNumThread
       !character(:,SKC)        , allocatable   :: plang
        integer(IK)             , allocatable   :: randomSeed
        character(:,SKC)        , allocatable   :: sysInfoFilePath
        real(RKH)               , allocatable   :: targetAcceptanceRate(:)
    end type

    !>  \brief
    !>  This is a derived type for constructing objects containing
    !>  the **optional** simulation properties of the ParaMonte library MCMC explorers and samplers.<br>
    !>
    !>  \details
    !>  Objects of this derived type are not meant to be used directly by the end users.<br>
    !>  Instead use one of the appropriate derived types:<br>
    !>  <ol>
    !>      <li>    [paradram_type](@ref pm_sampling::paradram_type)
    !>      <li>    [paranest_type](@ref pm_sampling::paranest_type)
    !>      <li>    [paradise_type](@ref pm_sampling::paradise_type)
    !>  </ol>
    !>
    !>  \note
    !>  For more information on any of the components of this derived type,
    !>  see the output `_report.txt` files from the corresponding sampler simulations.<br>
    !>  Alternatively, see the ParaMonte library cross-language sampler documentation
    !>  for detailed descriptions of the input simulation specifications at:<br>
    !>  <ol>
    !>      <li>    [ParaDRAM simulation options](https://www.cdslab.org/paramonte/notes/usage/paradram/specifications/)
    !>      <li>    [ParaDISE simulation options](https://www.cdslab.org/paramonte/notes/usage/paradise/specifications/)
    !>      <li>    [ParaNest simulation options](https://www.cdslab.org/paramonte/notes/usage/paranest/specifications/)
    !>  </ol>
    !>
    !>  \see
    !>  [paradram_type](@ref pm_sampling::paradram_type)<br>
    !>  [paranest_type](@ref pm_sampling::paranest_type)<br>
    !>  [paradise_type](@ref pm_sampling::paradise_type)<br>
    !>  [getErrSampling](@ref pm_sampling::getErrSampling)<br>
    !>
    !>  \finmain{paramcmc_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(sampler_type)                 :: paramcmc_type
        integer(IK)             , allocatable   :: outputChainSize
        integer(IK)             , allocatable   :: outputSampleRefinementCount
        character(:,SKC)        , allocatable   :: outputSampleRefinementMethod
        character(:,SKC)        , allocatable   :: proposal
        real(RKH)               , allocatable   :: proposalCorMat(:,:)
        real(RKH)               , allocatable   :: proposalCovMat(:,:)
        character(:,SKC)        , allocatable   :: proposalScaleFactor
        real(RKH)               , allocatable   :: proposalStart(:)
        real(RKH)               , allocatable   :: proposalStartDomainCubeLimitLower(:)
        real(RKH)               , allocatable   :: proposalStartDomainCubeLimitUpper(:)
        logical(LK)             , allocatable   :: proposalStartRandomized
        real(RKH)               , allocatable   :: proposalStdVec(:)
    end type

    !>  \brief
    !>  This is a derived type for constructing objects containing the **optional** simulation
    !>  properties of the ParaMonte library Delayed-Rejection Adaptive Metropolis (DRAM) MCMC explorers and samplers.<br>
    !>
    !>  \note
    !>  For more information on any of the components of this derived type,
    !>  see the output `_report.txt` files from the corresponding sampler simulations.<br>
    !>  Alternatively, see the ParaMonte library cross-language sampler documentation
    !>  for detailed descriptions of the input simulation specifications at:<br>
    !>  <ol>
    !>      <li>    [ParaDRAM simulation options](https://www.cdslab.org/paramonte/notes/usage/paradram/specifications/)
    !>      <li>    [ParaDISE simulation options](https://www.cdslab.org/paramonte/notes/usage/paradise/specifications/)
    !>      <li>    [ParaNest simulation options](https://www.cdslab.org/paramonte/notes/usage/paranest/specifications/)
    !>  </ol>
    !>
    !>  \see
    !>  [paradram_type](@ref pm_sampling::paradram_type)<br>
    !>  [paranest_type](@ref pm_sampling::paranest_type)<br>
    !>  [paradise_type](@ref pm_sampling::paradise_type)<br>
    !>  [getErrSampling](@ref pm_sampling::getErrSampling)<br>
    !>
    !>  \finmain{paradram_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(paramcmc_type)                :: paradram_type
        real(RKH)               , allocatable   :: burninAdaptationMeasure
        integer(IK)             , allocatable   :: proposalAdaptationCount
        integer(IK)             , allocatable   :: proposalAdaptationCountGreedy
        integer(IK)             , allocatable   :: proposalAdaptationPeriod
        integer(IK)             , allocatable   :: proposalDelayedRejectionCount
        real(RKH)               , allocatable   :: proposalDelayedRejectionScaleFactor(:)
    end type

    !>  \brief
    !>  This is a derived type for constructing objects containing the **optional** simulation
    !>  properties of the ParaMonte library DRAM-enhanced MCMC explorers and samplers.<br>
    !>
    !>  \note
    !>  For more information on any of the components of this derived type,
    !>  see the output `_report.txt` files from the corresponding sampler simulations.<br>
    !>  Alternatively, see the ParaMonte library cross-language sampler documentation
    !>  for detailed descriptions of the input simulation specifications at:<br>
    !>  <ol>
    !>      <li>    [ParaDRAM simulation options](https://www.cdslab.org/paramonte/notes/usage/paradram/specifications/)
    !>      <li>    [ParaDISE simulation options](https://www.cdslab.org/paramonte/notes/usage/paradise/specifications/)
    !>      <li>    [ParaNest simulation options](https://www.cdslab.org/paramonte/notes/usage/paranest/specifications/)
    !>  </ol>
    !>
    !>  \see
    !>  [paradram_type](@ref pm_sampling::paradram_type)<br>
    !>  [paranest_type](@ref pm_sampling::paranest_type)<br>
    !>  [paradise_type](@ref pm_sampling::paradise_type)<br>
    !>  [getErrSampling](@ref pm_sampling::getErrSampling)<br>
    !>
    !>  \finmain{paradise_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(paradram_type)                :: paradise_type
    end type

    !>  \brief
    !>  This is a derived type for constructing objects containing the **optional** simulation
    !>  properties of the ParaMonte library Nested explorers and samplers.<br>
    !>
    !>  \note
    !>  For more information on any of the components of this derived type,
    !>  see the output `_report.txt` files from the corresponding sampler simulations.<br>
    !>  Alternatively, see the ParaMonte library cross-language sampler documentation
    !>  for detailed descriptions of the input simulation specifications at:<br>
    !>  <ol>
    !>      <li>    [ParaDRAM simulation options](https://www.cdslab.org/paramonte/notes/usage/paradram/specifications/)
    !>      <li>    [ParaDISE simulation options](https://www.cdslab.org/paramonte/notes/usage/paradise/specifications/)
    !>      <li>    [ParaNest simulation options](https://www.cdslab.org/paramonte/notes/usage/paranest/specifications/)
    !>  </ol>
    !>
    !>  \see
    !>  [paradram_type](@ref pm_sampling::paradram_type)<br>
    !>  [paranest_type](@ref pm_sampling::paranest_type)<br>
    !>  [paradise_type](@ref pm_sampling::paradise_type)<br>
    !>  [getErrSampling](@ref pm_sampling::getErrSampling)<br>
    !>  [getErrSampling](@ref pm_sampling::getErrSampling)<br>
    !>
    !>  \finmain{paradise_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(sampler_type)                 :: paranest_type
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
        character(:,SKC)        , allocatable   :: domainPartitionMethod
        character(:,SKC)        , allocatable   :: domainPartitionObject
        logical(LK)             , allocatable   :: domainPartitionOptimizationScaleEnabled
        logical(LK)             , allocatable   :: domainPartitionOptimizationShapeEnabled
        logical(LK)             , allocatable   :: domainPartitionOptimizationShapeScaleEnabled
        real(RKH)               , allocatable   :: domainPartitionScaleFactor
        character(:,SKC)        , allocatable   :: domainSampler
        integer(IK)             , allocatable   :: liveSampleSize
        real(RKH)               , allocatable   :: tolerance
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !   This abstract interface defines the set of objective function interfaces
    !   that can be passed to the ParaMonte samplers (e.g., ParaDRAM, ParaNest, ...).
    !   These interfaces only have internal library specific use cases.
    !>  \cond excluded
    abstract interface

#if OMP_ENABLED && (MATLAB_ENABLED || PYTHON_ENABLED || R_ENABLED)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    function getLogFunc_proc_RK5(logFuncState, avgTimePerFunCall, avgCommPerFunCall) result(mold)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogFunc_proc_RK5
#endif
        use pm_kind, only: RKD, RKC => RK5
        real(RKC), intent(inout), contiguous :: logFuncState(0:,:)
        real(RKD), intent(inout) :: avgTimePerFunCall, avgCommPerFunCall
        real(RKC) :: mold
    end function
#endif

#if RK4_ENABLED
    function getLogFunc_proc_RK4(logFuncState, avgTimePerFunCall, avgCommPerFunCall) result(mold)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogFunc_proc_RK4
#endif
        use pm_kind, only: RKD, RKC => RK4
        real(RKC), intent(inout), contiguous :: logFuncState(0:,:)
        real(RKD), intent(inout) :: avgTimePerFunCall, avgCommPerFunCall
        real(RKC) :: mold
    end function
#endif

#if RK3_ENABLED
    function getLogFunc_proc_RK3(logFuncState, avgTimePerFunCall, avgCommPerFunCall) result(mold)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogFunc_proc_RK3
#endif
        use pm_kind, only: RKD, RKC => RK3
        real(RKC), intent(inout), contiguous :: logFuncState(0:,:)
        real(RKD), intent(inout) :: avgTimePerFunCall, avgCommPerFunCall
        real(RKC) :: mold
    end function
#endif

#if RK2_ENABLED
    function getLogFunc_proc_RK2(logFuncState, avgTimePerFunCall, avgCommPerFunCall) result(mold)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogFunc_proc_RK2
#endif
        use pm_kind, only: RKD, RKC => RK2
        real(RKC), intent(inout), contiguous :: logFuncState(0:,:)
        real(RKD), intent(inout) :: avgTimePerFunCall, avgCommPerFunCall
        real(RKC) :: mold
    end function
#endif

#if RK1_ENABLED
    function getLogFunc_proc_RK1(logFuncState, avgTimePerFunCall, avgCommPerFunCall) result(mold)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogFunc_proc_RK1
#endif
        use pm_kind, only: RKD, RKC => RK1
        real(RKC), intent(inout), contiguous :: logFuncState(0:,:)
        real(RKD), intent(inout) :: avgTimePerFunCall, avgCommPerFunCall
        real(RKC) :: mold
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#else

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    function getLogFunc_proc_RK5(state) result(logFunc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogFunc_proc_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC), intent(in), contiguous :: state(:)
        real(RKC) :: logFunc
    end function
#endif

#if RK4_ENABLED
    function getLogFunc_proc_RK4(state) result(logFunc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogFunc_proc_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC), intent(in), contiguous :: state(:)
        real(RKC) :: logFunc
    end function
#endif

#if RK3_ENABLED
    function getLogFunc_proc_RK3(state) result(logFunc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogFunc_proc_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC), intent(in), contiguous :: state(:)
        real(RKC) :: logFunc
    end function
#endif

#if RK2_ENABLED
    function getLogFunc_proc_RK2(state) result(logFunc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogFunc_proc_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC), intent(in), contiguous :: state(:)
        real(RKC) :: logFunc
    end function
#endif

#if RK1_ENABLED
    function getLogFunc_proc_RK1(state) result(logFunc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogFunc_proc_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC), intent(in), contiguous :: state(:)
        real(RKC) :: logFunc
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#endif

    end interface
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the procedure fails to fully accomplish the task of generating a
    !>  Monte Carlo sample of the specified input mathematical objective function, otherwise, return `.false.`.
    !>
    !>  \details
    !>  This generic interface is the entry point to all ParaMonte Monte Carlo samplers of mathematical density functions.<br>
    !>  Although the procedures of this generic interface return a single scalar object of type [err_type](@ref pm_err::err_type),
    !>  the procedures generate massive amounts of information about each simulation which are stored in appropriate external hard drive files.<br>
    !>
    !>  \param[in]  sampler     :   The input scalar object that can be of,<br>
    !>                              <ol>
    !>                                  <li>    type [paradram_type](@ref pm_sampling::paradram_type) indicating the use of the ParaDRAM sampler,<br>
    !>                                  <li>    type [paradise_type](@ref pm_sampling::paradise_type) indicating the use of the ParaDISE sampler,<br>
    !>                                  <li>    type [paranest_type](@ref pm_sampling::paranest_type) indicating the use of the ParaNest sampler,<br>
    !>                              </ol>
    !>  \param[in]  getLogFunc  :   The input user-specified `function` that takes an input `contiguous` vector of shape `state(1:ndim)` of type `real` of kind \RKALL.<br>
    !>                              On input, `state` represents a state (point) from within the domain of the user-specified target density function whose function value must be returned.<br>
    !>                              On output, the user-specified procedure `getLogFunc()` must return the function value corresponding to the input `state(1:ndim)`.<br>
    !>                              The following illustrate the generic interface of `getLogFunc(state)`,
    !>                              \code{.F90}
    !>                                  function getLogFunc(state) result(logFunc)
    !>                                      real(RKC), intent(in), contiguous :: state(:) ! (1 : ndim)
    !>                                      real(RKC) :: logFunc
    !>                                  end function
    !>                              \endcode
    !>                              where the condition `all(shape(state) == [ndim])` holds where `ndim` is the number of dimensions of the target density function.<br>
    !>                              **If you are an end user, this is where you stop reading about this input argument**.<br>
    !>                              If you are a developer, continue reading about the alternative internal interfaces of this argument for dynamic (MATLAB/Python/R) languages.<br>
    !>                              Setting the preprocessing flag `CFI_ENABLED` to a non-zero value at the time of library compilation
    !>                              activates the C-Fortran Interoperability (CFI) interface of this generic interface.<br>
    !>                              This requires `getLogFunc()` to have the following interfaces:<br>
    !>                              <ol>
    !>                                  <li>    When the library is built for OpenMP parallelism (i.e., `OMP_ENABLED=1` holds) for MATLAB/Python/R programming languages,
    !>                                          the input `getLogFunc()` must have the following interface:<br>
    !>                                          \code{.F90}
    !>                                              function getLogFunc(logFuncState, ndim, njob, avgTimePerFunCall, avgCommPerFunCall) result(mold)
    !>                                                  integer(IK), value :: ndim, njob
    !>                                                  real(RKC), intent(inout) :: logFuncState(ndim + 1, njob)
    !>                                                  real(RKD), intent(inout) :: avgTimePerFunCall, avgCommPerFunCall
    !>                                                  real(RKC) :: mold
    !>                                              end function
    !>                                          \endcode
    !>                                          where,
    !>                                          <ol>
    !>                                              <li>    the variable `ndim` is the number of dimensions of the target density function and
    !>                                                      the variable `njob` is the number of states for which the function must be evaluated concurrently.
    !>                                              <li>    the constant \RKD refers to the `real` kind `double precision` supported by the compiler.
    !>                                              <li>    the output `mold` merely dictates the `real` kind of the interface and its value is ignored.
    !>                                              <li>    on input, the subset `logFuncState(2:ndim+1, 1:njob)` represents the set of `njob` states (points)
    !>                                                      from within the domain of the user-specified target density function whose function values must be returned.<br>
    !>                                                      On output, the user-specified procedure `getLogFunc()` must return `njob` function values in elements `logFuncState(1, 1:njob)`
    !>                                                      corresponding to the input states `logFuncState(2:ndim+1, 1:njob)`.<br>
    !>                                              <li>    the two extra `intent(inout)` arguments `avgTimePerFunCall, avgCommPerFunCall` represent:<br>
    !>                                                      <ol>
    !>                                                          <li>    the average time (in seconds) it takes to compute the user-specified function for a **single** input state.
    !>                                                          <li>    the average time (in seconds) it takes to create images/processes/threads for parallel computation of
    !>                                                                  the user-specified function and communicate information among them, collectively called **parallelization overhead**.
    !>                                                      </ol>
    !>                                                      These two extra arguments are essential for accurate post-processing the performance of the parallel computations in the simulation.<br>
    !>                                                      If there is no parallelization within he user-specified `getLogFunc`, then the parallelization overhead can be set to zero.<br>
    !>                                          </ol>
    !>                                          This procedure interface is entirely different from the density function interfaces in the previous versions of the ParaMonte library.<br>
    !>                                          The new interface has been primarily developed to facilitate the **concurrent or parallel** evaluation
    !>                                          of target density function in higher-level languages for a set of `njob` input states.<br>
    !>                                          This interface is particularly vital in higher-level programming languages such as MATLAB, Python, and R where ParaMonte thread-level
    !>                                          parallelism is impossible due to the [global interpreter lock (GIL)](https://en.wikipedia.org/wiki/Global_interpreter_lock).<br>
    !>                                          This interface is currently activated **if and only if** the preprocessors `OMP_ENABLED=1` and either `MATLAB_ENABLED=1` or `PYTHON_ENABLED=1` or `R_ENABLED=1` are set.<br>
    !>                                  <li>    When the library is built for any build configuration other than the above with preprocessor `CDF_ENABLED=1`, particularly, for the C and C++ programming languages,
    !>                                          the input `getLogFunc()` must have the following interface:<br>
    !>                                          \code{.F90}
    !>                                              function getLogFunc(state, ndim) result(logFunc)
    !>                                                  integer(IK), value :: ndim
    !>                                                  real(RKC), intent(in) :: state(ndim)
    !>                                                  real(RKC) :: logFunc
    !>                                              end function
    !>                                          \endcode
    !>                                          where,
    !>                                          <ol>
    !>                                              <li>    the variable `ndim` is the number of dimensions of the target density function.
    !>                                              <li>    the variable `state` is an input point from within the domain of the density function.
    !>                                              <li>    the constant `RKC` refers to the `real` kind used for computations.
    !>                                              <li>    the output `logFunc` is the computed value of the density function at the specified input `state`.
    !>                                          </ol>
    !>                              </ol>
    !>  \param[in]  ndim        :   The input scalar `integer` of default kind \IK representing the number of
    !>                              dimensions of the domain of the objective function that is to be sampled.<br>
    !>
    !>  \return
    !>  `err`                   :   The output scalar of type [err_type](@ref pm_err::err_type) whose component `occur` is set to the `.true.` **if and only if** the
    !>                              procedure fails to accomplish the sampling task, in which case, the `msg` component of `err` is set to a description of the error origin.<br>
    !>                              In such a case, the description of the possible cause(s) of the error will be written to the output *report* file that is automatically generated for all simulations.<br>
    !>                              The primary reason for simulation failures is a wrong syntactic or semantic input value(s) for the simulation settings specified within the components of the input
    !>                              argument `sampler` or within the external output file specified by the `inputFile` component of the input `sampler`.<br>
    !>                              On output, the `msg` component of `err` is an empty string if the sampling completes successfully.<br>
    !>
    !>  \interface{getErrSampling}
    !>  \code{.F90}
    !>
    !>      use pm_sampling, only: getErrSampling
    !>
    !>      err = getErrSampling(sampler, getLogFunc, ndim)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < ndim` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \example{getErrSampling}
    !>  \include{lineno} example/pm_sampling/getErrSampling/main.F90
    !>  \compilef{getErrSampling}
    !>  \output{getErrSampling}
    !>  \include{lineno} example/pm_sampling/getErrSampling/main.out.F90
    !>
    !>  \example{mvn}
    !>  \include{lineno} example/pm_sampling/mvn/main.F90
    !>  \inputnml{mvn, input.nml}
    !>  \include{lineno} example/pm_sampling/mvn/input.nml
    !>  \compilef{mvn}
    !>  \postproc{mvn}
    !>  \include{lineno} example/pm_sampling/mvn/main.py
    !>  \vis{mvn}
    !>  \image html example/pm_sampling/mvn/mvn.traceplot.png width=700
    !>  \image html example/pm_sampling/mvn/mvn.scatterplot.png width=700
    !>  \image html example/pm_sampling/mvn/mvn.adaptationMeasure.png width=700
    !>
    !>  \example{himmelblau}
    !>  \include{lineno} example/pm_sampling/himmelblau/main.F90
    !>  \inputnml{himmelblau, input.nml}
    !>  \include{lineno} example/pm_sampling/himmelblau/input.nml
    !>  \compilef{himmelblau}
    !>  \postproc{himmelblau}
    !>  \include{lineno} example/pm_sampling/himmelblau/main.py
    !>  \vis{himmelblau}
    !>  \image html example/pm_sampling/himmelblau/himmelblau.traceplot.png width=700
    !>  \image html example/pm_sampling/himmelblau/himmelblau.scatterplot.png width=700
    !>  \image html example/pm_sampling/himmelblau/himmelblau.adaptationMeasure.png width=700
    !>
    !>  \example{mixture}
    !>  \include{lineno} example/pm_sampling/mixture/main.F90
    !>  \inputnml{mixture, input.nml}
    !>  \include{lineno} example/pm_sampling/mixture/input.nml
    !>  \compilef{mixture}
    !>  \postproc{mixture}
    !>  \include{lineno} example/pm_sampling/mixture/main.py
    !>  \vis{mixture}
    !>  \image html example/pm_sampling/mixture/mixLogNormLogNorm.traceplot.png width=700
    !>  \image html example/pm_sampling/mixture/mixLogNormLogNorm.scatterplot.png width=700
    !>  \image html example/pm_sampling/mixture/mixLogNormLogNorm.adaptationMeasure.png width=700
    !>  \image html example/pm_sampling/mixture/mixLogNormLogNorm.hist.png width=700
    !>  <br>
    !>  \image html example/pm_sampling/mixture/mixFlatPowetoFlatPowetTapered.traceplot.png width=700
    !>  \image html example/pm_sampling/mixture/mixFlatPowetoFlatPowetTapered.scatterplot.png width=700
    !>  \image html example/pm_sampling/mixture/mixFlatPowetoFlatPowetTapered.adaptationMeasure.png width=700
    !>  \image html example/pm_sampling/mixture/mixFlatPowetoFlatPowetTapered.hist.png width=700
    !>  <br>
    !>  \image html example/pm_sampling/mixture/mixLogNormFlatPowetTapered.traceplot.png width=700
    !>  \image html example/pm_sampling/mixture/mixLogNormFlatPowetTapered.scatterplot.png width=700
    !>  \image html example/pm_sampling/mixture/mixLogNormFlatPowetTapered.adaptationMeasure.png width=700
    !>  \image html example/pm_sampling/mixture/mixLogNormFlatPowetTapered.hist.png width=700
    !>  <br>
    !>  \image html example/pm_sampling/mixture/mixFlatPowetLogNorm.traceplot.png width=700
    !>  \image html example/pm_sampling/mixture/mixFlatPowetLogNorm.scatterplot.png width=700
    !>  \image html example/pm_sampling/mixture/mixFlatPowetLogNorm.adaptationMeasure.png width=700
    !>  \image html example/pm_sampling/mixture/mixFlatPowetLogNorm.hist.png width=700
    !>
    !>  \test
    !>  [test_pm_sampling](@ref test_pm_sampling)
    !>
    !>  \finmain{getErrSampling}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas at Austin
    interface getErrSampling

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getErrParaDRAM_RK5(sampler, getLogFunc, ndim) result(err)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getErrParaDRAM_RK5
#endif
        use pm_kind, only: RKC => RK5
        procedure(getLogFunc_proc_RK5) :: getLogFunc
        type(paradram_type) , intent(in) :: sampler
        integer(IK), intent(in) :: ndim
        type(err_type) :: err
    end function
#endif

#if RK4_ENABLED
    impure module function getErrParaDRAM_RK4(sampler, getLogFunc, ndim) result(err)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getErrParaDRAM_RK4
#endif
        use pm_kind, only: RKC => RK4
        procedure(getLogFunc_proc_RK4) :: getLogFunc
        type(paradram_type) , intent(in) :: sampler
        integer(IK), intent(in) :: ndim
        type(err_type) :: err
    end function
#endif

#if RK3_ENABLED
    impure module function getErrParaDRAM_RK3(sampler, getLogFunc, ndim) result(err)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getErrParaDRAM_RK3
#endif
        use pm_kind, only: RKC => RK3
        procedure(getLogFunc_proc_RK3) :: getLogFunc
        type(paradram_type) , intent(in) :: sampler
        integer(IK), intent(in) :: ndim
        type(err_type) :: err
    end function
#endif

#if RK2_ENABLED
    impure module function getErrParaDRAM_RK2(sampler, getLogFunc, ndim) result(err)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getErrParaDRAM_RK2
#endif
        use pm_kind, only: RKC => RK2
        procedure(getLogFunc_proc_RK2) :: getLogFunc
        type(paradram_type) , intent(in) :: sampler
        integer(IK), intent(in) :: ndim
        type(err_type) :: err
    end function
#endif

#if RK1_ENABLED
    impure module function getErrParaDRAM_RK1(sampler, getLogFunc, ndim) result(err)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getErrParaDRAM_RK1
#endif
        use pm_kind, only: RKC => RK1
        procedure(getLogFunc_proc_RK1) :: getLogFunc
        type(paradram_type) , intent(in) :: sampler
        integer(IK), intent(in) :: ndim
        type(err_type) :: err
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function runParaDRAML(getLogFuncPtr, ndim, input) result(failed) bind(C, name = "runParaDRAML")
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: runParaDRAML
#endif
        use, intrinsic :: iso_c_binding, only: SK => c_char, IK => c_int32_t, c_long_double
        integer, parameter :: RKC = merge(c_long_double, RKH, any(c_long_double == RKALL))
#define runParaDRAM_ENABLED 1
#include "pm_sampling@routines.inc.F90"
#undef  runParaDRAM_ENABLED
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function runParaDRAMD(getLogFuncPtr, ndim, input) result(failed) bind(C, name = "runParaDRAMD")
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: runParaDRAMD
#endif
        use, intrinsic :: iso_c_binding, only: SK => c_char, IK => c_int32_t, c_double
        integer, parameter :: RKC = merge(c_double, RKD, any(c_double == RKALL))
#define runParaDRAM_ENABLED 1
#include "pm_sampling@routines.inc.F90"
#undef  runParaDRAM_ENABLED
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function runParaDRAMF(getLogFuncPtr, ndim, input) result(failed) bind(C, name = "runParaDRAMF")
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: runParaDRAMF
#endif
        use pm_kind, only: RKALL
        use, intrinsic :: iso_c_binding, only: SK => c_char, IK => c_int32_t, c_float
        integer, parameter :: RKC = merge(c_float, RKL, any(c_float == RKALL))
#define runParaDRAM_ENABLED 1
#include "pm_sampling@routines.inc.F90"
#undef  runParaDRAM_ENABLED
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_sampling ! LCOV_EXCL_LINE