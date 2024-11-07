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

        ! \bug
        ! Avoid Intel ifort bug for too many `use` statements in a submodule by placing them all in the submodule header.
#if     CHECK_ENABLED
        use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG)call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif
        ! Define intel-specific shared property of files on Windows.
#if     INTEL_ENABLED && WINDOWS_ENABLED
#define SHARED, shared
#else
#define SHARED
#endif

    use pm_err, only: getFine
    use pm_err, only: err_type
    use pm_kind, only: SK, IK, LK, SKG => SK
    use pm_sysPath, only: PATHLEN => MAX_LEN_FILE_PATH
    use pm_strASCII, only: getStrLower, setStrLower
    use pm_arrayReplace, only: getReplaced
    use pm_io, only: setContentsFrom
    use pm_val2str, only: getStr

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_sampling_scio"

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !   The namelist specifications of this module are all internal ParaMonte libraries and not meant to be accessible to end users.
    !   This utter io mess is entirely due to shear inadequacy of the namelist feature and
    !   the standard committee inaction to expand and improve this valuable feature.
    !   Unfortunately, the following variables must remain public, because they are
    !   used by the following internal modules:
    !   pm_sampling_base
    !   pm_sampling_mcmc
    !   pm_sampling_dram
    !   pm_sampling_nest

    ! specbase
    !>  \cond excluded
    character(8191,SKG)                     :: description ! roughly 66Kb of memory
    character(63,SKG)                       :: domain
    character(63,SKG)       , allocatable   :: domainAxisName(:)
    real(RKG)               , allocatable   :: domainBallAvg(:)
    real(RKG)               , allocatable   :: domainBallCor(:,:)
    real(RKG)               , allocatable   :: domainBallCov(:,:)
    real(RKG)               , allocatable   :: domainBallStd(:)
    real(RKG)               , allocatable   :: domainCubeLimitLower(:)
    real(RKG)               , allocatable   :: domainCubeLimitUpper(:)
    integer(IK)                             :: domainErrCount
    integer(IK)                             :: domainErrCountMax
    logical(LK)                             :: inputFileHasPriority
    character(15,SKG)                       :: outputChainFileFormat
    integer(IK)                             :: outputColumnWidth
    character(PATHLEN,SKG)                  :: outputFileName
    character(15,SKG)                       :: outputStatus
    integer(IK)                             :: outputPrecision
    integer(IK)                             :: outputReportPeriod
    character(15,SKG)                       :: outputRestartFileFormat
    integer(IK)                             :: outputSampleSize
    character(63,SKG)                       :: outputSeparator
    character(63,SKG)                       :: outputSplashMode
    character(63,SKG)                       :: parallelism
    logical(LK)                             :: parallelismMpiFinalizeEnabled
    integer(IK)                             :: parallelismNumThread
   !character(511,SKG)                      :: plang
    integer(IK)                             :: randomSeed
    character(PATHLEN,SKG)                  :: sysInfoFilePath
    real(RKG)                               :: targetAcceptanceRate(2)
    !>  \endcond excluded

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! specmcmc
    !>  \cond excluded
    integer(IK)                             :: outputChainSize
    integer(IK)                             :: outputSampleRefinementCount
    character(63,SKG)                       :: outputSampleRefinementMethod
    character(63,SKG)                       :: proposal
    real(RKG)               , allocatable   :: proposalCor(:,:)
    real(RKG)               , allocatable   :: proposalCov(:,:)
    character(127,SKG)                      :: proposalScale
    real(RKG)               , allocatable   :: proposalStart(:)
    real(RKG)               , allocatable   :: proposalStartDomainCubeLimitLower(:)
    real(RKG)               , allocatable   :: proposalStartDomainCubeLimitUpper(:)
    logical(LK)                             :: proposalStartRandomized
    real(RKG)               , allocatable   :: proposalStd(:)
    !>  \endcond excluded

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! specdram
    !>  \cond excluded
    real(RKG)                               :: proposalAdaptationBurnin
    integer(IK)                             :: proposalAdaptationCount
    integer(IK)                             :: proposalAdaptationCountGreedy
    integer(IK)                             :: proposalAdaptationPeriod
    integer(IK)                             :: proposalDelayedRejectionCount
    real(RKG)               , allocatable   :: proposalDelayedRejectionScale(:)
    !>  \endcond excluded

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! specnest
    !>  \cond excluded
    integer(IK)                             :: domainPartitionAdaptationCount
    integer(IK)                             :: domainPartitionAdaptationPeriod
    logical(LK)                             :: domainPartitionBiasCorrectionEnabled
    integer(IK)                             :: domainPartitionCountMax
    real(RKG)                               :: domainPartitionFactorExpansion
    real(RKG)                               :: domainPartitionFactorShrinkage
    integer(IK)                             :: domainPartitionKmeansClusterCountMax
    integer(IK)                             :: domainPartitionKmeansClusterSizeMin
    logical(LK)                             :: domainPartitionKmeansNormalizationEnabled
    integer(IK)                             :: domainPartitionKmeansNumFailMax
    integer(IK)                             :: domainPartitionKmeansNumRecursionMax
    integer(IK)                             :: domainPartitionKmeansNumTry
    integer(IK)                             :: domainPartitionKvolumeNumRecursionMax
    real(RKG)                               :: domainPartitionKvolumeWeightExponent
    character(63,SKG)                       :: domainPartitionMethod
    character(63,SKG)                       :: domainPartitionObject
    logical(LK)                             :: domainPartitionOptimizationScaleEnabled
    logical(LK)                             :: domainPartitionOptimizationShapeEnabled
    logical(LK)                             :: domainPartitionOptimizationShapeScaleEnabled
    real(RKG)                               :: domainPartitionScaleFactor
    character(63,SKG)                       :: domainSampler
    integer(IK)                             :: liveSampleSize
    real(RKG)                               :: tolerance
    !>  \endcond excluded

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define namelist_ENABLED 1

#define NAMELIST paramonte
#include "pm_sampling_scio.inc.F90"
#undef NAMELIST

#define NAMELIST paradram
#include "pm_sampling_scio.inc.F90"
#undef NAMELIST

#define NAMELIST paradise
#include "pm_sampling_scio.inc.F90"
#undef NAMELIST

#define NAMELIST paranest
#include "pm_sampling_scio.inc.F90"
#undef NAMELIST

#undef namelist_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type                                    :: cfhbase_type
        character(31,SKG)                   :: processID = "processID"
        character(31,SKG)                   :: meanAcceptanceRate = "meanAcceptanceRate"
        character(31,SKG)                   :: sampleLogFunc = "sampleLogFunc"
    end type
    type(cfhbase_type)      , parameter     :: cfhbase = cfhbase_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: cfhbase
#endif

    !>  \brief
    !>  This is the `abstract` derived type containing the components for storing the basic properties of output chain files by [ParaMonte samplers](@ref pm_sampling).<br>
    !>
    !>  \see
    !>  [getErrChainRead](@ref pm_sampling_scio::getErrChainRead)<br>
    !>  [cfcbase_type](@ref pm_sampling_scio::cfcbase_type)<br>
    !>
    !>  \test
    !>  [test_pm_sampling_scio](@ref test_pm_sampling_scio)<br>
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{13}
    !>  \desc
    !>  Unfortunately, the lack of a fully functional implementation of PDTs in \gfortran{13} prevents the declaration of this derived type as a PDT.<br>
    !>  \remedy
    !>  See the mess with kind-specific sampling modules.
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to a PDT as soon as PDTs are supported bug-free in \gfortran.<br>
    !>
    !>  \final{cfcbase_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    type, abstract                          :: cfcbase_type
        integer(IK)                         :: numDefCol                !<  \public The scalar of type `integer` of default kind \IK. containing the number of the set of default columns that always appear in the chain file before columns corresponding to `sampleState` components appear.
        integer(IK)                         :: nsam                     !<  \public The scalar of type `integer` of default kind \IK containing the chain size which may be less than the size of the chain components below (the number of (potentially weighted) sampled states).<br>
        character(:,SKG)    , allocatable   :: sep                      !<  \public The scalar of type `character` of default kind \SK, containing the field separator used in the chain file.<br>
        character(:,SKG)    , allocatable   :: header                   !<  \public The scalar of type `character` of default kind \SK, containing the delimited chain file column names.<br>
       !type(css_type)      , allocatable   :: colname(:)               !<  \public The vector of type [css_type](@ref pm_container::css_type), containing all chain file column names, such that the column name `sampleLogFunc` correspond to index `0` of colnames.<br>
        integer(IK)         , allocatable   :: processID(:)             !<  \public The `allocatable` vector of type `integer` of default kind \IK, containing the IDs of the processes contributing the corresponding accepted state in the chain file.<br>
        real(RKG)           , allocatable   :: meanAcceptanceRate(:)    !<  \public The `allocatable` vector of type `real` of kind \RKALL, containing the vector of the average acceptance rates at the given point in the chain.<br>
        real(RKG)           , allocatable   :: sampleLogFunc(:)         !<  \public The `allocatable` vector of type `real` of kind \RKALL, containing the natural logarithm of the function value corresponding to the state in the same entry in the chain file.<br>
        real(RKG)           , allocatable   :: sampleState(:,:)         !<  \public The `allocatable` matrix of type `real` of kind \RKALL, each `(1 : ndim, i)` slice of which contains the accepted `ndim`-dimensional state in the chain at the entry `i` of the chain file.<br>
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type, extends(cfhbase_type)             :: cfhmcmc_type
        character(31,SKG)                   :: burninLocation = "burninLocation"
        character(31,SKG)                   :: sampleWeight = "sampleWeight"
    end type
    type(cfhmcmc_type)      , parameter     :: cfhmcmc = cfhmcmc_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: cfhmcmc
#endif

    !type :: size_type
    !    integer(IK) :: compact, verbose
    !end type

    !>  \brief
    !>  This is the `abstract` derived type containing the components for storing the MCMC properties of output chain files by [ParaMonte samplers](@ref pm_sampling).<br>
    !>
    !>  \see
    !>  [getErrChainRead](@ref pm_sampling_scio::getErrChainRead)<br>
    !>  [cfcbase_type](@ref pm_sampling_scio::cfcbase_type)<br>
    !>
    !>  \test
    !>  [test_pm_sampling_scio](@ref test_pm_sampling_scio)<br>
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{13}
    !>  \desc
    !>  Unfortunately, the lack of a fully functional implementation of PDTs in \gfortran{13} prevents the declaration of this derived type as a PDT.<br>
    !>  \remedy
    !>  See the mess with kind-specific sampling modules.
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to a PDT as soon as PDTs are supported bug-free in \gfortran.<br>
    !>
    !>  \final{cfcbase_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    type, extends(cfcbase_type)             :: cfcmcmc_type
        integer(IK)         , allocatable   :: burninLocation(:)    !<  \public The `allocatable` vector of type `integer` of default kind \IK, containing the burnin locations at the given locations in the chains.<br>
        integer(IK)         , allocatable   :: sampleWeight(:)      !<  \public The `allocatable` vector of type `integer` of default kind \IK, containing the weights of the accepted states.<br>
        integer(IK)                         :: sumw                 !<  \public The sum of all sample weights (corresponding to the verbose chain size).<br>
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type, extends(cfhmcmc_type)             :: cfhdram_type
        character(31,SKG)                   :: proposalAdaptation = "proposalAdaptation"
        character(31,SKG)                   :: delayedRejectionStage = "delayedRejectionStage"
    end type
    type(cfhdram_type)      , parameter     :: cfhdram = cfhdram_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: cfhdram
#endif

    character(31,SKG)       , parameter     :: chainFileColNameDRAM(*) =[ cfhdram%processID &
                                                                        , cfhdram%delayedRejectionStage &
                                                                        , cfhdram%meanAcceptanceRate &
                                                                        , cfhdram%proposalAdaptation &
                                                                        , cfhdram%burninLocation &
                                                                        , cfhdram%sampleWeight &
                                                                        , cfhdram%sampleLogFunc &
                                                                        ]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: chainFileColNameDRAM
#endif

    !>  \brief
    !>  This is the `abstract` derived type containing the components for storing the DRAM properties of output chain files by [ParaMonte samplers](@ref pm_sampling).<br>
    !>
    !>  \see
    !>  [getErrChainRead](@ref pm_sampling_scio::getErrChainRead)<br>
    !>  [cfcbase_type](@ref pm_sampling_scio::cfcbase_type)<br>
    !>
    !>  \test
    !>  [test_pm_sampling_scio](@ref test_pm_sampling_scio)<br>
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{13}
    !>  \desc
    !>  Unfortunately, the lack of a fully functional implementation of PDTs in \gfortran{13} prevents the declaration of this derived type as a PDT.<br>
    !>  \remedy
    !>  See the mess with kind-specific sampling modules.
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to a PDT as soon as PDTs are supported bug-free in \gfortran.<br>
    !>
    !>  \final{cfcbase_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    type, extends(cfcmcmc_type)             :: cfcdram_type
        real(RKG)           , allocatable   :: proposalAdaptation(:)     !<  \public The `allocatable` vector of type `real` of kind \RKALL, containing the vector of the adaptation measures at the accepted states.<br>
        integer(IK)         , allocatable   :: delayedRejectionStage(:) !<  \public The `allocatable` vector of type `integer` of default kind \IK, containing the delayed rejection stages at which the proposed states were accepted.<br>
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type, extends(cfhdram_type)             :: cfhdise_type
    end type
    type(cfhdise_type)      , parameter     :: cfhdise = cfhdise_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: cfhdise
#endif

    character(31,SKG)       , parameter     :: chainFileColNameDISE(*) =[ cfhdise%processID &
                                                                        , cfhdise%delayedRejectionStage &
                                                                        , cfhdise%meanAcceptanceRate &
                                                                        , cfhdise%proposalAdaptation &
                                                                        , cfhdise%burninLocation &
                                                                        , cfhdise%sampleWeight &
                                                                        , cfhdise%sampleLogFunc &
                                                                        ]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: chainFileColNameDISE
#endif

    !>  \brief
    !>  This is the `abstract` derived type containing the components for storing the DISE properties of output chain files by [ParaMonte samplers](@ref pm_sampling).<br>
    !>
    !>  \see
    !>  [getErrChainRead](@ref pm_sampling_scio::getErrChainRead)<br>
    !>  [cfcbase_type](@ref pm_sampling_scio::cfcbase_type)<br>
    !>
    !>  \test
    !>  [test_pm_sampling_scio](@ref test_pm_sampling_scio)<br>
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{13}
    !>  \desc
    !>  Unfortunately, the lack of a fully functional implementation of PDTs in \gfortran{13} prevents the declaration of this derived type as a PDT.<br>
    !>  \remedy
    !>  See the mess with kind-specific sampling modules.
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to a PDT as soon as PDTs are supported bug-free in \gfortran.<br>
    !>
    !>  \final{cfcbase_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    type, extends(cfcdram_type)             :: cfcdise_type
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type, extends(cfhbase_type)             :: cfhnest_type
        character(31,SKG)                   :: domainLogVol = "domainLogVol"
        character(31,SKG)                   :: funcLogIntegral = "funcLogIntegral"
        character(31,SKG)                   :: sampleLogWeight = "sampleLogWeight"
    end type
    type(cfhnest_type)      , parameter     :: cfhnest = cfhnest_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: cfhnest
#endif

    character(31,SKG)       , parameter     :: chainFileColNameNest(*) =[ cfhnest%processID &
                                                                        , cfhnest%meanAcceptanceRate &
                                                                        , cfhnest%domainLogVol &
                                                                        , cfhnest%funcLogIntegral &
                                                                        , cfhnest%sampleLogWeight &
                                                                        , cfhnest%sampleLogFunc &
                                                                        ]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: chainFileColNameNest
#endif

    !>  \brief
    !>  This is the derived type containing `allocatable` storage for columns of data in output chain files resulting from [ParaNest](@ref pm_sampling::paranest_type) simulations.
    !>
    !>  \details
    !>  Objects of this derived type can be passed to the generic interface [getErrChainRead](@ref pm_sampling_scio::getErrChainRead)
    !>  to read the contents of a given [ParaNest](@ref pm_sampling::paranest_type) simulation chain file.<br>
    !>  For more details, see the documentation of the components of the derived type.<br>
    !>
    !>  \note
    !>  This derived type does not have a custom constructor as of the most recent version of the ParaMonte library.<br>
    !>
    !>  \see
    !>  [getErrChainRead](@ref pm_sampling_scio::getErrChainRead)<br>
    !>  [cfcbase_type](@ref pm_sampling_scio::cfcbase_type)<br>
    !>
    !>  \test
    !>  [test_pm_sampling_scio](@ref test_pm_sampling_scio)<br>
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{13}
    !>  \desc
    !>  Unfortunately, the lack of a fully functional implementation of PDTs in \gfortran{13} prevents the declaration of this derived type as a PDT.<br>
    !>  \remedy
    !>  See the mess with kind-specific sampling modules.
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to a PDT as soon as PDTs are supported bug-free in \gfortran.<br>
    !>
    !>  \final{chainParaNest_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    type, extends(cfcbase_type)             :: cfcnest_type
        real(RKG)           , allocatable   :: domainLogVol(:)      !<  \public The `allocatable` vector of type `real` of kind \RKALL, containing the natural logarithm of the size of the domain of integration at the corresponding stage in the chain of sampling/integration.<br>
        real(RKG)           , allocatable   :: logMaxRelErr(:)      !<  \public The `allocatable` vector of type `real` of kind \RKALL, containing the natural logarithm of the size of the domain of integration at the corresponding stage in the chain of sampling/integration.<br>
        real(RKG)           , allocatable   :: funcLogIntegral(:)   !<  \public The `allocatable` vector of type `real` of kind \RKALL, containing the natural logarithm of the size of the integral of the objective function at the corresponding stage in the chain.<br>
        real(RKG)           , allocatable   :: sampleLogWeight(:)   !<  \public The `allocatable` vector of type `real` of kind \RKALL, containing the natural logarithm of the contribution of the corresponding state in chain to the integral of the objective function.<br>
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Set the module namelist variables declared in the parent module of this generic interface from the specified `inputFile`.<br>
    !>
    !>  \details
    !>  This is a low-level ParaMonte library interface that is meant to be directly accessible to the end users.<br>
    !>  This routines of this generic interface are used by various explorers and samplers of the ParaMonte library
    !>  to parse the user-specified simulation specifications.<br>
    !>
    !>  \param[in]      method      :   The input scalar object that can be
    !>                                  <ol>
    !>                                      <li>    The constant `"ParaDISE"` indicating the use of the ParaDISE sampler,<br>
    !>                                      <li>    The constant `"ParaDRAM"` indicating the use of the ParaDRAM sampler,<br>
    !>                                      <li>    The constant `"ParaNest"` indicating the use of the ParaNest sampler,<br>
    !>                                  </ol>
    !>  \param[in]      inputFile   :   The input scalar of type `character` of default kind \SK, containing either
    !>                                  <ol>
    !>                                      <li>    The path to the external input file containing a ParaMonte simulation specifications.<br>
    !>                                      <li>    A string containing the simulation specifications as if it contains the contents of
    !>                                              an external input simulation specification file in its entirety.
    !>                                  </ol>
    !>                                  The simulation specifications within the input file or string
    !>                                  must follow the rules of the Fortran language `namelist` IO.<br>
    !>                                  Four namelist grouping names are currently supported:<br>
    !>                                  <ol>
    !>                                      <li>    The grouping name `&ParaDISE` (case-<b>insensitive</b>) signifies
    !>                                              a specification grouping for the ParaDISE sampler and optimizer.<br>
    !>                                              When the user-specified `method` is `"ParaDISE"`,
    !>                                              this will be the default grouping name that will be search within the input file.<br>
    !>                                      <li>    The grouping name `&ParaDRAM` (case-<b>insensitive</b>) signifies
    !>                                              a specification grouping for the ParaDRAM sampler and integrator.<br>
    !>                                              When the user-specified `method` is `"ParaDRAM"`,
    !>                                              this will be the default grouping name that will be search within the input file.<br>
    !>                                      <li>    The grouping name `&ParaNest` (case-<b>insensitive</b>) signifies
    !>                                              a specification grouping for the ParaNest sampler and integrator.<br>
    !>                                              When the user-specified `method` is `"ParaNest"`,
    !>                                              this will be the default grouping name that will be search within the input file.<br>
    !>                                      <li>    The grouping name `ParaMonte` (case-<b>insensitive</b>) specifies a specification grouping
    !>                                              for any simulator/explorer supported by the ParaMonte library.<br>
    !>                                              This group naming is used as the last resort to read the simulation specifications
    !>                                              only if all of the above simulation-specific namelist groups are missing in the input file.<br>
    !>                                  </ol>
    !>  \param[inout]   err         :   The input/output scalar of type [err_type](@ref pm_err::err_type) containing
    !>                                  the relevant error information if any error occurs the IO within the routine.<br>
    !>                                  On output, if any error has occurred, an error message will be appended to the
    !>                                  `msg` component of `err` and the `occurred` component is set to `.true.`.<br>
    !>                                  Otherwise, the input contents of `err` remain intact on output.<br>
    !>
    !>  \interface{setSpecFromInput}
    !>  \code{.F90}
    !>
    !>      use pm_sampling_scio, only: setSpecFromInput
    !>
    !>      call setSpecFromInput(method, inputFile, err)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [pm_sampling](@ref pm_sampling)<br>
    !>  [pm_sampling_base](@ref pm_sampling_base)<br>
    !>  [pm_sampling_mcmc](@ref pm_sampling_mcmc)<br>
    !>  [pm_sampling_dram](@ref pm_sampling_dram)<br>
    !>  [pm_sampling_nest](@ref pm_sampling_nest)<br>
    !>
    !>  \final{setSpecFromInput}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 02:15 AM, September 27, 2021, Dallas, TX<br>
    subroutine setSpecFromInput(method, inputFile, err)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSpecFromInput
#endif
        character(*,SKG)    , intent(in)    :: inputFile
        character(*,SKG)    , intent(in)    :: method
        type(err_type)      , intent(inout) :: err
        logical(LK) :: opened
        integer(IK) :: unit
        integer :: iostat
        !logical(LK) :: missing
        !logical(LK) :: specIsNameList
        logical(LK) :: exist
        character(255, SK) :: iomsg
        character(:, SK), allocatable :: lcmethod
        character(:, SK), allocatable :: contents
        character(*, SK), parameter :: PROCEDURE_NAME = MODULE_NAME//SKG_"@setSpecFromInput()"
#define RETURN_IF_ERRED(LINE,FAILED,MSG) \
if (FAILED) then; \
err%msg = PROCEDURE_NAME//getFine(__FILE__, LINE)//SK_": "//trim(MSG); \
err%occurred = .true._LK; \
err%stat = iostat; \
return; \
end if;
        ! Determine if the file is external vs. internal.

        exist = .false.
        inquire(file = inputFile, exist = exist, iostat = iostat, opened = opened, number = unit)
        exist = exist .and. iostat == 0
        ! Intel compiler throws error for `inputFile` longer than 4096 characters.
        !return_if_erred(__LINE__, iostat /= 0, SK_": Failed to inquire the properties of the user-specified input file: '"//trim(inputFile)//SK_"'. "//trim(adjustl(iomsg)))

        ! Read the contents of the input file.

        if (exist) then
            ! Close the input file and read its contents.
            if (opened) then
                close(unit, iostat = iostat, iomsg = iomsg)
                RETURN_IF_ERRED(__LINE__, iostat /= 0, SK_"Failed to close the user-specified input file: '"//trim(inputFile)//SK_"'. "//trim(adjustl(iomsg)))
            end if
            call setContentsFrom(inputFile, contents, iostat, iomsg)
            !open(newunit = unit, file = inputFile, form = "formatted", access = "sequential", status = "old", action = "read", iostat = iostat, iomsg = iomsg SHARED)
            RETURN_IF_ERRED(__LINE__, iostat /= 0, SK_"Failed to open the user-specified input file: '"//trim(inputFile)//SK_"'. "//trim(adjustl(iomsg)))
        else
            contents = inputFile
        end if

        !iostat = index(getReplaced(contents, SK_" ", SK_""), SK_"&"//lcmethod, kind = IK)
        !RETURN_IF_ERRED(__LINE__, iostat == 0, SK_"Failed to find a `&"//lcmethod//SK_" namelist group in user-specified input file: '"//trim(inputFile)//SK_"'. ")

        ! Now write the contents to a temporary file for the processor to parse it.
        ! The rewrite is important to avoid gfortran end of file IO errors.

        open(newunit = unit, form = "formatted", access = "sequential", status = "scratch", action = "readwrite", iostat = iostat, iomsg = iomsg)
        RETURN_IF_ERRED(__LINE__, iostat /= 0, SK_"Failed to open a temporary input file to hold user-specified simulation specifications: "//trim(adjustl(iomsg)))
        write(unit, "(A)", iostat = iostat, iomsg = iomsg) contents//new_line(contents) ! The new line is essential for successful namelist content parsing by the GNU compiler.
        RETURN_IF_ERRED(__LINE__, iostat /= 0, SK_"Failed to write the user-specified simulation specifications to a temporary input file: "//trim(adjustl(iomsg)))
        rewind(unit, iostat = iostat, iomsg = iomsg)
        RETURN_IF_ERRED(__LINE__, iostat /= 0, SK_"Failed to rewind the temporary input file containing the user-specified simulation specifications: "//trim(adjustl(iomsg)))
        lcmethod = getStrLower(method)
        if (lcmethod == SK_"paradise") then
            read(unit, nml = paradise, iostat = iostat, iomsg = iomsg)
        elseif (lcmethod == SK_"paradram") then
            read(unit, nml = paradram, iostat = iostat, iomsg = iomsg)
        elseif (lcmethod == SK_"paranest") then
            read(unit, nml = paranest, iostat = iostat, iomsg = iomsg)
        else
            error stop PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SK_": Internal error occurred. Unrecognized input `method` = "//method
        end if
        close(unit, status = "delete", iostat = iostat)

        !if (iostat /= 0) then
        !    err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SK_": Failed to read the sampler-specific namelist group `&"//method//SK_"`: "//trim(adjustl(iomsg))
        !    ! Try the generic namelist group name: ParaMonte
        !    read(INPUT_FILE, nml = paramonte, iostat = iostat, iomsg = iomsg)
        !end if

        if (iostat /= 0) then
            err%occurred = .true._LK
            ! Check if the any namelist group exists at all.
            if (index(getStrLower(contents), lcmethod, kind = IK) == 0_IK) then! .and. index(contents, SK_"paramonte", kind = IK) == 0_IK) then
                err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SK_": Failed to detect any namelist groups named `&"//method//SK_"` in the specified input namelist file: "//trim(adjustl(iomsg))//SK_": """//contents//SK_""""
            else
                err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SK_": Failed to read the simulation specifications from the specified input string/file: "//trim(adjustl(iomsg))//SK_": """//contents//SK_""""
            end if
        end if

    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ParaDISE_ENABLED 1
#include "pm_sampling_scio.inc.F90"
#undef ParaDISE_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ParaDRAM_ENABLED 1
#include "pm_sampling_scio.inc.F90"
#undef ParaDRAM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ParaNest_ENABLED 1
#include "pm_sampling_scio.inc.F90"
#undef ParaNest_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! check if is end of record or file error condition.
    pure function iserf(stat) result(ended)
        integer, intent(in) :: stat
        logical :: ended
        ended = is_iostat_eor(stat) .or. is_iostat_end(stat)
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RETURN_IF_ERRED
#undef SHARED
