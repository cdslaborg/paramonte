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

        ! Define intel-specific shared property of files on Windows.
#if     INTEL_ENABLED && WINDOWS_ENABLED
#define SHARED, shared
#else
#define SHARED
#endif
    use pm_kind, only: SKG => SK, SK, IK, LK, RKD, modelr_type
   !use pm_io, only: FORMAT_GENERIC_BLANK_TABBED
   use pm_sysPath, only: MAX_LEN_FILE_PATH
    use pm_mathNumSys, only: getCountDigit
    use pm_paramonte, only: envis_type
    use pm_paramonte, only: envname, envis
    use pm_parallelism, only: image_type
    use pm_arrayRemove, only: getRemoved
    use pm_arrayResize, only: setResized
    use pm_matrixInit, only: getMatInit
    use pm_matrixInit, only: uppLowDia
    use pm_strASCII, only: getStrLower
    use pm_val2str, only: getStr
    use pm_except, only: setNAN
    use pm_except, only: isNAN
    use pm_except, only: isInf
    use pm_str, only: UNDEFINED
    use pm_io, only: field_type
    use pm_err, only: getFine
    use pm_io, only: filext
    use pm_io, only: TABEQV
    use pm_io, only: openArg_type
    use pm_distUnif, only: xoshiro256ssw_type
    use pm_io, only: display_type, mark_type, note_type, warn_type, stop_type, wrap_type
    use pm_strASCII, only: SUB ! Equivalent to the `Substitute` ASCII character used in the early days of communication.
    use pm_container, only: css_type
    use pm_timer, only: timer_type

    implicit none
#if OMP_ENABLED
    character(63,SKG)                       :: co_outputFileNameDef
#endif
    character(*,SKG)        , parameter     :: MODULE_NAME = SK_"@pm_sampling_base"
   !character(*,SKG)        , parameter     :: format = FORMAT_GENERIC_BLANK_TABBED ! format with which all simulation specs are reported to the output file.
    character(*,SKG)        , parameter     :: NL1 = new_line(SKG_"a")
    character(*,SKG)        , parameter     :: NL2 = NL1//NL1
   !character(*,SKG)        , parameter     :: RED = SKG_"ES" ! real edit descriptor.

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! simulation declarations.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type :: sfcbase_type
        integer(IK)                         :: nref = 0_IK              !<  \public The number of refinements, zero if sample size is prescribed by the user.
        integer(IK)         , allocatable   :: nsam(:)                  !<  \public The array of shape `(1:nref)` containing the number of unique samples (compact sample size) at each refinement stage.
        integer(IK)         , allocatable   :: sumw(:)                  !<  \public The array of shape `(1:nref)` containing the sum of the sample weights (verbose sample size) at each refinement stage.
        integer(IK)         , allocatable   :: sampleWeight(:)          !<  \public The array of shape `(1:nsam)` containing the final refined sample weights.
        real(RKG)           , allocatable   :: sampleLogFuncState(:,:)  !<  \public The array of shape `(1:ndim, nsam)` containing the final refined logFunc and sample weights.
        type(css_type)      , allocatable   :: colname(:)
        character(:,SKG)    , allocatable   :: header
    end type

    !> The Quantile derived type containing the distribution quantiles.
    type                                    :: quantile_type
        real(RKG)           , allocatable   :: quan(:,:) ! 9 by ndim.
        real(RKG)                           :: prob(9) = [0._RKG, 0.05_RKG, 0.10_RKG, 0.25_RKG, 0.50_RKG, 0.75_RKG, 0.90_RKG, 0.95_RKG, 1.0_RKG]
        character(4,SKG)                    :: name(9) = ["  Q0","  Q5"," Q10"," Q25"," Q50"," Q75"," Q90"," Q95","Q100"]
    end type

    type                                    :: mode_type
        integer(IK)                         :: loc
        real(RKG)                           :: val
        real(RKG)           , allocatable   :: crd(:)
    end type

    type                                    :: statistics_type
        integer(IK)                         :: size
        type(mode_type)                     :: mode
        type(quantile_type)                 :: quantile
        real(RKG)           , allocatable   :: cor(:,:)
        real(RKG)           , allocatable   :: cov(:,:)
        real(RKG)           , allocatable   :: avg(:)
    end type

    type                                    :: progress_type
        integer(IK)                         :: counterPRP = 0_IK                            !<  \public Counter for progress report file. The initialization is important.
        real(RKD)                           :: timeElapsedSinceStartInSeconds = 0._RKG      !<  \public Used in progress report. The initialization is important.
        real(RKD)                           :: clock = -1._RKG                              !<  \public Used in progress report. The initialization is important.
    end type

    type, abstract                          :: astatbase_type
        integer(IK)                         :: ess                                          !<  \public The effective sample size.
        class(timer_type)   , allocatable   :: timer                                        !<  \public timer.
        type(statistics_type)               :: chain                                        !<  \public chain statistics.
        type(statistics_type)               :: sample                                       !<  \public sample statistics.
        type(progress_type)                 :: progress                                     !<  \public sampling progress info.
        real(RKD)                           :: avgTimePerFunCall = 0._RKD                   !<  \public Average time per objective function call (in seconds).
        real(RKD)                           :: avgCommPerFunCall = 0._RKD                   !<  \public Average inter-process communication time per function call (in seconds).
        integer(IK)                         :: numFunCallAccepted = 0_IK                    !<  \public The number of accepted function calls during the sampling.
        integer(IK)                         :: numFunCallAcceptedRejected = 0_IK            !<  \public The number of accepted or rejected function calls during the sampling.
        integer(IK)                         :: numFunCallAcceptedLastReport = 0_IK          !<  \public Used in progress report.
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! auxiliary simulation information.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type                                    :: runis_type
        logical(LK)                         :: dry          !<  \public The scalar `logical` of default kind \LK, that is .true`. if and only if the simulation restarts from an existing incomplete simulation.
        logical(LK)                         :: new          !<  \public The scalar `logical` of default kind \LK, that is .true`. if and only if the simulation is fresh from the beginning.
    end type

    type                                    :: run_type
        integer(IK)                         :: id           !<  \public The scalar `integer` of default kind \IK, containing the current simulation run number (always starting at `1`).
        type(runis_type)                    :: is
    end type

    type                                    :: format_type
        character(:,SKG)    , allocatable   :: generic
        character(:,SKG)    , allocatable   :: integer
        character(:,SKG)    , allocatable   :: allreal
        character(:,SKG)    , allocatable   :: intreal
        character(:,SKG)    , allocatable   :: strreal
        character(:,SKG)    , allocatable   :: fixform
    end type

    type                                    :: tabularFileFormat_type
        character(:,SKG)    , allocatable   :: header
        character(:,SKG)    , allocatable   :: rows
    end type

    type, extends(openArg_type)             :: samplerFile_type
        logical(LK)                         :: extant = .false._LK
        character(:,SKG)    , allocatable   :: suffix
        character(:,SKG)    , allocatable   :: kind
        character(:,SKG)    , allocatable   :: ext
        type(css_type)      , allocatable   :: list(:)
    end type

    type, extends(samplerFile_type)         :: chainFile_type
        type(tabularFileFormat_type)        :: format
    end type

    type, extends(samplerFile_type)         :: progressFile_type
        type(tabularFileFormat_type)        :: format
    end type

    type, extends(samplerFile_type)         :: reportFile_type
        character(len(TABEQV, IK),SKG)      :: indent = TABEQV
        type(format_type)                   :: format
    end type

    type, extends(samplerFile_type)         :: restartFile_type
        integer(IK)                         :: counter = 0_IK
        character(:,SKG)    , allocatable   :: format
    end type

    type, extends(samplerFile_type)         :: sampleFile_type
        type(tabularFileFormat_type)        :: format
    end type

    type                                    :: ndim_type
        real(RKG)                           :: invhalf
        real(RKG)                           :: inv
        integer(IK)                         :: val
        integer(IK)                         :: vp1
        character(:,SKG)    , allocatable   :: str
        character(:,SKG)    , allocatable   :: desc
    end type

    type                                    :: method_type
        logical(LK)                         :: isParaDISE = .false._LK
        logical(LK)                         :: isParaDRAM = .false._LK
        logical(LK)                         :: isParaNest = .false._LK
        character(:,SKG)    , allocatable   :: val
        character(:,SKG)    , allocatable   :: desc
    end type

    type, extends(modelr_type)              :: real_type
        real(RKG)                           :: large = -huge(0._RKG)    !<  \public The largest working real number without undue Euclidean distance overflow.
        character(:,SKG)    , allocatable   :: desc
        character(:,SKG)    , allocatable   :: ed                       !<  \public real edit descriptor: ES
        character(:,SKG)    , allocatable   :: ex                       !<  \public real exponent descriptor: E0
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! specification declarations.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   !character(8192,SKG)                     :: description ! namelist input
    type                                    :: description_type
        character(:,SKG)    , allocatable   :: val
        character(:,SKG)    , allocatable   :: def
        character(:,SKG)    , allocatable   :: null
        character(:,SKG)    , allocatable   :: desc
    end type

   !character(:,SKG)        , allocatable   :: domain ! namelist input
    type                                    :: domain_type
        real(RKG)                           :: logVol = log(huge(0._RKG)) ! log-volume of the domain.
        logical(LK)                         :: isFinite ! paranest: always true. dram/dise: true only if domain is user-set.
        logical(LK)                         :: isCube
        logical(LK)                         :: isBall
        character(4,SKG)                    :: ball
        character(4,SKG)                    :: cube
        character(:,SKG)    , allocatable   :: val
        character(:,SKG)    , allocatable   :: def
        character(:,SKG)    , allocatable   :: null
        character(:,SKG)    , allocatable   :: desc
    end type

   !character(63,SKG)       , allocatable   :: domainAxisName(:) ! namelist input
    type                                    :: domainAxisName_type
        character(63,SKG)                   :: null
        character(63,SKG)   , allocatable   :: val(:)
        character(63,SKG)   , allocatable   :: def(:)
        character(:,SKG)    , allocatable   :: desc
        character(:,SKG)    , allocatable   :: prefix
    end type

   !real(RKG)               , allocatable   :: domainBallAvg(:) ! namelist input
    type                                    :: domainBallAvg_type
        real(RKG)           , allocatable   :: val(:)
        real(RKG)                           :: def
        character(:,SKG)    , allocatable   :: desc
    end type

   !real(RKG)               , allocatable   :: domainBallCor(:,:) ! namelist input
    type                                    :: domainBallCor_type
        real(RKG)           , allocatable   :: val(:,:)
        real(RKG)           , allocatable   :: def(:,:)
        character(:,SKG)    , allocatable   :: desc
    end type

   !real(RKG)               , allocatable   :: domainBallCov(:,:) ! namelist input
    type                                    :: domainBallCov_type
        logical(LK)                         :: isUserSet
        real(RKG)           , allocatable   :: def(:,:)
        real(RKG)           , allocatable   :: val(:,:)
        real(RKG)           , allocatable   :: inv(:,:)
        character(:,SKG)    , allocatable   :: desc
    end type

   !real(RKG)               , allocatable   :: domainBallStd(:) ! namelist input
    type                                    :: domainBallStd_type
        real(RKG)           , allocatable   :: val(:)
        real(RKG)           , allocatable   :: def
        character(:,SKG)    , allocatable   :: desc
    end type

   !real(RKG)               , allocatable   :: domainCubeLimitLower(:) ! namelist input
    type                                    :: domainCubeLimitLower_type
        real(RKG)                           :: def
        real(RKG)           , allocatable   :: val(:)
        character(:,SKG)    , allocatable   :: desc
    end type

   !real(RKG)               , allocatable   :: domainCubeLimitUpper(:) ! namelist input
    type                                    :: domainCubeLimitUpper_type
        real(RKG)                           :: def
        real(RKG)           , allocatable   :: val(:)
        character(:,SKG)    , allocatable   :: desc
    end type

   !integer(IK)                             :: domainErrCount ! namelist input
    type                                    :: domainErrCount_type
        integer(IK)                         :: def
        integer(IK)                         :: val
        integer(IK)                         :: null
        character(:,SKG)    , allocatable   :: str
        character(:,SKG)    , allocatable   :: desc
    end type

   !integer(IK)                             :: domainErrCountMax ! namelist input
    type                                    :: domainErrCountMax_type
        integer(IK)                         :: val
        integer(IK)                         :: def
        integer(IK)                         :: null
        character(:,SKG)    , allocatable   :: str
        character(:,SKG)    , allocatable   :: desc
    end type

   !logical(LK)                             :: inputFileHasPriority
    type                                    :: inputFileHasPriority_type
        logical(LK)                         :: val
        logical(LK)                         :: def
        character(:,SKG)    , allocatable   :: desc
    end type

   !character(:,SKG)        , allocatable   :: outputChainFileFormat
    type                                    :: outputChainFileFormat_type
        logical(LK)                         :: isCompact
        logical(LK)                         :: isVerbose
        logical(LK)                         :: isBinary
        character(7,SKG)                    :: compact
        character(7,SKG)                    :: verbose
        character(6,SKG)                    :: binary
        character(:,SKG)    , allocatable   :: def
        character(:,SKG)    , allocatable   :: val
        character(:,SKG)    , allocatable   :: null
        character(:,SKG)    , allocatable   :: desc
    end type

   !integer(IK)                             :: outputColumnWidth ! namelist input
    type                                    :: outputColumnWidth_type
        integer(IK)                         :: null
        integer(IK)                         :: def
        integer(IK)                         :: val
        character(:,SKG)    , allocatable   :: str
        character(:,SKG)    , allocatable   :: max
        character(:,SKG)    , allocatable   :: desc
    end type

   !character(:,SKG)        , allocatable   :: outputFileName ! namelist input
    type                                    :: outputFileName_type
        character(:,SKG)    , allocatable   :: full
        character(:,SKG)    , allocatable   :: base
        character(:,SKG)    , allocatable   :: seps
        character(:,SKG)    , allocatable   :: sep
        character(:,SKG)    , allocatable   :: dir
        character(:,SKG)    , allocatable   :: val
        character(:,SKG)    , allocatable   :: def
        character(:,SKG)    , allocatable   :: null
        character(:,SKG)    , allocatable   :: desc
    end type

   !character(:,SKG)        , allocatable   :: outputStatus
   type                                     :: outputStatusIs_type
        logical(LK)                         :: extend
        logical(LK)                         :: repeat
        logical(LK)                         :: retry
   end type
    type                                    :: outputStatus_type
        type(outputStatusIs_type)           :: is
        character(:,SKG)    , allocatable   :: def
        character(:,SKG)    , allocatable   :: val
        character(:,SKG)    , allocatable   :: null
        character(:,SKG)    , allocatable   :: desc
    end type

   !integer(IK)                             :: outputPrecision ! namelist input
    type                                    :: outputPrecision_type
        character(:,SKG)    , allocatable   :: str
        integer(IK)                         :: val
        integer(IK)                         :: def
        integer(IK)                         :: null
        character(:,SKG)    , allocatable   :: desc
    end type

   !integer(IK)                             :: outputReportPeriod ! namelist input
    type                                    :: outputReportPeriod_type
        real(RKG)                           :: inv
        integer(IK)                         :: val
        integer(IK)                         :: def
        integer(IK)                         :: null
        character(:,SKG)    , allocatable   :: desc
    end type

   !character(:,SKG)        , allocatable   :: outputRestartFileFormat ! namelist input
    type                                    :: outputRestartFileFormat_type
        logical(LK)                         :: isBinary
        logical(LK)                         :: isAscii
        character(6,SKG)                    :: binary
        character(5,SKG)                    :: ascii
        character(:,SKG)    , allocatable   :: def
        character(:,SKG)    , allocatable   :: val
        character(:,SKG)    , allocatable   :: null
        character(:,SKG)    , allocatable   :: desc
    end type

   !integer(IK)                             :: outputSampleSize ! namelist input
    type                                    :: outputSampleSize_type
        character(:,SKG)    , allocatable   :: str
        integer(IK)                         :: val
        integer(IK)                         :: def
        integer(IK)                         :: null
        character(:,SKG)    , allocatable   :: desc
    end type

   !character(:,SKG)        , allocatable   :: outputSeparator ! namelist input
    type                                    :: outputSeparator_type
        character(:,SKG)    , allocatable   :: val
        character(:,SKG)    , allocatable   :: def
        character(:,SKG)    , allocatable   :: null
        character(:,SKG)    , allocatable   :: desc
    end type

   !character(:,SKG)        , allocatable   :: outputSplashMode ! namelist input
    type                                    :: outputSplashModeIs_type
        logical(LK)                         :: normal = .false._LK
        logical(LK)                         :: silent = .false._LK
        logical(LK)                         :: quiet = .false._LK
    end type
    type                                    :: outputSplashMode_type
        type(outputSplashModeIs_type)       :: is
        character(6,SKG)                    :: normal = SKG_"normal"
        character(6,SKG)                    :: silent = SKG_"silent"
        character(5,SKG)                    :: quiet = SKG_"quiet"
        character(:,SKG)    , allocatable   :: null
        character(:,SKG)    , allocatable   :: def
        character(:,SKG)    , allocatable   :: val
        character(:,SKG)    , allocatable   :: desc
    end type

   !character(:,SKG)        , allocatable   :: parallelism
    type                                    :: parallelismis_type
        logical(LK)                         :: forkJoin = .false._LK
        logical(LK)                         :: multiChain = .false._LK
        logical(LK)                         :: singleChain = .false._LK
    end type
    type                                    :: parallelism_type
        type(parallelismis_type)            :: is
        character(11,SKG)                   :: singleChain
        character(10,SKG)                   :: multiChain
        character(:,SKG)    , allocatable   :: def
        character(:,SKG)    , allocatable   :: val
        character(:,SKG)    , allocatable   :: null
        character(:,SKG)    , allocatable   :: desc
    end type

   !logical(LK)                             :: parallelismMpiFinalizeEnabled ! namelist input
    type                                    :: parallelismMpiFinalizeEnabled_type
        logical(LK)                         :: val
        logical(LK)                         :: def
        character(:,SKG)    , allocatable   :: desc
    end type

    type                                    :: parallelismNumThread_type
        integer(IK)                         :: val
        integer(IK)                         :: def
        integer(IK)                         :: null
        character(:,SKG)    , allocatable   :: desc
    end type

   !character(:,SKG)        , allocatable   :: plang ! namelist input
    type                                    :: plang_type ! programming language environment.
        character(:,SKG)    , allocatable   :: val
        character(:,SKG)    , allocatable   :: def
        character(:,SKG)    , allocatable   :: null
       !character(:,SKG)    , allocatable   :: desc
       !type(envis_type)                    :: is = envis
    end type

   !integer(IK)                             :: randomSeed ! namelist input
    type                                    :: randomSeed_type
        integer(IK)                         :: val
        integer(IK)                         :: null
       !integer(IK)                         :: imageID
       !integer(IK)                         :: sizeSeed
       !integer(IK)                         :: imageCount
       !integer(IK)                         :: processID
       !integer(IK)         , allocatable   :: Seed(:,:)
        character(:,SKG)    , allocatable   :: desc
    end type

   !character(:,SKG)        , allocatable   :: sysInfoFilePath ! namelist input
    type                                    :: sysInfoFilePath_type
        character(:,SKG)    , allocatable   :: def
        character(:,SKG)    , allocatable   :: val
        character(:,SKG)    , allocatable   :: null
    end type

   !real(RKG)                               :: targetAcceptanceRate(2) ! namelist input
    type                                    :: targetAcceptanceRate_type
        logical(LK)                         :: enabled
        real(RKG)                           :: val(2)
        real(RKG)                           :: def(2)
        real(RKG)                           :: aim ! The target efficiency (depends on the sampler).
        character(:,SKG)    , allocatable   :: desc
    end type

    !>  \brief
    !>  This is the derived type containing the base properties of a ParaMonte sampler.
    !>
    !>  \details
    !>  This derived type contains (and shall contain) only the base sampler specifications
    !>  that are meant to be set and frozen before the start of the simulation kernel.<br>
    !>  Any component that changes over the course of simulation should become
    !>  part of [statbase_type](@ref pm_sampling_base::statbase_type).<br>
    !>
    type :: specbase_type
        type(xoshiro256ssw_type)                    :: rng                      !<  \public Random number generator.
        type(run_type)                              :: run                      !<  \public The simulation run ID and information.
        logical(LK)                                 :: overridable              !<  \public The scalar `logical` of default kind \LK that is `.true.` if and only if the namelist argument `inputFileHasPriority` is `.false.`.
        character(:,SKG), allocatable               :: msg                      !<  \public The scalar `character` containing the messages common between stdout and the report file.
        type(ndim_type)                             :: ndim
        type(image_type)                            :: image
        type(method_type)                           :: method
        type(real_type)                             :: real
        type(display_type)                          :: disp
        type(chainFile_type)                        :: chainFile
        type(progressFile_type)                     :: progressFile
        type(reportFile_type)                       :: reportFile
        type(restartFile_type)                      :: restartFile
        type(sampleFile_type)                       :: sampleFile
        type(description_type                   )   :: description
        type(domain_type                        )   :: domain
        type(domainAxisName_type                )   :: domainAxisName
        type(domainBallAvg_type                 )   :: domainBallAvg
        type(domainBallCor_type                 )   :: domainBallCor
        type(domainBallCov_type                 )   :: domainBallCov
        type(domainBallStd_type                 )   :: domainBallStd
        type(domainCubeLimitLower_type          )   :: domainCubeLimitLower
        type(domainCubeLimitUpper_type          )   :: domainCubeLimitUpper
        type(domainErrCount_type                )   :: domainErrCount
        type(domainErrCountMax_type             )   :: domainErrCountMax
        type(inputFileHasPriority_type          )   :: inputFileHasPriority
        type(outputChainFileFormat_type         )   :: outputChainFileFormat
        type(outputColumnWidth_type             )   :: outputColumnWidth
        type(outputFileName_type                )   :: outputFileName
        type(outputPrecision_type               )   :: outputPrecision
        type(outputReportPeriod_type            )   :: outputReportPeriod
        type(outputRestartFileFormat_type       )   :: outputRestartFileFormat
        type(outputSampleSize_type              )   :: outputSampleSize
        type(outputSeparator_type               )   :: outputSeparator
        type(outputSplashMode_type              )   :: outputSplashMode
        type(outputStatus_type                  )   :: outputStatus
        type(parallelism_type                   )   :: parallelism
        type(parallelismMpiFinalizeEnabled_type )   :: parallelismMpiFinalizeEnabled
        type(parallelismNumThread_type          )   :: parallelismNumThread
        type(plang_type                         )   :: plang
        type(randomSeed_type                    )   :: randomSeed
        type(sysInfoFilePath_type               )   :: sysInfoFilePath
        type(targetAcceptanceRate_type          )   :: targetAcceptanceRate
    contains
        procedure, pass, public                     :: set
        procedure, pass, private                    :: isFailedResizeList
        procedure, pass, private                    :: openFiles
        procedure, pass, private                    :: openFile
        procedure, pass, private                    :: sanitize
        procedure, pass, private                    :: report
    end type

    interface specbase_type
        module procedure :: specbase_typer
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine killMeAlreadyCMake1_RK5(); use pm_sampling_scio_RK5, only: RKG; end subroutine
    subroutine killMeAlreadyCMake1_RK4(); use pm_sampling_scio_RK4, only: RKG; end subroutine
    subroutine killMeAlreadyCMake1_RK3(); use pm_sampling_scio_RK3, only: RKG; end subroutine
    subroutine killMeAlreadyCMake1_RK2(); use pm_sampling_scio_RK2, only: RKG; end subroutine
    subroutine killMeAlreadyCMake1_RK1(); use pm_sampling_scio_RK1, only: RKG; end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function specbase_typer(modelr, method, ndim) result(spec)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: specbase_typer
#endif
        type(modelr_type), intent(in) :: modelr
        character(*,SKG), intent(in) :: method
        integer(IK), intent(in) :: ndim
        type(specbase_type) :: spec

        block
            integer(IK) :: width
            character(:,SKG), allocatable :: prefix
            prefix = spec%reportFile%indent//method
            width = 115_IK - len(prefix, IK)
            spec%disp = display_type( tmsize = 0_IK, bmsize = 1_IK & ! LCOV_EXCL_LINE
                                    , mark = mark_type(prefix = spec%reportFile%indent, tmsize = 0_IK, bmsize = 1_IK, width = 99999_IK) & ! LCOV_EXCL_LINE
                                    , note = note_type(prefix = prefix, tmsize = 0_IK, bmsize = 1_IK, width = width) & ! LCOV_EXCL_LINE
                                    , warn = warn_type(prefix = prefix, tmsize = 0_IK, bmsize = 1_IK, width = width) & ! LCOV_EXCL_LINE
                                    , stop = stop_type(prefix = prefix, tmsize = 1_IK, bmsize = 2_IK, width = width) & ! LCOV_EXCL_LINE
                                    )
        end block
        spec%disp%width = max(80_IK, min(132_IK, abs(spec%disp%width))) - 4_IK ! laptops have tiny screens - 2 * margin width (= 2).
        spec%disp%text%width = spec%disp%width
        spec%real%large = real(aint(sqrt(modelr%huge) / ndim), RKG)
        spec%image = image_type()

        ndim_block: block
            spec%ndim%val = ndim
            spec%ndim%vp1 = ndim + 1_IK
            spec%ndim%inv = 1._RKG / real(ndim, RKG)
            spec%ndim%invhalf = .5_RKG * spec%ndim%inv
           !spec%ndim%sqpndim = ndim**2 + ndim
            spec%ndim%str = getStr(ndim)
            spec%ndim%desc = &
            SKG_"The simulation specification `ndim` is a positive scalar of type `integer` of kind 32-bit, &
                &representing the number of dimensions of the domain of the objective function. &
                &Unlike most other simulation specifications, `ndim` cannot be set from within an external input function &
                &and the user must always specify it along with the objective function at the time of calling the exploration routine &
                &separately from the rest of the (optional) simulation specifications. &
                &As such, `ndim` has no default value and is always user-supplied."
        end block ndim_block

        method_block: block
            character(:,SKG), allocatable :: lowval
            spec%method%val = method
            lowval = getStrLower(spec%method%val)
            spec%method%isParaDISE = logical(lowval == SKG_"paradise", LK)
            spec%method%isParaDRAM = logical(lowval == SKG_"paradram", LK)
            spec%method%isParaNest = logical(lowval == SKG_"paranest", LK)
            spec%method%desc = &
            SKG_"The simulation specification `method` is a scalar string, &
                &representing the method of exploring the objective function. &
                &Unlike most other simulation specifications, `method` cannot be set from within an external input function &
                &and the user always specifies it by calling the corresponding objective function exploration routine. &
                &As such, `method` has no default value and is always user-supplied. &
                &The following density function exploration methods are currently or will be supported by the ParaMonte library:"//NL2//&
            SKG_"   ParaDRAM"//NL2//&
            SKG_"       The ParaDRAM acronym stands for Parallel Delayed-Rejection Adaptive Metropolis Markov Chain Monte Carlo. &
                        &The ParaDRAM sampler can optimize and explore the landscape of objective functions stochastically. &
                        &For more details of the ParaDRAM sampling scheme, refer to the references section in the documentation of the ParaMonte library:"//NL2//&
            SKG_"           https://www.cdslab.org/paramonte"
        end block method_block

        modelr_block: block
            spec%real%modelr_type = modelr
            spec%real%ed = SKG_"ES"
            spec%real%ex = SKG_"E0"
            spec%real%desc = &
            SKG_"The simulation specification `modelr` is an object whose components contain the &
                &characteristics of `real` (float) data type model (hence the name `modelr`) used for (most) computations during the simulation. &
                &The real/float type used within the simulation is the same as that of the input and output values to and from the objective function. &
                &It is set by the user when the sampler is called. As such, unlike most other simulation specifications, `modelr` cannot be set from &
                &within an external input function, the user always specifies it when calling the corresponding sampler. &
                &Therefore, `modelr` has no default value and is always user-supplied when calling the sampler. &
                &The components of `modelr` are important for future determination of the type and kind of the real/float used in the simulation. &
                &The `modelr` object contains the following real/float characteristics:"//NL2//&
            SKG_"+   `digits`"//NL2//&
            SKG_"    The number of significant bit digits in the real/float model used in the simulations."//NL2//&
            SKG_"+   `kind`"//NL2//&
            SKG_"    The processor-dependent kind type-parameter of the real/float model used in the current simulation."//NL2//&
            SKG_"+   `maxexponent`"//NL2//&
            SKG_"    The maximum exponent for the real/float model used in the current simulation."//NL2//&
            SKG_"+   `minexponent`"//NL2//&
            SKG_"    The minimum exponent for the real/float model used in the current simulation."//NL2//&
            SKG_"+   `precision`"//NL2//&
            SKG_"    The equivalent decimal precision (number of digits after the decimal point) of the real/float model used in the current simulation."//NL2//&
            SKG_"+   `radix`"//NL2//&
            SKG_"    The base in the real/float model used in the current simulation."//NL2//&
            SKG_"+   `range`"//NL2//&
            SKG_"    The decimal exponent range in the real/float model used in the current simulation. &
                     &It is the value corresponding to `int(min(log10(huge), -log10(tiny)))` where `huge` &
                     &and `tiny` are the largest and smallest positive numbers in the real/float model"//NL2//&
            SKG_"+   `storage_size`"//NL2//&
            SKG_"    The size, in bits, that would be taken in memory by a scalar of the real/float model used in the current simulation."//NL2//&
            SKG_"+   `epsilon`"//NL2//&
            SKG_"    The smallest normal value representable by the real/float model used in the current simulation such that `1 < 1 + epsilon`."//NL2//&
            SKG_"+   `huge`"//NL2//&
            SKG_"    The largest normal value representable by the real/float model used in the current simulation (without causing overflow)."//NL2//&
            SKG_"+   `tiny`"//NL2//&
            SKG_"    The smallest positive normal value representable by the real/float model used in the current simulation (without causing underflow)."
        end block modelr_block

        description_block: block
            use pm_sampling_scio, only: description
            spec%description%null = repeat(SUB, len(description, IK))
            spec%description%def = UNDEFINED
            spec%description%desc = &
            SKG_"The simulation specification `description` is a scalar string of maximum length `"//getStr(len(description, IK))//&
            SKG_"` characters containing general information about the simulation to be performed. It has no effects &
            &on the simulation and serves only as a general description for future reference. &
            &The sampler parser automatically recognizes the C-style `'\n'` escape sequence as the new-line character &
            &and `'\\'` as the backslash character `'\'` if used in the contents of `description`. For example, `'\\n'` &
            &will be converted to `'\n'` on the output, while `'\n'` translates to the new-line character. &
            &The default value for `description` is `'"//spec%description%def//SKG_"'`."
            !$omp master
            description = spec%description%null
            !$omp end master
        end block description_block

        domain_block: block
            use pm_sampling_scio, only: domain
            spec%domain%isFinite = spec%method%isParaNest
            spec%domain%isCube = .false._LK
            spec%domain%isBall = .false._LK
            spec%domain%cube   = SKG_"cube"
            spec%domain%ball   = SKG_"ball"
            spec%domain%null   = repeat(SUB, len(domain, IK))
            spec%domain%def    = spec%domain%cube
            spec%domain%desc   = &
            SKG_"The simulation specification `domain` is a scalar string of maximum length `"//getStr(len(domain, IK))//SKG_"`, &
                &containing the model's name that defines the objective function's domain. &
                &The string value must be enclosed by either single or double quotation marks when provided as input. &
                &The following domain models are currently supported:"//NL2//&
            SKG_"+   `domain = '"//spec%domain%cube//SKG_"'`"//NL2//&
            SKG_"    This is equivalent to a `ndim`-dimensional hyper-cube (n-cube) whose upper and lower bounds are specified &
                     &by the input simulation specifications `domainCubeLimitUpper` and `domainCubeLimitLower` respectively."//NL2//&
            SKG_"+   `domain = '"//spec%domain%ball//SKG_"'`"//NL2//&
            SKG_"    This is equivalent to a `ndim`-dimensional hyper-ellipsoid (n-ball) whose center and covariance matrix can be &
                     &specified by the input simulation specification `domainBallAvg` and `domainBallCov` respectively. Alternatively, &
                     &the user can let the ParaMonte samplers construct the covariance matrix of the ellipsoidal domain from the input &
                     &values for the `domainBallCor` and `domainBallStd` simulation specifications. Note that a spherical &
                     &domain can be defined by dropping the `domainBallCov` and `domainBallCor` specifications from the input &
                     &and setting all elements of `domainBallStd` to the desired radius of the domain."//NL2//&
            SKG_"The default value for `domain` is an infinite cube for the ParaDRAM and ParaDISE samplers, and a unit-sized cube for the ParaNest sampler."
            !$omp master
            domain = spec%domain%null
            !$omp end master
        end block domain_block

        domainAxisName_block: block
            use pm_sampling_scio, only: domainAxisName
            integer(IK) :: idim
            spec%domainAxisName%prefix = SKG_"sampleState"
            spec%domainAxisName%null = repeat(SUB, len(domainAxisName, IK))
            call setResized(spec%domainAxisName%def, ndim)
            do idim = 1, ndim
                spec%domainAxisName%def(idim) = spec%domainAxisName%prefix//getStr(idim)
            end do
            spec%domainAxisName%desc = &
            SKG_"The simulation specification `domainAxisName` is a vector of scalar string values, each element of which contains &
                &the names of the corresponding axis of the density function domain (to be explored/sampled). &
                &It is used to construct the headers of the simulation output files. &
                &Any element of `domainAxisName` not set by the user will automatically be assigned a default name. &
                &If all elements of `domainAxisName` are set to the same value, then a number will be suffixed to each element &
                &representing the ID of the corresponding dimension of the domain of the density function. &
                &The default value for `domainAxisName` is '"//spec%domainAxisName%prefix//SKG_"i', where the integer `'i'` at the &
                &end of the name is replaced by the index of the corresponding domain axis."
            !$omp master
            call setResized(domainAxisName, ndim)
            domainAxisName = spec%domainAxisName%null
            !$omp end master
        end block domainAxisName_block

        domainBallAvg_block: block
            use pm_sampling_scio, only: domainBallAvg
            spec%domainBallAvg%def = 0._RKG
            spec%domainBallAvg%desc = &
            SKG_"The simulation specification `domainBallAvg` is a vector of type `real` of the highest precision supported by the ParaMonte library, &
                &of size `ndim` containing the coordinates of the center of the hyper-ellipsoidal (or spherical) `ndim`-dimensional domain of the objective function, &
                &such that all states `X(1 : ndim)` visited by the sampler would obey the following inequality:"//NL2//&
            SKG_"$$"//NL1//&
            SKG_"(X - \mu)^T \Sigma^{-1} (X - \mu) \leq 1."//NL1//&
            SKG_"$$"//NL2//&
            SKG_"where $$\mu$$ represents the specification `domainBallAvg` and $$\Sigma^{-1}$$ is the inverse of the specification `domainBallCov`."//NL2//&
            SKG_"When passed to the sampler from within an external input file, every missing element of `domainBallAvg` will be set to the origin (zero). &
                &Together with `domainBallCov`, or with `domainBallCor` and `domainBallStd`, it forms a hyper-ellipsoidal &
                &or hyper-spherical domain for the ParaMonte samplers. Note that an ellipsoidal/spherical domain is used if and only if the &
                &input simulation specification `domain` is set to 'ellipsoid' or 'sphere' or `ball`. Otherwise, a cubical domain will be used. &
                &The key `Avg` (standing for `Average`) is used in the name of this variable to denote the center of the hyper-ellipsoid. &
                &The default value for `domainBallAvg` is the origin (i.e., a zero-valued vector of size `ndim`)."
            !$omp master
            call setResized(domainBallAvg, ndim)
            call setNAN(domainBallAvg(1 : ndim))
            !$omp end master
        end block domainBallAvg_block

        domainBallCor_block: block
            use pm_sampling_scio, only: domainBallCor
            spec%domainBallCor%def = getMatInit([ndim, ndim], uppLowDia, vupp = 0._RKG, vlow = 0._RKG, vdia = 1._RKG)
            spec%domainBallCor%desc = &
            SKG_"The simulation specification `domainBallCor` is a positive-definite matrix of type `real` of the highest precision available in the ParaMonte Library &
                &of size `(ndim, ndim)` representing the correlation matrix of the domain of the objective function, where `ndim` is the dimension of the domain, &
                &such that all states `X(1 : ndim)` visited by the sampler would obey the following inequality:"//NL2//&
            SKG_"$$"//NL1//&
            SKG_"(X - \mu)^T \Sigma^{-1} (X - \mu) \leq 1."//NL1//&
            SKG_"$$"//NL2//&
            SKG_"where $$\mu$$ represents the specification `domainBallAvg` and $$\Sigma^{-1}$$ is the inverse of $$\Sigma$$ computed as, "//NL2//&
            SKG_"$$"//NL1//&
            SKG_"\Sigma = \mathrm{eye}(V) ~ \rho ~ \mathrm{eye}(V)."//NL1//&
            SKG_"$$"//NL2//&
            SKG_"where $$\rho$$ stands for specification `domainBallCor` and $$\mathrm{eye}(V)$$ stands &
                &for diagonal matrix whose diagonals are set to the specification `domainBallStd`."//NL2//&
            SKG_"Combined with the input simulation specification `domainBallStd` it defines the objective function's hyper-ellipsoidal (or spherical) domain. &
                &If the input simulation specification `domainBallCov` is provided by the user, then any values set &
                &for `domainBallCor` and `domainBallStd` will be automatically ignored. The input specification `domainBallCor` &
                &along with `domainBallStd` are especially useful when covariance matrix computation is non-trivial. &
                &When passed to the sampler from within an external input sampler specification file, any missing element of `domainBallCov` &
                &will be set to the appropriate default value. The default value for `domainBallCor` is a `ndim`-by-`ndim` Identity matrix."
            !$omp master
            call setResized(domainBallCor, [ndim, ndim])
            call setNAN(domainBallCor)
            !$omp end master
        end block domainBallCor_block

        domainBallCov_block: block
            use pm_sampling_scio, only: domainBallCov
            spec%domainBallCov%isUserSet = .false._LK ! This must be set here. It is important for the proper setting from inputFile and inputArg.
            if (spec%method%isParaDRAM .or. spec%method%isParaDISE) then
                spec%domainBallCov%def = getMatInit([ndim, ndim], uppLowDia, 0._RKG, 0._RKG, spec%real%large**2)
            elseif (spec%method%isParaNest) then
                spec%domainBallCov%def = getMatInit([ndim, ndim], uppLowDia, 0._RKG, 0._RKG, 1._RKG)
            end if
            call setResized(spec%domainBallCov%def, [ndim, ndim])
            call setResized(spec%domainBallCov%val, [ndim, ndim])
            spec%domainBallCov%desc = &
            SKG_"The simulation specification `domainBallCov` is a positive-definite matrix of type `real` of the highest precision available in the ParaMonte Library &
                &of size `(ndim, ndim)` representing the Gramian matrix of the domain of the objective function, where `ndim` is the dimension of the domain, &
                &such that all states `X(1 : ndim)` visited by the sampler would obey the following inequality:"//NL2//&
            SKG_"$$"//NL1//&
            SKG_"(X - \mu)^T \Sigma^{-1} (X - \mu) \leq 1."//NL1//&
            SKG_"$$"//NL2//&
            SKG_"where $$\mu$$ represents the specification `domainBallAvg` and $$\Sigma^{-1}$$ is the inverse of the specification `domainBallCov`."//NL2//&
            SKG_"If the user provides this input simulation specification, any values set for the input simulation &
                &specifications `domainBallCor` and `domainBallStd` will be automatically ignored. When set from inside an external &
                &input ParaMonte specification file, any missing element of `domainBallCov` will be set to the appropriate default value. &
                &To specify a `ndim`-dimensional spherical domain, set `domainBallCov` to the identity matrix whose diagonal elements are &
                &radius-squared of the desired hyper-sphere (n-ball). The default value for `domainBallCov` is a `ndim`-by-`ndim` Identity matrix &
                &for simulations (such as the ParaNest integrator) that require a finite domain and a `ndim`-by-`ndim` diagonal matrix whose diagonals &
                &are practically set to infinity for simulations that do not require a finite domain (such as the ParaDRAM and ParaDISE MCMC samplers)."//NL2//&
            SKG_"> **Note**: The use of `Cov` in the name of this simulation specification is theoretically incorrect as the domain &
                &of the objective function is not a distribution. Even if it is considered a hyper-ellipsoidal uniform distribution &
                &this specification would still not represent its covariance matrix because it represents the Gramian matrix. &
                &However, the decision was made to name this specification as `Cov` because of its nice fit to the rest of &
                &the relevant simulation specifications and user familiarity with keywords."
            !$omp master
            call setResized(domainBallCov, [ndim, ndim])
            call setNAN(domainBallCov)
            !$omp end master
        end block domainBallCov_block

        domainBallStd_block: block
            use pm_sampling_scio, only: domainBallStd
            if (spec%method%isParaDRAM .or. spec%method%isParaDISE) then
                spec%domainBallStd%def = spec%real%large
            elseif (spec%method%isParaNest) then
                spec%domainBallStd%def = 1._RKG
            end if
            spec%domainBallStd%desc = &
            SKG_"The simulation specification `domainBallStd` is a positive-valued vector of type `real` of the highest precision available &
                &within the ParaMonte library, of size `ndim`, where `ndim` is the dimension of the domain of the objective function. &
                &It represents the square roots of the diagonal elements of the covariance matrix of the domain of the objective function, &
                &such that all states `X(1 : ndim)` visited by the sampler would obey the following inequality:"//NL2//&
            SKG_"$$"//NL1//&
            SKG_"(X - \mu)^T \Sigma^{-1} (X - \mu) \leq 1."//NL1//&
            SKG_"$$"//NL2//&
            SKG_"where $$\mu$$ represents the specification `domainBallAvg` and $$\Sigma^{-1}$$ is the inverse of $$\Sigma$$ computed as, "//NL2//&
            SKG_"$$"//NL1//&
            SKG_"\Sigma = \mathrm{eye}(V) ~ \rho ~ \mathrm{eye}(V)."//NL1//&
            SKG_"$$"//NL2//&
            SKG_"where $$\rho$$ stands for specification `domainBallCor` and $$\mathrm{eye}(V)$$ stands &
                &for diagonal matrix whose diagonals are set to the specification `domainBallStd`."//NL2//&
            SKG_"If the covariance matrix of the ellipsoidal/spherical domain (`domainBallCov`) is missing from the input specifications &
                &to the sampler, then `domainBallStd` (along with the input specification `domainBallCor`) will be used to &
                &construct the covariance matrix of the domain of the objective function. However, if `domainBallCov` is present &
                &among the input specifications to the sampler, then the input values for `domainBallStd` and `domainBallCor` &
                &will be ignored, and `domainBallCov` will be used to construct the domain of the user-specified objective function. &
                &To specify a `ndim`-dimensional spherical domain, drop `domainBallCov` and `domainBallCor` from the input and &
                &set all elements of `domainBallStd` to the desired radius of the hyper-spherical domain. The default value for &
                &any missing elements of `domainBallStd` is `1` for simulations requiring a finite domain (such as the ParaNest &
                &integrator) and `+Infinity` for simulations not requiring a finite domain (such as the ParaDRAM and ParaDISE samplers)."//NL2//&
            SKG_"> **Note**: The use of `Std` in the name of this simulation specification is theoretically incorrect as the domain &
                &of the objective function is not a distribution. Even if it is considered a hyper-ellipsoidal uniform distribution, &
                &this specification would still not represent its vector of standard deviations. &
                &However, the decision was made to name this specification as `Std` because of its nice fit to the rest of &
                &the relevant simulation specifications and user familiarity with keywords."
            !$omp master
            call setResized(domainBallStd, ndim)
            call setNAN(domainBallStd)
            !$omp end master
        end block domainBallStd_block

        domainCubeLimitLower_block: block
            use pm_sampling_scio, only: domainCubeLimitLower
            if (spec%method%isParaDRAM .or. spec%method%isParaDISE) then
                spec%domainCubeLimitLower%def = -spec%real%large
            elseif (spec%method%isParaNest) then
                spec%domainCubeLimitLower%def = -1._RKG
            end if
            spec%domainCubeLimitLower%desc = &
            SKG_"The simulation specification `domainCubeLimitLower` is a vector of type `real` of the highest precision available within &
                &the ParaMonte library, of size `ndim` where `ndim` is the number of dimensions of the domain of the target density function. &
                &It contains the lower boundaries of the cubical domain of the objective function to be sampled. &
                &When `domainCubeLimitLower` is specified inside an external input file supplied to the sampler, &
                &it is also possible to assign only select values of `domainCubeLimitLower` &
                &and leave the rest of the components to be assigned the default value. &
                &For example, having the following inside the input file, "//NL2//&
            SKG_"+   `domainCubeLimitLower(3:5) = -100`"//NL2//&
            SKG_"    will only set the lower limits of the third, fourth, and fifth dimensions to `-100`, or,"//NL2//&
            SKG_"+   `domainCubeLimitLower(1) = -100, domainCubeLimitLower(2) = -1.e6`"//NL2//&
            SKG_"    will set the lower limit on the first dimension to `-100`, and `1.e6` on the second dimension, or,"//NL2//&
            SKG_"+   `domainCubeLimitLower = 3*-2.5e100`"//NL2//&
            SKG_"    will only set the lower limits on the first, second, and third dimensions to `-2.5*10^100`, while the &
                     &remaining lower limits for the missing dimensions will be automatically set to the default value."//NL2//&
            SKG_"The default value for all elements of `domainCubeLimitLower` is `"//getStr(spec%domainCubeLimitLower%def)//SKG_"`."//NL1//&
            SKG_"Beware that some ParaMonte samplers, such as ParaNest, require the user to specify the domain boundaries explicitly."
            !$omp master
            call setResized(domainCubeLimitLower, ndim)
            call setNAN(domainCubeLimitLower)
            !$omp end master
        end block domainCubeLimitLower_block

        domainCubeLimitUpper_block: block
            use pm_sampling_scio, only: domainCubeLimitUpper
            if (spec%method%isParaDRAM .or. spec%method%isParaDISE) then
                spec%domainCubeLimitUpper%def = +spec%real%large
            elseif (spec%method%isParaNest) then
                spec%domainCubeLimitUpper%def = +1._RKG
            end if
            spec%domainCubeLimitUpper%desc = &
            SKG_"The simulation specification `domainCubeLimitUpper` is a vector of type `real` of the highest precision available within &
                &by the ParaMonte library, of size `ndim` where `ndim` is the number of dimensions of the domain of the target density function. &
                &It contains the upper boundaries of the cubical domain of the objective function to be sampled. &
                &When `domainCubeLimitUpper` is specified inside an external input file supplied to the sampler, &
                &it is also possible to assign only select values of `domainCubeLimitUpper` &
                &and leave the rest of the components to be assigned the default value. &
                &For example,"//NL2//&
            SKG_"+   `domainCubeLimitUpper(3:5) = 100`"//NL2//&
            SKG_"    will only set the upper limits of the third, fourth, and fifth dimensions to 100, or,"//NL2//&
            SKG_"+   `domainCubeLimitUpper(1) = 100, domainCubeLimitUpper(2) = 1.e6`"//NL2//&
            SKG_"    will set the upper limit on the first dimension to `100`, and 1.e6 on the second dimension, or,"//NL2//&
            SKG_"+   `domainCubeLimitUpper = 3 * 2.5e100`"//NL2//&
            SKG_"    will only set the upper limits on the first, second, and third dimensions to `2.5*10^100`, while the &
                     &remaining upper limits for the missing dimensions will be automatically set to the default value."//NL2//&
            SKG_"The default value for all elements of `domainCubeLimitUpper` is `"//getStr(spec%domainCubeLimitUpper%def)//SKG_"`."//NL1//&
            SKG_"Beware that some ParaMonte samplers, such as ParaNest, require the user to specify the domain boundaries explicitly."
            !$omp master
            call setResized(domainCubeLimitUpper, ndim)
            call setNAN(domainCubeLimitUpper)
            !$omp end master
        end block domainCubeLimitUpper_block

        domainErrCount_block: block
            use pm_sampling_scio, only: domainErrCount
            spec%domainErrCount%null = -huge(0_IK)
            spec%domainErrCount%def = 10000_IK
            spec%domainErrCount%desc = &
            SKG_"The simulation specification `domainErrCount` is a scalar of type `integer` beyond which the user will be warned &
                &about the newly proposed points excessively falling outside the domain of the objective function. For every `domainErrCount` &
                &consecutively-proposed new points that fall outside the domain of the objective function, the user will be warned until &
                &`domainErrCount = domainErrCountMax`, in which case the sampler returns a fatal error and the program stops globally. &
                &The counter for this warning is reset after a proposal sample from within the domain of the objective function is obtained. &
                &When out-of-domain sampling happens frequently, it strongly indicates something fundamentally wrong in the simulation. &
                &It is, therefore, important to closely inspect and monitor for frequent out-of-domain samplings. &
                &This can be done by setting `domainErrCount` to an appropriate value determined by the user. &
                &The default value for `domainErrCount` is `"//getStr(spec%domainErrCount%def)//SKG_"`."
            !$omp master
            domainErrCount = spec%domainErrCount%null
            !$omp end master
        end block domainErrCount_block

        domainErrCountMax_block: block
            use pm_sampling_scio, only: domainErrCountMax
            spec%domainErrCountMax%null = -huge(0_IK)
            spec%domainErrCountMax%def = 100000_IK
            spec%domainErrCountMax%desc = &
            SKG_"The simulation specification `domainErrCountMax` is a scalar of type `integer` beyond which the program will stop globally &
                &with a fatal error message declaring that the maximum number of proposal-out-of-domain-bounds has reached. &
                &The counter for this global-stop request is reset after a proposal is accepted as a sample from within the domain of the objective function. &
                &When out-of-domain sampling happens frequently, it strongly indicates something fundamentally wrong in the simulation. &
                &It is, therefore, important to closely inspect and monitor for frequent out-of-domain samplings. &
                &This can be done by setting `domainErrCountMax` to an appropriate value determined by the user. &
                &The default value for `domainErrCountMax` is `"//getStr(spec%domainErrCountMax%def)//SKG_"`."
            !$omp master
            domainErrCountMax = spec%domainErrCountMax%null
            !$omp end master
        end block domainErrCountMax_block

        inputFileHasPriority_block: block
            use pm_sampling_scio, only: inputFileHasPriority
            spec%inputFileHasPriority%def = .false._LK
            spec%inputFileHasPriority%desc = &
            SKG_"The simulation specification `inputFileHasPriority` is a scalar of type `logical` (Boolean) of default kind. &
            &If it is set to the logical true value (e.g., `.true.`, `true` or `.t.` from within an external input file), &
            &then, the input specifications of the simulation will be read from the user-specified input specification file &
            &and the simulation specification assignments from within the programming language environment (if any are made) &
            &will be completely ignored. If `inputFileHasPriority` is set to the logical false value, then all simulation &
            &specifications that are taken from the user-specified input file will be overwritten by their corresponding &
            &user-set input values from within the user programming environment (if such specifications are set). &
            &This feature is useful when, for example, some simulation specifications have to be computed and specified &
            &at runtime and, therefore, cannot be specified before the program execution. Currently, this functionality &
            &(i.e., prioritizing the input file values to input-procedure-argument values) is available only in the Fortran &
            &interface to the ParaMonte library routines. It can be set exclusively within an external input file, and its &
            &value is ignored in non-Fortran programming language simulation environments. &
            &The default value for `inputFileHasPriority` is `"//getStr(spec%inputFileHasPriority%def)//SKG_"`."
            !$omp master
            inputFileHasPriority = spec%inputFileHasPriority%def
            !$omp end master
        end block inputFileHasPriority_block

        outputChainFileFormat_block: block
            use pm_sampling_scio, only: outputChainFileFormat
            spec%outputChainFileFormat%compact = SKG_"compact"
            spec%outputChainFileFormat%verbose = SKG_"verbose"
            spec%outputChainFileFormat%binary = SKG_"binary"
            spec%outputChainFileFormat%isCompact = .false._LK
            spec%outputChainFileFormat%isVerbose = .false._LK
            spec%outputChainFileFormat%isBinary  = .false._LK
            spec%outputChainFileFormat%def = spec%outputChainFileFormat%compact
            spec%outputChainFileFormat%null = repeat(SUB, len(outputChainFileFormat, IK))
            spec%outputChainFileFormat%desc = &
            SKG_"The simulation specification `outputChainFileFormat` is a scalar string of maximum length `"//getStr(len(outputChainFileFormat, IK))//SKG_"` characters &
                &representing the sampler output chain file(s) format. If specified within an external input file, it must be singly or doubly quoted. &
                &Three values are possible:"//NL2//&
            SKG_"+   `outputChainFileFormat = 'compact'` or `outputChainFileFormat = 'ascii'`"//NL2//&
            SKG_"    This is the ASCII (text) file format, which is human-readable but does not preserve the full accuracy of the output values. &
                     &It is also a significantly slower chain file generation mode than the binary file format (see below). &
                     &If the compact format is specified, each of the repeating visited states in the simulation will be condensed into a single entry (row) &
                     &in the output simulation chain file. Each entry will be then assigned a sample weight that is equal to the number of repetitions of &
                     &that state in the chain. Thus, each row in the output chain file will represent a unique sample from the objective function. &
                     &This will lead to a significantly smaller ASCII chain file and faster output size than the verbose chain file format (see below)."//NL2//&
            SKG_"+   `outputChainFileFormat = 'verbose'`"//NL2//&
            SKG_"    This is the ASCII (text) file format, which is human-readable but does not preserve the full accuracy of the output values. &
                     &It is also a significantly slower chain file generation mode compared to compact and binary chain file formats (see above and below). &
                     &If the verbose format is specified, all visited states will have equal sampling weights of `1` in the output simulation chain file. &
                     &The verbose format can lead to much larger chain file sizes than the compact and binary file formats. &
                     &This is especially true if the target objective function has a high-dimensional domain."//NL2//&
            SKG_"+   `outputChainFileFormat = '"//spec%outputChainFileFormat%binary//SKG_"'`"//NL2//&
            SKG_"    This is the binary file format, which is not human-readable but preserves the exact values in the output chain file. &
                     &It is also often the fastest mode of chain file generation. If the binary file format is chosen, the chain will be &
                     &automatically output in the compact format (but as binary) to ensure the production of the smallest possible output chain file. &
                     &Binary chain files will have the "//filext%binary//SKG_" file extensions. Use the binary format if you need full accuracy &
                     &representation of the output values while having the smallest-size output chain file in the shortest time possible."//NL2//&
            SKG_"The default value for `outputChainFileFormat` is `'"//spec%outputChainFileFormat%def//SKG_"'`, which provides a reasonable trade-off &
                &between speed and output file size for the specified simulation task. The input values are case-INsensitive."
            !$omp master
            outputChainFileFormat = spec%outputChainFileFormat%null
            !$omp end master
        end block outputChainFileFormat_block

        outputColumnWidth_block: block
            use pm_sampling_scio, only: outputColumnWidth
            spec%outputColumnWidth%null = -huge(0)
            spec%outputColumnWidth%def = 0_IK
            spec%outputColumnWidth%desc = &
            SKG_"The simulation specification `outputColumnWidth` is a non-negative scalar of type `integer` that &
                &determines the width of the data columns in the formatted tabular files outputted by the sampler. &
                &If it is set to zero, the sampler will ensure that the width of each output element is &
                &the minimum possible value without losing the requested output precision. In other words, &
                &setting `outputColumnWidth = 0` will result in the smallest size for the formatted output files in ASCII format. &
                &The default value for `outputColumnWidth` is `"//getStr(spec%outputColumnWidth%def)//SKG_"`."
            !$omp master
            outputColumnWidth = spec%outputColumnWidth%null
            !$omp end master
        end block outputColumnWidth_block

        outputFileName_block: block
            use pm_sampling_scio, only: outputFileName
            character(8,SKG):: date
            character(10,SKG):: time
            call date_and_time(date,time)
            spec%outputFileName%def = method//SKG_"_"//date//SKG_"_"//time(1:6)//SKG_"_"//time(8:10)
            ! Set the default outputFileName to the same on all parallel images.
#if         CAF_ENABLED
            block
                character(63,SKG), save :: co_outputFileNameDef[*]
                if (this_image() == 1) then
                    co_outputFileNameDef = spec%outputFileName%def
                    sync images(*)
                else
                    sync images(1)
                    spec%outputFileName%def = trim(adjustl(co_outputFileNameDef[1]))
                end if
            end block
#elif       MPI_ENABLED
            block
                use mpi !mpi_f08, only: mpi_character, mpi_comm_world, mpi_bcast
                integer :: ierrMPI
                character(63) :: co_outputFileNameDef ! This must be the default kind.
                co_outputFileNameDef = spec%outputFileName%def
                ! bcast co_outputFileNameDef from image one to all others.
                call mpi_bcast  ( co_outputFileNameDef  & ! buffer
                                , 63                    & ! count
                                , mpi_character         & ! datatype
                                , 0                     & ! root
                                , mpi_comm_world        & ! comm
                                , ierrMPI               & ! ierr
                                )
                spec%outputFileName%def = trim(adjustl(co_outputFileNameDef))
            end block
!#elif       OMP_ENABLED
!            !$omp master
!            co_outputFileNameDef = spec%outputFileName%def
!            !$omp end master
!            !$omp barrier
!            spec%outputFileName%def = trim(adjustl(co_outputFileNameDef))
#endif
            spec%outputFileName%null = repeat(SUB, len(outputFileName, IK))
            spec%outputFileName%desc = &
            SKG_"The simulation specification `outputFileName` is a scalar string of maximum length `"//getStr(len(outputFileName, IK))//SKG_"` &
                &containing the path and the base of the filename for the sampler output files. If not provided by the user, &
                &the default `outputFileName` is constructed from the current date and time:"//NL2//&
            SKG_"    sampler_yyyymmdd_hhmmss_mmm"//NL2//&
            SKG_"where `sampler` is replaced with the name of the ParaMonte sampler invoked, and `yyyy`, `mm`, `dd`, `hh`, `mm`, `ss`, &
                &and `mmm` are replaced with the current year, month, day, hour, minute, second, and millisecond. &
                &In such a case, the output files' default directory will be the sampler's current working directory. &
                &If `outputFileName` is provided but ends with a separator character '/' or '\' (as in Unix or Windows OS), &
                &then its value will be used as the directory to which the sampler output files will be written. &
                &In this case, the default output file naming convention described above will be used. &
                &The specified directory will be automatically created if it does not exist already. &
                &Note that the specified path is left-adjusted and right-padded, erasing all trailing whitespace characters. &
                &The value of `outputFileName` is always automatically suffixed with `_run<i>_pid<j>_<type>.<ext>` where "//NL2//&
            SKG_"+   `<ext>` is replaced with an appropriate file extension, typically `.txt` or `.bin`, depending on the type of the file contents,"//NL2//&
            SKG_"+   `<type>` is replaced with the simulation file type (e.g., `chain`, `report`, `sample`, `restart`, etc),"//NL2//&
            SKG_"+   `<j>` is replaced with the process (image/thread) ID (PID) that generates the current simulation file,"//NL2//&
            SKG_"+   `<i>` is replaced with the simulation run number which depends on the existence of previous simulation files &
                     &with the same file name prefix and the specified value for the simulation specification `outputStatus`."
            !$omp master
            outputFileName = spec%outputFileName%null
            !$omp end master
        end block outputFileName_block

        outputPrecision_block: block
            use pm_sampling_scio, only: outputPrecision
            spec%outputPrecision%def = nint(modelr%precision * 1.2, IK)
            spec%outputPrecision%null = -huge(0_IK)
            spec%outputPrecision%desc = &
            SKG_"The simulation specification `outputPrecision` is a scalar of type `integer` representing the precision &
                &(i.e., the number of significant digits) of the real and complex numbers in the output simulation files. &
                &Any positive integer is acceptable as the input value of `outputPrecision`. However, any digits of the &
                &output real numbers beyond the actual accuracy of floating-point numbers (e.g., ~16 digits of significance for 64-bit `real`) &
                &will be meaningless and random. Set this variable to the precision of the requested floating point precision in the simulation &
                &(or to larger values) if full reproducibility of the simulation is needed in the future. However, keep in mind that higher &
                &precisions result in larger output files. This variable is ignored for binary output (if any occurs during the simulation). &
                &The binary output files preserve the full precision of numbers. &
                &The default value for `outputPrecision` depends on the `real` precision, e.g., `"//getStr(spec%outputPrecision%def)//SKG_"`."
            !$omp master
            outputPrecision = spec%outputPrecision%null
            !$omp end master
        end block outputPrecision_block

        outputReportPeriod_block: block
            use pm_sampling_scio, only: outputReportPeriod
            spec%outputReportPeriod%null = -huge(0_IK)
            spec%outputReportPeriod%def = 1000_IK
            spec%outputReportPeriod%desc = &
            SKG_"The simulation specification `outputReportPeriod` is a positive-valued scalar of type `integer`. &
                &Every `outputReportPeriod` calls to the objective function, sampling progress will be reported to the progress file. &
                &The default value for `outputReportPeriod` is `"//getStr(spec%outputReportPeriod%def)//SKG_"`."
            !$omp master
            outputReportPeriod = spec%outputReportPeriod%null
            !$omp end master
        end block outputReportPeriod_block

        outputRestartFileFormat_block: block
            use pm_sampling_scio, only: outputRestartFileFormat
            spec%outputRestartFileFormat%binary = SKG_"binary"
            spec%outputRestartFileFormat%ascii = SKG_"ascii"
            spec%outputRestartFileFormat%isBinary = .false._LK
            spec%outputRestartFileFormat%isAscii = .false._LK
            spec%outputRestartFileFormat%def = spec%outputRestartFileFormat%binary
            spec%outputRestartFileFormat%null = repeat(SUB, len(outputRestartFileFormat, IK))
            spec%outputRestartFileFormat%desc = &
            SKG_"The simulation specification `outputRestartFileFormat` is a scalar string of maximum length `"//getStr(len(outputRestartFileFormat, IK))//SKG_"` &
                &representing the output restart file(s) format used to restart an interrupted simulation. &
                &Two values are possible:"//NL2//&
            SKG_"+   `outputRestartFileFormat = '"//spec%outputRestartFileFormat%binary//SKG_"'`"//NL2//&
            SKG_"    This is the binary file format which is not human-readable but preserves the exact values of the &
                     &specification variables required for the simulation restart. This full accuracy representation is required &
                     &to exactly reproduce an interrupted simulation. The binary format is also normally the fastest mode of restart file &
                     &generation. Binary restart files will have the `"//filext%binary//SKG_"` file extensions."//NL2//&
            SKG_"+   `outputRestartFileFormat = '"//spec%outputRestartFileFormat%ascii//SKG_"'`"//NL2//&
            SKG_"    This is the ASCII (text) file format, which is human-readable but does not preserve the full accuracy of &
                     &the specification variables required for the simulation restart. It is also a significantly slower mode of &
                     &restart file generation, compared to the binary format. Therefore, its usage should be limited to situations where &
                     &the user wants to track the dynamics of simulation specifications throughout the simulation time. &
                     &ASCII restart file(s) will have the `"//filext%ascii//SKG_"` file extensions."//NL2//&
            SKG_"The default value for `outputRestartFileFormat` is `'"//spec%outputRestartFileFormat%def//SKG_"'`. &
                &Note that the input values are case-INsensitive."
            !$omp master
            outputRestartFileFormat = spec%outputRestartFileFormat%null
            !$omp end master
        end block outputRestartFileFormat_block

        outputSampleSize_block: block
            use pm_sampling_scio, only: outputSampleSize
            spec%outputSampleSize%null = -huge(0_IK)
            spec%outputSampleSize%def  = -1_IK
            spec%outputSampleSize%desc = &
            SKG_"The simulation specification `outputSampleSize` is a non-zero scalar of type `integer` whose value dictates the number of &
                &(hopefully, independent and identically distributed [i.i.d.]) samples to be drawn from the user-provided objective function. &
                &Three ranges of values are possible. If"//NL2//&
            SKG_"+   `outputSampleSize < 0`,"//NL2//&
            SKG_"    then, the absolute value of `outputSampleSize` dictates the sample size in units of the effective sample size. &
                     &The effective sample is, by definition, i.i.d., free from duplicates and residual autocorrelation. The &
                     &effective sample size is automatically determined by the sampler toward the end of the simulation. &
                     &For example:"//NL2//&
            SKG_"    +   `outputSampleSize = -1` yields the effective i.i.d. sample drawn from the objective function."//NL2//&
            SKG_"    +   `outputSampleSize = -2` yields a (potentially non-i.i.d.) sample twice as big as the effective sample."//NL2//&
            SKG_"+   `outputSampleSize > 0`,"//NL2//&
            SKG_"    then, the sample size is assumed to be in units of the number of points to be sampled. &
                     &If outputSampleSize turns out to be less than effectiveSampleSize, the resulting sample will be i.i.d.. &
                     &If outputSampleSize turns out to be larger than effectiveSampleSize, the resulting sample will be &
                     &potentially non-i.i.d.. The larger this difference, the more non-i.i.d. the resulting &
                     &final refined sample will be. For example:"//NL2//&
            SKG_"    +  `outputSampleSize = 1000` yields a `1000`-points sample from the objective function."//NL2//&
           !SKG_"    outputSampleSize = 0,"//NL2//&
           !SKG_"            then, no sample file will be generated."//NL2//&
            SKG_"The default value for `outputSampleSize` is `"//getStr(spec%outputSampleSize%def)//SKG_"`."
            !$omp master
            outputSampleSize = spec%outputSampleSize%null
            !$omp end master
        end block outputSampleSize_block

        outputSeparator_block: block
            use pm_sampling_scio, only: outputSeparator
            spec%outputSeparator%null = repeat(SUB, len(outputSeparator, IK))
            spec%outputSeparator%def = SKG_","
            spec%outputSeparator%desc = &
            SKG_"The simulation specification `outputSeparator` is a scalar string of maximum length `"//getStr(len(outputSeparator, IK))//SKG_"` &
                &containing a sequence of one or more allowed characters used to separate fields within records of tabular contents &
                &in the simulation output files. Digits, the period symbol `'.'`, and the addition and subtraction operators: `'+'` and `'-'`) are not allowed. &
                &To output in Comma-Separated-Values (CSV) format, set `outputSeparator = ','`. If the input value is not provided, &
                &the default separator `'"//spec%outputSeparator%def//SKG_"'` will be used when input `outputColumnWidth = 0`, and a single &
                &space character, '"//spec%outputSeparator%def//SKG_"' will be used when input `outputColumnWidth > 0`. &
                &A value of `'\t'` is interpreted as the TAB character. To avoid this interpretation, &
                &use '\\\t' to yield '\t' without being interpreted as the TAB character. &
                &The default value for `outputSeparator` is `'"//spec%outputSeparator%def//SKG_"'`."
            !$omp master
            outputSeparator = spec%outputSeparator%null
            !$omp end master
        end block outputSeparator_block

        outputSplashMode_block: block
            use pm_sampling_scio, only: outputSplashMode
            spec%outputSplashMode%null = repeat(SUB, len(outputSplashMode, IK))
            if (envis%c .or. envis%cpp .or. envis%fortran) then
                spec%outputSplashMode%def = spec%outputSplashMode%normal
            else
                spec%outputSplashMode%def = spec%outputSplashMode%quiet
            end if
            spec%outputSplashMode%desc = &
            SKG_"The simulation specification `outputSplashMode` is a scalar string of maximum length `"//getStr(len(spec%outputSplashMode%null, IK))//SKG_"` &
                &representing the level of information output on screen while running or postprocessing the simulation. &
                &Three values are possible:"//NL2//&
            SKG_"+   `outputSplashMode = '"//spec%outputSplashMode%normal//SKG_"'`"//NL2//&
            SKG_"    Under this option, the simulation splash and progress bar will be shown on screen in addition to all post-processing details. &
                     &This is the default behavior in compiled language environments (e.g., C, C++, Fortran, ...)."//NL2//&
            SKG_"+   `outputSplashMode = '"//spec%outputSplashMode%quiet//SKG_"'`"//NL2//&
            SKG_"    Under this option, the splash screen will be hidden from the standard output but other information will be displayed as usual. &
                     &This is the default behavior in dynamic language environments (e.g., MATLAB, Python, R, ...)."//NL2//&
            SKG_"+   `outputSplashMode = '"//spec%outputSplashMode%silent//SKG_"'`"//NL2//&
            SKG_"    Under this option, no information about the simulation will be shown on screen. &
                     &Use this option if the simulations are expected to be short and straightforward &
                     &or if the amount of text allowed for display in standard output is limited. &
                     &This situation happens, for example, in online code coverage and CI platforms."//NL2//&
            SKG_"The default value for `outputSplashMode` is `'normal'` in compiled programming language environments &
                &and `'quiet'` in dynamic programming language environments. Note that the input values are case-INsensitive."
            !$omp master
            outputSplashMode = spec%outputSplashMode%null
            !$omp end master
        end block outputSplashMode_block

        outputStatus_block: block
            use pm_sampling_scio, only: outputStatus
            spec%outputStatus%is%extend = .false._LK
            spec%outputStatus%is%repeat = .false._LK
            spec%outputStatus%is%retry = .false._LK
            spec%outputStatus%def = SKG_"extend"
            spec%outputStatus%null = repeat(SUB, len(outputStatus, IK))
            spec%outputStatus%desc = &
            SKG_"The simulation specification `outputStatus` is a scalar string with a maximum length of `"//getStr(len(outputStatus, IK))//SKG_"` characters, &
                &whose value describes the protocol for dealing with and handling the simulation output files concerning potentially existing past simulations. &
                &The string value must be enclosed by single or double quotation marks when provided as input in an external input file. &
                &Three values are possible:"//NL2//&
            SKG_"+   `outputStatus = 'extend'`"//NL2//&
            SKG_"    This is the default behavior where the sampler will search for any prior simulation output files &
                     &with the same user-specified file name prefix in the working directory to restart the simulation. &
                     &If an old interrupted set of simulation output files exists in the working directory, &
                     &the sampler will attempt to restart the simulation from the last recorded simulation state. &
                     &The restart operation may fail if the user has modified or tampered with the old simulation output files. &
                     &If prior simulation files exist and represent a complete simulation, &
                     &a new simulation run will be performed with a new set of output files starting from the last successful run. &
                     &The parameters of the new simulation are initialized based on the output of the most recent successful simulation. &
                     &For MCMC simulations, this means starting from the last sampled point in the output sample file of the previous &
                     &simulation and using an initial proposal covariance matrix constructed from the output sample of the previous simulation. &
                     &A new simulation will start if the sampler does not find any prior simulations with the same output file names. &
                     &This default behavior allows seamless restart functionality while ensuring old potentially valuable computationally &
                     &expensive simulations are not inadvertently erased and replaced by the new simulation output files."//NL2//&
            SKG_"+   `outputStatus = 'repeat'`"//NL2//&
            SKG_"    This option is nearly identical to `'extend'` except that the new simulation specifications &
                     &are not initialized from the specifications of the last successful simulation run (if any exist). &
                     &Instead, a new set of simulation files will be generated as if the last simulation run is replicated. &
                     &If the simulation configuration has not changed since the last successful simulation run, then the new simulation &
                     &output sample, chain, and restart files will be identical to those of the last successful simulation. &
                     &This outputting is primarily useful for cross-platform or cross-compiler testing and development."//NL2//&
            SKG_"+   `outputStatus = 'retry'`"//NL2//&
            SKG_"    This option is nearly identical to `'repeat'` except that the new simulation starts afresh and overwrites &
                     &any potentially existing output files from the most recent simulation with the same names without ever using them. &
                     &There is no restart functionality with this option. The most recent simulation files are deleted regardless of &
                     &completion status. This option is effectively equivalent to deleting the set of output files from the last simulation &
                     &run and rerunning the simulation with the default value `'extend'` for the specification `outputStatus`. &
                     &Use this option for quick tests or small exploratory problems where lots of quick runs must be performed."//NL2//&
            SKG_"The default value for `outputStatus` is `'"//spec%outputStatus%def//SKG_"'`. The input values are case-INsensitive."
            !$omp master
            outputStatus = spec%outputStatus%null
            !$omp end master
        end block outputStatus_block
            !SKG_"    outputStatus = 'copy-extend'"//NL2//&
            !SKG_"            This option is nearly identical to the default behavior 'extend' except for the fact that the most &
            !                &recent existing simulation output are copied into the new simulation file before starting the new simulation. &
            !                &If no simulation exists, the sampler will start a new simulation as in the default behavior 'extend'. &
            !                &If the most recent existing simulation is incomplete, the sampler will &
            !                &complete the interrupted simulation as in the default behavior 'extend'. &
            !                &Unlike 'extend', however, if the sampler detects any old completed simulation files, it will duplicate &
            !                &the most recent existing simulation under new unique output file names for the new simulation run. &
            !                &Then the sampler will use the last simulation state from the duplicated simulation &
            !                &to extend the most recently completed simulation with a brand new simulation. &
            !                &By convention, the sampler suffixes the user-specified `outputFileName` value for extensible &
            !                &simulations with `_run<i>` where `<i>` is replaced by the simulation extension attempt number. &
            !                &The 'extend' behavior is particularly useful for situations where the lack of convergence requires &
            !                &the new simulation to be extended from the last state and with the updated parameters of an older simulation. &
            !                &Successful extension of a series of simulations requires specifying `outputStatus = 'extend'` &
            !                &from the first (unextended) simulation to the last (extended) simulation. &
            !                &Beware that the repeated copy-extension of a simulation will lead to increasingly larger simulation output files, &
            !                &because each new simulation copies the entire past efforts into the new simulation files. &
            !                &This option has been added upon request by the GitHub user: https://github.com/Peku995"//NL2//&
            !SKG_"    outputStatus = 'copy-retry'"//NL2//&
            !SKG_"            This option is nearly identical to 'copy-extend' except for the fact that the most recent files are first &
            !                &deleted and the new simulation starts as if 'copy-extend' has been specified as the value of `outputStatus`."//NL2//&

        parallelism_block: block
            use pm_sampling_scio, only: parallelism
            spec%parallelism%is%multiChain = .false._LK
            spec%parallelism%is%singleChain = .false._LK
            spec%parallelism%singleChain = "singleChain"
            spec%parallelism%multiChain = "multiChain"
            spec%parallelism%def = spec%parallelism%singleChain
            spec%parallelism%null = repeat(SUB, len(parallelism, IK))
            spec%parallelism%desc = &
            SKG_"The simulation specification `parallelism` is a scalar string of maximum length `"//getStr(len(spec%parallelism%null, IK))//SKG_"` &
                &that represents the parallelization method to be used in the simulation. &
                &The string value must be enclosed by single or double quotation marks when provided in an external input file. &
                &Two options are currently supported:"//NL2//&
            SKG_"+   `parallelism = '"//spec%parallelism%multiChain//SKG_"'`"//NL2//&
            SKG_"    This method uses the Perfect Parallelism scheme, in which multiple MCMC chains are generated &
                     &independently of each other. In this case, multiple output MCMC chain files will also be generated. &
                     &The Perfect parallelism option is available only in Coarray/MPI-enabled parallel simulations (not in OpenMP or other shared-memory). &
                     &However, it can be readily emulated by running multiple independent simulations concurrently in any programming environment."//NL2//&
            SKG_"+   `parallelism = '"//spec%parallelism%singleChain//SKG_"'`"//NL2//&
            SKG_"    This method uses the fork-style parallelization scheme. In this case, a single MCMC chain file will be generated. &
                     &At each MCMC step, multiple proposal steps will be checked in parallel until one proposal is accepted. &
                     &This is the default for all parallelism paradigms supported by the ParaMonte library and the only option for shared memory parallelism."//NL2//&
            SKG_"Note that in serial mode, there is no parallelism. Therefore, this option does not affect non-parallel simulations and ignores its value. &
                &The serial mode is equivalent to either of the parallelism methods with only one simulation image (processor, core, or thread). &
                &The default value for `parallelism` is `'"//spec%parallelism%def//SKG_"'`. &
                &Note that the input values are case-INsensitive and whitespace characters are ignored."
            !$omp master
            parallelism = spec%parallelism%null
            !$omp end master
        end block parallelism_block

        parallelismMpiFinalizeEnabled_block: block
            use pm_sampling_scio, only: parallelismMpiFinalizeEnabled
            spec%parallelismMpiFinalizeEnabled%def = .true._LK
            spec%parallelismMpiFinalizeEnabled%desc = &
            SKG_"The simulation specification `parallelismMpiFinalizeEnabled` is a scalar of type `logical` (Boolean). &
                &In MPI parallel simulations, if `parallelismMpiFinalizeEnabled` is set to the logical/Boolean true value, &
                &then a call will be made to the `MPI_Finalize()` routine from inside the ParaMonte routine at the end of the &
                &simulation to finalize the MPI communications. Set this variable to the logical/Boolean if you do not want the &
                &ParaMonte library to finalize the MPI communications for you. When this specification is set within an external input file, &
                &the values `F`, `False`, `false`, `FALSE`, and `.false.` all represent the logical true value, and &
                &the values `T`, `True`, `true`, `TRUE`, and `.true.` all represent the logical true value. &
                &This is a low-level simulation specification variable relevant to MPI parallelism simulations. &
                &If you do not use MPI-routine calls in your main program, you can safely ignore this variable with its default value. &
                &If you intend the ParaMonte samplers or other MPI-enabled ParaMonte routines repeatedly in one run, &
                &then you will have to `parallelismMpiFinalizeEnabled` to the logical `false` value to prevent early finalization of the MPI-library. &
                &Note that in non-MPI-enabled simulations, such as serial and Coarray-enabled simulations, &
                &the value of this variable is completely ignored. The default value for `parallelismMpiFinalizeEnabled` is `"//getStr(spec%parallelismMpiFinalizeEnabled%def)//SKG_"`."
            !$omp master
            parallelismMpiFinalizeEnabled = spec%parallelismMpiFinalizeEnabled%def
            !$omp end master
        end block parallelismMpiFinalizeEnabled_block

        parallelismNumThread_block: block
            use pm_sampling_scio, only: parallelismNumThread
            spec%parallelismNumThread%null = -huge(spec%parallelismNumThread%null)
            spec%parallelismNumThread%def = 0_IK
            spec%parallelismNumThread%desc = &
            SKG_"The simulation specification `parallelismNumThread` is a non-negative scalar of type `integer` of kind 32-bit, &
                &representing the number of parallel user-specified objective function evaluations in a Fork-Join shared-memory parallelism. &
                &Such parallelism paradigms include OpenMP-enabled shared-memory parallel simulations in C, C++, and Fortran or shared-memory &
                &simulations in higher-level programming language environments such as MATLAB, Python, and R. &
                &This specification is currently relevant to only OpenMP-enabled parallel ParaMonte library builds or in the context of dynamic &
                &interpreted programming languages such as those mentioned above. As such, its value or presence is ignored in serial simulations &
                &or Coarray/MPI-enabled parallel simulations. Specifying `0` leads to using all available CPU threads for the requested simulation. &
                &The default value for `parallelismNumThread` is `0`, which implies using the maximum number of available threads in concurrent or &
                &OpenMP-enabled builds of the ParaMonte library."
            !$omp master
            parallelismNumThread = spec%parallelismNumThread%null
            !$omp end master
        end block parallelismNumThread_block

        plang_block: block
            use pm_sampling_scio, only: plang
            spec%plang%null = repeat(SUB, len(plang, IK))
            spec%plang%def = SKG_"The "//envname//SKG_"."
            !spec%plang%desc = &
            !SKG_"The simulation specification `plang` is a scalar string of maximum length `"//getStr(spec%plang%null)//&
            !SKG_"`. It is an internal ParaMonte variable used to provide information about other languages interface with the ParaMonte routines."
            !$omp master
            plang = spec%plang%null
            !$omp end master
        end block plang_block

        randomSeed_block: block
            use pm_sampling_scio, only: randomSeed
            spec%randomSeed%null = -huge(0_IK)
            spec%randomSeed%desc = &
            SKG_"The simulation specification `randomSeed` is a positive scalar of type `integer` of kind 32-bit &
                &whose value serves as the seed of the random number generator. When specified by the user, &
                &the seed of the simulation random number generator will be set in a specific deterministic manner to enable future &
                &replications of the simulation with the same configuration and input specifications. The default for `randomSeed` &
                &is a processor-dependent random value to ensure complete simulation randomness at every run. &
                &In parallel simulations, the seed, whether default or user-specified, is uniquely set on &
                &all processors, threads, cores, or images in all circumstances to ensure simulation randomness."
            !$omp master
            randomSeed = spec%randomSeed%null
            !$omp end master
        end block randomSeed_block

        sysInfoFilePath_block: block
            use pm_dateTime, only: getDateTime
            use pm_sampling_scio, only: sysInfoFilePath
            character(8,SKG) :: date
            call date_and_time(date)
            spec%sysInfoFilePath%def = SKG_".sysinfo."//getDateTime(format = SKG_"%Y%m")//SKG_".cache"
            spec%sysInfoFilePath%null = repeat(SUB, len(sysInfoFilePath, IK))
            !$omp master
            sysInfoFilePath = spec%sysInfoFilePath%null
            !$omp end master
        end block sysInfoFilePath_block

        targetAcceptanceRate_block: block
            use pm_sampling_scio, only: targetAcceptanceRate
            spec%targetAcceptanceRate%enabled = .false._LK
            spec%targetAcceptanceRate%def = [0._RKG, 1._RKG]
            spec%targetAcceptanceRate%desc = &
            SKG_"The simulation specification `targetAcceptanceRate` is a vector of type `real` of size `2` of the highest precision &
                &available within the ParaMonte library whose values (in the range `[0, 1]`) determine the optimal target range for the &
                &simulation efficiency, defined as the ratio of the number of accepted objective function calls to the total number of &
                &function calls made through the simulation. The first and the second elements of `targetAcceptanceRate` determine the &
                &lower and upper bounds of the desired acceptance rate, respectively. When the acceptance rate of the sampler is outside the &
                &specified limits, the simulation settings will be automatically adjusted to bring the overall acceptance rate to within the &
                &user-specified limits given by `targetAcceptanceRate`. When assigned from within a dynamic-language programming environment &
                &(such as MATLAB, Python, or R) or within an input file, `targetAcceptanceRate` can also be a scalar real number in the &
                &range `[0, 1]`. In such a case, the sampler will constantly attempt (albeit with no guarantee of success) to bring the average &
                &acceptance ratio of the sampler as close to the user-provided target ratio as possible. Specifically, the success of the adaptive &
                &MCMC samplers (e.g., ParaDRAM) in keeping the average acceptance ratio close to the requested target value depends heavily on:"//NL2//&
            SKG_"+   the specified value of `proposalAdaptationPeriod`; the larger, the easier."//NL1//&
            SKG_"+   the specified value of `proposalAdaptationCount`; the larger, the easier."//NL2//&
            SKG_"Note that the acceptance ratio adjustments will only occur in every `proposalAdaptationPeriod` &
                &sampling step for a total number of `proposalAdaptationCount` in adaptive MCMC samplings. &
                &The default value for `targetAcceptanceRate` is the range `[0, 1]`."
            !$omp master
            call setNAN(targetAcceptanceRate)
            !$omp end master
        end block targetAcceptanceRate_block

        !$omp barrier

    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function set(spec, sampler) result(err)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: set
#endif
        use pm_err, only: err_type
        use pm_sampling, only: sampler_type
        use pm_paramonte, only: getParaMonteSplash
        class(specbase_type), intent(inout) :: spec
        class(sampler_type), intent(in), optional :: sampler
        type(err_type) :: err

        integer(IK), parameter :: PROGRESS_BAR_WIDTH = 85_IK

        err%occurred = .false._LK
        err%msg = SKG_""

        ! Setup the splash screen.

        if (spec%image%is%first) then
            ! These messages are to be displayed only if outputSplashMode is set to normal.
            spec%msg = SKG_"Setting up the "//envname//SKG_" environment for a "//spec%method%val//SKG_" simulation..."
            if (allocated(sampler%inputFile)) then
                block
                    logical :: exist; integer :: iostat; exist = .false.
                    ! The following test is necessary since Intel compiler `inquire()`
                    ! implementation fails with larger than maximum file path length.
                    exist = len(sampler%inputFile, IK) <= MAX_LEN_FILE_PATH
                    if (exist) then
                        inquire(file = sampler%inputFile, exist = exist, iostat = iostat)
                        exist = exist .and. iostat == 0
                    end if
                    if (exist) then
                        spec%msg = spec%msg//NL2//SKG_"The user-specified input file for "//spec%method%val//SKG_" specifications detected."//NL1& ! LCOV_EXCL_LINE
                        //SKG_"All "//spec%method%val//SKG_" simulation specifications will be read from:"//NL1//SKG_""""//sampler%inputFile//SKG_""""
                    else
                        spec%msg = spec%msg//NL2//SKG_"The user-specified internal input file for "//spec%method%val//SKG_" specifications detected."//NL1& ! LCOV_EXCL_LINE
                        //SKG_"All "//spec%method%val//SKG_" simulation specifications will be read from the specified internal namelist file..."
                    end if
                end block
            else
                spec%msg = spec%msg//NL2//SKG_"No user-specified input file for the "//spec%method%val//SKG_" simulation specifications detected."//NL1& ! LCOV_EXCL_LINE
                //SKG_"The sampler will use the default simulation settings where needed."
            end if
        end if

        ! Read the namelist from the `inputFile`.

        if (allocated(sampler%inputFile)) then
            block
                use pm_parallelism, only: isFailedImage
                use pm_sampling_scio, only: setSpecFromInput
                !$omp master
                call setSpecFromInput(spec%method%val, sampler%inputFile, err)
                !$omp end master
                if (isFailedImage(err%occurred)) then
                    if (spec%image%is%first) then
                        call spec%disp%show(getParaMonteSplash(width = max(PROGRESS_BAR_WIDTH, spec%disp%width)))
                        call spec%disp%note%show(spec%msg)
                    end if
                    return
                end if
                !$omp barrier
            end block
        end if

        block
            use pm_sampling_scio, only: inputFileHasPriority
            spec%overridable = .not. inputFileHasPriority
        end block

        description_block: block
            use pm_sampling_scio, only: description
            use pm_strASCII, only: getAsciiFromEscaped
            if (spec%overridable .and. allocated(sampler%description)) then
                spec%description%val = getAsciiFromEscaped(trim(adjustl(sampler%description)))
            else
                spec%description%val = getAsciiFromEscaped(trim(adjustl(description)))
            end if
            if (spec%description%val == trim(adjustl(spec%description%null))) spec%description%val = trim(adjustl(spec%description%def))
        end block description_block

        domain_block: block
            use pm_sampling_scio, only: domain
            if (spec%overridable .and. allocated(sampler%domain)) then
                spec%domain%val = getStrLower(trim(adjustl(sampler%domain)))
            else
                spec%domain%val = getStrLower(trim(adjustl(domain)))
            end if
            if (spec%domain%val == spec%domain%null) spec%domain%val = spec%domain%def
            spec%domain%isBall = spec%domain%val == spec%domain%ball .or. spec%domain%val == SKG_"ball"
            spec%domain%isCube = spec%domain%val == spec%domain%cube
        end block domain_block

        domainAxisName_block: block
            use pm_sampling_scio, only: domainAxisName
            integer(IK) :: idim
            if (spec%overridable .and. allocated(sampler%domainAxisName)) then
                spec%domainAxisName%val = sampler%domainAxisName
            else
                spec%domainAxisName%val = spec%domainAxisName%def
            end if
            do concurrent(idim = 1 : size(spec%domainAxisName%val, 1, IK))
                if (trim(adjustl(domainAxisName(idim))) /= trim(adjustl(spec%domainAxisName%null))) spec%domainAxisName%val(idim) = domainAxisName(idim)
            end do
        end block domainAxisName_block

        domainBallAvg_block: block
            use pm_sampling_scio, only: domainBallAvg
            if (spec%overridable .and. allocated(sampler%domainBallAvg)) then
                spec%domainBallAvg%val = real(sampler%domainBallAvg, RKG)
            else
                spec%domainBallAvg%val = domainBallAvg
            end if
            where (isNAN(spec%domainBallAvg%val))
                spec%domainBallAvg%val = spec%domainBallAvg%def
            end where
        end block domainBallAvg_block

        domainBallCor_block: block
            use pm_sampling_scio, only: domainBallCor
            integer(IK) :: idim, jdim
            if (spec%overridable .and. allocated(sampler%domainBallCor)) then
                spec%domainBallCor%val = real(sampler%domainBallCor, RKG)
            else
                spec%domainBallCor%val = domainBallCor
            end if
            do jdim = 1, spec%ndim%val
                do idim = 1, spec%ndim%val
                    if (isNAN(spec%domainBallCor%val(idim, jdim))) then
                        spec%domainBallCor%val(idim, jdim) = spec%domainBallCor%def(idim, jdim)
                    else
                        spec%domain%isFinite = .true._LK
                    end if
                end do
            end do
        end block domainBallCor_block

        domainBallCov_block: block
            use pm_sampling_scio, only: domainBallCov
            integer(IK) :: idim, jdim
            if (spec%overridable .and. allocated(sampler%domainBallCov)) then
                spec%domainBallCov%val = real(sampler%domainBallCov, RKG)
            else
                spec%domainBallCov%val = domainBallCov
            end if
            do jdim = 1, spec%ndim%val
                do idim = 1, spec%ndim%val
                    if (isNAN(spec%domainBallCov%val(idim, jdim))) then
                        spec%domainBallCov%val(idim, jdim) = spec%domainBallCov%def(idim, jdim)
                    else
                        spec%domainBallCov%isUserSet = .true._LK
                    end if
                end do
            end do
            spec%domain%isFinite = spec%domain%isFinite .or. spec%domainBallCov%isUserSet
        end block domainBallCov_block

        domainBallStd_block: block
            use pm_sampling_scio, only: domainBallStd
            integer(IK) :: idim
            if (spec%overridable .and. allocated(sampler%domainBallStd)) then
                spec%domainBallStd%val = real(sampler%domainBallStd, RKG)
            else
                spec%domainBallStd%val = domainBallStd
            end if
            do idim = 1, spec%ndim%val
                if (isNAN(spec%domainBallStd%val(idim)) .or. isInf(spec%domainBallStd%val(idim))) then
                    spec%domainBallStd%val(idim) = spec%domainBallStd%def
                else
                    spec%domain%isFinite = .true._LK
                end if
            end do
        end block domainBallStd_block

        domainCubeLimitLower_block: block
            use pm_sampling_scio, only: domainCubeLimitLower
            integer(IK) :: idim
            if (spec%overridable .and. allocated(sampler%domainCubeLimitLower)) then
                spec%domainCubeLimitLower%val = real(sampler%domainCubeLimitLower, RKG)
            else
                spec%domainCubeLimitLower%val = domainCubeLimitLower
            end if
            do idim = 1, spec%ndim%val
                if (isNAN(spec%domainCubeLimitLower%val(idim)) .or. isInf(spec%domainCubeLimitLower%val(idim))) then
                    spec%domainCubeLimitLower%val(idim) = spec%domainCubeLimitLower%def
                else
                    spec%domain%isFinite = .true._LK
                end if
            end do
        end block domainCubeLimitLower_block

        domainCubeLimitUpper_block: block
            use pm_sampling_scio, only: domainCubeLimitUpper
            integer(IK) :: idim
            if (spec%overridable .and. allocated(sampler%domainCubeLimitUpper)) then
                spec%domainCubeLimitUpper%val = real(sampler%domainCubeLimitUpper, RKG)
            else
                spec%domainCubeLimitUpper%val = domainCubeLimitUpper
            end if
            do idim = 1, spec%ndim%val
                if (isNAN(spec%domainCubeLimitUpper%val(idim)) .or. isInf(spec%domainCubeLimitUpper%val(idim))) then
                    spec%domainCubeLimitUpper%val(idim) = spec%domainCubeLimitUpper%def
                else
                    spec%domain%isFinite = .true._LK
                end if
            end do
        end block domainCubeLimitUpper_block

        domainErrCount_block: block
            use pm_sampling_scio, only: domainErrCount
            if (spec%overridable .and. allocated(sampler%domainErrCount)) then
                spec%domainErrCount%val = sampler%domainErrCount
            else
                spec%domainErrCount%val = domainErrCount
            end if
            if (spec%domainErrCount%val == spec%domainErrCount%null) spec%domainErrCount%val = spec%domainErrCount%def
            spec%domainErrCount%str = getStr(spec%domainErrCount%val)
        end block domainErrCount_block

        domainErrCountMax_block: block
            use pm_sampling_scio, only: domainErrCountMax
            if (spec%overridable .and. allocated(sampler%domainErrCountMax)) then
                spec%domainErrCountMax%val = sampler%domainErrCountMax
            else
                spec%domainErrCountMax%val = domainErrCountMax
            end if
            if (spec%domainErrCountMax%val == spec%domainErrCountMax%null) spec%domainErrCountMax%val = spec%domainErrCountMax%def
            spec%domainErrCountMax%str = getStr(spec%domainErrCountMax%val)
        end block domainErrCountMax_block

        inputFileHasPriority_block: block
            use pm_sampling_scio, only: inputFileHasPriority
            !if (spec%overridable .and. allocated(sampler%inputFileHasPriority)) then
            !    spec%inputFileHasPriority%val = sampler%inputFileHasPriority
            !else
                spec%inputFileHasPriority%val = inputFileHasPriority
            !end if
        end block inputFileHasPriority_block

        outputChainFileFormat_block: block
            use pm_sampling_scio, only: outputChainFileFormat
            character(:,SKG), allocatable :: lowerCaseVal
            if (spec%overridable .and. allocated(sampler%outputChainFileFormat)) then
                spec%outputChainFileFormat%val = trim(adjustl(sampler%outputChainFileFormat))
            else
                spec%outputChainFileFormat%val = trim(adjustl(outputChainFileFormat))
            end if
            if (spec%outputChainFileFormat%val == trim(adjustl(spec%outputChainFileFormat%null))) spec%outputChainFileFormat%val = trim(adjustl(spec%outputChainFileFormat%def))
            lowerCaseVal = getStrLower(spec%outputChainFileFormat%val)
            spec%outputChainFileFormat%iscompact = lowerCaseVal == getStrLower(spec%outputChainFileFormat%compact) .or. lowerCaseVal == getStrLower(SKG_"ascii")
            spec%outputChainFileFormat%isverbose = lowerCaseVal == getStrLower(spec%outputChainFileFormat%verbose)
            spec%outputChainFileFormat%isBinary  = lowerCaseVal == getStrLower(spec%outputChainFileFormat%binary)
        end block outputChainFileFormat_block

        outputColumnWidth_block: block
            use pm_sampling_scio, only: outputColumnWidth
            if (spec%overridable .and. allocated(sampler%outputColumnWidth)) then
                spec%outputColumnWidth%val = sampler%outputColumnWidth
            else
                spec%outputColumnWidth%val = outputColumnWidth
            end if
            if (spec%outputColumnWidth%val == spec%outputColumnWidth%null) spec%outputColumnWidth%val = spec%outputColumnWidth%def
            spec%outputColumnWidth%str = getStr(spec%outputColumnWidth%val)
        end block outputColumnWidth_block

        outputFileName_block: block
            use pm_sampling_scio, only: outputFileName
            if (spec%overridable .and. allocated(sampler%outputFileName)) then
                spec%outputFileName%val = trim(adjustl(sampler%outputFileName))
            else
                spec%outputFileName%val = trim(adjustl(outputFileName))
            end if
            if (spec%outputFileName%val == trim(adjustl(spec%outputFileName%null))) spec%outputFileName%val = spec%outputFileName%def
        end block outputFileName_block

        outputStatus_block: block
            use pm_sampling_scio, only: outputStatus
            character(:,SKG), allocatable :: lowerCaseVal
            if (spec%overridable .and. allocated(sampler%outputStatus)) then
                spec%outputStatus%val = trim(adjustl(sampler%outputStatus))
            else
                spec%outputStatus%val = trim(adjustl(outputStatus))
            end if
            if (spec%outputStatus%val == trim(adjustl(spec%outputStatus%null))) spec%outputStatus%val = trim(adjustl(spec%outputStatus%def))
            lowerCaseVal = getRemoved(getRemoved(getStrLower(spec%outputStatus%val), SKG_" "), SKG_"-")
            spec%outputStatus%is%extend = index(lowerCaseVal, SKG_"extend") /= 0
            spec%outputStatus%is%repeat = index(lowerCaseVal, SKG_"repeat") /= 0
            spec%outputStatus%is%retry = index(lowerCaseVal, SKG_"retry") /= 0
        end block outputStatus_block

        outputPrecision_block: block
            use pm_sampling_scio, only: outputPrecision
            if (spec%overridable .and. allocated(sampler%outputPrecision)) then
                spec%outputPrecision%val = sampler%outputPrecision
            else
                spec%outputPrecision%val = outputPrecision
            end if
            if (spec%outputPrecision%val == spec%outputPrecision%null) spec%outputPrecision%val = spec%outputPrecision%def
            spec%outputPrecision%str = getStr(spec%outputPrecision%val)
        end block outputPrecision_block

        outputReportPeriod_block: block
            use pm_sampling_scio, only: outputReportPeriod
            if (spec%overridable .and. allocated(sampler%outputReportPeriod)) then
                spec%outputReportPeriod%val = sampler%outputReportPeriod
            else
                spec%outputReportPeriod%val = outputReportPeriod
            end if
            if (spec%outputReportPeriod%val == spec%outputReportPeriod%null) spec%outputReportPeriod%val = spec%outputReportPeriod%def
        end block outputReportPeriod_block

        outputRestartFileFormat_block: block
            use pm_sampling_scio, only: outputRestartFileFormat
            character(:,SKG), allocatable :: lowerCaseVal
            if (spec%overridable .and. allocated(sampler%outputRestartFileFormat)) then
                spec%outputRestartFileFormat%val = trim(adjustl(sampler%outputRestartFileFormat))
            else
                spec%outputRestartFileFormat%val = trim(adjustl(outputRestartFileFormat))
            end if
            if (spec%outputRestartFileFormat%val == trim(adjustl(spec%outputRestartFileFormat%null))) spec%outputRestartFileFormat%val = trim(adjustl(spec%outputRestartFileFormat%def))
            lowerCaseVal = getStrLower(spec%outputRestartFileFormat%val)
            spec%outputRestartFileFormat%isAscii = logical(lowerCaseVal == getStrLower(spec%outputRestartFileFormat%ascii), LK)
            spec%outputRestartFileFormat%isBinary = logical(lowerCaseVal == getStrLower(spec%outputRestartFileFormat%binary), LK)
        end block outputRestartFileFormat_block

        outputSampleSize_block: block
            use pm_sampling_scio, only: outputSampleSize
            if (spec%overridable .and. allocated(sampler%outputSampleSize)) then
                spec%outputSampleSize%val = sampler%outputSampleSize
            else
                spec%outputSampleSize%val = outputSampleSize
            end if
            if (spec%outputSampleSize%val == spec%outputSampleSize%null) spec%outputSampleSize%val = spec%outputSampleSize%def
            spec%outputSampleSize%str = getStr(spec%outputSampleSize%val)
        end block outputSampleSize_block

        outputSeparator_block: block
            use pm_strASCII, only: HT
            use pm_sampling_scio, only: outputSeparator
            if (spec%overridable .and. allocated(sampler%outputSeparator)) then
                spec%outputSeparator%val = trim(adjustl(sampler%outputSeparator))
            else
                spec%outputSeparator%val = trim(adjustl(outputSeparator))
            end if
            if (spec%outputSeparator%val == spec%outputSeparator%null) then
                if (allocated(spec%outputSeparator%val)) deallocate(spec%outputSeparator%val)
                spec%outputSeparator%val = spec%outputSeparator%def
                !if (outputColumnWidth == 0_IK) then
                !    spec%outputSeparator%val = spec%outputSeparator%def
                !else
                !    spec%outputSeparator%val = SKG_" "
                !end if
            elseif (spec%outputSeparator%val == SKG_"") then
                spec%outputSeparator%val = SKG_" "
            elseif (spec%outputSeparator%val == SKG_"\t") then
                spec%outputSeparator%val = HT
            elseif (spec%outputSeparator%val == SKG_"\\t") then
                spec%outputSeparator%val = SKG_"\t"
            end if
        end block outputSeparator_block

        outputSplashMode_block: block
            use pm_sampling_scio, only: outputSplashMode
            character(:,SKG), allocatable :: lowerCaseVal
            if (spec%overridable .and. allocated(sampler%outputSplashMode)) then
                spec%outputSplashMode%val = trim(adjustl(sampler%outputSplashMode))
            else
                spec%outputSplashMode%val = trim(adjustl(outputSplashMode))
            end if
            if (spec%outputSplashMode%val == trim(adjustl(spec%outputSplashMode%null))) spec%outputSplashMode%val = trim(adjustl(spec%outputSplashMode%def))
            lowerCaseVal = getStrLower(spec%outputSplashMode%val)
            spec%outputSplashMode%is%normal = logical(lowerCaseVal == getStrLower(spec%outputSplashMode%normal), LK)
            spec%outputSplashMode%is%silent = logical(lowerCaseVal == getStrLower(spec%outputSplashMode%silent), LK)
            spec%outputSplashMode%is%quiet = logical(lowerCaseVal == getStrLower(spec%outputSplashMode%quiet), LK)
            if (spec%image%is%first) then
                if (spec%outputSplashMode%is%normal) call spec%disp%show(getParaMonteSplash(width = max(PROGRESS_BAR_WIDTH, spec%disp%width)))
                if (.not. spec%outputSplashMode%is%silent) call spec%disp%note%show(spec%msg)
            end if
        end block outputSplashMode_block

        parallelism_block: block
            use pm_sampling_scio, only: parallelism
            character(:,SKG), allocatable :: lowerCaseVal
            if (spec%overridable .and. allocated(sampler%parallelism)) then
                spec%parallelism%val = sampler%parallelism
            else
                spec%parallelism%val = parallelism
            end if
            spec%parallelism%val = trim(adjustl(getRemoved(spec%parallelism%val, SKG_" ")))
            if (spec%parallelism%val == trim(adjustl(spec%parallelism%null))) spec%parallelism%val = trim(adjustl(spec%parallelism%def))
            lowerCaseVal = getStrLower(spec%parallelism%val)
#if         CAF_ENABLED || MPI_ENABLED
            spec%parallelism%is%singleChain = lowerCaseVal == getStrLower(spec%parallelism%singleChain)
            spec%parallelism%is%multiChain = lowerCaseVal == getStrLower(spec%parallelism%multiChain)
#else
            ! This is critical for serial, concurrent and openmp applications.
            spec%parallelism%val = spec%parallelism%singleChain
            spec%parallelism%is%singleChain = .true._LK
            spec%parallelism%is%multiChain = .false._LK
#endif
        end block parallelism_block

        parallelismMpiFinalizeEnabled_block: block
            use pm_sampling_scio, only: parallelismMpiFinalizeEnabled
            if (spec%overridable .and. allocated(sampler%parallelismMpiFinalizeEnabled)) then
                spec%parallelismMpiFinalizeEnabled%val = sampler%parallelismMpiFinalizeEnabled
            else
                spec%parallelismMpiFinalizeEnabled%val = parallelismMpiFinalizeEnabled
            end if
        end block parallelismMpiFinalizeEnabled_block

        parallelismNumThread_block: block
            use pm_sampling_scio, only: parallelismNumThread
            if (spec%overridable .and. allocated(sampler%parallelismNumThread)) then
                !print *, sampler%parallelismNumThread
                spec%parallelismNumThread%val = sampler%parallelismNumThread
            else
                spec%parallelismNumThread%val = parallelismNumThread
            end if
            if (spec%parallelismNumThread%val == spec%parallelismNumThread%null) spec%parallelismNumThread%val = spec%parallelismNumThread%def
#if         OMP_ENABLED
#if         MATLAB_ENABLED || PYTHON_ENABLED || R_ENABLED
            if (spec%parallelismNumThread%val < 0) then
                spec%image%count = abs(spec%parallelismNumThread%val)
                spec%parallelismNumThread%val = 0_IK
            end if
#endif
            if (0 < spec%parallelismNumThread%val) spec%image%count = spec%parallelismNumThread%val
#endif
        end block parallelismNumThread_block

        ! This is an exceptional internal namelist entry that end users must not access.
        plang_block: block
            use pm_sampling_scio, only: plang
            spec%plang%val = trim(adjustl(plang))
            if (spec%plang%val == trim(adjustl(spec%plang%null))) then
                spec%plang%val = spec%plang%def
            else
                spec%plang%val = spec%plang%def//SKG_" "//spec%plang%val
            end if
        end block plang_block

        randomSeed_block: block
            use pm_kind, only: IKD, RKD
            use pm_distUnif, only: setUnifRand
            use pm_sampling_scio, only: randomSeed
            use pm_distUnif, only: getUnifRand, splitmix64_type
            integer(IK), parameter :: HALF_HUGE = int(real(huge(0_IK)) * .5, IK)
            type(splitmix64_type) :: rng
            if (spec%overridable .and. allocated(sampler%randomSeed)) then
                spec%randomSeed%val = sampler%randomSeed
            else
                spec%randomSeed%val = randomSeed
            end if
            if (randomSeed == spec%randomSeed%null) then
                rng = splitmix64_type()
                spec%randomSeed%val = getUnifRand(rng, int(sqrt(real(HALF_HUGE, RKD)), IK), HALF_HUGE)
            end if
            spec%rng = xoshiro256ssw_type(seed = int(spec%randomSeed%val, IKD), imageID = spec%image%id)
            !#if     __GFORTRAN__ && CAF_ENABLED
            !        ! opencoarrays crashes without this, by somehow setting comv_randomSeed(1)%err%occurred = TRUE
            !        ! likely a result of memory corruption
            !        !if (comv_randomSeed(1)%err%occurred) write(*,*) ""
            !#endif
            !        if (comv_randomSeed(1)%err%occurred) then
            !        ! LCOV_EXCL_START
            !            err%occurred = .true._LK
            !            err%msg = err%msg//PROCEDURE_NAME//comv_randomSeed(1)%err%msg
            !            return
            !        ! LCOV_EXCL_STOP
            !        end if
            !        call comv_randomSeed(1)%get()
            !        spec%randomSeed%Seed(:,spec%randomSeed%imageID) = comv_randomSeed(1)%val(:)
            !#if     CAF_ENABLED
            !        sync all    ! allow all images to set the seed first, then fetch the values
            !        do imageID = 1, spec%randomSeed%imageCount
            !            if (imageID/=spec%randomSeed%imageID) spec%randomSeed%Seed(:,imageID) = comv_randomSeed(1)[imageID]%val(:)
            !        end do
            !#elif   MPI_ENABLED
            !        allocate(Seed(spec%randomSeed%sizeSeed,spec%randomSeed%imageCount))
            !        call mpi_barrier(mpi_comm_world,ierrMPI) ! allow all images to set the seed first, then fetch the values
            !        call mpi_allgather  ( spec%randomSeed%Seed(:,spec%randomSeed%imageID)   &   ! LCOV_EXCL_LINE : send buffer
            !                            , spec%randomSeed%sizeSeed                        &   ! LCOV_EXCL_LINE : send count
            !                            , mpi_integer                                   &   ! LCOV_EXCL_LINE : send datatype
            !                            , Seed(:,:)                                     &   ! LCOV_EXCL_LINE : receive buffer
            !                            , spec%randomSeed%sizeSeed                        &   ! LCOV_EXCL_LINE : receive count
            !                            , mpi_integer                                   &   ! LCOV_EXCL_LINE : receive datatype
            !                            , mpi_comm_world                                &   ! LCOV_EXCL_LINE : comm
            !                            , ierrMPI                                       &   ! LCOV_EXCL_LINE : ierr
            !                            )
            !        spec%randomSeed%Seed(:,:) = Seed
            !        deallocate(Seed)
            !#endif
            !        deallocate(comv_randomSeed)
        end block randomSeed_block

        sysInfoFilePath_block: block
            use pm_sampling_scio, only: sysInfoFilePath
            if (spec%overridable .and. allocated(sampler%sysInfoFilePath)) then
                spec%sysInfoFilePath%val = trim(adjustl(sampler%sysInfoFilePath))
            else
                spec%sysInfoFilePath%val = trim(adjustl(sysInfoFilePath))
            end if
            if (spec%sysInfoFilePath%val == spec%sysInfoFilePath%null) spec%sysInfoFilePath%val = spec%sysInfoFilePath%def
            if (allocated(spec%sysInfoFilePath%null)) deallocate(spec%sysInfoFilePath%null)
        end block sysInfoFilePath_block

        targetAcceptanceRate_block: block
            use pm_sampling_scio, only: targetAcceptanceRate
            logical(LK) :: lowerUpperSet(2)
            if (spec%overridable .and. allocated(sampler%targetAcceptanceRate)) then
                spec%targetAcceptanceRate%val = real(sampler%targetAcceptanceRate, RKG)
            else
                spec%targetAcceptanceRate%val = targetAcceptanceRate
            end if
            lowerUpperSet = .not. isNAN(spec%targetAcceptanceRate%val)
            if (lowerUpperSet(1) .and. .not. lowerUpperSet(2)) spec%targetAcceptanceRate%val(2) = spec%targetAcceptanceRate%val(1)
            if (lowerUpperSet(2) .and. .not. lowerUpperSet(1)) spec%targetAcceptanceRate%val(1) = spec%targetAcceptanceRate%val(2)
            if (.not.(lowerUpperSet(1) .or. lowerUpperSet(2))) spec%targetAcceptanceRate%val(:) = spec%targetAcceptanceRate%def
            spec%targetAcceptanceRate%enabled = logical(any(spec%targetAcceptanceRate%val /= spec%targetAcceptanceRate%def), LK)
            if (spec%targetAcceptanceRate%enabled) then
                spec%targetAcceptanceRate%aim = sum(spec%targetAcceptanceRate%val) * .5_RKG
            elseif (spec%method%isParaDISE .or. spec%method%isParaDRAM) then
                spec%targetAcceptanceRate%aim = .234_RKG
            elseif (spec%method%isParaNest) then
                spec%targetAcceptanceRate%aim = .2_RKG
            else
                error stop "This cannot happen." ! LCOV_EXCL_LINE
            end if
        end block targetAcceptanceRate_block

        ! Resolve the conflicting cases.

        block
            use pm_sampleCov, only: getCov, uppDia
            if (.not. spec%domainBallCov%isUserSet) spec%domainBallCov%val = getCov(spec%domainBallCor%val, uppDia, spec%domainBallStd%val)
        end block

        if (spec%outputColumnWidth%val /= 0_IK) spec%outputSeparator%val = SKG_" "

        ! Setup the parallel images info.

        spec%image%is%leader = spec%image%is%first .or. spec%parallelism%is%multiChain
        spec%image%is%rooter = .not. spec%image%is%leader

        ! Set the auxiliary dependent variables.

        spec%outputColumnWidth%max      = getStr(max(spec%outputPrecision%val, spec%outputColumnWidth%val, maxval(len_trim(spec%domainAxisName%val, IK))) + getCountDigit(spec%real%minexponent) + 6)
        spec%reportFile%format%allreal  = SKG_"('"//spec%reportFile%indent//SKG_"',*("//spec%real%ed//spec%outputColumnWidth%max//SKG_"."//spec%outputPrecision%str//spec%real%ex//SKG_",:,' '))"
        spec%reportFile%format%intreal  = SKG_"('"//spec%reportFile%indent//SKG_"',1I"//spec%outputColumnWidth%max//SKG_",' ',*("//spec%real%ed//spec%outputColumnWidth%max//SKG_"."//spec%outputPrecision%str//spec%real%ex//SKG_",:,' '))"
        spec%reportFile%format%strreal  = SKG_"('"//spec%reportFile%indent//SKG_"',1A"//spec%outputColumnWidth%max//SKG_",' ',*("//spec%real%ed//spec%outputColumnWidth%max//SKG_"."//spec%outputPrecision%str//spec%real%ex//SKG_",:,' '))"
        spec%reportFile%format%fixform  = SKG_"('"//spec%reportFile%indent//SKG_"',*(g"//spec%outputColumnWidth%max//SKG_"."//spec%outputPrecision%str//spec%real%ex//SKG_",:,' '))"
        spec%reportFile%format%integer  = SKG_"('"//spec%reportFile%indent//SKG_"',*(I"//spec%outputColumnWidth%max//SKG_",:,' '))"
        spec%reportFile%format%generic  = SKG_"('"//spec%reportFile%indent//SKG_"',*(g0,:,' '))"

        spec%restartFile%format         = SKG_"(*(g0,:,'"//NL1//SKG_"'))"
        spec%progressFile%format%header = SKG_"(*(g"//spec%outputColumnWidth%str//SKG_"."//spec%outputPrecision%str//spec%real%ex//SKG_",:,'"//spec%outputSeparator%val//SKG_"'))"
        spec%progressFile%format%rows   = spec%progressFile%format%header!SKG_"(2I"//spec%outputColumnWidth%str//SKG_"*("//spec%real%ed//spec%outputColumnWidth%str//SKG_"."//spec%outputPrecision%str//spec%real%ex//SKG_",:,'"//spec%outputSeparator%val//SKG_"'))"
        spec%sampleFile%format%header   = SKG_"(*(g"//spec%outputColumnWidth%str//SKG_"."//spec%outputPrecision%str//spec%real%ex//SKG_",:,'"//spec%outputSeparator%val//SKG_"'))"
        spec%chainFile%format%header    = spec%sampleFile%format%header

        ! open output files, report and sanitize.

        call spec%openFiles(err)
        if (err%occurred) return
        if (spec%image%is%leader) call spec%report()
        call spec%sanitize(err)

        ! The format setup of setupOutputFiles() uses the generic g0 edit descriptor.
        ! Here the format is revised to be more specific.
        ! g0 edit descriptor format is slightly more arbitrary and compiler-dependent.
        ! These must be set here (and not in `spec%openFile()`) because the required info needs
        ! to be validated first in `spec%sanitize()` after being reported to the output files.

        if (0_IK < spec%outputColumnWidth%val) then
            associate(colWidth => spec%outputColumnWidth%str, precision => spec%outputPrecision%str, sep => spec%outputSeparator%val)
                spec%sampleFile%format%rows = SKG_"("//getStr(1 + spec%ndim%val)//SKG_"("//spec%real%ed//colWidth//SKG_"."//precision//spec%real%ex//SKG_",:,'"//sep//SKG_"'))"
                if (spec%method%isParaDISE .or. spec%method%isParaDRAM) then
                    spec%chainFile%format%rows =    SKG_"("// &
                                                    SKG_"2(I"//colWidth//SKG_",'"//sep//SKG_"')"// &
                                                    SKG_","// &
                                                    SKG_"2("//spec%real%ed//colWidth//SKG_"."//precision//spec%real%ex//SKG_",'"//sep//SKG_"')"// &
                                                    SKG_","// &
                                                    SKG_"2(I"//colWidth//SKG_",'"//sep//SKG_"')"// &
                                                    SKG_","// &
                                                    spec%sampleFile%format%rows// &
                                                    SKG_")"
                elseif (spec%method%isParaNest) then
                    spec%chainFile%format%rows =    SKG_"("// &
                                                    SKG_"1(I"//colWidth//SKG_",'"//sep//SKG_"')"// &
                                                    SKG_","// &
                                                    getStr(spec%ndim%val + 5_IK)//spec%sampleFile%format%rows// &
                                                    !SKG_","// &
                                                    !SKG_"1(A1"//colWidth//SKG_",'"//sep//SKG_"')"// & ! what is this? apparently, a commented-out linefeed.
                                                    SKG_")"
                else
                    error stop getStr(__FILE__)//SK_"@"//getFine(__FILE__, __LINE__)//SK_": Internal library error occurred. This serious error must be reported to the developers." ! LCOV_EXCL_LINE
                end if
            end associate
        else
            spec%chainFile%format%rows = spec%chainFile%format%header
            spec%sampleFile%format%rows = spec%sampleFile%format%header
        end if

    end function set

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine openFiles(spec, err)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: openFiles
#endif
        use pm_err, only: err_type
        use pm_sysPath, only: isDir
        use pm_sysPath, only: isFile
        use pm_sysPath, only: getDirSep
        use pm_sysInfo, only: getSysInfo
        use pm_sysPath, only: getDirSeps
        use pm_sysPath, only: getDirName
        use pm_sysPath, only: getPathAbs
        use pm_sysPath, only: getBaseName
        use pm_sysPath, only: getDirCurrent
        use pm_sysPath, only: isFailedRemove
        use pm_sysPath, only: isFailedMakeDir
        use pm_paramonte, only: PARAMONTE_COMPILER_OPTIONS
        use pm_paramonte, only: PARAMONTE_COMPILER_VERSION
        use pm_parallelism, only: PARALLELIZATION_MODE
        use pm_arrayStrip, only: getStripped, right
        use pm_io, only: setContentsFrom, LEN_IOMSG
        use pm_paramonte, only: getParaMonteSplash
        use pm_sysShell, only: isShellWindows
        use pm_str, only: getStrWrapped
        use pm_sysPath, only: verbatim
        use pm_str, only: isEndedWith

        type(err_type), intent(inout) :: err
        class(specbase_type), intent(inout) :: spec
        character(*,SKG), parameter :: PROCEDURE_NAME = MODULE_NAME//SKG_"@openFiles()"
        character(*,SKG), parameter :: PARALLELIZATION_MODE_SKG = PARALLELIZATION_MODE
        character(:,SKG), allocatable :: strtemp
        character(LEN_IOMSG,SKG) :: errmsg
        integer(IK) :: iostat, iell
        logical(LK) :: endsWithDirSep
        logical(LK) :: failed

        ! get the directory separators.

        errmsg = repeat(" ", len(errmsg))
        spec%outputFileName%sep = getDirSep(failed, errmsg)
        if (failed) then
            err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Failed to infer directory separator: "//trim(adjustl(errmsg))//NL2 ! LCOV_EXCL_LINE
            return ! LCOV_EXCL_LINE
        end if

        spec%outputFileName%seps = getDirSeps(failed, errmsg)
        if (failed) then
            err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Failed to infer directory separator: "//trim(adjustl(errmsg))//NL2 ! LCOV_EXCL_LINE
            return ! LCOV_EXCL_LINE
        end if

        if (spec%outputFileName%val == spec%outputFileName%def) then
            spec%msg =  SKG_"No user-specified `outputFileName` for "//spec%method%val//SKG_" output files detected."//NL1//&
                        SKG_"Generating appropriate filenames for "//spec%method%val//SKG_" output files from the current date and time..."//NL2
                        spec%outputFileName%base = getBaseName(spec%outputFileName%val, spec%outputFileName%seps, verbatim)
        else
            spec%msg =  SKG_"Variable `outputFileName` detected among the user-supplied specifications of the "//spec%method%val//SKG_" sampler:"//NL1//SKG_""""//spec%outputFileName%val//SKG_""""//NL2
        end if

        ! get the current working directory, just FYI.

        strtemp = getDirCurrent(failed)
        spec%msg = spec%msg//SKG_"Path to the current working directory:"//NL1//SKG_""""//strtemp//SKG_""""//NL2

        ! For some reason, getPathAbs() (which relies on Intel specialized routine, returns a path with newline character at the end in OpenMP-parallel mode.
        spec%outputFileName%dir = getPathAbs(getDirName(spec%outputFileName%val, spec%outputFileName%seps, verbatim), failed, errmsg)
        if (failed) then
            err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred while inferring the directory name from the specified `outputFileName`. "//trim(adjustl(errmsg))//NL2 ! LCOV_EXCL_LINE
            return ! LCOV_EXCL_LINE
        end if
        if (len_trim(adjustl(spec%outputFileName%dir)) == 0_IK) then
            spec%outputFileName%dir = trim(adjustl(strtemp))!//spec%outputFileName%sep
            spec%msg = spec%msg//SKG_"All output files will be written to the current working directory:"//NL1//SKG_""""//spec%outputFileName%dir//SKG_""""//NL2
        else
            spec%msg = spec%msg//SKG_"Generating the requested directory for the "//spec%method%val//SKG_" simulation output files:"//NL1//SKG_""""//spec%outputFileName%dir//SKG_""""//NL2
        end if

        endsWithDirSep = any(isEndedWith(css_type(spec%outputFileName%val), css_type([(spec%outputFileName%seps(iell : iell), iell = 1, len(spec%outputFileName%seps))])))
        if (endsWithDirSep) spec%msg = spec%msg//SKG_"The specified `outputFileName` does not contain a filename prefix for "//spec%method%val//SKG_" output files."//NL1//&
        SKG_"Generating appropriate filenames for "//spec%method%val//SKG_" output files from the current date and time..."//NL2
        if (endsWithDirSep .or. spec%outputFileName%val == spec%outputFileName%def) then
            spec%outputFileName%base = spec%outputFileName%def
        else
            spec%outputFileName%base = getBaseName(spec%outputFileName%val, spec%outputFileName%seps, verbatim)
        end if

        ! construct the semi full name.

        spec%outputFileName%full = spec%outputFileName%dir
        endsWithDirSep = any(isEndedWith(css_type(spec%outputFileName%full), css_type([(spec%outputFileName%seps(iell : iell), iell = 1, len(spec%outputFileName%seps))])))
        if (.not. endsWithDirSep) spec%outputFileName%full = spec%outputFileName%full//spec%outputFileName%sep
        spec%outputFileName%full = spec%outputFileName%full//spec%outputFileName%base//SKG_"_run"

        ! Generate the output files directory.

        if (spec%image%is%first) then
            failed = isFailedMakeDir(spec%outputFileName%dir)
            call spec%image%sync()
        else
            call spec%image%sync()
            failed = .not. isDir(spec%outputFileName%dir)
        end if

        ! in parallel mode, ensure the directory exists before moving on.

        if (failed) then
            if (spec%image%is%first .and. .not. spec%outputSplashMode%is%silent) call spec%disp%note%show(spec%msg)
            err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred while making directory '"//spec%outputFileName%dir//SKG_"'."//NL1//trim(errmsg)//NL2
            err%occurred = .true._LK ! LCOV_EXCL_LINE
            return ! LCOV_EXCL_LINE
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Only leader images: Generate the output filenames, search for pre-existing runs, and open the report file.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        spec%chainFile      %openArg_type = openArg_type(status = SKG_"replace")
        spec%sampleFile     %openArg_type = openArg_type(status = SKG_"replace")
        spec%reportFile     %openArg_type = openArg_type(status = SKG_"replace")
        spec%restartFile    %openArg_type = openArg_type(status = SKG_"replace")
        spec%progressFile   %openArg_type = openArg_type(status = SKG_"replace")

        spec%chainFile      %kind = SKG_"chain"
        spec%sampleFile     %kind = SKG_"sample"
        spec%reportFile     %kind = SKG_"report"
        spec%restartFile    %kind = SKG_"restart"
        spec%progressFile   %kind = SKG_"progress"

        spec%chainFile      %ext = filext%ascii
        spec%sampleFile     %ext = filext%ascii
        spec%reportFile     %ext = filext%ascii
        spec%restartFile    %ext = filext%ascii
        spec%progressFile   %ext = filext%ascii

        ! Reset the file extensions as necessary.

        if (spec%outputChainFileFormat%isBinary) then
            spec%chainFile%form = SKG_"unformatted"
            spec%chainFile%access = SKG_"stream"
            spec%chainFile%ext = filext%binary
        end if

        if (spec%outputRestartFileFormat%isBinary) then
            spec%restartFile%form = SKG_"unformatted"
            spec%restartFile%access = SKG_"stream"
            spec%restartFile%ext = filext%binary
        end if

        ! Define the file full suffix.

        strtemp = SKG_"_pid"//getStr(merge(spec%image%id, 1_IK, spec%parallelism%is%multiChain))//SKG_"_"
        spec%chainFile      %suffix = strtemp//spec%chainFile    %kind//spec%chainFile    %ext
        spec%sampleFile     %suffix = strtemp//spec%sampleFile   %kind//spec%sampleFile   %ext
        spec%reportFile     %suffix = strtemp//spec%reportFile   %kind//spec%reportFile   %ext
        spec%restartFile    %suffix = strtemp//spec%restartFile  %kind//spec%restartFile  %ext
        spec%progressFile   %suffix = strtemp//spec%progressFile %kind//spec%progressFile %ext


        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Search for older simulations with the same naming prefix and determine simulation start mode (restart vs. new).
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (spec%isFailedResizeList(1_IK, errmsg)) then
            err%occurred = .true._LK
            err%msg = PROCEDURE_NAME//getFine(__FILE__, __LINE__)//trim(errmsg)
            return
        end if

        spec%run%id = 0_IK
        loopListFiles: do

            spec%run%id = spec%run%id + 1_IK
            if (size(spec%chainFile%list, 1, IK) < spec%run%id) then
                if (spec%isFailedResizeList(spec%run%id * 2_IK, errmsg)) then
                    err%occurred = .true._LK
                    err%msg = PROCEDURE_NAME//getFine(__FILE__, __LINE__)//trim(errmsg)
                    return
                end if
            end if

            strtemp = getStr(spec%run%id)
            spec%chainFile      %list(spec%run%id)%val = spec%outputFileName%full//strtemp//spec%chainFile    %suffix
            spec%sampleFile     %list(spec%run%id)%val = spec%outputFileName%full//strtemp//spec%sampleFile   %suffix
            spec%reportFile     %list(spec%run%id)%val = spec%outputFileName%full//strtemp//spec%reportFile   %suffix
            spec%restartFile    %list(spec%run%id)%val = spec%outputFileName%full//strtemp//spec%restartFile  %suffix
            spec%progressFile   %list(spec%run%id)%val = spec%outputFileName%full//strtemp//spec%progressFile %suffix
            spec%chainFile      %extant = isFile(spec%chainFile      %list(spec%run%id)%val)
            spec%sampleFile     %extant = isFile(spec%sampleFile     %list(spec%run%id)%val)
            spec%reportFile     %extant = isFile(spec%reportFile     %list(spec%run%id)%val)
            spec%restartFile    %extant = isFile(spec%restartFile    %list(spec%run%id)%val)
            spec%progressFile   %extant = isFile(spec%progressFile   %list(spec%run%id)%val)

            ! exceptional case: chainFile

            if (.not. spec%chainFile%extant) then
                if (spec%outputChainFileFormat%isBinary) then
                    strtemp = getStripped(spec%chainFile%list(spec%run%id)%val, filext%binary, right)//filext%ascii
                else
                    strtemp = getStripped(spec%chainFile%list(spec%run%id)%val, filext%ascii, right)//filext%binary
                end if
                spec%chainFile%extant = isFile(strtemp)
                if (spec%chainFile%extant) spec%chainFile%list(spec%run%id)%val = strtemp
            end if

            ! exceptional case: restartFile

            if (.not. spec%restartFile%extant) then
                if (spec%outputChainFileFormat%isBinary) then
                    strtemp = getStripped(spec%restartFile%list(spec%run%id)%val, filext%binary, right)//filext%ascii
                else
                    strtemp = getStripped(spec%restartFile%list(spec%run%id)%val, filext%ascii, right)//filext%binary
                end if
                spec%restartFile%extant = isFile(strtemp)
                if (spec%restartFile%extant) spec%restartFile%list(spec%run%id)%val = strtemp
            end if

            ! If all files exist, search for the next batch.

            if (spec%chainFile      %extant .and. & ! LCOV_EXCL_LINE
                spec%sampleFile     %extant .and. & ! LCOV_EXCL_LINE
                spec%reportFile     %extant .and. & ! LCOV_EXCL_LINE
                spec%restartFile    %extant .and. & ! LCOV_EXCL_LINE
                spec%progressFile   %extant) cycle
            exit ! only if no files or only some files exist.

        end do loopListFiles

        spec%run%is%new = .not. (spec%reportFile%extant .or. spec%progressFile%extant .or. spec%restartFile%extant .or. spec%chainFile%extant .or. spec%sampleFile%extant) ! fresh (no restart) simulation if no file exists.
        spec%msg = spec%msg//getStr(spec%run%id - 1)//trim(merge(SKG_" count ", SKG_" counts", spec%run%id <= 2))//SKG_" of preexisting complete simulation runs with the same name prefix were detected."//NL2

        if (spec%outputStatus%is%retry) then

            ! The following is possible only if the `outpuStatus` is in `retry` mode and the most recent simulation is complete, in which case, the file names must be regenerated.
            ! This ensures `run%id` is set to the latest existing simulation run (complete or incomplete) which will be overwritten by new fresh simulation.
            ! The latest incomplete or complete simulation files will be deleted (replaced) in retry mode.
            spec%outputStatus%is%extend = spec%outputStatus%is%retry
            spec%outputStatus%is%retry = .false._LK
            if (spec%run%is%new) then
                if (1_IK < spec%run%id) spec%run%id = spec%run%id - 1_IK
                spec%msg = spec%msg//SKG_"The last complete simulation run #"//getStr(spec%run%id)//" will be overwritten as requested."//NL1
            else
                spec%msg = spec%msg//SKG_"The last incomplete simulation run #"//getStr(spec%run%id)//" will be overwritten as requested."//NL1
                spec%run%is%new = .true._LK
            end if

            ! Delete all the last batch of existing files in retry mode.

            if (spec%image%is%first .or. spec%parallelism%is%multiChain) then
                if (isFailedRemove(spec%chainFile   %list(spec%run%id)%val, forced = .true._LK, ntry = 10_IK, errmsg = errmsg)) then; failed = .true._LK; err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//trim(errmsg); end if
                if (isFailedRemove(spec%sampleFile  %list(spec%run%id)%val, forced = .true._LK, ntry = 10_IK, errmsg = errmsg)) then; failed = .true._LK; err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//trim(errmsg); end if
                if (isFailedRemove(spec%reportFile  %list(spec%run%id)%val, forced = .true._LK, ntry = 10_IK, errmsg = errmsg)) then; failed = .true._LK; err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//trim(errmsg); end if
                if (isFailedRemove(spec%restartFile %list(spec%run%id)%val, forced = .true._LK, ntry = 10_IK, errmsg = errmsg)) then; failed = .true._LK; err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//trim(errmsg); end if
                if (isFailedRemove(spec%progressFile%list(spec%run%id)%val, forced = .true._LK, ntry = 10_IK, errmsg = errmsg)) then; failed = .true._LK; err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//trim(errmsg); end if
                if (failed) then
                    call spec%disp%note%show(spec%msg//err%msg)
                    err%occurred = .true._LK
                    return
                end if
            end if

        else

            spec%msg = spec%msg//SKG_"Checking for simulation restart possibility..."//NL2

        end if

        spec%chainFile      %file = spec%chainFile      %list(spec%run%id)%val
        spec%sampleFile     %file = spec%sampleFile     %list(spec%run%id)%val
        spec%reportFile     %file = spec%reportFile     %list(spec%run%id)%val
        spec%restartFile    %file = spec%restartFile    %list(spec%run%id)%val
        spec%progressFile   %file = spec%progressFile   %list(spec%run%id)%val
        !print *, spec%chainFile      %file
        !print *, spec%sampleFile     %file
        !print *, spec%reportFile     %file
        !print *, spec%restartFile    %file
        !print *, spec%progressFile   %file

        if (spec%sampleFile%extant) then
            ! compromised simulation because a file other than the output sample is missing.
            err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. All output files are essential for a successful simulation restart."//NL1//SKG_"List of missing simulation output files:"
            if (.not. spec%chainFile    %extant) err%msg = err%msg//NL1//SKG_""""//spec%chainFile    %file//SKG_""""
            if (.not. spec%reportFile   %extant) err%msg = err%msg//NL1//SKG_""""//spec%reportFile   %file//SKG_""""
            if (.not. spec%restartFile  %extant) err%msg = err%msg//NL1//SKG_""""//spec%restartFile  %file//SKG_""""
            if (.not. spec%progressFile %extant) err%msg = err%msg//NL1//SKG_""""//spec%progressFile %file//SKG_""""
            if (spec%image%is%first .and. .not. spec%outputSplashMode%is%silent) call spec%disp%note%show(spec%msg)
            err%occurred = .true._LK
            return
        elseif (spec%run%is%new) then
            if (1_IK < spec%run%id) then
                strtemp = SKG_" from the most recent completed simulation"
            else
                strtemp = SKG_""
            end if
            spec%msg = spec%msg//SKG_"Starting a fresh simulation run #"//getStr(spec%run%id)//strtemp//SKG_"..."//NL2
            spec%chainFile      %status = SKG_"new"
            spec%sampleFile     %status = SKG_"new"
            spec%reportFile     %status = SKG_"new"
            spec%restartFile    %status = SKG_"new"
            spec%progressFile   %status = SKG_"new"
            !block
            !    use pm_paramonte, only: PARAMONTE_WEB_ISSUES
            !    if (spec%outputStatus%isCopy) error stop SK_"copy-extension facility is not yet implemented. Please report this to the ParaMonte developers at: "//PARAMONTE_WEB_ISSUES
            !end block
        else ! legible restart (dry) mode.
            spec%msg = spec%msg//SKG_"Restarting the existing incomplete simulation run #"//getStr(spec%run%id)//NL2
            spec%chainFile      %status = SKG_"old"
            spec%sampleFile     %status = SKG_"new"
            spec%reportFile     %status = SKG_"old"
            spec%restartFile    %status = SKG_"old"
            spec%progressFile   %status = SKG_"old"
            spec%reportFile     %position = SKG_"append"
        end if
        spec%run%is%dry = .not. spec%run%is%new

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! print the stdout message for generating / appending the output report file(s):
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        spec%msg = spec%msg//SKG_"Running the simulation in "//PARALLELIZATION_MODE_SKG//SKG_" on "//getStr(spec%image%count)//trim(merge(SKG_" process...  ", SKG_" processes...", spec%image%count == 1_IK))//NL2

        ! ensure all images sync here to avoid wrong inquire() result for the existence of the files.

        call spec%image%sync()

        ! open the output files

        blockLeaderFileSetup: if (spec%image%is%leader) then

            call spec%openFile(spec%reportFile, err, spec%msg)
            if (err%occurred) return

            spec%msg = spec%msg//spec%reportFile%file//NL2
            spec%disp%unit = spec%reportFile%unit
            spec%disp%mark%unit = spec%disp%unit
            spec%disp%note%unit = spec%disp%unit
            spec%disp%warn%unit = spec%disp%unit
            spec%disp%stop%unit = spec%disp%unit
            spec%disp%text%unit = spec%disp%unit
            spec%disp%text%width = 128_IK
            spec%disp%width = 128_IK

            ! rewrite the accumulated info to the stdout and all report files.

            !if (spec%run%is%new) then
                if (isFile(spec%sysInfoFilePath%val)) then
                    call setContentsFrom(spec%sysInfoFilePath%val, strtemp, iostat, errmsg)
                    if (iostat /= 0_IK) strtemp = UNDEFINED ! LCOV_EXCL_LINE
                else
                    ! This can take a long time, for example, about 0.75 seconds to execute on Stampede Login nodes.
                    ! Get system info via only the first image.
                    ! On many parallel processors via singleChain this leads to
                    ! the creation of thousands of files on the system, simultaneously.
                    ! This is not needed by any process other than the leader images.
                    strtemp = getStripped(getSysInfo(failed), NL1)
                end if
                call spec%disp%show(getParaMonteSplash())
                call spec%disp%text%wrap(NL1//SKG_"ParaMonte.library.interface.specifications"//NL1)
                call spec%disp%show(getStrWrapped(spec%plang%val))
                call spec%disp%text%wrap(NL1//SKG_"ParaMonte.library.compiler.version"//NL1)
                call spec%disp%show(getStrWrapped(PARAMONTE_COMPILER_VERSION))
                call spec%disp%text%wrap(NL1//SKG_"ParaMonte.library.compiler.options"//NL1)
                call spec%disp%show(getStrWrapped(PARAMONTE_COMPILER_OPTIONS))
                call spec%disp%text%wrap(NL1//SKG_"ParaMonte.runtime.platform.specifications"//NL1)
                call spec%disp%show(getStrWrapped(strtemp))
                call spec%disp%text%wrap(NL1//spec%method%val//".simulation.environment.setup"//NL1)
            !end if


            ! The sample file will be opened by the samplers at the time of writing to the file. Its non-existence is indicative of restart mode.
            !call spec%openFile(spec%sampleFile, err, spec%msg); if (err%occurred) return ! open/append the output sampleFile.
            call spec%openFile(spec%progressFile, err, spec%msg); if (err%occurred) return ! open/append the output progressFile.
            call spec%openFile(spec%restartFile, err, spec%msg); if (err%occurred) return ! open/append the output restartFile.
            call spec%openFile(spec%chainFile, err, spec%msg); if (err%occurred) return ! open/append the output chainFile.

            call spec%disp%note%show(spec%msg//SKG_"Done.")

            ! read the header line of the time file, only by leader images

            iostat = 0
            if (spec%run%is%dry) read(spec%progressFile%unit, *, iostat = iostat)
            if (spec%run%is%new .or. is_iostat_end(iostat)) then
                write(spec%progressFile%unit,spec%progressFile%format%header) "numFuncCallTotal" &
                                                                            , "numFuncCallAccepted" &
                                                                            , "meanAcceptanceRateSinceStart" &
                                                                            , "meanAcceptanceRateSinceLastReport" &
                                                                            , "timeElapsedSinceLastReportInSeconds" &
                                                                            , "timeElapsedSinceStartInSeconds" &
                                                                            , "timeRemainedToFinishInSeconds"
                flush(spec%progressFile%unit)
            end if

        end if blockLeaderFileSetup

    end subroutine openFiles

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine openFile(spec, samplerFile, err, msg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: openFile
#endif
        use pm_arrayStrip, only: getStripped, right
        !use pm_str, only: isEndedWith
        use pm_err, only: err_type
        use pm_io, only: getFileUnit
        type(err_type), intent(inout) :: err
        class(specbase_type), intent(inout) :: spec
        character(:,SKG), intent(inout), allocatable, optional :: msg
        class(samplerFile_type), intent(inout) :: samplerFile
        character(*,SKG), parameter :: PROCEDURE_NAME = MODULE_NAME//SKG_"@openFiles()"
        ! open/append the output files:
        if (present(msg)) msg = msg//merge(SKG_"Generating the new output ", SKG_"Appending to the existing ", spec%run%is%new)//samplerFile%kind//SKG_" file:"//NL1
        ! for some unknown reason, if newunit is used, GFortran opens the file as an internal file. Therefore, do not use `newunit`.
        !samplerFile%unit = getFileUnit()
        ! The Intel ifort SHARED attribute is essential for file unlocking on Windows OS.
        if (samplerFile%access /= SKG_"direct") then
            open(newunit = samplerFile%unit, file = samplerFile%file, form = samplerFile%form, access = samplerFile%access, status = samplerFile%status, iostat = samplerFile%iostat, iomsg = samplerFile%iomsg SHARED, position = samplerFile%position)
        else
            open(newunit = samplerFile%unit, file = samplerFile%file, form = samplerFile%form, access = samplerFile%access, status = samplerFile%status, iostat = samplerFile%iostat, iomsg = samplerFile%iomsg SHARED)
        end if
        !block
        !use pm_io, only: isOpen
        !integer(IK) :: lenrec
        !character(:), allocatable :: record
        !select type (samplerFile)
        !type is (chainFile_type)
        !    associate(thisFile => samplerFile)
        !    print *, "isOpen(thisFile%unit)", isOpen(thisFile%unit)
        !    if (isOpen(thisFile%unit)) close(thisFile%unit)
        !    open(newunit = thisFile%unit, file = thisFile%file, form = "unformatted", access = "stream", status = "old")
        !    call setResized(record, 162_IK)
        !    !read(thisFile%unit) record
        !    !print *, record(1:1) == achar(0)
        !    !print *, """"//record//""""
        !    lenrec = 0_IK
        !    do
        !        lenrec = lenrec + 1_IK
        !        print *, lenrec
        !        if (len(record, IK) < lenrec) call setResized(record)
        !        read(thisFile%unit) record(lenrec : lenrec)
        !        if (record(lenrec : lenrec) /= achar(0)) cycle
        !        lenrec = lenrec - 1_IK
        !        exit
        !    end do
        !end associate
        !end select
        !end block
        if (samplerFile%iostat /= 0_IK) then
            ! LCOV_EXCL_START
            if (spec%image%is%first .and. same_type_as(samplerFile, spec%reportFile) .and. present(msg) .and. .not. spec%outputSplashMode%is%silent) call spec%disp%note%show(getStripped(msg, NL1, right))
            err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred while opening the "//spec%method%val//SKG_" *"//samplerFile%suffix//SKG_" file="""//samplerFile%file//SKG_""". "
            if (scan(SKG_" ", samplerFile%file) /= 0_IK) then
                err%msg = err%msg//& ! LCOV_EXCL_LINE
                SKG_"It appears that absolute path used for the output files contains whitespace characters. This could be one potential cause of the simulation failure. &
                &The whitespace characters are often problematic in paths. Ensure the path used for the output files does not contain whitespace characters."
            end if
            err%msg = err%msg//trim(samplerFile%iomsg)
            err%occurred = .true._LK
            ! LCOV_EXCL_STOP
        elseif (spec%image%is%first .and. same_type_as(samplerFile, spec%reportFile)) Then
            if (present(msg) .and. .not. spec%outputSplashMode%is%silent) call spec%disp%note%show(msg, bmsize = 0_IK)
            block
                use pm_arrayReplace, only: getReplaced
                character(:,SKG), allocatable :: filelist
                integer(IK) :: imageID
                filelist = SKG_""
                do imageID = 1, merge(spec%image%count, 1_IK, spec%parallelism%is%multiChain)
                    filelist = filelist//SKG_""""//getReplaced(spec%reportFile%file, spec%reportFile%suffix, replacement = getReplaced(spec%reportFile%suffix, SKG_"pid1", SKG_"pid"//getStr(imageID)))//SKG_""""//NL1
                end do
                if (.not. spec%outputSplashMode%is%silent) call spec%disp%note%show(filelist//NL1//SKG_"Please see the output *"//spec%reportFile%suffix//SKG_" and *"//spec%progressFile%suffix//SKG_" files for further realtime simulation details...", tmsize = 0_IK, bmsize = 1_IK)
            end block
        else
            if (present(msg)) msg = msg//SKG_""""//samplerFile%file//SKG_""""//NL2
        end if
    end subroutine openFile

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function isFailedResizeList(spec, newsize, errmsg) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedResizeList
#endif
        class(specbase_type), intent(inout) :: spec
        character(*,SKG), intent(out) :: errmsg
        integer(IK), intent(in) :: newsize
        logical(LK) :: failed
        call setResized(spec%chainFile      %list, newsize, failed = failed, errmsg = errmsg); if (failed) then; errmsg = MODULE_NAME//SK_"@isFailedResizeList(): "//getFine(__FILE__, __LINE__)//trim(errmsg); return; end if
        call setResized(spec%sampleFile     %list, newsize, failed = failed, errmsg = errmsg); if (failed) then; errmsg = MODULE_NAME//SK_"@isFailedResizeList(): "//getFine(__FILE__, __LINE__)//trim(errmsg); return; end if
        call setResized(spec%reportFile     %list, newsize, failed = failed, errmsg = errmsg); if (failed) then; errmsg = MODULE_NAME//SK_"@isFailedResizeList(): "//getFine(__FILE__, __LINE__)//trim(errmsg); return; end if
        call setResized(spec%restartFile    %list, newsize, failed = failed, errmsg = errmsg); if (failed) then; errmsg = MODULE_NAME//SK_"@isFailedResizeList(): "//getFine(__FILE__, __LINE__)//trim(errmsg); return; end if
        call setResized(spec%progressFile   %list, newsize, failed = failed, errmsg = errmsg); if (failed) then; errmsg = MODULE_NAME//SK_"@isFailedResizeList(): "//getFine(__FILE__, __LINE__)//trim(errmsg); return; end if
    end function isFailedResizeList

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! The report method must exclusively contain reposting and no other tasks.
    subroutine report(spec)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: report
#endif
        use pm_kind, only: SK
        class(specbase_type), intent(inout) :: spec

        call spec%disp%text%wrap(NL1//spec%method%val//SKG_".simulation.specifications.base"//NL1)

        associate(ndim => spec%ndim%val, format => spec%reportFile%format%generic)

            call spec%disp%show("ndim")
            call spec%disp%show(spec%ndim%val, format = format)
            call spec%disp%note%show(spec%ndim%desc)

            call spec%disp%show("method")
            call spec%disp%show(spec%method%val, format = format)
            call spec%disp%note%show(spec%method%desc)

            call spec%disp%show("modelr")
            call spec%disp%show(SK_"      digits: "//getStr(spec%real%      digits          , signed = .true._LK), format = format, tmsize = 0_IK, bmsize = 0_IK)
            call spec%disp%show(SK_"        kind: "//getStr(spec%real%        kind          , signed = .true._LK), format = format, tmsize = 0_IK, bmsize = 0_IK)
            call spec%disp%show(SK_" maxexponent: "//getStr(spec%real% maxexponent          , signed = .true._LK), format = format, tmsize = 0_IK, bmsize = 0_IK)
            call spec%disp%show(SK_" minexponent: "//getStr(spec%real% minexponent          , signed = .true._LK), format = format, tmsize = 0_IK, bmsize = 0_IK)
            call spec%disp%show(SK_"   precision: "//getStr(spec%real%   precision          , signed = .true._LK), format = format, tmsize = 0_IK, bmsize = 0_IK)
            call spec%disp%show(SK_"       radix: "//getStr(spec%real%       radix          , signed = .true._LK), format = format, tmsize = 0_IK, bmsize = 0_IK)
            call spec%disp%show(SK_"       range: "//getStr(spec%real%       range          , signed = .true._LK), format = format, tmsize = 0_IK, bmsize = 0_IK)
            call spec%disp%show(SK_"storage_size: "//getStr(spec%real%storage_size          , signed = .true._LK), format = format, tmsize = 0_IK, bmsize = 0_IK)
            call spec%disp%show(SK_"     epsilon: "//getStr(real(spec%real% epsilon, RKG)   , signed = .true._LK), format = format, tmsize = 0_IK, bmsize = 0_IK)
            call spec%disp%show(SK_"        huge: "//getStr(real(spec%real%    huge, RKG)   , signed = .true._LK), format = format, tmsize = 0_IK, bmsize = 0_IK)
            call spec%disp%show(SK_"        tiny: "//getStr(real(spec%real%    tiny, RKG)   , signed = .true._LK), format = format, tmsize = 0_IK)
            !associate(modelr => spec%real)
            !call spec%disp%show("modelr%digits", format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show( modelr%digits , format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show("modelr%kind", format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show( modelr%kind , format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show("modelr%maxexponent", format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show( modelr%maxexponent , format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show("modelr%minexponent", format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show( modelr%minexponent , format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show("modelr%precision", format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show( modelr%precision , format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show("modelr%radix", format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show( modelr%radix , format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show("modelr%range", format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show( modelr%range , format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show("modelr%storage_size", format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show( modelr%storage_size , format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show("modelr%epsilon", format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show( modelr%epsilon , format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show("modelr%huge", format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show( modelr%huge , format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show("modelr%tiny", format = format, tmsize = 0_IK, bmsize = 0_IK)
            !call spec%disp%show( modelr%tiny , format = format, tmsize = 0_IK, bmsize = 0_IK)
            !end associate
            call spec%disp%note%show(spec%real%desc)

            call spec%disp%show("description")
           !call spec%disp%show(spec%description%val, format = format)
            call spec%disp%mark%show(spec%description%val)
            call spec%disp%note%show(spec%description%desc)

            call spec%disp%show("domain")
            call spec%disp%show(spec%domain%val, format = format)
            call spec%disp%note%show(spec%domain%desc)

            call spec%disp%show("domainAxisName")
            call spec%disp%show(reshape(spec%domainAxisName%val, [ndim, 1_IK]), format = format)
            call spec%disp%note%show(spec%domainAxisName%desc)

            call spec%disp%show("domainBallAvg")
            call spec%disp%show(reshape(spec%domainBallAvg%val, [ndim, 1_IK]), format = format)
            call spec%disp%note%show(spec%domainBallAvg%desc)

            call spec%disp%show("domainBallCor")
            call spec%disp%show(spec%domainBallCor%val, format = format)
            call spec%disp%note%show(spec%domainBallCor%desc)

            call spec%disp%show("domainBallCov")
            call spec%disp%show(spec%domainBallCov%val, format = format)
            call spec%disp%note%show(spec%domainBallCov%desc)

            call spec%disp%show("domainBallStd")
            call spec%disp%show(reshape(spec%domainBallStd%val, [ndim, 1_IK]), format = format)
            call spec%disp%note%show(spec%domainBallStd%desc)

            call spec%disp%show("domainCubeLimitLower")
            call spec%disp%show(reshape(spec%domainCubeLimitLower%val, [ndim, 1_IK]), format = format)
            call spec%disp%note%show(spec%domainCubeLimitLower%desc)

            call spec%disp%show("domainCubeLimitUpper")
            call spec%disp%show(reshape(spec%domainCubeLimitUpper%val, [ndim, 1_IK]), format = format)
            call spec%disp%note%show(spec%domainCubeLimitUpper%desc)

            call spec%disp%show("domainErrCount")
            call spec%disp%show(spec%domainErrCount%val, format = format)
            call spec%disp%note%show(spec%domainErrCount%desc)

            call spec%disp%show("domainErrCountMax")
            call spec%disp%show(spec%domainErrCountMax%val, format = format)
            call spec%disp%note%show(spec%domainErrCountMax%desc)

            call spec%disp%show("inputFileHasPriority")
            call spec%disp%show(spec%inputFileHasPriority%val, format = format)
            call spec%disp%note%show(spec%inputFileHasPriority%desc)

            call spec%disp%show("outputChainFileFormat")
            call spec%disp%show(spec%outputChainFileFormat%val, format = format)
            call spec%disp%note%show(spec%outputChainFileFormat%desc)

            call spec%disp%show("outputColumnWidth")
            call spec%disp%show(spec%outputColumnWidth%val, format = format)
            call spec%disp%note%show(spec%outputColumnWidth%desc)

            call spec%disp%show("outputFileName")
            call spec%disp%show(spec%outputFileName%val, format = format)
            call spec%disp%note%show(spec%outputFileName%desc)

            call spec%disp%show("outputPrecision")
            call spec%disp%show(spec%outputPrecision%val, format = format)
            call spec%disp%note%show(spec%outputPrecision%desc)

            call spec%disp%show("outputReportPeriod")
            call spec%disp%show(spec%outputReportPeriod%val, format = format)
            call spec%disp%note%show(spec%outputReportPeriod%desc)

            call spec%disp%show("outputRestartFileFormat")
            call spec%disp%show(spec%outputRestartFileFormat%val, format = format)
            call spec%disp%note%show(spec%outputRestartFileFormat%desc)

            call spec%disp%show("outputSampleSize")
            call spec%disp%show(spec%outputSampleSize%val, format = format)
            call spec%disp%note%show(spec%outputSampleSize%desc)

            call spec%disp%show("outputSeparator")
            call spec%disp%show(spec%outputSeparator%val, format = format)
            call spec%disp%note%show(spec%outputSeparator%desc)

            call spec%disp%show("outputSplashMode")
            call spec%disp%show(spec%outputSplashMode%val, format = format)
            call spec%disp%note%show(spec%outputSplashMode%desc)

            call spec%disp%show("outputStatus")
            call spec%disp%show(spec%outputStatus%val, format = format)
            call spec%disp%note%show(spec%outputStatus%desc)

            call spec%disp%show("parallelism")
            call spec%disp%show(spec%parallelism%val, format = format)
            call spec%disp%note%show(spec%parallelism%desc)

            call spec%disp%show("parallelismMpiFinalizeEnabled")
            call spec%disp%show(spec%parallelismMpiFinalizeEnabled%val, format = format)
            call spec%disp%note%show(spec%parallelismMpiFinalizeEnabled%desc)

            call spec%disp%show("parallelismNumThread")
            call spec%disp%show(spec%parallelismNumThread%val, format = format)
            call spec%disp%note%show(spec%parallelismNumThread%desc)

            !call spec%disp%show("plang")
            !call spec%disp%show(spec%plang%val, format = format)
            !call spec%disp%note%show(spec%plang%desc)

            call spec%disp%show("randomSeed")
            call spec%disp%show(spec%randomSeed%val, format = format)
            call spec%disp%note%show(spec%randomSeed%desc)

            !call spec%disp%show("sysInfoFilePath")
            !call spec%disp%show(spec%sysInfoFilePath%val, format = format)
            !call spec%disp%note%show(spec%sysInfoFilePath%desc)

            call spec%disp%show("targetAcceptanceRate")
            call spec%disp%show(reshape(spec%targetAcceptanceRate%val, [2, 1]), format = format)
            call spec%disp%note%show(spec%targetAcceptanceRate%desc)

        end associate

    end subroutine report

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine sanitize(spec, err)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: sanitize
#endif
        use pm_err, only: err_type
        type(err_type), intent(inout) :: err
        class(specbase_type), intent(inout) :: spec
        character(*,SKG), parameter :: PROCEDURE_NAME = MODULE_NAME//SKG_"@sanitize()"

        ndim_block: block
            if (spec%ndim%val < 1_IK) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. &
                            &The input value for the specification `ndim` ("//spec%ndim%str//SKG_") must be a positive integer. &
                            &The specification `ndim` represents the number of dimensions of the domain of the objective function &
                            &to be explored via the user-specified exploration routine."
            end if
        end block ndim_block

        method_block: block
            if (.not. (spec%method%isParaDRAM .or. spec%method%isParaDISE .or. spec%method%isParaNest)) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Internal Error occurred. &
                            &The specified method name '"//spec%method%val//SKG_"' for the exploration routine does not match the internal list of method names."
            end if
        end block method_block

        description_block: block
            if (spec%description%val == spec%description%null) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. Internal library error detected. &
                            &The contents of `description` cannot be null-valued."
            end if
        end block description_block

        domain_block: block
            if (.not. (spec%domain%isBall .or. spec%domain%isCube)) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. &
                            &Invalid requested value for the sampler domain. &
                            &The input requested domain model ("//spec%domain%val//SKG_") is not supported. &
                            &The input specification domain cannot be set to anything other than '"//spec%domain%cube//SKG_"', or '"//spec%domain%ball//SKG_"'"
            end if
        end block domain_block

        domainAxisName_block: block
            character(len(spec%domainAxisName%val, IK),SKG) :: newval(size(spec%domainAxisName%val, 1, IK))
            logical(LK) :: redefined
            integer(IK) :: idim
            redefined = .false._LK
            do idim = 1, size(spec%domainAxisName%val, 1, IK)
                if (any(spec%domainAxisName%val(idim) == spec%domainAxisName%val(1 : idim - 1)) .or. any(spec%domainAxisName%val(idim) == spec%domainAxisName%val(idim + 1:))) then
                    newval(idim) = trim(adjustl(spec%domainAxisName%val(idim)))//getStr(idim)
                    redefined = .true._LK
                else
                    newval(idim) = spec%domainAxisName%val(idim)
                end if
            end do
            if (redefined) spec%domainAxisName%val(:) = newval
        end block domainAxisName_block

        domainBallAvg_block: block
            if (.not. all(-huge(0._RKG) < spec%domainBallAvg%val .and. spec%domainBallAvg%val < huge(0._RKG))) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. &
                            &The specified values for `domainBallAvg` must be all finite: "//getStr(spec%domainBallAvg%val)//SKG_"'"
            end if
        end block domainBallAvg_block

        domainBallCor_block: block
            use pm_matrixClass, only: posdefmat
            use pm_matrixClass, only: isMatClass
            !   There is no need to check for eyeness of the input correlation matrix. Only positive definiteness is enough.
            !   If the input correlation matrix is problematic, it will eventually lead to a non-positive-definite covariance matrix.
            if (.not. isMatClass(spec%domainBallCor%val, posdefmat)) then
                err%occurred = .true._LK ! This must be set only when .true.
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. The input requested `domainBallCor` for defining the &
                            &domain of the objective function passed to the sampler is not a positive-definite matrix: "//getStr(spec%domainBallCor%val)
            end if
        end block domainBallCor_block

        domainBallCov_block: block
            use pm_matrixInv, only: setMatInv, choUpp
            use pm_matrixTrace, only: getMatMulTraceLog
            use pm_matrixChol, only: setMatChol, transHerm, lowDia
            real(RKG), allocatable :: chol(:,:)
            integer(IK) :: info
            if (spec%domain%isFinite) then
                call setResized(chol, [spec%ndim%val, spec%ndim%val])
                call setResized(spec%domainBallCov%inv, [spec%ndim%val, spec%ndim%val])
                call setMatChol(spec%domainBallCov%val, lowDia, info, chol, transHerm)
                if (info /= 0_IK) then
                    err%occurred = .true._LK ! This must be set only when .true.
                    err%msg =   err%msg//&
                                PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. The input requested `domainBallCov` for defining the domain &
                                &of the objective function passed to the sampler is not a positive-definite matrix :"//getStr(spec%domainBallCov%val)
                end if
                spec%domain%logVol = getMatMulTraceLog(chol)
                call setMatInv(spec%domainBallCov%inv, chol, choUpp)
            end if
        end block domainBallCov_block

        domainBallStd_block: block
            integer(IK) :: idim
            do idim = 1, size(spec%domainBallStd%val, 1, IK)
                if (spec%domainBallStd%val(idim) <= 0._RKG) then
                    err%occurred = .true._LK
                    err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. &
                                &The input requested value ("//getStr(spec%domainBallStd%val(idim))//SKG_") for the component `"//getStr(idim)//&
                                SKG_"` of the input specification `domainBallStd` of the simulation must be a positive real number."
                end if
            end do
        end block domainBallStd_block

        domainCubeLimitLower_block: block
            integer(IK) :: idim
            if (.not. all(-huge(0._RKG) < spec%domainCubeLimitLower%val .and. spec%domainCubeLimitLower%val < huge(0._RKG))) then
                err%occurred = .true._LK
                err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. &
                &The specified values for `domainCubeLimitLower` must be all finite: "//getStr(spec%domainCubeLimitLower%val)//SKG_"'"
            end if
            ! domainCubeLimitLower must be all less than domainCubeLimitUpper.
            do idim = 1, size(spec%domainCubeLimitUpper%val(:), 1, IK)
                if (spec%domainCubeLimitUpper%val(idim) <= spec%domainCubeLimitLower%val(idim)) then
                    err%occurred = .true._LK
                    err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. &
                                &The input value for the upper limit of the component "//getStr(idim)//SKG_" of the &
                                    &variable domainCubeLimitUpper cannot be smaller than or equal to the input value &
                                    &for the lower limit of the corresponding dimension as given by domainCubeLimitLower:"//NL1//&
                                SKG_"    domainCubeLimitLower("//getStr(idim)//SKG_") = "//getStr(spec%domainCubeLimitLower%val(idim))//NL1//&
                                SKG_"    domainCubeLimitUpper("//getStr(idim)//SKG_") = "//getStr(spec%domainCubeLimitUpper%val(idim))
                end if
            end do
        end block domainCubeLimitLower_block

        domainCubeLimitUpper_block: block
            if (.not. all(-huge(0._RKG) < spec%domainCubeLimitUpper%val .and. spec%domainCubeLimitUpper%val < huge(0._RKG))) then
                err%occurred = .true._LK
                err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. &
                &The specified values for `domainCubeLimitUpper` must be all finite: "//getStr(spec%domainCubeLimitUpper%val)//SKG_"'"
            end if
        end block domainCubeLimitUpper_block

        domainErrCount_block: block
            if (spec%domainErrCount%val < 1_IK) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. &
                            &The input value for variable `domainErrCount` must be a positive integer. If you are unsure &
                            &about the appropriate value for this variable, simply drop it from the input. &
                            &The sampler will automatically assign an appropriate value to it."
            end if
        end block domainErrCount_block

        domainErrCountMax_block: block
            if (spec%domainErrCountMax%val < 1_IK) then
                err%occurred = .true._LK
                err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. &
                &The input value for variable `domainErrCountMax` must be a positive integer. &
                &If you are unsure about the appropriate value for this variable, simply drop it from the input. &
                &The sampler will automatically assign an appropriate value to it."
            end if
        end block domainErrCountMax_block

        !inputFileHasPriority_block: block
        !end block inputFileHasPriority_block

        outputChainFileFormat_block: block
            if (.not.(spec%outputChainFileFormat%isCompact .or. spec%outputChainFileFormat%isVerbose .or. spec%outputChainFileFormat%isBinary)) then
                err%occurred = .true._LK
                err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. The input requested chain file format ('"//spec%outputChainFileFormat%val//&
                SKG_"') contained within the simulation specification `outputChainFileFormat` cannot be anything other than '"//&
                spec%outputChainFileFormat%compact//SKG_"' or '"//spec%outputChainFileFormat%verbose//SKG_"' or '"//spec%outputChainFileFormat%binary//SKG_"'. If you do not &
                &know an appropriate value for outputChainFileFormat, drop it from the input list. &
                &The sampler will automatically assign an appropriate value to it."
            end if
        end block outputChainFileFormat_block

        outputColumnWidth_block: block
            if (spec%outputColumnWidth%val < 0_IK) then
                err%occurred = .true._LK
                err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. &
                &The input value for variable `outputColumnWidth` must be a non-negative integer. &
                &If you are unsure about the appropriate value for this variable, simply drop it from the input. &
                &The sampler will automatically assign an appropriate value to it."
            elseif (0_IK < spec%outputColumnWidth%val .and. spec%outputColumnWidth%val < spec%outputPrecision%val + 7_IK) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. The input value for variable &
                            &`outputColumnWidth` must be equal to or greater than the input value for `outputPrecision + 7`. &
                            &If you are unsure about the appropriate value for this variable, either set it to zero on input, &
                            &or simply drop it from the input. The sampler will automatically assign an appropriate value to it."
            end if
        end block outputColumnWidth_block

        outputFileName_block: block
            ! Ideally there should be a test here to ensure the prefix contains only valid characters.
            if (spec%outputFileName%val == SKG_"") then
                err%occurred = .true._LK
                err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. &
                &The input value for variable outputFileName cannot be empty."
            end if
        end block outputFileName_block

        outputStatus_block: block
            if (.not. (spec%outputStatus%is%extend .or. spec%outputStatus%is%retry)) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. The specified chain file format ('"//spec%outputStatus%val//&
                            SKG_"') via `outputStatus` cannot be anything other than 'retry', 'extend', 'copy-retry', or 'copy-extend'. &
                            &If you do not know an appropriate value for `outputStatus`, drop it from the input list. &
                            &The sampler will automatically assign an appropriate value to it."
            end if
        end block outputStatus_block

        outputPrecision_block: block
            if (spec%outputPrecision%val < 1_IK) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. &
                            &The input value for variable `outputPrecision` must be a positive integer. &
                            &If you are unsure about the appropriate value for this variable, simply drop it from the input. &
                            &The sampler will automatically assign an appropriate value to it."
            end if
        end block outputPrecision_block

        outputReportPeriod_block: block
            if (spec%outputReportPeriod%val < 1_IK) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. &
                            &The input value for variable outputReportPeriod must be a positive integer value. &
                            &If you are unsure about the appropriate value for this variable, simply drop it from the input. &
                            &The sampler will automatically assign an appropriate value to it."
            else
                spec%outputReportPeriod%inv = 1._RKG / spec%outputReportPeriod%val
            end if
        end block outputReportPeriod_block

        outputRestartFileFormat_block: block
            if (.not. (spec%outputRestartFileFormat%isBinary .or. spec%outputRestartFileFormat%isAscii)) then
                err%occurred = .true._LK
                err%msg = err%msg//PROCEDURE_NAME//SKG_"@sanitize(): Error occurred. &
                &The input requested restart file format ('"//spec%outputRestartFileFormat%val//&
                SKG_"') represented by the variable `outputRestartFileFormat` cannot be anything other than '"//&
                spec%outputRestartFileFormat%binary//SKG_"' or '"//spec%outputRestartFileFormat%ascii//SKG_"'. &
                &If you are unsure about an appropriate value for `outputRestartFileFormat`, drop it from the input list. &
                &The sampler will automatically assign an appropriate value to it."
            end if
        end block outputRestartFileFormat_block

        outputSampleSize_block: block
            if (.not. abs(spec%outputSampleSize%val) < huge(0_IK)) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. &
                            &The input value for `outputSampleSize` cannot be larger than the largest representable integer. &
                            &If you are unsure of an appropriate value for this variable, simply drop it from the input. &
                            &The sampler will automatically assign an appropriate value to it."
            end if
        end block outputSampleSize_block

        outputSeparator_block: block
            use pm_strASCII, only: isCharDigit
            character(:,SKG), allocatable :: separator
            integer(IK) :: separatorLen, i
            separator = trim(adjustl(spec%outputSeparator%val))
            separatorLen = len(separator, IK)
            do i = 1_IK, separatorLen
                if (isCharDigit(separator(i:i)) .or. separator(i:i) == SKG_"." .or. separator(i:i) == SKG_"-" .or. separator(i:i) == SKG_"+") then
                    err%occurred = .true._LK
                    err%msg =   err%msg//NL2//PROCEDURE_NAME//SKG_"@sanitize(): Error occurred. &
                                &The input value for variable outputSeparator cannot contain any digits or the period symbol '.' or '-' &
                                &or '+'. If you are unsure about the appropriate value for this variable, simply drop it from the input. &
                                &The sampler will automatically assign an appropriate value to it."
                    exit
                end if
            end do
        end block outputSeparator_block

        outputSplashMode_block: block
            if (.not. (spec%outputSplashMode%is%normal .or. spec%outputSplashMode%is%silent .or. spec%outputSplashMode%is%quiet)) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//SKG_"@sanitize(): Error occurred. &
                            &The input requested splash mode ('"//spec%outputSplashMode%val//&
                            SKG_"') represented by the variable `outputSplashMode` cannot be anything other than '"//&
                            spec%outputSplashMode%normal//SKG_"' or '"//spec%outputSplashMode%quiet//SKG_"' or '"//spec%outputSplashMode%silent//SKG_"'. &
                            &If you are unsure about an appropriate value for `outputSplashMode`, drop it from the input list. &
                            &The sampler will automatically assign an appropriate value to it."
            end if
        end block outputSplashMode_block

        parallelism_block: block
            if (.not.(spec%parallelism%is%singleChain .or. spec%parallelism%is%multiChain)) then
                err%occurred = .true._LK
                err%msg = err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. &
                &The input requested parallelization method ("//spec%parallelism%val//&
                SKG_") represented by variable `parallelism` cannot be anything other than &
                &'"//spec%parallelism%singleChain//SKG_"' or '"//spec%parallelism%multiChain//SKG_"'. If you are unsure about &
                &the appropriate value for `parallelism`, drop it from the input list. &
                &The sampler will automatically assign an appropriate value to it."
            end if
            spec%parallelism%is%forkJoin = spec%parallelism%is%singleChain .and. 1_IK < spec%image%count
        end block parallelism_block

        parallelismMpiFinalizeEnabled_block: block
        end block parallelismMpiFinalizeEnabled_block

        parallelismNumThread_block: block
#if         OMP_ENABLED
            if (spec%parallelismNumThread%val < 0_IK) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. &
                            &The input value for the specification `parallelismNumThread` ("//getStr(spec%parallelismNumThread%val)//SKG_") must be a non-negative integer. &
                            &The specification `parallelismNumThread` represents the number of shared-memory CPU threads over which multiple parallel evaluations of the &
                            &objective function must be done. This specification is relevant only to OpenMP-enabled parallel simulations in C, C++, and Fortran, &
                            &and shared-memory non-MPI parallel simulations in higher-level programming languages such as MATLAB, Python, or R. &
                            &Specifying `0` leads to using all available CPU threads for the requested simulation task."
            end if
#endif
        end block parallelismNumThread_block

        !plang_block: block
        !end block plang_block

        randomSeed_block: block
            if (spec%randomSeed%val < 1_IK) then
                err%occurred = .true._LK
                err%msg = err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. &
                &The input value for `randomSeed` ("//getStr(spec%randomSeed%val)//SKG_") must be a positive integer. &
                &If you are unsure about the appropriate value for this variable, simply drop it from the input. &
                &The sampler will automatically assign an appropriate value to it."
            end if
        end block randomSeed_block

        sysInfoFilePath_block: block
        end block sysInfoFilePath_block

        targetAcceptanceRate_block: block
            if (any(spec%targetAcceptanceRate%val < 0._RKG) .or. any(1._RKG < spec%targetAcceptanceRate%val)) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. The target acceptance ratio limits `targetAcceptanceRate` ["//&
                            getStr(spec%targetAcceptanceRate%val(1))//SKG_", "//getStr(spec%targetAcceptanceRate%val(2))//SKG_"] cannot be less than 0 or larger than 1."
            end if
            if (all(spec%targetAcceptanceRate%val == 0._RKG) .or. all(spec%targetAcceptanceRate%val == 1._RKG)) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. The target acceptance ratio limits targetAcceptanceRate ["//&
                            getStr(spec%targetAcceptanceRate%val(1))//SKG_", "//getStr(spec%targetAcceptanceRate%val(2))//SKG_"] cannot be both 0 or both 1."
            end if
            if (spec%targetAcceptanceRate%val(2) < spec%targetAcceptanceRate%val(1)) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKG_": Error occurred. The the lower limit of the input specification targetAcceptanceRate ["//&
                            getStr(spec%targetAcceptanceRate%val(1))//SKG_","//getStr(spec%targetAcceptanceRate%val(2))//SKG_"] cannot be larger than the specified upper limit."
            end if
        end block targetAcceptanceRate_block

        if (spec%domain%isCube .and. spec%domain%isFinite) spec%domain%logVol = sum(log(spec%domainCubeLimitUpper%val - spec%domainCubeLimitLower%val))

    end subroutine sanitize

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SHARED
