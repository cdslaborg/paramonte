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

    use pm_strASCII, only: SUB
    use pm_except, only: isNAN
    use pm_except, only: setNAN
    use pm_val2str, only: getStr
    use pm_strASCII, only: getStrLower
    use pm_kind, only: SKC => SK, SK, IK, LK
    use pm_sampling_base, only: specbase_type, astatbase_type, sfcbase_type, NL2, NL1
    use pm_sampling_scio, only: cfcnest_type

    implicit none

    character(*,SKC)    , parameter         :: MODULE_NAME = SK_"@pm_sampling_nest"

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! simulation declarations.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! type for output sample file contents.
    type, extends(sfcbase_type)             :: sfcnest_type
    end type

    type, abstract, extends(astatbase_type) :: astatnest_type
    end type

    type, abstract, extends(astatnest_type) :: statnest_type
        type(cfcnest_type)                  :: cfc
        type(sfcnest_type)                  :: sfc
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! specification declarations.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type                                :: domainPartitionAdaptationCount_type
        integer(IK)                     :: val
        integer(IK)                     :: def
        integer(IK)                     :: null
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainPartitionAdaptationPeriod_type
        integer(IK)                     :: val
        integer(IK)                     :: def
        integer(IK)                     :: null
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainPartitionBiasCorrectionEnabled_type
        logical(LK)                     :: val
        logical(LK)                     :: def
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainPartitionCountMax_type
        integer(IK)                     :: val
        integer(IK)                     :: def
        integer(IK)                     :: null
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainPartitionFactorExpansion_type
        real(RKC)                       :: val
        real(RKC)                       :: def
       !real(RKC)                       :: null
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainPartitionFactorShrinkage_type
        real(RKC)                       :: val
        real(RKC)                       :: def
       !real(RKC)                       :: null
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainPartitionKmeansClusterCountMax_type
        integer(IK)                     :: val
        integer(IK)                     :: def
        integer(IK)                     :: null
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainPartitionKmeansClusterSizeMin_type
        integer(IK)                     :: val
        integer(IK)                     :: def
        integer(IK)                     :: null
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainPartitionKmeansNormalizationEnabled_type
        logical(LK)                     :: val
        logical(LK)                     :: def
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainPartitionKmeansNumFailMax_type
        integer(IK)                     :: val
        integer(IK)                     :: def
        integer(IK)                     :: null
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainPartitionKmeansNumRecursionMax_type
        integer(IK)                     :: val
        integer(IK)                     :: def
        integer(IK)                     :: null
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainPartitionKmeansNumTry_type
        integer(IK)                     :: val
        integer(IK)                     :: def
        integer(IK)                     :: null
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainPartitionKvolumeNumRecursionMax_type
        integer(IK)                     :: val
        integer(IK)                     :: def
        integer(IK)                     :: null
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainPartitionKvolumeWeightExponent_type
        real(RKC)                       :: val
        real(RKC)                       :: def
       !real(RKC)                       :: null
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainPartitionMethod_type
        logical(LK)                     :: isMultiNest = .false._LK
        logical(LK)                     :: isDynesty = .false._LK
        logical(LK)                     :: isMinVol = .false._LK
        logical(LK)                     :: isMaxDen = .false._LK
        character(9,SKC)                :: multinest = SKC_"multinest"
        character(7,SKC)                :: dynesty = SKC_"dynesty"
        character(7,SKC)                :: minvol = SKC_"minvol"
        character(7,SKC)                :: maxden = SKC_"maxden"
        character(:,SKC), allocatable   :: val
        character(:,SKC), allocatable   :: def
        character(:,SKC), allocatable   :: null
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainPartitionObject_type
        logical(LK)                     :: isBall = .false._LK
        character(:,SKC), allocatable   :: val
        character(:,SKC), allocatable   :: def
        character(:,SKC), allocatable   :: null
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainPartitionOptimizationScaleEnabled_type
        logical(LK)                     :: val
        logical(LK)                     :: def
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainPartitionOptimizationShapeEnabled_type
        logical(LK)                     :: val
        logical(LK)                     :: def
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainPartitionOptimizationShapeScaleEnabled_type
        logical(LK)                     :: val
        logical(LK)                     :: def
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainPartitionScaleFactor_type
        real(RKC)                       :: val
        real(RKC)                       :: def
       !real(RKC)                       :: null
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: domainSampler_type
        logical(LK)                     :: isRejection
        character(:,SKC), allocatable   :: val
        character(:,SKC), allocatable   :: def
        character(:,SKC), allocatable   :: null
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: liveSampleSize_type
        real(RKC)                       :: log
        integer(IK)                     :: val
        integer(IK)                     :: def
        integer(IK)                     :: null
        character(:,SKC), allocatable   :: desc
    end type

    type                                :: tolerance_type
        real(RKC)                       :: val
        real(RKC)                       :: def
       !real(RKC)                       :: null
        character(:,SKC), allocatable   :: desc
    end type

    type, extends(specbase_type)                                    :: specnest_type
        type(domainPartitionAdaptationCount_type                )   :: domainPartitionAdaptationCount
        type(domainPartitionAdaptationPeriod_type               )   :: domainPartitionAdaptationPeriod
        type(domainPartitionBiasCorrectionEnabled_type          )   :: domainPartitionBiasCorrectionEnabled
        type(domainPartitionCountMax_type                       )   :: domainPartitionCountMax
        type(domainPartitionFactorExpansion_type                )   :: domainPartitionFactorExpansion
        type(domainPartitionFactorShrinkage_type                )   :: domainPartitionFactorShrinkage
        type(domainPartitionKmeansClusterCountMax_type          )   :: domainPartitionKmeansClusterCountMax
        type(domainPartitionKmeansClusterSizeMin_type           )   :: domainPartitionKmeansClusterSizeMin
        type(domainPartitionKmeansNormalizationEnabled_type     )   :: domainPartitionKmeansNormalizationEnabled
        type(domainPartitionKmeansNumFailMax_type               )   :: domainPartitionKmeansNumFailMax
        type(domainPartitionKmeansNumRecursionMax_type          )   :: domainPartitionKmeansNumRecursionMax
        type(domainPartitionKmeansNumTry_type                   )   :: domainPartitionKmeansNumTry
        type(domainPartitionKvolumeNumRecursionMax_type         )   :: domainPartitionKvolumeNumRecursionMax
        type(domainPartitionKvolumeWeightExponent_type          )   :: domainPartitionKvolumeWeightExponent
        type(domainPartitionMethod_type                         )   :: domainPartitionMethod
        type(domainPartitionObject_type                         )   :: domainPartitionObject
        type(domainPartitionOptimizationScaleEnabled_type       )   :: domainPartitionOptimizationScaleEnabled
        type(domainPartitionOptimizationShapeEnabled_type       )   :: domainPartitionOptimizationShapeEnabled
        type(domainPartitionOptimizationShapeScaleEnabled_type  )   :: domainPartitionOptimizationShapeScaleEnabled
        type(domainPartitionScaleFactor_type                    )   :: domainPartitionScaleFactor
        type(domainSampler_type                                 )   :: domainSampler
        type(liveSampleSize_type                                )   :: liveSampleSize
        type(tolerance_type                                     )   :: tolerance
    contains
        procedure, pass, private                                    :: sanitize
        procedure, pass, private                                    :: report
        procedure, pass, public                                     :: set
    end type

    interface specnest_type
        module procedure :: construct
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine killMeAlreadyCMake1_RK5(); use pm_sampling_scio_RK5, only: RKC; end subroutine
    subroutine killMeAlreadyCMake1_RK4(); use pm_sampling_scio_RK4, only: RKC; end subroutine
    subroutine killMeAlreadyCMake1_RK3(); use pm_sampling_scio_RK3, only: RKC; end subroutine
    subroutine killMeAlreadyCMake1_RK2(); use pm_sampling_scio_RK2, only: RKC; end subroutine
    subroutine killMeAlreadyCMake1_RK1(); use pm_sampling_scio_RK1, only: RKC; end subroutine

    subroutine killMeAlreadyCMake2_RK5(); use pm_sampling_base_RK5, only: RKC; end subroutine
    subroutine killMeAlreadyCMake2_RK4(); use pm_sampling_base_RK4, only: RKC; end subroutine
    subroutine killMeAlreadyCMake2_RK3(); use pm_sampling_base_RK3, only: RKC; end subroutine
    subroutine killMeAlreadyCMake2_RK2(); use pm_sampling_base_RK2, only: RKC; end subroutine
    subroutine killMeAlreadyCMake2_RK1(); use pm_sampling_base_RK1, only: RKC; end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function construct(modelr, method, ndim) result(spec)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct
#endif
        use pm_kind, only: modelr_type
        type(modelr_type), intent(in) :: modelr
        character(*,SKC), intent(in) :: method
        integer(IK), intent(in) :: ndim
        type(specnest_type) :: spec

        spec%specbase_type = specbase_type(modelr, method, ndim)

        domainPartitionAdaptationCount_block: block
            use pm_sampling_scio, only: domainPartitionAdaptationCount
            spec%domainPartitionAdaptationCount%def = huge(0_IK)
            spec%domainPartitionAdaptationCount%null = -huge(0_IK)
            spec%domainPartitionAdaptationCount%desc = &
            SKC_"The simulation specification `domainPartitionAdaptationCount` is a positive-valued scalar of type `integer` of default kind representing the total number of adaptive &
                &updates to make to the parameters of the domain sampler to increase the efficiency of the sampler thus increasing the overall sampling efficiency. &
                &Every `domainPartitionAdaptationPeriod` number of calls to the objective function, the parameters of the domain sampler &
                &will be updated until either the total number of adaptive updates reaches the value of `domainPartitionAdaptationCount`. &
                &If the condition `domainPartitionAdaptationCount == 0`, then the domain sampler parameters will be fixed to &
                &the initial input values throughout the entire sampling/integration. &
                &The default value is "//getStr(spec%domainPartitionAdaptationCount%def)//SKC_"."
            !$omp master
            domainPartitionAdaptationCount = spec%domainPartitionAdaptationCount%null
            !$omp end master
        end block domainPartitionAdaptationCount_block

        domainPartitionAdaptationPeriod_block: block
            use pm_sampling_scio, only: domainPartitionAdaptationPeriod
            spec%domainPartitionAdaptationPeriod%def = 500_IK ! + 1_IK ! max(nd+1_IK,100_IK)
            spec%domainPartitionAdaptationPeriod%null = -huge(0_IK)
            spec%domainPartitionAdaptationPeriod%desc = &
            SKC_"The simulation specification `domainPartitionAdaptationPeriod` is a positive-valued scalar of type `integer` of default kind. &
                &Every `domainPartitionAdaptationPeriod` calls to the objective function, the approximated shape of the domain of the objective function will be updated. &
                &The smaller the value of `domainPartitionAdaptationPeriod`, the easier it will be for the integrator to adapt the domain sampler to the constrained &
                &shape of the domain of the objective function. However, this will happen at the expense of slower simulation runtime because the adaptation &
                &process is generally computationally expensive, in particular, for very high dimensional objective functions (`ndim >> 1`). &
                &The larger the value of `domainPartitionAdaptationPeriod`, the easier it will be for the integrator kernel to keep the &
                &sampling efficiency close to the requested target acceptance rate range (specified via the `targetAcceptanceRate`). &
                &However, too large values for `domainPartitionAdaptationPeriod` will only delay the adaptation of the domain &
                &sampler to the global structure of the domain of the objective function that is being sampled. &
                &The default value is "//getStr(spec%domainPartitionAdaptationPeriod%def)//SKC_"."
            !$omp master
            domainPartitionAdaptationPeriod = spec%domainPartitionAdaptationPeriod%null
            !$omp end master
        end block domainPartitionAdaptationPeriod_block

        domainPartitionBiasCorrectionEnabled_block: block
            use pm_sampling_scio, only: domainPartitionBiasCorrectionEnabled
            spec%domainPartitionBiasCorrectionEnabled%def = .true._LK
            spec%domainPartitionBiasCorrectionEnabled%desc = &
            SKC_"The simulation specification `domainPartitionBiasCorrectionEnabled` is a scalar of type `logical` (Boolean) of default kind. &
                &When `domainPartitionBiasCorrectionEnabled` is set to the logical true value (equivalent to `.true.`, `.t.`, or `true`, all case-INsensitive, &
                &in an external input file, or True from within Python environment or true from within a MATLAB session), the final scales of the partitions of &
                &the live points of the sampler will be optimized. This will minimize the inaccuracies in the partitioning of the domain of the objective &
                &function, although, it comes with an extra cost in terms of sampler runtime efficiency. The default value is "//getStr(spec%domainPartitionBiasCorrectionEnabled%def)//SKC_"."
            !$omp master
            domainPartitionBiasCorrectionEnabled = spec%domainPartitionBiasCorrectionEnabled%def
            !$omp end master
        end block domainPartitionBiasCorrectionEnabled_block

        domainPartitionCountMax_block: block
            use pm_sampling_scio, only: domainPartitionCountMax
            spec%domainPartitionCountMax%def  = huge(0_IK)
            spec%domainPartitionCountMax%null = -huge(0_IK)
            spec%domainPartitionCountMax%desc = &
            SKC_"The simulation specification `domainPartitionCountMax` is a positive-valued scalar of type `integer` of default kind, &
                &representing the maximum number of partitions to be used in constraining the domain of the objective function as defined &
                &by the live (active) points at any stage of the simulation. When the condition `domainPartitionCountMax == 1` holds, &
                &a single partition (e.g., ellipsoid) will be used to constrain the domain of the objective function. &
                &Note that the condition `1 =< domainPartitionCountMax =< liveSampleSize / (ndim + 1)` must hold in all simulations where `ndim` &
                &is the number of dimensions of the domain of the density function to be sampled. Do not set this simulation specification to small &
                &numbers unless you know the consequences. Doing so, may cause the sampler to become highly inefficient. The default value is "//getStr(spec%domainPartitionCountMax%def)//SKC_"."
            !$omp master
            domainPartitionCountMax = spec%domainPartitionCountMax%null
            !$omp end master
        end block domainPartitionCountMax_block

        domainPartitionFactorExpansion_block: block
            use pm_sampling_scio, only: domainPartitionFactorExpansion
            spec%domainPartitionFactorExpansion%def = 1._RKC !1.1_RKC
            spec%domainPartitionFactorExpansion%desc = &
            SKC_"The simulation specification `domainPartitionFactorExpansion` is a positive-valued scalar of type `real` of the highest precision &
                &available within the ParaMonte library, representing the factor by which the volume of the bounding regions of the live (active) &
                &points are enlarged after the partitioning of the domain of the objective function. &
                &Setting `domainPartitionFactorExpansion << 1` may lead to biased and unreliable results. &
                &Conversely, setting `domainPartitionFactorExpansion >> 1` will lead to conservative inefficient &
                &partitioning of the domain and slow sampling and integration of the objective function. &
                &The default value is "//getStr(spec%domainPartitionFactorExpansion%def)//SKC_"."
                !&If ellipsoidalBoundingVolume < domainPartitionFactorExpansion * actualLivePointVolume, then, further partitioning of the live &
                !&points will be performed. Conversely, if ellipsoidalBoundingVolume >= domainPartitionFactorExpansion * actualLivePointVolume, &
                !&then, no further partitioning of the live points will be performed and the estimated ellipsoidalBoundingVolume will be enlarged &
                !&to match domainPartitionFactorExpansion * actualLivePointVolume.
            !$omp master
            call setNAN(domainPartitionFactorExpansion)
            !$omp end master
        end block domainPartitionFactorExpansion_block

        domainPartitionFactorShrinkage_block: block
            use pm_sampling_scio, only: domainPartitionFactorShrinkage
            spec%domainPartitionFactorShrinkage%def = 1._RKC !1.1_RKC
            spec%domainPartitionFactorShrinkage%desc = &
            SKC_"The simulation specification `domainPartitionFactorShrinkage` is a positive-valued scalar of type `real` of the highest precision &
                &available within the ParaMonte library, representing the factor by which the volume of the bounding regions of the live (active) &
                &points is multiplied before its usage in ellipsoidal partitioning. If the condition &
                &`ellipsoidalBoundingVolume < domainPartitionFactorShrinkage * actualLivePointVolume` holds, then, further partitioning of the live points &
                &will be performed. Conversely, if the condition `ellipsoidalBoundingVolume >= domainPartitionFactorShrinkage * actualLivePointVolume` holds, &
                & no further partitioning of the live points will be performed and the estimated `ellipsoidalBoundingVolume` will be enlarged &
                &to match `domainPartitionFactorShrinkage * actualLivePointVolume`. This is a low-level simulation specification. &
                &Setting `domainPartitionFactorShrinkage << 1` may lead to biased and unreliable results. &
                &Conversely, setting `domainPartitionFactorShrinkage >> 1` will lead to conservative inefficient partitioning of &
                &the domain and slow sampling and integration of the objective function. The default value is "//getStr(spec%domainPartitionFactorShrinkage%def)//SKC_"."
            !$omp master
            call setNAN(domainPartitionFactorShrinkage)
            !$omp end master
        end block domainPartitionFactorShrinkage_block

        domainPartitionKmeansClusterCountMax_block: block
            use pm_sampling_scio, only: domainPartitionKmeansClusterCountMax
            spec%domainPartitionKmeansClusterCountMax%def  = 2_IK ! 10000_IK
            spec%domainPartitionKmeansClusterCountMax%null = -huge(0_IK)
            spec%domainPartitionKmeansClusterCountMax%desc = &
            SKC_"The simulation specification `domainPartitionKmeansClusterCountMax` is a scalar positive `integer` of default kind that sets &
                &the maximum number of clusters to be found by the Kmeans clustering algorithm at any level of partitioning of the live sample. &
                &Note that the condition `1 <= domainPartitionKmeansClusterCountMax <= domainPartitionCountMax` must hold. &
                &The default value is `min("//getStr(spec%domainPartitionKmeansClusterCountMax%def)//SKC_", domainPartitionCountMax)`."
            !$omp master
            domainPartitionKmeansClusterCountMax = spec%domainPartitionKmeansClusterCountMax%null
            !$omp end master
        end block domainPartitionKmeansClusterCountMax_block

        domainPartitionKmeansClusterSizeMin_block: block
            use pm_sampling_scio, only: domainPartitionKmeansClusterSizeMin
            spec%domainPartitionKmeansClusterSizeMin%def  = 0_IK
            spec%domainPartitionKmeansClusterSizeMin%null = -huge(0_IK)
            spec%domainPartitionKmeansClusterSizeMin%desc = &
            SKC_"The simulation specification `domainPartitionKmeansClusterSizeMin` in a scalar non-negative `integer` of default kind representing the minimum &
            &size of the clusters allowed to be identified by the Kmeans clustering algorithms when partitioning any subsets of the live sample of the sampler. &
            &Note that the condition `0 <= domainPartitionKmeansClusterSizeMin <= liveSampleSize` must hold. Do not set this simulation specification to values &
            &significantly larger than `ndim + 1` where `ndim` is the number of dimensions of the domain of the objective function. Doing so may cause the &
            &sampler to become highly inefficient. The default value is "//getStr(spec%domainPartitionKmeansClusterSizeMin%def)//SKC_"."
            !$omp master
            domainPartitionKmeansClusterSizeMin = spec%domainPartitionKmeansClusterSizeMin%null
            !$omp end master
        end block domainPartitionKmeansClusterSizeMin_block

        domainPartitionKmeansNormalizationEnabled_block: block
            use pm_sampling_scio, only: domainPartitionKmeansNormalizationEnabled
            spec%domainPartitionKmeansNormalizationEnabled%def = .true._LK
            spec%domainPartitionKmeansNormalizationEnabled%desc = &
            SKC_"The simulation specification `domainPartitionKmeansNormalizationEnabled` is a scalar of type `logical` (Boolean) of default kind. &
                &When `domainPartitionKmeansNormalizationEnabled` is set to the logical (Boolean) true value (equivalent to `.true.`, `.t.`, or `true`, all case-INsensitive, &
                &from inside the input file or `True` from inside Python or `true` from inside MATLAB), any subset of live points will be standardized before being passed &
                &to the Kmeans algorithm. Enabling this option will lead to better behavior of the Kmeans algorithm and potentially better &
                &identification of the subclusters in data. This improvement, however, comes with an extra cost in terms of runtime &
                &efficiency of the "//spec%method%val//SKC_" sampler. The default value is "//getStr(spec%domainPartitionKmeansNormalizationEnabled%def)//SKC_"."
            !$omp master
            domainPartitionKmeansNormalizationEnabled = spec%domainPartitionKmeansNormalizationEnabled%def
            !$omp end master
        end block domainPartitionKmeansNormalizationEnabled_block

        domainPartitionKmeansNumFailMax_block: block
            use pm_sampling_scio, only: domainPartitionKmeansNumFailMax
            spec%domainPartitionKmeansNumFailMax%def  = 10_IK ! 10000_IK
            spec%domainPartitionKmeansNumFailMax%null = -huge(0_IK)
            spec%domainPartitionKmeansNumFailMax%desc = &
            SKC_"The simulation specification `domainPartitionKmeansNumFailMax` is a positive-valued scalar of type `integer` of default kind &
                &that sets the maximum allowed number of Kmeans clustering failures at each clustering try within simulations that employ the rejection sampling method. &
                &This is a low-level simulation specification whose value should not be changed from the default unless the consequences are well understood. &
                &Setting this simulation property to a small value, for example, on the order of ten, may lead to highly inefficient simulations. &
                &The default value is "//getStr(spec%domainPartitionKmeansNumFailMax%def)//SKC_"."
            !$omp master
            domainPartitionKmeansNumFailMax = spec%domainPartitionKmeansNumFailMax%null
            !$omp end master
        end block domainPartitionKmeansNumFailMax_block

        domainPartitionKmeansNumRecursionMax_block: block
            use pm_sampling_scio, only: domainPartitionKmeansNumRecursionMax
            spec%domainPartitionKmeansNumRecursionMax%def  = 10_IK ! 10000_IK
            spec%domainPartitionKmeansNumRecursionMax%null = -huge(0_IK)
            spec%domainPartitionKmeansNumRecursionMax%desc = &
            SKC_"The simulation specification `domainPartitionKmeansNumRecursionMax` is a positive scalar of type `integer` of default kind &
                &that sets the maximum allowed number of Kmeans convergence iterations in simulations that employ the rejection sampling method. &
                &This is a low-level simulation specification whose value should not be changed from the default unless the consequences are well understood. &
                &Setting this simulation property to a small value, for example, on the order of ten, may lead to highly inefficient "//spec%method%val//SKC_" simulations. &
                &The default value is "//getStr(spec%domainPartitionKmeansNumRecursionMax%def)//SKC_"."
            !$omp master
            domainPartitionKmeansNumRecursionMax = spec%domainPartitionKmeansNumRecursionMax%null
            !$omp end master
        end block domainPartitionKmeansNumRecursionMax_block

        domainPartitionKmeansNumTry_block: block
            use pm_sampling_scio, only: domainPartitionKmeansNumTry
            spec%domainPartitionKmeansNumTry%def  = 10_IK ! 10000_IK
            spec%domainPartitionKmeansNumTry%null = -huge(0_IK)
            spec%domainPartitionKmeansNumTry%desc = &
            SKC_"The simulation specification `domainPartitionKmeansNumTry` is a positive scalar of type `integer` of default kind &
                &that sets the number of Kmeans clustering tries at each level of partitioning of the set of live points during the simulation. &
                &This simulation specification should be not be set to very large values (e.g., >> 10) because the simulation will significantly &
                &slow down. The default value is "//getStr(spec%domainPartitionKmeansNumTry%def)//SKC_"."
            !$omp master
            domainPartitionKmeansNumTry = spec%domainPartitionKmeansNumTry%null
            !$omp end master
        end block domainPartitionKmeansNumTry_block

        domainPartitionKvolumeNumRecursionMax_block: block
            use pm_sampling_scio, only: domainPartitionKvolumeNumRecursionMax
            spec%domainPartitionKvolumeNumRecursionMax%def  = 300_IK
            spec%domainPartitionKvolumeNumRecursionMax%null = -huge(0_IK)
            spec%domainPartitionKvolumeNumRecursionMax%desc = &
            "The simulation specification `domainPartitionKvolumeNumRecursionMax` is a positive scalar of type `integer` of default kind &
            &that sets the maximum allowed number of volume minimization convergence iterations in simulations that employ the rejection sampling method. &
            &This is a low-level simulation specification whose value should not be changed from the default unless the consequences are well understood. &
            &Setting this simulation property to a small value, for example, on the order of ten, may lead to highly inefficient "//spec%method%val//SKC_" simulations. &
            &The default value is "//getStr(spec%domainPartitionKvolumeNumRecursionMax%def)//SKC_"."
            !$omp master
            domainPartitionKvolumeNumRecursionMax = spec%domainPartitionKvolumeNumRecursionMax%null
            !$omp end master
        end block domainPartitionKvolumeNumRecursionMax_block

        domainPartitionKvolumeWeightExponent_block: block
            use pm_sampling_scio, only: domainPartitionKvolumeWeightExponent
            spec%domainPartitionKvolumeWeightExponent%def = 0._RKC
           !spec%domainPartitionKvolumeWeightExponent%null = NULL_RKC
            spec%domainPartitionKvolumeWeightExponent%desc = &
            SKC_"The simulation specification `domainPartitionKvolumeWeightExponent` is a positive scalar of type `real` of the highest precision available &
                &within the ParaMonte library, representing the exponent of the power-law density weights `domainPartitionKvolumeWeightExponent` weights of &
                &the Mahalanobis distances of live (active) points from their corresponding cluster centers when the rejection sampling method is employed. &
                &Reasonable values include `0` or `1`. This is a low-level simulation specification whose value should not be changed from the default &
                &unless the consequences are well understood. Setting this simulation property to any values other than the default value can drastically &
                &affect the performance or even convergence and success of the simulations when the rejection sampling method is used. &
                &The default value for `domainPartitionKvolumeWeightExponent` is "//getStr(spec%domainPartitionKvolumeWeightExponent%def)//SKC_"."
            !$omp master
            call setNAN(domainPartitionKvolumeWeightExponent)
            !$omp end master
        end block domainPartitionKvolumeWeightExponent_block

        domainPartitionMethod_block: block
            use pm_sampling_scio, only: domainPartitionMethod
            spec%domainPartitionMethod%def = "ParaMonte-OptDen"
            spec%domainPartitionMethod%null = repeat(SUB, len(domainPartitionMethod))
            spec%domainPartitionMethod%desc = &
            SKC_"The simulation specification `domainPartitionMethod` is a scalar string of default kind of maximum length "//getStr(domainPartitionMethod)//SKC_" containing &
                &the name of the method to be used for the domain partitioning within the Nested sampler. The string value must be enclosed by either single or double &
                &quotation marks when provided as input. Options that are currently supported include:"//NL2//&
            SKC_"+   domainPartitionMethod = 'ParaMonte' or 'ParaMonte optimal', all case-INsensitive."//NL2//&
            SKC_"    This is equivalent to the constrained rejection sampling via partitions &
                     &formed by bounding objects around the set of live (active) points using a complex method that optimizes the shapes &
                     &and count of bounding objects for the most-likely (optimal) density of live points within each bounding object."//NL2//&
            SKC_"+   domainPartitionMethod = 'ParaMonte MinVol' or 'ParaMonte Minimum Volume', all case-INsensitive."//NL2//&
            SKC_"    This is equivalent to the constrained rejection sampling via partitions &
                     &formed by bounding objects around the set of live (active) points using a complex method that optimizes the shapes &
                     &and count of the bounding objects to achieve a minimum overall domain volume."//NL2//&
            SKC_"+   domainPartitionMethod = 'MultiNest', all case-INsensitive."//NL2//&
            SKC_"    in which case, the method of partitioning based on volume-minimization used by the MultiNest Fortran &
                     &library will be used to find the optimal partitioning of the domain of the objective function. &
                     &The 'MultiNest' method of partitioning is similar to 'ParaMonte-MinVol' in the sense that both methods &
                     &minimize the collective volumes of the bounding objects. However, the minimization approach taken by MultiNest is much &
                     &more aggressive and often results in systematically biased results, especially in high-dimensional problems."//NL2//&
            SKC_"+   domainPartitionMethod = 'Dynesty', all case-INsensitive."//NL2//&
            SKC_"    in which case, the method of partitioning based on volume-minimization used by the Dynesty Python &
                     &library will be used to find the optimal partitioning of the domain of the objective function. &
                     &The 'often' method of partitioning is similar to 'ParaMonte-MinVol' in the sense that both methods &
                     &minimize the collective volumes of the bounding objects. However, the minimization approach taken by often is less &
                     &aggressive and often results in lower efficiencies for the Nested Sampler."//NL2//&
            SKC_"Note that all values are case-INsensitive and all hyphens (dashes, -), white-space, or other separators are ignored.&
                &The default value is '"//spec%domainPartitionMethod%def//SKC_"'. Choose the other non-default methods ONLY if you know the implications or &
                &if you want to perform benchmarks."
            !$omp master
            domainPartitionMethod = spec%domainPartitionMethod%null
            !$omp end master
        end block domainPartitionMethod_block

        domainPartitionObject_block: block
            use pm_sampling_scio, only: domainPartitionObject
            spec%domainPartitionObject%def = SKC_"hyper-ellipsoid"
            spec%domainPartitionObject%null = repeat(SUB, len(domainPartitionObject))
            spec%domainPartitionObject%desc = &
            SKC_"The simulation specification `domainPartitionObject` is a scalar string of default kind of maximum length "//getStr(domainPartitionObject)//SKC_" containing &
                &the name of the ndim-dimensional geometrical object to be used for partitioning of the domain of the objective function. The string value must be enclosed &
                &by either single or double quotation marks when specified from within an external input file. Options that are currently supported include:"//NL2//&
            SKC_"    domainPartitionObject = 'ball'"//NL2//&
            SKC_"            This is equivalent to using hyper-ellipsoids to partition the domain of the objective function."//NL2//&
            SKC_"Note that all values are case-INsensitive and all hyphens (dashes, -) and white-space characters are ignored.&
                &The default value is '"//spec%domainPartitionObject%def//SKC_"'."
            !$omp master
            domainPartitionObject = spec%domainPartitionObject%null
            !$omp end master
        end block domainPartitionObject_block

        domainPartitionOptimizationScaleEnabled_block: block
            use pm_sampling_scio, only: domainPartitionOptimizationScaleEnabled
            spec%domainPartitionOptimizationScaleEnabled%def = .true._LK
            spec%domainPartitionOptimizationScaleEnabled%desc = &
            SKC_"The simulation specification `domainPartitionOptimizationScaleEnabled` is a scalar of type `logical` (Boolean). &
                &When `domainPartitionOptimizationScaleEnabled` is set to the logical true value (equivalent to `.true.`, `.t.` or `true`, all case-INsensitive, &
                &from within an external input file or `True` from inside Python or `true` from inside MATLAB, the final scales of the partitions of the live &
                &points of the sampler will be optimized. This will minimize the inaccuracies in the partitioning of the domain of the objective &
                &function, although, it comes with an extra cost in terms of runtime efficiency of the "//spec%method%val//SKC_" sampler. &
                &The default value is "//getStr(spec%domainPartitionOptimizationScaleEnabled%def)//SKC_"."
            !$omp master
            domainPartitionOptimizationScaleEnabled = spec%domainPartitionOptimizationScaleEnabled%def
            !$omp end master
        end block domainPartitionOptimizationScaleEnabled_block

        domainPartitionOptimizationShapeEnabled_block: block
            use pm_sampling_scio, only: domainPartitionOptimizationShapeEnabled
            spec%domainPartitionOptimizationShapeEnabled%def = .true._LK
            spec%domainPartitionOptimizationShapeEnabled%desc = &
            SKC_"The simulation specification `domainPartitionOptimizationShapeEnabled` is a scalar of type `logical` (Boolean). &
                &When domainPartitionOptimizationShapeEnabled = .true. (or .t. or true, all case-INsensitive, from inside the input file) &
                &or True (from inside Python) or true (from inside MATLAB), the final centers and shapes and sizes of the partitions of the &
                &live points of the sampler will be optimized. This will minimize the inaccuracies in the partitioning of the domain of &
                &the objective function, although, it comes with an extra cost in terms of runtime efficiency of the "//spec%method%val//SKC_" sampler. &
                &The default value is "//getStr(spec%domainPartitionOptimizationShapeEnabled%def)//SKC_"."
            !$omp master
            domainPartitionOptimizationShapeEnabled = spec%domainPartitionOptimizationShapeEnabled%def
            !$omp end master
        end block domainPartitionOptimizationShapeEnabled_block

        domainPartitionOptimizationShapeScaleEnabled_block: block
            use pm_sampling_scio, only: domainPartitionOptimizationShapeScaleEnabled
            spec%domainPartitionOptimizationShapeScaleEnabled%def = .true._LK
            spec%domainPartitionOptimizationShapeScaleEnabled%desc = &
            SKC_"The simulation specification `domainPartitionOptimizationShapeScaleEnabled` = .true. (or .t. or true, all case-INsensitive, from inside the input file) &
                &When `domainPartitionOptimizationShapeScaleEnabled` is set to the logical true value (equivalent to `.true.`, `.t.` or `true`, all case-INsensitive, &
                &from within an external input file or `True` from inside Python or `true` from inside MATLAB, the final centers and shapes and sizes of the partitions &
                &of the live live points of the sampler will be enlarged for as much as possible. This will minimize the inaccuracies in the &
                &partitioning of the domain of the objective function, although, it comes with an extra cost in terms of runtime efficiency &
                &of the "//spec%method%val//SKC_" sampler. The default value is "//getStr(spec%domainPartitionOptimizationShapeScaleEnabled%def)//SKC_"."
            !$omp master
            domainPartitionOptimizationShapeScaleEnabled = spec%domainPartitionOptimizationShapeScaleEnabled%def
            !$omp end master
        end block domainPartitionOptimizationShapeScaleEnabled_block

        domainPartitionScaleFactor_block: block
            use pm_sampling_scio, only: domainPartitionScaleFactor
            spec%domainPartitionScaleFactor%def = 1.5_RKC
            spec%domainPartitionScaleFactor%desc = &
            SKC_"The simulation specification `domainPartitionScaleFactor` is a positive-valued scalar of type `real` of the highest precision available within the &
                &ParaMonte library representing the factor by which the effective radius of the rejection sampling region is enlarged to ensure a minimally-biased integration. &
                &Values smaller than 1 will lead to highly biased simulation results whereas values significantly larger than 1 will lead to highly inefficient simulations. &
                &Relevant `domainPartitionScaleFactor` values are generally bound to the range `1 < domainPartitionScaleFactor < 2`. &
                &The default value for `domainPartitionScaleFactor` is "//getStr(spec%domainPartitionScaleFactor%def)//SKC_"."
            !$omp master
            call setNAN(domainPartitionScaleFactor)
            !$omp end master
        end block domainPartitionScaleFactor_block

        domainSampler_block: block
            use pm_sampling_scio, only: domainSampler
            spec%domainSampler%def = SKC_"rejection" ! spec%%rejection//SKC_"-"//spec%%ellipsoidal
            spec%domainSampler%null = repeat(SUB, len(domainSampler))
            spec%domainSampler%desc = &
            SKC_"The simulation specification `domainSampler` is a scalar string of default kind of maximum length "//getStr(domainSampler)//SKC_" containing the name of &
                &the domain sampler for the Nested sampling integration. The string value must be singly or doubly quoted when specified within an external input file. &
                &Options that are currently supported include:"//NL2//&
            SKC_"    domainSampler = 'rejection'"//NL2//&
            SKC_"            This is equivalent to the constrained rejection sampling via partitions formed by bounding objects around the set of live (active) points."//NL2//&
            SKC_"Note that all values are case-INsensitive and all hyphens (dashes, -) and white-space characters are ignored. The default value for `domainSampler` is '"//spec%domainSampler%def//SKC_"'."
            !$omp master
            domainSampler = spec%domainSampler%null
            !$omp end master
        end block domainSampler_block

        liveSampleSize_block: block
            use pm_sampling_scio, only: liveSampleSize
            spec%liveSampleSize%def  = 1000_IK ! @todo: perhaps an intelligent setting strategy here would be useful.
            spec%liveSampleSize%null = -huge(0_IK)
            spec%liveSampleSize%desc = &
            SKC_"The simulation specification `liveSampleSize` is a positive scalar of type `integer` of default kind representing the number of live (active) &
                &points that are initially sampled uniformly from the domain of the objective function. New points will be subsequently added to this sample &
                &throughout the simulation and the lowest-value will be removed such that the size of the live (active) sample of points at any stage during &
                &the simulation remains fixed as specified by the input simulation specification `liveSampleSize`. &
                &The condition `ndim < liveSampleSize` must hold at all times, where `ndim` is number of dimensions of the domain of the objective function. &
                &The specified value for `liveSampleSize` should be preferably an integer on the order of hundreds or thousands. The larger the value &
                &of `liveSampleSize` is, the better the odds of capturing all modes of the objective function will be. &
                &You can think of `liveSampleSize` (the number of live (active) points) as a random meshing of the integration domain of objective function. &
                &Therefore, a larger `liveSampleSize` means a finer mesh for the integration, thus yielding more accurate results. &
                &In general, higher-dimensional domains require larger the values specified for `liveSampleSize`. &
                &The default value for `liveSampleSize` is "//getStr(spec%liveSampleSize%def)//SKC_"."
            !$omp master
            liveSampleSize = spec%liveSampleSize%null
            !$omp end master
        end block liveSampleSize_block

        tolerance_block: block
            use pm_sampling_scio, only: tolerance
            spec%tolerance%def = 0.01_RKC
           !spec%tolerance%null = NULL_RKC
            spec%tolerance%desc = &
            SKC_"The simulation specification `tolerance` is a positive scalar of type `real` of the highest precision available within the ParaMonte library &
                &representing the threshold below which the integration of the objective function by the "//spec%method%val//" sampler is assumed to have converged. &
                &Typical values are around `0.1 - 1`. Avoid setting this number to extremely tiny values on the order of `1e-6` or smaller, otherwise the integration &
                &result is not guaranteed to converge and the simulation might never end. The default value for `tolerance` is "//getStr(spec%tolerance%def)//SKC_"."
            !$omp master
            call setNAN(tolerance)
            !$omp end master
        end block tolerance_block

    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function set(spec, sampler) result(err)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: set
#endif
        use pm_err, only: err_type
        use pm_sampling, only: paranest_type
        use pm_sampling, only: sampler_type
        class(specnest_type), intent(inout) :: spec
        class(sampler_type), intent(in), optional :: sampler
        type(err_type) :: err

        select type(sampler)
        type is (paranest_type)

            err = spec%specbase_type%set(sampler%sampler_type)
            if (err%occurred) return

            domainPartitionAdaptationCount_block: block
                use pm_sampling_scio, only: domainPartitionAdaptationCount
                if (spec%overridable .and. allocated(sampler%domainPartitionAdaptationCount)) then
                    spec%domainPartitionAdaptationCount%val = sampler%domainPartitionAdaptationCount
                else
                    spec%domainPartitionAdaptationCount%val = domainPartitionAdaptationCount
                end if
                if (spec%domainPartitionAdaptationCount%val == spec%domainPartitionAdaptationCount%null) spec%domainPartitionAdaptationCount%val = spec%domainPartitionAdaptationCount%def
            end block domainPartitionAdaptationCount_block

            domainPartitionAdaptationPeriod_block: block
                use pm_sampling_scio, only: domainPartitionAdaptationPeriod
                if (spec%overridable .and. allocated(sampler%domainPartitionAdaptationPeriod)) then
                    spec%domainPartitionAdaptationPeriod%val = sampler%domainPartitionAdaptationPeriod
                else
                    spec%domainPartitionAdaptationPeriod%val = domainPartitionAdaptationPeriod
                end if
                if (spec%domainPartitionAdaptationPeriod%val == spec%domainPartitionAdaptationPeriod%null) spec%domainPartitionAdaptationPeriod%val = spec%domainPartitionAdaptationPeriod%def
            end block domainPartitionAdaptationPeriod_block

            domainPartitionBiasCorrectionEnabled_block: block
                use pm_sampling_scio, only: domainPartitionBiasCorrectionEnabled
                if (spec%overridable .and. allocated(sampler%domainPartitionBiasCorrectionEnabled)) then
                    spec%domainPartitionBiasCorrectionEnabled%val = sampler%domainPartitionBiasCorrectionEnabled
                else
                    spec%domainPartitionBiasCorrectionEnabled%val = domainPartitionBiasCorrectionEnabled
                end if
            end block domainPartitionBiasCorrectionEnabled_block

            domainPartitionCountMax_block: block
                use pm_sampling_scio, only: domainPartitionCountMax
                if (spec%overridable .and. allocated(sampler%domainPartitionCountMax)) then
                    spec%domainPartitionCountMax%val = sampler%domainPartitionCountMax
                else
                    spec%domainPartitionCountMax%val = domainPartitionCountMax
                end if
                if (spec%domainPartitionCountMax%val == spec%domainPartitionCountMax%null) spec%domainPartitionCountMax%val = spec%domainPartitionCountMax%def
            end block domainPartitionCountMax_block

            domainPartitionFactorExpansion_block: block
                use pm_sampling_scio, only: domainPartitionFactorExpansion
                if (spec%overridable .and. allocated(sampler%domainPartitionFactorExpansion)) then
                    spec%domainPartitionFactorExpansion%val = real(sampler%domainPartitionFactorExpansion, RKC)
                else
                    spec%domainPartitionFactorExpansion%val = domainPartitionFactorExpansion
                end if
                if (isNAN(spec%domainPartitionFactorExpansion%val)) spec%domainPartitionFactorExpansion%val = spec%domainPartitionFactorExpansion%def
            end block domainPartitionFactorExpansion_block

            domainPartitionFactorShrinkage_block: block
                use pm_sampling_scio, only: domainPartitionFactorShrinkage
                if (spec%overridable .and. allocated(sampler%domainPartitionFactorShrinkage)) then
                    spec%domainPartitionFactorShrinkage%val = real(sampler%domainPartitionFactorShrinkage, RKC)
                else
                    spec%domainPartitionFactorShrinkage%val = domainPartitionFactorShrinkage
                end if
                if (isNAN(spec%domainPartitionFactorShrinkage%val)) spec%domainPartitionFactorShrinkage%val = spec%domainPartitionFactorShrinkage%def
            end block domainPartitionFactorShrinkage_block

            domainPartitionKmeansClusterCountMax_block: block
                use pm_sampling_scio, only: domainPartitionKmeansClusterCountMax
                if (spec%overridable .and. allocated(sampler%domainPartitionKmeansClusterCountMax)) then
                    spec%domainPartitionKmeansClusterCountMax%val = sampler%domainPartitionKmeansClusterCountMax
                else
                    spec%domainPartitionKmeansClusterCountMax%val = domainPartitionKmeansClusterCountMax
                end if
                if (spec%domainPartitionKmeansClusterCountMax%val == spec%domainPartitionKmeansClusterCountMax%null) spec%domainPartitionKmeansClusterCountMax%val = spec%domainPartitionKmeansClusterCountMax%def
            end block domainPartitionKmeansClusterCountMax_block

            domainPartitionKmeansClusterSizeMin_block: block
                use pm_sampling_scio, only: domainPartitionKmeansClusterSizeMin
                if (spec%overridable .and. allocated(sampler%domainPartitionKmeansClusterSizeMin)) then
                    spec%domainPartitionKmeansClusterSizeMin%val = sampler%domainPartitionKmeansClusterSizeMin
                else
                    spec%domainPartitionKmeansClusterSizeMin%val = domainPartitionKmeansClusterSizeMin
                end if
                if (spec%domainPartitionKmeansClusterSizeMin%val == spec%domainPartitionKmeansClusterSizeMin%null) spec%domainPartitionKmeansClusterSizeMin%val = spec%domainPartitionKmeansClusterSizeMin%def
            end block domainPartitionKmeansClusterSizeMin_block

            domainPartitionKmeansNormalizationEnabled_block: block
                use pm_sampling_scio, only: domainPartitionKmeansNormalizationEnabled
                if (spec%overridable .and. allocated(sampler%domainPartitionKmeansNormalizationEnabled)) then
                    spec%domainPartitionKmeansNormalizationEnabled%val = sampler%domainPartitionKmeansNormalizationEnabled
                else
                    spec%domainPartitionKmeansNormalizationEnabled%val = domainPartitionKmeansNormalizationEnabled
                end if
            end block domainPartitionKmeansNormalizationEnabled_block

            domainPartitionKmeansNumFailMax_block: block
                use pm_sampling_scio, only: domainPartitionKmeansNumFailMax
                if (spec%overridable .and. allocated(sampler%domainPartitionKmeansNumFailMax)) then
                    spec%domainPartitionKmeansNumFailMax%val = sampler%domainPartitionKmeansNumFailMax
                else
                    spec%domainPartitionKmeansNumFailMax%val = domainPartitionKmeansNumFailMax
                end if
                if (spec%domainPartitionKmeansNumFailMax%val == spec%domainPartitionKmeansNumFailMax%null) spec%domainPartitionKmeansNumFailMax%val = spec%domainPartitionKmeansNumFailMax%def
            end block domainPartitionKmeansNumFailMax_block

            domainPartitionKmeansNumRecursionMax_block: block
                use pm_sampling_scio, only: domainPartitionKmeansNumRecursionMax
                if (spec%overridable .and. allocated(sampler%domainPartitionKmeansNumRecursionMax)) then
                    spec%domainPartitionKmeansNumRecursionMax%val = sampler%domainPartitionKmeansNumRecursionMax
                else
                    spec%domainPartitionKmeansNumRecursionMax%val = domainPartitionKmeansNumRecursionMax
                end if
                if (spec%domainPartitionKmeansNumRecursionMax%val == spec%domainPartitionKmeansNumRecursionMax%null) spec%domainPartitionKmeansNumRecursionMax%val = spec%domainPartitionKmeansNumRecursionMax%def
            end block domainPartitionKmeansNumRecursionMax_block

            domainPartitionKmeansNumTry_block: block
                use pm_sampling_scio, only: domainPartitionKmeansNumTry
                if (spec%overridable .and. allocated(sampler%domainPartitionKmeansNumTry)) then
                    spec%domainPartitionKmeansNumTry%val = sampler%domainPartitionKmeansNumTry
                else
                    spec%domainPartitionKmeansNumTry%val = domainPartitionKmeansNumTry
                end if
                if (spec%domainPartitionKmeansNumTry%val == spec%domainPartitionKmeansNumTry%null) spec%domainPartitionKmeansNumTry%val = spec%domainPartitionKmeansNumTry%def
            end block domainPartitionKmeansNumTry_block

            domainPartitionKvolumeNumRecursionMax_block: block
                use pm_sampling_scio, only: domainPartitionKvolumeNumRecursionMax
                if (spec%overridable .and. allocated(sampler%domainPartitionKvolumeNumRecursionMax)) then
                    spec%domainPartitionKvolumeNumRecursionMax%val = sampler%domainPartitionKvolumeNumRecursionMax
                else
                    spec%domainPartitionKvolumeNumRecursionMax%val = domainPartitionKvolumeNumRecursionMax
                end if
                if (spec%domainPartitionKvolumeNumRecursionMax%val == spec%domainPartitionKvolumeNumRecursionMax%null) spec%domainPartitionKvolumeNumRecursionMax%val = spec%domainPartitionKvolumeNumRecursionMax%def
            end block domainPartitionKvolumeNumRecursionMax_block

            domainPartitionKvolumeWeightExponent_block: block
                use pm_sampling_scio, only: domainPartitionKvolumeWeightExponent
                if (spec%overridable .and. allocated(sampler%domainPartitionKvolumeWeightExponent)) then
                    spec%domainPartitionKvolumeWeightExponent%val = real(sampler%domainPartitionKvolumeWeightExponent, RKC)
                else
                    spec%domainPartitionKvolumeWeightExponent%val = domainPartitionKvolumeWeightExponent
                end if
                if (isNAN(spec%domainPartitionKvolumeWeightExponent%val)) spec%domainPartitionKvolumeWeightExponent%val = spec%domainPartitionKvolumeWeightExponent%def
            end block domainPartitionKvolumeWeightExponent_block

            domainPartitionMethod_block: block
                use pm_sampling_scio, only: domainPartitionMethod
                if (spec%overridable .and. allocated(sampler%domainPartitionMethod)) then
                    spec%domainPartitionMethod%val = sampler%domainPartitionMethod
                else
                    spec%domainPartitionMethod%val = getStrLower(trim(adjustl(domainPartitionMethod)))
                end if
                if (spec%domainPartitionMethod%val == spec%domainPartitionMethod%null) spec%domainPartitionMethod%val = getStrLower(spec%domainPartitionMethod%def)
                spec%domainPartitionMethod%isMinVol = index(spec%domainPartitionMethod%val, SKC_"minvol") > 0
                spec%domainPartitionMethod%isMaxDen = index(spec%domainPartitionMethod%val, SKC_"maxden") > 0
                spec%domainPartitionMethod%isDynesty = index(spec%domainPartitionMethod%val, SKC_"dynesty") > 0
                spec%domainPartitionMethod%isMultiNest = index(spec%domainPartitionMethod%val, SKC_"multinest") > 0
            end block domainPartitionMethod_block

            domainPartitionObject_block: block
                use pm_sampling_scio, only: domainPartitionObject
                if (spec%overridable .and. allocated(sampler%domainPartitionObject)) then
                    spec%domainPartitionObject%val = sampler%domainPartitionObject
                else
                    spec%domainPartitionObject%val = getStrLower(trim(adjustl(domainPartitionObject)))
                end if
                if (spec%domainPartitionObject%val == spec%domainPartitionObject%null) spec%domainPartitionObject%val = spec%domainPartitionObject%def
                spec%domainPartitionObject%isBall = index(spec%domainPartitionObject%val, SKC_"ellipsoid") > 0 .or. index(spec%domainPartitionObject%val, SKC_"ball") > 0
            end block domainPartitionObject_block

            domainPartitionOptimizationScaleEnabled_block: block
                use pm_sampling_scio, only: domainPartitionOptimizationScaleEnabled
                if (spec%overridable .and. allocated(sampler%domainPartitionOptimizationScaleEnabled)) then
                    spec%domainPartitionOptimizationScaleEnabled%val = sampler%domainPartitionOptimizationScaleEnabled
                else
                    spec%domainPartitionOptimizationScaleEnabled%val = domainPartitionOptimizationScaleEnabled
                end if
            end block domainPartitionOptimizationScaleEnabled_block

            domainPartitionOptimizationShapeEnabled_block: block
                use pm_sampling_scio, only: domainPartitionOptimizationShapeEnabled
                if (spec%overridable .and. allocated(sampler%domainPartitionOptimizationShapeEnabled)) then
                    spec%domainPartitionOptimizationShapeEnabled%val = sampler%domainPartitionOptimizationShapeEnabled
                else
                    spec%domainPartitionOptimizationShapeEnabled%val = domainPartitionOptimizationShapeEnabled
                end if
            end block domainPartitionOptimizationShapeEnabled_block

            domainPartitionOptimizationShapeScaleEnabled_block: block
                use pm_sampling_scio, only: domainPartitionOptimizationShapeScaleEnabled
                if (spec%overridable .and. allocated(sampler%domainPartitionOptimizationShapeScaleEnabled)) then
                    spec%domainPartitionOptimizationShapeScaleEnabled%val = sampler%domainPartitionOptimizationShapeScaleEnabled
                else
                    spec%domainPartitionOptimizationShapeScaleEnabled%val = domainPartitionOptimizationShapeScaleEnabled
                end if
            end block domainPartitionOptimizationShapeScaleEnabled_block

            domainPartitionScaleFactor_block: block
                use pm_sampling_scio, only: domainPartitionScaleFactor
                if (spec%overridable .and. allocated(sampler%domainPartitionScaleFactor)) then
                    spec%domainPartitionScaleFactor%val = real(sampler%domainPartitionScaleFactor, RKC)
                else
                    spec%domainPartitionScaleFactor%val = domainPartitionScaleFactor
                end if
                if (isNAN(spec%domainPartitionScaleFactor%val)) spec%domainPartitionScaleFactor%val = spec%domainPartitionScaleFactor%def
            end block domainPartitionScaleFactor_block

            domainSampler_block: block
                use pm_sampling_scio, only: domainSampler
                if (spec%overridable .and. allocated(sampler%domainSampler)) then
                    spec%domainSampler%val = sampler%domainSampler
                else
                    spec%domainSampler%val = getStrLower(trim(adjustl(domainSampler)))
                end if
                if (spec%domainSampler%val == spec%domainSampler%null) spec%domainSampler%val = spec%domainSampler%def
                spec%domainSampler%isRejection = index(spec%domainSampler%val, SKC_"rejection") > 0
            end block domainSampler_block

            liveSampleSize_block: block
                use pm_sampling_scio, only: liveSampleSize
                if (spec%overridable .and. allocated(sampler%liveSampleSize)) then
                    spec%liveSampleSize%val = sampler%liveSampleSize
                else
                    spec%liveSampleSize%val = liveSampleSize
                end if
                if (spec%liveSampleSize%val == spec%liveSampleSize%null) spec%liveSampleSize%val = spec%liveSampleSize%def
            end block liveSampleSize_block

            tolerance_block: block
                use pm_sampling_scio, only: tolerance
                if (spec%overridable .and. allocated(sampler%tolerance)) then
                    spec%tolerance%val = real(sampler%tolerance, RKC)
                else
                    spec%tolerance%val = tolerance
                end if
                if (isNAN(spec%tolerance%val)) spec%tolerance%val = spec%tolerance%def
            end block tolerance_block

            !$omp barrier

            ! Now resolve the conflicts.

            spec%domainPartitionKmeansClusterCountMax%val = min(spec%domainPartitionKmeansClusterCountMax%def, spec%domainPartitionCountMax%val)

            ! open output files, report and sanitize.

            if (spec%image%is%leader) call spec%report() ! if (spec%run%is%new)
            call spec%sanitize(err)

        class default
            error stop "The input `sampler` must be of type `paranest_type`." ! LCOV_EXCL_LINE
        end select

    end function set

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine report(spec)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: report
#endif
        class(specnest_type), intent(inout) :: spec

        call spec%disp%text%wrap(NL1//spec%method%val//SKC_" simulation Nest specifications"//NL1)

        associate(ndim => spec%ndim%val, format => spec%reportFile%format%generic)

            call spec%disp%show("domainPartitionAdaptationCount")
            call spec%disp%show(spec%domainPartitionAdaptationCount%val, format = format)
            call spec%disp%note%show(spec%domainPartitionAdaptationCount%desc)

            call spec%disp%show("domainPartitionAdaptationPeriod")
            call spec%disp%show(spec%domainPartitionAdaptationPeriod%val, format = format)
            call spec%disp%note%show(spec%domainPartitionAdaptationPeriod%desc)

            call spec%disp%show("domainPartitionBiasCorrectionEnabled")
            call spec%disp%show(spec%domainPartitionBiasCorrectionEnabled%val, format = format)
            call spec%disp%note%show(spec%domainPartitionBiasCorrectionEnabled%desc)

            call spec%disp%show("domainPartitionCountMax")
            call spec%disp%show(spec%domainPartitionCountMax%val, format = format)
            call spec%disp%note%show(spec%domainPartitionCountMax%desc)

            call spec%disp%show("domainPartitionFactorExpansion")
            call spec%disp%show(spec%domainPartitionFactorExpansion%val, format = format)
            call spec%disp%note%show(spec%domainPartitionFactorExpansion%desc)

            call spec%disp%show("domainPartitionFactorShrinkage")
            call spec%disp%show(spec%domainPartitionFactorShrinkage%val, format = format)
            call spec%disp%note%show(spec%domainPartitionFactorShrinkage%desc)

            call spec%disp%show("domainPartitionKmeansClusterCountMax")
            call spec%disp%show(spec%domainPartitionKmeansClusterCountMax%val, format = format)
            call spec%disp%note%show(spec%domainPartitionKmeansClusterCountMax%desc)

            call spec%disp%show("domainPartitionKmeansClusterSizeMin")
            call spec%disp%show(spec%domainPartitionKmeansClusterSizeMin%val, format = format)
            call spec%disp%note%show(spec%domainPartitionKmeansClusterSizeMin%desc)

            call spec%disp%show("domainPartitionKmeansNormalizationEnabled")
            call spec%disp%show(spec%domainPartitionKmeansNormalizationEnabled%val, format = format)
            call spec%disp%note%show(spec%domainPartitionKmeansNormalizationEnabled%desc)

            call spec%disp%show("domainPartitionKmeansNumFailMax")
            call spec%disp%show(spec%domainPartitionKmeansNumFailMax%val, format = format)
            call spec%disp%note%show(spec%domainPartitionKmeansNumFailMax%desc)

            call spec%disp%show("domainPartitionKmeansNumRecursionMax")
            call spec%disp%show(spec%domainPartitionKmeansNumRecursionMax%val, format = format)
            call spec%disp%note%show(spec%domainPartitionKmeansNumRecursionMax%desc)

            call spec%disp%show("domainPartitionKmeansNumTry")
            call spec%disp%show(spec%domainPartitionKmeansNumTry%val, format = format)
            call spec%disp%note%show(spec%domainPartitionKmeansNumTry%desc)

            call spec%disp%show("domainPartitionKvolumeNumRecursionMax")
            call spec%disp%show(spec%domainPartitionKvolumeNumRecursionMax%val, format = format)
            call spec%disp%note%show(spec%domainPartitionKvolumeNumRecursionMax%desc)

            call spec%disp%show("domainPartitionKvolumeWeightExponent")
            call spec%disp%show(spec%domainPartitionKvolumeWeightExponent%val, format = format)
            call spec%disp%note%show(spec%domainPartitionKvolumeWeightExponent%desc)

            call spec%disp%show("domainPartitionMethod")
            call spec%disp%show(spec%domainPartitionMethod%val, format = format)
            call spec%disp%note%show(spec%domainPartitionMethod%desc)

            call spec%disp%show("domainPartitionObject")
            call spec%disp%show(spec%domainPartitionObject%val, format = format)
            call spec%disp%note%show(spec%domainPartitionObject%desc)

            call spec%disp%show("domainPartitionOptimizationScaleEnabled")
            call spec%disp%show(spec%domainPartitionOptimizationScaleEnabled%val, format = format)
            call spec%disp%note%show(spec%domainPartitionOptimizationScaleEnabled%desc)

            call spec%disp%show("domainPartitionOptimizationShapeEnabled")
            call spec%disp%show(spec%domainPartitionOptimizationShapeEnabled%val, format = format)
            call spec%disp%note%show(spec%domainPartitionOptimizationShapeEnabled%desc)

            call spec%disp%show("domainPartitionOptimizationShapeScaleEnabled")
            call spec%disp%show(spec%domainPartitionOptimizationShapeScaleEnabled%val, format = format)
            call spec%disp%note%show(spec%domainPartitionOptimizationShapeScaleEnabled%desc)

            call spec%disp%show("domainPartitionScaleFactor")
            call spec%disp%show(spec%domainPartitionScaleFactor%val, format = format)
            call spec%disp%note%show(spec%domainPartitionScaleFactor%desc)

            call spec%disp%show("domainSampler")
            call spec%disp%show(spec%domainSampler%val, format = format)
            call spec%disp%note%show(spec%domainSampler%desc)

            call spec%disp%show("liveSampleSize")
            call spec%disp%show(spec%liveSampleSize%val, format = format)
            call spec%disp%note%show(spec%liveSampleSize%desc)

            call spec%disp%show("tolerance")
            call spec%disp%show(spec%tolerance%val, format = format)
            call spec%disp%note%show(spec%tolerance%desc)

        end associate

    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine sanitize(spec, err)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: sanitize
#endif
        use pm_err, only: err_type, getFine
        type(err_type), intent(inout) :: err
        class(specnest_type), intent(inout) :: spec
        character(*,SKC), parameter :: PROCEDURE_NAME = MODULE_NAME//SKC_"@sanitize()"

        call spec%disp%text%wrap(NL1//spec%method%val//SKC_".simulation.specifications.nest"//NL1)

        domainPartitionAdaptationCount_block: block
            if (spec%domainPartitionAdaptationCount%val < 0_IK) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The specified value for `domainPartitionAdaptationCount` ("//getStr(spec%domainPartitionAdaptationCount%val)//SKC_") cannot be negative. &
                            &If you are unsure of the appropriate value for `domainPartitionAdaptationCount`, drop it from the input list. The "//&
                            spec%method%val//SKC_" sampler will automatically assign an appropriate value to it."
            end if
        end block domainPartitionAdaptationCount_block

        domainPartitionAdaptationPeriod_block: block
            if (spec%domainPartitionAdaptationPeriod%val<1_IK) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &Invalid requested value for `domainPartitionAdaptationPeriod`. &
                            &The specified value for `domainPartitionAdaptationPeriod` ("//getStr(spec%domainPartitionAdaptationPeriod%val)//SKC_") cannot be less than 1. &
                            &If you are unsure of the appropriate value for `domainPartitionAdaptationPeriod`, drop it from the input list. The "//&
                            spec%method%val//SKC_" sample will automatically assign an appropriate value to it."
            end if
        end block domainPartitionAdaptationPeriod_block

        domainPartitionBiasCorrectionEnabled_block: block
        end block domainPartitionBiasCorrectionEnabled_block

        domainPartitionCountMax_block: block
            if (spec%domainPartitionCountMax%val < 1_IK) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The specified value for `domainPartitionCountMax` ("//getStr(spec%domainPartitionCountMax%val)//SKC_") &
                            &cannot be less than `1` or larger than `liveSampleSize / (ndim + 1)`, where `ndim` is the number of dimensions of &
                            &the domain of the objective function. If you are unsure of the appropriate value for `domainPartitionCountMax`, drop it &
                            &from the input list. The "//spec%method%val//SKC_" sampler will automatically assign an appropriate value to it."
            end if
        end block domainPartitionCountMax_block

        domainPartitionFactorExpansion_block: block
            if (spec%domainPartitionFactorExpansion%val <= 0._RKC) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The input variable `domainPartitionFactorExpansion` ("//getStr(spec%domainPartitionFactorExpansion%val)//") must be positive. &
                            &If you are unsure of the appropriate value for `domainPartitionFactorExpansion`, drop it from the input list. &
                            &The "//spec%method%val//SKC_" sampler will automatically assign an appropriate value to it."
            end if
        end block domainPartitionFactorExpansion_block

        domainPartitionFactorShrinkage_block: block
            if (spec%domainPartitionFactorShrinkage%val <= 0._RKC) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The input variable `domainPartitionFactorShrinkage` ("//getStr(spec%domainPartitionFactorShrinkage%val)//SKC_") must be positive. &
                            &If you are unsure of the appropriate value for `domainPartitionFactorExpansion`, drop it from the input list. &
                            &The "//spec%method%val//SKC_" sampler will automatically assign an appropriate value to it."
            end if
        end block domainPartitionFactorShrinkage_block

        domainPartitionKmeansClusterCountMax_block: block
            if (spec%domainPartitionKmeansClusterCountMax%val < 1_IK .or. spec%domainPartitionCountMax%val < spec%domainPartitionKmeansClusterCountMax%val) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The specified value for `domainPartitionKmeansClusterCountMax` ("//getStr(spec%domainPartitionKmeansClusterCountMax%val)//SKC_") &
                            &can not be less than `1` or larger than the value of `domainPartitionCountMax` ("//getStr(spec%domainPartitionCountMax%val)//SKC_"). &
                            &If you are unsure about the appropriate value for domainPartitionKmeansClusterCountMax, drop it from the input list. &
                            &The "//spec%method%val//SKC_" sampler will automatically assign an appropriate value to it."
            end if
        end block domainPartitionKmeansClusterCountMax_block

        domainPartitionKmeansClusterSizeMin_block: block
            if (spec%domainPartitionKmeansClusterSizeMin%val < 0_IK) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The specified value for `domainPartitionKmeansClusterSizeMin` ("//getStr(spec%domainPartitionKmeansClusterSizeMin%val)//SKC_") &
                            &cannot be less than 0 or larger than liveSampleSize. If you are unsure of the appropriate value for &
                            &`domainPartitionKmeansClusterSizeMin`, drop it from the input list. The "//spec%method%val//SKC_" &
                            &sampler will automatically assign an appropriate value to it."
            end if
        end block domainPartitionKmeansClusterSizeMin_block

        domainPartitionKmeansNormalizationEnabled_block: block
        end block domainPartitionKmeansNormalizationEnabled_block

        domainPartitionKmeansNumFailMax_block: block
            if (spec%domainPartitionKmeansNumFailMax%val < 0_IK) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The specified value for `domainPartitionKmeansNumFailMax` ("//getStr(spec%domainPartitionKmeansNumFailMax%val)//SKC_") &
                            &cannot be negative. If you are unsure of the appropriate value for `domainPartitionKmeansNumFailMax`, drop it &
                            &from the input list. "//spec%method%val//SKC_" will automatically assign an appropriate value to it."
            end if
        end block domainPartitionKmeansNumFailMax_block

        domainPartitionKmeansNumRecursionMax_block: block
            if (spec%domainPartitionKmeansNumRecursionMax%val < 0_IK) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The specified value for `domainPartitionKmeansNumRecursionMax` ("//getStr(spec%domainPartitionKmeansNumRecursionMax%val)//SKC_") &
                            &cannot be negative. If you are unsure of the appropriate value for `domainPartitionKmeansNumRecursionMax`, drop it &
                            &from the input list. "//spec%method%val//SKC_" will automatically assign an appropriate value to it."
            end if
        end block domainPartitionKmeansNumRecursionMax_block

        domainPartitionKmeansNumTry_block: block
            if (spec%domainPartitionKmeansNumTry%val < 1_IK) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The specified value for `domainPartitionKmeansNumTry` ("//getStr(spec%domainPartitionKmeansNumTry%val)//SKC_") &
                            &can not be less than 1. If you are unsure about the appropriate value for `domainPartitionKmeansNumTry`, drop it &
                            &from the input list. The "//spec%method%val//SKC_" sampler will automatically assign an appropriate value to it."
            end if
        end block domainPartitionKmeansNumTry_block

        domainPartitionKvolumeNumRecursionMax_block: block
            if (spec%domainPartitionKvolumeNumRecursionMax%val < 0_IK) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The specified value for domainPartitionKvolumeNumRecursionMax ("//getStr(spec%domainPartitionKvolumeNumRecursionMax%val)//SKC_") &
                            &cannot be negative. If you are unsure of the appropriate value for domainPartitionKvolumeNumRecursionMax, drop it &
                            &from the input list. "//spec%method%val//SKC_" will automatically assign an appropriate value to it."
            end if
        end block domainPartitionKvolumeNumRecursionMax_block

        domainPartitionKvolumeWeightExponent_block: block
            !if (spec%domainPartitionKvolumeWeightExponent%val <= 0._RKC) then
            !    err%occurred = .true._LK
            !    err%msg =   err%msg//NL2//&
            !                PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
            !                &The input variable `domainPartitionKvolumeWeightExponent` ("//getStr(spec%domainPartitionKvolumeWeightExponent%val)//&
            !                ") cannot be less than 0. If you are unsure of the appropriate value for `domainPartitionKvolumeWeightExponent`, drop it &
            !                &from the input list. "//spec%method%val//SKC_" will automatically assign an appropriate value to it."
            !end if
        end block domainPartitionKvolumeWeightExponent_block

        domainPartitionMethod_block: block
            if (count(  [ spec%domainPartitionMethod%isMinVol & ! LCOV_EXCL_LINE
                        , spec%domainPartitionMethod%isMaxDen & ! LCOV_EXCL_LINE
                        , spec%domainPartitionMethod%isDynesty & ! LCOV_EXCL_LINE
                        , spec%domainPartitionMethod%isMultiNest & ! LCOV_EXCL_LINE
                        ]) /= 1_IK) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The specified value for `domainPartitionMethod` ("//spec%domainPartitionMethod%val//SKC_") is not supported. &
                            &The variable `domainPartitionMethod` cannot be set to anything other than the possible values described &
                            &in the description of the simulation specification `domainPartitionMethod`. &
                            &Please specify only one partitioning method."
            end if
        end block domainPartitionMethod_block

        domainPartitionObject_block: block
            if  (.not. spec%domainPartitionObject%isBall) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The specified value for the `domainPartitionObject` ("//spec%domainPartitionObject%val//SKC_") is not supported. &
                            &The variable `domainPartitionObject` cannot be set to anything other than the possible values described &
                            &in the description of the simulation specification `domainPartitionObject`. &
                            &Please specify only one partitioning object."
            end if
        end block domainPartitionObject_block

        domainPartitionOptimizationScaleEnabled_block: block
        end block domainPartitionOptimizationScaleEnabled_block

        domainPartitionOptimizationShapeEnabled_block: block
        end block domainPartitionOptimizationShapeEnabled_block

        domainPartitionOptimizationShapeScaleEnabled_block: block
        end block domainPartitionOptimizationShapeScaleEnabled_block

        domainPartitionScaleFactor_block: block
            if (spec%domainPartitionScaleFactor%val <= 0._RKC) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The specified value for `domainPartitionScaleFactor` ("//getStr(spec%domainPartitionScaleFactor%val)//SKC_") cannot be less than 0. &
                            &If you are unsure of the appropriate value for `domainPartitionScaleFactor`, drop it from the input list. &
                            &The "//spec%method%val//SKC_" sampler will automatically assign an appropriate value to it."
            end if
        end block domainPartitionScaleFactor_block

        domainSampler_block: block
            if (.not. spec%domainSampler%isRejection) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The specified value for `domainSampler` ("//spec%domainSampler%val//SKC_") is unsupported. &
                            &The specification `domainSampler` cannot be set to anything other than the values described &
                            &in the description of the simulation specification `domainSampler`."
            end if
        end block domainSampler_block

        liveSampleSize_block: block
            if (spec%liveSampleSize%val < spec%ndim%val + 1_IK) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The specified value for liveSampleSize ("//getStr(spec%liveSampleSize%val)//SKC_") &
                            &cannot be less than or equal to the number of dimensions of the domain of the objective function. &
                            &If you are unsure of the appropriate value for liveSampleSize, drop it from the input list&. &
                            &The "//spec%method%val//SKC_" sampler will automatically assign an appropriate value to it."
            else
                spec%liveSampleSize%log = log(real(spec%liveSampleSize%val, RKC))
            end if
        end block liveSampleSize_block

        tolerance_block: block
            if (spec%tolerance%val <= 0._RKC) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The specified value for `tolerance` ("//getStr(spec%tolerance%val)//SKC_") must be positive. &
                            &If you are unsure of the appropriate value for tolerance, drop it from the input list. &
                            &The "//spec%method%val//SKC_" will automatically assign an appropriate value to it."
            end if
        end block tolerance_block

    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
