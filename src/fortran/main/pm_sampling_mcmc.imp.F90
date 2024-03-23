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

#if CHECK_ENABLED
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif
    use pm_kind, only: SKC => SK, SK, IK, LK
    use pm_matrixClass, only: isMatClass, posdefmat
    use pm_sampling_base, only: specbase_type, astatbase_type, sfcbase_type, NL2, NL1
    use pm_arrayRemove, only: getRemoved
    use pm_arrayResize, only: setResized
    use pm_strASCII, only: getStrLower
    use pm_val2str, only: getStr
    use pm_except, only: setNAN
    use pm_except, only: isNAN
    use pm_strASCII, only: SUB
    use pm_err, only: err_type
    use pm_err, only: getFine

    implicit none

    character(*,SKC)        , parameter     :: MODULE_NAME = SK_"@pm_sampling_mcmc"

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! simulation declarations.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type :: burninLoc_type
        integer(IK) :: compact
        integer(IK) :: verbose
    end type

    ! type for chain refinements
    type, extends(sfcbase_type)             :: sfcmcmc_type
        real(RKC)           , allocatable   :: act(:,:) !<  \public The array of shape `(1:ndim, 0:nref)` containing he Integrated AutoCorrelation Time at each refinement stage.
    contains
        procedure           , pass          :: getErrRefinement
    end type

    type, abstract, extends(astatbase_type) :: astatmcmc_type
        type(burninLoc_type)                :: burninLocMCMC
        type(sfcmcmc_type)                  :: sfc
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! specification declarations.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   !integer(IK)                             :: outputChainSize ! namelist input
    type                                    :: outputChainSize_type
        integer(IK)                         :: val
        integer(IK)                         :: def
        integer(IK)                         :: null
        character(:,SKC)    , allocatable   :: desc
    end type

   !integer(IK)                             :: outputSampleRefinementCount ! namelist input
    type                                    :: outputSampleRefinementCount_type
        character(:,SKC)    , allocatable   :: str
        integer(IK)                         :: val
        integer(IK)                         :: def
        integer(IK)                         :: null
        character(:,SKC)    , allocatable   :: desc
    end type

   !character(63,SKC)                       :: outputSampleRefinementMethod ! namelist input
    type                                    :: outputSampleRefinementMethod_type
        logical(LK)                         :: isBatchMeansCompact = .false._LK
        logical(LK)                         :: isBatchMeansVerbose = .false._LK
        logical(LK)                         :: isMaxCumSumACF = .false._LK
        logical(LK)                         :: isBatchMeans = .false._LK
        character(12,SKC)                   :: maxCumSumACF = "maxCumSumACF"
        character(10,SKC)                   :: batchMeans = "BatchMeans"
        character(9,SKC)                    :: cutoffACF = "cutoffACF"
        character(:,SKC)    , allocatable   :: def
        character(:,SKC)    , allocatable   :: val
        character(:,SKC)    , allocatable   :: null
        character(:,SKC)    , allocatable   :: desc
    end type

   !character(63,SKC)                       :: proposal ! namelist input
    type                                    :: proposalis_type
        logical(LK)                         :: normal = .false._LK
        logical(LK)                         :: uniform = .false._LK
    end type
    type                                    :: proposal_type
        type(proposalis_type)               :: is
        character(7,SKC)                    :: uniform
        character(6,SKC)                    :: normal
        character(:,SKC)    , allocatable   :: val
        character(:,SKC)    , allocatable   :: def
        character(:,SKC)    , allocatable   :: null
        character(:,SKC)    , allocatable   :: desc
    end type

   !real(RKC)               , allocatable   :: proposalCorMat(:,:) ! namelist input
    type                                    :: proposalCorMat_type
        real(RKC)           , allocatable   :: val(:,:)
        real(RKC)           , allocatable   :: def(:,:)
        character(:,SKC)    , allocatable   :: desc
    end type

   !real(RKC)               , allocatable   :: proposalCovMat(:,:) ! namelist input
    type                                    :: proposalCovMat_type
        logical(LK)                         :: isUserSet
        real(RKC)           , allocatable   :: def(:,:)
        real(RKC)           , allocatable   :: val(:,:)
        character(:,SKC)    , allocatable   :: desc
    end type

   !character(127,SKC)                      :: proposalScaleFactor
    type                                    :: proposalScaleFactor_type
        real(RKC)                           :: val, valdef
        character(:,SKC)    , allocatable   :: str, strdef, null, desc
    end type

   !real(RKC)               , allocatable   :: proposalStart(:) ! namelist input
    type                                    :: proposalStart_type
        real(RKC)           , allocatable   :: val(:)
        character(:,SKC)    , allocatable   :: desc
    end type

   !real(RKC)               , allocatable   :: proposalStartDomainCubeLimitLower(:) ! namelist input
    type                                    :: proposalStartDomainCubeLimitLower_type
        real(RKC)                           :: null
        real(RKC)                           :: def
        real(RKC)           , allocatable   :: val(:)
        character(:,SKC)    , allocatable   :: desc
    end type

   !real(RKC)               , allocatable   :: proposalStartDomainCubeLimitUpper(:) ! namelist input
    type                                    :: proposalStartDomainCubeLimitUpper_type
        real(RKC)                           :: def
        real(RKC)           , allocatable   :: val(:)
        character(:,SKC)    , allocatable   :: desc
    end type

   !logical(LK)                             :: proposalStartRandomized ! namelist input
    type                                    :: proposalStartRandomized_type
        logical(LK)                         :: val
        logical(LK)                         :: def
        character(:,SKC)    , allocatable   :: desc
    end type

   !real(RKC)               , allocatable   :: proposalStdVec(:) ! namelist input
    type                                    :: proposalStdVec_type
        real(RKC)           , allocatable   :: val(:)
        real(RKC)           , allocatable   :: def(:)
        character(:,SKC)    , allocatable   :: desc
    end type

    type, extends(specbase_type)                        :: specmcmc_type
        type(outputChainSize_type)                      :: outputChainSize
        type(proposal_type)                             :: proposal
        type(ProposalScaleFactor_type)                  :: ProposalScaleFactor
        type(proposalStart_type)                        :: proposalStart
        type(proposalStdVec_type)                       :: proposalStdVec
        type(proposalCorMat_type)                       :: proposalCorMat
        type(proposalCovMat_type)                       :: proposalCovMat
        type(outputSampleRefinementCount_type)          :: outputSampleRefinementCount
        type(outputSampleRefinementMethod_type)         :: outputSampleRefinementMethod
        type(proposalStartRandomized_type)              :: proposalStartRandomized
        type(proposalStartDomainCubeLimitLower_type)    :: proposalStartDomainCubeLimitLower
        type(proposalStartDomainCubeLimitUpper_type)    :: proposalStartDomainCubeLimitUpper
    contains
        procedure, pass, private                        :: sanitize
        procedure, pass, private                        :: report
        procedure, pass, public                         :: set
    end type

    interface specmcmc_type
        module procedure                                :: construct
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type, private   :: method_type
        logical(LK) :: isMax = .false._LK
        logical(LK) :: isMed = .false._LK
        logical(LK) :: isMin = .false._LK
        logical(LK) :: isAvg = .false._LK
        logical(LK) :: isCumSum = .false._LK
        logical(LK) :: isCumSumMax = .false._LK
        logical(LK) :: isBatchMeans = .false._LK
        logical(LK) :: isBatchMeansMax = .false._LK
        logical(LK) :: isViaCompactChain = .true._LK
        logical(LK) :: isViaVerboseChain = .true._LK
    end type

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
        type(specmcmc_type) :: spec

        spec%specbase_type = specbase_type(modelr, method, ndim)

        outputChainSize_block: block
            use pm_sampling_scio, only: outputChainSize
            spec%outputChainSize%null   = -huge(0_IK)
            spec%outputChainSize%def    = 100000_IK
            spec%outputChainSize%desc   = &
            SKC_"The simulation specification `outputChainSize` is a positive scalar of type `integer` whose value determines the number &
                &of non-refined, potentially auto-correlated, but unique samples drawn by the MCMC sampler before stopping the sampling. &
                &For example, if `outputChainSize = 10000`, then `10000` unique sample points (with no duplicates) will be drawn from the &
                &target objective function that the user has provided. &
                &The input value for `outputChainSize` must be a positive integer of a minimum value `ndim + 1` or larger, &
                &where `ndim` is the number of dimensions of the domain of the objective function to be sampled. &
                &Note that `outputChainSize` is different from and always smaller than the length of the constructed MCMC chain. &
                &The default value is `outputChainSize = "//getStr(spec%outputChainSize%def)//SKC_"`."
            !!$omp master
            outputChainSize = spec%outputChainSize%null
            !!$omp end master
        end block outputChainSize_block

        outputSampleRefinementCount_block: block
            use pm_sampling_scio, only: outputSampleRefinementCount
            spec%outputSampleRefinementCount%null = -huge(0_IK)
            spec%outputSampleRefinementCount%def  = huge(0_IK)
            spec%outputSampleRefinementCount%desc = &
            SKC_"The simulation specification `outputSampleRefinementCount` is a positive-valued scalar of type `integer`. &
                &When `outputSampleSize < 0`, the value of `outputSampleRefinementCount` dictates the maximum number of times &
                &the MCMC chain will be refined to remove the autocorrelation within the output MCMC sample. For example,"//NL2//&
            SKC_"+   if `outputSampleRefinementCount = 0`,"//NL2//&
            SKC_"    no refinement of the output MCMC chain will be performed. The resulting MCMC sample will correspond &
                     &to the full MCMC chain in verbose format (i.e., each sampled state has a weight of one)."//NL2//&
            SKC_"+   if `outputSampleRefinementCount = 1`,"//NL2//&
            SKC_"    the refinement of the output MCMC chain will be done only once if needed, and no more, &
                     &even though some residual autocorrelation in the output MCMC sample may still exist. &
                     &In practice, only one refinement of the final output MCMC chain should be enough to remove &
                     &the existing autocorrelations in the final output sample. Exceptions occur when the integrated &
                     &Autocorrelation (ACT) of the output MCMC chain is comparable to or larger than the length of the chain. &
                     &In such cases, neither the BatchMeans method nor any other method of ACT computation will be able to &
                     &accurately compute the ACT. Consequently, the samples generated based on the computed ACT values will &
                     &likely not be i.i.d. and will still be significantly autocorrelated. In such scenarios, more than &
                     &one refinement of the MCMC chain will be necessary. Very small sample size resulting from multiple &
                     &refinements of the sample could be a strong indication of the bad mixing of the MCMC chain and &
                     &the lack of convergence to the target objective function."//NL2//&
            SKC_"+   if `outputSampleRefinementCount > 1`,"//NL2//&
            SKC_"    the refinement of the output MCMC chain will be done for a maximum `outputSampleRefinementCount` number of &
                     &times, even though some residual autocorrelation in the final output MCMC sample may still exist."//NL2//&
            SKC_"+   if `outputSampleRefinementCount >> 1` (e.g., comparable to or larger than the length of the MCMC chain),"//NL2//&
            SKC_"    the refinement of the output MCMC chain will continue until the integrated autocorrelation of the resulting &
                     &final sample is less than 2, virtually implying that an independent identically-distributed (i.i.d.) sample &
                     &from the target objective function has finally been obtained."//NL2//&
            SKC_"Note that to obtain i.i.d. samples from a multidimensional chain, the sampler will, by default, use the maximum of &
                &integrated Autocorrelation (ACT) among all chain dimensions to refine the chain. &
                &Note that the value specified for `outputSampleRefinementCount` is used only when the variable outputSampleSize < 0, &
                &otherwise, it will be ignored. The default value is `outputSampleRefinementCount = "//getStr(spec%outputSampleRefinementCount%def)//SKC_"`."
            !!$omp master
            outputSampleRefinementCount = spec%outputSampleRefinementCount%null
            !!$omp end master
        end block outputSampleRefinementCount_block

        outputSampleRefinementMethod_block: block
            use pm_sampling_scio, only: outputSampleRefinementMethod
            spec%outputSampleRefinementMethod%def = spec%outputSampleRefinementMethod%batchMeans
            spec%outputSampleRefinementMethod%null = repeat(SUB, len(outputSampleRefinementMethod, IK))
            spec%outputSampleRefinementMethod%desc = &
            SKC_"The simulation specification `outputSampleRefinementMethod` is a scalar string of maximum length "//getStr(len(outputSampleRefinementMethod, IK))//&
            SKC_" representing the method of computing the (integrated) AutoCorrelation Time (ACT) to be used in the simulation for refining &
                &the final output MCMC chain and sample. If specified within an external input file, it must be either singly or doubly quoted. &
                &Methods that are currently supported include:"//NL2//&
            SKC_"+  `outputSampleRefinementMethod = '"//spec%outputSampleRefinementMethod%batchMeans//SKC_"'`"//NL2//&
            SKC_"    This method of computing the integrated Autocorrelation Time is based on the approach described in &
                     &**SCHMEISER, B., 1982, Batch size effects in the analysis of simulation output, Oper. Res. 30 556-568**. The &
                     &batch sizes in the BatchMeans method are chosen to be `int(N^(2/3))`, where `N` is the length of the MCMC chain. &
                     &As long as the batch size is larger than the ACT of the chain and there are significantly more than 10 &
                     &batches, the BatchMeans method will provide reliable estimates of the ACT. &
                     &Note that the refinement strategy involves two separate phases of sample decorrelation. At the first stage, &
                     &the Markov chain is decorrelated recursively (for as long as needed) based on the ACT of its compact format, &
                     &where only the uniquely visited states are kept in the (compact) chain. Once the Markov chain is refined &
                     &such that its compact format is fully decorrelated, the second phase of the decorrelation begins, during which &
                     &the Markov chain is decorrelated based on the ACT of the chain in its verbose (Markov) format. This process &
                     &is repeated recursively for as long as residual autocorrelation exists in the refined sample."//NL2//&
            SKC_"+   `outputSampleRefinementMethod = '"//spec%outputSampleRefinementMethod%batchMeans//SKC_"-compact'`"//NL2//&
            SKC_"    This is the same as the first case in the above, except that only the first phase of the sample refinement &
                     &described in the above will be performed; that is, the (verbose) Markov chain is refined only based on the &
                     &ACT computed from the compact format of the Markov chain. This will lead to a larger, final, refined sample. &
                     &However, the final sample will likely not be fully decorrelated."//NL2//&
            SKC_"+   `outputSampleRefinementMethod = '"//spec%outputSampleRefinementMethod%batchMeans//SKC_"-verbose'`"//NL2//&
            SKC_"    This is the same as the first case in the above, except that only the second phase of the sample refinement &
                     &described in the above will be performed; that is, the (verbose) Markov chain is refined only based on the ACT &
                     &computed from the verbose format of the Markov chain. While the resulting refined sample will be fully decorrelated, &
                     &the size of the refined sample may be smaller than the default choice in the first case in the above."//NL2//&
            SKC_"Note that to obtain i.i.d. samples from a multidimensional chain, the sampler will use the average of the &
                &ACT among all dimensions of the chain to refine the chain. If the maximum, minimum, or the median of IACs is preferred &
                &add `'-max'` (or `'-maximum'`), `'-min'` (or `'-minimum'`), `'-med'` (or `'-median'`), respectively, to the value of &
                &`outputSampleRefinementMethod`. For example, "//NL2//&
            SKC_"+   `outputSampleRefinementMethod = '"//spec%outputSampleRefinementMethod%batchMeans//SKC_"-max'`"//NL2//&
            SKC_"or, "//NL2//&
            SKC_"+   `outputSampleRefinementMethod = '"//spec%outputSampleRefinementMethod%batchMeans//SKC_"-compact-max'`"//NL2//&
            SKC_"or, "//NL2//&
            SKC_"    `outputSampleRefinementMethod = '"//spec%outputSampleRefinementMethod%batchMeans//SKC_"-max-compact'`"//NL2//&
            SKC_"Note that the specified `outputSampleRefinementCount` is used only when the condition `outputSampleSize < 0` holds. &
                &Otherwise, it is ignored. The default value is `outputSampleRefinementMethod = '"//spec%outputSampleRefinementMethod%def//SKC_"'`. &
                &Note that the input values are case-INsensitive and white-space characters are ignored."
            !!$omp master
            outputSampleRefinementMethod = spec%outputSampleRefinementMethod%null
            !!$omp end master
        end block outputSampleRefinementMethod_block

        proposal_block: block
            use pm_sampling_scio, only: proposal
            spec%proposal%is%uniform = .false._LK
            spec%proposal%is%normal = .false._LK
            spec%proposal%uniform = "uniform"
            spec%proposal%normal = "normal"
            spec%proposal%def = spec%proposal%normal
            spec%proposal%null = repeat(SUB, len(proposal, IK))
            spec%proposal%desc = &
            SKC_"The simulation specification `proposal` is a scalar string of maximum length "//getStr(len(proposal, IK))//SKC_" containing &
                &the name of the proposal distribution for the MCMC sampler. When specified from within an external input file, it must &
                &be singly or doubly quoted. Options that are currently supported include:"//NL2//&
            SKC_"+   `proposal = '"//spec%proposal%normal//SKC_"'`"//NL2//&
            SKC_"    This is equivalent to the multivariate normal distribution, &
                     &which is the most widely-used proposal model along with MCMC samplers."//NL2//&
            SKC_"+   `proposal = '"//spec%proposal%uniform//SKC_"'`"//NL2//&
            SKC_"    The proposals will be drawn uniformly from within a ndim-dimensional ellipsoid whose covariance matrix &
                     &and scale are initialized by the user and optionally adaptively updated throughout the simulation."//NL2//&
            SKC_"The default value is `'"//spec%proposal%def//SKC_"'`."
            !!$omp master
            proposal = spec%proposal%null
            !!$omp end master
        end block proposal_block

        proposalCorMat_block: block
            use pm_except, only: setNAN
            use pm_sampling_scio, only: proposalCorMat
            use pm_matrixInit, only: getMatInit, uppLowDia
            spec%proposalCorMat%def = getMatInit([ndim, ndim], uppLowDia, 0._RKC, 0._RKC, 1._RKC)
            spec%proposalCorMat%desc = &
            SKC_"The simulation specification `proposalCorMat` is a positive-definite square matrix of type `real` of the highest precision &
                &available within the ParaMonte library matrix of size `(ndim, ndim)`, where `ndim` is the dimension of the sampling space. &
                &It serves as the best-guess starting correlation matrix of the proposal distribution used by the sampler. &
                &It is used (along with the input vector `proposalStdVec`) to construct the covariance matrix of the proposal &
                &distribution when the input covariance matrix is missing in the input list of variables. &
                &If the covariance matrix is specified as input to the sampler, any input values for `proposalCorMat`, &
                &and `proposalStdVec` will be automatically ignored. Specifying `proposalCorMat` along with `proposalStdVec` is especially &
                &useful when obtaining the best-guess covariance matrix is not trivial. &
                &The default value of `proposalCorMat` is a square Identity matrix of rank `ndim`."
            !!$omp master
            if (allocated(proposalCorMat)) deallocate(proposalCorMat)
            allocate(proposalCorMat(ndim, ndim))
            call setNAN(proposalCorMat)
            !!$omp end master
        end block proposalCorMat_block

        proposalCovMat_block: block
            use pm_except, only: setNAN
            use pm_sampling_scio, only: proposalCovMat
            use pm_matrixInit, only: getMatInit, uppLowDia
            spec%proposalCovMat%def = getMatInit([ndim, ndim], uppLowDia, 0._RKC, 0._RKC, 1._RKC)
            ! This must be set here. It is important for the proper setting from inputFile and inputArg.
            spec%proposalCovMat%isUserSet = .false._LK
            spec%proposalCovMat%desc = &
            SKC_"The simulation specification `proposalCovMat` is a square positive-definite matrix of type `real` of the highest precision &
                &available within the ParaMonte library, of shape `(ndim, ndim)`, where `ndim` is the number of dimensions of the sampling space. &
                &It serves as the best-guess starting covariance matrix of the proposal distribution. &
                &To bring the sampling efficiency of the sampler to within the desired requested range, the covariance matrix will &
                &be adaptively updated throughout the simulation, according to the user-specified schedule. If `proposalCovMat` &
                &is not provided by the user or it is completely missing from the input file, its value will be automatically &
                &computed via the input variables `proposalCorMat` and `proposalStdVec` (or via their default values, if not provided). &
                &If the simulation specification `outputStatus` is set to ""extend"" and a successful prior simulation run exists, &
                &then `proposalCovMat` will be set to the covariance matrix of the output sample from the most recent simulation run. &
                &In this case, the computed `proposalCovMat` will override any user-specified value. &
                &Otherwise, the default value of `proposalCovMat` is a square Identity matrix of rank `ndim`."
            !!$omp master
            if (allocated(proposalCovMat)) deallocate(proposalCovMat)
            allocate(proposalCovMat(ndim, ndim))
            call setNAN(proposalCovMat)
            !!$omp end master
        end block proposalCovMat_block

        proposalScaleFactor_block: block
            use pm_sampling_scio, only: proposalScaleFactor
            spec%proposalScaleFactor%strdef = SKC_"gelman"
            spec%proposalScaleFactor%valdef = 2.38_RKC / sqrt(real(ndim, RKC)) ! Gelman, Roberts, Gilks (1996): Efficient Metropolis Jumping Rules.
            spec%proposalScaleFactor%null = repeat(SUB, len(proposalScaleFactor, IK))
            spec%proposalScaleFactor%desc = &
            SKC_"The simulation specification `proposalScaleFactor` is a scalar string of maximum length "//getStr(len(proposalScaleFactor, IK))//SKC_" &
                &containing a positive real-valued number whose square will be multiplied with the covariance matrix of the proposal distribution &
                &of the MCMC sampler to shrink or enlarge it. In other words, the proposal distribution will be scaled in every direction &
                &by the specified numeric value of `proposalScaleFactor`. It can also be given in units of the string keyword 'gelman' &
                &(which is case-INsensitive) after the paper:"//NL2//&
            SKC_"    Gelman, Roberts, and Gilks (1996): 'Efficient Metropolis Jumping Rules'."//NL2//&
            SKC_"The paper finds that the optimal scaling factor for a Multivariate Gaussian proposal distribution for the Metropolis-Hastings &
                &Markov Chain Monte Carlo sampling of a target Multivariate Normal Distribution of dimension `ndim` is given by:"//NL2//&
            SKC_"    proposalScaleFactor = 2.38 / sqrt(ndim)  ,  in the limit of ndim -> Infinity."//NL2//&
            SKC_"Multiples of the Gelman scale factors are also acceptable as input and can be specified like the following examples:"//NL2//&
            SKC_"+   `proposalScaleFactor = '1'`"//NL2//&
            SKC_"    multiplies the ndim-dimensional proposal covariance matrix by 1, &
                     &essentially no change occurs to the covariance matrix."//NL2//&
            SKC_"+   `proposalScaleFactor = ""1""`"//NL2//&
            SKC_"    same as the previous example. The double-quotation marks act the same way as single-quotation marks."//NL2//&
            SKC_"+   `proposalScaleFactor = '2.5'`"//NL2//&
            SKC_"    multiplies the ndim-dimensional proposal covariance matrix by 2.5."//NL2//&
            SKC_"+   `proposalScaleFactor = '2.5*Gelman'`"//NL2//&
            SKC_"    multiplies the `ndim`-dimensional proposal covariance matrix by 2.5 * 2.38/sqrt(ndim)."//NL2//&
            SKC_"+   `proposalScaleFactor = ""2.5 * gelman""`"//NL2//&
            SKC_"    same as the previous example but with double-quotation marks. space characters are ignored."//NL2//&
            SKC_"+   `proposalScaleFactor = ""2.5 * gelman*gelman*2""`"//NL2//&
            SKC_"    equivalent to gelmanFactor-squared multiplied by `5`."//NL2//&
            SKC_"Note, however, that the result of Gelman et al. paper applies only to multivariate normal proposal distributions, in the limit &
                &of infinite dimensions. Therefore, care must be taken when using Gelman's scaling factor with non-Gaussian proposals and target &
                &objective functions. Only the product symbol `*` can be parsed in the string value of `proposalScaleFactor`. &
                &The presence of other mathematical symbols or multiple appearances of the product symbol will lead to a simulation crash. &
                &Also, note that the prescription of an acceptance range specified by the input variable `targetAcceptanceRate` will lead &
                &to dynamic modification of the initial input value of `proposalScaleFactor` throughout sampling for `proposalAdaptationCount` times. &
                &The default string value for `proposalScaleFactor` is ""gelman"" (for all proposal distributions), &
                &which is subsequently converted to `2.38 / sqrt(ndim)`."
            !!$omp master
            proposalScaleFactor = spec%proposalScaleFactor%null
            !!$omp end master
        end block proposalScaleFactor_block

        proposalStart_block: block
            use pm_sampling_scio, only: proposalStart
            spec%proposalStart%desc = &
            SKC_"The simulation specification `proposalStart` is a vector type of `real` of the highest precision available &
                &within the ParaMonte library of length `ndim` where `ndim` is the dimension of the domain of the objective function. &
                &For every element of `proposalStart` that is not provided as input, the default value will be the center of the sampling domain &
                &as determined by `domainCubeLimitLower` and `domainCubeLimitUpper` input specifications. &
                &If the condition `proposalStartRandomized` is set to the logical/Boolean true value, then the missing &
                &elements of `proposalStart` will be initialized to values drawn randomly from within the corresponding &
                &ranges specified by the input variables `proposalStartDomainCubeLimitLower` and `proposalStartDomainCubeLimitUpper`. &
                &If the simulation specification `outputStatus` is set to ""extend"" and a prior successful simulation run exists, &
                &then `proposalStart` will be set to the average of the sampled states of the most recent successful run. &
                &In this case, any user-specified value will be overridden by the computed `proposalStart`."
            !!$omp master
            if (allocated(proposalStart)) deallocate(proposalStart)
            allocate(proposalStart(ndim))
            call setNAN(proposalStart)
            !!$omp end master
        end block proposalStart_block

        proposalStartDomainCubeLimitLower_block: block
            use pm_sampling_scio, only: proposalStartDomainCubeLimitLower
            spec%proposalStartDomainCubeLimitLower%def = spec%domainCubeLimitLower%def
            spec%proposalStartDomainCubeLimitLower%desc = &
            SKC_"The simulation specification `proposalStartDomainCubeLimitLower` is a vector of type `real` of the highest precision available &
                &in the ParaMonet library of size `ndim` is the number of dimensions of the domain of the objective function. It contains the &
                &lower boundaries of the cubical domain from which the starting point(s) of the MCMC chain(s) will be initialized randomly &
                &(only if requested via the input variable `proposalStartRandomized`). &
                &This happens only when some or all of the input specification `proposalStart` elements are missing. &
                &In such cases, every missing value of the input `proposalStart` will be set to the center point between &
                &`proposalStartDomainCubeLimitLower` and `proposalStartDomainCubeLimitUpper` in the corresponding dimension. &
                &If `proposalStartRandomized` is set to the logical/Boolean true value, then the missing elements of &
                &`proposalStart` will be initialized to values drawn randomly from within the corresponding ranges &
                &whose lower limits are specified by the input `proposalStartDomainCubeLimitLower`. &
                &When specified from within an external input file to the sampler, it is also possible to assign only select &
                &values of `proposalStartDomainCubeLimitLower` and leave the rest of the components to be assigned the default value. &
                &For example, having the following inside the input file, "//NL2//&
            SKC_"+   `proposalStartDomainCubeLimitLower(3:5) = -100`"//NL2//&
            SKC_"    will only set the lower limits of the third, fourth, and the fifth dimensions to -100, or,"//NL2//&
            SKC_"+   `proposalStartDomainCubeLimitLower(1) = -100, proposalStartDomainCubeLimitLower(2) = -1.e6`"//NL2//&
            SKC_"    will set the lower limit on the first dimension to -100, and 1.e6 on the second dimension, or,"//NL2//&
            SKC_"+   `proposalStartDomainCubeLimitLower = 3*-2.5e100`"//NL2//&
            SKC_"    will only set the lower limits on the first, second, and the third dimensions to `-2.5*10^100`, &
                     &while the rest of the lower limits for the missing dimensions will be automatically set to the default value."//NL2//&
            SKC_"The default for all `proposalStartDomainCubeLimitLower` elements are taken from the corresponding elements of `domainCubeLimitLower`."
            !!$omp master
            call setResized(proposalStartDomainCubeLimitLower, ndim)
            call setNAN(proposalStartDomainCubeLimitLower)
            !!$omp end master
        end block proposalStartDomainCubeLimitLower_block

        proposalStartDomainCubeLimitUpper_block: block
            use pm_sampling_scio, only: proposalStartDomainCubeLimitUpper
            spec%proposalStartDomainCubeLimitUpper%def = spec%domainCubeLimitUpper%def
            spec%proposalStartDomainCubeLimitUpper%desc = &
            SKC_"The simulation specification `proposalStartDomainCubeLimitUpper` is a vector of type `real` of the highest precision available &
                &in the ParaMonet library of size `ndim` is the number of dimensions of the domain of the objective function. It contains the &
                &upper boundaries of the cubical domain from which the starting point(s) of the MCMC chain(s) will be initialized randomly &
                &(only if requested via the input variable `proposalStartRandomized`). &
                &This happens only when some or all of the input specification `proposalStart` elements are missing. &
                &In such cases, every missing value of the input `proposalStart` will be set to the center point between &
                &`proposalStartDomainCubeLimitLower` and `proposalStartDomainCubeLimitUpper` in the corresponding dimension. &
                &If `proposalStartRandomized` is set to the logical/Boolean true value, then the missing elements of &
                &`proposalStart` will be initialized to values drawn randomly from within the corresponding ranges &
                &whose upper limits are specified by the input `proposalStartDomainCubeLimitUpper`. &
                &When specified from within an external input file to the sampler, it is also possible to assign only select &
                &values of `proposalStartDomainCubeLimitUpper` and leave the rest of the components to be assigned the default value. &
                &For example, having the following inside the input file, "//NL2//&
            SKC_"+   `proposalStartDomainCubeLimitUpper(3:5) = -100`"//NL2//&
            SKC_"    will only set the upper limits of the third, fourth, and the fifth dimensions to `-100`, or,"//NL2//&
            SKC_"+   `proposalStartDomainCubeLimitUpper(1) = -100, proposalStartDomainCubeLimitUpper(2) = -1.e6`"//NL2//&
            SKC_"    will set the upper limit on the first dimension to -100, and 1.e6 on the second dimension, or,"//NL2//&
            SKC_"+   `proposalStartDomainCubeLimitUpper = 3*-2.5e100`"//NL2//&
            SKC_"    will only set the upper limits on the first, second, and the third dimensions to `-2.5*10**100`, while &
                     &the rest of the upper limits for the missing dimensions will be automatically set to the default value."//NL2//&
            SKC_"The default values for all elements of proposalStartDomainCubeLimitUpper are &
                &taken from the corresponding values in the input variable `domainCubeLimitUpper`."
            !!$omp master
            call setResized(proposalStartDomainCubeLimitUpper, ndim)
            call setNAN(proposalStartDomainCubeLimitUpper)
            !!$omp end master
        end block proposalStartDomainCubeLimitUpper_block

        proposalStartRandomized_block: block
            use pm_sampling_scio, only: proposalStartRandomized
            spec%proposalStartRandomized%def = .false._LK
            spec%proposalStartRandomized%desc = &
            SKC_"The simulation specification `proposalStartRandomized` is scalar of type `logical` (Boolean). &
                &If `true` (or `.true.` or `TRUE` or `.t.` from within an external input file), then the variable `proposalStart` &
                &will be initialized randomly for each MCMC chain that is to be generated by the sampler. The random values will be &
                &drawn from the specified or the default domain of `proposalStart`, given by `proposalStartDomainCubeLimitLower` and &
                &`proposalStartDomainCubeLimitUpper` variable. Note that the value of `proposalStart`, if provided, has precedence over &
                &random initialization. In other words, only uninitialized elements of `proposalStart` will be randomly initialized only &
                &if `proposalStartRandomized` is set to the logical true value. Note that even if `proposalStart` is randomly initialized, &
                &its random value will be deterministic between different independent simulation runs if the input variable `randomSeed` &
                &is specified by the user. The default value is `"//getStr(spec%proposalStartRandomized%def)//SKC_"`."
            !!$omp master
            proposalStartRandomized = spec%proposalStartRandomized%def
            !!$omp end master
        end block proposalStartRandomized_block

        proposalStdVec_block: block
            use pm_arrayFill, only: getFilled
            use pm_sampling_scio, only: proposalStdVec
            spec%proposalStdVec%def = getFilled(1._RKC, ndim)
            spec%proposalStdVec%desc = &
            SKC_"The simulation specification `proposalStdVec` is a positive-valued vector of type `real` of the highest precision available &
                &within the ParaMonte library, of size `ndim`, where `ndim` is the dimension of the domain of the objective function. &
                &It serves as the best-guess starting standard deviation for each component of the proposal distribution. &
                &If the initial covariance matrix (`proposalCovMat`) is missing as an input specification to the sampler, &
                &then `proposalStdVec` (along with the specified `proposalCorMat`) will be used to construct &
                &the initial covariance matrix of the proposal distribution of the MCMC sampler. &
                &However, if `proposalCovMat` is specified for the sampler, then the input `proposalStdVec` and `proposalCorMat` will &
                &be completely ignored, and the input value for `proposalCovMat` will be used to construct the initial covariance matrix &
                &of the proposal distribution. The default value of `proposalStdVec` is a vector of size `ndim` of unit values (i.e., ones)."
            !!$omp master
            call setResized(proposalStdVec, ndim)
            call setNAN(proposalStdVec)
            !!$omp end master
        end block proposalStdVec_block

        !!$omp barrier

    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function set(spec, sampler) result(err)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: set
#endif
        use pm_err, only: err_type
        use pm_sampling, only: paramcmc_type
        use pm_sampling, only: sampler_type
        class(specmcmc_type), intent(inout) :: spec
        class(sampler_type), intent(in), optional :: sampler
        type(err_type) :: err

        character(*, SK), parameter :: PROCEDURE_NAME = MODULE_NAME//"@set()"

        select type(sampler)
        type is (paramcmc_type)

            err = spec%specbase_type%set(sampler%sampler_type)
            if (err%occurred) return

            outputChainSize_block: block
                use pm_sampling_scio, only: outputChainSize
                if (spec%overridable .and. allocated(sampler%outputChainSize)) then
                    spec%outputChainSize%val = sampler%outputChainSize
                else
                    spec%outputChainSize%val = outputChainSize
                end if
                if (spec%outputChainSize%val == spec%outputChainSize%null) spec%outputChainSize%val = spec%outputChainSize%def
            end block outputChainSize_block

            outputSampleRefinementCount_block: block
                use pm_sampling_scio, only: outputSampleRefinementCount
                if (spec%overridable .and. allocated(sampler%outputSampleRefinementCount)) then
                    spec%outputSampleRefinementCount%val = sampler%outputSampleRefinementCount
                else
                    spec%outputSampleRefinementCount%val = outputSampleRefinementCount
                end if
                if (spec%outputSampleRefinementCount%val == spec%outputSampleRefinementCount%null) spec%outputSampleRefinementCount%val = spec%outputSampleRefinementCount%def
                spec%outputSampleRefinementCount%str = getStr(spec%outputSampleRefinementCount%val)
            end block outputSampleRefinementCount_block

            outputSampleRefinementMethod_block: block
                use pm_sampling_scio, only: outputSampleRefinementMethod
                character(:,SKC), allocatable :: lowval
                if (spec%overridable .and. allocated(sampler%outputSampleRefinementMethod)) then
                    spec%outputSampleRefinementMethod%val = trim(adjustl(getRemoved(sampler%outputSampleRefinementMethod, SKC_" ")))
                else
                    spec%outputSampleRefinementMethod%val = trim(adjustl(getRemoved(outputSampleRefinementMethod, SKC_" ")))
                end if
                if (spec%outputSampleRefinementMethod%val == trim(adjustl(spec%outputSampleRefinementMethod%null))) spec%outputSampleRefinementMethod%val = spec%outputSampleRefinementMethod%def
                lowval = getStrLower(spec%outputSampleRefinementMethod%val)
                spec%outputSampleRefinementMethod%isMaxCumSumACF = logical(lowval == SKC_"maxcumsumacf", LK)
                spec%outputSampleRefinementMethod%isBatchMeans = logical(0 < index(lowval, SKC_"batchmeans"), LK)
                spec%outputSampleRefinementMethod%isBatchMeansCompact = spec%outputSampleRefinementMethod%isBatchMeans .and. logical(0 < index(lowval, SKC_"compact"), LK)
                spec%outputSampleRefinementMethod%isBatchMeansVerbose = spec%outputSampleRefinementMethod%isBatchMeans .and. logical(0 < index(lowval, SKC_"verbose"), LK)
            end block outputSampleRefinementMethod_block

            proposal_block: block
                use pm_sampling_scio, only: proposal
                if (spec%overridable .and. allocated(sampler%proposal)) then
                    spec%proposal%val = getStrLower(trim(adjustl(sampler%proposal)))
                else
                    spec%proposal%val = getStrLower(trim(adjustl(proposal)))
                end if
                if (spec%proposal%val == spec%proposal%null) spec%proposal%val = spec%proposal%def
                spec%proposal%is%uniform = spec%proposal%val == spec%proposal%uniform
                spec%proposal%is%normal = spec%proposal%val == spec%proposal%normal
            end block proposal_block

            proposalCorMat_block: block
                use pm_sampling_scio, only: proposalCorMat
                if (spec%overridable .and. allocated(sampler%proposalCorMat)) then
                    spec%proposalCorMat%val = real(sampler%proposalCorMat, RKC)
                else
                    spec%proposalCorMat%val = proposalCorMat
                end if
                where (isNAN(spec%proposalCorMat%val))
                    spec%proposalCorMat%val = spec%proposalCorMat%def
                end where
            end block proposalCorMat_block

            proposalCovMat_block: block
                use pm_sampling_scio, only: proposalCovMat
                integer(IK) :: idim, jdim
                if (spec%overridable .and. allocated(sampler%proposalCovMat)) then
                    spec%proposalCovMat%val = real(sampler%proposalCovMat, RKC)
                else
                    spec%proposalCovMat%val = proposalCovMat
                end if
                do jdim = 1, spec%ndim%val
                    do idim = 1, spec%ndim%val
                        if (isNAN(spec%proposalCovMat%val(idim, jdim))) then
                            spec%proposalCovMat%val(idim, jdim) = spec%proposalCovMat%def(idim, jdim)
                        else
                            spec%proposalCovMat%isUserSet = .true._LK
                        end if
                    end do
                end do
            end block proposalCovMat_block

            proposalScaleFactor_block: block
                use pm_sampling_scio, only: proposalScaleFactor
                if (spec%overridable .and. allocated(sampler%proposalScaleFactor)) then
                    spec%proposalScaleFactor%str = trim(adjustl(sampler%proposalScaleFactor))
                else
                    spec%proposalScaleFactor%str = trim(adjustl(proposalScaleFactor))
                end if
                if (spec%proposalScaleFactor%str == spec%proposalScaleFactor%null) then
                    spec%proposalScaleFactor%val = spec%proposalScaleFactor%valdef
                    spec%proposalScaleFactor%str = spec%proposalScaleFactor%strdef
                end if
            end block proposalScaleFactor_block

            proposalStart_block: block
                use pm_sampling_scio, only: proposalStart
                if (spec%overridable .and. allocated(sampler%proposalStart)) then
                    spec%proposalStart%val = real(sampler%proposalStart, RKC)
                else
                    spec%proposalStart%val = proposalStart
                end if
            end block proposalStart_block

            proposalStartDomainCubeLimitLower_block: block
                use pm_sampling_scio, only: proposalStartDomainCubeLimitLower
                if (spec%overridable .and. allocated(sampler%proposalStartDomainCubeLimitLower)) then
                    spec%proposalStartDomainCubeLimitLower%val = real(sampler%proposalStartDomainCubeLimitLower, RKC)
                else
                    spec%proposalStartDomainCubeLimitLower%val = proposalStartDomainCubeLimitLower
                end if
                where(isNAN(spec%proposalStartDomainCubeLimitLower%val))
                    spec%proposalStartDomainCubeLimitLower%val = spec%domainCubeLimitLower%val
                end where
            end block proposalStartDomainCubeLimitLower_block

            proposalStartDomainCubeLimitUpper_block: block
                use pm_sampling_scio, only: proposalStartDomainCubeLimitUpper
                if (spec%overridable .and. allocated(sampler%proposalStartDomainCubeLimitUpper)) then
                    spec%proposalStartDomainCubeLimitUpper%val = real(sampler%proposalStartDomainCubeLimitUpper, RKC)
                else
                    spec%proposalStartDomainCubeLimitUpper%val = proposalStartDomainCubeLimitUpper
                end if
                where(isNAN(spec%proposalStartDomainCubeLimitUpper%val))
                    spec%proposalStartDomainCubeLimitUpper%val = spec%domainCubeLimitUpper%val
                end where
            end block proposalStartDomainCubeLimitUpper_block

            proposalStartRandomized_block: block
                use pm_sampling_scio, only: proposalStartRandomized
                if (spec%overridable .and. allocated(sampler%proposalStartRandomized)) then
                    spec%proposalStartRandomized%val = sampler%proposalStartRandomized
                else
                    spec%proposalStartRandomized%val = proposalStartRandomized
                end if
            end block proposalStartRandomized_block

            proposalStdVec_block: block
                use pm_sampling_scio, only: proposalStdVec
                if (spec%overridable .and. allocated(sampler%proposalStdVec)) then
                    spec%proposalStdVec%val = real(sampler%proposalStdVec, RKC)
                else
                    spec%proposalStdVec%val = proposalStdVec
                end if
                where (isNAN(spec%proposalStdVec%val))
                    spec%proposalStdVec%val = spec%proposalStdVec%def
                end where
            end block proposalStdVec_block

            ! Resolve the conflicting cases.

            if (spec%outputStatus%is%extend .and. 1 < spec%run%id) then
                block
                    use pm_sampleMean, only: setMean
                    use pm_sampleCov, only: setCov, uppDia
                    use pm_io, only: getErrTableRead, LEN_IOMSG
                    use pm_matrixCopy, only: setMatCopy, rdpack, upp, transHerm
                    real(RKC), allocatable :: logFuncState(:,:)
                    call setResized(err%msg, LEN_IOMSG)
                    err%stat = getErrTableRead(spec%sampleFile%list(spec%run%id - 1)%val, logFuncState, sep = spec%outputSeparator%val, roff = 1_IK, iomsg = err%msg)
                    err%occurred = err%stat /= 0
                    if (err%occurred) then
                        err%msg = PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SK_": Failed to read the sample file """//spec%sampleFile%list(spec%run%id - 1)%val//& ! LCOV_EXCL_LINE
                        SK_""" contents from the previous simulation run. "//trim(err%msg) ! LCOV_EXCL_LINE
                        return ! LCOV_EXCL_LINE
                    end if
                    call setMean(spec%proposalStart%val, logFuncState(:, 2 : spec%ndim%val + 1), dim = 1_IK)
                    call setCov(spec%proposalCovMat%val, uppDia, logFuncState(:, 2 : spec%ndim%val + 1), dim = 1_IK)
                    call setMatCopy(spec%proposalCovMat%val, rdpack, spec%proposalCovMat%val, rdpack, upp, transHerm)
                end block
            else
                block
                    use pm_sampleCov, only: getCov, uppDia
                    use pm_distUnif, only: setUnifRand
                    integer(IK) :: idim
                    if (.not. spec%proposalCovMat%isUserSet) spec%proposalCovMat%val = getCov(spec%proposalCorMat%val, uppDia, spec%proposalStdVec%val)
                    do idim = 1, size(spec%proposalStart%val, 1, IK)
                        if (isNAN(spec%proposalStart%val(idim))) then
                            if (spec%proposalStartRandomized%val) then
                                call setUnifRand(spec%rng, spec%proposalStart%val(idim), spec%proposalStartDomainCubeLimitLower%val(idim), spec%proposalStartDomainCubeLimitUpper%val(idim))
                            else
                                spec%proposalStart%val(idim) = 0.5_RKC * spec%proposalStartDomainCubeLimitLower%val(idim) + 0.5_RKC * spec%proposalStartDomainCubeLimitUpper%val(idim)
                            end if
                        end if
                    end do
                end block
            end if

            ! open output files, report and sanitize.

            if (spec%image%is%leader) call spec%report() ! if (spec%run%is%new)
            call spec%sanitize(err)

        class default
            error stop "The input `sampler` must be of type `paramcmc_type`." ! LCOV_EXCL_LINE
        end select

    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine report(spec)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: report
#endif
        class(specmcmc_type), intent(inout) :: spec

        call spec%disp%text%wrap(NL1//spec%method%val//SKC_".simulation.specifications.mcmc"//NL1)

        associate(ndim => spec%ndim%val, format => spec%reportFile%format%generic)

            call spec%disp%show("outputChainSize")
            call spec%disp%show(spec%outputChainSize%val, format = format)
            call spec%disp%note%show(spec%outputChainSize%desc)

            call spec%disp%show("outputSampleRefinementCount")
            call spec%disp%show(spec%outputSampleRefinementCount%val, format = format)
            call spec%disp%note%show(spec%outputSampleRefinementCount%desc)

            call spec%disp%show("outputSampleRefinementMethod")
            call spec%disp%show(spec%outputSampleRefinementMethod%val, format = format)
            call spec%disp%note%show(spec%outputSampleRefinementMethod%desc)

            call spec%disp%show("proposal")
            call spec%disp%show(spec%proposal%val, format = format)
            call spec%disp%note%show(spec%proposal%desc)

            call spec%disp%show("proposalCorMat")
            call spec%disp%show(spec%proposalCorMat%val, format = format)
            call spec%disp%note%show(spec%proposalCorMat%desc)

            call spec%disp%show("proposalCovMat")
            call spec%disp%show(spec%proposalCovMat%val, format = format)
            call spec%disp%note%show(spec%proposalCovMat%desc)

            call spec%disp%show("proposalScaleFactor")
            call spec%disp%show(spec%proposalScaleFactor%str//SKC_" (gelman(ndim = "//getStr(spec%ndim%val)//SKC_") = "//getStr(spec%proposalScaleFactor%valdef)//SKC_")", format = format)
            call spec%disp%note%show(spec%proposalScaleFactor%desc)

            call spec%disp%show("proposalStart")
            call spec%disp%show(reshape(spec%proposalStart%val, [ndim, 1_IK]), format = format)
            call spec%disp%note%show(spec%proposalStart%desc)

            call spec%disp%show("proposalStartDomainCubeLimitLower")
            call spec%disp%show(reshape(spec%proposalStartDomainCubeLimitLower%val, [ndim, 1_IK]), format = format)
            call spec%disp%note%show(spec%proposalStartDomainCubeLimitLower%desc)

            call spec%disp%show("proposalStartDomainCubeLimitUpper")
            call spec%disp%show(reshape(spec%proposalStartDomainCubeLimitUpper%val, [ndim, 1_IK]), format = format)
            call spec%disp%note%show(spec%proposalStartDomainCubeLimitUpper%desc)

            call spec%disp%show("proposalStartRandomized")
            call spec%disp%show(spec%proposalStartRandomized%val, format = format)
            call spec%disp%note%show(spec%proposalStartRandomized%desc)

            call spec%disp%show("proposalStdVec")
            call spec%disp%show(reshape(spec%proposalStdVec%val, [ndim, 1_IK]), format = format)
            call spec%disp%note%show(spec%proposalStdVec%desc)

        end associate

    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine sanitize(spec, err)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: sanitize
#endif
        use pm_err, only: err_type
        type(err_type), intent(inout) :: err
        class(specmcmc_type), intent(inout) :: spec
        character(*,SKC), parameter :: PROCEDURE_NAME = MODULE_NAME//SKC_"@sanitize()"

        outputChainSize_block: block
            if (spec%outputChainSize%val < spec%ndim%val + 1_IK) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. The input requested value for `outputChainSize` ("//getStr(spec%outputChainSize%val)//&
                            SKC_") can neither be negative nor smaller than `ndim + 1`, where `ndim` represents the dimension of the domain of the objective &
                            &function, here `ndim = "//getStr(spec%ndim%val)//SKC_". If you do not know an appropriate value for `outputChainSize`, &
                            &drop it from the input list. The MCMC sampler will automatically assign an appropriate value to it."
            end if
        end block outputChainSize_block

        outputSampleRefinementCount_block: block
            if (spec%outputSampleRefinementCount%val < 0_IK) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The input value for variable `outputSampleRefinementCount` must be a non-negative integer. &
                            &If you are unsure about the appropriate value for this variable, simply drop it from the input. &
                            &The sampler will automatically assign an appropriate value to it."
            end if
        end block outputSampleRefinementCount_block

        outputSampleRefinementMethod_block: block
            character(:,SKC), allocatable :: lowerCaseVal
            lowerCaseVal = getStrLower(spec%outputSampleRefinementMethod%val)
            if  (& ! LCOV_EXCL_LINE
                index(lowerCaseVal, getStrLower(spec%outputSampleRefinementMethod%maxCumSumACF)) == 0_IK & ! LCOV_EXCL_LINE
                .and. & ! LCOV_EXCL_LINE
                index(lowerCaseVal, getStrLower(spec%outputSampleRefinementMethod%batchMeans)) == 0_IK & ! LCOV_EXCL_LINE
                .and. & ! LCOV_EXCL_LINE
                index(lowerCaseVal, getStrLower(spec%outputSampleRefinementMethod%cutoffACF)) == 0_IK & ! LCOV_EXCL_LINE
               ) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The input requested method for the computation of the integrated AutoCorrelation Time ("//spec%outputSampleRefinementMethod%val//&
                            SKC_") assigned to the variable `outputSampleRefinementMethod` cannot be anything other than "//spec%outputSampleRefinementMethod%batchMeans//SKC_". "//&
                            ! " or "//spec%maxCumSumACF//". &
                            SKC_"If you are unsure of the appropriate value for `outputSampleRefinementMethod`, &
                            &drop it from the input list. The sampler will automatically assign an appropriate value to it."
            end if
        end block outputSampleRefinementMethod_block

        proposal_block: block
            if (.not. (spec%proposal%is%normal .or. spec%proposal%is%uniform)) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The input requested proposal model ("//spec%proposal%val//") is not supported. &
                            &The variable proposal cannot be set to anything other than '"//&
                            spec%proposal%normal//SKC_"', or '"//&
                            spec%proposal%uniform//SKC_"'."
            end if
        end block proposal_block

        proposalCorMat_block: block
            if (.not. isMatClass(spec%proposalCorMat%val, posdefmat)) then
                err%occurred = .true._LK
                err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. The specified `proposalCorMat` is not positive-definite: "//getStr(spec%proposalCorMat%val)//SKC_""
            end if
        end block proposalCorMat_block

        proposalCovMat_block: block
            if (.not. isMatClass(spec%proposalCovMat%val, posdefmat)) then
                err%occurred = .true._LK
                err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. The specified `proposalCovMat` is not positive-definite: "//getStr(spec%proposalCovMat%val)//SKC_""
            end if
        end block proposalCovMat_block

        proposalScaleFactor_block: block
            use pm_val2real, only: setReal
            use pm_arraySplit, only: setSplit
            integer(IK), allocatable :: sindex(:,:)
            character(:,SKC), allocatable :: str
            real(RKC) :: conversion
            integer(IK) :: ipart
            ! First convert the proposalScaleFactor string to real value:
            str = getRemoved(spec%proposalScaleFactor%str, SKC_" ") ! remove the white spaces.
            if (len_trim(adjustl(str)) == 0_IK) then
                err%occurred = .true._LK
                err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                            &The input string value ("//spec%proposalScaleFactor%str//SKC_") for the variable `proposalScaleFactor` is empty. &
                            &Ensure the specified string follows the syntax rules of the sampler for this variable. Otherwise drop &
                            &it from the input list. The sampler will automatically assign an appropriate value to it."
            end if
            ! Now split the string by "*" to real coefficient and character (gelman) parts for further evaluations.
            call setSplit(sindex, str, sep = SKC_"*")
            spec%proposalScaleFactor%val = 1._RKC
            do ipart = 1, size(sindex, 2, IK)
                if (getStrLower(str(sindex(1, ipart) : sindex(2, ipart))) == SKC_"gelman") then
                    spec%proposalScaleFactor%val = spec%proposalScaleFactor%val * spec%proposalScaleFactor%valdef
                else
                    call setReal(conversion, str(sindex(1, ipart) : sindex(2, ipart)), iostat = err%stat)
                    if (err%stat /= 0_IK) then
                        err%occurred = .true._LK
                        err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred while reading real number. &
                        &The input string value for the variable `proposalScaleFactor` ("//spec%proposalScaleFactor%str//SKC_") does not appear to follow &
                        &the standard syntax rules of the sampler for this variable. The slice `"//str(sindex(1, ipart) : sindex(1, ipart))//& ! LCOV_EXCL_LINE
                        SKC_"` cannot be parsed into any meaningful token. Please correct the input value, or drop it from the input list &
                        &in which case, the sampler will automatically assign an appropriate value to it."
                    else
                        spec%proposalScaleFactor%val = spec%proposalScaleFactor%val * conversion
                    end if
                end if
            end do
            ! Now check if the real value is positive
            if (spec%proposalScaleFactor%val <= 0_IK) then
                err%occurred = .true._LK
                err%msg = err%msg//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                &The input string value ("""//spec%proposalScaleFactor%str//""") translates to a non-positive real value: "//getStr(spec%proposalScaleFactor%val)//SKC_". &
                &Make sure the input string follows the syntax rules of the sampler for this variable. &
                &Otherwise drop it from the input list. The sampler will automatically assign an appropriate value to it."
            end if
        end block proposalScaleFactor_block

        proposalStart_block: block
            use pm_sampling_scio, only: proposalStart
            !!$omp barrier
            !!$omp master
            if (allocated(proposalStart)) deallocate(proposalStart)
            !!$omp end master
        end block proposalStart_block

        proposalStartDomainCubeLimitLower_block: block
            use pm_sampling_scio, only: proposalStartDomainCubeLimitLower
            !!$omp barrier
            !!$omp master
            if (allocated(proposalStartDomainCubeLimitLower)) deallocate(proposalStartDomainCubeLimitLower)
            !!$omp end master
        end block proposalStartDomainCubeLimitLower_block

        proposalStartDomainCubeLimitUpper_block: block
            use pm_sampling_scio, only: proposalStartDomainCubeLimitUpper
            integer(IK) :: idim
            do idim = 1, size(spec%proposalStartDomainCubeLimitUpper%val, kind = IK)
                ! Check if the domain is set when random start point is requested.
                if (spec%proposalStartRandomized%val .and. spec%proposalStartDomainCubeLimitUpper%val(idim) == spec%domainCubeLimitUpper%def) then
                    err%occurred = .true._LK
                    err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                                &You have requested a random start point by setting `proposalStartRandomized` to TRUE while the &
                                &element #"//getStr(idim)//SKC_" of `proposalStartDomainCubeLimitLower` has not been preset to a finite value. &
                                &This information is essential otherwise, how could the sampler draw points randomly from within an unspecified domain?"
                end if
                ! The upper boundary of the domain of random-start-point must be smaller than the upper boundary of the target's domain.
                if (spec%domainCubeLimitUpper%val(idim) < spec%proposalStartDomainCubeLimitUpper%val(idim)) then
                    err%occurred = .true._LK
                    err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. The component "//getStr(idim)//& ! LCOV_EXCL_LINE
                                SKC_" of the variable `proposalStartDomainCubeLimitUpper` ("//getStr(spec%proposalStartDomainCubeLimitUpper%val(idim))//SKC_") cannot be "//& ! LCOV_EXCL_LINE
                                SKC_"larger than the corresponding component of the variable domainCubeLimitUpper ("//& ! LCOV_EXCL_LINE
                                    getStr(spec%domainCubeLimitUpper%val(idim))//SKC_"). If you do not know an appropriate &
                                    &value to set for `proposalStartDomainCubeLimitUpper`, drop it from the input list. &
                                    &The sampler will automatically assign an appropriate value to it."
                end if
                ! The upper boundary of the domain of random-start-point must be smaller than the corresponding lower boundary.
                if (spec%proposalStartDomainCubeLimitUpper%val(idim) <= spec%proposalStartDomainCubeLimitLower%val(idim)) then
                    err%occurred = .true._LK
                    err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. The input upper limit value in the component "//getStr(idim)//&
                                SKC_" of the variable `proposalStartDomainCubeLimitUpper` cannot be smaller than or equal to the corresponding input &
                                    &lower limit value in `proposalStartDomainCubeLimitLower`:"//NL1//&
                                SKC_"    proposalStartDomainCubeLimitLower("//getStr(idim)//SKC_") = "//getStr(spec%proposalStartDomainCubeLimitLower%val(idim))//NL1//&
                                SKC_"    proposalStartDomainCubeLimitUpper("//getStr(idim)//SKC_") = "//getStr(spec%proposalStartDomainCubeLimitUpper%val(idim))//SKC_""
                end if

            end do
            !!$omp barrier
            !!$omp master
            if (allocated(proposalStartDomainCubeLimitUpper)) deallocate(proposalStartDomainCubeLimitUpper)
            !!$omp end master
        end block proposalStartDomainCubeLimitUpper_block

        proposalStartRandomized_block: block
        end block proposalStartRandomized_block

        proposalStdVec_block: block
            integer(IK) :: idim
            do idim = 1, spec%ndim%val
                if (spec%proposalStdVec%val(idim) <= 0._RKC) then
                    err%occurred = .true._LK
                    err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                                &The input requested value ("//getStr(spec%proposalStdVec%val(idim))//SKC_") for the component "//getStr(idim)//&
                                SKC_" of the variable proposalStdVec for the proposal distribution of the sampler must be a positive real number."
                end if
            end do
        end block proposalStdVec_block

        ! Resolve conflicts.

        block
            ! proposalStart must be within the domain.
            integer(IK) :: idim
            do idim = 1, size(spec%proposalStart%val, 1, IK)
                if (spec%proposalStart%val(idim) < spec%domainCubeLimitLower%val(idim) .or. spec%domainCubeLimitUpper%val(idim) < spec%proposalStart%val(idim)) then
                    err%occurred = .true._LK
                    err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                                &The input requested value for the component "//getStr(idim)//SKC_" of the vector `proposalStart` ("//& ! LCOV_EXCL_LINE
                                getStr(spec%proposalStart%val(idim))//SKC_") must be within the range of the sampling domain defined in the program: ["//& ! LCOV_EXCL_LINE
                                getStr([spec%domainCubeLimitLower%val(idim), spec%domainCubeLimitUpper%val(idim)])//SKC_"]. &
                                &If you do not know an appropriate value for proposalStart, drop it from the input list. &
                                &The sampler will automatically assign an appropriate value to it."
                end if
            end do
        end block

        block
            ! `proposalStartDomainCubeLimitLower` must be larger than `domainCubeLimitLower` when start point is randomized.
            integer(IK) :: idim
            do idim = 1, size(spec%proposalStartDomainCubeLimitLower%val, kind = IK)
                ! check if the domain is set when random start point is requested.
                if (spec%proposalStartRandomized%val .and. spec%proposalStartDomainCubeLimitLower%val(idim) == spec%domainCubeLimitLower%def) then
                    err%occurred = .true._LK
                    err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. &
                                &You have requested a random start point by setting `proposalStartRandomized` to the logical/Boolean true value &
                                &while the element #"//getStr(idim)//SKC_" of `proposalStartDomainCubeLimitLower` has not been preset to a finite value. &
                                &This information is essential otherwise, how could the sampler draw points randomly from within an unspecified domain?"
                end if
                ! check if the random start point domain is within the boundaries of the domain of the target.
                if (spec%proposalStartDomainCubeLimitLower%val(idim) < spec%domainCubeLimitLower%val(idim)) then
                    err%occurred = .true._LK
                    err%msg =   err%msg//NL2//PROCEDURE_NAME//getFine(__FILE__, __LINE__)//SKC_": Error occurred. The component "//getStr(idim)//& ! LCOV_EXCL_LINE
                                SKC_" of the variable `proposalStartDomainCubeLimitLower` ("//getStr(spec%proposalStartDomainCubeLimitLower%val(idim))//& ! LCOV_EXCL_LINE
                                SKC_") cannot be smaller than the corresponding component of the variable `domainCubeLimitLower` ("//& ! LCOV_EXCL_LINE
                                getStr(spec%domainCubeLimitLower%val(idim))//"). If you do not know &an appropriate value to &
                                &set for `proposalStartDomainCubeLimitLower`, drop it from the input list. The sampler will &
                                &automatically assign an appropriate value to it."
                end if

            end do
        end block

    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  Return the refined Markov chain, given the input Markov chain and its specifications.
    !>  This procedure is a method of the [sfc_type](@ref sfc_type) class.
    !>
    !> \param[inout] sfc                            :   An object of class [sfc_type](@ref sfc_type).
    !> \param[in]    outputSampleSize               :   The requested refined sample size (**optional**). If the size of the refined sample is given as input,
    !!                                                  then the requested sample is directly generated based on the input size.
    !> \param[in]    outputSampleRefinementCount    :   The maximum number of times the sample can be refined (**optional**, default = `Infinity`).
    !!                                              :   For example, if set to 1, then only one round of refinement will be performed on the Markov chain.
    !> \param[in]    outputSampleRefinementMethod   :   The requested method of refining the sample (**optional**, default = "BatchMeans").
    function getErrRefinement(sfc, sampleWeight, sampleLogFunc, sampleState, outputSampleRefinementCount, outputSampleRefinementMethod) result(err)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getErrRefinement
#endif
        use, intrinsic :: iso_fortran_env, only: output_unit
        use pm_arrayRefine, only: setRefined
        use pm_arrayRebind, only: setRebound
        use pm_arrayResize, only: setResized
        use pm_sampleACT, only: getACT, batchMeans, batchMeansMax, cumSum, cumSumMax
        use pm_sampleQuan, only: getQuan, piwilin

        type(err_type)                                      :: err
        class(sfcmcmc_type) , intent(inout)                 :: sfc
        integer(IK)         , intent(in)    , contiguous    :: sampleWeight(:)
        real(RKC)           , intent(in)    , contiguous    :: sampleLogFunc(:), sampleState(:,:)
        integer(IK)         , intent(in)    , optional      :: outputSampleRefinementCount
        character(*, SK)    , intent(in)    , optional      :: outputSampleRefinementMethod
        real(RKC)           , allocatable                   :: logFuncStateAxis(:) ! ensures contiguity.
        character(*, SK)    , parameter                     :: PROCEDURE_NAME = SK_"@getErrRefinement()"
        integer(IK)                                         :: idim, ndim, nsam, ndimp1, nrefAlloc, outputSampleRefinementCount_def
        real(RKC)                                           :: act, ndimp1inv
        type(method_type)                                   :: method

        CHECK_ASSERTION(__LINE__, 0_IK < size(sampleState, 1, IK), SK_"The condition `0 < size(sampleState, 1)` must hold. shape(sampleState) = "//getStr(shape(sampleState, IK)))
        CHECK_ASSERTION(__LINE__, size(sampleWeight, 1, IK) == size(sampleState, 2, IK), SK_"The condition `size(sampleWeight) == size(sampleState, 2)` must hold. size(sampleWeight), shape(sampleState) = "//getStr([size(sampleWeight, 1, IK), shape(sampleState, IK)]))
        CHECK_ASSERTION(__LINE__, size(sampleWeight, 1, IK) == size(sampleLogFunc, 1, IK), SK_"The condition `size(sampleWeight) == size(sampleLogFunc)` must hold. size(sampleWeight), size(sampleLogFunc) = "//getStr([size(sampleWeight, 1, IK), size(sampleLogFunc, 1, IK)]))
        err = getErrRefinementMethod(method, outputSampleRefinementMethod)
        ndim = size(sampleState, 1, IK)
        nsam = size(sampleState, 2, IK)
        ndimp1 = ndim + 1_IK

        ! This is to avoid memory overflow due to extremely large `outputSampleRefinementCount` requested by the user.

        ndimp1inv = 1._RKC / ndimp1
        outputSampleRefinementCount_def = 20_IK  ! this is a temporary maximum value, to be increased later if needed
        if (present(outputSampleRefinementCount)) outputSampleRefinementCount_def = outputSampleRefinementCount
        nrefAlloc = min(2_IK, outputSampleRefinementCount_def)

        ! Allocate components.

        call setRebound(sfc%act, [0_IK, 0_IK], [ndim, nrefAlloc])
        call setRebound(sfc%nsam, 0_IK, nrefAlloc)
        call setRebound(sfc%sumw, 0_IK, nrefAlloc)

        sfc%nref = 0_IK
        sfc%nsam(sfc%nref) = size(sampleWeight, 1, IK)
        call setResized(logFuncStateAxis, sfc%nsam(sfc%nref))
        call setRebound(sfc%sampleLogFuncState, [0_IK, 1_IK], [ndim, nsam])
        sfc%sampleLogFuncState(1 : ndim, :) = sampleState
        sfc%sampleLogFuncState(0, :) = sampleLogFunc
        sfc%sampleWeight = sampleWeight

        ! Perform chain refinement.

        blockSufficientChainSize: if (ndim < nsam) then

            loopRefinement: do

                sfc%sumw(sfc%nref) = sum(sfc%sampleWeight(1 : sfc%nsam(sfc%nref)))

                ! Obtain the ACT for each individual axes.
                block
                    do idim = 0, ndim
                        logFuncStateAxis(1 : sfc%nsam(sfc%nref)) = sfc%sampleLogFuncState(idim, 1 : sfc%nsam(sfc%nref))
                        if (method%isViaCompactChain) then
                            if (method%isBatchMeansMax) then
                                sfc%act(idim, sfc%nref) = getACT(logFuncStateAxis(1 : sfc%nsam(sfc%nref)), batchMeansMax)
                            elseif (method%isBatchMeans) then
                                sfc%act(idim, sfc%nref) = getACT(logFuncStateAxis(1 : sfc%nsam(sfc%nref)), batchMeans)
                            elseif (method%isCumSumMax) then
                                sfc%act(idim, sfc%nref) = getACT(logFuncStateAxis(1 : sfc%nsam(sfc%nref)), cumSumMax)
                            elseif (method%isCumSum) then
                                sfc%act(idim, sfc%nref) = getACT(logFuncStateAxis(1 : sfc%nsam(sfc%nref)), cumSum)
                            else
                                error stop PROCEDURE_NAME//SK_"Internal Error: Unrecognized ACT computation method." ! LCOV_EXCL_LINE
                            end if
                        elseif (method%isViaVerboseChain) then
                            if (method%isBatchMeansMax) then
                                sfc%act(idim, sfc%nref) = getACT(logFuncStateAxis(1 : sfc%nsam(sfc%nref)), sfc%sampleWeight(1 : sfc%nsam(sfc%nref)), batchMeansMax)
                            elseif (method%isBatchMeans) then
                                sfc%act(idim, sfc%nref) = getACT(logFuncStateAxis(1 : sfc%nsam(sfc%nref)), sfc%sampleWeight(1 : sfc%nsam(sfc%nref)), batchMeans)
                            elseif (method%isCumSumMax) then
                                sfc%act(idim, sfc%nref) = getACT(logFuncStateAxis(1 : sfc%nsam(sfc%nref)), sfc%sampleWeight(1 : sfc%nsam(sfc%nref)), cumSumMax)
                            elseif (method%isCumSum) then
                                sfc%act(idim, sfc%nref) = getACT(logFuncStateAxis(1 : sfc%nsam(sfc%nref)), sfc%sampleWeight(1 : sfc%nsam(sfc%nref)), cumSum)
                            else
                                error stop PROCEDURE_NAME//SK_"Internal Error: Unrecognized ACT computation method." ! LCOV_EXCL_LINE
                            end if
                        else
                            error stop PROCEDURE_NAME//SK_": Internal error: Unrecognized chain format." ! LCOV_EXCL_LINE
                        endif
                    end do
                end block

                ! act is guaranteed to be larger than `1` by definition.

                if (method%isAvg) then
                    act = sum(sfc%act(0 : ndim, sfc%nref)) * ndimp1inv
                elseif (method%isMed) then
                    act = getQuan(piwilin, .5_RKC, sfc%act(0 : ndim, sfc%nref))
                elseif (method%isMax) then
                    act = maxval(sfc%act(0 : ndim, sfc%nref))
                elseif (method%isMin) then
                    act = minval(sfc%act(0 : ndim, sfc%nref))
                else
                    error stop PROCEDURE_NAME//SK_": Internal error: Unrecognized refinement method." ! LCOV_EXCL_LINE
                end if

                ! So far, we have computed the max ACT of the sample, but no refinement. Refine the sample only if needed.
                ! The output sample size can be less than the requested outputSampleSize if nsam < outputSampleSize.
                ! This is the right place to quit if needed, as all components for the current stage of the refinement are set.

                if (sfc%nref == outputSampleRefinementCount_def) exit loopRefinement
                if (act < 2._RKC) then
                    if (method%isViaCompactChain .and. method%isViaVerboseChain) then
                        method%isViaCompactChain = .false._LK
                        cycle loopRefinement
                    else
                        exit loopRefinement
                    end if
                end if

                ! Generate the refined sample, then put it back into sfc to start over again.

                sfc%nref = sfc%nref + 1_IK

                ! Reallocate to bigger array if needed.

                if (nrefAlloc < sfc%nref) then
                    nrefAlloc = min(nrefAlloc * 2, outputSampleRefinementCount_def)
                    call setRebound(sfc%act, [0_IK, 0_IK], [ndim, nrefAlloc])
                    call setRebound(sfc%nsam, 0_IK, nrefAlloc)
                    call setRebound(sfc%sumw, 0_IK, nrefAlloc)
                end if
                !block
                !   use pm_io
                !   call disp%show("act")
                !   call disp%show( act )
                !   print *, sfc%nref
                !   print *, sfc%nsam(0 : sfc%nref - 1)
                !   !call disp%show("sfc%nref")
                !   !call disp%show( sfc%nref )
                !   call disp%show("sfc%nsam(0 : sfc%nref - 1)")
                !   call disp%show( sfc%nsam(0 : sfc%nref - 1) )
                !end block
                if (act < 2._RKC) cycle loopRefinement ! no need for refinement. should happen only when transitioning from compact to verbose.
                call setRefined(sfc%sampleLogFuncState(:, 1 : sfc%nsam(sfc%nref - 1)), 2_IK, sfc%sampleWeight(1 : sfc%nsam(sfc%nref - 1)), skip = int(act, IK), rsize = sfc%nsam(sfc%nref))

            end do loopRefinement

        elseif (0_IK < nsam) then blockSufficientChainSize

            ! Return only the last sampled state if the input sample is too small.

            sfc%nsam(sfc%nref) = size(sampleWeight, 1, IK)
            call setResized(logFuncStateAxis, sfc%nsam(sfc%nref))
            call setRebound(sfc%sampleLogFuncState, [0_IK, 1_IK], [ndim, nsam])
            sfc%sampleLogFuncState(1 : ndim, :) = sampleState
            sfc%sampleLogFuncState(0, :) = sampleLogFunc
            sfc%sampleWeight = sampleWeight

        end if blockSufficientChainSize

        call setRebound(sfc%sampleLogFuncState, [0_IK, 1_IK], [ndim, sfc%nsam(sfc%nref)])
        call setResized(sfc%sampleWeight, sfc%nsam(sfc%nref))

    end function getErrRefinement

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getErrRefinementMethod(method, outputSampleRefinementMethod) result(err)

        use pm_val2int, only: getInt
        use pm_strASCII, only: getStrLower
        use pm_arrayReplace, only: getReplaced
        type(method_type), intent(inout) :: method
        character(*, SK), parameter :: PROCEDURE_NAME = SK_"@getErrRefinementMethod()"
        character(*, SK), intent(in), optional :: outputSampleRefinementMethod
        character(:, SK), allocatable :: outputSampleRefinementMethod_def
        type(err_type) :: err
        integer(IK) :: idim

        err%msg = ""
        err%occurred = .false._LK

        ! Define the ACT computation method.

        if (present(outputSampleRefinementMethod)) then
            outputSampleRefinementMethod_def = getStrLower(outputSampleRefinementMethod)
            if (0 < index(outputSampleRefinementMethod_def, SK_"batchmeansmax")) then
                method%isBatchMeansMax = .true._LK
            elseif (0 < index(outputSampleRefinementMethod_def, SK_"batchmeans")) then
                method%isBatchMeans = .true._LK
            elseif (0 < index(outputSampleRefinementMethod_def, SK_"cumsummax")) then
                method%isCumSumMax = .true._LK
            elseif (0 < index(outputSampleRefinementMethod_def, SK_"cumsum")) then
                method%isCumSum = .true._LK
            else
                err%occurred = .true._LK ! LCOV_EXCL_LINE
                err%msg = getFine(__FILE__, __LINE__)//PROCEDURE_NAME//SK_": Unknown unsupported ACT computation method name: "//outputSampleRefinementMethod_def ! LCOV_EXCL_LINE
                return ! LCOV_EXCL_LINE
            end if
            method%isViaCompactChain = 0 < index(outputSampleRefinementMethod_def, SK_"compact")
            method%isViaVerboseChain = 0 < index(outputSampleRefinementMethod_def, SK_"verbose")
        else
            outputSampleRefinementMethod_def = SK_"batchmeans"
            method%isBatchMeans = .true._LK
        end if

        ! Define the chain types to use for the ACT computation.

        if (.not. (method%isViaCompactChain .or. method%isViaVerboseChain)) then
            method%isViaCompactChain = .true._LK
            method%isViaVerboseChain = .true._LK
        end if

        ! Define the statistic to use for the collective ACT computation.

        outputSampleRefinementMethod_def = getReplaced(outputSampleRefinementMethod_def, SK_" ", SK_"-")
        method%isMed = index(outputSampleRefinementMethod_def, SK_"-med") > 0 .or. index(outputSampleRefinementMethod_def, SK_"-median") > 0
        method%isAvg = index(outputSampleRefinementMethod_def, SK_"-avg") > 0 .or. index(outputSampleRefinementMethod_def, SK_"-average") > 0
        method%isMin = index(outputSampleRefinementMethod_def, SK_"-min") > 0 .or. index(outputSampleRefinementMethod_def, SK_"-minimum") > 0
        method%isMax = index(outputSampleRefinementMethod_def, SK_"-max") > 0 .or. index(outputSampleRefinementMethod_def, SK_"-maximum") > 0
        idim = sum(getInt([method%isAvg, method%isMed, method%isMax, method%isMin]))
        err%occurred = 1_IK < idim
        if (err%occurred) then
            err%msg = getFine(__FILE__, __LINE__)//PROCEDURE_NAME//SK_": avg, med, min, max cannot be simultaneously specified in outputSampleRefinementMethod: "//outputSampleRefinementMethod_def ! LCOV_EXCL_LINE
            return ! LCOV_EXCL_LINE
        end if
        if (idim == 0_IK) method%isAvg = .true._LK ! default method of ACT summarization.

    end function getErrRefinementMethod

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
