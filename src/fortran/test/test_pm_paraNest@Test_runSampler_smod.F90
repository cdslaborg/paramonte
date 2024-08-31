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

!>  \brief These submodules contain tests of the modules [pm_sampling](@ref pm_sampling).
!>  \author Amir Shahmoradi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_paraNest) Test_runSampler_smod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test the ParaNest sampler with no input arguments or input file.
    module function test_runSampler_1() result(assertion)
        use pm_distMultiNormGrid, only: distMultiNormGrid_type
        use pm_distNormShellMix, only: DistMultiNormShellMix_type
        use pm_distMultiNormExpGammaMix, only: distMultiNormExpGammaMix_type
        implicit none
        logical(LK)                         :: assertion
        integer(IK), parameter              :: NDIM = 2_IK
        type(distMultiNormGrid_type)        :: DistMultiNormGrid
        type(DistMultiNormShellMix_type)    :: DistMultiNormShellMix
        type(distMultiNormExpGammaMix_type) :: DistMultiNormExpGammaMix
        type(paranest_type)                 :: self
        integer(IK)                         :: i

        assertion = .true._LK

!#if CODECOV_ENABLED || SAMPLER_TEST_ENABLED

        DistMultiNormGrid = distMultiNormGrid_type(NDIM, nmode = 5_IK, elongation = 1000._RK, domainMargin = 4.5_RK)
        DistMultiNormShellMix = DistMultiNormShellMix_type(NDIM)
        DistMultiNormExpGammaMix = distMultiNormExpGammaMix_type(NDIM)

        call self%runSampler( ndim = NDIM & ! LCOV_EXCL_LINE

                            ! MVN
                            !, getLogFunc = getLogFuncMVN & ! LCOV_EXCL_LINE
                            !, domainCubeLimitLower = [(-10._RK,i = 1, NDIM)] & ! LCOV_EXCL_LINE
                            !, domainCubeLimitUpper = [(+10._RK,i = 1, NDIM)] & ! LCOV_EXCL_LINE
                            !, domainBallAvg = [(0._RK,i = 1, NDIM)] & ! LCOV_EXCL_LINE
                            !, domainBallStd = [(5._RK,i = 1, NDIM)] & ! LCOV_EXCL_LINE
                            !, domain = "ball" & ! LCOV_EXCL_LINE

                            ! NormGrid
                            !, getLogFunc = getLogFuncDistMultiNormGrid & ! LCOV_EXCL_LINE
                            !, domainCubeLimitLower = DistMultiNormGrid%LowerLim & ! LCOV_EXCL_LINE
                            !, domainCubeLimitUpper = DistMultiNormGrid%UpperLim & ! LCOV_EXCL_LINE

                            ! NormShellMix
                            !, getLogFunc = getLogFuncDistMultiNormShellMix & ! LCOV_EXCL_LINE
                            !, domainCubeLimitLower = DistMultiNormShellMix%Domain%Cube%Limit%Lower & ! LCOV_EXCL_LINE
                            !, domainCubeLimitUpper = DistMultiNormShellMix%Domain%Cube%Limit%Upper & ! LCOV_EXCL_LINE
                            !, domainBallAvg = DistMultiNormShellMix%Domain%Ellipsoid%Center & ! LCOV_EXCL_LINE
                            !, domainBallCov = DistMultiNormShellMix%Domain%Ellipsoid%choLowCovUpp & ! LCOV_EXCL_LINE
                            !, domain = "ball" & ! LCOV_EXCL_LINE

                            ! NormLogGammaMix
                            !, getLogFunc = getLogFuncDistMultiNormExpGammaMix & ! LCOV_EXCL_LINE
                            !, domainCubeLimitLower = DistMultiNormExpGammaMix%LowerLim & ! LCOV_EXCL_LINE
                            !, domainCubeLimitUpper = DistMultiNormExpGammaMix%UpperLim & ! LCOV_EXCL_LINE

                            ! 2D Banana
                            , getLogFunc = getLogFuncBanana2D & ! LCOV_EXCL_LINE
                            , domainCubeLimitLower = [(-7._RK,i = 1, NDIM)] & ! LCOV_EXCL_LINE
                            , domainCubeLimitUpper = [(+7._RK,i = 1, NDIM)] & ! LCOV_EXCL_LINE

                            ! 2D EggBox
                            !, getLogFunc = getLogFuncEggBox2D & ! LCOV_EXCL_LINE
                            !, domainCubeLimitLower = [(0._RK,i = 1, NDIM)] & ! LCOV_EXCL_LINE
                            !, domainCubeLimitUpper = [(1._RK,i = 1, NDIM)] & ! LCOV_EXCL_LINE

                            , parallelismMpiFinalizeEnabled = .false._LK & ! LCOV_EXCL_LINE
                            , outputFileName = test%dir%out//"/"//MODULE_NAME//SK_"/test_runSampler_1" & ! LCOV_EXCL_LINE
                            , outputRestartFileFormat = "ascii" & ! LCOV_EXCL_LINE
                            , outputReportPeriod = 2000_IK & ! LCOV_EXCL_LINE
                            !, domain = "cube" & ! LCOV_EXCL_LINE
                            !, domainErrCount = huge(1_IK) & ! LCOV_EXCL_LINE
                            !, domainErrCountMax = huge(1_IK) & ! LCOV_EXCL_LINE
                            !, domainPartitionAdaptationCount = 0_IK & ! LCOV_EXCL_LINE
                            , domainPartitionAdaptationPeriod = 1000_IK & ! LCOV_EXCL_LINE
                            , liveSampleSize = 2000_IK & ! LCOV_EXCL_LINE
                            , proposal = "partition" & ! LCOV_EXCL_LINE
                            , domainPartitionObject = "ball" & ! LCOV_EXCL_LINE
                            !, domainPartitionMethod = "paramonte-maxden" & ! LCOV_EXCL_LINE
                            !, domainPartitionMethod = "dynesty" & ! LCOV_EXCL_LINE
                            !, domainPartitionMethod = "multinest" & ! LCOV_EXCL_LINE
                            , randomSeed = 1357_IK & ! LCOV_EXCL_LINE
                            !, domainPartitionKvolumeNumRecursionMax = 0_IK & ! LCOV_EXCL_LINE
                            !, domainPartitionBiasCorrectionEnabled = .false._LK  & ! LCOV_EXCL_LINE
                            !, domainPartitionOptimizationShapeEnabled = .false._LK  & ! LCOV_EXCL_LINE
                            !, domainPartitionOptimizationScaleEnabled = .false._LK  & ! LCOV_EXCL_LINE
                            !, domainPartitionOptimizationShapeScaleEnabled = .false._LK  & ! LCOV_EXCL_LINE
                            !, domainPartitionFactorExpansion = 77._RK & ! LCOV_EXCL_LINE
                            !, domainPartitionCountMax = 1_IK & ! LCOV_EXCL_LINE
                            )
        assertion = assertion .and. .not. self%err%occurred
!#endif

    contains

#if     CFI_ENABLED
        function getLogFuncDistMultiNormGrid(ndim,Point) result(logFunc) bind(C)
            integer(IK) , intent(in), value :: ndim
#else
        function getLogFuncDistMultiNormGrid(ndim,Point) result(logFunc)
            integer(IK) , intent(in)        :: ndim
#endif
            real(RK)    , intent(in)        :: Point(ndim)
            real(RK)                        :: logFunc
            logFunc = DistMultiNormGrid%logpdf(Point)
        end function getLogFuncDistMultiNormGrid

#if     CFI_ENABLED
        function getLogFuncDistMultiNormShellMix(ndim,Point) result(logFunc) bind(C)
            integer(IK) , intent(in), value :: ndim
#else
        function getLogFuncDistMultiNormShellMix(ndim,Point) result(logFunc)
            integer(IK) , intent(in)        :: ndim
#endif
            real(RK)    , intent(in)        :: Point(ndim)
            real(RK)                        :: logFunc
            logFunc = DistMultiNormShellMix%logpdf(Point)
        end function getLogFuncDistMultiNormShellMix

#if     CFI_ENABLED
        function getLogFuncDistMultiNormExpGammaMix(ndim,Point) result(logFunc) bind(C)
            integer(IK) , intent(in), value :: ndim
#else
        function getLogFuncDistMultiNormExpGammaMix(ndim,Point) result(logFunc)
            integer(IK) , intent(in)        :: ndim
#endif
            real(RK)    , intent(in)        :: Point(ndim)
            real(RK)                        :: logFunc
            logFunc = DistMultiNormExpGammaMix%logpdf(Point)
        end function getLogFuncDistMultiNormExpGammaMix

    end function test_runSampler_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule Test_runSampler_smod
