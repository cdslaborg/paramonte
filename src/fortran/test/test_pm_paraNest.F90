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

!>  \brief This module contains tests of the module [pm_sampling](@ref pm_sampling).
!>  \author Amir Shahmoradi

module test_pm_paraNest

    use pm_kind, only: IK, RK
    use pm_test, only: test_type, LK, getLogFuncMVN, getLogFuncBanana2D, getLogFuncEggBox2D
    use pm_sampling

    implicit none

    private
    public :: setTest
    type(test_type) :: test

    interface
        module function test_runSampler_1 () result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_2 () result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_3 () result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_4 () result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_5 () result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_6 () result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_7 () result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_8 () result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_9 () result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_10() result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_11() result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_12() result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_13() result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_14() result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_15() result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_16() result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_17() result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_18() result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_19() result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_20() result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_21() result(assertion); logical(LK) :: assertion; end
        !module function test_runSampler_22() result(assertion); logical(LK) :: assertion; end
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none

        test = test_type(MODULE_NAME)

        call test%run(test_runSampler_1, SK_"test_runSampler_1")
        !call test%run(test_runSampler_2 , SK_"test_runSampler_2")
        !call test%run(test_runSampler_3 , SK_"test_runSampler_3")
        !call test%run(test_runSampler_4 , SK_"test_runSampler_4")
        !call test%run(test_runSampler_5 , SK_"test_runSampler_5")
        !call test%run(test_runSampler_6 , SK_"test_runSampler_6")
        !call test%run(test_runSampler_7 , SK_"test_runSampler_7")
        !call test%run(test_runSampler_8 , SK_"test_runSampler_8")
        !call test%run(test_runSampler_9 , SK_"test_runSampler_9")
        !call test%run(test_runSampler_10, SK_"test_runSampler_10")
        !call test%run(test_runSampler_11, SK_"test_runSampler_11")
        !call test%run(test_runSampler_12, SK_"test_runSampler_12")
        !call test%run(test_runSampler_13, SK_"test_runSampler_13")
        !call test%run(test_runSampler_14, SK_"test_runSampler_14")
        !call test%run(test_runSampler_15, SK_"test_runSampler_15")
        !call test%run(test_runSampler_16, SK_"test_runSampler_16")
        !call test%run(test_runSampler_17, SK_"test_runSampler_17")
        !call test%run(test_runSampler_18, SK_"test_runSampler_18")
        !call test%run(test_runSampler_19, SK_"test_runSampler_19")
        !call test%run(test_runSampler_20, SK_"test_runSampler_20")
        !call test%run(test_runSampler_21, SK_"test_runSampler_21")
        !call test%run(test_runSampler_22, SK_"test_runSampler_22")

        !call test%run(test_specbase_randomSeed_type_1, SK_"test_specbase_randomSeed_type_1")
        !call test%run(test_specbase_randomSeed_type_2, SK_"test_specbase_randomSeed_type_2")
        !call test%run(test_specbase_outputSampleSize_type_1, SK_"test_specbase_outputSampleSize_type_1")
        !call test%run(test_specbase_outputSampleSize_type_2, SK_"test_specbase_outputSampleSize_type_2")
        !call test%run(test_specbase_outputSampleSize_type_3, SK_"test_specbase_outputSampleSize_type_3")
        !call test%run(test_specbase_outputSampleSize_type_4, SK_"test_specbase_outputSampleSize_type_4")
        !call test%run(test_specbase_outputSeparator_type_1, SK_"test_specbase_outputSeparator_type_1")
        !call test%run(test_specbase_outputSeparator_type_2, SK_"test_specbase_outputSeparator_type_2")
        !call test%run(test_specbase_outputSeparator_type_3, SK_"test_specbase_outputSeparator_type_3")
        !call test%run(test_specbase_outputSeparator_type_4, SK_"test_specbase_outputSeparator_type_4")
        !call test%run(test_specbase_outputSeparator_type_5, SK_"test_specbase_outputSeparator_type_5")
        !call test%run(test_specbase_outputSeparator_type_6, SK_"test_specbase_outputSeparator_type_6")
        !call test%run(test_specbase_outputChainFileFormat_type_1, SK_"test_specbase_outputChainFileFormat_type_1")
        !call test%run(test_specbase_outputChainFileFormat_type_2, SK_"test_specbase_outputChainFileFormat_type_2")
        !call test%run(test_specbase_outputChainFileFormat_type_3, SK_"test_specbase_outputChainFileFormat_type_3")
        !call test%run(test_specbase_outputChainFileFormat_type_4, SK_"test_specbase_outputChainFileFormat_type_4")
        !call test%run(test_specbase_outputChainFileFormat_type_5, SK_"test_specbase_outputChainFileFormat_type_5")
        !call test%run(test_specbase_outputColumnWidth_type_1, SK_"test_specbase_outputColumnWidth_type_1")
        !call test%run(test_specbase_outputColumnWidth_type_2, SK_"test_specbase_outputColumnWidth_type_2")
        !call test%run(test_specbase_outputColumnWidth_type_3, SK_"test_specbase_outputColumnWidth_type_3")
        !call test%run(test_specbase_outputColumnWidth_type_4, SK_"test_specbase_outputColumnWidth_type_4")
        !call test%run(test_specbase_outputRestartFileFormat_type_1, SK_"test_specbase_outputRestartFileFormat_type_1")
        !call test%run(test_specbase_outputRestartFileFormat_type_2, SK_"test_specbase_outputRestartFileFormat_type_2")
        !call test%run(test_specbase_outputRestartFileFormat_type_3, SK_"test_specbase_outputRestartFileFormat_type_3")
        !call test%run(test_specbase_outputStatus_type_1, SK_"test_specbase_outputStatus_type_1")
        !call test%run(test_specbase_outputStatus_type_2, SK_"test_specbase_outputStatus_type_2")
        !call test%run(test_specbase_domainCubeLimitLower_type_1, SK_"test_specbase_domainCubeLimitLower_type_1")
        !call test%run(test_specbase_domainCubeLimitLower_type_2, SK_"test_specbase_domainCubeLimitLower_type_2")
        !call test%run(test_specbase_domainCubeLimitUpper_type_1, SK_"test_specbase_domainCubeLimitUpper_type_1")
        !call test%run(test_specbase_domainCubeLimitUpper_type_2, SK_"test_specbase_domainCubeLimitUpper_type_2")
        !call test%run(test_specbase_domainCubeLimitUpper_type_3, SK_"test_specbase_domainCubeLimitUpper_type_3")
        !call test%run(test_specbase_domainCubeLimitUpper_type_4, SK_"test_specbase_domainCubeLimitUpper_type_4")
        !call test%run(test_specbase_outputPrecision_type_1, SK_"test_specbase_outputPrecision_type_1")
        !call test%run(test_specbase_outputPrecision_type_2, SK_"test_specbase_outputPrecision_type_2")
        !call test%run(test_specbase_outputReportPeriod_type_1, SK_"test_specbase_outputReportPeriod_type_1")
        !call test%run(test_specbase_outputReportPeriod_type_2, SK_"test_specbase_outputReportPeriod_type_2")
        !call test%run(test_specbase_parallelism_type_1, SK_"test_specbase_parallelism_type_1")
        !call test%run(test_specbase_parallelism_type_2, SK_"test_specbase_parallelism_type_2")
        !call test%run(test_specbase_parallelism_type_3, SK_"test_specbase_parallelism_type_3")
        !call test%run(test_specbase_targetAcceptanceRate_type_1, SK_"test_specbase_targetAcceptanceRate_type_1")
        !call test%run(test_specbase_targetAcceptanceRate_type_2, SK_"test_specbase_targetAcceptanceRate_type_2")
        !call test%run(test_specbase_targetAcceptanceRate_type_3, SK_"test_specbase_targetAcceptanceRate_type_3")
        !call test%run(test_specbase_targetAcceptanceRate_type_4, SK_"test_specbase_targetAcceptanceRate_type_4")
        !call test%run(test_specbase_targetAcceptanceRate_type_5, SK_"test_specbase_targetAcceptanceRate_type_5")
        !call test%run(test_specbase_targetAcceptanceRate_type_6, SK_"test_specbase_targetAcceptanceRate_type_6")
        !call test%run(test_specbase_domainErrCount_type_1, SK_"test_specbase_domainErrCount_type_1")
        !call test%run(test_specbase_domainErrCount_type_2, SK_"test_specbase_domainErrCount_type_2")
        !call test%run(test_specbase_domainErrCount_type_3, SK_"test_specbase_domainErrCount_type_3")
        !call test%run(test_specbase_domainErrCountMax_type_1, SK_"test_specbase_domainErrCountMax_type_1")
        !call test%run(test_specbase_domainErrCountMax_type_2, SK_"test_specbase_domainErrCountMax_type_2")
        !call test%run(test_specbase_domainErrCountMax_type_3, SK_"test_specbase_domainErrCountMax_type_3")

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_paraNest
