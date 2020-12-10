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
!!!!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This module contains tests of the module [ParaDRAM_mod](@ref paradram_mod).
!>  \author Amir Shahmoradi

! module Test_ParaDXXX_mod

    use Constants_mod, only: IK, RK
    use Test_mod, only: Test_type, getLogFuncMVN
    use ParaDXXX_mod

    !use Statistics_mod, only: paradramPrintEnabled
    implicit none
    !paradramPrintEnabled = .true.

    private
    public :: test_ParaDXXX

    type(Test_type) :: Test

    interface
        module function test_runSampler_1 () result(assertion); logical :: assertion; end
        module function test_runSampler_2 () result(assertion); logical :: assertion; end
        module function test_runSampler_3 () result(assertion); logical :: assertion; end
        module function test_runSampler_4 () result(assertion); logical :: assertion; end
        module function test_runSampler_5 () result(assertion); logical :: assertion; end
        module function test_runSampler_6 () result(assertion); logical :: assertion; end
        module function test_runSampler_7 () result(assertion); logical :: assertion; end
        module function test_runSampler_8 () result(assertion); logical :: assertion; end
        module function test_runSampler_9 () result(assertion); logical :: assertion; end
        module function test_runSampler_10() result(assertion); logical :: assertion; end
        module function test_runSampler_11() result(assertion); logical :: assertion; end
        module function test_runSampler_12() result(assertion); logical :: assertion; end
        module function test_runSampler_13() result(assertion); logical :: assertion; end
        module function test_runSampler_14() result(assertion); logical :: assertion; end
        module function test_runSampler_15() result(assertion); logical :: assertion; end
        module function test_runSampler_16() result(assertion); logical :: assertion; end
        module function test_runSampler_17() result(assertion); logical :: assertion; end
        module function test_runSampler_18() result(assertion); logical :: assertion; end
        module function test_runSampler_19() result(assertion); logical :: assertion; end
    end interface

    interface
        module function test_SpecBase_RandomSeed_type_1             () result(assertion); logical :: assertion; end
        module function test_SpecBase_RandomSeed_type_2             () result(assertion); logical :: assertion; end
        module function test_SpecBase_SampleSize_type_1             () result(assertion); logical :: assertion; end
        module function test_SpecBase_SampleSize_type_2             () result(assertion); logical :: assertion; end
        module function test_SpecBase_SampleSize_type_3             () result(assertion); logical :: assertion; end
        module function test_SpecBase_SampleSize_type_4             () result(assertion); logical :: assertion; end
        module function test_SpecBase_OutputDelimiter_type_1        () result(assertion); logical :: assertion; end
        module function test_SpecBase_OutputDelimiter_type_2        () result(assertion); logical :: assertion; end
        module function test_SpecBase_OutputDelimiter_type_3        () result(assertion); logical :: assertion; end
        module function test_SpecBase_OutputDelimiter_type_4        () result(assertion); logical :: assertion; end
        module function test_SpecBase_OutputDelimiter_type_5        () result(assertion); logical :: assertion; end
        module function test_SpecBase_OutputDelimiter_type_6        () result(assertion); logical :: assertion; end
        module function test_SpecBase_ChainFileFormat_type_1        () result(assertion); logical :: assertion; end
        module function test_SpecBase_ChainFileFormat_type_2        () result(assertion); logical :: assertion; end
        module function test_SpecBase_ChainFileFormat_type_3        () result(assertion); logical :: assertion; end
        module function test_SpecBase_ChainFileFormat_type_4        () result(assertion); logical :: assertion; end
        module function test_SpecBase_ChainFileFormat_type_5        () result(assertion); logical :: assertion; end
        module function test_SpecBase_OutputColumnWidth_type_1      () result(assertion); logical :: assertion; end
        module function test_SpecBase_OutputColumnWidth_type_2      () result(assertion); logical :: assertion; end
        module function test_SpecBase_OutputColumnWidth_type_3      () result(assertion); logical :: assertion; end
        module function test_SpecBase_OutputColumnWidth_type_4      () result(assertion); logical :: assertion; end
        module function test_SpecBase_RestartFileFormat_type_1      () result(assertion); logical :: assertion; end
        module function test_SpecBase_RestartFileFormat_type_2      () result(assertion); logical :: assertion; end
        module function test_SpecBase_RestartFileFormat_type_3      () result(assertion); logical :: assertion; end
        module function test_SpecBase_OverwriteRequested_type_1     () result(assertion); logical :: assertion; end
        module function test_SpecBase_OverwriteRequested_type_2     () result(assertion); logical :: assertion; end
        module function test_SpecBase_DomainLowerLimitVec_type_1    () result(assertion); logical :: assertion; end
        module function test_SpecBase_DomainLowerLimitVec_type_2    () result(assertion); logical :: assertion; end
        module function test_SpecBase_DomainUpperLimitVec_type_1    () result(assertion); logical :: assertion; end
        module function test_SpecBase_DomainUpperLimitVec_type_2    () result(assertion); logical :: assertion; end
        module function test_SpecBase_DomainUpperLimitVec_type_3    () result(assertion); logical :: assertion; end
        module function test_SpecBase_DomainUpperLimitVec_type_4    () result(assertion); logical :: assertion; end
        module function test_SpecBase_OutputRealPrecision_type_1    () result(assertion); logical :: assertion; end
        module function test_SpecBase_OutputRealPrecision_type_2    () result(assertion); logical :: assertion; end
        module function test_SpecBase_ProgressReportPeriod_type_1   () result(assertion); logical :: assertion; end
        module function test_SpecBase_ProgressReportPeriod_type_2   () result(assertion); logical :: assertion; end
        module function test_SpecBase_ParallelizationModel_type_1   () result(assertion); logical :: assertion; end
        module function test_SpecBase_ParallelizationModel_type_2   () result(assertion); logical :: assertion; end
        module function test_SpecBase_ParallelizationModel_type_3   () result(assertion); logical :: assertion; end
        module function test_SpecBase_TargetAcceptanceRate_type_1   () result(assertion); logical :: assertion; end
        module function test_SpecBase_TargetAcceptanceRate_type_2   () result(assertion); logical :: assertion; end
        module function test_SpecBase_TargetAcceptanceRate_type_3   () result(assertion); logical :: assertion; end
        module function test_SpecBase_TargetAcceptanceRate_type_4   () result(assertion); logical :: assertion; end
        module function test_SpecBase_TargetAcceptanceRate_type_5   () result(assertion); logical :: assertion; end
        module function test_SpecBase_TargetAcceptanceRate_type_6   () result(assertion); logical :: assertion; end
        module function test_SpecBase_MaxNumDomainCheckToWarn_type_1() result(assertion); logical :: assertion; end
        module function test_SpecBase_MaxNumDomainCheckToWarn_type_2() result(assertion); logical :: assertion; end
        module function test_SpecBase_MaxNumDomainCheckToWarn_type_3() result(assertion); logical :: assertion; end
        module function test_SpecBase_MaxNumDomainCheckToStop_type_1() result(assertion); logical :: assertion; end
        module function test_SpecBase_MaxNumDomainCheckToStop_type_2() result(assertion); logical :: assertion; end
        module function test_SpecBase_MaxNumDomainCheckToStop_type_3() result(assertion); logical :: assertion; end
    end interface

    interface
        module function test_SpecMCMC_ChainSize_type_1                  () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ChainSize_type_2                  () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ScaleFactor_type_1                () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ScaleFactor_type_2                () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ScaleFactor_type_3                () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ScaleFactor_type_4                () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ScaleFactor_type_5                () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ScaleFactor_type_6                () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ScaleFactor_type_7                () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ScaleFactor_type_8                () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ProposalModel_type_1              () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ProposalModel_type_2              () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ProposalModel_type_3              () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ProposalModel_type_4              () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_StartPointVec_type_1              () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_StartPointVec_type_2              () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_StartPointVec_type_3              () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_StartPointVec_type_4              () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_StartPointVec_type_5              () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ProposalStartCovMat_type_1        () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ProposalStartCovMat_type_2        () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ProposalStartCovMat_type_3        () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ProposalStartCovMat_type_4        () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ProposalStartCovMat_type_5        () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ProposalStartCorMat_type_1        () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ProposalStartCorMat_type_2        () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ProposalStartCorMat_type_3        () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ProposalStartCorMat_type_4        () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ProposalStartCorMat_type_5        () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ProposalStartStdVec_type_1        () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ProposalStartStdVec_type_2        () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ProposalStartStdVec_type_3        () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ProposalStartStdVec_type_4        () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_ProposalStartStdVec_type_5        () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_SampleRefinementCount_type_1      () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_SampleRefinementCount_type_2      () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_SampleRefinementCount_type_3      () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_SampleRefinementMethod_type_1     () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_SampleRefinementMethod_type_2     () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_SampleRefinementMethod_type_3     () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_SampleRefinementMethod_type_4     () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_SampleRefinementMethod_type_5     () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_SampleRefinementMethod_type_6     () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_SampleRefinementMethod_type_7     () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_SampleRefinementMethod_type_8     () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_SampleRefinementMethod_type_9     () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_SampleRefinementMethod_type_10    () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_SampleRefinementMethod_type_11    () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_SampleRefinementMethod_type_12    () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_SampleRefinementMethod_type_13    () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_RandomStartPointRequested_type_1  () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_RandomStartPointRequested_type_2  () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_RandomStartPointRequested_type_3  () result(assertion); logical :: assertion; end
        module function test_SpecMCMC_RandomStartPointRequested_type_4  () result(assertion); logical :: assertion; end
        module function test_RSPDLowerLimitVec_type_1                   () result(assertion); logical :: assertion; end
        module function test_RSPDLowerLimitVec_type_2                   () result(assertion); logical :: assertion; end
        module function test_RSPDLowerLimitVec_type_3                   () result(assertion); logical :: assertion; end
        module function test_RSPDUpperLimitVec_type_1                   () result(assertion); logical :: assertion; end
        module function test_RSPDUpperLimitVec_type_2                   () result(assertion); logical :: assertion; end
        module function test_RSPDUpperLimitVec_type_3                   () result(assertion); logical :: assertion; end
        module function test_RSPDUpperLimitVec_type_4                   () result(assertion); logical :: assertion; end
    end interface

    interface
        module function test_SpecDRAM_AdaptiveUpdateCount_type_1            () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_AdaptiveUpdateCount_type_2            () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_AdaptiveUpdateCount_type_3            () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_AdaptiveUpdatePeriod_type_1           () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_AdaptiveUpdatePeriod_type_2           () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_AdaptiveUpdatePeriod_type_3           () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_DelayedRejectionCount_type_1          () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_DelayedRejectionCount_type_2          () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_DelayedRejectionCount_type_3          () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_DelayedRejectionCount_type_4          () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_GreedyAdaptationCount_type_1          () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_GreedyAdaptationCount_type_2          () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_GreedyAdaptationCount_type_3          () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_BurninAdaptationMeasure_type_1        () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_BurninAdaptationMeasure_type_2        () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_BurninAdaptationMeasure_type_3        () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_BurninAdaptationMeasure_type_4        () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_DelayedRejectionScaleFactorVec_type_1 () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_DelayedRejectionScaleFactorVec_type_2 () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_DelayedRejectionScaleFactorVec_type_3 () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_DelayedRejectionScaleFactorVec_type_4 () result(assertion); logical :: assertion; end
        module function test_SpecDRAM_DelayedRejectionScaleFactorVec_type_5 () result(assertion); logical :: assertion; end
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_ParaDXXX()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)

        call Test%run(test_runSampler_1 , "test_runSampler_1")
        call Test%run(test_runSampler_2 , "test_runSampler_2")
        call Test%run(test_runSampler_3 , "test_runSampler_3")
        call Test%run(test_runSampler_4 , "test_runSampler_4")
        call Test%run(test_runSampler_5 , "test_runSampler_5")
        call Test%run(test_runSampler_6 , "test_runSampler_6")
        call Test%run(test_runSampler_7 , "test_runSampler_7")
        call Test%run(test_runSampler_8 , "test_runSampler_8")
        call Test%run(test_runSampler_9 , "test_runSampler_9")
        call Test%run(test_runSampler_10, "test_runSampler_10")
        call Test%run(test_runSampler_11, "test_runSampler_11")
        call Test%run(test_runSampler_12, "test_runSampler_12")
        call Test%run(test_runSampler_13, "test_runSampler_13")
        call Test%run(test_runSampler_14, "test_runSampler_14")
        call Test%run(test_runSampler_15, "test_runSampler_15")
        call Test%run(test_runSampler_16, "test_runSampler_16")
        call Test%run(test_runSampler_17, "test_runSampler_17")
        call Test%run(test_runSampler_18, "test_runSampler_18")
        call Test%run(test_runSampler_19, "test_runSampler_19")

        call Test%run(test_SpecBase_RandomSeed_type_1, "test_SpecBase_RandomSeed_type_1")
        call Test%run(test_SpecBase_RandomSeed_type_2, "test_SpecBase_RandomSeed_type_2")
        call Test%run(test_SpecBase_SampleSize_type_1, "test_SpecBase_SampleSize_type_1")
        call Test%run(test_SpecBase_SampleSize_type_2, "test_SpecBase_SampleSize_type_2")
        call Test%run(test_SpecBase_SampleSize_type_3, "test_SpecBase_SampleSize_type_3")
        call Test%run(test_SpecBase_SampleSize_type_4, "test_SpecBase_SampleSize_type_4")
        call Test%run(test_SpecBase_OutputDelimiter_type_1, "test_SpecBase_OutputDelimiter_type_1")
        call Test%run(test_SpecBase_OutputDelimiter_type_2, "test_SpecBase_OutputDelimiter_type_2")
        call Test%run(test_SpecBase_OutputDelimiter_type_3, "test_SpecBase_OutputDelimiter_type_3")
        call Test%run(test_SpecBase_OutputDelimiter_type_4, "test_SpecBase_OutputDelimiter_type_4")
        call Test%run(test_SpecBase_OutputDelimiter_type_5, "test_SpecBase_OutputDelimiter_type_5")
        call Test%run(test_SpecBase_OutputDelimiter_type_6, "test_SpecBase_OutputDelimiter_type_6")
        call Test%run(test_SpecBase_ChainFileFormat_type_1, "test_SpecBase_ChainFileFormat_type_1")
        call Test%run(test_SpecBase_ChainFileFormat_type_2, "test_SpecBase_ChainFileFormat_type_2")
        call Test%run(test_SpecBase_ChainFileFormat_type_3, "test_SpecBase_ChainFileFormat_type_3")
        call Test%run(test_SpecBase_ChainFileFormat_type_4, "test_SpecBase_ChainFileFormat_type_4")
        call Test%run(test_SpecBase_ChainFileFormat_type_5, "test_SpecBase_ChainFileFormat_type_5")
        call Test%run(test_SpecBase_OutputColumnWidth_type_1, "test_SpecBase_OutputColumnWidth_type_1")
        call Test%run(test_SpecBase_OutputColumnWidth_type_2, "test_SpecBase_OutputColumnWidth_type_2")
        call Test%run(test_SpecBase_OutputColumnWidth_type_3, "test_SpecBase_OutputColumnWidth_type_3")
        call Test%run(test_SpecBase_OutputColumnWidth_type_4, "test_SpecBase_OutputColumnWidth_type_4")
        call Test%run(test_SpecBase_RestartFileFormat_type_1, "test_SpecBase_RestartFileFormat_type_1")
        call Test%run(test_SpecBase_RestartFileFormat_type_2, "test_SpecBase_RestartFileFormat_type_2")
        call Test%run(test_SpecBase_RestartFileFormat_type_3, "test_SpecBase_RestartFileFormat_type_3")
        call Test%run(test_SpecBase_OverwriteRequested_type_1, "test_SpecBase_OverwriteRequested_type_1")
        call Test%run(test_SpecBase_OverwriteRequested_type_2, "test_SpecBase_OverwriteRequested_type_2")
        call Test%run(test_SpecBase_DomainLowerLimitVec_type_1, "test_SpecBase_DomainLowerLimitVec_type_1")
        call Test%run(test_SpecBase_DomainLowerLimitVec_type_2, "test_SpecBase_DomainLowerLimitVec_type_2")
        call Test%run(test_SpecBase_DomainUpperLimitVec_type_1, "test_SpecBase_DomainUpperLimitVec_type_1")
        call Test%run(test_SpecBase_DomainUpperLimitVec_type_2, "test_SpecBase_DomainUpperLimitVec_type_2")
        call Test%run(test_SpecBase_DomainUpperLimitVec_type_3, "test_SpecBase_DomainUpperLimitVec_type_3")
        call Test%run(test_SpecBase_DomainUpperLimitVec_type_4, "test_SpecBase_DomainUpperLimitVec_type_4")
        call Test%run(test_SpecBase_OutputRealPrecision_type_1, "test_SpecBase_OutputRealPrecision_type_1")
        call Test%run(test_SpecBase_OutputRealPrecision_type_2, "test_SpecBase_OutputRealPrecision_type_2")
        call Test%run(test_SpecBase_ProgressReportPeriod_type_1, "test_SpecBase_ProgressReportPeriod_type_1")
        call Test%run(test_SpecBase_ProgressReportPeriod_type_2, "test_SpecBase_ProgressReportPeriod_type_2")
        call Test%run(test_SpecBase_ParallelizationModel_type_1, "test_SpecBase_ParallelizationModel_type_1")
        call Test%run(test_SpecBase_ParallelizationModel_type_2, "test_SpecBase_ParallelizationModel_type_2")
        call Test%run(test_SpecBase_ParallelizationModel_type_3, "test_SpecBase_ParallelizationModel_type_3")
        call Test%run(test_SpecBase_TargetAcceptanceRate_type_1, "test_SpecBase_TargetAcceptanceRate_type_1")
        call Test%run(test_SpecBase_TargetAcceptanceRate_type_2, "test_SpecBase_TargetAcceptanceRate_type_2")
        call Test%run(test_SpecBase_TargetAcceptanceRate_type_3, "test_SpecBase_TargetAcceptanceRate_type_3")
        call Test%run(test_SpecBase_TargetAcceptanceRate_type_4, "test_SpecBase_TargetAcceptanceRate_type_4")
        call Test%run(test_SpecBase_TargetAcceptanceRate_type_5, "test_SpecBase_TargetAcceptanceRate_type_5")
        call Test%run(test_SpecBase_TargetAcceptanceRate_type_6, "test_SpecBase_TargetAcceptanceRate_type_6")
        call Test%run(test_SpecBase_MaxNumDomainCheckToWarn_type_1, "test_SpecBase_MaxNumDomainCheckToWarn_type_1")
        call Test%run(test_SpecBase_MaxNumDomainCheckToWarn_type_2, "test_SpecBase_MaxNumDomainCheckToWarn_type_2")
        call Test%run(test_SpecBase_MaxNumDomainCheckToWarn_type_3, "test_SpecBase_MaxNumDomainCheckToWarn_type_3")
        call Test%run(test_SpecBase_MaxNumDomainCheckToStop_type_1, "test_SpecBase_MaxNumDomainCheckToStop_type_1")
        call Test%run(test_SpecBase_MaxNumDomainCheckToStop_type_2, "test_SpecBase_MaxNumDomainCheckToStop_type_2")
        call Test%run(test_SpecBase_MaxNumDomainCheckToStop_type_3, "test_SpecBase_MaxNumDomainCheckToStop_type_3")

        call Test%run(test_SpecMCMC_ChainSize_type_1, "test_SpecMCMC_ChainSize_type_1")
        call Test%run(test_SpecMCMC_ChainSize_type_2, "test_SpecMCMC_ChainSize_type_2")
        call Test%run(test_SpecMCMC_ScaleFactor_type_1, "test_SpecMCMC_ScaleFactor_type_1")
        call Test%run(test_SpecMCMC_ScaleFactor_type_2, "test_SpecMCMC_ScaleFactor_type_2")
        call Test%run(test_SpecMCMC_ScaleFactor_type_3, "test_SpecMCMC_ScaleFactor_type_3")
        call Test%run(test_SpecMCMC_ScaleFactor_type_4, "test_SpecMCMC_ScaleFactor_type_4")
        call Test%run(test_SpecMCMC_ScaleFactor_type_5, "test_SpecMCMC_ScaleFactor_type_5")
        call Test%run(test_SpecMCMC_ScaleFactor_type_6, "test_SpecMCMC_ScaleFactor_type_6")
        call Test%run(test_SpecMCMC_ScaleFactor_type_7, "test_SpecMCMC_ScaleFactor_type_7")
        call Test%run(test_SpecMCMC_ScaleFactor_type_8, "test_SpecMCMC_ScaleFactor_type_8")
        call Test%run(test_SpecMCMC_ProposalModel_type_1, "test_SpecMCMC_ProposalModel_type_1")
        call Test%run(test_SpecMCMC_ProposalModel_type_2, "test_SpecMCMC_ProposalModel_type_2")
        call Test%run(test_SpecMCMC_ProposalModel_type_3, "test_SpecMCMC_ProposalModel_type_3")
        call Test%run(test_SpecMCMC_ProposalModel_type_4, "test_SpecMCMC_ProposalModel_type_4")
        call Test%run(test_SpecMCMC_StartPointVec_type_1, "test_SpecMCMC_StartPointVec_type_1")
        call Test%run(test_SpecMCMC_StartPointVec_type_2, "test_SpecMCMC_StartPointVec_type_2")
        call Test%run(test_SpecMCMC_StartPointVec_type_3, "test_SpecMCMC_StartPointVec_type_3")
        call Test%run(test_SpecMCMC_StartPointVec_type_4, "test_SpecMCMC_StartPointVec_type_4")
        call Test%run(test_SpecMCMC_StartPointVec_type_5, "test_SpecMCMC_StartPointVec_type_5")
        call Test%run(test_SpecMCMC_ProposalStartCovMat_type_1, "test_SpecMCMC_ProposalStartCovMat_type_1")
        call Test%run(test_SpecMCMC_ProposalStartCovMat_type_2, "test_SpecMCMC_ProposalStartCovMat_type_2")
        call Test%run(test_SpecMCMC_ProposalStartCovMat_type_3, "test_SpecMCMC_ProposalStartCovMat_type_3")
        call Test%run(test_SpecMCMC_ProposalStartCovMat_type_4, "test_SpecMCMC_ProposalStartCovMat_type_4")
        call Test%run(test_SpecMCMC_ProposalStartCovMat_type_5, "test_SpecMCMC_ProposalStartCovMat_type_5")
        call Test%run(test_SpecMCMC_ProposalStartCorMat_type_1, "test_SpecMCMC_ProposalStartCorMat_type_1")
        call Test%run(test_SpecMCMC_ProposalStartCorMat_type_2, "test_SpecMCMC_ProposalStartCorMat_type_2")
        call Test%run(test_SpecMCMC_ProposalStartCorMat_type_3, "test_SpecMCMC_ProposalStartCorMat_type_3")
        call Test%run(test_SpecMCMC_ProposalStartCorMat_type_4, "test_SpecMCMC_ProposalStartCorMat_type_4")
        call Test%run(test_SpecMCMC_ProposalStartCorMat_type_5, "test_SpecMCMC_ProposalStartCorMat_type_5")
        call Test%run(test_SpecMCMC_ProposalStartStdVec_type_1, "test_SpecMCMC_ProposalStartStdVec_type_1")
        call Test%run(test_SpecMCMC_ProposalStartStdVec_type_2, "test_SpecMCMC_ProposalStartStdVec_type_2")
        call Test%run(test_SpecMCMC_ProposalStartStdVec_type_3, "test_SpecMCMC_ProposalStartStdVec_type_3")
        call Test%run(test_SpecMCMC_ProposalStartStdVec_type_4, "test_SpecMCMC_ProposalStartStdVec_type_4")
        call Test%run(test_SpecMCMC_ProposalStartStdVec_type_5, "test_SpecMCMC_ProposalStartStdVec_type_5")
        call Test%run(test_SpecMCMC_SampleRefinementCount_type_1, "test_SpecMCMC_SampleRefinementCount_type_1")
        call Test%run(test_SpecMCMC_SampleRefinementCount_type_2, "test_SpecMCMC_SampleRefinementCount_type_2")
        call Test%run(test_SpecMCMC_SampleRefinementCount_type_3, "test_SpecMCMC_SampleRefinementCount_type_3")
        call Test%run(test_SpecMCMC_SampleRefinementMethod_type_1, "test_SpecMCMC_SampleRefinementMethod_type_1")
        call Test%run(test_SpecMCMC_SampleRefinementMethod_type_2, "test_SpecMCMC_SampleRefinementMethod_type_2")
        call Test%run(test_SpecMCMC_SampleRefinementMethod_type_3, "test_SpecMCMC_SampleRefinementMethod_type_3")
        call Test%run(test_SpecMCMC_SampleRefinementMethod_type_4, "test_SpecMCMC_SampleRefinementMethod_type_4")
        call Test%run(test_SpecMCMC_SampleRefinementMethod_type_5, "test_SpecMCMC_SampleRefinementMethod_type_5")
        call Test%run(test_SpecMCMC_SampleRefinementMethod_type_6, "test_SpecMCMC_SampleRefinementMethod_type_6")
        call Test%run(test_SpecMCMC_SampleRefinementMethod_type_7, "test_SpecMCMC_SampleRefinementMethod_type_7")
        call Test%run(test_SpecMCMC_SampleRefinementMethod_type_8, "test_SpecMCMC_SampleRefinementMethod_type_8")
        call Test%run(test_SpecMCMC_SampleRefinementMethod_type_9, "test_SpecMCMC_SampleRefinementMethod_type_9")
        call Test%run(test_SpecMCMC_SampleRefinementMethod_type_10, "test_SpecMCMC_SampleRefinementMethod_type_10")
        call Test%run(test_SpecMCMC_SampleRefinementMethod_type_11, "test_SpecMCMC_SampleRefinementMethod_type_11")
        call Test%run(test_SpecMCMC_SampleRefinementMethod_type_12, "test_SpecMCMC_SampleRefinementMethod_type_12")
        call Test%run(test_SpecMCMC_SampleRefinementMethod_type_13, "test_SpecMCMC_SampleRefinementMethod_type_13")
        call Test%run(test_SpecMCMC_RandomStartPointRequested_type_1, "test_SpecMCMC_RandomStartPointRequested_type_1")
        call Test%run(test_SpecMCMC_RandomStartPointRequested_type_2, "test_SpecMCMC_RandomStartPointRequested_type_2")
        call Test%run(test_SpecMCMC_RandomStartPointRequested_type_3, "test_SpecMCMC_RandomStartPointRequested_type_3")
        call Test%run(test_SpecMCMC_RandomStartPointRequested_type_4, "test_SpecMCMC_RandomStartPointRequested_type_4")
        call Test%run(test_RSPDLowerLimitVec_type_1, "test_RSPDLowerLimitVec_type_1") ! global names shortened to accommodate the file length limits of Intel ifort on windows.
        call Test%run(test_RSPDLowerLimitVec_type_2, "test_RSPDLowerLimitVec_type_2") ! global names shortened to accommodate the file length limits of Intel ifort on windows.
        call Test%run(test_RSPDLowerLimitVec_type_3, "test_RSPDLowerLimitVec_type_3") ! global names shortened to accommodate the file length limits of Intel ifort on windows.
        call Test%run(test_RSPDUpperLimitVec_type_1, "test_RSPDUpperLimitVec_type_1") ! global names shortened to accommodate the file length limits of Intel ifort on windows.
        call Test%run(test_RSPDUpperLimitVec_type_2, "test_RSPDUpperLimitVec_type_2") ! global names shortened to accommodate the file length limits of Intel ifort on windows.
        call Test%run(test_RSPDUpperLimitVec_type_3, "test_RSPDUpperLimitVec_type_3") ! global names shortened to accommodate the file length limits of Intel ifort on windows.
        call Test%run(test_RSPDUpperLimitVec_type_4, "test_RSPDUpperLimitVec_type_4") ! global names shortened to accommodate the file length limits of Intel ifort on windows.

        call Test%run(test_SpecDRAM_AdaptiveUpdateCount_type_1, "test_SpecDRAM_AdaptiveUpdateCount_type_1")
        call Test%run(test_SpecDRAM_AdaptiveUpdateCount_type_2, "test_SpecDRAM_AdaptiveUpdateCount_type_2")
        call Test%run(test_SpecDRAM_AdaptiveUpdateCount_type_3, "test_SpecDRAM_AdaptiveUpdateCount_type_3")
        call Test%run(test_SpecDRAM_AdaptiveUpdatePeriod_type_1, "test_SpecDRAM_AdaptiveUpdatePeriod_type_1")
        call Test%run(test_SpecDRAM_AdaptiveUpdatePeriod_type_2, "test_SpecDRAM_AdaptiveUpdatePeriod_type_2")
        call Test%run(test_SpecDRAM_AdaptiveUpdatePeriod_type_3, "test_SpecDRAM_AdaptiveUpdatePeriod_type_3")
        call Test%run(test_SpecDRAM_DelayedRejectionCount_type_1, "test_SpecDRAM_DelayedRejectionCount_type_1")
        call Test%run(test_SpecDRAM_DelayedRejectionCount_type_2, "test_SpecDRAM_DelayedRejectionCount_type_2")
        call Test%run(test_SpecDRAM_DelayedRejectionCount_type_3, "test_SpecDRAM_DelayedRejectionCount_type_3")
        call Test%run(test_SpecDRAM_DelayedRejectionCount_type_4, "test_SpecDRAM_DelayedRejectionCount_type_4")
        call Test%run(test_SpecDRAM_GreedyAdaptationCount_type_1, "test_SpecDRAM_GreedyAdaptationCount_type_1")
        call Test%run(test_SpecDRAM_GreedyAdaptationCount_type_2, "test_SpecDRAM_GreedyAdaptationCount_type_2")
        call Test%run(test_SpecDRAM_GreedyAdaptationCount_type_3, "test_SpecDRAM_GreedyAdaptationCount_type_3")
        call Test%run(test_SpecDRAM_BurninAdaptationMeasure_type_1, "test_SpecDRAM_BurninAdaptationMeasure_type_1")
        call Test%run(test_SpecDRAM_BurninAdaptationMeasure_type_2, "test_SpecDRAM_BurninAdaptationMeasure_type_2")
        call Test%run(test_SpecDRAM_BurninAdaptationMeasure_type_3, "test_SpecDRAM_BurninAdaptationMeasure_type_3")
        call Test%run(test_SpecDRAM_BurninAdaptationMeasure_type_4, "test_SpecDRAM_BurninAdaptationMeasure_type_4")
        call Test%run(test_SpecDRAM_DelayedRejectionScaleFactorVec_type_1, "test_SpecDRAM_DelayedRejectionScaleFactorVec_type_1")
        call Test%run(test_SpecDRAM_DelayedRejectionScaleFactorVec_type_2, "test_SpecDRAM_DelayedRejectionScaleFactorVec_type_2")
        call Test%run(test_SpecDRAM_DelayedRejectionScaleFactorVec_type_3, "test_SpecDRAM_DelayedRejectionScaleFactorVec_type_3")
        call Test%run(test_SpecDRAM_DelayedRejectionScaleFactorVec_type_4, "test_SpecDRAM_DelayedRejectionScaleFactorVec_type_4")
        call Test%run(test_SpecDRAM_DelayedRejectionScaleFactorVec_type_5, "test_SpecDRAM_DelayedRejectionScaleFactorVec_type_5")

        call Test%finalize()
    end subroutine test_ParaDXXX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! end module Test_ParaDXXX_mod