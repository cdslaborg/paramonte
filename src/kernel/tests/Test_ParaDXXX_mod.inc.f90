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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_ParaDXXX()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)

        call Test%run(test_runSampler_1, "test_runSampler_1")
        call Test%run(test_runSampler_2, "test_runSampler_2")
        call Test%run(test_runSampler_3, "test_runSampler_3")
        call Test%run(test_runSampler_4, "test_runSampler_4")
        call Test%run(test_runSampler_5, "test_runSampler_5")
        call Test%run(test_runSampler_6, "test_runSampler_6")
        call Test%run(test_runSampler_7, "test_runSampler_7")
        call Test%run(test_runSampler_8, "test_runSampler_8")
        call Test%run(test_runSampler_9, "test_runSampler_9")
        call Test%run(test_runSampler_10, "test_runSampler_10")
        call Test%run(test_runSampler_11, "test_runSampler_11")
        call Test%run(test_runSampler_12, "test_runSampler_12")
        call Test%run(test_runSampler_13, "test_runSampler_13")
        call Test%run(test_runSampler_14, "test_runSampler_14")
        call Test%run(test_runSampler_15, "test_runSampler_15")

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
        call Test%run(test_RSPDLowerLimitVec_type_1, "test_RSPDLowerLimitVec_type_1")
        call Test%run(test_RSPDLowerLimitVec_type_2, "test_RSPDLowerLimitVec_type_2")
        call Test%run(test_RSPDLowerLimitVec_type_3, "test_RSPDLowerLimitVec_type_3")
        call Test%run(test_RSPDUpperLimitVec_type_1, "test_RSPDUpperLimitVec_type_1")
        call Test%run(test_RSPDUpperLimitVec_type_2, "test_RSPDUpperLimitVec_type_2")
        call Test%run(test_RSPDUpperLimitVec_type_3, "test_RSPDUpperLimitVec_type_3")
        call Test%run(test_RSPDUpperLimitVec_type_4, "test_RSPDUpperLimitVec_type_4")

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

#include "Test_ParaDXXX_mod@SpecBase_mod.inc.f90"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "Test_ParaDXXX_mod@SpecMCMC_mod.inc.f90"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "Test_ParaDXXX_mod@SpecDRAM_mod.inc.f90"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with no input arguments or input file.
    function test_runSampler_1() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_1" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_runSampler_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with an internal input file.
    function test_runSampler_2() result(assertion)
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: chainSize_ref = 100_IK
        integer(IK) , parameter :: userSeed_ref = 1111_IK

        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_2" &
                            , inputFile = ParaDXXX_NML//" randomSeed = "//num2str(userSeed_ref)//" chainSize = "//num2str(chainSize_ref)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. PD%SpecBase%RandomSeed%userSeed==userSeed_ref .and. PD%SpecMCMC%ChainSize%val==chainSize_ref
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "PD%Err%occurred                   = ", PD%Err%occurred
            write(Test%outputUnit,"(*(g0))") "PD%Err%msg                        = ", PD%Err%msg
            write(Test%outputUnit,"(*(g0))") "userSeed_ref                      = ", userSeed_ref
            write(Test%outputUnit,"(*(g0))") "PD%SpecBase%RandomSeed%userSeed   = ", PD%SpecBase%RandomSeed%userSeed
            write(Test%outputUnit,"(*(g0))") "chainSize_ref                     = ", chainSize_ref
            write(Test%outputUnit,"(*(g0))") "PD%SpecMCMC%ChainSize%val         = ", PD%SpecMCMC%ChainSize%val
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
#endif
    end function test_runSampler_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with an external input file.
    function test_runSampler_3() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , inputFile = Test%inDir//"/Test_ParaDXXX_mod@test_runSampler_3.in" &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_3" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_runSampler_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with a path to a non-existing external input file.
    function test_runSampler_4() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_4" &
                            , inputFile = " " &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_runSampler_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with wrong SpecBase input values.
    function test_runSampler_5() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_5" &
                            , inputFile = Test%inDir//"/Test_ParaDXXX_mod@test_runSampler_5.in" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_runSampler_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler to restart a complete simulation without an output sample file.
    function test_runSampler_6() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDXXX_type) :: PD1, PD2
        assertion = .true.
#if defined CODECOV_ENABLED

        ! Run the fresh simulation

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_6" &
                            , outputRealPrecision = 16_IK &
                            , randomSeed = 12345_IK &
                            , sampleSize = 500_IK &
                            , chainSize = 1000_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        ! delete the sample file

        open(newunit = PD1%SampleFile%unit, file = PD1%SampleFile%Path%original, status = "replace")
        close(PD1%SampleFile%unit, status = "delete")

        ! restart the simulation with the same configuration

        call PD2%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_6" &
                            , outputRealPrecision = 16_IK &
                            , randomSeed = 12345_IK &
                            , sampleSize = 500_IK &
                            , chainSize = 1000_IK &
                            )
        assertion = assertion .and. .not. PD2%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        assertion = assertion .and. all( abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState) < 1.e-12_RK )
#endif
    end function test_runSampler_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with a wrong internal input namelist group:
    !> +    Infinity values for `domainLowerLimitVec` and `domainUpperLimitVec`.
    function test_runSampler_7() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , inputFile = "&ParaDXXX randomSeed = 1111 /" &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecBase/test_runSampler_7" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_runSampler_7

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with a wrong internal input namelist group:
    !> +    Infinity values for `domainLowerLimitVec` and `domainUpperLimitVec`.
    function test_runSampler_8() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , inputFile = Test%inDir//"/Test_ParaDXXX_mod@test_runSampler_8.in" &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecBase/test_runSampler_8" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_runSampler_8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler to restart a complete simulation without an output sample file.
    function test_runSampler_9() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDXXX_type) :: PD1, PD2
        assertion = .true.
#if defined CODECOV_ENABLED

        ! Run the fresh simulation as a reference run

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_9.ref" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 500_IK &
                            , chainSize = 1000_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        ! Run the fresh simulation

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_9" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 500_IK &
                            , chainSize = 1000_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        ! delete the sample file

        open(newunit = PD1%SampleFile%unit, file = PD1%SampleFile%Path%original, status = "replace")
        close(PD1%SampleFile%unit, status = "delete")

        ! restart the simulation with the same configuration

        call PD2%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_9" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 500_IK &
                            , chainSize = 1000_IK &
                            )
        assertion = assertion .and. .not. PD2%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        assertion = assertion .and. all( abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState) < 1.e-6_RK ) ! by default, the output precision is only 8 digits
#endif
    end function test_runSampler_9

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler to restart a complete simulation without an output sample file, with an ASCII output file.
    function test_runSampler_10() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDXXX_type) :: PD1, PD2
        assertion = .true.
#if defined CODECOV_ENABLED

        ! Run the fresh simulation as a reference run

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_10.ref" &
                            , restartFileFormat = "ascii" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 500_IK &
                            , chainSize = 1000_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        ! Run the fresh simulation

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_10" &
                            , restartFileFormat = "ascii" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 500_IK &
                            , chainSize = 1000_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        ! delete the sample file

        open(newunit = PD1%SampleFile%unit, file = PD1%SampleFile%Path%original, status = "replace")
        close(PD1%SampleFile%unit, status = "delete")

        ! restart the simulation with the same configuration

        call PD2%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_10" &
                            , restartFileFormat = "ascii" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 500_IK &
                            , chainSize = 1000_IK &
                            )
        assertion = assertion .and. .not. PD2%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        assertion = assertion .and. all( abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState) < 1.e-6_RK ) ! by default, the output precision is only 8 digits
#endif
    end function test_runSampler_10

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler to restart a complete simulation without an 
    !> output sample file, with an ASCII output file but with binary chain file.
    function test_runSampler_11() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDXXX_type) :: PD1, PD2
        assertion = .true.
#if defined CODECOV_ENABLED

        ! Run the fresh simulation as a reference run

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_11.ref" &
                            , restartFileFormat = "ascii" &
                            , chainFileFormat = "verbose" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 500_IK &
                            , chainSize = 1000_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        ! Run the fresh simulation

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_11" &
                            , restartFileFormat = "ascii" &
                            , chainFileFormat = "verbose" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 500_IK &
                            , chainSize = 1000_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        ! delete the sample file

        open(newunit = PD1%SampleFile%unit, file = PD1%SampleFile%Path%original, status = "replace")
        close(PD1%SampleFile%unit, status = "delete")

        ! restart the simulation with the same configuration

        call PD2%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_11" &
                            , restartFileFormat = "ascii" &
                            , chainFileFormat = "verbose" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 500_IK &
                            , chainSize = 1000_IK &
                            )
        assertion = assertion .and. .not. PD2%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        assertion = assertion .and. all( abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState) < 1.e-6_RK ) ! by default, the output precision is only 8 digits
#endif
    end function test_runSampler_11

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler to restart a complete simulation without an 
    !> output sample file, with an ASCII output file but with binary chain file.
    function test_runSampler_12() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDXXX_type) :: PD1, PD2
        assertion = .true.
#if defined CODECOV_ENABLED

        ! Run the fresh simulation as a reference run

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_12.ref" &
                            , restartFileFormat = "ascii" &
                            , chainFileFormat = "binary" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 500_IK &
                            , chainSize = 1000_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        ! Run the fresh simulation

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_12" &
                            , restartFileFormat = "ascii" &
                            , chainFileFormat = "binary" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 500_IK &
                            , chainSize = 1000_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        ! delete the sample file

        open(newunit = PD1%SampleFile%unit, file = PD1%SampleFile%Path%original, status = "replace")
        close(PD1%SampleFile%unit, status = "delete")

        ! restart the simulation with the same configuration

        call PD2%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_12" &
                            , restartFileFormat = "ascii" &
                            , chainFileFormat = "binary" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 500_IK &
                            , chainSize = 1000_IK &
                            )
        assertion = assertion .and. .not. PD2%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        assertion = assertion .and. all( abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState) < 1.e-6_RK ) ! by default, the output precision is only 8 digits
#endif
    end function test_runSampler_12

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with an empty external input file. 
    !> This should first search for the specific namelist, then for the generic namelist, 
    !> then issues a warning and runs the simulation without setting the simulation specifications.
    function test_runSampler_13() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , inputFile = Test%inDir//"/Test_ParaDXXX_mod@test_runSampler_13.in" &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_13" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_runSampler_13

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with an empty external input file. 
    !> This should first search for the specific namelist, then for the generic namelist, 
    !> then issues a warning and runs the simulation without setting the simulation specifications.
    function test_runSampler_14() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , inputFile = "&nonsense /" &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_14" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_runSampler_14

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the read method of the `ParaMCMCRefinedChain_type` class can successfully read an external sample file. 
    function test_runSampler_15() result(assertion)
        use ParaMCMCRefinedChain_mod, only: readRefinedChain, RefinedChain_type
        implicit none
        logical                 :: assertion
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        integer(IK) , parameter :: NDIM = 2_IK
        character(*), parameter :: DELIM = "delim"
        type(ParaDXXX_type)     :: PD
        type(RefinedChain_type) :: RefinedChain
        real(RK), allocatable   :: Difference(:,:)
        assertion = .true.
#if defined CODECOV_ENABLED

        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_15" &
                            , mpiFinalizeRequested = .false. &
                            , outputRealPrecision = 15_IK &
                            , outputDelimiter = DELIM &
                            , sampleSize = 10_IK &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        RefinedChain = readRefinedChain( sampleFilePath=PD%SampleFile%Path%original, delimiter = PD%SpecBase%OutputDelimiter%val, ndim = PD%nd%val )

        assertion = assertion .and. RefinedChain%numRefinement == 0_IK

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(*,"(10(g0,:,', '))")
            write(*,"(10(g0,:,', '))") "RefinedChain%numRefinement", RefinedChain%numRefinement
            write(*,"(10(g0,:,', '))")
            return
        end if
        ! LCOV_EXCL_STOP

        assertion = assertion .and. all(shape(RefinedChain%LogFuncState) == shape(PD%RefinedChain%LogFuncState))

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(*,"(10(g0,:,', '))")
            write(*,"(10(g0,:,', '))") "shape(PD%RefinedChain%LogFuncState) ", shape(PD%RefinedChain%LogFuncState)
            write(*,"(10(g0,:,', '))") "shape(RefinedChain%LogFuncState)    ", shape(RefinedChain%LogFuncState)
            write(*,"(10(g0,:,', '))")
            return
        end if
        ! LCOV_EXCL_STOP

        Difference = abs(RefinedChain%LogFuncState-PD%RefinedChain%LogFuncState)
        assertion = assertion .and. all(Difference < tolerance)

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(*,"(10(g0,:,', '))")
            write(*,"(10(g0,:,', '))") "PD%RefinedChain%LogFuncState", PD%RefinedChain%LogFuncState
            write(*,"(10(g0,:,', '))") "RefinedChain%LogFuncState   ", RefinedChain%LogFuncState
            write(*,"(10(g0,:,', '))") "Difference                  ", Difference
            write(*,"(10(g0,:,', '))")
            return
        end if
        ! LCOV_EXCL_STOP

#endif
    end function test_runSampler_15

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! end module Test_ParaDXXX_mod