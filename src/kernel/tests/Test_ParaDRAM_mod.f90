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
!>  @author Amir Shahmoradi

module Test_ParaDRAM_mod

    use Constants_mod, only: IK, RK
    use Test_mod, only: Test_type
    use ParaDRAM_mod

    !use Statistics_mod, only: paradramPrintEnabled
    implicit none
    !paradramPrintEnabled = .true.

    private
    public :: test_ParaDRAM

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_ParaDRAM()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)

        call Test%run(test_runSampler_1, "test_runSampler_1")
        call Test%run(test_runSampler_2, "test_runSampler_2")
        call Test%run(test_runSampler_3, "test_runSampler_3")
        call Test%run(test_runSampler_4, "test_runSampler_4")
        call Test%run(test_runSampler_5, "test_runSampler_5")
        call Test%run(test_runSampler_7, "test_runSampler_7")
        call Test%run(test_runSampler_8, "test_runSampler_8")

        call Test%run(test_SpecBase_RandomSeed_type_1, "test_SpecBase_RandomSeed_type_1")
        call Test%run(test_SpecBase_RandomSeed_type_2, "test_SpecBase_RandomSeed_type_2")
        call Test%run(test_SpecBase_OutputDelimiter_type_1, "test_SpecBase_OutputDelimiter_type_1")
        call Test%run(test_SpecBase_OutputDelimiter_type_2, "test_SpecBase_OutputDelimiter_type_2")
        call Test%run(test_SpecBase_OutputDelimiter_type_3, "test_SpecBase_OutputDelimiter_type_3")
        call Test%run(test_SpecBase_OutputDelimiter_type_4, "test_SpecBase_OutputDelimiter_type_4")
        call Test%run(test_SpecBase_OutputDelimiter_type_5, "test_SpecBase_OutputDelimiter_type_5")
        call Test%run(test_SpecBase_OutputDelimiter_type_6, "test_SpecBase_OutputDelimiter_type_6")
        call Test%run(test_SpecBase_ChainFileFormat_type_1, "test_SpecBase_ChainFileFormat_type_1")
        call Test%run(test_SpecBase_ChainFileFormat_type_2, "test_SpecBase_ChainFileFormat_type_2")
        call Test%run(test_SpecBase_OutputColumnWidth_type_1, "test_SpecBase_OutputColumnWidth_type_1")
        call Test%run(test_SpecBase_OutputColumnWidth_type_2, "test_SpecBase_OutputColumnWidth_type_2")
        call Test%run(test_SpecBase_OutputColumnWidth_type_3, "test_SpecBase_OutputColumnWidth_type_3")
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
        call Test%run(test_SpecMCMC_RandomStartPointRequested_type_1, "test_SpecMCMC_RandomStartPointRequested_type_1")
        call Test%run(test_SpecMCMC_RandomStartPointRequested_type_2, "test_SpecMCMC_RandomStartPointRequested_type_2")
        call Test%run(test_SpecMCMC_RandomStartPointRequested_type_3, "test_SpecMCMC_RandomStartPointRequested_type_3")
        call Test%run(test_SpecMCMC_RandomStartPointRequested_type_4, "test_SpecMCMC_RandomStartPointRequested_type_4")
        call Test%run(test_SpecMCMC_RandomStartPointDomainLowerLimitVec_type_1, "test_SpecMCMC_RandomStartPointDomainLowerLimitVec_type_1")
        call Test%run(test_SpecMCMC_RandomStartPointDomainLowerLimitVec_type_2, "test_SpecMCMC_RandomStartPointDomainLowerLimitVec_type_2")
        call Test%run(test_SpecMCMC_RandomStartPointDomainLowerLimitVec_type_3, "test_SpecMCMC_RandomStartPointDomainLowerLimitVec_type_3")
        call Test%run(test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_1, "test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_1")
        call Test%run(test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_2, "test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_2")
        call Test%run(test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_3, "test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_3")
        call Test%run(test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_4, "test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_4")

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
    end subroutine test_ParaDRAM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "Test_ParaDRAM_mod@SpecBase_mod.inc.f90"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "Test_ParaDRAM_mod@SpecMCMC_mod.inc.f90"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "Test_ParaDRAM_mod@SpecDRAM_mod.inc.f90"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with no input arguments or input file.
    function test_runSampler_1() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_runSampler_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with an internal input file.
    function test_runSampler_2() result(assertion)
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: chainSize_ref = 100_IK
        integer(IK) , parameter :: userSeed_ref = 1111_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , inputFile = "&ParaDRAM randomSeed = "//num2str(userSeed_ref)//" chainSize = "//num2str(chainSize_ref)//" /" &
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
    !> Test the ParaDRAM sampler with an external input file.
    function test_runSampler_3() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , inputFile = Test%inDir//"/Test_ParaDRAM_mod@test_runSampler_3.in" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_runSampler_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a path to a non-existing external input file.
    function test_runSampler_4() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , inputFile = " " &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_runSampler_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with wrong SpecBase input values.
    function test_runSampler_5() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , inputFile = Test%inDir//"/Test_ParaDRAM_mod@test_runSampler_5.in" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_runSampler_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the density function value of the uncorrelated Multivariate Normal distribution for the given input
    !> location specified by the vector `Point` of length `ndim`.
#if defined CFI_ENABLED
    function getLogFuncMVN(ndim,Point) result(logFunc) bind(C)
#else
    function getLogFuncMVN(ndim,Point) result(logFunc)
#endif
        ! This function returns the probability density function of the standard multivariate normal distribution of ndim dimensions.
        !use Statistics_mod, only: getLogProbMixMVN
        use Constants_mod, only : IK, RK
        implicit none

        !! Standard MultiVariate Normal (SMVN) specifications: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
        !real(RK), parameter :: LOG_SMVN_COEF = NDIM*log(1._RK/sqrt(2._RK*acos(-1._RK))) ! log(1/sqrt(2*Pi)^ndim)

#if defined CFI_ENABLED
        integer(IK), intent(in), value  :: ndim
#else
        integer(IK), intent(in)         :: ndim
#endif
        real(RK), intent(in)            :: Point(ndim)
        real(RK)                        :: logFunc

        !block
        !    use System_mod, only: sleep
        !    use Err_mod, only: Err_type
        !    type(Err_type) :: Err
        !    call sleep(seconds=5000.e-6_RK,Err=Err)
        !end block

        !block
        !    real(RK), allocatable :: unifrnd(:,:)
        !    allocate(unifrnd(200,20))
        !    call random_number(unifrnd)
        !    logFunc = sum(unifrnd) - 0.5_RK * sum(Point**2) - sum(unifrnd)
        !    deallocate(unifrnd)
        !end block

        logFunc = -sum(Point**2)

       !block
       !    integer(IK), parameter :: nmode = 2_IK
       !    real(RK) :: LogAmplitude(nmode), MeanVec(nmode), InvCovMat(nmode), LogSqrtDetInvCovMat(nmode)
       !    LogAmplitude        = [1._RK, 1._RK]
       !    MeanVec             = [0._RK, 7._RK]
       !    InvCovMat           = [1._RK,1._RK]
       !    LogSqrtDetInvCovMat = [1._RK,1._RK]
       !    logFunc = getLogProbGausMix ( nmode = 2_IK &
       !                                , nd = 1_IK &
       !                                , np = 1_IK &
       !                                , LogAmplitude = LogAmplitude &
       !                                , MeanVec = MeanVec &
       !                                , InvCovMat = InvCovMat &
       !                                , LogSqrtDetInvCovMat = LogSqrtDetInvCovMat &
       !                                , Point = Point(1) &
       !                                )
       !end block
    end function getLogFuncMVN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_ParaDRAM_mod