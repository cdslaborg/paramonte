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
!!!!       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This include file contains the body of the submodules 
!>  [ParaDRAM_mod@Test_SpecDRAM_smod](@ref paradram_mod@test_specdram_smod) and 
!>  [ParaDISE_mod@@Test_SpecDRAM_smod](@ref paradise_mod@@test_specdram_smod).
!>  \author Amir Shahmoradi

#if defined PARADRAM

#define ParaXXXX_mod ParaDRAM_mod
#define ParaXXXX_type ParaDRAM_type
#define test_ParaXXXX test_ParaDRAM
#define ParaXXXX_NML "&ParaDRAM"
#define ParaXXXX ParaDRAM
#define ParaXXXX_RefinedChain_mod ParaDRAM_RefinedChain_mod

#elif defined PARADISE

#define ParaXXXX_mod ParaDISE_mod
#define ParaXXXX_type ParaDISE_type
#define test_ParaXXXX test_ParaDISE
#define ParaXXXX_NML "&ParaDISE"
#define ParaXXXX ParaDISE
#define ParaXXXX_RefinedChain_mod ParaDISE_RefinedChain_mod

#endif

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a wrong input value for `daptiveUpdateCount < 0`.
    module function test_SpecDRAM_AdaptiveUpdateCount_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type) :: PD
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_AdaptiveUpdateCount_type_1" &
                            , adaptiveUpdateCount = -1_IK &
                            )
        assertion = assertion .and. PD%Err%occurred
        end block
#endif
    end function test_SpecDRAM_AdaptiveUpdateCount_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a valid input value for `daptiveUpdateCount`.
    module function test_SpecDRAM_AdaptiveUpdateCount_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type) :: PD
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_AdaptiveUpdateCount_type_2" &
                            , adaptiveUpdateCount = 0_IK &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
        end block
#endif
    end function test_SpecDRAM_AdaptiveUpdateCount_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a valid input value for `daptiveUpdateCount`.
    module function test_SpecDRAM_AdaptiveUpdateCount_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type) :: PD
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_AdaptiveUpdateCount_type_3" &
                            , inputFile = ParaXXXX_NML//" adaptiveUpdateCount = 0 /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
        end block
#endif
    end function test_SpecDRAM_AdaptiveUpdateCount_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a wrong input value for `adaptiveUpdatePeriod < 0`.
    module function test_SpecDRAM_AdaptiveUpdatePeriod_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type) :: PD
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_AdaptiveUpdatePeriod_type_1" &
                            , adaptiveUpdatePeriod = 0_IK &
                            )
        assertion = assertion .and. PD%Err%occurred
        end block
#endif
    end function test_SpecDRAM_AdaptiveUpdatePeriod_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a valid input value for `adaptiveUpdatePeriod`.
    module function test_SpecDRAM_AdaptiveUpdatePeriod_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type) :: PD
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_AdaptiveUpdatePeriod_type_2" &
                            , adaptiveUpdatePeriod = 1_IK &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
        end block
#endif
    end function test_SpecDRAM_AdaptiveUpdatePeriod_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a valid input value for `adaptiveUpdatePeriod`.
    module function test_SpecDRAM_AdaptiveUpdatePeriod_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type) :: PD
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_AdaptiveUpdatePeriod_type_3" &
                            , inputFile = ParaXXXX_NML//" adaptiveUpdatePeriod = 1 /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
        end block
#endif
    end function test_SpecDRAM_AdaptiveUpdatePeriod_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a wrong input value for `burninAdaptationMeasure < 0.`.
    module function test_SpecDRAM_BurninAdaptationMeasure_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type) :: PD
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_BurninAdaptationMeasure_type_1" &
                            , burninAdaptationMeasure = -0.1_RK &
                            )
        assertion = assertion .and. PD%Err%occurred
        end block
#endif
    end function test_SpecDRAM_BurninAdaptationMeasure_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a wrong input value for `burninAdaptationMeasure > 1.`.
    module function test_SpecDRAM_BurninAdaptationMeasure_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type) :: PD
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_BurninAdaptationMeasure_type_2" &
                            , burninAdaptationMeasure = 1.1_RK &
                            )
        assertion = assertion .and. PD%Err%occurred
        end block
#endif
    end function test_SpecDRAM_BurninAdaptationMeasure_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a valid input value for `adaptiveUpdatePeriod`.
    module function test_SpecDRAM_BurninAdaptationMeasure_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type) :: PD
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_BurninAdaptationMeasure_type_3" &
                            , burninAdaptationMeasure = 1._RK &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
        end block
#endif
    end function test_SpecDRAM_BurninAdaptationMeasure_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a valid input value for `adaptiveUpdatePeriod`.
    module function test_SpecDRAM_BurninAdaptationMeasure_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type) :: PD
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_BurninAdaptationMeasure_type_4" &
                            , inputFile = ParaXXXX_NML//" burninAdaptationMeasure = 0.5 /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
        end block
#endif
    end function test_SpecDRAM_BurninAdaptationMeasure_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a wrong input value for `burninAdaptationMeasure < 0.`.
    module function test_SpecDRAM_DelayedRejectionCount_type_1() result(assertion)
        use SpecDRAM_DelayedRejectionCount_mod, only: MIN_DELAYED_REJECTION_COUNT
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type) :: PD
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_DelayedRejectionCount_type_1" &
                            , delayedRejectionCount = MIN_DELAYED_REJECTION_COUNT - 1_IK &
                            )
        assertion = assertion .and. PD%Err%occurred
        end block
#endif
    end function test_SpecDRAM_DelayedRejectionCount_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a wrong input value for `DelayedRejectionCount > 1.`.
    module function test_SpecDRAM_DelayedRejectionCount_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use SpecDRAM_DelayedRejectionCount_mod, only: MAX_DELAYED_REJECTION_COUNT
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type) :: PD
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_DelayedRejectionCount_type_2" &
                            , delayedRejectionCount = MAX_DELAYED_REJECTION_COUNT + 1_IK &
                            )
        assertion = assertion .and. PD%Err%occurred
        end block
#endif
    end function test_SpecDRAM_DelayedRejectionCount_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a valid input value for `adaptiveUpdatePeriod`.
    module function test_SpecDRAM_DelayedRejectionCount_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type) :: PD
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_DelayedRejectionCount_type_3" &
                            , delayedRejectionCount = 0_IK &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
        end block
#endif
    end function test_SpecDRAM_DelayedRejectionCount_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a valid input value for `adaptiveUpdatePeriod`.
    module function test_SpecDRAM_DelayedRejectionCount_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type) :: PD
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_DelayedRejectionCount_type_4" &
                            , inputFile = ParaXXXX_NML//" delayedRejectionCount = 3 /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
        end block
#endif
    end function test_SpecDRAM_DelayedRejectionCount_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a valid input value for `DelayedRejectionScaleFactorVec`.
    module function test_SpecDRAM_DelayedRejectionScaleFactorVec_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type) :: PD
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_DelayedRejectionScaleFactorVec_type_1" &
                            , delayedRejectionCount = 3_IK &
                            , inputFile = ParaXXXX_NML//" DelayedRejectionScaleFactorVec = 3., 2., 1., /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
        end block
#endif
    end function test_SpecDRAM_DelayedRejectionScaleFactorVec_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a valid input value for `delayedRejectionCount`, which should lead to default values for
    !> the vector `DelayedRejectionScaleFactorVec`.
    module function test_SpecDRAM_DelayedRejectionScaleFactorVec_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        integer(IK), parameter  :: NDIM = 2_IK
        type(ParaXXXX_type)     :: PD
        call PD%runSampler  ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_DelayedRejectionScaleFactorVec_type_2" &
                            , inputFile = ParaXXXX_NML//" delayedRejectionCount = 3 /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(abs(PD%SpecDRAM%DelayedRejectionScaleFactorVec%Val-0.5_RK**(1._RK/real(NDIM,kind=RK)))<1.e-12_RK)
        end block
#endif
    end function test_SpecDRAM_DelayedRejectionScaleFactorVec_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a valid input value for `DelayedRejectionScaleFactorVec`.
    module function test_SpecDRAM_DelayedRejectionScaleFactorVec_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type)     :: PD
        real(RK), parameter     :: DelayedRejectionScaleFactorVec(*) = [3._RK, 2._RK]
        integer(IK), parameter  :: delayedRejectionCount = 2_IK
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_DelayedRejectionScaleFactorVec_type_3" &
                            , DelayedRejectionScaleFactorVec = DelayedRejectionScaleFactorVec &
                            , inputFile = ParaXXXX_NML//" delayedRejectionCount = "//num2str(delayedRejectionCount)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%SpecDRAM%DelayedRejectionScaleFactorVec%Val==DelayedRejectionScaleFactorVec(1:delayedRejectionCount))
        end block
#endif
    end function test_SpecDRAM_DelayedRejectionScaleFactorVec_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a valid input value for `DelayedRejectionScaleFactorVec`.
    module function test_SpecDRAM_DelayedRejectionScaleFactorVec_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type) :: PD
        real(RK), parameter :: DelayedRejectionScaleFactorVec(*) = [-1._RK, 2._RK]
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_DelayedRejectionScaleFactorVec_type_4" &
                            , DelayedRejectionScaleFactorVec = DelayedRejectionScaleFactorVec &
                            , inputFile = ParaXXXX_NML//" delayedRejectionCount = 2 /" &
                            )
        assertion = assertion .and. PD%Err%occurred
        end block
#endif
    end function test_SpecDRAM_DelayedRejectionScaleFactorVec_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a valid input value for `DelayedRejectionScaleFactorVec`.
    module function test_SpecDRAM_DelayedRejectionScaleFactorVec_type_5() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type) :: PD
        real(RK), parameter :: DelayedRejectionScaleFactorVec(*) = [2._RK, 0._RK, 1._RK]
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_DelayedRejectionScaleFactorVec_type_5" &
                            , DelayedRejectionScaleFactorVec = DelayedRejectionScaleFactorVec &
                            , inputFile = ParaXXXX_NML//" delayedRejectionCount = 2 /" &
                            )
        assertion = assertion .and. PD%Err%occurred
        end block
#endif
    end function test_SpecDRAM_DelayedRejectionScaleFactorVec_type_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a valid input value for `greedyAdaptationCount >= 0`.
    module function test_SpecDRAM_GreedyAdaptationCount_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type)     :: PD
        integer(IK), parameter  :: greedyAdaptationCount = 10_IK
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_GreedyAdaptationCount_type_1" &
                            , greedyAdaptationCount = greedyAdaptationCount &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. PD%SpecDRAM%greedyAdaptationCount%val==greedyAdaptationCount
        end block
#endif
    end function test_SpecDRAM_GreedyAdaptationCount_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a valid input value for `greedyAdaptationCount >= 0`.
    module function test_SpecDRAM_GreedyAdaptationCount_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type)     :: PD
        integer(IK), parameter  :: greedyAdaptationCount = 0_IK
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_GreedyAdaptationCount_type_2" &
                            , inputFile = ParaXXXX_NML//" greedyAdaptationCount = "//num2str(greedyAdaptationCount)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. PD%SpecDRAM%greedyAdaptationCount%val==greedyAdaptationCount
        end block
#endif
    end function test_SpecDRAM_GreedyAdaptationCount_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with an invalid input value for `greedyAdaptationCount >= 0`.
    module function test_SpecDRAM_GreedyAdaptationCount_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block
        type(ParaXXXX_type)     :: PD
        integer(IK), parameter  :: greedyAdaptationCount = -1_IK
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecDRAM/test_SpecDRAM_GreedyAdaptationCount_type_3" &
                            , inputFile = ParaXXXX_NML//" greedyAdaptationCount = "//num2str(greedyAdaptationCount)//" /" &
                            )
        assertion = assertion .and. PD%Err%occurred
        end block
#endif
    end function test_SpecDRAM_GreedyAdaptationCount_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ParaXXXX_mod
#undef ParaXXXX_type
#undef test_ParaXXXX
#undef ParaXXXX_NML
#undef ParaXXXX
#undef ParaXXXX_RefinedChain_mod
