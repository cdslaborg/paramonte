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

!>  \brief This include file contains tests of the module [SpecBase_mod](@ref specbase_mod) when accessed by the ParaDRAM sampler.
!>  @author Amir Shahmoradi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a wrong input value for `chainFileFormat`.
    function test_SpecBase_ChainFileFormat_type_1() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , chainFileFormat = "nonsense" &
                            , outputFileName = "test_SpecBase_ChainFileFormat_type_1" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_ChainFileFormat_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a wrong input value for `chainFileFormat`.
    function test_SpecBase_ChainFileFormat_type_2() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , inputFile = "&ParaDRAM chainFileFormat = 'nonsense' /" &
                            , outputFileName = "test_SpecBase_ChainFileFormat_type_2" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_ChainFileFormat_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a too-small input value for `domainLowerLimitVec`.
    function test_SpecBase_DomainLowerLimitVec_type_1() result(assertion)
        use Constants_mod, only: RK, HUGE_RK
        implicit none
        logical             :: assertion
        real(RK), parameter :: domainLowerLimitVec(*) = [-HUGE_RK/2._RK] ! NOTE: HUGE_RK is the null value.
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , domainLowerLimitVec = domainLowerLimitVec &
                            , outputFileName = "test_SpecBase_DomainLowerLimitVec_type_1" &
                            )
        assertion = PD%Err%occurred
#endif
    end function test_SpecBase_DomainLowerLimitVec_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a too-small input value for `domainLowerLimitVec`.
    function test_SpecBase_DomainLowerLimitVec_type_2() result(assertion)
        use Constants_mod, only: RK, HUGE_RK
        use String_mod, only: num2str
        implicit none
        logical             :: assertion
        real(RK), parameter :: domainLowerLimitVec(*) = [-HUGE_RK/2._RK] ! NOTE: HUGE_RK is the null value.
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , inputFile = "&ParaDRAM domainLowerLimitVec = "//num2str(domainLowerLimitVec)//" /" &
                            , outputFileName = "test_SpecBase_DomainLowerLimitVec_type_2" &
                            )
        assertion = PD%Err%occurred
#endif
    end function test_SpecBase_DomainLowerLimitVec_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a too-large input value for `domainUpperLimitVec`.
    function test_SpecBase_DomainUpperLimitVec_type_1() result(assertion)
        use Constants_mod, only: RK, HUGE_RK
        implicit none
        logical             :: assertion
        real(RK), parameter :: domainUpperLimitVec(*) = [+HUGE_RK/2._RK] ! NOTE: HUGE_RK is the null value.
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , domainUpperLimitVec = domainUpperLimitVec &
                            , outputFileName = "test_SpecBase_DomainUpperLimitVec_type_1" &
                            )
        assertion = PD%Err%occurred
#endif
    end function test_SpecBase_DomainUpperLimitVec_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a too-large input value for `domainUpperLimitVec`.
    function test_SpecBase_DomainUpperLimitVec_type_2() result(assertion)
        use Constants_mod, only: RK, HUGE_RK
        use String_mod, only: num2str
        implicit none
        logical             :: assertion
        real(RK), parameter :: domainUpperLimitVec(*) = [+HUGE_RK/2._RK] ! NOTE: HUGE_RK is the null value.
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , inputFile = "&ParaDRAM domainUpperLimitVec = "//num2str(domainUpperLimitVec)//" /" &
                            , outputFileName = "test_SpecBase_DomainUpperLimitVec_type_2" &
                            )
        assertion = PD%Err%occurred
#endif
    end function test_SpecBase_DomainUpperLimitVec_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with an input value for `domainUpperLimitVec` 
    !> that is smaller than the input value for `domainLowerLimitVec`.
    function test_SpecBase_DomainUpperLimitVec_type_3() result(assertion)
        use Constants_mod, only: RK, HUGE_RK
        implicit none
        logical             :: assertion
        real(RK), parameter :: domainLowerLimitVec(*) = [+2._RK] ! NOTE: HUGE_RK is the null value.
        real(RK), parameter :: domainUpperLimitVec(*) = [-2._RK] ! NOTE: HUGE_RK is the null value.
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , domainLowerLimitVec = domainLowerLimitVec &
                            , domainUpperLimitVec = domainUpperLimitVec &
                            , outputFileName = "test_SpecBase_DomainUpperLimitVec_type_3" &
                            )
        assertion = PD%Err%occurred
#endif
    end function test_SpecBase_DomainUpperLimitVec_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with an input value for `domainUpperLimitVec` 
    !> that is equal to the input value for `domainLowerLimitVec`.
    function test_SpecBase_DomainUpperLimitVec_type_4() result(assertion)
        use Constants_mod, only: RK, HUGE_RK
        implicit none
        logical             :: assertion
        real(RK), parameter :: domainLowerLimitVec(*) = [+2._RK] ! NOTE: HUGE_RK is the null value.
        real(RK), parameter :: domainUpperLimitVec(*) = [+2._RK] ! NOTE: HUGE_RK is the null value.
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , domainLowerLimitVec = domainLowerLimitVec &
                            , domainUpperLimitVec = domainUpperLimitVec &
                            , outputFileName = "test_SpecBase_DomainUpperLimitVec_type_4" &
                            )
        assertion = PD%Err%occurred
#endif
    end function test_SpecBase_DomainUpperLimitVec_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a wrong internal input namelist group:
    !> +    Infinity values for `domainLowerLimitVec` and `domainUpperLimitVec`.
    function test_runSampler_7() result(assertion)
        use Constants_mod, only: RK, HUGE_RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , inputFile = "&ParaDXXX randomSeed = 1111 /" &
                            , outputFileName = "test_runSampler_7" &
                            )
        assertion = .not. PD%Err%occurred
#endif
    end function test_runSampler_7

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a wrong internal input namelist group:
    !> +    Infinity values for `domainLowerLimitVec` and `domainUpperLimitVec`.
    function test_runSampler_8() result(assertion)
        use Constants_mod, only: RK, HUGE_RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , inputFile = Test%inDir//"/Test_ParaDRAM_mod@test_runSampler_8.in" &
                            , outputFileName = "test_runSampler_8" &
                            )
        assertion = .not. PD%Err%occurred
#endif
    end function test_runSampler_8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler quits with an error message when `maxNumDomainCheckToWarn` has reached.
    function test_runSampler_9() result(assertion)
        use Constants_mod, only: RK, HUGE_RK
        implicit none
        logical             :: assertion
        real(RK), parameter :: domainLowerLimitVec(*) = [-1.e-1_RK] ! NOTE: HUGE_RK is the null value.
        real(RK), parameter :: domainUpperLimitVec(*) = [+1.e+1_RK] ! NOTE: HUGE_RK is the null value.
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_runSampler_9" &
                            , domainLowerLimitVec = domainLowerLimitVec &
                            , domainUpperLimitVec = domainUpperLimitVec &
                            , maxNumDomainCheckToWarn = 1_IK &
                            )
        assertion = PD%Err%occurred
#endif
    end function test_runSampler_9

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler quits with an error message when `maxNumDomainCheckToStop` has reached.
    function test_runSampler_10() result(assertion)
        use Constants_mod, only: RK, HUGE_RK
        implicit none
        logical             :: assertion
        real(RK), parameter :: domainLowerLimitVec(*) = [-1.e-1_RK] ! NOTE: HUGE_RK is the null value.
        real(RK), parameter :: domainUpperLimitVec(*) = [+1.e+1_RK] ! NOTE: HUGE_RK is the null value.
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_runSampler_10" &
                            , domainLowerLimitVec = domainLowerLimitVec &
                            , domainUpperLimitVec = domainUpperLimitVec &
                            , maxNumDomainCheckToStop = 1_IK &
                            )
        assertion = PD%Err%occurred
#endif
    end function test_runSampler_10

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
