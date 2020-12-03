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
        use Constants_mod, only: IK, RK, HUGE_RK
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
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_DomainLowerLimitVec_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a too-small input value for `domainLowerLimitVec`.
    function test_SpecBase_DomainLowerLimitVec_type_2() result(assertion)
        use Constants_mod, only: IK, RK, HUGE_RK
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
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_DomainLowerLimitVec_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a too-large input value for `domainUpperLimitVec`.
    function test_SpecBase_DomainUpperLimitVec_type_1() result(assertion)
        use Constants_mod, only: IK, RK, HUGE_RK
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
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_DomainUpperLimitVec_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a too-large input value for `domainUpperLimitVec`.
    function test_SpecBase_DomainUpperLimitVec_type_2() result(assertion)
        use Constants_mod, only: IK, RK, HUGE_RK
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
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_DomainUpperLimitVec_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with an input value for `domainUpperLimitVec` 
    !> that is smaller than the input value for `domainLowerLimitVec`.
    function test_SpecBase_DomainUpperLimitVec_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        real(RK), parameter :: domainLowerLimitVec(*) = [+2._RK]
        real(RK), parameter :: domainUpperLimitVec(*) = [-2._RK]
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
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_DomainUpperLimitVec_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with an input value for `domainUpperLimitVec` 
    !> that is equal to the input value for `domainLowerLimitVec`.
    function test_SpecBase_DomainUpperLimitVec_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        real(RK), parameter :: domainLowerLimitVec(*) = [+2._RK]
        real(RK), parameter :: domainUpperLimitVec(*) = [+2._RK]
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
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_DomainUpperLimitVec_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler quits with an error message when `maxNumDomainCheckToWarn < 1`.
    function test_SpecBase_MaxNumDomainCheckToWarn_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        real(RK), parameter :: domainLowerLimitVec(*) = [-1.e-1_RK]
        real(RK), parameter :: domainUpperLimitVec(*) = [+1.e+1_RK]
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_MaxNumDomainCheckToWarn_type_1" &
                            , domainLowerLimitVec = domainLowerLimitVec &
                            , domainUpperLimitVec = domainUpperLimitVec &
                            , maxNumDomainCheckToWarn = 0_IK &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_MaxNumDomainCheckToWarn_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler quits with an error message when `maxNumDomainCheckToWarn < 1`.
    function test_SpecBase_MaxNumDomainCheckToWarn_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        real(RK), parameter :: domainLowerLimitVec(*) = [-1.e-1_RK]
        real(RK), parameter :: domainUpperLimitVec(*) = [+1.e+1_RK]
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_MaxNumDomainCheckToWarn_type_2" &
                            , domainLowerLimitVec = domainLowerLimitVec &
                            , domainUpperLimitVec = domainUpperLimitVec &
                            , inputFile = "&ParaDRAM maxNumDomainCheckToWarn = 0 /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_MaxNumDomainCheckToWarn_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler quits with an error message when `maxNumDomainCheckToStop < 1`.
    function test_SpecBase_MaxNumDomainCheckToStop_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        real(RK), parameter :: domainLowerLimitVec(*) = [-1.e-1_RK]
        real(RK), parameter :: domainUpperLimitVec(*) = [+1.e+1_RK]
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_MaxNumDomainCheckToStop_type_1" &
                            , domainLowerLimitVec = domainLowerLimitVec &
                            , domainUpperLimitVec = domainUpperLimitVec &
                            , maxNumDomainCheckToStop = 0_IK &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_MaxNumDomainCheckToStop_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler quits with an error message when `maxNumDomainCheckToStop` has reached.
    function test_SpecBase_MaxNumDomainCheckToStop_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        real(RK), parameter :: domainLowerLimitVec(*) = [-1.e-1_RK]
        real(RK), parameter :: domainUpperLimitVec(*) = [+1.e+1_RK]
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_MaxNumDomainCheckToStop_type_2" &
                            , domainLowerLimitVec = domainLowerLimitVec &
                            , domainUpperLimitVec = domainUpperLimitVec &
                            , inputFile = "&ParaDRAM  maxNumDomainCheckToStop = 0 /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_MaxNumDomainCheckToStop_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler quits with an error message when `maxNumDomainCheckToStop` has reached.
    function test_SpecBase_MaxNumDomainCheckToStop_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        real(RK), parameter :: domainLowerLimitVec(*) = [-1.e-1_RK]
        real(RK), parameter :: domainUpperLimitVec(*) = [+1.e+1_RK]
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_MaxNumDomainCheckToStop_type_3" &
                            , domainLowerLimitVec = domainLowerLimitVec &
                            , domainUpperLimitVec = domainUpperLimitVec &
                            , maxNumDomainCheckToStop = 1_IK &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_MaxNumDomainCheckToStop_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler quits with an error message when `outputColumnWidth < 0`.
    function test_SpecBase_OutputColumnWidth_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_OutputColumnWidth_type_1" &
                            , outputColumnWidth = -1_IK &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_OutputColumnWidth_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler quits with an error message when `outputColumnWidth < 0`.
    function test_SpecBase_OutputColumnWidth_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_OutputColumnWidth_type_2" &
                            , inputFile = "&ParaDRAM outputColumnWidth = -1 /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_OutputColumnWidth_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler quits with an error message when `outputColumnWidth < outputRealPrecision + 7`.
    function test_SpecBase_OutputColumnWidth_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_OutputColumnWidth_type_3" &
                            , outputRealPrecision = 8_IK &
                            , outputColumnWidth = 14_IK &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_OutputColumnWidth_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler quits with an error message when `outputDelimiter` contains `.`.
    function test_SpecBase_OutputDelimiter_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_OutputDelimiter_type_1" &
                            , outputDelimiter = "this.that" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_OutputDelimiter_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler quits with an error message when `outputDelimiter` contains `+`.
    function test_SpecBase_OutputDelimiter_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_OutputDelimiter_type_2" &
                            , inputFile = "&ParaDRAM outputDelimiter = 'this+that' /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_OutputDelimiter_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler quits with an error message when `outputDelimiter` contains digit.
    function test_SpecBase_OutputDelimiter_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_OutputDelimiter_type_3" &
                            , inputFile = "&ParaDRAM outputDelimiter = 'this1234that' /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_OutputDelimiter_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler converts `\t` to tab character when `outputDelimiter = "\t"`.
    function test_SpecBase_OutputDelimiter_type_4() result(assertion)
        use Constants_mod, only: IK, RK, TAB
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_OutputDelimiter_type_4" &
                            , outputDelimiter = "\t" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. PD%SpecBase%OutputDelimiter%val == TAB
#endif
    end function test_SpecBase_OutputDelimiter_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler converts `\\t` to `\t` when `outputDelimiter = "\\t"`.
    function test_SpecBase_OutputDelimiter_type_5() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_OutputDelimiter_type_5" &
                            , outputDelimiter = "\\t" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. PD%SpecBase%OutputDelimiter%val == "\t"
#endif
    end function test_SpecBase_OutputDelimiter_type_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler converts `""` to `" "` when `outputDelimiter = ""` contains digit.
    function test_SpecBase_OutputDelimiter_type_6() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_OutputDelimiter_type_6" &
                            , outputDelimiter = "" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. PD%SpecBase%OutputDelimiter%val == " " .and. len(PD%SpecBase%OutputDelimiter%val) == 1
#endif
    end function test_SpecBase_OutputDelimiter_type_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler quits with an error message when `outputRealPrecision < 1`.
    function test_SpecBase_OutputRealPrecision_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_OutputRealPrecision_type_1" &
                            , outputRealPrecision = 0_IK &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_OutputRealPrecision_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler quits with an error message when `outputRealPrecision < 1`.
    function test_SpecBase_OutputRealPrecision_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_OutputRealPrecision_type_1" &
                            , inputFile = "&ParaDRAM outputRealPrecision = 0 /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_OutputRealPrecision_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler overwrites existing simulation files when `overwriteRequested = .true.`.
    function test_SpecBase_OverwriteRequested_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_OverwriteRequested_type_1" &
                            , overwriteRequested = .true. &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_OverwriteRequested_type_1" &
                            , inputFile = "&ParaDRAM overwriteRequested = true /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecBase_OverwriteRequested_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler exits with an error message when attempting to overwrite an existing simulation with `overwriteRequested = .false.`.
    function test_SpecBase_OverwriteRequested_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_OverwriteRequested_type_2" &
                            , overwriteRequested = .false. &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_OverwriteRequested_type_2" &
                            , inputFile = "&ParaDRAM overwriteRequested = false /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_OverwriteRequested_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler handles properly the `singleChain` parallelization model.
    function test_SpecBase_ParallelizationModel_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_ParallelizationModel_type_1" &
                            , parallelizationModel = "SinGleChAin" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecBase_ParallelizationModel_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler handles properly the `multiChain` parallelization model.
    function test_SpecBase_ParallelizationModel_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_ParallelizationModel_type_2" &
                            , inputFile = "&ParaDRAM parallelizationModel = 'MULTI ChAin' /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecBase_ParallelizationModel_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler exits with an error message when the specified value for `parallelizationModel` is not recognized.
    function test_SpecBase_ParallelizationModel_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_ParallelizationModel_type_3" &
                            , inputFile = "&ParaDRAM parallelizationModel = 'nonsense' /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_ParallelizationModel_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler exits with an error message when the `progressReportPeriod < 1`.
    function test_SpecBase_ProgressReportPeriod_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_ProgressReportPeriod_type_1" &
                            , progressReportPeriod = -1_IK &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_ProgressReportPeriod_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler exits with an error message when the `progressReportPeriod < 1`.
    function test_SpecBase_ProgressReportPeriod_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_ProgressReportPeriod_type_2" &
                            , inputFile = "&ParaDRAM progressReportPeriod = 0 /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_ProgressReportPeriod_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler appropriately considers the user's input random seed value.
    function test_SpecBase_RandomSeed_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_RandomSeed_type_1" &
                            , randomSeed = -12345_IK &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecBase_RandomSeed_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler appropriately considers the user's input random seed value.
    function test_SpecBase_RandomSeed_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_RandomSeed_type_2" &
                            , inputFile = "&ParaDRAM randomSeed = 12345 /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecBase_RandomSeed_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler appropriately considers a binary restart file format.
    function test_SpecBase_RestartFileFormat_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_RestartFileFormat_type_1" &
                            , restartFileFormat = "BINARY" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecBase_RestartFileFormat_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler appropriately considers an ascii restart file format.
    function test_SpecBase_RestartFileFormat_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_RestartFileFormat_type_2" &
                            , inputFile = "&ParaDRAM restartFileFormat = 'asCII' /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecBase_RestartFileFormat_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler exits with an error message when the specified restart file format is not recognized.
    function test_SpecBase_RestartFileFormat_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_RestartFileFormat_type_3" &
                            , inputFile = "&ParaDRAM restartFileFormat = 'nonsense' /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_RestartFileFormat_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler exits with an error message when the specified `targetAcceptanceRate` is out of range.
    function test_SpecBase_TargetAcceptanceRate_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        real(RK), parameter :: targetAcceptanceRate(*) = [-1._RK, 0.5_RK]
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_TargetAcceptanceRate_type_1" &
                            , targetAcceptanceRate = targetAcceptanceRate &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_TargetAcceptanceRate_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler exits with an error message when the specified `targetAcceptanceRate` is out of range.
    function test_SpecBase_TargetAcceptanceRate_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_TargetAcceptanceRate_type_2" &
                            , inputFile = "&ParaDRAM targetAcceptanceRate = +2. /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_TargetAcceptanceRate_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler exits with an error message when the specified `targetAcceptanceRate` is out of range.
    function test_SpecBase_TargetAcceptanceRate_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_TargetAcceptanceRate_type_3" &
                            , inputFile = "&ParaDRAM targetAcceptanceRate = +0.5, 0.2 /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_TargetAcceptanceRate_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler exits with an error message when the specified `targetAcceptanceRate` is out of range.
    function test_SpecBase_TargetAcceptanceRate_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_TargetAcceptanceRate_type_4" &
                            , inputFile = "&ParaDRAM targetAcceptanceRate(1) = 1. /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_TargetAcceptanceRate_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler exits with an error message when the specified `targetAcceptanceRate` is out of range.
    function test_SpecBase_TargetAcceptanceRate_type_5() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_TargetAcceptanceRate_type_5" &
                            , inputFile = "&ParaDRAM targetAcceptanceRate(2) = 0. /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecBase_TargetAcceptanceRate_type_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler does not do scaling when the `targetAcceptanceRate = [0., 1.]`.
    function test_SpecBase_TargetAcceptanceRate_type_6() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        real(RK), parameter :: targetAcceptanceRate(*) = [0._RK, 1._RK]
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = "test_SpecBase_TargetAcceptanceRate_type_6" &
                            , targetAcceptanceRate = targetAcceptanceRate &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. .not. PD%SpecBase%TargetAcceptanceRate%scalingRequested
#endif
    end function test_SpecBase_TargetAcceptanceRate_type_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a wrong internal input namelist group:
    !> +    Infinity values for `domainLowerLimitVec` and `domainUpperLimitVec`.
    function test_runSampler_7() result(assertion)
        use Constants_mod, only: IK, RK
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
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_runSampler_7

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a wrong internal input namelist group:
    !> +    Infinity values for `domainLowerLimitVec` and `domainUpperLimitVec`.
    function test_runSampler_8() result(assertion)
        use Constants_mod, only: IK, RK
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
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_runSampler_8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%