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
    !> Test the ParaDRAM sampler with a wrong input value for `chainSize < ndim + 1`.
    function test_SpecMCMC_ChainSize_type_1() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ChainSize_type_1" &
                            , chainSize = 0_IK &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ChainSize_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a wrong input value for `chainSize < ndim + 1`.
    function test_SpecMCMC_ChainSize_type_2() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ChainSize_type_2" &
                            , inputFile = "&ParaDRAM chainSize = 1 /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ChainSize_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a uniform proposal model.
    function test_SpecMCMC_ProposalModel_type_1() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalModel_type_1" &
                            , proposalModel = "UniForm" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalModel_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a uniform proposal model.
    function test_SpecMCMC_ProposalModel_type_2() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalModel_type_2" &
                            , inputFile = "&ParaDRAM proposalModel = 'UNIFORM' /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalModel_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a uniform proposal model.
    function test_SpecMCMC_ProposalModel_type_3() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalModel_type_3" &
                            , proposalModel = "Normal" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalModel_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a wrong unrecognized proposal model.
    function test_SpecMCMC_ProposalModel_type_4() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDRAM_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalModel_type_4" &
                            , inputFile = "&ParaDRAM proposalModel = 'nonsense' /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalModel_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a valid unidimensional `ProposalStartCorMat`.
    function test_SpecMCMC_ProposalStartCorMat_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 1_IK
        real(RK)    , parameter :: ProposalStartCorMat(NDIM,NDIM) = reshape([1._RK], shape = shape(ProposalStartCorMat))
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartCorMat_type_1" &
                            , ProposalStartCorMat = ProposalStartCorMat &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalStartCorMat_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a valid unidimensional `ProposalStartCorMat`.
    function test_SpecMCMC_ProposalStartCorMat_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 1_IK
        real(RK)    , parameter :: ProposalStartCorMat(NDIM,NDIM) = reshape([1._RK], shape = shape(ProposalStartCorMat))
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartCorMat_type_2" &
                            , inputFile = "&ParaDRAM ProposalStartCorMat = "//num2str(ProposalStartCorMat)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalStartCorMat_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a wrong value for the unidimensional `ProposalStartCorMat`.
    function test_SpecMCMC_ProposalStartCorMat_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 1_IK
        real(RK)    , parameter :: ProposalStartCorMat(NDIM,NDIM) = reshape([0._RK], shape = shape(ProposalStartCorMat))
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartCorMat_type_3" &
                            , inputFile = "&ParaDRAM ProposalStartCorMat = "//num2str(ProposalStartCorMat)//" /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalStartCorMat_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a multidimensional `ProposalStartCorMat`.
    function test_SpecMCMC_ProposalStartCorMat_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: ProposalStartCorMat(NDIM,NDIM) = reshape([1._RK, 0.5_RK, 0.5_RK, 1._RK], shape = shape(ProposalStartCorMat))
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartCorMat_type_4" &
                            , inputFile = "&ParaDRAM ProposalStartCorMat = "//num2str(ProposalStartCorMat)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%SpecMCMC%ProposalStartCorMat%Val==PD%SpecMCMC%ProposalStartCovMat%Val)
#endif
    end function test_SpecMCMC_ProposalStartCorMat_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a wrong multidimensional `ProposalStartCorMat`, 
    !> which must be fine with the sampler, as long as it leads to a correct positive-definite covariance matrix.
    function test_SpecMCMC_ProposalStartCorMat_type_5() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: ProposalStartCorMat(NDIM,NDIM) = reshape([2._RK, 0.5_RK, 0.5_RK, 2._RK], shape = shape(ProposalStartCorMat))
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartCorMat_type_5" &
                            , inputFile = "&ParaDRAM ProposalStartCorMat = "//num2str(ProposalStartCorMat)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalStartCorMat_type_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a valid unidimensional `ProposalStartCorMat`.
    function test_SpecMCMC_ProposalStartCovMat_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 1_IK
        real(RK)    , parameter :: ProposalStartCovMat(NDIM,NDIM) = reshape([2._RK], shape = shape(ProposalStartCovMat))
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartCovMat_type_1" &
                            , ProposalStartCovMat = ProposalStartCovMat &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalStartCovMat_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a valid unidimensional `ProposalStartCovMat`.
    function test_SpecMCMC_ProposalStartCovMat_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 1_IK
        real(RK)    , parameter :: ProposalStartCovMat(NDIM,NDIM) = reshape([1._RK], shape = shape(ProposalStartCovMat))
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartCovMat_type_2" &
                            , inputFile = "&ParaDRAM ProposalStartCovMat = "//num2str(ProposalStartCovMat)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalStartCovMat_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a wrong value for the unidimensional `ProposalStartCovMat`.
    function test_SpecMCMC_ProposalStartCovMat_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 1_IK
        real(RK)    , parameter :: ProposalStartCovMat(NDIM,NDIM) = reshape([0._RK], shape = shape(ProposalStartCovMat))
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartCovMat_type_3" &
                            , inputFile = "&ParaDRAM ProposalStartCovMat = "//num2str(ProposalStartCovMat)//" /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalStartCovMat_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a multidimensional `ProposalStartCovMat`.
    function test_SpecMCMC_ProposalStartCovMat_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: ProposalStartCovMat(NDIM,NDIM) = reshape([1._RK, 0.5_RK, 0.5_RK, 1._RK], shape = shape(ProposalStartCovMat))
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartCovMat_type_4" &
                            , inputFile = "&ParaDRAM ProposalStartCovMat = "//num2str(ProposalStartCovMat)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%SpecMCMC%ProposalStartCovMat%Val==PD%SpecMCMC%ProposalStartCovMat%Val)
#endif
    end function test_SpecMCMC_ProposalStartCovMat_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a multidimensional `ProposalStartCovMat`, in the presence of `ProposalStartCorMat` 
    !> and `ProposalStartStdVec`, in which case, the `ProposalStartCovMat` must be preferred.
    function test_SpecMCMC_ProposalStartCovMat_type_5() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: ProposalStartCovMat(NDIM,NDIM) = reshape([2._RK, 0.5_RK, 0.5_RK, 2._RK], shape = shape(ProposalStartCovMat))
        real(RK)    , parameter :: ProposalStartCorMat(NDIM,NDIM) = reshape([1._RK, 0.0_RK, 0.0_RK, 1._RK], shape = shape(ProposalStartCovMat))
        real(RK)    , parameter :: ProposalStartStdVec(NDIM) = [1._RK, 1._RK]
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartCovMat_type_5" &
                            , ProposalStartCovMat = ProposalStartCovMat &
                            , ProposalStartCorMat = ProposalStartCorMat &
                            , ProposalStartStdVec = ProposalStartStdVec &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalStartCovMat_type_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a valid unidimensional `ProposalStartCorMat`.
    function test_SpecMCMC_ProposalStartStdVec_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 1_IK
        real(RK)    , parameter :: ProposalStartStdVec(NDIM) = [2._RK]
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartStdVec_type_1" &
                            , ProposalStartStdVec = ProposalStartStdVec &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalStartStdVec_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a valid unidimensional `ProposalStartStdVec`.
    function test_SpecMCMC_ProposalStartStdVec_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 1_IK
        real(RK)    , parameter :: ProposalStartStdVec(NDIM) = [1._RK]
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartStdVec_type_2" &
                            , inputFile = "&ParaDRAM ProposalStartStdVec = "//num2str(ProposalStartStdVec)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalStartStdVec_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a wrong value for the unidimensional `ProposalStartStdVec`.
    function test_SpecMCMC_ProposalStartStdVec_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 1_IK
        real(RK)    , parameter :: ProposalStartStdVec(NDIM) = [0._RK]
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartStdVec_type_3" &
                            , inputFile = "&ParaDRAM ProposalStartStdVec = "//num2str(ProposalStartStdVec)//" /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalStartStdVec_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a multidimensional `ProposalStartStdVec`.
    function test_SpecMCMC_ProposalStartStdVec_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: ProposalStartStdVec(NDIM) = reshape([1._RK, 0.5_RK], shape = shape(ProposalStartStdVec))
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartStdVec_type_4" &
                            , inputFile = "&ParaDRAM ProposalStartStdVec = "//num2str(ProposalStartStdVec)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%SpecMCMC%ProposalStartStdVec%Val==PD%SpecMCMC%ProposalStartStdVec%Val)
#endif
    end function test_SpecMCMC_ProposalStartStdVec_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDRAM sampler with a multidimensional `ProposalStartStdVec`, in the presence of `ProposalStartCorMat`,
    !> in which case, the `ProposalStartStdVec` must be correctly computed.
    function test_SpecMCMC_ProposalStartStdVec_type_5() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: ProposalStartStdVec(NDIM) = reshape([1._RK, 2._RK], shape = shape(ProposalStartStdVec))
        real(RK)    , parameter :: ProposalStartCorMat(NDIM,NDIM) = reshape([1._RK, 0.5_RK, 0.5_RK, 1._RK], shape = shape(ProposalStartCorMat))
        real(RK)    , parameter :: ProposalStartCovMat(NDIM,NDIM) = reshape([1._RK, 1.0_RK, 1.0_RK, 4._RK], shape = shape(ProposalStartCovMat))
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartStdVec_type_5" &
                            , ProposalStartStdVec = ProposalStartStdVec &
                            , ProposalStartCorMat = ProposalStartCorMat &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(abs(PD%SpecMCMC%ProposalStartCovMat%Val-ProposalStartCovMat) <= tolerance)
#endif
    end function test_SpecMCMC_ProposalStartStdVec_type_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
