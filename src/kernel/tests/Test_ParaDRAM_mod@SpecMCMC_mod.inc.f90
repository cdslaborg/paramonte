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

    !> \brief
    !> Test whether the ParaDRAM sampler can correctly set the value of `RandomStartPointDomainLowerLimitVec` from input argument.
    function test_SpecMCMC_RandomStartPointDomainLowerLimitVec_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: RandomStartPointDomainLowerLimitVec(NDIM) = [1._RK, 2._RK]
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_RandomStartPointDomainLowerLimitVec_type_1" &
                            , RandomStartPointDomainLowerLimitVec = RandomStartPointDomainLowerLimitVec &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%SpecMCMC%RandomStartPointDomainLowerLimitVec%Val == RandomStartPointDomainLowerLimitVec)
#endif
    end function test_SpecMCMC_RandomStartPointDomainLowerLimitVec_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler can correctly set the value of
    !> `RandomStartPointDomainLowerLimitVec` for two consecutive simulations.
    function test_SpecMCMC_RandomStartPointDomainLowerLimitVec_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: RandomStartPointDomainLowerLimitVec(NDIM) = [1._RK, 2._RK]
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_RandomStartPointDomainLowerLimitVec_type_2A" &
                            , inputFile = "&ParaDRAM RandomStartPointDomainLowerLimitVec = "//num2str(-RandomStartPointDomainLowerLimitVec)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%SpecMCMC%RandomStartPointDomainLowerLimitVec%Val == -RandomStartPointDomainLowerLimitVec)
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_RandomStartPointDomainLowerLimitVec_type_2B" &
                            , inputFile = "&ParaDRAM RandomStartPointDomainLowerLimitVec = "//num2str(RandomStartPointDomainLowerLimitVec)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%SpecMCMC%RandomStartPointDomainLowerLimitVec%Val == RandomStartPointDomainLowerLimitVec)
#endif
    end function test_SpecMCMC_RandomStartPointDomainLowerLimitVec_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler returns with an error message when
    !> `RandomStartPointDomainLowerLimitVec` goes below the limits of `DomainLowerLimitVec`.
    function test_SpecMCMC_RandomStartPointDomainLowerLimitVec_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: RandomStartPointDomainLowerLimitVec(NDIM) = [+1._RK, +2._RK]
        real(RK)    , parameter :: DomainLowerLimitVec(NDIM) = [-1._RK, +3._RK]
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_RandomStartPointDomainLowerLimitVec_type_3" &
                            , RandomStartPointDomainLowerLimitVec = RandomStartPointDomainLowerLimitVec &
                            , DomainLowerLimitVec = DomainLowerLimitVec &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_RandomStartPointDomainLowerLimitVec_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler can correctly set the value of `RandomStartPointDomainLowerLimitVec` from input argument.
    function test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: RandomStartPointDomainUpperLimitVec(NDIM) = [1._RK, 2._RK]
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_1" &
                            , RandomStartPointDomainUpperLimitVec = RandomStartPointDomainUpperLimitVec &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%SpecMCMC%RandomStartPointDomainUpperLimitVec%Val == RandomStartPointDomainUpperLimitVec)
#endif
    end function test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler can correctly set the value of
    !> `RandomStartPointDomainUpperLimitVec` for two consecutive simulations.
    function test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: RandomStartPointDomainUpperLimitVec(NDIM) = [1._RK, 2._RK]
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_2A" &
                            , inputFile = "&ParaDRAM RandomStartPointDomainUpperLimitVec = "//num2str(-RandomStartPointDomainUpperLimitVec)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%SpecMCMC%RandomStartPointDomainUpperLimitVec%Val == -RandomStartPointDomainUpperLimitVec)
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_2B" &
                            , inputFile = "&ParaDRAM RandomStartPointDomainUpperLimitVec = "//num2str(RandomStartPointDomainUpperLimitVec)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%SpecMCMC%RandomStartPointDomainUpperLimitVec%Val == RandomStartPointDomainUpperLimitVec)
#endif
    end function test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler returns with an error message when
    !> `RandomStartPointDomainUpperLimitVec <= RandomStartPointDomainLowerLimitVec`.
    function test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: RandomStartPointDomainLowerLimitVec(NDIM) = [1._RK, 2._RK]
        real(RK)    , parameter :: RandomStartPointDomainUpperLimitVec(NDIM) = [-1._RK, -2._RK]
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_3" &
                            , RandomStartPointDomainLowerLimitVec = RandomStartPointDomainLowerLimitVec &
                            , RandomStartPointDomainUpperLimitVec = RandomStartPointDomainUpperLimitVec &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler returns with an error message when
    !> `RandomStartPointDomainUpperLimitVec` goes beyond the limits of `DomainUpperLimitVec`.
    function test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: RandomStartPointDomainUpperLimitVec(NDIM) = [+1._RK, +2._RK]
        real(RK)    , parameter :: DomainUpperLimitVec(NDIM) = [-1._RK, +3._RK]
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_4" &
                            , RandomStartPointDomainUpperLimitVec = RandomStartPointDomainUpperLimitVec &
                            , DomainUpperLimitVec = DomainUpperLimitVec &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_RandomStartPointDomainUpperLimitVec_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler correctly randomizes the start point when requested.
    function test_SpecMCMC_RandomStartPointRequested_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_RandomStartPointRequested_type_1" &
                            , randomStartPointRequested = .false. &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. .not. PD%SpecMCMC%RandomStartPointRequested%val
#endif
    end function test_SpecMCMC_RandomStartPointRequested_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test `randomStartPointRequested = true`, the use must also specify domain of either the target or the random start point.
    function test_SpecMCMC_RandomStartPointRequested_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_RandomStartPointRequested_type_2" &
                            , inputFile = "&paradram randomStartPointRequested = true /" &
                            )
        assertion = assertion .and. PD%Err%occurred .and. PD%SpecMCMC%RandomStartPointRequested%val
#endif
    end function test_SpecMCMC_RandomStartPointRequested_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test `randomStartPointRequested = true`, the use must also specify domain of either the target or the random start point.
    function test_SpecMCMC_RandomStartPointRequested_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: RandomStartPointDomainLowerLimitVec(NDIM) = [-1.e0_RK, +1.e1_RK]
        real(RK)    , parameter :: RandomStartPointDomainUpperLimitVec(NDIM) = [+2.e0_RK, +2.e1_RK]
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_RandomStartPointRequested_type_3" &
                            , inputFile = "&paradram randomStartPointRequested = true /" &
                            , RandomStartPointDomainLowerLimitVec = RandomStartPointDomainLowerLimitVec &
                            , RandomStartPointDomainUpperLimitVec = RandomStartPointDomainUpperLimitVec &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. PD%SpecMCMC%RandomStartPointRequested%val
#endif
    end function test_SpecMCMC_RandomStartPointRequested_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test `randomStartPointRequested = true`, the use must also specify domain of either the target or the random start point.
    function test_SpecMCMC_RandomStartPointRequested_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: DomainLowerLimitVec(NDIM) = [-1.e0_RK, +1.e2_RK]
        real(RK)    , parameter :: DomainUpperLimitVec(NDIM) = [+2.e0_RK, +2.e2_RK]
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_RandomStartPointRequested_type_4" &
                            , inputFile = "&paradram randomStartPointRequested = true /" &
                            , DomainLowerLimitVec = DomainLowerLimitVec &
                            , DomainUpperLimitVec = DomainUpperLimitVec &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. PD%SpecMCMC%RandomStartPointRequested%val
#endif
    end function test_SpecMCMC_RandomStartPointRequested_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler returns an error message when `sampleRefinementCount < 0`.
    function test_SpecMCMC_SampleRefinementCount_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementCount_type_1" &
                            , sampleRefinementCount = -1_IK &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementCount_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler returns an error message when `sampleRefinementCount < 0`.
    function test_SpecMCMC_SampleRefinementCount_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementCount_type_2" &
                            , inputFile = "&ParaDRAM sampleRefinementCount = -1 /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementCount_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler succeeds with a valid `sampleRefinementCount >= 0`.
    function test_SpecMCMC_SampleRefinementCount_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementCount_type_3" &
                            , inputFile = "&ParaDRAM sampleRefinementCount = 0 /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementCount_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler returns an error message with an unrecognized value for `sampleRefinementMethod`.
    function test_SpecMCMC_SampleRefinementMethod_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementMethod_type_1" &
                            , inputFile = "&ParaDRAM sampleRefinementMethod = 'nonsense' /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementMethod_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler succeeds with a valid value for `sampleRefinementMethod`.
    function test_SpecMCMC_SampleRefinementMethod_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementMethod_type_2" &
                            , sampleRefinementMethod = "Batch  Means" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementMethod_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler succeeds with a valid value for `sampleRefinementMethod`.
    function test_SpecMCMC_SampleRefinementMethod_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementMethod_type_3" &
                            , inputFile = "&ParaDRAM sampleRefinementMethod = 'CutOffAutoCorr' /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementMethod_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler succeeds with a valid value for `sampleRefinementMethod`.
    function test_SpecMCMC_SampleRefinementMethod_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementMethod_type_4" &
                            , inputFile = "&ParaDRAM sampleRefinementMethod = 'CutOf  fAuto Corr' /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementMethod_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler succeeds with a valid value for `sampleRefinementMethod`.
    function test_SpecMCMC_SampleRefinementMethod_type_5() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementMethod_type_5" &
                            , inputFile = "&ParaDRAM sampleRefinementMethod = 'MaxCum SumAut oCorr' /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementMethod_type_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler returns with an error message for a wrong input value for `scaleFactor`.
    function test_SpecMCMC_ScaleFactor_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ScaleFactor_type_1" &
                            , inputFile = "&ParaDRAM scaleFactor = 'nonsense' /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ScaleFactor_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler returns with an error message for an empty input value for `scaleFactor`.
    function test_SpecMCMC_ScaleFactor_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ScaleFactor_type_2" &
                            , scaleFactor = " " &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ScaleFactor_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler returns with an error message for an empty input value for `scaleFactor`.
    function test_SpecMCMC_ScaleFactor_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ScaleFactor_type_3" &
                            , scaleFactor = " " &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ScaleFactor_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler returns with an error message for a wrong input value for `scaleFactor`.
    function test_SpecMCMC_ScaleFactor_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ScaleFactor_type_4" &
                            , scaleFactor = "Gelman / 2" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ScaleFactor_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler returns successfully for a valid input value for `scaleFactor`.
    function test_SpecMCMC_ScaleFactor_type_5() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ScaleFactor_type_5" &
                            , inputFile = "&ParaDRAM scaleFactor = '2 * Gelman'" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_ScaleFactor_type_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler returns with an error message for a negative input value for `scaleFactor`.
    function test_SpecMCMC_ScaleFactor_type_6() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ScaleFactor_type_6" &
                            , inputFile = "&ParaDRAM scaleFactor = '-0.5 * Gelman'" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ScaleFactor_type_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler returns with an error message for a wrong input value for `scaleFactor`.
    function test_SpecMCMC_ScaleFactor_type_7() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ScaleFactor_type_7" &
                            , inputFile = "&ParaDRAM scaleFactor = '*GELMAN'" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ScaleFactor_type_7

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler returns successfully for a valid input value for `scaleFactor`.
    function test_SpecMCMC_ScaleFactor_type_8() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ScaleFactor_type_8" &
                            , inputFile = "&ParaDRAM scaleFactor = '2. * 0.5 * Gelman'" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. abs(PD%SpecMCMC%ScaleFactor%val-2.38_RK/sqrt(real(NDIM,kind=RK)))<1.e-12_RK
#endif
    end function test_SpecMCMC_ScaleFactor_type_8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler returns successfully for a valid input value for `StartPointVec`.
    function test_SpecMCMC_StartPointVec_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: StartPointVec(NDIM) = [1._RK, -10._RK]
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_StartPointVec_type_1" &
                            , StartPointVec = StartPointVec &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. abs(PD%SpecMCMC%ScaleFactor%val-2.38_RK/sqrt(real(NDIM,kind=RK)))<1.e-12_RK
#endif
    end function test_SpecMCMC_StartPointVec_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler returns successfully for a valid input value for `StartPointVec`.
    function test_SpecMCMC_StartPointVec_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_StartPointVec_type_2" &
                            , inputFile = "&ParaDRAM StartPointVec = 1., -10., /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_StartPointVec_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler returns successfully and correctly computes `StartPointVec` when not provided.
    function test_SpecMCMC_StartPointVec_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_StartPointVec_type_3" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all( abs(PD%SpecMCMC%StartPointVec%Val) < 1.e-12_RK )
#endif
    end function test_SpecMCMC_StartPointVec_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDRAM sampler returns successfully and correctly computes `StartPointVec` when not provided.
    function test_SpecMCMC_StartPointVec_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDRAM_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: DomainLowerLimitVec(NDIM) = [-1.e0_RK, +1.e1_RK]
        real(RK)    , parameter :: DomainUpperLimitVec(NDIM) = [+2.e0_RK, +2.e1_RK]
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_StartPointVec_type_4" &
                            , DomainLowerLimitVec = DomainLowerLimitVec &
                            , DomainUpperLimitVec = DomainUpperLimitVec &
                            , chainSize = 100_IK &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all( abs(PD%SpecMCMC%StartPointVec%Val-0.5_RK*(DomainLowerLimitVec+DomainUpperLimitVec)) < 1.e-12_RK )
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "PD%Err%occurred", PD%Err%occurred
            write(Test%outputUnit,"(*(g0,:,' '))") "PD%SpecMCMC%StartPointVec%Val", PD%SpecMCMC%StartPointVec%Val
            write(Test%outputUnit,"(*(g0,:,' '))") "DomainUpperLimitVec", DomainUpperLimitVec
            write(Test%outputUnit,"(*(g0,:,' '))") "DomainLowerLimitVec", DomainLowerLimitVec
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
#endif
    end function test_SpecMCMC_StartPointVec_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%