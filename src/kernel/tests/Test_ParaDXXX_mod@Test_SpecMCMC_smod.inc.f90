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

!>  \brief This include file contains the body of the submodules 
!>  [ParaDRAM_mod@Test_SpecMCMC_smod](@ref paradram_mod@test_specmcmc_smod) and 
!>  [ParaDISE_mod@@Test_SpecMCMC_smod](@ref paradise_mod@@test_specmcmc_smod).
!>  \author Amir Shahmoradi

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with a wrong input value for `chainSize < ndim + 1`.
    module function test_SpecMCMC_ChainSize_type_1() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
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
    !> Test the ParaDXXX sampler with a wrong input value for `chainSize < ndim + 1`.
    module function test_SpecMCMC_ChainSize_type_2() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ChainSize_type_2" &
                            , inputFile = ParaDXXX_NML//" chainSize = 1 /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ChainSize_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with a uniform proposal model.
    module function test_SpecMCMC_ProposalModel_type_1() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
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
    !> Test the ParaDXXX sampler with a uniform proposal model.
    module function test_SpecMCMC_ProposalModel_type_2() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalModel_type_2" &
                            , inputFile = ParaDXXX_NML//" proposalModel = 'UNIFORM' /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalModel_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with a uniform proposal model.
    module function test_SpecMCMC_ProposalModel_type_3() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
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
    !> Test the ParaDXXX sampler with a wrong unrecognized proposal model.
    module function test_SpecMCMC_ProposalModel_type_4() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaDXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalModel_type_4" &
                            , inputFile = ParaDXXX_NML//" proposalModel = 'nonsense' /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalModel_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with a valid unidimensional `ProposalStartCorMat`.
    module function test_SpecMCMC_ProposalStartCorMat_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 1_IK
        real(RK)    , parameter :: ProposalStartCorMat(NDIM,NDIM) = reshape([1._RK], shape = shape(ProposalStartCorMat))
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
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
    !> Test the ParaDXXX sampler with a valid unidimensional `ProposalStartCorMat`.
    module function test_SpecMCMC_ProposalStartCorMat_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 1_IK
        real(RK)    , parameter :: ProposalStartCorMat(NDIM,NDIM) = reshape([1._RK], shape = shape(ProposalStartCorMat))
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartCorMat_type_2" &
                            , inputFile = ParaDXXX_NML//" ProposalStartCorMat = "//num2str(ProposalStartCorMat)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalStartCorMat_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with a wrong value for the unidimensional `ProposalStartCorMat`.
    module function test_SpecMCMC_ProposalStartCorMat_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 1_IK
        real(RK)    , parameter :: ProposalStartCorMat(NDIM,NDIM) = reshape([0._RK], shape = shape(ProposalStartCorMat))
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartCorMat_type_3" &
                            , inputFile = ParaDXXX_NML//" ProposalStartCorMat = "//num2str(ProposalStartCorMat)//" /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalStartCorMat_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with a multidimensional `ProposalStartCorMat`.
    module function test_SpecMCMC_ProposalStartCorMat_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: ProposalStartCorMat(NDIM,NDIM) = reshape([1._RK, 0.5_RK, 0.5_RK, 1._RK], shape = shape(ProposalStartCorMat))
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartCorMat_type_4" &
                            , inputFile = ParaDXXX_NML//" ProposalStartCorMat = "//num2str(ProposalStartCorMat)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%SpecMCMC%ProposalStartCorMat%Val==PD%SpecMCMC%ProposalStartCovMat%Val)
#endif
    end function test_SpecMCMC_ProposalStartCorMat_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with a wrong multidimensional `ProposalStartCorMat`,
    !> which must be fine with the sampler, as long as it leads to a correct positive-definite covariance matrix.
    module function test_SpecMCMC_ProposalStartCorMat_type_5() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: ProposalStartCorMat(NDIM,NDIM) = reshape([2._RK, 0.5_RK, 0.5_RK, 2._RK], shape = shape(ProposalStartCorMat))
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartCorMat_type_5" &
                            , inputFile = ParaDXXX_NML//" ProposalStartCorMat = "//num2str(ProposalStartCorMat)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalStartCorMat_type_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with a valid unidimensional `ProposalStartCorMat`.
    module function test_SpecMCMC_ProposalStartCovMat_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 1_IK
        real(RK)    , parameter :: ProposalStartCovMat(NDIM,NDIM) = reshape([2._RK], shape = shape(ProposalStartCovMat))
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
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
    !> Test the ParaDXXX sampler with a valid unidimensional `ProposalStartCovMat`.
    module function test_SpecMCMC_ProposalStartCovMat_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 1_IK
        real(RK)    , parameter :: ProposalStartCovMat(NDIM,NDIM) = reshape([1._RK], shape = shape(ProposalStartCovMat))
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartCovMat_type_2" &
                            , inputFile = ParaDXXX_NML//" ProposalStartCovMat = "//num2str(ProposalStartCovMat)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalStartCovMat_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with a wrong value for the unidimensional `ProposalStartCovMat`.
    module function test_SpecMCMC_ProposalStartCovMat_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 1_IK
        real(RK)    , parameter :: ProposalStartCovMat(NDIM,NDIM) = reshape([0._RK], shape = shape(ProposalStartCovMat))
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartCovMat_type_3" &
                            , inputFile = ParaDXXX_NML//" ProposalStartCovMat = "//num2str(ProposalStartCovMat)//" /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalStartCovMat_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with a multidimensional `ProposalStartCovMat`.
    module function test_SpecMCMC_ProposalStartCovMat_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: ProposalStartCovMat(NDIM,NDIM) = reshape([1._RK, 0.5_RK, 0.5_RK, 1._RK], shape = shape(ProposalStartCovMat))
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartCovMat_type_4" &
                            , inputFile = ParaDXXX_NML//" ProposalStartCovMat = "//num2str(ProposalStartCovMat)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%SpecMCMC%ProposalStartCovMat%Val==PD%SpecMCMC%ProposalStartCovMat%Val)
#endif
    end function test_SpecMCMC_ProposalStartCovMat_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with a multidimensional `ProposalStartCovMat`, in the presence of `ProposalStartCorMat`
    !> and `ProposalStartStdVec`, in which case, the `ProposalStartCovMat` must be preferred.
    module function test_SpecMCMC_ProposalStartCovMat_type_5() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: ProposalStartCovMat(NDIM,NDIM) = reshape([2._RK, 0.5_RK, 0.5_RK, 2._RK], shape = shape(ProposalStartCovMat))
        real(RK)    , parameter :: ProposalStartCorMat(NDIM,NDIM) = reshape([1._RK, 0.0_RK, 0.0_RK, 1._RK], shape = shape(ProposalStartCovMat))
        real(RK)    , parameter :: ProposalStartStdVec(NDIM) = [1._RK, 1._RK]
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
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
    !> Test the ParaDXXX sampler with a valid unidimensional `ProposalStartCorMat`.
    module function test_SpecMCMC_ProposalStartStdVec_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 1_IK
        real(RK)    , parameter :: ProposalStartStdVec(NDIM) = [2._RK]
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
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
    !> Test the ParaDXXX sampler with a valid unidimensional `ProposalStartStdVec`.
    module function test_SpecMCMC_ProposalStartStdVec_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 1_IK
        real(RK)    , parameter :: ProposalStartStdVec(NDIM) = [1._RK]
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartStdVec_type_2" &
                            , inputFile = ParaDXXX_NML//" ProposalStartStdVec = "//num2str(ProposalStartStdVec)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalStartStdVec_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with a wrong value for the unidimensional `ProposalStartStdVec`.
    module function test_SpecMCMC_ProposalStartStdVec_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 1_IK
        real(RK)    , parameter :: ProposalStartStdVec(NDIM) = [0._RK]
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartStdVec_type_3" &
                            , inputFile = ParaDXXX_NML//" ProposalStartStdVec = "//num2str(ProposalStartStdVec)//" /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ProposalStartStdVec_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with a multidimensional `ProposalStartStdVec`.
    module function test_SpecMCMC_ProposalStartStdVec_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: ProposalStartStdVec(NDIM) = reshape([1._RK, 0.5_RK], shape = shape(ProposalStartStdVec))
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ProposalStartStdVec_type_4" &
                            , inputFile = ParaDXXX_NML//" ProposalStartStdVec = "//num2str(ProposalStartStdVec)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%SpecMCMC%ProposalStartStdVec%Val==PD%SpecMCMC%ProposalStartStdVec%Val)
#endif
    end function test_SpecMCMC_ProposalStartStdVec_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaDXXX sampler with a multidimensional `ProposalStartStdVec`, in the presence of `ProposalStartCorMat`,
    !> in which case, the `ProposalStartStdVec` must be correctly computed.
    module function test_SpecMCMC_ProposalStartStdVec_type_5() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: ProposalStartStdVec(NDIM) = reshape([1._RK, 2._RK], shape = shape(ProposalStartStdVec))
        real(RK)    , parameter :: ProposalStartCorMat(NDIM,NDIM) = reshape([1._RK, 0.5_RK, 0.5_RK, 1._RK], shape = shape(ProposalStartCorMat))
        real(RK)    , parameter :: ProposalStartCovMat(NDIM,NDIM) = reshape([1._RK, 1.0_RK, 1.0_RK, 4._RK], shape = shape(ProposalStartCovMat))
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
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
    !> Test whether the ParaDXXX sampler can correctly set the value of `RandomStartPointDomainLowerLimitVec` from input argument.
    module function test_RSPDLowerLimitVec_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: RandomStartPointDomainLowerLimitVec(NDIM) = [1._RK, 2._RK]
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_RSPDLowerLimitVec_type_1" &
                            , RandomStartPointDomainLowerLimitVec = RandomStartPointDomainLowerLimitVec &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%SpecMCMC%RandomStartPointDomainLowerLimitVec%Val == RandomStartPointDomainLowerLimitVec)
#endif
    end function test_RSPDLowerLimitVec_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler can correctly set the value of
    !> `RandomStartPointDomainLowerLimitVec` for two consecutive simulations.
    module function test_RSPDLowerLimitVec_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: RandomStartPointDomainLowerLimitVec(NDIM) = [1._RK, 2._RK]
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_RSPDLowerLimitVec_type_2A" &
                            , inputFile = ParaDXXX_NML//" RandomStartPointDomainLowerLimitVec = "//num2str(-RandomStartPointDomainLowerLimitVec)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%SpecMCMC%RandomStartPointDomainLowerLimitVec%Val == -RandomStartPointDomainLowerLimitVec)
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_RSPDLowerLimitVec_type_2B" &
                            , inputFile = ParaDXXX_NML//" RandomStartPointDomainLowerLimitVec = "//num2str(RandomStartPointDomainLowerLimitVec)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%SpecMCMC%RandomStartPointDomainLowerLimitVec%Val == RandomStartPointDomainLowerLimitVec)
#endif
    end function test_RSPDLowerLimitVec_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler returns with an error message when
    !> `RandomStartPointDomainLowerLimitVec` goes below the limits of `DomainLowerLimitVec`.
    module function test_RSPDLowerLimitVec_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: RandomStartPointDomainLowerLimitVec(NDIM) = [+1._RK, +2._RK]
        real(RK)    , parameter :: DomainLowerLimitVec(NDIM) = [-1._RK, +3._RK]
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_RSPDLowerLimitVec_type_3" &
                            , RandomStartPointDomainLowerLimitVec = RandomStartPointDomainLowerLimitVec &
                            , DomainLowerLimitVec = DomainLowerLimitVec &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_RSPDLowerLimitVec_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler can correctly set the value of `RandomStartPointDomainLowerLimitVec` from input argument.
    module function test_RSPDUpperLimitVec_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: RandomStartPointDomainUpperLimitVec(NDIM) = [1._RK, 2._RK]
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_RSPDUpperLimitVec_type_1" &
                            , RandomStartPointDomainUpperLimitVec = RandomStartPointDomainUpperLimitVec &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%SpecMCMC%RandomStartPointDomainUpperLimitVec%Val == RandomStartPointDomainUpperLimitVec)
#endif
    end function test_RSPDUpperLimitVec_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler can correctly set the value of
    !> `RandomStartPointDomainUpperLimitVec` for two consecutive simulations.
    module function test_RSPDUpperLimitVec_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: RandomStartPointDomainUpperLimitVec(NDIM) = [1._RK, 2._RK]
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_RSPDUpperLimitVec_type_2A" &
                            , inputFile = ParaDXXX_NML//" RandomStartPointDomainUpperLimitVec = "//num2str(-RandomStartPointDomainUpperLimitVec)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%SpecMCMC%RandomStartPointDomainUpperLimitVec%Val == -RandomStartPointDomainUpperLimitVec)
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_RSPDUpperLimitVec_type_2B" &
                            , inputFile = ParaDXXX_NML//" RandomStartPointDomainUpperLimitVec = "//num2str(RandomStartPointDomainUpperLimitVec)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%SpecMCMC%RandomStartPointDomainUpperLimitVec%Val == RandomStartPointDomainUpperLimitVec)
#endif
    end function test_RSPDUpperLimitVec_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler returns with an error message when
    !> `RandomStartPointDomainUpperLimitVec <= RandomStartPointDomainLowerLimitVec`.
    module function test_RSPDUpperLimitVec_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: RandomStartPointDomainLowerLimitVec(NDIM) = [1._RK, 2._RK]
        real(RK)    , parameter :: RandomStartPointDomainUpperLimitVec(NDIM) = [-1._RK, -2._RK]
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_RSPDUpperLimitVec_type_3" &
                            , RandomStartPointDomainLowerLimitVec = RandomStartPointDomainLowerLimitVec &
                            , RandomStartPointDomainUpperLimitVec = RandomStartPointDomainUpperLimitVec &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_RSPDUpperLimitVec_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler returns with an error message when
    !> `RandomStartPointDomainUpperLimitVec` goes beyond the limits of `DomainUpperLimitVec`.
    module function test_RSPDUpperLimitVec_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: RandomStartPointDomainUpperLimitVec(NDIM) = [+1._RK, +2._RK]
        real(RK)    , parameter :: DomainUpperLimitVec(NDIM) = [-1._RK, +3._RK]
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_RSPDUpperLimitVec_type_4" &
                            , RandomStartPointDomainUpperLimitVec = RandomStartPointDomainUpperLimitVec &
                            , DomainUpperLimitVec = DomainUpperLimitVec &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_RSPDUpperLimitVec_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler correctly randomizes the start point when requested.
    module function test_SpecMCMC_RandomStartPointRequested_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
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
    module function test_SpecMCMC_RandomStartPointRequested_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_RandomStartPointRequested_type_2" &
                            , inputFile = ParaDXXX_NML//" randomStartPointRequested = true /" &
                            )
        assertion = assertion .and. PD%Err%occurred .and. PD%SpecMCMC%RandomStartPointRequested%val
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "PD%Err%occurred :", PD%Err%occurred
            write(Test%outputUnit,"(*(g0,:,', '))") "PD%SpecMCMC%RandomStartPointRequested%val :", PD%SpecMCMC%RandomStartPointRequested%val
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
#endif
    end function test_SpecMCMC_RandomStartPointRequested_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> When `randomStartPointRequested = true`, the user must also specify domain of either the target or the random start point.
    module function test_SpecMCMC_RandomStartPointRequested_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: RandomStartPointDomainLowerLimitVec(NDIM) = [-1.e0_RK, +1.e1_RK]
        real(RK)    , parameter :: RandomStartPointDomainUpperLimitVec(NDIM) = [+2.e0_RK, +2.e1_RK]
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_RandomStartPointRequested_type_3" &
                            , inputFile = ParaDXXX_NML//" randomStartPointRequested = true /" &
                            , RandomStartPointDomainLowerLimitVec = RandomStartPointDomainLowerLimitVec &
                            , RandomStartPointDomainUpperLimitVec = RandomStartPointDomainUpperLimitVec &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. PD%SpecMCMC%RandomStartPointRequested%val
#endif
    end function test_SpecMCMC_RandomStartPointRequested_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> When `randomStartPointRequested = true`, the user must also specify domain of either the target or the random start point.
    module function test_SpecMCMC_RandomStartPointRequested_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: DomainLowerLimitVec(NDIM) = [-1.e0_RK, +1.e2_RK]
        real(RK)    , parameter :: DomainUpperLimitVec(NDIM) = [+2.e0_RK, +2.e2_RK]
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_RandomStartPointRequested_type_4" &
                            , inputFile = ParaDXXX_NML//" randomStartPointRequested = true /" &
                            , DomainLowerLimitVec = DomainLowerLimitVec &
                            , DomainUpperLimitVec = DomainUpperLimitVec &
                            , chainSize = 100_IK &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. PD%SpecMCMC%RandomStartPointRequested%val
#endif
    end function test_SpecMCMC_RandomStartPointRequested_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler returns an error message when `sampleRefinementCount < 0`.
    module function test_SpecMCMC_SampleRefinementCount_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
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
    !> Test whether the ParaDXXX sampler returns an error message when `sampleRefinementCount < 0`.
    module function test_SpecMCMC_SampleRefinementCount_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementCount_type_2" &
                            , inputFile = ParaDXXX_NML//" sampleRefinementCount = -1 /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementCount_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler succeeds with a valid `sampleRefinementCount >= 0`.
    module function test_SpecMCMC_SampleRefinementCount_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementCount_type_3" &
                            , inputFile = ParaDXXX_NML//" sampleRefinementCount = 0 /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementCount_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler returns an error message with an unrecognized value for `sampleRefinementMethod`.
    module function test_SpecMCMC_SampleRefinementMethod_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementMethod_type_1" &
                            , inputFile = ParaDXXX_NML//" sampleRefinementMethod = 'nonsense' /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementMethod_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler succeeds with a valid value for `sampleRefinementMethod`.
    module function test_SpecMCMC_SampleRefinementMethod_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
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
    !> Test whether the ParaDXXX sampler succeeds with a valid value for `sampleRefinementMethod`.
    module function test_SpecMCMC_SampleRefinementMethod_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementMethod_type_3" &
                            , inputFile = ParaDXXX_NML//" sampleRefinementMethod = 'CutOffAutoCorr' /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementMethod_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler succeeds with a valid value for `sampleRefinementMethod`.
    module function test_SpecMCMC_SampleRefinementMethod_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementMethod_type_4" &
                            , inputFile = ParaDXXX_NML//" sampleRefinementMethod = 'CutOf  fAuto Corr' /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementMethod_type_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler succeeds with a valid value for `sampleRefinementMethod`.
    module function test_SpecMCMC_SampleRefinementMethod_type_5() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementMethod_type_5" &
                            , inputFile = ParaDXXX_NML//" sampleRefinementMethod = 'MaxCum SumAut oCorr' /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementMethod_type_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler succeeds with a valid value for `sampleRefinementMethod`.
    module function test_SpecMCMC_SampleRefinementMethod_type_6() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementMethod_type_6" &
                            , inputFile = ParaDXXX_NML//" sampleRefinementMethod = 'batchmeans-avg' /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementMethod_type_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler succeeds with a valid value for `sampleRefinementMethod`.
    module function test_SpecMCMC_SampleRefinementMethod_type_7() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementMethod_type_7" &
                            , inputFile = ParaDXXX_NML//" sampleRefinementMethod = 'batchmeans-average' /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementMethod_type_7

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler succeeds with a valid value for `sampleRefinementMethod`.
    module function test_SpecMCMC_SampleRefinementMethod_type_8() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementMethod_type_8" &
                            , inputFile = ParaDXXX_NML//" sampleRefinementMethod = 'batchmeans-med' /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementMethod_type_8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler succeeds with a valid value for `sampleRefinementMethod`.
    module function test_SpecMCMC_SampleRefinementMethod_type_9() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementMethod_type_9" &
                            , inputFile = ParaDXXX_NML//" sampleRefinementMethod = 'batchmeans-median' /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementMethod_type_9

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler succeeds with a valid value for `sampleRefinementMethod`.
    module function test_SpecMCMC_SampleRefinementMethod_type_10() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementMethod_type_10" &
                            , inputFile = ParaDXXX_NML//" sampleRefinementMethod = 'batchmeans-min' /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementMethod_type_10

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler succeeds with a valid value for `sampleRefinementMethod`.
    module function test_SpecMCMC_SampleRefinementMethod_type_11() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementMethod_type_11" &
                            , inputFile = ParaDXXX_NML//" sampleRefinementMethod = 'batchmeans-minimum' /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementMethod_type_11

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler succeeds with a valid value for `sampleRefinementMethod`.
    module function test_SpecMCMC_SampleRefinementMethod_type_12() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementMethod_type_12" &
                            , inputFile = ParaDXXX_NML//" sampleRefinementMethod = 'batchmeans-max' /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementMethod_type_12

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler succeeds with a valid value for `sampleRefinementMethod`.
    module function test_SpecMCMC_SampleRefinementMethod_type_13() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_SampleRefinementMethod_type_13" &
                            , inputFile = ParaDXXX_NML//" sampleRefinementMethod = 'batchmeans-maximum' /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_SampleRefinementMethod_type_13

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler returns with an error message for a wrong input value for `scaleFactor`.
    module function test_SpecMCMC_ScaleFactor_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ScaleFactor_type_1" &
                            , inputFile = ParaDXXX_NML//" scaleFactor = 'nonsense' /" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ScaleFactor_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler returns with an error message for an empty input value for `scaleFactor`.
    module function test_SpecMCMC_ScaleFactor_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
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
    !> Test whether the ParaDXXX sampler returns with an error message for an empty input value for `scaleFactor`.
    module function test_SpecMCMC_ScaleFactor_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
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
    !> Test whether the ParaDXXX sampler returns with an error message for a wrong input value for `scaleFactor`.
    module function test_SpecMCMC_ScaleFactor_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
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
    !> Test whether the ParaDXXX sampler returns successfully for a valid input value for `scaleFactor`.
    module function test_SpecMCMC_ScaleFactor_type_5() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ScaleFactor_type_5" &
                            , inputFile = ParaDXXX_NML//" scaleFactor = '2 * Gelman'" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_ScaleFactor_type_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler returns with an error message for a negative input value for `scaleFactor`.
    module function test_SpecMCMC_ScaleFactor_type_6() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ScaleFactor_type_6" &
                            , inputFile = ParaDXXX_NML//" scaleFactor = '-0.5 * Gelman'" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ScaleFactor_type_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler returns with an error message for a wrong input value for `scaleFactor`.
    module function test_SpecMCMC_ScaleFactor_type_7() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ScaleFactor_type_7" &
                            , inputFile = ParaDXXX_NML//" scaleFactor = '*GELMAN'" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_SpecMCMC_ScaleFactor_type_7

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler returns successfully for a valid input value for `scaleFactor`.
    module function test_SpecMCMC_ScaleFactor_type_8() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_ScaleFactor_type_8" &
                            , inputFile = ParaDXXX_NML//" scaleFactor = '2. * 0.5 * Gelman'" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. abs(PD%SpecMCMC%ScaleFactor%val-2.38_RK/sqrt(real(NDIM,kind=RK)))<1.e-12_RK
#endif
    end function test_SpecMCMC_ScaleFactor_type_8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler returns successfully for a valid input value for `StartPointVec`.
    module function test_SpecMCMC_StartPointVec_type_1() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: StartPointVec(NDIM) = [1._RK, -10._RK]
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
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
    !> Test whether the ParaDXXX sampler returns successfully for a valid input value for `StartPointVec`.
    module function test_SpecMCMC_StartPointVec_type_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_StartPointVec_type_2" &
                            , inputFile = ParaDXXX_NML//" StartPointVec = 1., -10., /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_SpecMCMC_StartPointVec_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the ParaDXXX sampler returns successfully and correctly computes `StartPointVec` when not provided.
    module function test_SpecMCMC_StartPointVec_type_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
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
    !> Test whether the ParaDXXX sampler returns successfully and correctly computes `StartPointVec` when not provided.
    module function test_SpecMCMC_StartPointVec_type_4() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: DomainLowerLimitVec(NDIM) = [-1.e0_RK, +1.e1_RK]
        real(RK)    , parameter :: DomainUpperLimitVec(NDIM) = [+2.e0_RK, +2.e1_RK]
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
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

    !> \brief
    !> Test whether the ParaDXXX sampler returns with an error message when `StartPointVec` is out of domain boundary.
    module function test_SpecMCMC_StartPointVec_type_5() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaDXXX_type)     :: PD
        integer(IK) , parameter :: NDIM = 2_IK
        real(RK)    , parameter :: DomainLowerLimitVec(NDIM) = [-1.e0_RK, +1.e1_RK]
        real(RK)    , parameter :: DomainUpperLimitVec(NDIM) = [+2.e0_RK, +2.e1_RK]
        real(RK)    , parameter :: StartPointVec(NDIM) = [-2.e1_RK, +2.e2_RK]
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecMCMC/test_SpecMCMC_StartPointVec_type_5" &
                            , DomainLowerLimitVec = DomainLowerLimitVec &
                            , DomainUpperLimitVec = DomainUpperLimitVec &
                            , StartPointVec = StartPointVec &
                            , chainSize = 100_IK &
                            )
        assertion = assertion .and. PD%Err%occurred .and. all(PD%SpecMCMC%StartPointVec%Val==StartPointVec)
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
    end function test_SpecMCMC_StartPointVec_type_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
