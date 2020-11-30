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

!>  \brief This module contains tests of the module [CrossCorr_mod](@ref crosscorr_mod).
!>  @author Amir Shahmoradi

module Test_CrossCorr_mod

    use Test_mod, only: Test_type
    use Constants_mod, only: IK, RK
    use CrossCorr_mod
    implicit none

    private
    public :: Test_CrossCorr

    type(Test_type) :: Test

    ! input Test data

    type :: WeightedData_type
        integer(IK)                 :: nd = 1_IK
        integer(IK)                 :: np = 9985_IK
        integer(IK) , allocatable   :: Weight(:)
        real(RK)    , allocatable   :: Data(:,:)
        real(RK)    , allocatable   :: NormedData(:,:)
        real(RK)    , allocatable   :: ref_AutoCorr(:,:)
        real(RK)    , allocatable   :: InverseSumNormedDataSq(:)
        real(RK)    , allocatable   :: InverseSumNormedDataSq_ref(:)
    contains
        procedure   , pass          :: read => readData
    end type WeightedData_type
    type(WeightedData_type)         :: WeightedData !< An object of class WeightedData_type containing the input Compact-Weighted Data (WCD).

    ! Computed Autocorrelation

    type :: AutoCorr_type
        real(RK)    , allocatable   :: NormedDataFFT1(:,:), NormedDataFFT2(:,:)
        real(RK)    , allocatable   :: AutoCorrWeightedFFT(:,:)
        real(RK)    , allocatable   :: AutoCorrDirect_ref(:,:)
        real(RK)    , allocatable   :: AutoCorrFFT_ref(:,:)
        real(RK)    , allocatable   :: AutoCorrDirect(:,:)
        real(RK)    , allocatable   :: AutoCorrFFT(:,:)
        integer(IK) , allocatable   :: Lag_ref(:)
        integer(IK) , allocatable   :: Lag(:)
        integer(IK)                 :: nlag
        integer(IK)                 :: paddedLen
    end type AutoCorr_type
    type(AutoCorr_type)             :: AutoCorr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine Test_CrossCorr()

        use Constants_mod, only: IK, RK
        implicit none

        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(Test_padZero_1, "Test_padZero_1")
        call Test%run(Test_padZero_2, "Test_padZero_2")
        call Test%run(Test_getCumSumIAC_1, "Test_getCumSumIAC_1")
        call Test%run(Test_getCumSumIAC_2, "Test_getCumSumIAC_2")
        call Test%run(Test_getCumSumIAC_3, "Test_getCumSumIAC_3")
        call Test%run(Test_getCrossCorrFFT, "Test_getCrossCorrFFT")
        call Test%run(Test_getMaxCumSumIAC_1, "Test_getMaxCumSumIAC_1")
        call Test%run(Test_getMaxCumSumIAC_2, "Test_getMaxCumSumIAC_2")
        call Test%run(Test_getPaddedLen_IK_1, "Test_getPaddedLen_IK_1")
        call Test%run(Test_getPaddedLen_RK_2, "Test_getPaddedLen_RK_2")
        call Test%run(Test_getNextExponent_1, "Test_getNextExponent_1")
        call Test%run(Test_getNextExponent_2, "Test_getNextExponent_2")
        call Test%run(Test_getNextExponent_3, "Test_getNextExponent_3")
        call Test%run(Test_getBatchMeansIAC_1, "Test_getBatchMeansIAC_1")
        call Test%run(Test_getBatchMeansIAC_2, "Test_getBatchMeansIAC_2")
        call Test%run(Test_getBatchMeansIAC_3, "Test_getBatchMeansIAC_3")
        call Test%run(Test_getAutoCorrDirect_1, "Test_getAutoCorrDirect_1")
        call Test%run(Test_getAutoCorrDirect_2, "Test_getAutoCorrDirect_2")
        call Test%run(Test_getPreviousExponent_1, "Test_getPreviousExponent_1")
        call Test%run(Test_getPreviousExponent_2, "Test_getPreviousExponent_2")
        call Test%run(Test_getInverseSumNormedDataSq_1, "Test_getInverseSumNormedDataSq_1")
        call Test%run(Test_getCrossCorrWeightedFFT_1, "Test_getCrossCorrWeightedFFT_1")
        call Test%run(Test_getCrossCorrWeightedFFT_2, "Test_getCrossCorrWeightedFFT_2")
        call Test%run(Test_getCrossCorrWeightedFFT_3, "Test_getCrossCorrWeightedFFT_3")
        call Test%finalize()

    end subroutine Test_CrossCorr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_padZero_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                     :: assertion
        real(RK)    , parameter     :: tolerance = 1.e-14_RK
        real(RK)    , parameter     :: CurrentArray(*) = [1._RK, 1._RK, 1._RK, 1._RK]
        real(RK)    , parameter     :: PaddedArray_ref(*) = [1._RK, 1._RK, 1._RK, 1._RK, 0._RK, 0._RK, 0._RK, 0._RK, 0._RK, 0._RK]
        integer(IK) , parameter     :: lenCurrentArray = size(CurrentArray)
        integer(IK) , parameter     :: paddedLen = size(PaddedArray_ref)
        real(RK)                    :: PaddedArray(size(PaddedArray_ref)) ! Gfortran 7.1 fails to automatically allocate an allocatable version of these arrays
        real(RK)                    :: Difference(size(PaddedArray_ref))  ! Gfortran 7.1 fails to automatically allocate an allocatable version of these arrays
        PaddedArray = padZero(currentLen = lenCurrentArray, Array = CurrentArray, paddedLen = paddedLen)
        Difference = abs(PaddedArray - PaddedArray_ref)
        assertion = all(Difference < tolerance)
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "CurrentArray      =", CurrentArray
            write(Test%outputUnit,"(*(g0.15,:,' '))") "PaddedArray_ref   =", PaddedArray_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "PaddedArray       =", PaddedArray
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Difference        =", Difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function Test_padZero_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_padZero_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                     :: assertion
        real(RK)    , parameter     :: tolerance = 1.e-14_RK
        real(RK)    , parameter     :: CurrentArray(*) = [1._RK, 1._RK, 1._RK, 1._RK]
        real(RK)    , parameter     :: PaddedArray_ref(*) = [1._RK, 1._RK, 1._RK, 1._RK, 0._RK, 0._RK, 0._RK, 0._RK]
        integer(IK) , parameter     :: lenCurrentArray = size(CurrentArray)
        integer(IK) , parameter     :: paddedLen = size(PaddedArray_ref)
        real(RK)                    :: PaddedArray(size(PaddedArray_ref)) ! Gfortran 7.1 fails to automatically allocate an allocatable version of these arrays
        real(RK)                    :: Difference(size(PaddedArray_ref))  ! Gfortran 7.1 fails to automatically allocate an allocatable version of these arrays
        PaddedArray = padZero(currentLen = lenCurrentArray, Array = CurrentArray)
        Difference = abs(PaddedArray - PaddedArray_ref)
        assertion = all(Difference < tolerance)
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "CurrentArray      =", CurrentArray
            write(Test%outputUnit,"(*(g0.15,:,' '))") "PaddedArray_ref   =", PaddedArray_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "PaddedArray       =", PaddedArray
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Difference        =", Difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function Test_padZero_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_getPaddedLen_IK_1() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        integer(IK) , parameter     :: actualLen = 4_IK
        integer(IK) , parameter     :: paddedLen_ref = 27_IK
        integer(IK) , parameter     :: base = 3_IK
        integer(IK)                 :: paddedLen
        integer(IK)                 :: difference
        paddedLen = getPaddedLen(actualLen = actualLen, base = base)
        difference = abs(paddedLen - paddedLen_ref)
        assertion = difference == 0_IK
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "actualLen     =", actualLen
            write(Test%outputUnit,"(*(g0.15,:,' '))") "paddedLen_ref =", paddedLen_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "paddedLen     =", paddedLen
            write(Test%outputUnit,"(*(g0.15,:,' '))") "difference    =", difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function Test_getPaddedLen_IK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_getPaddedLen_IK_2() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        integer(IK) , parameter     :: actualLen = 4_IK
        integer(IK) , parameter     :: paddedLen_ref = 8_IK
        integer(IK)                 :: paddedLen
        integer(IK)                 :: difference
        paddedLen = getPaddedLen(actualLen = actualLen)
        difference = abs(paddedLen - paddedLen_ref)
        assertion = difference == 0_IK
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "actualLen     =", actualLen
            write(Test%outputUnit,"(*(g0.15,:,' '))") "paddedLen_ref =", paddedLen_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "paddedLen     =", paddedLen
            write(Test%outputUnit,"(*(g0.15,:,' '))") "difference    =", difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function Test_getPaddedLen_IK_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_getPaddedLen_RK_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                     :: assertion
        real(RK)    , parameter     :: actualLen = 4._RK
        integer(IK) , parameter     :: paddedLen_ref = 43_IK
        real(RK)    , parameter     :: base = 3.5_RK
        integer(IK)                 :: paddedLen
        integer(IK)                 :: difference
        paddedLen = getPaddedLen(actualLen = actualLen, base = base)
        difference = abs(paddedLen - paddedLen_ref)
        assertion = difference == 0_IK
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "actualLen     =", actualLen
            write(Test%outputUnit,"(*(g0.15,:,' '))") "paddedLen_ref =", paddedLen_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "paddedLen     =", paddedLen
            write(Test%outputUnit,"(*(g0.15,:,' '))") "difference    =", difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function Test_getPaddedLen_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_getPaddedLen_RK_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                     :: assertion
        real(RK)    , parameter     :: actualLen = 4._RK
        integer(IK) , parameter     :: paddedLen_ref = 8_IK
        real(RK)    , parameter     :: base = 2._RK
        integer(IK)                 :: paddedLen
        integer(IK)                 :: difference
        paddedLen = getPaddedLen(actualLen = actualLen, base = base)
        difference = abs(paddedLen - paddedLen_ref)
        assertion = difference == 0_IK
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "actualLen     =", actualLen
            write(Test%outputUnit,"(*(g0.15,:,' '))") "paddedLen_ref =", paddedLen_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "paddedLen     =", paddedLen
            write(Test%outputUnit,"(*(g0.15,:,' '))") "difference    =", difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function Test_getPaddedLen_RK_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_getNextExponent_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK)             :: difference
        integer(IK)             :: nextExponent
        integer(IK) , parameter :: nextExponent_ref = 4_IK
        real(RK)    , parameter :: absoluteValue = 10._RK
        real(RK)    , parameter :: base = 2._RK
        nextExponent = getNextExponent(absoluteValue = absoluteValue, base = base)
        difference = abs(nextExponent - nextExponent_ref)
        assertion = difference == 0_IK
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "nextExponent_ref      =", nextExponent_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Computed nextExponent =", nextExponent
            write(Test%outputUnit,"(*(g0.15,:,' '))") "difference            =", difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function Test_getNextExponent_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_getNextExponent_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK)             :: difference
        integer(IK)             :: nextExponent
        integer(IK) , parameter :: nextExponent_ref = 4_IK
        real(RK)    , parameter :: absoluteValue = 16._RK
        real(RK)    , parameter :: base = 2.5_RK
        nextExponent = getNextExponent(absoluteValue = absoluteValue, base = base)
        difference = abs(nextExponent - nextExponent_ref)
        assertion = difference == 0_IK
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "nextExponent_ref      =", nextExponent_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Computed nextExponent =", nextExponent
            write(Test%outputUnit,"(*(g0.15,:,' '))") "difference            =", difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function Test_getNextExponent_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_getNextExponent_3() result(assertion)
        use Constants_mod, only: RK
        implicit none
        logical                 :: assertion
        integer(IK)             :: difference
        integer(IK)             :: nextExponent
        integer(IK) , parameter :: nextExponent_ref = 4_IK
        real(RK)    , parameter :: absoluteValue = 16._RK
        nextExponent = getNextExponent(absoluteValue = absoluteValue)
        difference = abs(nextExponent - nextExponent_ref)
        assertion = difference == 0_IK
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "nextExponent_ref      =", nextExponent_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Computed nextExponent =", nextExponent
            write(Test%outputUnit,"(*(g0.15,:,' '))") "difference            =", difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function Test_getNextExponent_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_getPreviousExponent_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK)             :: difference
        integer(IK)             :: previousExponent
        integer(IK) , parameter :: previousExponent_ref = 13_IK
        real(RK)    , parameter :: absoluteValue = 9985._RK
        previousExponent = getPreviousExponent(absoluteValue = absoluteValue)
        difference = abs(previousExponent - previousExponent_ref)
        assertion = difference == 0_IK
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "previousExponent_ref    =", previousExponent_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "previousExponent        =", previousExponent
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Difference              =", difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function Test_getPreviousExponent_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_getPreviousExponent_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK)             :: difference
        integer(IK)             :: previousExponent
        integer(IK) , parameter :: previousExponent_ref = 7_IK
        real(RK)    , parameter :: absoluteValue = 9985._RK
        real(RK)    , parameter :: base = 3.5_RK
        previousExponent = getPreviousExponent(absoluteValue = absoluteValue, base = base)
        difference = abs(previousExponent - previousExponent_ref)
        assertion = difference == 0_IK
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "previousExponent_ref    =", previousExponent_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "previousExponent        =", previousExponent
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Difference              =", difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function Test_getPreviousExponent_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine readData(WeightedData)

        use Statistics_mod, only: getNormData
        implicit none
        class(WeightedData_type), intent(inout) :: WeightedData
        integer(IK) :: ip, fileUnit

        ! read the input data required for other Tests

        WeightedData%nd = 1_IK
        WeightedData%np = 9985_IK
        if (allocated(WeightedData%Weight)) deallocate(WeightedData%Weight); allocate(WeightedData%Weight(WeightedData%np))
        if (allocated(WeightedData%Data)) deallocate(WeightedData%Data); allocate(WeightedData%Data(WeightedData%nd,WeightedData%np))

        open(newunit=fileUnit,file=Test%inDir//"Test_CrossCorr_mod@WeightedData.txt",status="old")
        do ip = 1, WeightedData%np
            read(fileUnit,*) WeightedData%Weight(ip), WeightedData%Data(1:WeightedData%nd,ip)
        end do
        close(fileUnit)

        ! normalize data

        WeightedData%NormedData = getNormData(WeightedData%nd, WeightedData%np, WeightedData%Data)
        WeightedData%InverseSumNormedDataSq_ref = [ 5.935290321338481E-004_RK ]

    end subroutine readData

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_getInverseSumNormedDataSq_1() result(assertion)

        use Constants_mod, only: RK
        implicit none
        logical                 :: assertion
        real(RK), allocatable   :: Difference(:)
        real(RK), parameter     :: tolerance = 1.e-12_RK
        call WeightedData%read()

        WeightedData%InverseSumNormedDataSq = getInverseSumNormedDataSq(1_IK, WeightedData%np, WeightedData%NormedData)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(Difference)) deallocate(Difference); allocate(Difference, mold = WeightedData%InverseSumNormedDataSq)

        Difference = abs( (WeightedData%InverseSumNormedDataSq - WeightedData%InverseSumNormedDataSq_ref) / WeightedData%InverseSumNormedDataSq_ref)
        assertion = all( Difference < tolerance )
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "InverseSumNormedDataSq_ref  =", WeightedData%InverseSumNormedDataSq_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "InverseSumNormedDataSq      =", WeightedData%InverseSumNormedDataSq
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Difference                  =", Difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function Test_getInverseSumNormedDataSq_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! compute AutoCorr for the given lags using the classical definition.
    function Test_getAutoCorrDirect_1() result(assertion)

        use Constants_mod, only: RK
        use String_mod, only: num2str
        implicit none
        logical                     :: assertion
        integer(IK) , allocatable   :: DifferenceLag(:)
        real(RK)    , allocatable   :: DifferenceAutoCorrDirect(:,:)
        real(RK)    , parameter     :: tolerance = 1.e-12_RK
        integer(IK)                 :: ilag, fileUnit

        call WeightedData%read()

        ! Generate and verify Lags

        AutoCorr%Lag_ref = [ 0_IK, 1_IK, 2_IK, 4_IK, 8_IK, 16_IK, 32_IK, 64_IK, 128_IK, 256_IK, 512_IK, 1024_IK, 2048_IK, 4096_IK, 8192_IK ]
        AutoCorr%AutoCorrDirect_ref = reshape(  [ 1._RK &
                                                , .8985120580850618_RK &
                                                , .8142306385443645_RK &
                                                , .6793649391868671_RK &
                                                , .5018573940223264_RK &
                                                , .2963258939172878_RK &
                                                , .2018267302958381_RK &
                                                , .1415530903106628_RK &
                                                , .5232924941049354E-01_RK &
                                                , .4982051840374713E-01_RK &
                                                , .5464128642962141E-01_RK &
                                                , -.2355552565791153E-01_RK &
                                                , .2584696857029676E-01_RK &
                                                , -.4164933936801392E-01_RK &
                                                , .4757873176344895E-02_RK &
                                                ], shape = [ 1, size(AutoCorr%Lag_ref) ] )
        AutoCorr%nlag = getPreviousExponent( real(WeightedData%np, kind=RK) ) + 1
        AutoCorr%Lag = [ 0, ( 2_IK**(ilag-1), ilag = 1, AutoCorr%nlag ) ]

        ! Gfortran 7.1 fails to automatically allocate an allocatable version of these arrays
        if (allocated(DifferenceLag)) deallocate(DifferenceLag); allocate(DifferenceLag, mold = AutoCorr%Lag)

        DifferenceLag = abs( AutoCorr%Lag - AutoCorr%Lag_ref )
        assertion = all( DifferenceLag == 0_IK )

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "AutoCorr%Lag_ref    =", AutoCorr%Lag_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "AutoCorr%Lag        =", AutoCorr%Lag
            write(Test%outputUnit,"(*(g0.15,:,' '))") "nlag                =", AutoCorr%nlag
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

        ! Generate and verify AutoCorrs

        if (allocated(AutoCorr%AutoCorrDirect)) deallocate(AutoCorr%AutoCorrDirect); allocate( AutoCorr%AutoCorrDirect(WeightedData%nd,AutoCorr%nlag + 1_IK) )
        call getAutoCorrDirect  ( nd = WeightedData%nd &
                                , np = WeightedData%np &
                                , NormedData = WeightedData%NormedData(1:WeightedData%nd,1:WeightedData%np) &
                                , nlag = AutoCorr%nlag + 1_IK &
                                , Lag = AutoCorr%Lag &
                                , AutoCorr = AutoCorr%AutoCorrDirect &
                                , InverseSumNormedDataSq = WeightedData%InverseSumNormedDataSq &
                                )

        ! Gfortran 7.1 fails to automatically allocate an reallocatable version of these arrays
        if (allocated(DifferenceAutoCorrDirect)) deallocate(DifferenceAutoCorrDirect); allocate(DifferenceAutoCorrDirect, mold = AutoCorr%AutoCorrDirect)

        DifferenceAutoCorrDirect = abs( AutoCorr%AutoCorrDirect - AutoCorr%AutoCorrDirect_ref )

        assertion = assertion .and. all( DifferenceAutoCorrDirect < tolerance )

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then

            ! write data to output for further investigation

            open(newunit=fileUnit,file=Test%outDir//"Test_CrossCorr_mod@WeightedData@AutoCorr@getAutoCorrDirect."//num2str(Test%Image%id)//".txt",status="replace")
            do ilag = 1, AutoCorr%nlag + 1_IK
                write(fileUnit,"(*(g0.15,:,' '))") AutoCorr%Lag(ilag), AutoCorr%AutoCorrDirect(1:WeightedData%nd,ilag)
            end do
            close(fileUnit)

            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "AutoCorr%AutoCorrDirect_ref =", AutoCorr%AutoCorrDirect_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "AutoCorr%AutoCorrDirect     =", AutoCorr%AutoCorrDirect
            write(Test%outputUnit,"(*(g0.15,:,' '))") "DifferenceAutoCorrDirect    =", DifferenceAutoCorrDirect
            write(Test%outputUnit,"(*(g0.15,:,' '))")

        end if
        ! LCOV_EXCL_STOP

    end function Test_getAutoCorrDirect_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! compute AutoCorr for the given lags using the classical definition, but without passing `InverseSumNormedDataSq`.
    function Test_getAutoCorrDirect_2() result(assertion)

        use Constants_mod, only: RK
        use String_mod, only: num2str
        implicit none
        logical                     :: assertion
        integer(IK) , allocatable   :: DifferenceLag(:)
        real(RK)    , allocatable   :: DifferenceAutoCorrDirect(:,:)
        real(RK)    , parameter     :: tolerance = 1.e-12_RK
        integer(IK)                 :: ilag, fileUnit

        call WeightedData%read()

        ! Generate and verify Lags

        AutoCorr%Lag_ref = [ 0_IK, 1_IK, 2_IK, 4_IK, 8_IK, 16_IK, 32_IK, 64_IK, 128_IK, 256_IK, 512_IK, 1024_IK, 2048_IK, 4096_IK, 8192_IK ]
        AutoCorr%AutoCorrDirect_ref = reshape(  [ 1._RK &
                                                , .8985120580850618_RK &
                                                , .8142306385443645_RK &
                                                , .6793649391868671_RK &
                                                , .5018573940223264_RK &
                                                , .2963258939172878_RK &
                                                , .2018267302958381_RK &
                                                , .1415530903106628_RK &
                                                , .5232924941049354E-01_RK &
                                                , .4982051840374713E-01_RK &
                                                , .5464128642962141E-01_RK &
                                                , -.2355552565791153E-01_RK &
                                                , .2584696857029676E-01_RK &
                                                , -.4164933936801392E-01_RK &
                                                , .4757873176344895E-02_RK &
                                                ], shape = [ 1, size(AutoCorr%Lag_ref) ] )
        AutoCorr%nlag = getPreviousExponent( real(WeightedData%np, kind=RK) ) + 1
        AutoCorr%Lag = [ 0, ( 2_IK**(ilag-1), ilag = 1, AutoCorr%nlag ) ]

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(DifferenceLag)) deallocate(DifferenceLag); allocate(DifferenceLag, mold = AutoCorr%Lag)

        DifferenceLag = abs( AutoCorr%Lag - AutoCorr%Lag_ref )
        assertion = all( DifferenceLag == 0_IK )

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "AutoCorr%Lag_ref    =", AutoCorr%Lag_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "AutoCorr%Lag        =", AutoCorr%Lag
            write(Test%outputUnit,"(*(g0.15,:,' '))") "nlag                =", AutoCorr%nlag
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

        ! Generate and verify AutoCorrs

        if (allocated(AutoCorr%AutoCorrDirect)) deallocate(AutoCorr%AutoCorrDirect); allocate( AutoCorr%AutoCorrDirect(WeightedData%nd,AutoCorr%nlag+1_IK) )
        call getAutoCorrDirect  ( nd = WeightedData%nd &
                                , np = WeightedData%np &
                                , NormedData = WeightedData%NormedData(1:WeightedData%nd,1:WeightedData%np) &
                                , nlag = AutoCorr%nlag + 1_IK &
                                , Lag = AutoCorr%Lag &
                                , AutoCorr = AutoCorr%AutoCorrDirect &
                                )

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(DifferenceAutoCorrDirect)) deallocate(DifferenceAutoCorrDirect); allocate(DifferenceAutoCorrDirect, mold = AutoCorr%AutoCorrDirect)

        DifferenceAutoCorrDirect = abs( AutoCorr%AutoCorrDirect - AutoCorr%AutoCorrDirect_ref )

        assertion = assertion .and. all( DifferenceAutoCorrDirect <= tolerance )

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then

            ! write data to output for further investigation

            open(newunit=fileUnit,file=Test%outDir//"WeightedDataAutoCorrDirect.Without.InverseSumNormedDataSq."//num2str(Test%Image%id)//".txt",status="replace")
            do ilag = 1, AutoCorr%nlag + 1_IK
                write(fileUnit,"(*(g0.15,:,' '))") AutoCorr%Lag(ilag), AutoCorr%AutoCorrDirect(1:WeightedData%nd,ilag)
            end do
            close(fileUnit)

            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "AutoCorr%AutoCorrDirect_ref =", AutoCorr%AutoCorrDirect_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "AutoCorr%AutoCorrDirect     =", AutoCorr%AutoCorrDirect
            write(Test%outputUnit,"(*(g0.15,:,' '))") "DifferenceAutoCorrDirect    =", DifferenceAutoCorrDirect
            write(Test%outputUnit,"(*(g0.15,:,' '))")

        end if
        ! LCOV_EXCL_STOP

    end function Test_getAutoCorrDirect_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Test AutoCorr computation for Compact Data using FFT.
    function Test_getCrossCorrFFT() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                     :: assertion, assertionCurrent
        real(RK)    , allocatable   :: Difference(:,:)
        real(RK)    , allocatable   :: DummyArray(:)
        real(RK)    , parameter     :: tolerance = 1.e-12_RK
        integer(IK)                 :: id, ip, fileUnit, ilag

        assertion = .true.

        call WeightedData%read()

        AutoCorr%paddedLen = getPaddedLen(WeightedData%np)
        allocate( AutoCorr%AutoCorrFFT_ref(AutoCorr%paddedLen,WeightedData%nd) &
                , AutoCorr%NormedDataFFT1(AutoCorr%paddedLen,WeightedData%nd) &
                , AutoCorr%NormedDataFFT2(AutoCorr%paddedLen,WeightedData%nd) &
                , AutoCorr%AutoCorrFFT(AutoCorr%paddedLen,WeightedData%nd) &
                )
        do id = 1, WeightedData%nd
            DummyArray = WeightedData%NormedData(id,1:WeightedData%np)
            AutoCorr%NormedDataFFT1(1:AutoCorr%paddedLen,id) = padZero(WeightedData%np, DummyArray, AutoCorr%paddedLen)
        end do
        AutoCorr%NormedDataFFT2 = AutoCorr%NormedDataFFT1

        do id = 1, WeightedData%nd

            AutoCorr%AutoCorrFFT(1:AutoCorr%paddedLen,id) = getCrossCorrFFT ( AutoCorr%paddedLen &
                                                                            , AutoCorr%NormedDataFFT1(1:AutoCorr%paddedLen,id) &
                                                                            , AutoCorr%NormedDataFFT1(1:AutoCorr%paddedLen,id) )
            AutoCorr%AutoCorrFFT(1:AutoCorr%paddedLen,id) = AutoCorr%AutoCorrFFT(1:AutoCorr%paddedLen,id) / AutoCorr%AutoCorrFFT(1,id)

            ! ensure NormedDataFFT1 does not change upon entering and exiting getCrossCorrFFT()

            do ip = 1, AutoCorr%paddedLen
                assertionCurrent = AutoCorr%NormedDataFFT2(ip,id) == AutoCorr%NormedDataFFT1(ip,id)
                assertion = assertion .and. assertionCurrent
                if (Test%isDebugMode .and. .not. assertionCurrent) then
                    write(Test%outputUnit,"(*(g0.15,:,' '))") "Nonsense. These two must be equal: ", AutoCorr%NormedDataFFT1(ip,id), AutoCorr%NormedDataFFT2(ip,id)
                end if
            end do

        end do

        ! read reference AutoCorrFFT_ref from the input data

        open(newunit=fileUnit,file=Test%inDir//"Test_CrossCorr_mod@WeightedData@AutoCorrFFT@Compact.txt",status="old")
        do ip = 1, AutoCorr%paddedLen
            read(fileUnit,*) ilag, AutoCorr%AutoCorrFFT_ref(ip,1)
        end do
        close(fileUnit)

        ! Gfortran 7.1 fails to automatically allocate an allocatable version of these arrays
        if (allocated(Difference)) deallocate(Difference); allocate(Difference, mold = AutoCorr%AutoCorrFFT)

        Difference = abs( AutoCorr%AutoCorrFFT - AutoCorr%AutoCorrFFT_ref )
        assertion = assertion .and. all( Difference <= tolerance )

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then

            ! write data to output for further investigation

            open(newunit=fileUnit,file=Test%outDir//"Test_CrossCorr_mod@WeightedData@AutoCorrFFT@Compact@getCrossCorrFFT."//num2str(Test%Image%id)//".out",status="replace")
            do ip = 1, AutoCorr%paddedLen
                write(fileUnit,"(*(g0.15,:,' '))") ip-1, AutoCorr%AutoCorrFFT(ip,1)
            end do
            close(fileUnit)

        end if
        ! LCOV_EXCL_STOP

    end function Test_getCrossCorrFFT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Test the AutoCorr computation using FFT with missing Weighted Data, that is, for the compact data.
    function Test_getCrossCorrWeightedFFT_1() result(assertion)

        use Statistics_mod, only: getNormData
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                     :: assertion
        real(RK)    , allocatable   :: Difference(:,:)
        real(RK)    , allocatable   :: DummyArray(:)
        real(RK)    , parameter     :: tolerance = 1.e-12_RK
        integer(IK)                 :: id, ip, fileUnit
        integer(IK) , allocatable   :: DiffMaxLoc(:)

        assertion = .true.

        call WeightedData%read()

        ! open(newunit=fileUnit,file=Test%outDir//"NormedData.compact."//num2str(Test%Image%id)//".txt",status="replace")
        ! write(fileUnit,"(*(g0.15,:,/))") WeightedData%NormedData
        ! close(fileUnit)

        if (allocated(AutoCorr%AutoCorrWeightedFFT)) deallocate(AutoCorr%AutoCorrWeightedFFT)
        allocate( AutoCorr%AutoCorrWeightedFFT(AutoCorr%paddedLen, WeightedData%nd) )

        do id = 1, WeightedData%nd
            DummyArray = WeightedData%NormedData(id,1:WeightedData%np)
            AutoCorr%AutoCorrWeightedFFT(1:AutoCorr%paddedLen,id) = getCrossCorrWeightedFFT ( lenCompactData1 = WeightedData%np &
                                                                                            , lenCompactData2 = WeightedData%np &
                                                                                            , paddedLen = AutoCorr%paddedLen &
                                                                                            , CompactData1 = DummyArray &
                                                                                            , CompactData2 = DummyArray &
                                                                                            )
            AutoCorr%AutoCorrWeightedFFT(1:AutoCorr%paddedLen,id) = AutoCorr%AutoCorrWeightedFFT(1:AutoCorr%paddedLen,id) / AutoCorr%AutoCorrWeightedFFT(1,id)
        end do

        ! Gfortran 7.1 fails to automatically allocate an allocatable version of these arrays
        if (allocated(Difference)) deallocate(Difference); allocate(Difference, mold = AutoCorr%AutoCorrWeightedFFT)

        Difference = abs( AutoCorr%AutoCorrWeightedFFT - AutoCorr%AutoCorrFFT_ref )
        assertion = assertion .and. all( Difference <= tolerance )

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then

            ! write data to output for further investigation

            open(newunit=fileUnit,file=Test%outDir//"Test_CrossCorr_mod@WeightedData@AutoCorrFFT@Compact@getCrossCorrWeightedFFT."//num2str(Test%Image%id)//".out",status="replace")
            do ip = 1, AutoCorr%paddedLen
                write(fileUnit,"(*(g0.15,:,' '))") ip-1, AutoCorr%AutoCorrWeightedFFT(ip,1)
            end do
            close(fileUnit)

            !allocate(DiffMaxLoc(WeightedData%nd))
            DiffMaxLoc = maxloc(Difference(1:AutoCorr%paddedLen,1:WeightedData%nd))

            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Stats at max difference:"
            write(Test%outputUnit,"(*(g0.15,:,' '))") "AutoCorr%AutoCorrFFT_ref        =", AutoCorr%AutoCorrFFT_ref(DiffMaxLoc(1),1)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "AutoCorr%AutoCorrWeightedFFT    =", AutoCorr%AutoCorrWeightedFFT(DiffMaxLoc(1),1)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "maxDifference                   =", Difference(DiffMaxLoc(1),1)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "location                        =", DiffMaxLoc(1)
            write(Test%outputUnit,"(*(g0.15,:,' '))")

        end if
        ! LCOV_EXCL_STOP

    end function Test_getCrossCorrWeightedFFT_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Test AutoCorr computation for the verbose data, stored Compact Weighted format, using FFT.
    function Test_getCrossCorrWeightedFFT_2() result(assertion)

        use Statistics_mod, only: getNormData
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                     :: assertion, assertionCurrent
        real(RK)    , allocatable   :: Difference(:,:)
        real(RK)    , parameter     :: tolerance = 1.e-12_RK
        integer(IK)                 :: id, ip, fileUnit, ilag

        assertion = .true.

        call WeightedData%read()

        ! normalize data according to weights

        WeightedData%NormedData = getNormData(WeightedData%nd, WeightedData%np, WeightedData%Data, WeightedData%Weight, tenabled = .true.)
        AutoCorr%NormedDataFFT1 = WeightedData%NormedData

        AutoCorr%paddedLen = getPaddedLen( sum(WeightedData%Weight) )
        if (allocated(AutoCorr%AutoCorrWeightedFFT)) deallocate(AutoCorr%AutoCorrWeightedFFT)
        allocate( AutoCorr%AutoCorrWeightedFFT(AutoCorr%paddedLen, WeightedData%nd) )
        if (allocated(AutoCorr%AutoCorrFFT_ref)) deallocate(AutoCorr%AutoCorrFFT_ref)
        allocate( AutoCorr%AutoCorrFFT_ref(AutoCorr%paddedLen, WeightedData%nd) )
        do id = 1, WeightedData%nd

            AutoCorr%AutoCorrWeightedFFT(1:AutoCorr%paddedLen,id) = getCrossCorrWeightedFFT ( lenCompactData1 = WeightedData%np &
                                                                                            , lenCompactData2 = WeightedData%np &
                                                                                            , paddedLen = AutoCorr%paddedLen &
                                                                                            , CompactData1 = WeightedData%NormedData(1:WeightedData%np,id) &
                                                                                            , CompactData2 = WeightedData%NormedData(1:WeightedData%np,id) &
                                                                                            , Weight1 = WeightedData%Weight(1:WeightedData%np) &
                                                                                            , Weight2 = WeightedData%Weight(1:WeightedData%np) &
                                                                                            )
            AutoCorr%AutoCorrWeightedFFT(1:AutoCorr%paddedLen,id) = AutoCorr%AutoCorrWeightedFFT(1:AutoCorr%paddedLen,id) / AutoCorr%AutoCorrWeightedFFT(1,id)

            ! ensure NormedData does not change upon entering and exiting getCrossCorrFFT()

            do ip = 1, WeightedData%np
                assertionCurrent = AutoCorr%NormedDataFFT1(ip,id) == WeightedData%NormedData(ip,id)
                assertion = assertion .and. assertionCurrent
                ! LCOV_EXCL_START
                if (Test%isDebugMode .and. .not. assertionCurrent) then
                    write(Test%outputUnit,"(*(g0.15,:,' '))") "Nonsense. These two must be equal ", ip, id, WeightedData%NormedData(ip,id), AutoCorr%NormedDataFFT1(ip,id)
                end if
                ! LCOV_EXCL_STOP
            end do

        end do

        ! read reference AutoCorrFFT_ref from the input data

        open(newunit=fileUnit,file=Test%inDir//"Test_CrossCorr_mod@WeightedData@AutoCorrFFT@Verbose.txt",status="old")
        do ip = 1, AutoCorr%paddedLen
            read(fileUnit,*) ilag, AutoCorr%AutoCorrFFT_ref(ip,1)
        end do
        close(fileUnit)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(Difference)) deallocate(Difference); allocate(Difference, mold = AutoCorr%AutoCorrWeightedFFT)

        Difference = abs( AutoCorr%AutoCorrWeightedFFT - AutoCorr%AutoCorrFFT_ref )
        assertion = assertion .and. all( Difference <= tolerance )

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then

            ! write data to output for further investigation

            open(newunit=fileUnit,file=Test%outDir//"Test_CrossCorr_mod@WeightedData@AutoCorrFFT@Verbose@getCrossCorrWeightedFFT."//num2str(Test%Image%id)//".out",status="replace")
            do ip = 1, AutoCorr%paddedLen
                write(fileUnit,"(*(g0.15,:,' '))") ip-1, AutoCorr%AutoCorrWeightedFFT(ip,1)
            end do
            close(fileUnit)

        end if
        ! LCOV_EXCL_STOP

    end function Test_getCrossCorrWeightedFFT_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Test the AutoCorr computation using FFT for an input identity data vector with a weight vector of ones.
    function Test_getCrossCorrWeightedFFT_3() result(assertion)

        use Statistics_mod, only: getNormData
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                     :: assertion, assertionCurrent
        real(RK)    , allocatable   :: Difference(:,:)
        real(RK)    , parameter     :: tolerance = 1.e-12_RK
        integer(IK)                 :: id, ip, fileUnit, ilag

        assertion = .true.

        call WeightedData%read()

        ! normalize data

        WeightedData%nd = 1
        WeightedData%np = 100
        if (allocated(WeightedData%Data)) deallocate(WeightedData%Data)
        if (allocated(WeightedData%Weight)) deallocate(WeightedData%Weight)
        if (allocated(WeightedData%NormedData)) deallocate(WeightedData%NormedData)
        if (allocated(AutoCorr%NormedDataFFT1)) deallocate(AutoCorr%NormedDataFFT1)
        allocate( WeightedData%Data(WeightedData%nd,WeightedData%np), WeightedData%Weight(WeightedData%np) )
        allocate( WeightedData%NormedData(WeightedData%nd,WeightedData%np), AutoCorr%NormedDataFFT1(WeightedData%np,WeightedData%nd) )
        WeightedData%Data = 1._RK
        WeightedData%Data(1,WeightedData%np) = .999999_RK
        WeightedData%Weight = 1_IK
        WeightedData%NormedData = getNormData(WeightedData%nd, WeightedData%np, WeightedData%Data, WeightedData%Weight, tenabled = .true.)
        AutoCorr%NormedDataFFT1 = WeightedData%NormedData

        AutoCorr%paddedLen = getPaddedLen(sum(WeightedData%Weight))
        if (allocated(AutoCorr%AutoCorrWeightedFFT)) deallocate(AutoCorr%AutoCorrWeightedFFT)
        allocate(AutoCorr%AutoCorrWeightedFFT(AutoCorr%paddedLen,WeightedData%nd))
        do id = 1, WeightedData%nd

            AutoCorr%AutoCorrWeightedFFT(1:AutoCorr%paddedLen,id) = getCrossCorrWeightedFFT ( lenCompactData1 = WeightedData%np                             &
                                                                                            , lenCompactData2 = WeightedData%np                             &
                                                                                            , paddedLen = AutoCorr%paddedLen                                &
                                                                                            , CompactData1 = WeightedData%NormedData(1:WeightedData%np,id)  &
                                                                                            , CompactData2 = WeightedData%NormedData(1:WeightedData%np,id)  &
                                                                                            , Weight1 = WeightedData%Weight(1:WeightedData%np)              &
                                                                                            , Weight2 = WeightedData%Weight(1:WeightedData%np)              &
                                                                                            )
            AutoCorr%AutoCorrWeightedFFT(1:AutoCorr%paddedLen,id) = AutoCorr%AutoCorrWeightedFFT(1:AutoCorr%paddedLen,id) / AutoCorr%AutoCorrWeightedFFT(1,id)

            ! ensure NormedData does not change upon entering and exiting getCrossCorrFFT()

            do ip = 1, WeightedData%np
                assertionCurrent = AutoCorr%NormedDataFFT1(ip,id) == WeightedData%NormedData(ip,id)
                assertion = assertion .and. assertionCurrent
                ! LCOV_EXCL_START
                if (Test%isDebugMode .and. .not. assertionCurrent) then
                    write(Test%outputUnit,"(*(g0.15,:,' '))") "Nonsense. These two must be equal ", ip, id, WeightedData%NormedData(ip,id), AutoCorr%NormedDataFFT1(ip,id)
                end if
                ! LCOV_EXCL_STOP
            end do

        end do

        ! read reference AutoCorrFFT_ref from the input data

        if (allocated(AutoCorr%AutoCorrFFT_ref)) deallocate(AutoCorr%AutoCorrFFT_ref)
        allocate( AutoCorr%AutoCorrFFT_ref(AutoCorr%paddedLen, WeightedData%nd) )
        open(newunit=fileUnit,file=Test%inDir//"Test_CrossCorr_mod@AutoCorrFFT@Identity.txt",status="old")
        do ip = 1, AutoCorr%paddedLen
            read(fileUnit,*) ilag, AutoCorr%AutoCorrFFT_ref(ip,1)
        end do
        close(fileUnit)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(Difference)) deallocate(Difference); allocate(Difference, mold = AutoCorr%AutoCorrWeightedFFT)

        Difference = abs( AutoCorr%AutoCorrWeightedFFT - AutoCorr%AutoCorrFFT_ref )
        assertion = assertion .and. all( Difference <= tolerance )

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then

            ! write data to output for further investigation

            open(newunit=fileUnit,file=Test%outDir//"Test_CrossCorr_mod@AutoCorrFFT@Identity@getCrossCorrWeightedFFT."//num2str(Test%Image%id)//".out",status="replace")
            do ip = 1, AutoCorr%paddedLen
                write(fileUnit,"(*(g0.15,:,' '))") ip-1, AutoCorr%AutoCorrWeightedFFT(ip,1)
            end do
            close(fileUnit)

        end if
        ! LCOV_EXCL_STOP

    end function Test_getCrossCorrWeightedFFT_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_getCumSumIAC_1() result(assertion)
        use Constants_mod, only: RK
        implicit none
        logical                 :: assertion
        real(RK)                :: difference
        real(RK)                :: cumSumIAC
        real(RK), parameter     :: cumSumIAC_ref = 45.0461450858553_RK
        real(RK), parameter     :: tolerance = 1.e-12_RK
        call WeightedData%read()
        cumSumIAC = getCumSumIAC(np = WeightedData%np, Point = WeightedData%Data)
        difference = abs( (cumSumIAC - cumSumIAC_ref) / cumSumIAC_ref)
        assertion = difference < tolerance
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "cumSumIAC_ref    =", cumSumIAC_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "cumSumIAC        =", cumSumIAC
            write(Test%outputUnit,"(*(g0.15,:,' '))") "difference       =", difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function Test_getCumSumIAC_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_getCumSumIAC_2() result(assertion)
        use Constants_mod, only: RK
        implicit none
        logical                 :: assertion
        real(RK)                :: difference
        real(RK)                :: cumSumIAC
        real(RK), parameter     :: cumSumIAC_ref = 240.805734399049_RK
        real(RK), parameter     :: tolerance = 1.e-12_RK
        call WeightedData%read()
        cumSumIAC = getCumSumIAC(np = WeightedData%np, Point = WeightedData%Data, Weight = WeightedData%Weight)
        difference = abs( (cumSumIAC - cumSumIAC_ref) / cumSumIAC_ref)
        assertion = difference < tolerance
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "cumSumIAC_ref    =", cumSumIAC_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "cumSumIAC        =", cumSumIAC
            write(Test%outputUnit,"(*(g0.15,:,' '))") "difference       =", difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function Test_getCumSumIAC_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_getCumSumIAC_3() result(assertion)
        use Constants_mod, only: RK
        implicit none
        logical                 :: assertion
        real(RK)                :: difference
        real(RK)                :: cumSumIAC
        real(RK), parameter     :: significance = 1._RK
        real(RK), parameter     :: cumSumIAC_ref = 240.864727585537_RK
        real(RK), parameter     :: tolerance = 1.e-12_RK
        call WeightedData%read()
        cumSumIAC = getCumSumIAC(np = WeightedData%np, Point = WeightedData%Data, Weight = WeightedData%Weight, significance = significance)
        difference = abs( (cumSumIAC - cumSumIAC_ref) / cumSumIAC_ref)
        assertion = difference < tolerance
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "cumSumIAC_ref    =", cumSumIAC_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "cumSumIAC        =", cumSumIAC
            write(Test%outputUnit,"(*(g0.15,:,' '))") "difference       =", difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function Test_getCumSumIAC_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_getMaxCumSumIAC_1() result(assertion)
        use Constants_mod, only: RK
        implicit none
        logical                 :: assertion
        real(RK)                :: difference
        real(RK)                :: maxCumSumIAC
        real(RK), parameter     :: maxCumSumIAC_ref = 276.9829480126210_RK
        real(RK), parameter     :: tolerance = 1.e-12_RK
        call WeightedData%read()
        maxCumSumIAC = getMaxCumSumIAC(np = WeightedData%np, Point = WeightedData%Data)
        difference = abs( (maxCumSumIAC - maxCumSumIAC_ref) / maxCumSumIAC_ref)
        assertion = difference < tolerance
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "maxCumSumIAC_ref    =", maxCumSumIAC_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "maxCumSumIAC        =", maxCumSumIAC
            write(Test%outputUnit,"(*(g0.15,:,' '))") "difference          =", difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function Test_getMaxCumSumIAC_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_getMaxCumSumIAC_2() result(assertion)
        use Constants_mod, only: RK
        implicit none
        logical                 :: assertion
        real(RK)                :: difference
        real(RK)                :: maxCumSumIAC
        real(RK), parameter     :: maxCumSumIAC_ref = 1111.683837194696_RK
        real(RK), parameter     :: tolerance = 1.e-12_RK
        call WeightedData%read()
        maxCumSumIAC = getMaxCumSumIAC(np = WeightedData%np, Point = WeightedData%Data, Weight = WeightedData%Weight)
        difference = abs( (maxCumSumIAC - maxCumSumIAC_ref) / maxCumSumIAC_ref)
        assertion = difference < tolerance
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "maxCumSumIAC_ref    =", maxCumSumIAC_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "maxCumSumIAC        =", maxCumSumIAC
            write(Test%outputUnit,"(*(g0.15,:,' '))") "difference          =", difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function Test_getMaxCumSumIAC_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_getBatchMeansIAC_1() result(assertion)
        use Constants_mod, only: RK
        implicit none
        logical                 :: assertion
        real(RK)                :: difference
        real(RK)                :: batchMeansIAC
        real(RK), parameter     :: batchMeansIAC_ref = 52.5214820575914_RK
        real(RK), parameter     :: tolerance = 1.e-12_RK
        call WeightedData%read()
        batchMeansIAC = getBatchMeansIAC(np = WeightedData%np, Point = WeightedData%Data)
        difference = abs( (batchMeansIAC - batchMeansIAC_ref) / batchMeansIAC_ref)
        assertion = difference < tolerance
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "batchMeansIAC_ref    =", batchMeansIAC_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "batchMeansIAC        =", batchMeansIAC
            write(Test%outputUnit,"(*(g0.15,:,' '))") "difference           =", difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function Test_getBatchMeansIAC_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_getBatchMeansIAC_2() result(assertion)
        use Constants_mod, only: RK
        implicit none
        logical                 :: assertion
        real(RK)                :: difference
        real(RK)                :: batchMeansIAC
        real(RK), parameter     :: batchMeansIAC_ref = 190.438982011941_RK
        real(RK), parameter     :: tolerance = 1.e-12_RK
        call WeightedData%read()
        batchMeansIAC = getBatchMeansIAC(np = WeightedData%np, Point = WeightedData%Data, Weight = WeightedData%Weight)
        difference = abs( (batchMeansIAC - batchMeansIAC_ref) / batchMeansIAC_ref)
        assertion = difference < tolerance
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "batchMeansIAC_ref    =", batchMeansIAC_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "batchMeansIAC        =", batchMeansIAC
            write(Test%outputUnit,"(*(g0.15,:,' '))") "difference           =", difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function Test_getBatchMeansIAC_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test_getBatchMeansIAC_3() result(assertion)
        use Constants_mod, only: RK
        implicit none
        logical                 :: assertion
        real(RK)                :: difference
        real(RK)                :: batchMeansIAC
        real(RK)    , parameter :: batchMeansIAC_ref = 66.8090158674070_RK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        integer(IK) , parameter :: batchSize = 100_IK
        call WeightedData%read()
        batchMeansIAC = getBatchMeansIAC(np = WeightedData%np, Point = WeightedData%Data, Weight = WeightedData%Weight, batchSize = batchSize)
        difference = abs( (batchMeansIAC - batchMeansIAC_ref) / batchMeansIAC_ref)
        assertion = difference < tolerance
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "batchMeansIAC_ref    =", batchMeansIAC_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "batchMeansIAC        =", batchMeansIAC
            write(Test%outputUnit,"(*(g0.15,:,' '))") "difference           =", difference
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function Test_getBatchMeansIAC_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_CrossCorr_mod ! LCOV_EXCL_LINE