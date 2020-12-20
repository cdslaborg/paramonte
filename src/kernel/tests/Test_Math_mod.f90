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

!>  \brief This module contains tests of the module [Math_mod](@ref math_mod).
!>  \author Amir Shahmoradi

module Test_Math_mod

    use Math_mod
    use Err_mod, only: Err_type
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_Math

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_Math()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_getCumSum_IK_1, "test_getCumSum_IK_1")
        call Test%run(test_getCumSum_RK_1, "test_getCumSum_RK_1")
        call Test%run(test_getFactorial_1, "test_getFactorial_1")
        call Test%run(test_getDistanceSq_1, "test_getDistanceSq_1")
        call Test%run(test_getEllVolCoef_1, "test_getEllVolCoef_1")
        call Test%run(test_getEllVolCoef_2, "test_getEllVolCoef_2")
        call Test%run(test_getLowerGamma_1, "test_getLowerGamma_1")
        call Test%run(test_getLowerGamma_2, "test_getLowerGamma_2")
        call Test%run(test_getLowerGamma_3, "test_getLowerGamma_3")
        call Test%run(test_getLowerGamma_4, "test_getLowerGamma_4")
        call Test%run(test_getLowerGamma_5, "test_getLowerGamma_5")
        call Test%run(test_getUpperGamma_1, "test_getUpperGamma_1")
        call Test%run(test_getUpperGamma_2, "test_getUpperGamma_2")
        call Test%run(test_getUpperGamma_3, "test_getUpperGamma_3")
        call Test%run(test_getUpperGamma_4, "test_getUpperGamma_4")
        call Test%run(test_getUpperGamma_5, "test_getUpperGamma_5")
        call Test%run(test_getGammaSeries_1, "test_getGammaSeries_1")
        call Test%run(test_getLogSubExp_RK_1, "test_getLogSubExp_RK_1")
        call Test%run(test_getLogSumExp_RK_1, "test_getLogSumExp_RK_1")
        call Test%run(test_getLogSumExp_RK_2, "test_getLogSumExp_RK_2")
        call Test%run(test_getLogSumExp_CK_1, "test_getLogSumExp_CK_1")
        call Test%run(test_getLogSumExp_CK_2, "test_getLogSumExp_CK_2")
        call Test%run(test_getLogFactorial_1, "test_getLogFactorial_1")
        call Test%run(test_getGammaHalfInt_1, "test_getGammaHalfInt_1")
        call Test%run(test_getGammaContFrac_1, "test_getGammaContFrac_1")
        call Test%run(test_getLogEllVolCoef_1, "test_getLogEllVolCoef_1")
        call Test%run(test_getLogEllVolCoef_2, "test_getLogEllVolCoef_2")
        call Test%run(test_getLogEggBoxSD_RK_1, "test_getLogEggBoxSD_RK_1")
        call Test%run(test_getLogEggBoxSD_CK_1, "test_getLogEggBoxSD_CK_1")
        call Test%run(test_getLogEggBoxMD_RK_1, "test_getLogEggBoxMD_RK_1")
        call Test%run(test_getLogEggBoxMD_CK_1, "test_getLogEggBoxMD_CK_1")
        call Test%run(test_getLogVolUnitBall_1, "test_getLogVolUnitBall_1")
        call Test%run(test_getLogVolUnitBall_2, "test_getLogVolUnitBall_2")
        call Test%run(test_getLogGammaHalfInt_1, "test_getLogGammaHalfInt_1")
        call Test%run(test_getLogVolEllipsoid_1, "test_getLogVolEllipsoid_1")
        call Test%run(test_getLogVolEllipsoids_1, "test_getLogVolEllipsoids_1")
        call Test%run(test_getCumSumReverse_IK_1, "test_getCumSumReverse_IK_1")
        call Test%run(test_getCumSumReverse_RK_1, "test_getCumSumReverse_RK_1")
        call Test%run(test_getCorCeofFromFisherTrans_1, "test_getCorCeofFromFisherTrans_1")
        call Test%run(test_getFisherTransFromcorCoef_1, "test_getFisherTransFromcorCoef_1")
        call Test%finalize()
    end subroutine test_Math

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getDistanceSq_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical             :: assertion
        real(RK), parameter :: Point1(*) = [0._RK, 1._RK, 2._RK, 3._RK, 4._RK]
        real(RK), parameter :: Point2(*) = Point1 + 1._RK
        real(RK), parameter :: tolerance = 1.e-10_RK
        real(RK), parameter :: distanceSq_ref = norm2(Point2-Point1)**2
        real(RK)            :: distanceSq
        real(RK)            :: difference
        distanceSq = getDistanceSq(nd = size(Point1), Point1 = Point1, Point2 = Point2)
        difference = abs(distanceSq - distanceSq_ref) / distanceSq_ref
        assertion  = difference < tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "distanceSq_ref    = ", distanceSq_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "distanceSq        = ", distanceSq
            write(Test%outputUnit,"(*(g0,:,' '))") "difference        = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getDistanceSq_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getCorCeofFromFisherTrans_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical             :: assertion
        real(RK), parameter :: fisherTrans = 1.5_RK
        real(RK), parameter :: corCoef_ref = 0.905148253644866_RK ! tanh(fisherTrans)
        real(RK), parameter :: tolerance = 1.e-10_RK
        real(RK)            :: corCoef
        real(RK)            :: difference
        corCoef = getCorCeofFromFisherTrans(fisherTrans = fisherTrans)
        difference = abs(corCoef - corCoef_ref) / corCoef_ref
        assertion  = difference < tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "corCoef_ref   = ", corCoef_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "corCoef       = ", corCoef
            write(Test%outputUnit,"(*(g0,:,' '))") "difference    = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getCorCeofFromFisherTrans_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getFisherTransFromcorCoef_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical             :: assertion
        real(RK), parameter :: corCoef = 0.905148253644866_RK ! tanh(fisherTrans)
        real(RK), parameter :: fisherTrans_ref = 1.5_RK
        real(RK), parameter :: tolerance = 1.e-10_RK
        real(RK)            :: fisherTrans
        real(RK)            :: difference
        fisherTrans = getFisherTransFromcorCoef(corCoef = corCoef)
        difference = abs(fisherTrans - fisherTrans_ref) / fisherTrans_ref
        assertion  = difference < tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "fisherTrans_ref   = ", fisherTrans_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "fisherTrans       = ", fisherTrans
            write(Test%outputUnit,"(*(g0,:,' '))") "difference        = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getFisherTransFromcorCoef_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getCumSum_IK_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical                     :: assertion
        integer(IK)                 :: i
        integer(IK), parameter      :: CumSum_ref(*) = [1_IK, 3_IK, 6_IK, 10_IK, 15_IK, 21_IK, 28_IK, 36_IK, 45_IK, 55_IK]
        integer(IK), parameter      :: Vector(*) = [(i,i=1,size(CumSum_ref))]
        integer(IK), allocatable    :: CumSum(:)
        integer(IK), allocatable    :: Difference(:)
        CumSum = getCumSum(vecLen = int(size(Vector),kind=IK), Vec = Vector)
        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(Difference)) deallocate(Difference); allocate(Difference, mold = CumSum)
        Difference = abs(CumSum - CumSum_ref)
        assertion  = all(Difference == 0_IK)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "CumSum_ref    = ", CumSum_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "CumSum        = ", CumSum
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference    = ", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getCumSum_IK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getCumSum_RK_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical                 :: assertion
        integer(IK)             :: i
        real(RK), parameter     :: CumSum_ref(*) = [1._RK, 3._RK, 6._RK, 10._RK, 15._RK, 21._RK, 28._RK, 36._RK, 45._RK, 55._RK]
        real(RK), parameter     :: Vector(*) = [(real(i,kind=RK),i=1,size(CumSum_ref))]
        real(RK), parameter     :: tolerance = 1.e-14_RK
        real(RK), allocatable   :: CumSum(:)
        real(RK), allocatable   :: Difference(:)
        CumSum = getCumSum(vecLen = int(size(Vector),kind=IK), Vec = Vector)
        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(Difference)) deallocate(Difference); allocate(Difference, mold = CumSum)
        Difference = abs(CumSum - CumSum_ref)
        assertion  = all(Difference < tolerance)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "CumSum_ref    = ", CumSum_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "CumSum        = ", CumSum
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference    = ", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getCumSum_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getCumSumReverse_IK_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical                     :: assertion
        integer(IK)                 :: i
        integer(IK), parameter      :: CumSumReverse_ref(*) = [10_IK, 19_IK, 27_IK, 34_IK, 40_IK, 45_IK, 49_IK, 52_IK, 54_IK, 55_IK]
        integer(IK), parameter      :: Vector(*) = [(i,i=1,size(CumSumReverse_ref),1)]
        integer(IK), allocatable    :: CumSumReverse(:)
        integer(IK), allocatable    :: Difference(:)
        CumSumReverse = getCumSumReverse_IK(vecLen = int(size(Vector),kind=IK), Vec = Vector)
        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(Difference)) deallocate(Difference); allocate(Difference, mold = CumSumReverse)
        Difference = abs(CumSumReverse - CumSumReverse_ref)
        assertion  = all(Difference == 0_IK)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector              = ", Vector
            write(Test%outputUnit,"(*(g0,:,' '))") "CumSumReverse_ref   = ", CumSumReverse_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "CumSumReverse       = ", CumSumReverse
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference          = ", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getCumSumReverse_IK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getCumSumReverse_RK_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical                 :: assertion
        integer(IK)             :: i
        real(RK), parameter     :: CumSumReverse_ref(*) = [10._RK, 19._RK, 27._RK, 34._RK, 40._RK, 45._RK, 49._RK, 52._RK, 54._RK, 55._RK]
        real(RK), parameter     :: Vector(*) = [(real(i,kind=RK),i=1,size(CumSumReverse_ref))]
        real(RK), parameter     :: tolerance = 1.e-14_RK
        real(RK), allocatable   :: CumSumReverse(:)
        real(RK), allocatable   :: Difference(:)
        CumSumReverse = getCumSumReverse_RK(vecLen = int(size(Vector),kind=IK), Vec = Vector)
        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(Difference)) deallocate(Difference); allocate(Difference, mold = CumSumReverse)
        Difference = abs(CumSumReverse - CumSumReverse_ref)
        assertion  = all(Difference < tolerance)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "Vector              = ", Vector
            write(Test%outputUnit,"(*(g0,:,' '))") "CumSumReverse_ref   = ", CumSumReverse_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "CumSumReverse       = ", CumSumReverse
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference          = ", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getCumSumReverse_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogSubExp_RK_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical                 :: assertion
        real(RK), parameter     :: logTiny2 = log(2*tiny(1._RK))
        real(RK), parameter     :: logTiny1 = log(2._RK) + logTiny2
        real(RK), parameter     :: tolerance = 1.e-10_RK
        real(RK), parameter     :: logSubExp_ref = -707.7032713517043_RK
        real(RK)                :: logSubExp
        real(RK)                :: difference
        logSubExp = getLogSubExp_RK(logValueLarger = logTiny1, logValueSamller = logTiny2)
        difference = abs(logSubExp - logSubExp_ref)
        assertion  = difference < tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logSubExp_ref   = ", logSubExp_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logSubExp       = ", logSubExp
            write(Test%outputUnit,"(*(g0,:,' '))") "difference      = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogSubExp_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogSumExp_RK_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical                 :: assertion
        real(RK), parameter     :: LogValue(*) = [ log(0.5*huge(1._RK)), log(0.9*huge(1._RK)), log(0.1*huge(1._RK)) ]
        real(RK), parameter     :: tolerance = 1.e-10_RK
        real(RK), parameter     :: logSumExp_ref = 710.1881779865910_RK
        real(RK)                :: logSumExp
        real(RK)                :: difference
        logSumExp = getLogSumExp_RK(lenLogValue = size(LogValue), LogValue = LogValue, maxLogValue = maxval(LogValue))
        difference = abs(logSumExp - logSumExp_ref)
        assertion  = difference < tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logSumExp_ref   = ", logSumExp_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logSumExp       = ", logSumExp
            write(Test%outputUnit,"(*(g0,:,' '))") "difference      = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogSumExp_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogSumExp_RK_2() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical                 :: assertion
        real(RK), parameter     :: LogValue(*) = [ log(0.5*huge(1._RK)), log(0.9*huge(1._RK)), log(0.1*huge(1._RK)) ]
        real(RK), parameter     :: tolerance = 1.e-10_RK
        real(RK), parameter     :: logSumExp_ref = 710.1881779865910_RK
        real(RK)                :: logSumExp
        real(RK)                :: difference
        logSumExp = getLogSumExp_RK(lenLogValue = size(LogValue), LogValue = LogValue)
        difference = abs(logSumExp - logSumExp_ref)
        assertion  = difference < tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logSumExp_ref   = ", logSumExp_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logSumExp       = ", logSumExp
            write(Test%outputUnit,"(*(g0,:,' '))") "difference      = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogSumExp_RK_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogSumExp_CK_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical                 :: assertion
        complex(RK), parameter  :: LogValue(*) =    [ cmplx( log(0.5*huge(1._RK)), 0._RK, kind = RK ) &
                                                    , cmplx( log(0.9*huge(1._RK)), 0._RK, kind = RK ) &
                                                    , cmplx( log(0.1*huge(1._RK)), 0._RK, kind = RK ) &
                                                    ]
        real(RK), parameter     :: logSumExp_ref = cmplx(710.1881779865910_RK, 0._RK, RK)
        complex(RK), parameter  :: tolerance = cmplx(1.e-10_RK, 0._RK, RK)
        complex(RK)             :: logSumExp
        complex(RK)             :: difference
        logSumExp = getLogSumExp_CK ( lenLogValue = int(size(LogValue),kind=IK) &
                                    , LogValue = LogValue &
                                    , maxLogValue = cmplx( maxval(real(LogValue,kind=RK)), kind = RK ) &
                                    )
        difference = abs(logSumExp - logSumExp_ref)
        assertion  = real(difference,RK) < real(tolerance,RK)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logSumExp_ref   = ", logSumExp_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logSumExp       = ", logSumExp
            write(Test%outputUnit,"(*(g0,:,' '))") "difference      = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogSumExp_CK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogSumExp_CK_2() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical                 :: assertion
        complex(RK), parameter  :: LogValue(*) =    [ cmplx( log(0.5*huge(1._RK)), 0._RK, kind = RK ) &
                                                    , cmplx( log(0.9*huge(1._RK)), 0._RK, kind = RK ) &
                                                    , cmplx( log(0.1*huge(1._RK)), 0._RK, kind = RK ) &
                                                    ]
        real(RK), parameter     :: logSumExp_ref = cmplx(710.1881779865910_RK, 0._RK, RK)
        complex(RK), parameter  :: tolerance = cmplx(1.e-10_RK, 0._RK, RK)
        complex(RK)             :: logSumExp
        complex(RK)             :: difference
        logSumExp = getLogSumExp_CK(lenLogValue = int(size(LogValue),IK), LogValue = LogValue)
        difference = abs(logSumExp - logSumExp_ref)
        assertion  = real(difference,RK) < real(tolerance,RK)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logSumExp_ref   = ", logSumExp_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logSumExp       = ", logSumExp
            write(Test%outputUnit,"(*(g0,:,' '))") "difference      = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogSumExp_CK_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogEggBoxSD_RK_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical                 :: assertion
        real(RK), parameter     :: coef = 1._RK
        real(RK), parameter     :: point = 1.e1_RK
        real(RK), parameter     :: constant = 1._RK
        real(RK), parameter     :: exponent = 1._RK
        real(RK), parameter     :: tolerance = 1.e-10_RK
        real(RK), parameter     :: logEggBox_ref = 0.16092847092354756_RK
        real(RK)                :: logEggBox
        real(RK)                :: difference
        logEggBox = getLogEggBoxSD_RK(constant = constant, exponent = exponent, coef = coef, point = point)
        difference = abs(logEggBox - logEggBox_ref)
        assertion  = difference < tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logEggBox_ref   = ", logEggBox_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logEggBox       = ", logEggBox
            write(Test%outputUnit,"(*(g0,:,' '))") "difference      = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogEggBoxSD_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogEggBoxSD_CK_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical                 :: assertion
        complex(RK), parameter  :: coef = cmplx(1._RK, 0._RK, RK)
        complex(RK), parameter  :: point = cmplx(1.e1_RK, 0._RK, RK)
        complex(RK), parameter  :: constant = cmplx(1._RK, 0._RK, RK)
        complex(RK), parameter  :: exponent = cmplx(1._RK, 0._RK, RK)
        complex(RK), parameter  :: tolerance = cmplx(1.e-10_RK, 0._RK, RK)
        complex(RK), parameter  :: logEggBox_ref = cmplx(0.16092847092354756_RK, 0._RK, RK)
        complex(RK)             :: logEggBox
        complex(RK)             :: difference
        logEggBox = getLogEggBoxSD_CK(constant = constant, exponent = exponent, coef = coef, point = point)
        difference = abs(logEggBox - logEggBox_ref)
        assertion  = real(difference,RK) < real(tolerance,RK)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logEggBox_ref   = ", logEggBox_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logEggBox       = ", logEggBox
            write(Test%outputUnit,"(*(g0,:,' '))") "difference      = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogEggBoxSD_CK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogEggBoxMD_RK_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: coef(nd) = [1._RK, 2._RK, 3._RK]
        real(RK)    , parameter :: point(nd) = [1.e1_RK, 1.e2_RK, 1.e3_RK]
        real(RK)    , parameter :: constant = 1._RK
        real(RK)    , parameter :: exponent = 1._RK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: logEggBox_ref = 1.3988445480199623_RK
        real(RK)                :: logEggBox
        real(RK)                :: difference
        logEggBox = getLogEggBoxMD_RK(nd = nd, constant = constant, exponent = exponent, coef = coef, point = point)
        difference = abs(logEggBox - logEggBox_ref)
        assertion  = difference < tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logEggBox_ref   = ", logEggBox_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logEggBox       = ", logEggBox
            write(Test%outputUnit,"(*(g0,:,' '))") "difference      = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogEggBoxMD_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogEggBoxMD_CK_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical                 :: assertion
        integer(IK), parameter  :: nd = 3_IK
        complex(RK), parameter  :: coef(nd) = [cmplx(1._RK, 0._RK, RK), cmplx(2._RK, 0._RK, RK), cmplx(3._RK, 0._RK, RK)]
        complex(RK), parameter  :: point(nd) = [cmplx(1.e1_RK, 0._RK, RK), cmplx(1.e2_RK, 0._RK, RK), cmplx(1.e3_RK, 0._RK, RK)]
        complex(RK), parameter  :: constant = cmplx(1._RK, 0._RK, RK)
        complex(RK), parameter  :: exponent = cmplx(1._RK, 0._RK, RK)
        complex(RK), parameter  :: tolerance = cmplx(1.e-10_RK, 0._RK, RK)
        complex(RK), parameter  :: logEggBox_ref = cmplx(1.3988445480199623_RK, 0._RK, RK)
        complex(RK)             :: logEggBox
        complex(RK)             :: difference
        logEggBox = getLogEggBoxMD_CK(nd = nd, constant = constant, exponent = exponent, coef = coef, point = point)
        difference = abs(logEggBox - logEggBox_ref)
        assertion  = real(difference,RK) < real(tolerance,RK)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logEggBox_ref   = ", logEggBox_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logEggBox       = ", logEggBox
            write(Test%outputUnit,"(*(g0,:,' '))") "difference      = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogEggBoxMD_CK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getFactorial_1() result(assertion)
        use Constants_mod, only: RK, IK ! LCOV_EXCL_LINE
        implicit none
        logical                 :: assertion
        integer(IK), parameter  :: factorial_ref = 3628800_IK
        integer(IK), parameter  :: positiveInteger = 10_IK
        real(RK)                :: difference
        real(RK)                :: factorial
        factorial = getFactorial(positiveInteger = positiveInteger)
        difference = abs(factorial - factorial_ref)
        assertion  = difference == 0_IK
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "factorial_ref   = ", factorial_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "factorial       = ", factorial
            write(Test%outputUnit,"(*(g0,:,' '))") "difference      = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getFactorial_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogFactorial_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: positiveInteger = 10_IK
        real(RK)    , parameter :: logFactorial_ref = 15.10441257307552_RK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)                :: logFactorial
        real(RK)                :: difference
        logFactorial = getLogFactorial(positiveInteger = positiveInteger)
        difference = abs(logFactorial - logFactorial_ref)
        assertion  = difference < tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logFactorial_ref    = ", logFactorial_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logFactorial        = ", logFactorial
            write(Test%outputUnit,"(*(g0,:,' '))") "difference          = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogFactorial_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getEllVolCoef_1() result(assertion)
        use Constants_mod, only: RK, IK, PI
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 10_IK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: ellVolCoef_ref = PI**(0.5_RK*nd) / gamma(0.5_RK*nd+1._RK) ! 2.550164039877345_RK
        real(RK)                :: ellVolCoef
        real(RK)                :: difference
        ellVolCoef = getEllVolCoef(nd = nd)
        difference = abs(ellVolCoef - ellVolCoef_ref)
        assertion  = difference < tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "ellVolCoef_ref  = ", ellVolCoef_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "ellVolCoef      = ", ellVolCoef
            write(Test%outputUnit,"(*(g0,:,' '))") "difference      = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getEllVolCoef_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getEllVolCoef_2() result(assertion)
        use Constants_mod, only: RK, IK, PI
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 11_IK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: ellVolCoef_ref = PI**(0.5_RK*nd) / gamma(0.5_RK*nd+1._RK) ! 1.884103879389900_RK
        real(RK)                :: ellVolCoef
        real(RK)                :: difference
        !integer(IK)             :: i
        !do i = 1, 10000000
        ellVolCoef = getEllVolCoef(nd = nd)
        !end do
        difference = abs(ellVolCoef - ellVolCoef_ref)
        assertion  = difference < tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "ellVolCoef_ref  = ", ellVolCoef_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "ellVolCoef      = ", ellVolCoef
            write(Test%outputUnit,"(*(g0,:,' '))") "difference      = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getEllVolCoef_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogEllVolCoef_1() result(assertion)
        use Constants_mod, only: RK, IK, PI
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 10_IK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: logEllVolCoef_ref = log( PI**(0.5_RK*nd) / gamma(0.5_RK*nd+1._RK) ) ! 0.936157686464955_RK
        real(RK)                :: logEllVolCoef
        real(RK)                :: difference
        logEllVolCoef = getLogEllVolCoef(nd = nd)
        difference = abs(logEllVolCoef - logEllVolCoef_ref)
        assertion  = difference < tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logEllVolCoef_ref   = ", logEllVolCoef_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logEllVolCoef       = ", logEllVolCoef
            write(Test%outputUnit,"(*(g0,:,' '))") "difference          = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogEllVolCoef_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogEllVolCoef_2() result(assertion)
        use Constants_mod, only: RK, IK, PI
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 11_IK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: logEllVolCoef_ref = log( PI**(0.5_RK*nd) / gamma(0.5_RK*nd+1._RK) ) ! 0.633452312314559_RK
        real(RK)                :: logEllVolCoef
        real(RK)                :: difference
        !integer(IK)             :: i
        !do i = 1, 10000000
        logEllVolCoef = getLogEllVolCoef(nd = nd)
        !end do
        difference = abs(logEllVolCoef - logEllVolCoef_ref)
        assertion  = difference < tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logEllVolCoef_ref   = ", logEllVolCoef_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logEllVolCoef       = ", logEllVolCoef
            write(Test%outputUnit,"(*(g0,:,' '))") "difference          = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogEllVolCoef_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogVolUnitBall_1() result(assertion)
        use Constants_mod, only: RK, IK, PI
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 10_IK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: logVolUnitBall_ref = log( PI**(0.5_RK*nd) / gamma(0.5_RK*nd+1._RK) ) ! .9361576864649548_RK
        real(RK)                :: logVolUnitBall
        real(RK)                :: difference
        logVolUnitBall = getLogVolUnitBall(nd = nd)
        difference = abs(logVolUnitBall - logVolUnitBall_ref)
        assertion  = difference < tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logVolUnitBall_ref  = ", logVolUnitBall_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logVolUnitBall      = ", logVolUnitBall
            write(Test%outputUnit,"(*(g0,:,' '))") "difference          = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogVolUnitBall_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogVolUnitBall_2() result(assertion)
        use Constants_mod, only: RK, IK, PI
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 11_IK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: logVolUnitBall_ref = log( PI**(0.5_RK*nd) / gamma(0.5_RK*nd+1._RK) ) ! .6334523123145592_RK
        real(RK)                :: logVolUnitBall
        real(RK)                :: difference
        !integer(IK)             :: i
        !do i = 1, 10000000
        logVolUnitBall = getLogVolUnitBall(nd = nd)
        !end do
        difference = abs(logVolUnitBall - logVolUnitBall_ref)
        assertion  = difference < tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logVolUnitBall_ref  = ", logVolUnitBall_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logVolUnitBall      = ", logVolUnitBall
            write(Test%outputUnit,"(*(g0,:,' '))") "difference          = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogVolUnitBall_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogVolEllipsoid_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 10_IK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: logSqrtDetCovMat = 2._RK
        real(RK)    , parameter :: logVolEllipsoid_ref = 2.936157686464955_RK
        real(RK)                :: logVolEllipsoid
        real(RK)                :: difference
        logVolEllipsoid = getLogVolEllipsoid(nd = nd, logSqrtDetCovMat = logSqrtDetCovMat)
        difference = abs(logVolEllipsoid - logVolEllipsoid_ref)
        assertion  = difference < tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logVolEllipsoid_ref = ", logVolEllipsoid_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logVolEllipsoid     = ", logVolEllipsoid
            write(Test%outputUnit,"(*(g0,:,' '))") "difference          = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogVolEllipsoid_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogVolEllipsoids_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK)             :: i
        logical                 :: assertion
        integer(IK) , parameter :: nEllipsoid = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: LogSqrtDetCovMat(nEllipsoid) = [ (log(real(i,RK)), i = 1, nEllipsoid) ]
        real(RK)    , parameter :: EllipsoidVolume_ref(*) = [ 1.144729885849400_RK, 1.837877066409345_RK ]
        real(RK), allocatable   :: EllipsoidVolume(:)
        real(RK), allocatable   :: Difference(:)
        EllipsoidVolume = getLogVolEllipsoids(nd = nd, nEllipsoid = nEllipsoid, LogSqrtDetCovMat = LogSqrtDetCovMat)
        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(Difference)) deallocate(Difference); allocate(Difference, mold = EllipsoidVolume)
        Difference = abs(EllipsoidVolume - EllipsoidVolume_ref)
        assertion  = all(Difference < tolerance)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "EllipsoidVolume_ref = ", EllipsoidVolume_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "EllipsoidVolume     = ", EllipsoidVolume
            write(Test%outputUnit,"(*(g0,:,' '))") "difference          = ", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogVolEllipsoids_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief 
    !> Test [getLowerGamma](@ref math_mod::getlowergamma) with a small `tolerance` input optional argument.
    function test_getLowerGamma_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK)             :: i
        logical                 :: assertion
        integer(IK) , parameter :: ntest = 3
        real(RK)    , parameter :: tolerance = 1.e-10_RK ! 1.e-7_RK
        real(RK)    , parameter :: Exponent(ntest) = [1.0_RK, 5.0_RK, 0.5_RK]
        real(RK)    , parameter :: UpperLim(ntest) = [1.0_RK, 3.0_RK, 7.0_RK]
        real(RK)    , parameter :: LowerGamma_ref(ntest) = [ 0.632120558828558_RK, 0.184736755476228_RK, 0.999817189367018_RK ]
        real(RK)                :: LowerGamma(ntest)
        real(RK)                :: difference
        do i = 1, ntest
            LowerGamma(i) = getLowerGamma   ( exponent = Exponent(i) &
                                            , logGammaExponent = log_gamma(Exponent(i)) &
                                            , upperLim = UpperLim(i) &
                                            , tolerance = tolerance &
                                            )
            difference = 2 * abs(LowerGamma(i) - LowerGamma_ref(i)) / LowerGamma_ref(i)
            assertion = difference < tolerance
            if (Test%isDebugMode .and. .not. assertion) then
            ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "Exponent, UpperLim, Reference LowerGamma, Computed LowerGamma, difference:"
                write(Test%outputUnit,"(*(g0,:,', '))") Exponent(i), UpperLim(i), LowerGamma(i), LowerGamma_ref(i), difference
                write(Test%outputUnit,"(*(g0,:,', '))")
            end if
            ! LCOV_EXCL_STOP
        end do
    end function test_getLowerGamma_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief 
    !> Test [getLowerGamma](@ref math_mod::getlowergamma) with a medium `tolerance` input optional argument.
    function test_getLowerGamma_2() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK)             :: i
        logical                 :: assertion
        integer(IK) , parameter :: ntest = 3
        real(RK)    , parameter :: tolerance = 1.e-7_RK
        real(RK)    , parameter :: Exponent(ntest) = [1.0_RK, 5.0_RK, 0.5_RK]
        real(RK)    , parameter :: UpperLim(ntest) = [1.0_RK, 3.0_RK, 7.0_RK]
        real(RK)    , parameter :: LowerGamma_ref(ntest) = [ 0.632120558828558_RK, 0.184736755476228_RK, 0.999817189367018_RK ]
        real(RK)                :: LowerGamma(ntest)
        real(RK)                :: difference
        do i = 1, ntest
            LowerGamma(i) = getLowerGamma   ( exponent = Exponent(i) &
                                            , logGammaExponent = log_gamma(Exponent(i)) &
                                            , upperLim = UpperLim(i) &
                                            , tolerance = tolerance &
                                            )
            difference = 2 * abs(LowerGamma(i) - LowerGamma_ref(i)) / LowerGamma_ref(i)
            assertion = difference < tolerance
            if (Test%isDebugMode .and. .not. assertion) then
            ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "Exponent, UpperLim, Reference LowerGamma, Computed LowerGamma, difference:"
                write(Test%outputUnit,"(*(g0,:,', '))") Exponent(i), UpperLim(i), LowerGamma(i), LowerGamma_ref(i), difference
                write(Test%outputUnit,"(*(g0,:,', '))")
            end if
            ! LCOV_EXCL_STOP
        end do
    end function test_getLowerGamma_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief 
    !> Test [getLowerGamma](@ref math_mod::getlowergamma) with a large `tolerance` input optional argument.
    function test_getLowerGamma_3() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK)             :: i
        logical                 :: assertion
        integer(IK) , parameter :: ntest = 3
        real(RK)    , parameter :: tolerance = 1.e-3_RK
        real(RK)    , parameter :: Exponent(ntest) = [1.0_RK, 5.0_RK, 0.5_RK]
        real(RK)    , parameter :: UpperLim(ntest) = [1.0_RK, 3.0_RK, 7.0_RK]
        real(RK)    , parameter :: LowerGamma_ref(ntest) = [ 0.632120558828558_RK, 0.184736755476228_RK, 0.999817189367018_RK ]
        real(RK)                :: LowerGamma(ntest)
        real(RK)                :: difference
        do i = 1, ntest
            LowerGamma(i) = getLowerGamma   ( exponent = Exponent(i) &
                                            , logGammaExponent = log_gamma(Exponent(i)) &
                                            , upperLim = UpperLim(i) &
                                            , tolerance = tolerance &
                                            )
            difference = 2 * abs(LowerGamma(i) - LowerGamma_ref(i)) / LowerGamma_ref(i)
            assertion = difference < tolerance
            if (Test%isDebugMode .and. .not. assertion) then
            ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "Exponent, UpperLim, Reference LowerGamma, Computed LowerGamma, difference:"
                write(Test%outputUnit,"(*(g0,:,', '))") Exponent(i), UpperLim(i), LowerGamma(i), LowerGamma_ref(i), difference
                write(Test%outputUnit,"(*(g0,:,', '))")
            end if
            ! LCOV_EXCL_STOP
        end do
    end function test_getLowerGamma_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief 
    !> Test [getLowerGamma](@ref math_mod::getlowergamma) without the `tolerance` input optional argument, 
    !> in which case, the procedure should default to `epsilon` for the value of tolerance.
    function test_getLowerGamma_4() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK)             :: i
        logical                 :: assertion
        integer(IK) , parameter :: ntest = 3
        real(RK)    , parameter :: tolerance = epsilon(1._RK)
        real(RK)    , parameter :: Exponent(ntest) = [1.0_RK, 5.0_RK, 0.5_RK]
        real(RK)    , parameter :: UpperLim(ntest) = [1.0_RK, 3.0_RK, 7.0_RK]
        real(RK)    , parameter :: LowerGamma_ref(ntest) = [ 0.632120558828558_RK, 0.184736755476228_RK, 0.999817189367018_RK ]
        real(RK)                :: LowerGamma(ntest)
        real(RK)                :: difference
        do i = 1, ntest
            LowerGamma(i) = getLowerGamma   ( exponent = Exponent(i) &
                                            , logGammaExponent = log_gamma(Exponent(i)) &
                                            , upperLim = UpperLim(i) &
                                            )
            difference = 2 * abs(LowerGamma(i) - LowerGamma_ref(i)) / LowerGamma_ref(i)
            assertion = difference < 1000 * tolerance
            if (Test%isDebugMode .and. .not. assertion) then
            ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "Exponent, UpperLim, Reference LowerGamma, Computed LowerGamma, difference:"
                write(Test%outputUnit,"(*(g0,:,', '))") Exponent(i), UpperLim(i), LowerGamma(i), LowerGamma_ref(i), difference
                write(Test%outputUnit,"(*(g0,:,', '))")
            end if
            ! LCOV_EXCL_STOP
        end do
    end function test_getLowerGamma_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief 
    !> Test [getLowerGamma](@ref math_mod::getlowergamma) with a wrong positive value for the input argument `upperLim`.
    function test_getLowerGamma_5() result(assertion)
        use Constants_mod, only: RK, IK, HUGE_RK
        implicit none
        logical                 :: assertion
        real(RK)    , parameter :: exponent = +1.0_RK
        real(RK)    , parameter :: upperlim = -1.0_RK
        real(RK)    , parameter :: lowerGamma_ref = -HUGE_RK
        real(RK)                :: lowerGamma
        lowerGamma = getLowerGamma  ( exponent = exponent &
                                    , logGammaExponent = log_gamma(exponent) &
                                    , upperLim = upperLim &
                                    )
        assertion = lowerGamma == lowerGamma_ref
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Exponent, UpperLim, Reference LowerGamma, Computed LowerGamma:"
            write(Test%outputUnit,"(*(g0,:,', '))") exponent, upperlim, lowerGamma, lowerGamma_ref
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLowerGamma_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief 
    !> Test [getUpperGamma](@ref math_mod::getuppergamma) with a small `tolerance` input optional argument.
    function test_getUpperGamma_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK)             :: i
        logical                 :: assertion
        integer(IK) , parameter :: ntest = 3
        real(RK)    , parameter :: tolerance = 1.e-10_RK ! 1.e-7_RK
        real(RK)    , parameter :: Exponent(ntest) = [1.0_RK, 5.0_RK, 0.5_RK]
        real(RK)    , parameter :: LowerLim(ntest) = [1.0_RK, 3.0_RK, 7.0_RK]
        real(RK)    , parameter :: UpperGamma_ref(ntest) = [0.367879441171442_RK, 0.815263244523772_RK, 1.828106329818355e-04_RK]
        real(RK)                :: UpperGamma(ntest)
        real(RK)                :: difference
        do i = 1, ntest
            UpperGamma(i) = getUpperGamma   ( exponent = Exponent(i) &
                                            , logGammaExponent = log_gamma(Exponent(i)) &
                                            , lowerLim = LowerLim(i) &
                                            , tolerance = tolerance &
                                            )
            difference = 2 * abs(UpperGamma(i) - UpperGamma_ref(i)) / UpperGamma_ref(i)
            assertion = difference < tolerance
            if (Test%isDebugMode .and. .not. assertion) then
            ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "Exponent, UpperLim, Reference UpperGamma, Computed UpperGamma, difference:"
                write(Test%outputUnit,"(*(g0,:,', '))") Exponent(i), LowerLim(i), UpperGamma(i), UpperGamma_ref(i), difference
                write(Test%outputUnit,"(*(g0,:,', '))")
            end if
            ! LCOV_EXCL_STOP
        end do
    end function test_getUpperGamma_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief 
    !> Test [getUpperGamma](@ref math_mod::getuppergamma) with a medium `tolerance` input optional argument.
    function test_getUpperGamma_2() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK)             :: i
        logical                 :: assertion
        integer(IK) , parameter :: ntest = 3
        real(RK)    , parameter :: tolerance = 1.e-7_RK
        real(RK)    , parameter :: Exponent(ntest) = [1.0_RK, 5.0_RK, 0.5_RK]
        real(RK)    , parameter :: LowerLim(ntest) = [1.0_RK, 3.0_RK, 7.0_RK]
        real(RK)    , parameter :: UpperGamma_ref(ntest) = [0.367879441171442_RK, 0.815263244523772_RK, 1.828106329818355e-04_RK]
        real(RK)                :: UpperGamma(ntest)
        real(RK)                :: difference
        do i = 1, ntest
            UpperGamma(i) = getUpperGamma   ( exponent = Exponent(i) &
                                            , logGammaExponent = log_gamma(Exponent(i)) &
                                            , lowerLim = LowerLim(i) &
                                            , tolerance = tolerance &
                                            )
            difference = 2 * abs(UpperGamma(i) - UpperGamma_ref(i)) / UpperGamma_ref(i)
            assertion = difference < tolerance
            if (Test%isDebugMode .and. .not. assertion) then
            ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "Exponent, UpperLim, Reference UpperGamma, Computed UpperGamma, difference:"
                write(Test%outputUnit,"(*(g0,:,', '))") Exponent(i), LowerLim(i), UpperGamma(i), UpperGamma_ref(i), difference
                write(Test%outputUnit,"(*(g0,:,', '))")
            end if
            ! LCOV_EXCL_STOP
        end do
    end function test_getUpperGamma_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief 
    !> Test [getUpperGamma](@ref math_mod::getuppergamma) with a coarse `tolerance` input optional argument.
    function test_getUpperGamma_3() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK)             :: i
        logical                 :: assertion
        integer(IK) , parameter :: ntest = 3
        real(RK)    , parameter :: tolerance = 1.e-3_RK
        real(RK)    , parameter :: Exponent(ntest) = [1.0_RK, 5.0_RK, 0.5_RK]
        real(RK)    , parameter :: LowerLim(ntest) = [1.0_RK, 3.0_RK, 7.0_RK]
        real(RK)    , parameter :: UpperGamma_ref(ntest) = [0.367879441171442_RK, 0.815263244523772_RK, 1.828106329818355e-04_RK]
        real(RK)                :: UpperGamma(ntest)
        real(RK)                :: difference
        do i = 1, ntest
            UpperGamma(i) = getUpperGamma   ( exponent = Exponent(i) &
                                            , logGammaExponent = log_gamma(Exponent(i)) &
                                            , lowerLim = LowerLim(i) &
                                            , tolerance = tolerance &
                                            )
            difference = 2 * abs(UpperGamma(i) - UpperGamma_ref(i)) / UpperGamma_ref(i)
            assertion = difference < tolerance
            if (Test%isDebugMode .and. .not. assertion) then
            ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "Exponent, UpperLim, Reference UpperGamma, Computed UpperGamma, difference:"
                write(Test%outputUnit,"(*(g0,:,', '))") Exponent(i), LowerLim(i), UpperGamma(i), UpperGamma_ref(i), difference
                write(Test%outputUnit,"(*(g0,:,', '))")
            end if
            ! LCOV_EXCL_STOP
        end do
    end function test_getUpperGamma_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief 
    !> Test [getUpperGamma](@ref math_mod::getuppergamma) without the `tolerance` input optional argument, in which case
    !> the procedure should default to `epsilon` for the value of tolerance.
    function test_getUpperGamma_4() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK)             :: i
        logical                 :: assertion
        integer(IK) , parameter :: ntest = 3
        real(RK)    , parameter :: tolerance = epsilon(1._RK)
        real(RK)    , parameter :: Exponent(ntest) = [1.0_RK, 5.0_RK, 0.5_RK]
        real(RK)    , parameter :: LowerLim(ntest) = [1.0_RK, 3.0_RK, 7.0_RK]
        real(RK)    , parameter :: UpperGamma_ref(ntest) = [0.367879441171442_RK, 0.815263244523772_RK, 1.828106329818355e-04_RK]
        real(RK)                :: UpperGamma(ntest)
        real(RK)                :: difference
        do i = 1, ntest
            UpperGamma(i) = getUpperGamma   ( exponent = Exponent(i) &
                                            , logGammaExponent = log_gamma(Exponent(i)) &
                                            , lowerLim = LowerLim(i) &
                                            )
            difference = 2 * abs(UpperGamma(i) - UpperGamma_ref(i)) / UpperGamma_ref(i)
            assertion = difference < 1000 * tolerance
            if (Test%isDebugMode .and. .not. assertion) then
            ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "Exponent, UpperLim, Reference UpperGamma, Computed UpperGamma, difference:"
                write(Test%outputUnit,"(*(g0,:,', '))") Exponent(i), LowerLim(i), UpperGamma(i), UpperGamma_ref(i), difference
                write(Test%outputUnit,"(*(g0,:,', '))")
            end if
            ! LCOV_EXCL_STOP
        end do
    end function test_getUpperGamma_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief 
    !> Test [getUpperGamma](@ref math_mod::getuppergamma) with a wrong positive value for the input argument `lowerLim`.
    function test_getUpperGamma_5() result(assertion)
        use Constants_mod, only: RK, IK, HUGE_RK
        implicit none
        logical                 :: assertion
        real(RK)    , parameter :: exponent = +1.0_RK
        real(RK)    , parameter :: lowerLim = -1.0_RK
        real(RK)    , parameter :: upperGamma_ref = -HUGE_RK
        real(RK)                :: upperGamma
        upperGamma = getUpperGamma  ( exponent = exponent &
                                    , logGammaExponent = log_gamma(exponent) &
                                    , lowerLim = lowerLim &
                                    )
        assertion = upperGamma == upperGamma_ref
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "exponent       ", exponent
            write(Test%outputUnit,"(*(g0,:,', '))") "lowerLim       ", lowerLim
            write(Test%outputUnit,"(*(g0,:,', '))") "upperGamma     ", upperGamma
            write(Test%outputUnit,"(*(g0,:,', '))") "upperGamma_ref ", upperGamma_ref
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getUpperGamma_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief 
    !> Test [getGammaSeries](@ref math_mod::getgammaseries) with a zero value for the input argument `upperLim`.
    function test_getGammaSeries_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical                 :: assertion
        real(RK)    , parameter :: exponent = +1.0_RK
        real(RK)    , parameter :: upperLim = 0._RK
        real(RK)    , parameter :: tolerance = 1.e-3_RK
        real(RK)    , parameter :: gammaSeries_ref = 0._RK
        real(RK)                :: gammaSeries
        gammaSeries = getGammaSeries( exponent = exponent &
                                    , logGammaExponent = log_gamma(exponent) &
                                    , upperLim = upperLim &
                                    , tolerance = tolerance &
                                    )
        assertion = gammaSeries == gammaSeries_ref
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "exponent       ", exponent
            write(Test%outputUnit,"(*(g0,:,', '))") "upperLim       ", upperLim
            write(Test%outputUnit,"(*(g0,:,', '))") "tolerance      ", tolerance
            write(Test%outputUnit,"(*(g0,:,', '))") "gammaSeries    ", gammaSeries
            write(Test%outputUnit,"(*(g0,:,', '))") "gammaSeries_ref", gammaSeries_ref
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getGammaSeries_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief 
    !> Test [getGammaContFrac](@ref math_mod::getgammacontfrac) with a zero value for the input argument `lowerLim`.
    function test_getGammaContFrac_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical                 :: assertion
        real(RK)    , parameter :: exponent = +1.0_RK
        real(RK)    , parameter :: lowerLim = 0._RK
        real(RK)    , parameter :: tolerance = 1.e-3_RK
        real(RK)    , parameter :: gammaContFrac_ref = 1._RK
        real(RK)                :: gammaContFrac
        gammaContFrac = getGammaContFrac( exponent = exponent &
                                        , logGammaExponent = log_gamma(exponent) &
                                        , lowerLim = lowerLim &
                                        , tolerance = tolerance &
                                        )
        assertion = gammaContFrac == gammaContFrac_ref
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "exponent           ", exponent
            write(Test%outputUnit,"(*(g0,:,', '))") "lowerLim           ", lowerLim
            write(Test%outputUnit,"(*(g0,:,', '))") "tolerance          ", tolerance
            write(Test%outputUnit,"(*(g0,:,', '))") "gammaContFrac      ", gammaContFrac
            write(Test%outputUnit,"(*(g0,:,', '))") "gammaContFrac_ref  ", gammaContFrac_ref
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getGammaContFrac_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief 
    !> Test the accuracy of [getLogGammaHalfInt](@ref math_mod::getloggammahalfint).
    function test_getGammaHalfInt_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical                 :: assertion
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: positiveHalfInteger = 0.5_RK
        real(RK)    , parameter :: gammaHalfInt_ref = gamma(positiveHalfInteger)
        real(RK)                :: gammaHalfInt
        real(RK)                :: difference
        gammaHalfInt = getGammaHalfInt(positiveHalfInteger)
        difference = abs(gammaHalfInt - gammaHalfInt_ref) / abs(gammaHalfInt_ref)
        assertion = difference < tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "gammaHalfInt_ref   ", gammaHalfInt_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "gammaHalfInt       ", gammaHalfInt
            write(Test%outputUnit,"(*(g0,:,', '))") "difference         ", difference
            write(Test%outputUnit,"(*(g0,:,', '))") "tolerance          ", tolerance
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getGammaHalfInt_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief 
    !> Test the accuracy of [getLogGammaHalfInt](@ref math_mod::getloggammahalfint).
    function test_getLogGammaHalfInt_1() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical                 :: assertion
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: positiveHalfInteger = 0.5_RK
        real(RK)    , parameter :: logGammaHalfInt_ref = log_gamma(positiveHalfInteger)
        real(RK)                :: logGammaHalfInt
        real(RK)                :: difference
        logGammaHalfInt = getLogGammaHalfInt(positiveHalfInteger)
        difference = abs(logGammaHalfInt - logGammaHalfInt_ref) / abs(logGammaHalfInt_ref)
        assertion = difference < tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "logGammaHalfInt_ref", logGammaHalfInt_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "logGammaHalfInt    ", logGammaHalfInt
            write(Test%outputUnit,"(*(g0,:,', '))") "difference         ", difference
            write(Test%outputUnit,"(*(g0,:,', '))") "tolerance          ", tolerance
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogGammaHalfInt_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_Math_mod ! LCOV_EXCL_LINE