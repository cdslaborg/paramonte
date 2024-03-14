!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This module contains tests of the module [pm_ellipsoid](@ref pm_math).
!>  \author Amir Shahmoradi

module test_pm_ellipsoid

    use pm_ellipsoid
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    use pm_math, only: PI_RK
    use pm_kind, only: LK

    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)
        call test%run(test_getEllVolCoef_1, SK_"test_getEllVolCoef_1")
        call test%run(test_getEllVolCoef_2, SK_"test_getEllVolCoef_2")
        call test%run(test_getLogEllVolCoef_1, SK_"test_getLogEllVolCoef_1")
        call test%run(test_getLogEllVolCoef_2, SK_"test_getLogEllVolCoef_2")
        call test%run(test_getLogVolUnitBall_1, SK_"test_getLogVolUnitBall_1")
        call test%run(test_getLogVolUnitBall_2, SK_"test_getLogVolUnitBall_2")
        call test%run(test_getLogVolUnitBall_3, SK_"test_getLogVolUnitBall_3")
        call test%run(test_getLogSurfUnitBall_1, SK_"test_getLogSurfUnitBall_1")
        call test%run(test_getLogVolEll_1, SK_"test_getLogVolEll_1")
        call test%run(test_getLogVolEll_2, SK_"test_getLogVolEll_2")
        call test%run(test_isInsideEllipsoid_1, SK_"test_isInsideEllipsoid_1")
        call test%run(test_getRandPointOnEllipsoid_1, SK_"test_getRandPointOnEllipsoid_1")
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getEllVolCoef_1() result(assertion)
        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 10_IK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: ellVolCoef_ref = PI_RK**(0.5_RK*nd) / gamma(0.5_RK*nd+1._RK) ! 2.550164039877345_RK
        real(RK)                :: ellVolCoef
        real(RK)                :: difference
        ellVolCoef = getEllVolCoef(nd = nd)
        difference = abs(ellVolCoef - ellVolCoef_ref)
        assertion  = difference < tolerance
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") "ellVolCoef_ref  = ", ellVolCoef_ref
            write(test%disp%unit,"(*(g0,:,' '))") "ellVolCoef      = ", ellVolCoef
            write(test%disp%unit,"(*(g0,:,' '))") "difference      = ", difference
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getEllVolCoef_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getEllVolCoef_2() result(assertion)
        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 11_IK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: ellVolCoef_ref = PI_RK**(0.5_RK*nd) / gamma(0.5_RK*nd+1._RK) ! 1.884103879389900_RK
        real(RK)                :: ellVolCoef
        real(RK)                :: difference
        !integer(IK)             :: i
        !do i = 1, 10000000
        ellVolCoef = getEllVolCoef(nd = nd)
        !end do
        difference = abs(ellVolCoef - ellVolCoef_ref)
        assertion  = difference < tolerance
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") "ellVolCoef_ref  = ", ellVolCoef_ref
            write(test%disp%unit,"(*(g0,:,' '))") "ellVolCoef      = ", ellVolCoef
            write(test%disp%unit,"(*(g0,:,' '))") "difference      = ", difference
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getEllVolCoef_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogEllVolCoef_1() result(assertion)
        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 10_IK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: logEllVolCoef_ref = log( PI_RK**(0.5_RK*nd) / gamma(0.5_RK*nd+1._RK) ) ! 0.936157686464955_RK
        real(RK)                :: logEllVolCoef
        real(RK)                :: difference
        logEllVolCoef = getLogEllVolCoef(nd = nd)
        difference = abs(logEllVolCoef - logEllVolCoef_ref)
        assertion  = difference < tolerance
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") "logEllVolCoef_ref   = ", logEllVolCoef_ref
            write(test%disp%unit,"(*(g0,:,' '))") "logEllVolCoef       = ", logEllVolCoef
            write(test%disp%unit,"(*(g0,:,' '))") "difference          = ", difference
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogEllVolCoef_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogEllVolCoef_2() result(assertion)
        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 11_IK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: logEllVolCoef_ref = log( PI_RK**(0.5_RK*nd) / gamma(0.5_RK*nd+1._RK) ) ! 0.633452312314559_RK
        real(RK)                :: logEllVolCoef
        real(RK)                :: difference
        !integer(IK)             :: i
        !do i = 1, 10000000
        logEllVolCoef = getLogEllVolCoef(nd = nd)
        !end do
        difference = abs(logEllVolCoef - logEllVolCoef_ref)
        assertion  = difference < tolerance
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") "logEllVolCoef_ref   = ", logEllVolCoef_ref
            write(test%disp%unit,"(*(g0,:,' '))") "logEllVolCoef       = ", logEllVolCoef
            write(test%disp%unit,"(*(g0,:,' '))") "difference          = ", difference
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogEllVolCoef_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogVolUnitBall_1() result(assertion)
        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 10_IK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: logVolUnitBall_ref = log( PI_RK**(0.5_RK*nd) / gamma(0.5_RK*nd+1._RK) ) ! .9361576864649548_RK
        real(RK)                :: logVolUnitBall
        real(RK)                :: difference
        logVolUnitBall = getLogVolUnitBall(nd = nd)
        difference = abs(logVolUnitBall - logVolUnitBall_ref)
        assertion  = difference < tolerance
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") "logVolUnitBall_ref  = ", logVolUnitBall_ref
            write(test%disp%unit,"(*(g0,:,' '))") "logVolUnitBall      = ", logVolUnitBall
            write(test%disp%unit,"(*(g0,:,' '))") "difference          = ", difference
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogVolUnitBall_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogVolUnitBall_2() result(assertion)
        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 11_IK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: logVolUnitBall_ref = log( PI_RK**(0.5_RK*nd) / gamma(0.5_RK*nd+1._RK) ) ! .6334523123145592_RK
        real(RK)                :: logVolUnitBall
        real(RK)                :: difference
        !integer(IK)             :: i
        !do i = 1, 10000000
        logVolUnitBall = getLogVolUnitBall(nd = nd)
        !end do
        difference = abs(logVolUnitBall - logVolUnitBall_ref)
        assertion  = difference < tolerance
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") "logVolUnitBall_ref  = ", logVolUnitBall_ref
            write(test%disp%unit,"(*(g0,:,' '))") "logVolUnitBall      = ", logVolUnitBall
            write(test%disp%unit,"(*(g0,:,' '))") "difference          = ", difference
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogVolUnitBall_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Test `getLogVolUnitBall()` for a range of values from an independent source.
    ! The values are taken from https://upload.wikimedia.org/wikipedia/commons/6/6c/Hypersphere_volume_and_surface_area_graphs.svg
    function test_getLogVolUnitBall_3() result(assertion)
        logical(LK)             :: assertion
        integer(IK)             :: nd
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: LogVolUnitBall_ref(0:6) = log([1._RK, 2._RK, PI_RK, 4*PI_RK/3, PI_RK**2/2, 8*PI_RK**2/15, PI_RK**3/6])
        real(RK)                :: logVolUnitBall
        real(RK)                :: difference
        do nd = 0_IK, 6_IK
            logVolUnitBall = getLogVolUnitBall(nd = nd)
            difference = abs(logVolUnitBall - LogVolUnitBall_ref(nd))
            assertion  = difference < tolerance
            if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,' '))")
                write(test%disp%unit,"(*(g0,:,' '))") "nd                      = ", nd
                write(test%disp%unit,"(*(g0,:,' '))") "LogVolUnitBall_ref(nd)  = ", LogVolUnitBall_ref(nd)
                write(test%disp%unit,"(*(g0,:,' '))") "logVolUnitBall          = ", logVolUnitBall
                write(test%disp%unit,"(*(g0,:,' '))") "difference              = ", difference
                write(test%disp%unit,"(*(g0,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end do
    end function test_getLogVolUnitBall_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Test `getLogSurfUnitBall()` for a range of values from an independent source.
    ! The values are taken from https://upload.wikimedia.org/wikipedia/commons/6/6c/Hypersphere_volume_and_surface_area_graphs.svg
    function test_getLogSurfUnitBall_1() result(assertion)
        logical(LK)             :: assertion
        integer(IK)             :: nd
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: SurfUnitBall_ref(0:8) = [0._RK, 2._RK, 2*PI_RK, 4*PI_RK, 2*PI_RK**2, 8*PI_RK**2/3, PI_RK**3, 16*PI_RK**3/15, PI_RK**4/3]
        real(RK)                :: surfUnitBall
        real(RK)                :: difference
        do nd = 0, 8
            surfUnitBall = exp(getLogSurfUnitBall(nd = nd))
            difference = abs(surfUnitBall - SurfUnitBall_ref(nd))
            assertion  = difference < tolerance
            if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,' '))")
                write(test%disp%unit,"(*(g0,:,' '))") "nd                      = ", nd
                write(test%disp%unit,"(*(g0,:,' '))") "SurfUnitBall_ref(nd)    = ", SurfUnitBall_ref(nd)
                write(test%disp%unit,"(*(g0,:,' '))") "surfUnitBall            = ", surfUnitBall
                write(test%disp%unit,"(*(g0,:,' '))") "difference              = ", difference
                write(test%disp%unit,"(*(g0,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end do
    end function test_getLogSurfUnitBall_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogVolEll_1() result(assertion)
        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 10_IK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        real(RK)    , parameter :: logSqrtDetCovMat = 2._RK
        real(RK)    , parameter :: logVolEll_ref = 2.936157686464955_RK
        real(RK)                :: logVolEll
        real(RK)                :: difference
        logVolEll = getLogVolEll(nd = nd, logSqrtDetCovMat = logSqrtDetCovMat)
        difference = abs(logVolEll - logVolEll_ref)
        assertion  = difference < tolerance
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") "logVolEll_ref = ", logVolEll_ref
            write(test%disp%unit,"(*(g0,:,' '))") "logVolEll     = ", logVolEll
            write(test%disp%unit,"(*(g0,:,' '))") "difference          = ", difference
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogVolEll_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test the accuracy of [pm_ellipsoid::getLogVolEll_D1()](@ref pm_ellipsoid::getLogVolEll_D1).
    function test_getLogVolEll_2() result(assertion)
        use pm_kind, only: RK, IK
        implicit none
        integer(IK)             :: i
        logical(LK)             :: assertion
        integer(IK) , parameter :: nEllipsoid = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-10_RK
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: LogSqrtDetCovMat(nEllipsoid) = [ (log(real(i,RK)), i = 1, nEllipsoid) ]
        real(RK)    , parameter :: EllipsoidVolume_ref(*) = [ 1.144729885849400_RK, 1.837877066409345_RK ]
        real(RK), allocatable   :: EllipsoidVolume(:)
        real(RK), allocatable   :: Difference(:)
        EllipsoidVolume = getLogVolEll(nd = nd, LogSqrtDetCovMat = LogSqrtDetCovMat)
        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(Difference)) deallocate(Difference); allocate(Difference, mold = EllipsoidVolume)
        Difference = abs(EllipsoidVolume - EllipsoidVolume_ref)
        assertion  = all(Difference < tolerance)
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") "EllipsoidVolume_ref = ", EllipsoidVolume_ref
            write(test%disp%unit,"(*(g0,:,' '))") "EllipsoidVolume     = ", EllipsoidVolume
            write(test%disp%unit,"(*(g0,:,' '))") "difference          = ", Difference
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogVolEll_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isInsideEllipsoid_1() result(assertion)
        use pm_domainBall, only: getUnifRand
        use pm_kind, only: IK, RK
        implicit none
        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: mean(nd) = [ 1._RK, 2._RK ]
        real(RK)    , parameter :: CovMat(nd,nd) = reshape( [ 1._RK, 0.5_RK, 0.5_RK, 1._RK ], shape = shape(CovMat) )
        real(RK)    , parameter :: choDia(nd) = [ 1._RK, 0.866025403784439_RK ]
        real(RK)    , parameter :: choLow(nd,nd) = reshape( [ 1._RK, 0.5_RK, 0.5_RK, 1._RK ], shape = shape(CovMat) )
        real(RK)    , parameter :: invCov(nd,nd) = reshape(  [ +1.333333333333333_RK, -0.666666666666667_RK &
                                                                , -0.666666666666667_RK, +1.333333333333333_RK ] &
                                                                , shape = shape(invCov) )
        real(RK)                :: X(nd), NormedPoint(nd)
        call getUnifRand(X, mean, choLow, choDia)
        NormedPoint = X - mean
        assertion = isInsideEllipsoid(nd, NormedPoint, invCov)
        NormedPoint = [-1.e2_RK, 1.e2_RK]
        assertion = assertion .and. .not. isInsideEllipsoid(nd, NormedPoint, invCov)
    end function test_isInsideEllipsoid_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRandPointOnEllipsoid_1() result(assertion)
        use pm_kind, only: IK, RK
        implicit none
        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: mean(nd) = [ 1._RK, 2._RK ]
        real(RK)    , parameter :: CovMat(nd,nd) = reshape( [ 1._RK, 0.5_RK, 0.5_RK, 1._RK ], shape = shape(CovMat) )
        real(RK)    , parameter :: choDia(nd) = [ 1._RK, 0.866025403784439_RK ]
        real(RK)    , parameter :: choLow(nd,nd) = reshape( [ 1._RK, 0.5_RK, 0.5_RK, 1._RK ], shape = shape(CovMat) )
        real(RK)    , parameter :: invCov(nd,nd) = reshape(  [ +1.333333333333333_RK, -0.666666666666667_RK &
                                                                , -0.666666666666667_RK, +1.333333333333333_RK ] &
                                                                , shape = shape(invCov) )
        real(RK)                :: X(nd), NormedPoint(nd)
        X = getRandPointOnEllipsoid(nd,mean,choLow,choDia)
        NormedPoint = X - mean
        assertion = dot_product(NormedPoint,matmul(invCov,NormedPoint)) - 1._RK < tolerance
        NormedPoint = [-1.e2_RK, 1.e2_RK]
        assertion = assertion .and. .not. isInsideEllipsoid(nd, NormedPoint, invCov)
        ! LCOV_EXCL_START
        if (test%traceable .and. .not. assertion) then
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "RandPointOnEllipsoid   ", X
            write(test%disp%unit,"(*(g0,:,', '))") "distance from center   ", dot_product(NormedPoint,matmul(invCov,NormedPoint))
            write(test%disp%unit,"(*(g0,:,', '))") "expected distance      ", 1._RK
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getRandPointOnEllipsoid_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_ellipsoid ! LCOV_EXCL_LINE