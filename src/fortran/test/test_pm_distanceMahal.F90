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

!>  \brief This module contains tests of the module [pm_distanceMahal](@ref pm_distanceMahal).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, March 22, 2012, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_distanceMahal

    use pm_distanceMahal
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    use pm_kind, only: LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED
    module function test_getDisMahalSq_CK3 () result(assertion); logical(LK) :: assertion; end function
#endif
#if CK2_ENABLED
    module function test_getDisMahalSq_CK2  () result(assertion); logical(LK) :: assertion; end function
#endif
#if CK1_ENABLED
    module function test_getDisMahalSq_CK1  () result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED
    module function test_getDisMahalSq_RK3 () result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function test_getDisMahalSq_RK2  () result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function test_getDisMahalSq_RK1  () result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED
        call test%run(test_getDisMahalSq_CK3, SK_"test_getDisMahalSq_CK3")
#endif
#if CK2_ENABLED
        call test%run(test_getDisMahalSq_CK2 , SK_"test_getDisMahalSq_CK2")
#endif
#if CK1_ENABLED
        call test%run(test_getDisMahalSq_CK1, SK_"test_getDisMahalSq_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED
        call test%run(test_getDisMahalSq_RK3, SK_"test_getDisMahalSq_RK3")
#endif
#if RK2_ENABLED
        call test%run(test_getDisMahalSq_RK2 , SK_"test_getDisMahalSq_RK2")
#endif
#if RK1_ENABLED
        call test%run(test_getDisMahalSq_RK1, SK_"test_getDisMahalSq_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !call test%run(test_getDisMahalSqSP_RK_1, SK_"test_getDisMahalSqSP_RK_1")
        !call test%run(test_getDisMahalSqMP_RK_1, SK_"test_getDisMahalSqMP_RK_1")
        !call test%run(test_getDisMahalSqSP_CK_1, SK_"test_getDisMahalSqSP_CK_1")
        !call test%run(test_getDisMahalSqMP_CK_1, SK_"test_getDisMahalSqMP_CK_1")
        call test%summarize()

    end subroutine setTest

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getDisMahalSqSP_RK_1() result(assertion)
!        use pm_kind, only: IK, RK
!        implicit none
!        integer(IK)             :: i
!        logical(LK)             :: assertion
!        integer(IK) , parameter :: nd = 3_IK
!        real(RK)    , parameter :: tolerance = 1.e-12_RK
!        real(RK)    , parameter :: mahalSq_ref = 180._RK
!        real(RK)    , parameter :: Point(nd) = [(real(i,RK),i=1,nd)]
!        real(RK)    , parameter :: mean(nd) = [(real(i**2+1._RK,RK),i=1,nd)]
!        real(RK)    , parameter :: invCov(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
!                                                                , 0._RK, 2._RK, 0._RK &
!                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(invCov) )
!        real(RK)                :: mahalSq
!        real(RK)                :: difference
!        mahalSq = getDisMahalSqSP_RK(nd = nd, mean = mean, invCov = invCov, Point = Point)
!        difference = abs(mahalSq - mahalSq_ref) / mahalSq_ref
!        assertion = difference <= tolerance
!
!        ! LCOV_EXCL_START
!        if (test%traceable .and. .not. assertion) then
!            write(test%disp%unit,"(*(g0,:,', '))")
!            write(test%disp%unit,"(*(g0,:,', '))") "mahalSq_ref    ", mahalSq_ref
!            write(test%disp%unit,"(*(g0,:,', '))") "mahalSq        ", mahalSq
!            write(test%disp%unit,"(*(g0,:,', '))") "difference     ", difference
!            write(test%disp%unit,"(*(g0,:,', '))")
!        end if
!        ! LCOV_EXCL_STOP
!
!    end function test_getDisMahalSqSP_RK_1
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getDisMahalSqMP_RK_1() result(assertion)
!        use pm_kind, only: IK, RK
!        implicit none
!        integer(IK)             :: i
!        logical(LK)             :: assertion
!        integer(IK) , parameter :: nd = 3_IK, np = 2_IK
!        real(RK)    , parameter :: tolerance = 1.e-12_RK
!        real(RK)    , parameter :: MahalSq_ref(np) = [180._RK, 36._RK]
!        real(RK)    , parameter :: Point(nd,np) = reshape([(real(i,RK),i=1,nd*np)], shape = shape(Point))
!        real(RK)    , parameter :: mean(nd) = [(real(i**2+1._RK,RK),i=1,nd)]
!        real(RK)    , parameter :: invCov(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
!                                                                , 0._RK, 2._RK, 0._RK &
!                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(invCov) )
!        real(RK)                :: mahalSq(np)
!        real(RK)                :: Difference(np)
!        mahalSq = getDisMahalSqMP_RK(nd = nd, np = np, mean = mean, invCov = invCov, Point = Point)
!        Difference = abs(mahalSq - MahalSq_ref) / MahalSq_ref
!        assertion = all(Difference <= tolerance)
!
!        ! LCOV_EXCL_START
!        if (test%traceable .and. .not. assertion) then
!            write(test%disp%unit,"(*(g0,:,', '))")
!            write(test%disp%unit,"(*(g0,:,', '))") "MahalSq_ref    ", MahalSq_ref
!            write(test%disp%unit,"(*(g0,:,', '))") "mahalSq        ", mahalSq
!            write(test%disp%unit,"(*(g0,:,', '))") "Difference     ", Difference
!            write(test%disp%unit,"(*(g0,:,', '))")
!        end if
!        ! LCOV_EXCL_STOP
!
!    end function test_getDisMahalSqMP_RK_1
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getDisMahalSqSP_CK_1() result(assertion)
!
!        use pm_kind, only: IK, RK, CK
!        implicit none
!        integer(IK)             :: i
!        logical(LK)             :: assertion
!        integer(IK) , parameter :: nd = 3_IK
!        real(RK)    , parameter :: tolerance = 1.e-12_RK
!        complex(CK) , parameter :: mahalSq_ref = 180._CK
!        complex(CK) , parameter :: Point(nd) = [(cmplx(i,0.,kind=RK),i=1,nd)]
!        complex(CK) , parameter :: mean(nd) = [(cmplx(i**2+1._RK,0.,kind=RK),i=1,nd)]
!        complex(CK) , parameter :: invCov(nd,nd) = cmplx(reshape([ 1._RK, 0._RK, 1._RK &
!                                                                    , 0._RK, 2._RK, 0._RK &
!                                                                    , 1._RK, 0._RK, 3._RK ], shape = shape(invCov) ), kind = RK )
!        complex(CK)             :: mahalSq
!        real(RK)                :: difference
!        mahalSq = getDisMahalSqSP_CK(nd = nd, mean = mean, invCov = invCov, Point = Point)
!        difference = abs(real(mahalSq - mahalSq_ref,RK)) / real(mahalSq_ref,RK)
!        assertion = difference <= tolerance
!
!        ! LCOV_EXCL_START
!        if (test%traceable .and. .not. assertion) then
!            write(test%disp%unit,"(*(g0,:,', '))")
!            write(test%disp%unit,"(*(g0,:,', '))") "mahalSq_ref    ", mahalSq_ref
!            write(test%disp%unit,"(*(g0,:,', '))") "mahalSq        ", mahalSq
!            write(test%disp%unit,"(*(g0,:,', '))") "difference     ", difference
!            write(test%disp%unit,"(*(g0,:,', '))")
!        end if
!        ! LCOV_EXCL_STOP
!
!    end function test_getDisMahalSqSP_CK_1
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getDisMahalSqMP_CK_1() result(assertion)
!        use pm_kind, only: IK, RK, CK
!        implicit none
!        integer(IK)             :: i
!        logical(LK)             :: assertion
!        real(RK)    , parameter :: tolerance = 1.e-12_RK
!        integer(IK) , parameter :: nd = 3_IK, np = 2_IK
!        complex(CK) , parameter :: MahalSq_ref(np) = [180._CK, 36._CK]
!        complex(CK) , parameter :: Point(nd,np) = reshape([(cmplx(i,0.,kind=RK),i=1,nd*np)], shape = shape(Point))
!        complex(CK) , parameter :: mean(nd) = [(cmplx(i**2+1._RK,0.,kind=RK),i=1,nd)]
!        complex(CK) , parameter :: invCov(nd,nd) = cmplx(reshape([ 1._RK, 0._RK, 1._RK &
!                                                                    , 0._RK, 2._RK, 0._RK &
!                                                                    , 1._RK, 0._RK, 3._RK ], shape = shape(invCov) ), kind=RK )
!
!        complex(CK)             :: mahalSq(np)
!        real(RK)                :: Difference(np)
!        mahalSq = getDisMahalSqMP_CK(nd = nd, np = np, mean = mean, invCov = invCov, Point = Point)
!        Difference = abs(real(mahalSq - MahalSq_ref,RK) / real(MahalSq_ref,RK))
!        assertion = all(Difference <= tolerance)
!
!        ! LCOV_EXCL_START
!        if (test%traceable .and. .not. assertion) then
!            write(test%disp%unit,"(*(g0,:,', '))")
!            write(test%disp%unit,"(*(g0,:,', '))") "MahalSq_ref    ", MahalSq_ref
!            write(test%disp%unit,"(*(g0,:,', '))") "mahalSq        ", mahalSq
!            write(test%disp%unit,"(*(g0,:,', '))") "Difference     ", Difference
!            write(test%disp%unit,"(*(g0,:,', '))")
!        end if
!        ! LCOV_EXCL_STOP
!
!    end function test_getDisMahalSqMP_CK_1
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_distanceMahal