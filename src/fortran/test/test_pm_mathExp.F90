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

!>  \brief
!>  This module contains tests of the module [pm_mathExp](@ref pm_mathExp).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 25, 2015, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_mathExp

    use pm_test, only: test_type
    use pm_mathExp
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     IK5_ENABLED
        module function test_getExpNext_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_getExpNext_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_getExpNext_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_getExpNext_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_getExpNext_IK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK5_ENABLED
        module function test_getExpNext_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getExpNext_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getExpNext_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getExpNext_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getExpNext_RK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     IK5_ENABLED
        module function test_getExpPrev_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_getExpPrev_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_getExpPrev_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_getExpPrev_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_getExpPrev_IK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK5_ENABLED
        module function test_getExpPrev_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getExpPrev_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getExpPrev_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getExpPrev_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getExpPrev_RK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none
        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_getExpNext_IK5, SK_"test_getExpNext_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_getExpNext_IK4, SK_"test_getExpNext_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_getExpNext_IK3, SK_"test_getExpNext_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_getExpNext_IK2, SK_"test_getExpNext_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_getExpNext_IK1, SK_"test_getExpNext_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getExpNext_RK5, SK_"test_getExpNext_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getExpNext_RK4, SK_"test_getExpNext_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getExpNext_RK3, SK_"test_getExpNext_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getExpNext_RK2, SK_"test_getExpNext_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getExpNext_RK1, SK_"test_getExpNext_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_getExpPrev_IK5, SK_"test_getExpPrev_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_getExpPrev_IK4, SK_"test_getExpPrev_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_getExpPrev_IK3, SK_"test_getExpPrev_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_getExpPrev_IK2, SK_"test_getExpPrev_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_getExpPrev_IK1, SK_"test_getExpPrev_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getExpPrev_RK5, SK_"test_getExpPrev_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getExpPrev_RK4, SK_"test_getExpPrev_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getExpPrev_RK3, SK_"test_getExpPrev_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getExpPrev_RK2, SK_"test_getExpPrev_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getExpPrev_RK1, SK_"test_getExpPrev_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_mathExp ! LCOV_EXCL_LINE