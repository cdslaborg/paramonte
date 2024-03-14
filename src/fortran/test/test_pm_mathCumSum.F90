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
!>  This module contains tests of the module [pm_mathCumSum](@ref pm_mathCumSum).
!>
!>  \fintest
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 2:06 AM, September 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_mathCumSum

    use pm_mathCumSum ! LCOV_EXCL_LINE
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    use pm_kind, only: LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     IK5_ENABLED
        module function test_getCumSum_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_getCumSum_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_getCumSum_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_getCumSum_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_getCumSum_IK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     CK5_ENABLED
        module function test_getCumSum_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_getCumSum_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_getCumSum_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_getCumSum_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_getCumSum_CK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK5_ENABLED
        module function test_getCumSum_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getCumSum_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getCumSum_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getCumSum_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getCumSum_RK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     IK5_ENABLED
        module function test_setCumSum_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setCumSum_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setCumSum_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setCumSum_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setCumSum_IK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     CK5_ENABLED
        module function test_setCumSum_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setCumSum_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setCumSum_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setCumSum_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setCumSum_CK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK5_ENABLED
        module function test_setCumSum_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setCumSum_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setCumSum_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setCumSum_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setCumSum_RK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_getCumSum_IK5, SK_"test_getCumSum_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_getCumSum_IK4, SK_"test_getCumSum_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_getCumSum_IK3, SK_"test_getCumSum_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_getCumSum_IK2, SK_"test_getCumSum_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_getCumSum_IK1, SK_"test_getCumSum_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_getCumSum_CK5, SK_"test_getCumSum_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_getCumSum_CK4, SK_"test_getCumSum_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_getCumSum_CK3, SK_"test_getCumSum_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_getCumSum_CK2, SK_"test_getCumSum_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_getCumSum_CK1, SK_"test_getCumSum_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getCumSum_RK5, SK_"test_getCumSum_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getCumSum_RK4, SK_"test_getCumSum_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getCumSum_RK3, SK_"test_getCumSum_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getCumSum_RK2, SK_"test_getCumSum_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getCumSum_RK1, SK_"test_getCumSum_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setCumSum_IK5, SK_"test_setCumSum_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_setCumSum_IK4, SK_"test_setCumSum_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_setCumSum_IK3, SK_"test_setCumSum_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_setCumSum_IK2, SK_"test_setCumSum_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_setCumSum_IK1, SK_"test_setCumSum_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setCumSum_CK5, SK_"test_setCumSum_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_setCumSum_CK4, SK_"test_setCumSum_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_setCumSum_CK3, SK_"test_setCumSum_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_setCumSum_CK2, SK_"test_setCumSum_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_setCumSum_CK1, SK_"test_setCumSum_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setCumSum_RK5, SK_"test_setCumSum_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setCumSum_RK4, SK_"test_setCumSum_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setCumSum_RK3, SK_"test_setCumSum_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setCumSum_RK2, SK_"test_setCumSum_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setCumSum_RK1, SK_"test_setCumSum_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%summarize()

    end subroutine setTest

end module test_pm_mathCumSum ! LCOV_EXCL_LINE