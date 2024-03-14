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
!>  This module contains tests of the module [pm_option](@ref pm_option).
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_option

    use pm_option
    use pm_test, only: test_type, LK

    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function test_getOption_CK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if CK4_ENABLED
    module function test_getOption_CK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if CK3_ENABLED
    module function test_getOption_CK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if CK2_ENABLED
    module function test_getOption_CK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if CK1_ENABLED
    module function test_getOption_CK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function test_getOption_RK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if RK4_ENABLED
    module function test_getOption_RK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if RK3_ENABLED
    module function test_getOption_RK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function test_getOption_RK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function test_getOption_RK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK4_ENABLED
    module function test_getOption_IK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if IK3_ENABLED
    module function test_getOption_IK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if IK2_ENABLED
    module function test_getOption_IK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if IK1_ENABLED
    module function test_getOption_IK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function test_getOption_LK_1 () result(assertion); logical(LK) :: assertion; end function
    module function test_getOption_SK_1 () result(assertion); logical(LK) :: assertion; end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
        call test%run(test_getOption_CK5_1, SK_"test_getOption_CK5_1")
#endif
#if CK4_ENABLED
        call test%run(test_getOption_CK4_1, SK_"test_getOption_CK4_1")
#endif
#if CK3_ENABLED
        call test%run(test_getOption_CK3_1, SK_"test_getOption_CK3_1")
#endif
#if CK2_ENABLED
        call test%run(test_getOption_CK2_1, SK_"test_getOption_CK2_1")
#endif
#if CK1_ENABLED
        call test%run(test_getOption_CK1_1, SK_"test_getOption_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
        call test%run(test_getOption_RK5_1, SK_"test_getOption_RK5_1")
#endif
#if RK4_ENABLED
        call test%run(test_getOption_RK4_1, SK_"test_getOption_RK4_1")
#endif
#if RK3_ENABLED
        call test%run(test_getOption_RK3_1, SK_"test_getOption_RK3_1")
#endif
#if RK2_ENABLED
        call test%run(test_getOption_RK2_1, SK_"test_getOption_RK2_1")
#endif
#if RK1_ENABLED
        call test%run(test_getOption_RK1_1, SK_"test_getOption_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK4_ENABLED
        call test%run(test_getOption_IK4_1, SK_"test_getOption_IK4_1")
#endif
#if IK3_ENABLED
        call test%run(test_getOption_IK3_1, SK_"test_getOption_IK3_1")
#endif
#if IK2_ENABLED
        call test%run(test_getOption_IK2_1, SK_"test_getOption_IK2_1")
#endif
#if IK1_ENABLED
        call test%run(test_getOption_IK1_1  , SK_"test_getOption_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%run(test_getOption_LK_1, SK_"test_getOption_LK_1")
        call test%run(test_getOption_SK_1, SK_"test_getOption_SK_1")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_option