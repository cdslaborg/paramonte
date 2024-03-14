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

!>  \brief This module contains tests of the module [pm_mathMinMax](@ref pm_mathMinMax).
!>  \author Amir Shahmoradi

module test_pm_mathMinMax

    use pm_mathMinMax ! LCOV_EXCL_LINE
    use pm_err, only: err_type
    use pm_test, only: test_type, LK

    implicit none

    private
    public :: setTest
    type(test_type) :: test

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if CK3_ENABLED
    module function test_getMinMax_CK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if CK2_ENABLED
    module function test_getMinMax_CK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if CK1_ENABLED
    module function test_getMinMax_CK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if RK3_ENABLED
    module function test_getMinMax_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function test_getMinMax_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function test_getMinMax_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if IK4_ENABLED
    module function test_getMinMax_IK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK3_ENABLED
    module function test_getMinMax_IK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK2_ENABLED
    module function test_getMinMax_IK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK1_ENABLED
    module function test_getMinMax_IK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
    module function test_getMinMax_SK_1() result(assertion); logical(LK) :: assertion; end function
    module function test_getMinMax_SSK_1() result(assertion); logical(LK) :: assertion; end function
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if CK3_ENABLED
    module function test_setMinMax_CK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if CK2_ENABLED
    module function test_setMinMax_CK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if CK1_ENABLED
    module function test_setMinMax_CK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if RK3_ENABLED
    module function test_setMinMax_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function test_setMinMax_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function test_setMinMax_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if IK4_ENABLED
    module function test_setMinMax_IK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK3_ENABLED
    module function test_setMinMax_IK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK2_ENABLED
    module function test_setMinMax_IK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK1_ENABLED
    module function test_setMinMax_IK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
    module function test_setMinMax_SK_1() result(assertion); logical(LK) :: assertion; end function
    module function test_setMinMax_SSK_1() result(assertion); logical(LK) :: assertion; end function
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED
        call test%run(test_getMinMax_CK3_1, SK_"test_getMinMax_CK3_1")
#endif
#if CK2_ENABLED
        call test%run(test_getMinMax_CK2_1, SK_"test_getMinMax_CK2_1")
#endif
#if CK1_ENABLED
        call test%run(test_getMinMax_CK1_1, SK_"test_getMinMax_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED
        call test%run(test_getMinMax_RK3_1, SK_"test_getMinMax_RK3_1")
#endif
#if RK2_ENABLED
        call test%run(test_getMinMax_RK2_1, SK_"test_getMinMax_RK2_1")
#endif
#if RK1_ENABLED
        call test%run(test_getMinMax_RK1_1, SK_"test_getMinMax_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK4_ENABLED
        call test%run(test_getMinMax_IK4_1, SK_"test_getMinMax_IK4_1")
#endif
#if IK3_ENABLED
        call test%run(test_getMinMax_IK3_1, SK_"test_getMinMax_IK3_1")
#endif
#if IK2_ENABLED
        call test%run(test_getMinMax_IK2_1, SK_"test_getMinMax_IK2_1")
#endif
#if IK1_ENABLED
        call test%run(test_getMinMax_IK1_1, SK_"test_getMinMax_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%run(test_getMinMax_SK_1, SK_"test_getMinMax_SK_1")
        call test%run(test_getMinMax_SSK_1, SK_"test_getMinMax_SSK_1")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED
        call test%run(test_setMinMax_CK3_1, SK_"test_setMinMax_CK3_1")
#endif
#if CK2_ENABLED
        call test%run(test_setMinMax_CK2_1, SK_"test_setMinMax_CK2_1")
#endif
#if CK1_ENABLED
        call test%run(test_setMinMax_CK1_1, SK_"test_setMinMax_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED
        call test%run(test_setMinMax_RK3_1, SK_"test_setMinMax_RK3_1")
#endif
#if RK2_ENABLED
        call test%run(test_setMinMax_RK2_1, SK_"test_setMinMax_RK2_1")
#endif
#if RK1_ENABLED
        call test%run(test_setMinMax_RK1_1, SK_"test_setMinMax_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK4_ENABLED
        call test%run(test_setMinMax_IK4_1, SK_"test_setMinMax_IK4_1")
#endif
#if IK3_ENABLED
        call test%run(test_setMinMax_IK3_1, SK_"test_setMinMax_IK3_1")
#endif
#if IK2_ENABLED
        call test%run(test_setMinMax_IK2_1, SK_"test_setMinMax_IK2_1")
#endif
#if IK1_ENABLED
        call test%run(test_setMinMax_IK1_1, SK_"test_setMinMax_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%run(test_setMinMax_SK_1, SK_"test_setMinMax_SK_1")
        call test%run(test_setMinMax_SSK_1, SK_"test_setMinMax_SSK_1")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%summarize()

    end subroutine setTest

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_mathMinMax ! LCOV_EXCL_LINE