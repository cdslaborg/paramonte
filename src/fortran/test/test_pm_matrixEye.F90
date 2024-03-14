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

!>  \brief This module contains tests of the module [pm_matrixInitDia](@ref pm_matrixInitDia).
!>  \author Amir Shahmoradi

module test_pm_matrixInitDia

    use pm_matrixInitDia ! LCOV_EXCL_LINE
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK4_ENABLED
    module function test_getMatEye_IK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK3_ENABLED
    module function test_getMatEye_IK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK2_ENABLED
    module function test_getMatEye_IK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK1_ENABLED
    module function test_getMatEye_IK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED
    module function test_getMatEye_CK3 () result(assertion); logical(LK) :: assertion; end function
#endif
#if CK2_ENABLED
    module function test_getMatEye_CK2  () result(assertion); logical(LK) :: assertion; end function
#endif
#if CK1_ENABLED
    module function test_getMatEye_CK1  () result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED
    module function test_getMatEye_RK3 () result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function test_getMatEye_RK2  () result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function test_getMatEye_RK1  () result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK4_ENABLED
    module function test_setMatEye_IK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK3_ENABLED
    module function test_setMatEye_IK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK2_ENABLED
    module function test_setMatEye_IK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK1_ENABLED
    module function test_setMatEye_IK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED
    module function test_setMatEye_CK3 () result(assertion); logical(LK) :: assertion; end function
#endif
#if CK2_ENABLED
    module function test_setMatEye_CK2  () result(assertion); logical(LK) :: assertion; end function
#endif
#if CK1_ENABLED
    module function test_setMatEye_CK1  () result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED
    module function test_setMatEye_RK3 () result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function test_setMatEye_RK2  () result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function test_setMatEye_RK1  () result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK4_ENABLED
        call test%run(test_getMatEye_IK4_1, SK_"test_getMatEye_IK4_1")
#endif
#if IK3_ENABLED
        call test%run(test_getMatEye_IK3_1, SK_"test_getMatEye_IK3_1")
#endif
#if IK2_ENABLED
        call test%run(test_getMatEye_IK2_1, SK_"test_getMatEye_IK2_1")
#endif
#if IK1_ENABLED
        call test%run(test_getMatEye_IK1_1, SK_"test_getMatEye_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED
        call test%run(test_getMatEye_CK3, SK_"test_getMatEye_CK3")
#endif
#if CK2_ENABLED
        call test%run(test_getMatEye_CK2 , SK_"test_getMatEye_CK2")
#endif
#if CK1_ENABLED
        call test%run(test_getMatEye_CK1, SK_"test_getMatEye_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED
        call test%run(test_getMatEye_RK3, SK_"test_getMatEye_RK3")
#endif
#if RK2_ENABLED
        call test%run(test_getMatEye_RK2 , SK_"test_getMatEye_RK2")
#endif
#if RK1_ENABLED
        call test%run(test_getMatEye_RK1, SK_"test_getMatEye_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK4_ENABLED
        call test%run(test_setMatEye_IK4_1, SK_"test_setMatEye_IK4_1")
#endif
#if IK3_ENABLED
        call test%run(test_setMatEye_IK3_1, SK_"test_setMatEye_IK3_1")
#endif
#if IK2_ENABLED
        call test%run(test_setMatEye_IK2_1, SK_"test_setMatEye_IK2_1")
#endif
#if IK1_ENABLED
        call test%run(test_setMatEye_IK1_1, SK_"test_setMatEye_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED
        call test%run(test_setMatEye_CK3, SK_"test_setMatEye_CK3")
#endif
#if CK2_ENABLED
        call test%run(test_setMatEye_CK2 , SK_"test_setMatEye_CK2")
#endif
#if CK1_ENABLED
        call test%run(test_setMatEye_CK1, SK_"test_setMatEye_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED
        call test%run(test_setMatEye_RK3, SK_"test_setMatEye_RK3")
#endif
#if RK2_ENABLED
        call test%run(test_setMatEye_RK2 , SK_"test_setMatEye_RK2")
#endif
#if RK1_ENABLED
        call test%run(test_setMatEye_RK1, SK_"test_setMatEye_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_matrixInitDia ! LCOV_EXCL_LINE