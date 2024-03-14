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

!>  \brief This module contains tests of the module [val2int_pmod](@ref pm_val2int).
!>  \author Amir Shahmoradi

module test_pm_val2int

    use pm_val2int
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if IK4_ENABLED
    module function test_getInt64_LK_1() result(assertion); logical(LK) :: assertion; end function
    module function test_getInt64_SK_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK3_ENABLED
    module function test_getInt32_LK_1() result(assertion); logical(LK) :: assertion; end function
    module function test_getInt32_SK_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK2_ENABLED
    module function test_getInt16_LK_1() result(assertion); logical(LK) :: assertion; end function
    module function test_getInt16_SK_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK1_ENABLED
    module function test_getInt8_LK_1() result(assertion); logical(LK) :: assertion; end function
    module function test_getInt8_SK_1() result(assertion); logical(LK) :: assertion; end function
#endif
    module function test_getInt_LK_1() result(assertion); logical(LK) :: assertion; end function
    module function test_getInt_SK_1() result(assertion); logical(LK) :: assertion; end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()
        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK4_ENABLED
        call test%run(test_getInt64_LK_1, SK_"test_getInt64_LK_1")
        call test%run(test_getInt64_SK_1, SK_"test_getInt64_SK_1")
#endif
#if IK3_ENABLED
        call test%run(test_getInt32_LK_1, SK_"test_getInt32_LK_1")
        call test%run(test_getInt32_SK_1, SK_"test_getInt32_SK_1")
#endif
#if IK2_ENABLED
        call test%run(test_getInt16_LK_1, SK_"test_getInt16_LK_1")
        call test%run(test_getInt16_SK_1, SK_"test_getInt16_SK_1")
#endif
#if IK1_ENABLED
        call test%run(test_getInt8_LK_1, SK_"test_getInt8_LK_1")
        call test%run(test_getInt8_SK_1, SK_"test_getInt8_SK_1")
#endif
        call test%run(test_getInt_LK_1, SK_"test_getInt_LK_1")
        call test%run(test_getInt_SK_1, SK_"test_getInt_SK_1")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_val2int ! LCOV_EXCL_LINE