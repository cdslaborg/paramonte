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

!>  \brief This module contains tests of the module [pm_mathLogSubExp](@ref pm_mathLogSubExp).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_mathLogSubExp

    use pm_mathLogSubExp ! LCOV_EXCL_LINE
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: test_MathLogSubExp

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if CK5_ENABLED
    module function test_getLogSubExp_CK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if CK4_ENABLED
    module function test_getLogSubExp_CK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if CK3_ENABLED
    module function test_getLogSubExp_CK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if CK2_ENABLED
    module function test_getLogSubExp_CK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if CK1_ENABLED
    module function test_getLogSubExp_CK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if RK5_ENABLED
    module function test_getLogSubExp_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK4_ENABLED
    module function test_getLogSubExp_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK3_ENABLED
    module function test_getLogSubExp_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function test_getLogSubExp_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function test_getLogSubExp_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_MathLogSubExp()
        implicit none
        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
        call test%run(test_getLogSubExp_CK5_1, SK_"test_getLogSubExp_CK5_1")
#endif
#if CK4_ENABLED
        call test%run(test_getLogSubExp_CK4_1, SK_"test_getLogSubExp_CK4_1")
#endif
#if CK3_ENABLED
        call test%run(test_getLogSubExp_CK3_1, SK_"test_getLogSubExp_CK3_1")
#endif
#if CK2_ENABLED
        call test%run(test_getLogSubExp_CK2_1, SK_"test_getLogSubExp_CK2_1")
#endif
#if CK1_ENABLED
        call test%run(test_getLogSubExp_CK1_1, SK_"test_getLogSubExp_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
        call test%run(test_getLogSubExp_RK5_1, SK_"test_getLogSubExp_RK5_1")
#endif
#if RK4_ENABLED
        call test%run(test_getLogSubExp_RK4_1, SK_"test_getLogSubExp_RK4_1")
#endif
#if RK3_ENABLED
        call test%run(test_getLogSubExp_RK3_1, SK_"test_getLogSubExp_RK3_1")
#endif
#if RK2_ENABLED
        call test%run(test_getLogSubExp_RK2_1, SK_"test_getLogSubExp_RK2_1")
#endif
#if RK1_ENABLED
        call test%run(test_getLogSubExp_RK1_1, SK_"test_getLogSubExp_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    end subroutine test_MathLogSubExp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_mathLogSubExp ! LCOV_EXCL_LINE