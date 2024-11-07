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
!>  This module contains tests of the module [pm_complexCompareAny](@ref pm_complexCompareAny).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_complexCompareAny

    use pm_complexCompareAny
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if CK5_ENABLED
    module function test_isanyless_CK5_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanyleq_CK5_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanyneq_CK5_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanyeq_CK5_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanymeq_CK5_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanymore_CK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if CK4_ENABLED
    module function test_isanyless_CK4_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanyleq_CK4_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanyneq_CK4_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanyeq_CK4_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanymeq_CK4_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanymore_CK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if CK3_ENABLED
    module function test_isanyless_CK3_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanyleq_CK3_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanyneq_CK3_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanyeq_CK3_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanymeq_CK3_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanymore_CK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if CK2_ENABLED
    module function test_isanyless_CK2_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanyleq_CK2_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanyneq_CK2_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanyeq_CK2_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanymeq_CK2_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanymore_CK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if CK1_ENABLED
    module function test_isanyless_CK1_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanyleq_CK1_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanyneq_CK1_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanyeq_CK1_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanymeq_CK1_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isanymore_CK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none
        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
        call test%run(test_isanyless_CK5_1, SK_"test_isanyless_CK5_1")
        call test%run(test_isanyleq_CK5_1, SK_"test_isanyleq_CK5_1")
        call test%run(test_isanyneq_CK5_1, SK_"test_isanyneq_CK5_1")
        call test%run(test_isanyeq_CK5_1, SK_"test_isanyeq_CK5_1")
        call test%run(test_isanymeq_CK5_1, SK_"test_isanymeq_CK5_1")
        call test%run(test_isanymore_CK5_1, SK_"test_isanymore_CK5_1")
#endif
#if CK4_ENABLED
        call test%run(test_isanyless_CK4_1, SK_"test_isanyless_CK4_1")
        call test%run(test_isanyleq_CK4_1, SK_"test_isanyleq_CK4_1")
        call test%run(test_isanyneq_CK4_1, SK_"test_isanyneq_CK4_1")
        call test%run(test_isanyeq_CK4_1, SK_"test_isanyeq_CK4_1")
        call test%run(test_isanymeq_CK4_1, SK_"test_isanymeq_CK4_1")
        call test%run(test_isanymore_CK4_1, SK_"test_isanymore_CK4_1")
#endif
#if CK3_ENABLED
        call test%run(test_isanyless_CK3_1, SK_"test_isanyless_CK3_1")
        call test%run(test_isanyleq_CK3_1, SK_"test_isanyleq_CK3_1")
        call test%run(test_isanyneq_CK3_1, SK_"test_isanyneq_CK3_1")
        call test%run(test_isanyeq_CK3_1, SK_"test_isanyeq_CK3_1")
        call test%run(test_isanymeq_CK3_1, SK_"test_isanymeq_CK3_1")
        call test%run(test_isanymore_CK3_1, SK_"test_isanymore_CK3_1")
#endif
#if CK2_ENABLED
        call test%run(test_isanyless_CK2_1, SK_"test_isanyless_CK2_1")
        call test%run(test_isanyleq_CK2_1, SK_"test_isanyleq_CK2_1")
        call test%run(test_isanyneq_CK2_1, SK_"test_isanyneq_CK2_1")
        call test%run(test_isanyeq_CK2_1, SK_"test_isanyeq_CK2_1")
        call test%run(test_isanymeq_CK2_1, SK_"test_isanymeq_CK2_1")
        call test%run(test_isanymore_CK2_1, SK_"test_isanymore_CK2_1")
#endif
#if CK1_ENABLED
        call test%run(test_isanyless_CK1_1, SK_"test_isanyless_CK1_1")
        call test%run(test_isanyleq_CK1_1, SK_"test_isanyleq_CK1_1")
        call test%run(test_isanyneq_CK1_1, SK_"test_isanyneq_CK1_1")
        call test%run(test_isanyeq_CK1_1, SK_"test_isanyeq_CK1_1")
        call test%run(test_isanymeq_CK1_1, SK_"test_isanymeq_CK1_1")
        call test%run(test_isanymore_CK1_1, SK_"test_isanymore_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_complexCompareAny ! LCOV_EXCL_LINE