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
!>  This module contains tests of the module [pm_complexCompareAll](@ref pm_complexCompareAll).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_complexCompareAll

    use pm_complexCompareAll
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if CK5_ENABLED
    module function test_isallless_CK5_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallleq_CK5_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallneq_CK5_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallmeq_CK5_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallmore_CK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if CK4_ENABLED
    module function test_isallless_CK4_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallleq_CK4_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallneq_CK4_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallmeq_CK4_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallmore_CK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if CK3_ENABLED
    module function test_isallless_CK3_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallleq_CK3_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallneq_CK3_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallmeq_CK3_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallmore_CK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if CK2_ENABLED
    module function test_isallless_CK2_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallleq_CK2_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallneq_CK2_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallmeq_CK2_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallmore_CK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if CK1_ENABLED
    module function test_isallless_CK1_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallleq_CK1_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallneq_CK1_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallmeq_CK1_1() result(assertion); logical(LK) :: assertion; end function
    module function test_isallmore_CK1_1() result(assertion); logical(LK) :: assertion; end function
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
        call test%run(test_isallless_CK5_1, SK_"test_isallless_CK5_1")
        call test%run(test_isallleq_CK5_1, SK_"test_isallleq_CK5_1")
        call test%run(test_isallneq_CK5_1, SK_"test_isallneq_CK5_1")
        call test%run(test_isallmeq_CK5_1, SK_"test_isallmeq_CK5_1")
        call test%run(test_isallmore_CK5_1, SK_"test_isallmore_CK5_1")
#endif
#if CK4_ENABLED
        call test%run(test_isallless_CK4_1, SK_"test_isallless_CK4_1")
        call test%run(test_isallleq_CK4_1, SK_"test_isallleq_CK4_1")
        call test%run(test_isallneq_CK4_1, SK_"test_isallneq_CK4_1")
        call test%run(test_isallmeq_CK4_1, SK_"test_isallmeq_CK4_1")
        call test%run(test_isallmore_CK4_1, SK_"test_isallmore_CK4_1")
#endif
#if CK3_ENABLED
        call test%run(test_isallless_CK3_1, SK_"test_isallless_CK3_1")
        call test%run(test_isallleq_CK3_1, SK_"test_isallleq_CK3_1")
        call test%run(test_isallneq_CK3_1, SK_"test_isallneq_CK3_1")
        call test%run(test_isallmeq_CK3_1, SK_"test_isallmeq_CK3_1")
        call test%run(test_isallmore_CK3_1, SK_"test_isallmore_CK3_1")
#endif
#if CK2_ENABLED
        call test%run(test_isallless_CK2_1, SK_"test_isallless_CK2_1")
        call test%run(test_isallleq_CK2_1, SK_"test_isallleq_CK2_1")
        call test%run(test_isallneq_CK2_1, SK_"test_isallneq_CK2_1")
        call test%run(test_isallmeq_CK2_1, SK_"test_isallmeq_CK2_1")
        call test%run(test_isallmore_CK2_1, SK_"test_isallmore_CK2_1")
#endif
#if CK1_ENABLED
        call test%run(test_isallless_CK1_1, SK_"test_isallless_CK1_1")
        call test%run(test_isallleq_CK1_1, SK_"test_isallleq_CK1_1")
        call test%run(test_isallneq_CK1_1, SK_"test_isallneq_CK1_1")
        call test%run(test_isallmeq_CK1_1, SK_"test_isallmeq_CK1_1")
        call test%run(test_isallmore_CK1_1, SK_"test_isallmore_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_complexCompareAll ! LCOV_EXCL_LINE