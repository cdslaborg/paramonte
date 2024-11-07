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
!>  This module contains tests of the module [pm_matrixChol](@ref pm_matrixChol).
!>
!>  \todo
!>  \pvhigh
!>  This test collection must be expanded to cover all interfaces.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_matrixChol

    use pm_matrixChol ! LCOV_EXCL_LINE
    use pm_test, only: test_type, LK

    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_setChoLow_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setChoLow_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setChoLow_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function  test_setChoLow_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function  test_setChoLow_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_getMatChol_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_getMatChol_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_getMatChol_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function  test_getMatChol_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function  test_getMatChol_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getMatChol_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getMatChol_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getMatChol_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function  test_getMatChol_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function  test_getMatChol_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setMatChol_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setMatChol_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setMatChol_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function  test_setMatChol_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function  test_setMatChol_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setMatChol_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setMatChol_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setMatChol_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function  test_setMatChol_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function  test_setMatChol_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setChoLow_RK5, SK_"test_setChoLow_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setChoLow_RK4, SK_"test_setChoLow_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setChoLow_RK3, SK_"test_setChoLow_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setChoLow_RK2, SK_"test_setChoLow_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setChoLow_RK1, SK_"test_setChoLow_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_getMatChol_CK5, SK_"test_getMatChol_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_getMatChol_CK4, SK_"test_getMatChol_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_getMatChol_CK3, SK_"test_getMatChol_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_getMatChol_CK2, SK_"test_getMatChol_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_getMatChol_CK1, SK_"test_getMatChol_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getMatChol_RK5, SK_"test_getMatChol_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getMatChol_RK4, SK_"test_getMatChol_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getMatChol_RK3, SK_"test_getMatChol_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getMatChol_RK2, SK_"test_getMatChol_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getMatChol_RK1, SK_"test_getMatChol_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setMatChol_CK5, SK_"test_setMatChol_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_setMatChol_CK4, SK_"test_setMatChol_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_setMatChol_CK3, SK_"test_setMatChol_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_setMatChol_CK2, SK_"test_setMatChol_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_setMatChol_CK1, SK_"test_setMatChol_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setMatChol_RK5, SK_"test_setMatChol_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setMatChol_RK4, SK_"test_setMatChol_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setMatChol_RK3, SK_"test_setMatChol_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setMatChol_RK2, SK_"test_setMatChol_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setMatChol_RK1, SK_"test_setMatChol_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_matrixChol ! LCOV_EXCL_LINE