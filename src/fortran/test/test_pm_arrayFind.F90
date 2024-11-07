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
!>  This module contains tests of the module [pm_arrayFind](@ref pm_arrayFind).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_arrayFind

    use pm_arrayFind
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_loc_D0_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_loc_D0_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_loc_D0_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_loc_D0_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_loc_D0_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_loc_D1_D0_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_loc_D1_D0_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_loc_D1_D0_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_loc_D1_D0_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_loc_D1_D0_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_loc_D1_D0_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_loc_D1_D0_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_loc_D1_D0_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_loc_D1_D0_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_loc_D1_D0_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_loc_D1_D0_LK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_loc_D1_D0_LK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_loc_D1_D0_LK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_loc_D1_D0_LK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_loc_D1_D0_LK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_loc_D1_D0_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_loc_D1_D0_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_loc_D1_D0_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_loc_D1_D0_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_loc_D1_D0_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_loc_D1_D0_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_loc_D1_D0_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_loc_D1_D0_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_loc_D1_D0_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_loc_D1_D0_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_loc_D1_D1_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_loc_D1_D1_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_loc_D1_D1_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_loc_D1_D1_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_loc_D1_D1_SK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_loc_D1_D1_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_loc_D1_D1_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_loc_D1_D1_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_loc_D1_D1_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_loc_D1_D1_IK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_loc_D1_D1_LK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_loc_D1_D1_LK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_loc_D1_D1_LK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_loc_D1_D1_LK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_loc_D1_D1_LK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_loc_D1_D1_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_loc_D1_D1_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_loc_D1_D1_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_loc_D1_D1_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_loc_D1_D1_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_loc_D1_D1_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_loc_D1_D1_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_loc_D1_D1_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_loc_D1_D1_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_loc_D1_D1_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_loc_D1_D0_SK5, SK_"test_loc_D1_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_loc_D0_D0_SK4, SK_"test_loc_D0_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_loc_D0_D0_SK3, SK_"test_loc_D0_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_loc_D0_D0_SK2, SK_"test_loc_D0_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_loc_D0_D0_SK1, SK_"test_loc_D0_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_loc_D1_D0_SK5, SK_"test_loc_D1_D0_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_loc_D1_D0_SK4, SK_"test_loc_D1_D0_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_loc_D1_D0_SK3, SK_"test_loc_D1_D0_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_loc_D1_D0_SK2, SK_"test_loc_D1_D0_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_loc_D1_D0_SK1, SK_"test_loc_D1_D0_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_loc_D1_D0_IK5, SK_"test_loc_D1_D0_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_loc_D1_D0_IK4, SK_"test_loc_D1_D0_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_loc_D1_D0_IK3, SK_"test_loc_D1_D0_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_loc_D1_D0_IK2, SK_"test_loc_D1_D0_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_loc_D1_D0_IK1, SK_"test_loc_D1_D0_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_loc_D1_D0_LK5, SK_"test_loc_D1_D0_LK5")
#endif
#if     LK4_ENABLED
        call test%run(test_loc_D1_D0_LK4, SK_"test_loc_D1_D0_LK4")
#endif
#if     LK3_ENABLED
        call test%run(test_loc_D1_D0_LK3, SK_"test_loc_D1_D0_LK3")
#endif
#if     LK2_ENABLED
        call test%run(test_loc_D1_D0_LK2, SK_"test_loc_D1_D0_LK2")
#endif
#if     LK1_ENABLED
        call test%run(test_loc_D1_D0_LK1, SK_"test_loc_D1_D0_LK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_loc_D1_D0_CK5, SK_"test_loc_D1_D0_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_loc_D1_D0_CK4, SK_"test_loc_D1_D0_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_loc_D1_D0_CK3, SK_"test_loc_D1_D0_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_loc_D1_D0_CK2, SK_"test_loc_D1_D0_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_loc_D1_D0_CK1, SK_"test_loc_D1_D0_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_loc_D1_D0_RK5, SK_"test_loc_D1_D0_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_loc_D1_D0_RK4, SK_"test_loc_D1_D0_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_loc_D1_D0_RK3, SK_"test_loc_D1_D0_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_loc_D1_D0_RK2, SK_"test_loc_D1_D0_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_loc_D1_D0_RK1, SK_"test_loc_D1_D0_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_loc_D1_D1_SK5, SK_"test_loc_D1_D1_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_loc_D1_D1_SK4, SK_"test_loc_D1_D1_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_loc_D1_D1_SK3, SK_"test_loc_D1_D1_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_loc_D1_D1_SK2, SK_"test_loc_D1_D1_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_loc_D1_D1_SK1, SK_"test_loc_D1_D1_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_loc_D1_D1_IK5, SK_"test_loc_D1_D1_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_loc_D1_D1_IK4, SK_"test_loc_D1_D1_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_loc_D1_D1_IK3, SK_"test_loc_D1_D1_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_loc_D1_D1_IK2, SK_"test_loc_D1_D1_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_loc_D1_D1_IK1, SK_"test_loc_D1_D1_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_loc_D1_D1_LK5, SK_"test_loc_D1_D1_LK5")
#endif
#if     LK4_ENABLED
        call test%run(test_loc_D1_D1_LK4, SK_"test_loc_D1_D1_LK4")
#endif
#if     LK3_ENABLED
        call test%run(test_loc_D1_D1_LK3, SK_"test_loc_D1_D1_LK3")
#endif
#if     LK2_ENABLED
        call test%run(test_loc_D1_D1_LK2, SK_"test_loc_D1_D1_LK2")
#endif
#if     LK1_ENABLED
        call test%run(test_loc_D1_D1_LK1, SK_"test_loc_D1_D1_LK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_loc_D1_D1_CK5, SK_"test_loc_D1_D1_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_loc_D1_D1_CK4, SK_"test_loc_D1_D1_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_loc_D1_D1_CK3, SK_"test_loc_D1_D1_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_loc_D1_D1_CK2, SK_"test_loc_D1_D1_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_loc_D1_D1_CK1, SK_"test_loc_D1_D1_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_loc_D1_D1_RK5, SK_"test_loc_D1_D1_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_loc_D1_D1_RK4, SK_"test_loc_D1_D1_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_loc_D1_D1_RK3, SK_"test_loc_D1_D1_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_loc_D1_D1_RK2, SK_"test_loc_D1_D1_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_loc_D1_D1_RK1, SK_"test_loc_D1_D1_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_arrayFind ! LCOV_EXCL_LINE