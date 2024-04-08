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
!>  This module contains tests of the module [pm_sampleShift](@ref pm_sampleShift).
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, Thursday 12:45 AM, August 20, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_sampleShift

    use pm_sampleShift
    use pm_kind, only: IK, RK
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_getShifted_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_getShifted_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_getShifted_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_getShifted_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_getShifted_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getShifted_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getShifted_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getShifted_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getShifted_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getShifted_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setShifted_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setShifted_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setShifted_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setShifted_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setShifted_CK1() result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setShifted_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setShifted_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setShifted_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setShifted_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setShifted_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

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

#if     CK5_ENABLED
        call test%run(test_getShifted_CK5, SK_"test_getShifted_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_getShifted_CK4, SK_"test_getShifted_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_getShifted_CK3, SK_"test_getShifted_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_getShifted_CK2, SK_"test_getShifted_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_getShifted_CK1, SK_"test_getShifted_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getShifted_RK5, SK_"test_getShifted_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getShifted_RK4, SK_"test_getShifted_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getShifted_RK3, SK_"test_getShifted_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getShifted_RK2, SK_"test_getShifted_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getShifted_RK1, SK_"test_getShifted_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setShifted_CK5, SK_"test_setShifted_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_setShifted_CK4, SK_"test_setShifted_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_setShifted_CK3, SK_"test_setShifted_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_setShifted_CK2, SK_"test_setShifted_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_setShifted_CK1, SK_"test_setShifted_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setShifted_RK5, SK_"test_setShifted_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setShifted_RK4, SK_"test_setShifted_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setShifted_RK3, SK_"test_setShifted_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setShifted_RK2, SK_"test_setShifted_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setShifted_RK1, SK_"test_setShifted_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_sampleShift