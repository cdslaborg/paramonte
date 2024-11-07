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

!>  \brief This module contains tests of the module [pm_arrayRange](@ref pm_arrayRange).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_arrayRange

    use pm_arrayRange
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if SK5_ENABLED
    module function test_getRange_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if SK4_ENABLED
    module function test_getRange_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if SK3_ENABLED
    module function test_getRange_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if SK2_ENABLED
    module function test_getRange_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if SK1_ENABLED
    module function test_getRange_SK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if IK5_ENABLED
    module function test_getRange_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK4_ENABLED
    module function test_getRange_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK3_ENABLED
    module function test_getRange_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK2_ENABLED
    module function test_getRange_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK1_ENABLED
    module function test_getRange_IK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if RK5_ENABLED
    module function test_getRange_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK4_ENABLED
    module function test_getRange_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK3_ENABLED
    module function test_getRange_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function test_getRange_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function test_getRange_RK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if SK5_ENABLED
    module function test_setRange_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if SK4_ENABLED
    module function test_setRange_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if SK3_ENABLED
    module function test_setRange_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if SK2_ENABLED
    module function test_setRange_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if SK1_ENABLED
    module function test_setRange_SK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if IK5_ENABLED
    module function test_setRange_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK4_ENABLED
    module function test_setRange_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK3_ENABLED
    module function test_setRange_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK2_ENABLED
    module function test_setRange_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if IK1_ENABLED
    module function test_setRange_IK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if RK5_ENABLED
    module function test_setRange_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK4_ENABLED
    module function test_setRange_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK3_ENABLED
    module function test_setRange_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function test_setRange_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function test_setRange_RK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none
        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
        call test%run(test_getRange_SK5, SK_"test_getRange_SK5")
#endif
#if SK4_ENABLED
        call test%run(test_getRange_SK4, SK_"test_getRange_SK4")
#endif
#if SK3_ENABLED
        call test%run(test_getRange_SK3, SK_"test_getRange_SK3")
#endif
#if SK2_ENABLED
        call test%run(test_getRange_SK2, SK_"test_getRange_SK2")
#endif
#if SK1_ENABLED
        call test%run(test_getRange_SK1, SK_"test_getRange_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
        call test%run(test_getRange_IK5, SK_"test_getRange_IK5")
#endif
#if IK4_ENABLED
        call test%run(test_getRange_IK4, SK_"test_getRange_IK4")
#endif
#if IK3_ENABLED
        call test%run(test_getRange_IK3, SK_"test_getRange_IK3")
#endif
#if IK2_ENABLED
        call test%run(test_getRange_IK2, SK_"test_getRange_IK2")
#endif
#if IK1_ENABLED
        call test%run(test_getRange_IK1, SK_"test_getRange_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
        call test%run(test_getRange_RK5, SK_"test_getRange_RK5")
#endif
#if RK4_ENABLED
        call test%run(test_getRange_RK4, SK_"test_getRange_RK4")
#endif
#if RK3_ENABLED
        call test%run(test_getRange_RK3, SK_"test_getRange_RK3")
#endif
#if RK2_ENABLED
        call test%run(test_getRange_RK2, SK_"test_getRange_RK2")
#endif
#if RK1_ENABLED
        call test%run(test_getRange_RK1, SK_"test_getRange_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
        call test%run(test_setRange_SK5, SK_"test_setRange_SK5")
#endif
#if SK4_ENABLED
        call test%run(test_setRange_SK4, SK_"test_setRange_SK4")
#endif
#if SK3_ENABLED
        call test%run(test_setRange_SK3, SK_"test_setRange_SK3")
#endif
#if SK2_ENABLED
        call test%run(test_setRange_SK2, SK_"test_setRange_SK2")
#endif
#if SK1_ENABLED
        call test%run(test_setRange_SK1, SK_"test_setRange_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
        call test%run(test_setRange_IK5, SK_"test_setRange_IK5")
#endif
#if IK4_ENABLED
        call test%run(test_setRange_IK4, SK_"test_setRange_IK4")
#endif
#if IK3_ENABLED
        call test%run(test_setRange_IK3, SK_"test_setRange_IK3")
#endif
#if IK2_ENABLED
        call test%run(test_setRange_IK2, SK_"test_setRange_IK2")
#endif
#if IK1_ENABLED
        call test%run(test_setRange_IK1, SK_"test_setRange_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
        call test%run(test_setRange_RK5, SK_"test_setRange_RK5")
#endif
#if RK4_ENABLED
        call test%run(test_setRange_RK4, SK_"test_setRange_RK4")
#endif
#if RK3_ENABLED
        call test%run(test_setRange_RK3, SK_"test_setRange_RK3")
#endif
#if RK2_ENABLED
        call test%run(test_setRange_RK2, SK_"test_setRange_RK2")
#endif
#if RK1_ENABLED
        call test%run(test_setRange_RK1, SK_"test_setRange_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_arrayRange ! LCOV_EXCL_LINE