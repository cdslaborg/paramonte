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

!>  \brief This module contains tests of the module [pm_cosmicRate](@ref pm_cosmicRate).
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_cosmicRate

    use pm_cosmicRate
    use pm_err, only: err_type
    use pm_test, only: test_type, IK, LK

    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getLogRateDensityH06_D0_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getLogRateDensityH06_D0_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getLogRateDensityH06_D0_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getLogRateDensityH06_D0_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getLogRateDensityH06_D0_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getLogRateDensityL08_D0_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getLogRateDensityL08_D0_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getLogRateDensityL08_D0_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getLogRateDensityL08_D0_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getLogRateDensityL08_D0_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getLogRateDensityB10_D0_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getLogRateDensityB10_D0_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getLogRateDensityB10_D0_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getLogRateDensityB10_D0_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getLogRateDensityB10_D0_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getLogRateDensityM14_D0_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getLogRateDensityM14_D0_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getLogRateDensityM14_D0_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getLogRateDensityM14_D0_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getLogRateDensityM14_D0_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getLogRateDensityP15_D0_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getLogRateDensityP15_D0_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getLogRateDensityP15_D0_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getLogRateDensityP15_D0_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getLogRateDensityP15_D0_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getLogRateDensityM17_D0_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getLogRateDensityM17_D0_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getLogRateDensityM17_D0_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getLogRateDensityM17_D0_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getLogRateDensityM17_D0_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getLogRateDensityF18_D0_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getLogRateDensityF18_D0_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getLogRateDensityF18_D0_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getLogRateDensityF18_D0_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getLogRateDensityF18_D0_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getLogRateDensityH06_D0_RK5, SK_"test_getLogRateDensityH06_D0_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getLogRateDensityH06_D0_RK4, SK_"test_getLogRateDensityH06_D0_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getLogRateDensityH06_D0_RK3, SK_"test_getLogRateDensityH06_D0_RK3")
#endif
#if     RK3_ENABLED
        call test%run(test_getLogRateDensityH06_D0_RK3, SK_"test_getLogRateDensityH06_D0_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getLogRateDensityH06_D0_RK2, SK_"test_getLogRateDensityH06_D0_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getLogRateDensityH06_D0_RK1, SK_"test_getLogRateDensityH06_D0_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getLogRateDensityL08_D0_RK5, SK_"test_getLogRateDensityL08_D0_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getLogRateDensityL08_D0_RK4, SK_"test_getLogRateDensityL08_D0_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getLogRateDensityL08_D0_RK3, SK_"test_getLogRateDensityL08_D0_RK3")
#endif
#if     RK3_ENABLED
        call test%run(test_getLogRateDensityL08_D0_RK3, SK_"test_getLogRateDensityL08_D0_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getLogRateDensityL08_D0_RK2, SK_"test_getLogRateDensityL08_D0_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getLogRateDensityL08_D0_RK1, SK_"test_getLogRateDensityL08_D0_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getLogRateDensityB10_D0_RK5, SK_"test_getLogRateDensityB10_D0_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getLogRateDensityB10_D0_RK4, SK_"test_getLogRateDensityB10_D0_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getLogRateDensityB10_D0_RK3, SK_"test_getLogRateDensityB10_D0_RK3")
#endif
#if     RK3_ENABLED
        call test%run(test_getLogRateDensityB10_D0_RK3, SK_"test_getLogRateDensityB10_D0_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getLogRateDensityB10_D0_RK2, SK_"test_getLogRateDensityB10_D0_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getLogRateDensityB10_D0_RK1, SK_"test_getLogRateDensityB10_D0_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getLogRateDensityM14_D0_RK5, SK_"test_getLogRateDensityM14_D0_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getLogRateDensityM14_D0_RK4, SK_"test_getLogRateDensityM14_D0_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getLogRateDensityM14_D0_RK3, SK_"test_getLogRateDensityM14_D0_RK3")
#endif
#if     RK3_ENABLED
        call test%run(test_getLogRateDensityM14_D0_RK3, SK_"test_getLogRateDensityM14_D0_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getLogRateDensityM14_D0_RK2, SK_"test_getLogRateDensityM14_D0_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getLogRateDensityM14_D0_RK1, SK_"test_getLogRateDensityM14_D0_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getLogRateDensityP15_D0_RK5, SK_"test_getLogRateDensityP15_D0_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getLogRateDensityP15_D0_RK4, SK_"test_getLogRateDensityP15_D0_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getLogRateDensityP15_D0_RK3, SK_"test_getLogRateDensityP15_D0_RK3")
#endif
#if     RK3_ENABLED
        call test%run(test_getLogRateDensityP15_D0_RK3, SK_"test_getLogRateDensityP15_D0_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getLogRateDensityP15_D0_RK2, SK_"test_getLogRateDensityP15_D0_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getLogRateDensityP15_D0_RK1, SK_"test_getLogRateDensityP15_D0_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getLogRateDensityM17_D0_RK5, SK_"test_getLogRateDensityM17_D0_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getLogRateDensityM17_D0_RK4, SK_"test_getLogRateDensityM17_D0_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getLogRateDensityM17_D0_RK3, SK_"test_getLogRateDensityM17_D0_RK3")
#endif
#if     RK3_ENABLED
        call test%run(test_getLogRateDensityM17_D0_RK3, SK_"test_getLogRateDensityM17_D0_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getLogRateDensityM17_D0_RK2, SK_"test_getLogRateDensityM17_D0_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getLogRateDensityM17_D0_RK1, SK_"test_getLogRateDensityM17_D0_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getLogRateDensityF18_D0_RK5, SK_"test_getLogRateDensityF18_D0_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getLogRateDensityF18_D0_RK4, SK_"test_getLogRateDensityF18_D0_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getLogRateDensityF18_D0_RK3, SK_"test_getLogRateDensityF18_D0_RK3")
#endif
#if     RK3_ENABLED
        call test%run(test_getLogRateDensityF18_D0_RK3, SK_"test_getLogRateDensityF18_D0_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getLogRateDensityF18_D0_RK2, SK_"test_getLogRateDensityF18_D0_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getLogRateDensityF18_D0_RK1, SK_"test_getLogRateDensityF18_D0_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_cosmicRate