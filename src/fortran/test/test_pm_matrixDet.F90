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

!>  \brief This module contains tests of the module [pm_matrixDet](@ref pm_matrixDet).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_matrixDet

    use pm_matrixDet ! LCOV_EXCL_LINE
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    use pm_kind, only: LK

    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if RK3_ENABLED
    module function test_getDet_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function  test_getDet_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function  test_getDet_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if RK3_ENABLED
    module function test_setDet_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function  test_setDet_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function  test_setDet_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if RK3_ENABLED
    module function test_getMatDetSqrt_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function  test_getMatDetSqrt_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function  test_getMatDetSqrt_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if RK3_ENABLED
    module function test_setMatDetSqrt_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function  test_setMatDetSqrt_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function  test_setMatDetSqrt_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if RK3_ENABLED
    module function test_getMatDetSqrtLog_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function  test_getMatDetSqrtLog_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function  test_getMatDetSqrtLog_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if RK3_ENABLED
    module function test_setMatDetSqrtLog_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function  test_setMatDetSqrtLog_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function  test_setMatDetSqrtLog_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED
        call test%run(test_getDet_RK3, SK_"test_getDet_RK3")
#endif
#if RK2_ENABLED
        call test%run(test_getDet_RK2, SK_"test_getDet_RK2")
#endif
#if RK1_ENABLED
        call test%run(test_getDet_RK1, SK_"test_getDet_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED
        call test%run(test_setDet_RK3, SK_"test_setDet_RK3")
#endif
#if RK2_ENABLED
        call test%run(test_setDet_RK2, SK_"test_setDet_RK2")
#endif
#if RK1_ENABLED
        call test%run(test_setDet_RK1, SK_"test_setDet_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED
        call test%run(test_getMatDetSqrt_RK3, SK_"test_getMatDetSqrt_RK3")
#endif
#if RK2_ENABLED
        call test%run(test_getMatDetSqrt_RK2, SK_"test_getMatDetSqrt_RK2")
#endif
#if RK1_ENABLED
        call test%run(test_getMatDetSqrt_RK1, SK_"test_getMatDetSqrt_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED
        call test%run(test_setMatDetSqrt_RK3, SK_"test_setMatDetSqrt_RK3")
#endif
#if RK2_ENABLED
        call test%run(test_setMatDetSqrt_RK2, SK_"test_setMatDetSqrt_RK2")
#endif
#if RK1_ENABLED
        call test%run(test_setMatDetSqrt_RK1, SK_"test_setMatDetSqrt_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED
        call test%run(test_getMatDetSqrtLog_RK3, SK_"test_getMatDetSqrtLog_RK3")
#endif
#if RK2_ENABLED
        call test%run(test_getMatDetSqrtLog_RK2, SK_"test_getMatDetSqrtLog_RK2")
#endif
#if RK1_ENABLED
        call test%run(test_getMatDetSqrtLog_RK1, SK_"test_getMatDetSqrtLog_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED
        call test%run(test_setMatDetSqrtLog_RK3, SK_"test_setMatDetSqrtLog_RK3")
#endif
#if RK2_ENABLED
        call test%run(test_setMatDetSqrtLog_RK2, SK_"test_setMatDetSqrtLog_RK2")
#endif
#if RK1_ENABLED
        call test%run(test_setMatDetSqrtLog_RK1, SK_"test_setMatDetSqrtLog_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_matrixDet ! LCOV_EXCL_LINE