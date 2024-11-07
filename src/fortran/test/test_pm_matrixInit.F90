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

!>  \brief This module contains tests of the module [pm_matrixInit](@ref pm_matrixInit).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_matrixInit

    use pm_matrixInit
    use pm_kind, only: IK, LK, RK
    use pm_test, only: test_type
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_getMatInitULD_D2_SK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_getMatInitULD_D2_SK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_getMatInitULD_D2_SK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_getMatInitULD_D2_SK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_getMatInitULD_D2_SK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_getMatInitULD_D2_IK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_getMatInitULD_D2_IK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_getMatInitULD_D2_IK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_getMatInitULD_D2_IK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_getMatInitULD_D2_IK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_getMatInitULD_D2_LK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_getMatInitULD_D2_LK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_getMatInitULD_D2_LK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_getMatInitULD_D2_LK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_getMatInitULD_D2_LK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_getMatInitULD_D2_CK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_getMatInitULD_D2_CK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_getMatInitULD_D2_CK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_getMatInitULD_D2_CK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_getMatInitULD_D2_CK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getMatInitULD_D2_RK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getMatInitULD_D2_RK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getMatInitULD_D2_RK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getMatInitULD_D2_RK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getMatInitULD_D2_RK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        module function test_setMatInitULD_D2_SK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setMatInitULD_D2_SK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setMatInitULD_D2_SK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setMatInitULD_D2_SK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setMatInitULD_D2_SK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        module function test_setMatInitULD_D2_IK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setMatInitULD_D2_IK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setMatInitULD_D2_IK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setMatInitULD_D2_IK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setMatInitULD_D2_IK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        module function test_setMatInitULD_D2_LK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_setMatInitULD_D2_LK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_setMatInitULD_D2_LK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_setMatInitULD_D2_LK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_setMatInitULD_D2_LK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        module function test_setMatInitULD_D2_CK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setMatInitULD_D2_CK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setMatInitULD_D2_CK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setMatInitULD_D2_CK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setMatInitULD_D2_CK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setMatInitULD_D2_RK5_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setMatInitULD_D2_RK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setMatInitULD_D2_RK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setMatInitULD_D2_RK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setMatInitULD_D2_RK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_getMatInitULD_D2_SK5_1, SK_"test_getMatInitULD_D2_SK5_1")
#endif
#if     SK4_ENABLED
        call test%run(test_getMatInitULD_D2_SK4_1, SK_"test_getMatInitULD_D2_SK4_1")
#endif
#if     SK3_ENABLED
        call test%run(test_getMatInitULD_D2_SK3_1, SK_"test_getMatInitULD_D2_SK3_1")
#endif
#if     SK2_ENABLED
        call test%run(test_getMatInitULD_D2_SK2_1, SK_"test_getMatInitULD_D2_SK2_1")
#endif
#if     SK1_ENABLED
        call test%run(test_getMatInitULD_D2_SK1_1, SK_"test_getMatInitULD_D2_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_getMatInitULD_D2_SK5_1, SK_"test_getMatInitULD_D2_SK5_1")
#endif
#if     SK4_ENABLED
        call test%run(test_getMatInitULD_D2_SK4_1, SK_"test_getMatInitULD_D2_SK4_1")
#endif
#if     SK3_ENABLED
        call test%run(test_getMatInitULD_D2_SK3_1, SK_"test_getMatInitULD_D2_SK3_1")
#endif
#if     SK2_ENABLED
        call test%run(test_getMatInitULD_D2_SK2_1, SK_"test_getMatInitULD_D2_SK2_1")
#endif
#if     SK1_ENABLED
        call test%run(test_getMatInitULD_D2_SK1_1, SK_"test_getMatInitULD_D2_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_getMatInitULD_D2_IK5_1, SK_"test_getMatInitULD_D2_IK5_1")
#endif
#if     IK4_ENABLED
        call test%run(test_getMatInitULD_D2_IK4_1, SK_"test_getMatInitULD_D2_IK4_1")
#endif
#if     IK3_ENABLED
        call test%run(test_getMatInitULD_D2_IK3_1, SK_"test_getMatInitULD_D2_IK3_1")
#endif
#if     IK2_ENABLED
        call test%run(test_getMatInitULD_D2_IK2_1, SK_"test_getMatInitULD_D2_IK2_1")
#endif
#if     IK1_ENABLED
        call test%run(test_getMatInitULD_D2_IK1_1, SK_"test_getMatInitULD_D2_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_getMatInitULD_D2_LK5_1, SK_"test_getMatInitULD_D2_LK5_1")
#endif
#if     LK4_ENABLED
        call test%run(test_getMatInitULD_D2_LK4_1, SK_"test_getMatInitULD_D2_LK4_1")
#endif
#if     LK3_ENABLED
        call test%run(test_getMatInitULD_D2_LK3_1, SK_"test_getMatInitULD_D2_LK3_1")
#endif
#if     LK2_ENABLED
        call test%run(test_getMatInitULD_D2_LK2_1, SK_"test_getMatInitULD_D2_LK2_1")
#endif
#if     LK1_ENABLED
        call test%run(test_getMatInitULD_D2_LK1_1, SK_"test_getMatInitULD_D2_LK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_getMatInitULD_D2_CK5_1, SK_"test_getMatInitULD_D2_CK5_1")
#endif
#if     CK4_ENABLED
        call test%run(test_getMatInitULD_D2_CK4_1, SK_"test_getMatInitULD_D2_CK4_1")
#endif
#if     CK3_ENABLED
        call test%run(test_getMatInitULD_D2_CK3_1, SK_"test_getMatInitULD_D2_CK3_1")
#endif
#if     CK2_ENABLED
        call test%run(test_getMatInitULD_D2_CK2_1, SK_"test_getMatInitULD_D2_CK2_1")
#endif
#if     CK1_ENABLED
        call test%run(test_getMatInitULD_D2_CK1_1, SK_"test_getMatInitULD_D2_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getMatInitULD_D2_RK5_1, SK_"test_getMatInitULD_D2_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_getMatInitULD_D2_RK4_1, SK_"test_getMatInitULD_D2_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_getMatInitULD_D2_RK3_1, SK_"test_getMatInitULD_D2_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_getMatInitULD_D2_RK2_1, SK_"test_getMatInitULD_D2_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_getMatInitULD_D2_RK1_1, SK_"test_getMatInitULD_D2_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setMatInitULD_D2_SK5_1, SK_"test_setMatInitULD_D2_SK5_1")
#endif
#if     SK4_ENABLED
        call test%run(test_setMatInitULD_D2_SK4_1, SK_"test_setMatInitULD_D2_SK4_1")
#endif
#if     SK3_ENABLED
        call test%run(test_setMatInitULD_D2_SK3_1, SK_"test_setMatInitULD_D2_SK3_1")
#endif
#if     SK2_ENABLED
        call test%run(test_setMatInitULD_D2_SK2_1, SK_"test_setMatInitULD_D2_SK2_1")
#endif
#if     SK1_ENABLED
        call test%run(test_setMatInitULD_D2_SK1_1, SK_"test_setMatInitULD_D2_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setMatInitULD_D2_SK5_1, SK_"test_setMatInitULD_D2_SK5_1")
#endif
#if     SK4_ENABLED
        call test%run(test_setMatInitULD_D2_SK4_1, SK_"test_setMatInitULD_D2_SK4_1")
#endif
#if     SK3_ENABLED
        call test%run(test_setMatInitULD_D2_SK3_1, SK_"test_setMatInitULD_D2_SK3_1")
#endif
#if     SK2_ENABLED
        call test%run(test_setMatInitULD_D2_SK2_1, SK_"test_setMatInitULD_D2_SK2_1")
#endif
#if     SK1_ENABLED
        call test%run(test_setMatInitULD_D2_SK1_1, SK_"test_setMatInitULD_D2_SK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setMatInitULD_D2_IK5_1, SK_"test_setMatInitULD_D2_IK5_1")
#endif
#if     IK4_ENABLED
        call test%run(test_setMatInitULD_D2_IK4_1, SK_"test_setMatInitULD_D2_IK4_1")
#endif
#if     IK3_ENABLED
        call test%run(test_setMatInitULD_D2_IK3_1, SK_"test_setMatInitULD_D2_IK3_1")
#endif
#if     IK2_ENABLED
        call test%run(test_setMatInitULD_D2_IK2_1, SK_"test_setMatInitULD_D2_IK2_1")
#endif
#if     IK1_ENABLED
        call test%run(test_setMatInitULD_D2_IK1_1, SK_"test_setMatInitULD_D2_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_setMatInitULD_D2_LK5_1, SK_"test_setMatInitULD_D2_LK5_1")
#endif
#if     LK4_ENABLED
        call test%run(test_setMatInitULD_D2_LK4_1, SK_"test_setMatInitULD_D2_LK4_1")
#endif
#if     LK3_ENABLED
        call test%run(test_setMatInitULD_D2_LK3_1, SK_"test_setMatInitULD_D2_LK3_1")
#endif
#if     LK2_ENABLED
        call test%run(test_setMatInitULD_D2_LK2_1, SK_"test_setMatInitULD_D2_LK2_1")
#endif
#if     LK1_ENABLED
        call test%run(test_setMatInitULD_D2_LK1_1, SK_"test_setMatInitULD_D2_LK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setMatInitULD_D2_CK5_1, SK_"test_setMatInitULD_D2_CK5_1")
#endif
#if     CK4_ENABLED
        call test%run(test_setMatInitULD_D2_CK4_1, SK_"test_setMatInitULD_D2_CK4_1")
#endif
#if     CK3_ENABLED
        call test%run(test_setMatInitULD_D2_CK3_1, SK_"test_setMatInitULD_D2_CK3_1")
#endif
#if     CK2_ENABLED
        call test%run(test_setMatInitULD_D2_CK2_1, SK_"test_setMatInitULD_D2_CK2_1")
#endif
#if     CK1_ENABLED
        call test%run(test_setMatInitULD_D2_CK1_1, SK_"test_setMatInitULD_D2_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setMatInitULD_D2_RK5_1, SK_"test_setMatInitULD_D2_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_setMatInitULD_D2_RK4_1, SK_"test_setMatInitULD_D2_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_setMatInitULD_D2_RK3_1, SK_"test_setMatInitULD_D2_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_setMatInitULD_D2_RK2_1, SK_"test_setMatInitULD_D2_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_setMatInitULD_D2_RK1_1, SK_"test_setMatInitULD_D2_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_matrixInit