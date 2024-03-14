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

!>  \brief This module contains tests of the module [pm_except](@ref pm_except).
!>
!>  \warning
!>  The tests of this module are bound to fail with non-IEEE-compliant processors or compiler optimization flags.
!>  For example, in `-fast-math` mode of gfortran, the regular tests of `NaN` value will fail.
!>
!>  \fintest
!>
!>  \author
!>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_except

    use pm_except
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     CK5_ENABLED
        module function test_getInfPos_CK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_getInfPos_CK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_getInfPos_CK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_getInfPos_CK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_getInfPos_CK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK5_ENABLED
        module function test_getInfPos_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getInfPos_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getInfPos_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getInfPos_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getInfPos_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     CK5_ENABLED
        module function test_setInfPos_CK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setInfPos_CK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setInfPos_CK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setInfPos_CK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setInfPos_CK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK5_ENABLED
        module function test_setInfPos_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setInfPos_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setInfPos_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setInfPos_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setInfPos_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     CK5_ENABLED
        module function test_getInfNeg_CK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_getInfNeg_CK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_getInfNeg_CK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_getInfNeg_CK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_getInfNeg_CK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK5_ENABLED
        module function test_getInfNeg_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getInfNeg_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getInfNeg_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getInfNeg_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getInfNeg_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     CK5_ENABLED
        module function test_setInfNeg_CK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setInfNeg_CK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setInfNeg_CK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setInfNeg_CK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setInfNeg_CK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK5_ENABLED
        module function test_setInfNeg_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setInfNeg_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setInfNeg_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setInfNeg_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setInfNeg_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     CK5_ENABLED
        module function test_getNAN_CK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_getNAN_CK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_getNAN_CK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_getNAN_CK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_getNAN_CK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK5_ENABLED
        module function test_getNAN_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getNAN_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getNAN_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getNAN_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getNAN_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     CK5_ENABLED
        module function test_setNAN_CK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setNAN_CK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setNAN_CK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setNAN_CK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setNAN_CK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK5_ENABLED
        module function test_setNAN_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setNAN_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setNAN_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setNAN_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setNAN_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none
        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_getInfPos_CK5_1, SK_"test_getInfPos_CK5_1")
#endif
#if     CK4_ENABLED
        call test%run(test_getInfPos_CK4_1, SK_"test_getInfPos_CK4_1")
#endif
#if     CK3_ENABLED
        call test%run(test_getInfPos_CK3_1, SK_"test_getInfPos_CK3_1")
#endif
#if     CK2_ENABLED
        call test%run(test_getInfPos_CK2_1, SK_"test_getInfPos_CK2_1")
#endif
#if     CK1_ENABLED
        call test%run(test_getInfPos_CK1_1, SK_"test_getInfPos_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getInfPos_RK5_1, SK_"test_getInfPos_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_getInfPos_RK4_1, SK_"test_getInfPos_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_getInfPos_RK3_1, SK_"test_getInfPos_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_getInfPos_RK2_1, SK_"test_getInfPos_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_getInfPos_RK1_1, SK_"test_getInfPos_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setInfPos_CK5_1, SK_"test_setInfPos_CK5_1")
#endif
#if     CK4_ENABLED
        call test%run(test_setInfPos_CK4_1, SK_"test_setInfPos_CK4_1")
#endif
#if     CK3_ENABLED
        call test%run(test_setInfPos_CK3_1, SK_"test_setInfPos_CK3_1")
#endif
#if     CK2_ENABLED
        call test%run(test_setInfPos_CK2_1, SK_"test_setInfPos_CK2_1")
#endif
#if     CK1_ENABLED
        call test%run(test_setInfPos_CK1_1, SK_"test_setInfPos_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setInfPos_RK5_1, SK_"test_setInfPos_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_setInfPos_RK4_1, SK_"test_setInfPos_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_setInfPos_RK3_1, SK_"test_setInfPos_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_setInfPos_RK2_1, SK_"test_setInfPos_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_setInfPos_RK1_1, SK_"test_setInfPos_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_getInfNeg_CK5_1, SK_"test_getInfNeg_CK5_1")
#endif
#if     CK4_ENABLED
        call test%run(test_getInfNeg_CK4_1, SK_"test_getInfNeg_CK4_1")
#endif
#if     CK3_ENABLED
        call test%run(test_getInfNeg_CK3_1, SK_"test_getInfNeg_CK3_1")
#endif
#if     CK2_ENABLED
        call test%run(test_getInfNeg_CK2_1, SK_"test_getInfNeg_CK2_1")
#endif
#if     CK1_ENABLED
        call test%run(test_getInfNeg_CK1_1, SK_"test_getInfNeg_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getInfNeg_RK5_1, SK_"test_getInfNeg_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_getInfNeg_RK4_1, SK_"test_getInfNeg_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_getInfNeg_RK3_1, SK_"test_getInfNeg_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_getInfNeg_RK2_1, SK_"test_getInfNeg_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_getInfNeg_RK1_1, SK_"test_getInfNeg_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setInfNeg_CK5_1, SK_"test_setInfNeg_CK5_1")
#endif
#if     CK4_ENABLED
        call test%run(test_setInfNeg_CK4_1, SK_"test_setInfNeg_CK4_1")
#endif
#if     CK3_ENABLED
        call test%run(test_setInfNeg_CK3_1, SK_"test_setInfNeg_CK3_1")
#endif
#if     CK2_ENABLED
        call test%run(test_setInfNeg_CK2_1, SK_"test_setInfNeg_CK2_1")
#endif
#if     CK1_ENABLED
        call test%run(test_setInfNeg_CK1_1, SK_"test_setInfNeg_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setInfNeg_RK5_1, SK_"test_setInfNeg_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_setInfNeg_RK4_1, SK_"test_setInfNeg_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_setInfNeg_RK3_1, SK_"test_setInfNeg_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_setInfNeg_RK2_1, SK_"test_setInfNeg_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_setInfNeg_RK1_1, SK_"test_setInfNeg_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_getNAN_CK5_1, SK_"test_getNAN_CK5_1")
#endif
#if     CK4_ENABLED
        call test%run(test_getNAN_CK4_1, SK_"test_getNAN_CK4_1")
#endif
#if     CK3_ENABLED
        call test%run(test_getNAN_CK3_1, SK_"test_getNAN_CK3_1")
#endif
#if     CK2_ENABLED
        call test%run(test_getNAN_CK2_1, SK_"test_getNAN_CK2_1")
#endif
#if     CK1_ENABLED
        call test%run(test_getNAN_CK1_1, SK_"test_getNAN_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getNAN_RK5_1, SK_"test_getNAN_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_getNAN_RK4_1, SK_"test_getNAN_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_getNAN_RK3_1, SK_"test_getNAN_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_getNAN_RK2_1, SK_"test_getNAN_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_getNAN_RK1_1, SK_"test_getNAN_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setNAN_CK5_1, SK_"test_setNAN_CK5_1")
#endif
#if     CK4_ENABLED
        call test%run(test_setNAN_CK4_1, SK_"test_setNAN_CK4_1")
#endif
#if     CK3_ENABLED
        call test%run(test_setNAN_CK3_1, SK_"test_setNAN_CK3_1")
#endif
#if     CK2_ENABLED
        call test%run(test_setNAN_CK2_1, SK_"test_setNAN_CK2_1")
#endif
#if     CK1_ENABLED
        call test%run(test_setNAN_CK1_1, SK_"test_setNAN_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setNAN_RK5_1, SK_"test_setNAN_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_setNAN_RK4_1, SK_"test_setNAN_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_setNAN_RK3_1, SK_"test_setNAN_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_setNAN_RK2_1, SK_"test_setNAN_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_setNAN_RK1_1, SK_"test_setNAN_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_except ! LCOV_EXCL_LINE