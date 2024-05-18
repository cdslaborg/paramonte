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
!>  This include file contains procedure implementations of the tests of [pm_distBern](@ref pm_distBern).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 2:06 AM, September 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_distBern

    use pm_distBern
    use pm_err, only: err_type
    use pm_test, only: test_type, LK

    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK5_ENABLED
        module function test_isHead_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_isHead_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_isHead_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_isHead_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_isHead_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK5_ENABLED
        module function test_getBernRand_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getBernRand_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getBernRand_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getBernRand_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getBernRand_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED && RK5_ENABLED
        module function test_setBernRand_IK5_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK5_ENABLED && RK4_ENABLED
        module function test_setBernRand_IK5_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK5_ENABLED && RK3_ENABLED
        module function test_setBernRand_IK5_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK5_ENABLED && RK2_ENABLED
        module function test_setBernRand_IK5_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK5_ENABLED && RK1_ENABLED
        module function test_setBernRand_IK5_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK4_ENABLED && RK5_ENABLED
        module function test_setBernRand_IK4_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED && RK4_ENABLED
        module function test_setBernRand_IK4_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED && RK3_ENABLED
        module function test_setBernRand_IK4_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED && RK2_ENABLED
        module function test_setBernRand_IK4_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED && RK1_ENABLED
        module function test_setBernRand_IK4_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK3_ENABLED && RK5_ENABLED
        module function test_setBernRand_IK3_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED && RK4_ENABLED
        module function test_setBernRand_IK3_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED && RK3_ENABLED
        module function test_setBernRand_IK3_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED && RK2_ENABLED
        module function test_setBernRand_IK3_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED && RK1_ENABLED
        module function test_setBernRand_IK3_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK2_ENABLED && RK5_ENABLED
        module function test_setBernRand_IK2_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED && RK4_ENABLED
        module function test_setBernRand_IK2_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED && RK3_ENABLED
        module function test_setBernRand_IK2_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED && RK2_ENABLED
        module function test_setBernRand_IK2_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED && RK1_ENABLED
        module function test_setBernRand_IK2_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK1_ENABLED && RK5_ENABLED
        module function test_setBernRand_IK1_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED && RK4_ENABLED
        module function test_setBernRand_IK1_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED && RK3_ENABLED
        module function test_setBernRand_IK1_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED && RK2_ENABLED
        module function test_setBernRand_IK1_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED && RK1_ENABLED
        module function test_setBernRand_IK1_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED && RK5_ENABLED
        module function test_setBernRand_LK5_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK5_ENABLED && RK4_ENABLED
        module function test_setBernRand_LK5_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK5_ENABLED && RK3_ENABLED
        module function test_setBernRand_LK5_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK5_ENABLED && RK2_ENABLED
        module function test_setBernRand_LK5_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK5_ENABLED && RK1_ENABLED
        module function test_setBernRand_LK5_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK4_ENABLED && RK5_ENABLED
        module function test_setBernRand_LK4_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED && RK4_ENABLED
        module function test_setBernRand_LK4_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED && RK3_ENABLED
        module function test_setBernRand_LK4_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED && RK2_ENABLED
        module function test_setBernRand_LK4_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED && RK1_ENABLED
        module function test_setBernRand_LK4_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK3_ENABLED && RK5_ENABLED
        module function test_setBernRand_LK3_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED && RK4_ENABLED
        module function test_setBernRand_LK3_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED && RK3_ENABLED
        module function test_setBernRand_LK3_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED && RK2_ENABLED
        module function test_setBernRand_LK3_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED && RK1_ENABLED
        module function test_setBernRand_LK3_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK2_ENABLED && RK5_ENABLED
        module function test_setBernRand_LK2_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED && RK4_ENABLED
        module function test_setBernRand_LK2_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED && RK3_ENABLED
        module function test_setBernRand_LK2_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED && RK2_ENABLED
        module function test_setBernRand_LK2_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED && RK1_ENABLED
        module function test_setBernRand_LK2_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK1_ENABLED && RK5_ENABLED
        module function test_setBernRand_LK1_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED && RK4_ENABLED
        module function test_setBernRand_LK1_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED && RK3_ENABLED
        module function test_setBernRand_LK1_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED && RK2_ENABLED
        module function test_setBernRand_LK1_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED && RK1_ENABLED
        module function test_setBernRand_LK1_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setBernRand_RK5_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setBernRand_RK4_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setBernRand_RK3_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setBernRand_RK2_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setBernRand_RK1_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_isHead_RK5_1, SK_"test_isHead_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_isHead_RK4_1, SK_"test_isHead_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_isHead_RK3_1, SK_"test_isHead_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_isHead_RK2_1, SK_"test_isHead_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_isHead_RK1_1, SK_"test_isHead_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getBernRand_RK5_1, SK_"test_getBernRand_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_getBernRand_RK4_1, SK_"test_getBernRand_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_getBernRand_RK3_1, SK_"test_getBernRand_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_getBernRand_RK2_1, SK_"test_getBernRand_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_getBernRand_RK1_1, SK_"test_getBernRand_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED && RK5_ENABLED
        call test%run(test_setBernRand_IK5_RK5_1, SK_"test_setBernRand_IK5_RK5_1")
#endif
#if     IK5_ENABLED && RK4_ENABLED
        call test%run(test_setBernRand_IK5_RK4_1, SK_"test_setBernRand_IK5_RK4_1")
#endif
#if     IK5_ENABLED && RK3_ENABLED
        call test%run(test_setBernRand_IK5_RK3_1, SK_"test_setBernRand_IK5_RK3_1")
#endif
#if     IK5_ENABLED && RK2_ENABLED
        call test%run(test_setBernRand_IK5_RK2_1, SK_"test_setBernRand_IK5_RK2_1")
#endif
#if     IK5_ENABLED && RK1_ENABLED
        call test%run(test_setBernRand_IK5_RK1_1, SK_"test_setBernRand_IK5_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK4_ENABLED && RK5_ENABLED
        call test%run(test_setBernRand_IK4_RK5_1, SK_"test_setBernRand_IK4_RK5_1")
#endif
#if     IK4_ENABLED && RK4_ENABLED
        call test%run(test_setBernRand_IK4_RK4_1, SK_"test_setBernRand_IK4_RK4_1")
#endif
#if     IK4_ENABLED && RK3_ENABLED
        call test%run(test_setBernRand_IK4_RK3_1, SK_"test_setBernRand_IK4_RK3_1")
#endif
#if     IK4_ENABLED && RK2_ENABLED
        call test%run(test_setBernRand_IK4_RK2_1, SK_"test_setBernRand_IK4_RK2_1")
#endif
#if     IK4_ENABLED && RK1_ENABLED
        call test%run(test_setBernRand_IK4_RK1_1, SK_"test_setBernRand_IK4_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK3_ENABLED && RK5_ENABLED
        call test%run(test_setBernRand_IK3_RK5_1, SK_"test_setBernRand_IK3_RK5_1")
#endif
#if     IK3_ENABLED && RK4_ENABLED
        call test%run(test_setBernRand_IK3_RK4_1, SK_"test_setBernRand_IK3_RK4_1")
#endif
#if     IK3_ENABLED && RK3_ENABLED
        call test%run(test_setBernRand_IK3_RK3_1, SK_"test_setBernRand_IK3_RK3_1")
#endif
#if     IK3_ENABLED && RK2_ENABLED
        call test%run(test_setBernRand_IK3_RK2_1, SK_"test_setBernRand_IK3_RK2_1")
#endif
#if     IK3_ENABLED && RK1_ENABLED
        call test%run(test_setBernRand_IK3_RK1_1, SK_"test_setBernRand_IK3_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK2_ENABLED && RK5_ENABLED
        call test%run(test_setBernRand_IK2_RK5_1, SK_"test_setBernRand_IK2_RK5_1")
#endif
#if     IK2_ENABLED && RK4_ENABLED
        call test%run(test_setBernRand_IK2_RK4_1, SK_"test_setBernRand_IK2_RK4_1")
#endif
#if     IK2_ENABLED && RK3_ENABLED
        call test%run(test_setBernRand_IK2_RK3_1, SK_"test_setBernRand_IK2_RK3_1")
#endif
#if     IK2_ENABLED && RK2_ENABLED
        call test%run(test_setBernRand_IK2_RK2_1, SK_"test_setBernRand_IK2_RK2_1")
#endif
#if     IK2_ENABLED && RK1_ENABLED
        call test%run(test_setBernRand_IK2_RK1_1, SK_"test_setBernRand_IK2_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK1_ENABLED && RK5_ENABLED
        call test%run(test_setBernRand_IK1_RK5_1, SK_"test_setBernRand_IK1_RK5_1")
#endif
#if     IK1_ENABLED && RK4_ENABLED
        call test%run(test_setBernRand_IK1_RK4_1, SK_"test_setBernRand_IK1_RK4_1")
#endif
#if     IK1_ENABLED && RK3_ENABLED
        call test%run(test_setBernRand_IK1_RK3_1, SK_"test_setBernRand_IK1_RK3_1")
#endif
#if     IK1_ENABLED && RK2_ENABLED
        call test%run(test_setBernRand_IK1_RK2_1, SK_"test_setBernRand_IK1_RK2_1")
#endif
#if     IK1_ENABLED && RK1_ENABLED
        call test%run(test_setBernRand_IK1_RK1_1, SK_"test_setBernRand_IK1_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED && RK5_ENABLED
        call test%run(test_setBernRand_LK5_RK5_1, SK_"test_setBernRand_LK5_RK5_1")
#endif
#if     LK5_ENABLED && RK4_ENABLED
        call test%run(test_setBernRand_LK5_RK4_1, SK_"test_setBernRand_LK5_RK4_1")
#endif
#if     LK5_ENABLED && RK3_ENABLED
        call test%run(test_setBernRand_LK5_RK3_1, SK_"test_setBernRand_LK5_RK3_1")
#endif
#if     LK5_ENABLED && RK2_ENABLED
        call test%run(test_setBernRand_LK5_RK2_1, SK_"test_setBernRand_LK5_RK2_1")
#endif
#if     LK5_ENABLED && RK1_ENABLED
        call test%run(test_setBernRand_LK5_RK1_1, SK_"test_setBernRand_LK5_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK4_ENABLED && RK5_ENABLED
        call test%run(test_setBernRand_LK4_RK5_1, SK_"test_setBernRand_LK4_RK5_1")
#endif
#if     LK4_ENABLED && RK4_ENABLED
        call test%run(test_setBernRand_LK4_RK4_1, SK_"test_setBernRand_LK4_RK4_1")
#endif
#if     LK4_ENABLED && RK3_ENABLED
        call test%run(test_setBernRand_LK4_RK3_1, SK_"test_setBernRand_LK4_RK3_1")
#endif
#if     LK4_ENABLED && RK2_ENABLED
        call test%run(test_setBernRand_LK4_RK2_1, SK_"test_setBernRand_LK4_RK2_1")
#endif
#if     LK4_ENABLED && RK1_ENABLED
        call test%run(test_setBernRand_LK4_RK1_1, SK_"test_setBernRand_LK4_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK3_ENABLED && RK5_ENABLED
        call test%run(test_setBernRand_LK3_RK5_1, SK_"test_setBernRand_LK3_RK5_1")
#endif
#if     LK3_ENABLED && RK4_ENABLED
        call test%run(test_setBernRand_LK3_RK4_1, SK_"test_setBernRand_LK3_RK4_1")
#endif
#if     LK3_ENABLED && RK3_ENABLED
        call test%run(test_setBernRand_LK3_RK3_1, SK_"test_setBernRand_LK3_RK3_1")
#endif
#if     LK3_ENABLED && RK2_ENABLED
        call test%run(test_setBernRand_LK3_RK2_1, SK_"test_setBernRand_LK3_RK2_1")
#endif
#if     LK3_ENABLED && RK1_ENABLED
        call test%run(test_setBernRand_LK3_RK1_1, SK_"test_setBernRand_LK3_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK2_ENABLED && RK5_ENABLED
        call test%run(test_setBernRand_LK2_RK5_1, SK_"test_setBernRand_LK2_RK5_1")
#endif
#if     LK2_ENABLED && RK4_ENABLED
        call test%run(test_setBernRand_LK2_RK4_1, SK_"test_setBernRand_LK2_RK4_1")
#endif
#if     LK2_ENABLED && RK3_ENABLED
        call test%run(test_setBernRand_LK2_RK3_1, SK_"test_setBernRand_LK2_RK3_1")
#endif
#if     LK2_ENABLED && RK2_ENABLED
        call test%run(test_setBernRand_LK2_RK2_1, SK_"test_setBernRand_LK2_RK2_1")
#endif
#if     LK2_ENABLED && RK1_ENABLED
        call test%run(test_setBernRand_LK2_RK1_1, SK_"test_setBernRand_LK2_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK1_ENABLED && RK5_ENABLED
        call test%run(test_setBernRand_LK1_RK5_1, SK_"test_setBernRand_LK1_RK5_1")
#endif
#if     LK1_ENABLED && RK4_ENABLED
        call test%run(test_setBernRand_LK1_RK4_1, SK_"test_setBernRand_LK1_RK4_1")
#endif
#if     LK1_ENABLED && RK3_ENABLED
        call test%run(test_setBernRand_LK1_RK3_1, SK_"test_setBernRand_LK1_RK3_1")
#endif
#if     LK1_ENABLED && RK2_ENABLED
        call test%run(test_setBernRand_LK1_RK2_1, SK_"test_setBernRand_LK1_RK2_1")
#endif
#if     LK1_ENABLED && RK1_ENABLED
        call test%run(test_setBernRand_LK1_RK1_1, SK_"test_setBernRand_LK1_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setBernRand_RK5_RK5_1, SK_"test_setBernRand_RK5_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_setBernRand_RK4_RK4_1, SK_"test_setBernRand_RK4_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_setBernRand_RK3_RK3_1, SK_"test_setBernRand_RK3_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_setBernRand_RK2_RK2_1, SK_"test_setBernRand_RK2_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_setBernRand_RK1_RK1_1, SK_"test_setBernRand_RK1_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_distBern