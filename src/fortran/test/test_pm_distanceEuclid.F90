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

!>  \brief This module contains tests of the module [pm_distanceEuclid](@ref pm_distanceEuclid).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, March 22, 2012, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_distanceEuclid

    use pm_distanceEuclid
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    use pm_kind, only: LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getDisEuclid_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getDisEuclid_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getDisEuclid_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getDisEuclid_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getDisEuclid_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setDisEuclid_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setDisEuclid_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setDisEuclid_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setDisEuclid_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setDisEuclid_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_getDisMatEuclid_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getDisMatEuclid_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getDisMatEuclid_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getDisMatEuclid_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getDisMatEuclid_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setDisMatEuclid_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setDisMatEuclid_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setDisMatEuclid_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setDisMatEuclid_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setDisMatEuclid_RK1() result(assertion); logical(LK) :: assertion; end function
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
        call test%run(test_getDisEuclid_RK5, SK_"test_getDisEuclid_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getDisEuclid_RK4, SK_"test_getDisEuclid_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getDisEuclid_RK3, SK_"test_getDisEuclid_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getDisEuclid_RK2, SK_"test_getDisEuclid_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getDisEuclid_RK1, SK_"test_getDisEuclid_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setDisEuclid_RK5, SK_"test_setDisEuclid_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setDisEuclid_RK4, SK_"test_setDisEuclid_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setDisEuclid_RK3, SK_"test_setDisEuclid_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setDisEuclid_RK2, SK_"test_setDisEuclid_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setDisEuclid_RK1, SK_"test_setDisEuclid_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getDisMatEuclid_RK5, SK_"test_getDisMatEuclid_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getDisMatEuclid_RK4, SK_"test_getDisMatEuclid_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getDisMatEuclid_RK3, SK_"test_getDisMatEuclid_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getDisMatEuclid_RK2, SK_"test_getDisMatEuclid_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getDisMatEuclid_RK1, SK_"test_getDisMatEuclid_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setDisMatEuclid_RK5, SK_"test_setDisMatEuclid_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setDisMatEuclid_RK4, SK_"test_setDisMatEuclid_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setDisMatEuclid_RK3, SK_"test_setDisMatEuclid_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setDisMatEuclid_RK2, SK_"test_setDisMatEuclid_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setDisMatEuclid_RK1, SK_"test_setDisMatEuclid_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_distanceEuclid