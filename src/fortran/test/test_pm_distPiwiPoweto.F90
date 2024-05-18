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

!>  \brief This module contains tests of the module [pm_distPiwiPoweto](@ref pm_distPiwiPoweto).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_distPiwiPoweto

    use pm_distPiwiPoweto
    use pm_err, only: err_type
    use pm_test, only: test_type, IK, LK

    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if RK5_ENABLED
    module function test_getPiwiPowetoLogPDF_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK4_ENABLED
    module function test_getPiwiPowetoLogPDF_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK3_ENABLED
    module function test_getPiwiPowetoLogPDF_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function test_getPiwiPowetoLogPDF_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function test_getPiwiPowetoLogPDF_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if RK5_ENABLED
    module function test_setPiwiPowetoLogPDF_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK4_ENABLED
    module function test_setPiwiPowetoLogPDF_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK3_ENABLED
    module function test_setPiwiPowetoLogPDF_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function test_setPiwiPowetoLogPDF_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function test_setPiwiPowetoLogPDF_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_getPiwiPowetoCDF_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getPiwiPowetoCDF_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getPiwiPowetoCDF_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getPiwiPowetoCDF_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getPiwiPowetoCDF_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_setPiwiPowetoCDF_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setPiwiPowetoCDF_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setPiwiPowetoCDF_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setPiwiPowetoCDF_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setPiwiPowetoCDF_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
        call test%run(test_getPiwiPowetoLogPDF_RK5_1, SK_"test_getPiwiPowetoLogPDF_RK5_1")
#endif
#if RK4_ENABLED
        call test%run(test_getPiwiPowetoLogPDF_RK4_1, SK_"test_getPiwiPowetoLogPDF_RK4_1")
#endif
#if RK3_ENABLED
        call test%run(test_getPiwiPowetoLogPDF_RK3_1, SK_"test_getPiwiPowetoLogPDF_RK3_1")
#endif
#if RK2_ENABLED
        call test%run(test_getPiwiPowetoLogPDF_RK2_1, SK_"test_getPiwiPowetoLogPDF_RK2_1")
#endif
#if RK1_ENABLED
        call test%run(test_getPiwiPowetoLogPDF_RK1_1, SK_"test_getPiwiPowetoLogPDF_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
        call test%run(test_setPiwiPowetoLogPDF_RK5_1, SK_"test_setPiwiPowetoLogPDF_RK5_1")
#endif
#if RK4_ENABLED
        call test%run(test_setPiwiPowetoLogPDF_RK4_1, SK_"test_setPiwiPowetoLogPDF_RK4_1")
#endif
#if RK3_ENABLED
        call test%run(test_setPiwiPowetoLogPDF_RK3_1, SK_"test_setPiwiPowetoLogPDF_RK3_1")
#endif
#if RK2_ENABLED
        call test%run(test_setPiwiPowetoLogPDF_RK2_1, SK_"test_setPiwiPowetoLogPDF_RK2_1")
#endif
#if RK1_ENABLED
        call test%run(test_setPiwiPowetoLogPDF_RK1_1, SK_"test_setPiwiPowetoLogPDF_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getPiwiPowetoCDF_RK5_1, SK_"test_getPiwiPowetoCDF_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_getPiwiPowetoCDF_RK4_1, SK_"test_getPiwiPowetoCDF_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_getPiwiPowetoCDF_RK3_1, SK_"test_getPiwiPowetoCDF_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_getPiwiPowetoCDF_RK2_1, SK_"test_getPiwiPowetoCDF_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_getPiwiPowetoCDF_RK1_1, SK_"test_getPiwiPowetoCDF_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setPiwiPowetoCDF_RK5_1, SK_"test_setPiwiPowetoCDF_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_setPiwiPowetoCDF_RK4_1, SK_"test_setPiwiPowetoCDF_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_setPiwiPowetoCDF_RK3_1, SK_"test_setPiwiPowetoCDF_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_setPiwiPowetoCDF_RK2_1, SK_"test_setPiwiPowetoCDF_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_setPiwiPowetoCDF_RK1_1, SK_"test_setPiwiPowetoCDF_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_distPiwiPoweto