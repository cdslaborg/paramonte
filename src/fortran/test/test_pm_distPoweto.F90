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

!>  \brief This module contains tests of the module [pm_distPoweto](@ref pm_distPoweto).
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_distPoweto

    use pm_distPoweto
    use pm_err, only: err_type
    use pm_test, only: test_type, IK, LK

    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if RK5_ENABLED
    module function test_getPowetoLogPDF_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK4_ENABLED
    module function test_getPowetoLogPDF_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK3_ENABLED
    module function test_getPowetoLogPDF_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function test_getPowetoLogPDF_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function test_getPowetoLogPDF_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if RK5_ENABLED
    module function test_setPowetoLogPDF_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK4_ENABLED
    module function test_setPowetoLogPDF_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK3_ENABLED
    module function test_setPowetoLogPDF_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function test_setPowetoLogPDF_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function test_setPowetoLogPDF_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_getPowetoCDF_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getPowetoCDF_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getPowetoCDF_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getPowetoCDF_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getPowetoCDF_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_setPowetoCDF_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setPowetoCDF_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setPowetoCDF_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setPowetoCDF_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setPowetoCDF_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
        call test%run(test_getPowetoLogPDF_RK5_1, SK_"test_getPowetoLogPDF_RK5_1")
#endif
#if RK4_ENABLED
        call test%run(test_getPowetoLogPDF_RK4_1, SK_"test_getPowetoLogPDF_RK4_1")
#endif
#if RK3_ENABLED
        call test%run(test_getPowetoLogPDF_RK3_1, SK_"test_getPowetoLogPDF_RK3_1")
#endif
#if RK2_ENABLED
        call test%run(test_getPowetoLogPDF_RK2_1, SK_"test_getPowetoLogPDF_RK2_1")
#endif
#if RK1_ENABLED
        call test%run(test_getPowetoLogPDF_RK1_1, SK_"test_getPowetoLogPDF_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
        call test%run(test_setPowetoLogPDF_RK5_1, SK_"test_setPowetoLogPDF_RK5_1")
#endif
#if RK4_ENABLED
        call test%run(test_setPowetoLogPDF_RK4_1, SK_"test_setPowetoLogPDF_RK4_1")
#endif
#if RK3_ENABLED
        call test%run(test_setPowetoLogPDF_RK3_1, SK_"test_setPowetoLogPDF_RK3_1")
#endif
#if RK2_ENABLED
        call test%run(test_setPowetoLogPDF_RK2_1, SK_"test_setPowetoLogPDF_RK2_1")
#endif
#if RK1_ENABLED
        call test%run(test_setPowetoLogPDF_RK1_1, SK_"test_setPowetoLogPDF_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getPowetoCDF_RK5_1, SK_"test_getPowetoCDF_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_getPowetoCDF_RK4_1, SK_"test_getPowetoCDF_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_getPowetoCDF_RK3_1, SK_"test_getPowetoCDF_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_getPowetoCDF_RK2_1, SK_"test_getPowetoCDF_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_getPowetoCDF_RK1_1, SK_"test_getPowetoCDF_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setPowetoCDF_RK5_1, SK_"test_setPowetoCDF_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_setPowetoCDF_RK4_1, SK_"test_setPowetoCDF_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_setPowetoCDF_RK3_1, SK_"test_setPowetoCDF_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_setPowetoCDF_RK2_1, SK_"test_setPowetoCDF_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_setPowetoCDF_RK1_1, SK_"test_setPowetoCDF_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_distPoweto