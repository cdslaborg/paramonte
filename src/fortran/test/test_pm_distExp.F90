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

!>  \brief This module contains tests of the module [pm_distExp](@ref pm_distExp).
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_distExp

    use pm_distExp
    use pm_err, only: err_type
    use pm_test, only: test_type, LK

    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_getExpLogPDF_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getExpLogPDF_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getExpLogPDF_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getExpLogPDF_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getExpLogPDF_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_setExpLogPDF_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setExpLogPDF_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setExpLogPDF_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setExpLogPDF_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setExpLogPDF_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_getExpCDF_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getExpCDF_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getExpCDF_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getExpCDF_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getExpCDF_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_setExpCDF_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setExpCDF_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setExpCDF_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setExpCDF_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setExpCDF_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK5_ENABLED
        module function test_getExpRand_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getExpRand_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getExpRand_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getExpRand_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getExpRand_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK5_ENABLED
        module function test_setExpRand_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setExpRand_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setExpRand_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setExpRand_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setExpRand_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getExpLogPDF_RK5, SK_"test_getExpLogPDF_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getExpLogPDF_RK4, SK_"test_getExpLogPDF_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getExpLogPDF_RK3, SK_"test_getExpLogPDF_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getExpLogPDF_RK2, SK_"test_getExpLogPDF_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getExpLogPDF_RK1, SK_"test_getExpLogPDF_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setExpLogPDF_RK5, SK_"test_setExpLogPDF_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setExpLogPDF_RK4, SK_"test_setExpLogPDF_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setExpLogPDF_RK3, SK_"test_setExpLogPDF_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setExpLogPDF_RK2, SK_"test_setExpLogPDF_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setExpLogPDF_RK1, SK_"test_setExpLogPDF_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getExpCDF_RK5, SK_"test_getExpCDF_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getExpCDF_RK4, SK_"test_getExpCDF_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getExpCDF_RK3, SK_"test_getExpCDF_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getExpCDF_RK2, SK_"test_getExpCDF_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getExpCDF_RK1, SK_"test_getExpCDF_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setExpCDF_RK5, SK_"test_setExpCDF_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setExpCDF_RK4, SK_"test_setExpCDF_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setExpCDF_RK3, SK_"test_setExpCDF_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setExpCDF_RK2, SK_"test_setExpCDF_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setExpCDF_RK1, SK_"test_setExpCDF_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
        call test%run(test_getExpRand_RK5_1, SK_"test_getExpRand_RK5_1")
#endif
#if RK4_ENABLED
        call test%run(test_getExpRand_RK4_1, SK_"test_getExpRand_RK4_1")
#endif
#if RK3_ENABLED
        call test%run(test_getExpRand_RK3_1, SK_"test_getExpRand_RK3_1")
#endif
#if RK2_ENABLED
        call test%run(test_getExpRand_RK2_1, SK_"test_getExpRand_RK2_1")
#endif
#if RK1_ENABLED
        call test%run(test_getExpRand_RK1_1, SK_"test_getExpRand_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
        call test%run(test_setExpRand_RK5_1, SK_"test_setExpRand_RK5_1")
#endif
#if RK4_ENABLED
        call test%run(test_setExpRand_RK4_1, SK_"test_setExpRand_RK4_1")
#endif
#if RK3_ENABLED
        call test%run(test_setExpRand_RK3_1, SK_"test_setExpRand_RK3_1")
#endif
#if RK2_ENABLED
        call test%run(test_setExpRand_RK2_1, SK_"test_setExpRand_RK2_1")
#endif
#if RK1_ENABLED
        call test%run(test_setExpRand_RK1_1, SK_"test_setExpRand_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_distExp