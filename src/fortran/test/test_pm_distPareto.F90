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

!>  \brief This module contains tests of the module [pm_distPareto](@ref pm_distPareto).<br>
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_distPareto

    use pm_distPareto
    use pm_err, only: err_type
    use pm_test, only: test_type, IK, LK

    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_getParetoLogPDF_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getParetoLogPDF_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getParetoLogPDF_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getParetoLogPDF_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getParetoLogPDF_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_setParetoLogPDF_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setParetoLogPDF_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setParetoLogPDF_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setParetoLogPDF_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setParetoLogPDF_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_getParetoLogCDF_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getParetoLogCDF_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getParetoLogCDF_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getParetoLogCDF_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getParetoLogCDF_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_setParetoLogCDF_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setParetoLogCDF_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setParetoLogCDF_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setParetoLogCDF_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setParetoLogCDF_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_getParetoLogQuan_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getParetoLogQuan_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getParetoLogQuan_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getParetoLogQuan_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getParetoLogQuan_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_setParetoLogQuan_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setParetoLogQuan_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setParetoLogQuan_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setParetoLogQuan_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setParetoLogQuan_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_getParetoLogRand_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getParetoLogRand_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getParetoLogRand_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getParetoLogRand_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getParetoLogRand_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_setParetoLogRand_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setParetoLogRand_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setParetoLogRand_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setParetoLogRand_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setParetoLogRand_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setParetoLogPDF_RK5_1, SK_"test_setParetoLogPDF_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_setParetoLogPDF_RK4_1, SK_"test_setParetoLogPDF_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_setParetoLogPDF_RK3_1, SK_"test_setParetoLogPDF_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_setParetoLogPDF_RK2_1, SK_"test_setParetoLogPDF_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_setParetoLogPDF_RK1_1, SK_"test_setParetoLogPDF_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getParetoLogPDF_RK5_1, SK_"test_getParetoLogPDF_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_getParetoLogPDF_RK4_1, SK_"test_getParetoLogPDF_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_getParetoLogPDF_RK3_1, SK_"test_getParetoLogPDF_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_getParetoLogPDF_RK2_1, SK_"test_getParetoLogPDF_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_getParetoLogPDF_RK1_1, SK_"test_getParetoLogPDF_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setParetoLogCDF_RK5_1, SK_"test_setParetoLogCDF_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_setParetoLogCDF_RK4_1, SK_"test_setParetoLogCDF_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_setParetoLogCDF_RK3_1, SK_"test_setParetoLogCDF_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_setParetoLogCDF_RK2_1, SK_"test_setParetoLogCDF_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_setParetoLogCDF_RK1_1, SK_"test_setParetoLogCDF_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getParetoLogCDF_RK5_1, SK_"test_getParetoLogCDF_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_getParetoLogCDF_RK4_1, SK_"test_getParetoLogCDF_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_getParetoLogCDF_RK3_1, SK_"test_getParetoLogCDF_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_getParetoLogCDF_RK2_1, SK_"test_getParetoLogCDF_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_getParetoLogCDF_RK1_1, SK_"test_getParetoLogCDF_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setParetoLogQuan_RK5_1, SK_"test_setParetoLogQuan_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_setParetoLogQuan_RK4_1, SK_"test_setParetoLogQuan_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_setParetoLogQuan_RK3_1, SK_"test_setParetoLogQuan_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_setParetoLogQuan_RK2_1, SK_"test_setParetoLogQuan_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_setParetoLogQuan_RK1_1, SK_"test_setParetoLogQuan_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getParetoLogQuan_RK5_1, SK_"test_getParetoLogQuan_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_getParetoLogQuan_RK4_1, SK_"test_getParetoLogQuan_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_getParetoLogQuan_RK3_1, SK_"test_getParetoLogQuan_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_getParetoLogQuan_RK2_1, SK_"test_getParetoLogQuan_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_getParetoLogQuan_RK1_1, SK_"test_getParetoLogQuan_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setParetoLogRand_RK5_1, SK_"test_setParetoLogRand_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_setParetoLogRand_RK4_1, SK_"test_setParetoLogRand_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_setParetoLogRand_RK3_1, SK_"test_setParetoLogRand_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_setParetoLogRand_RK2_1, SK_"test_setParetoLogRand_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_setParetoLogRand_RK1_1, SK_"test_setParetoLogRand_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getParetoLogRand_RK5_1, SK_"test_getParetoLogRand_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_getParetoLogRand_RK4_1, SK_"test_getParetoLogRand_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_getParetoLogRand_RK3_1, SK_"test_getParetoLogRand_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_getParetoLogRand_RK2_1, SK_"test_getParetoLogRand_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_getParetoLogRand_RK1_1, SK_"test_getParetoLogRand_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_distPareto