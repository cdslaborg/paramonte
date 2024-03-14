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

!>  \brief This module contains tests of the module [pm_distPower](@ref pm_distPower).
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_distPower

    use pm_distPower
    use pm_err, only: err_type
    use pm_test, only: test_type, IK, LK

    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_getPowerLogPDF_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getPowerLogPDF_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getPowerLogPDF_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getPowerLogPDF_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getPowerLogPDF_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_setPowerLogPDF_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setPowerLogPDF_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setPowerLogPDF_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setPowerLogPDF_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setPowerLogPDF_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_getPowerLogCDF_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getPowerLogCDF_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getPowerLogCDF_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getPowerLogCDF_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getPowerLogCDF_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_setPowerLogCDF_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setPowerLogCDF_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setPowerLogCDF_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setPowerLogCDF_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setPowerLogCDF_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_getPowerLogQuan_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getPowerLogQuan_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getPowerLogQuan_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getPowerLogQuan_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getPowerLogQuan_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_setPowerLogQuan_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setPowerLogQuan_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setPowerLogQuan_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setPowerLogQuan_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setPowerLogQuan_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_getPowerLogRand_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getPowerLogRand_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getPowerLogRand_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getPowerLogRand_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getPowerLogRand_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

#if     RK5_ENABLED
        module function test_setPowerLogRand_RK5_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setPowerLogRand_RK4_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setPowerLogRand_RK3_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setPowerLogRand_RK2_1() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setPowerLogRand_RK1_1() result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setPowerLogPDF_RK5_1, SK_"test_setPowerLogPDF_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_setPowerLogPDF_RK4_1, SK_"test_setPowerLogPDF_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_setPowerLogPDF_RK3_1, SK_"test_setPowerLogPDF_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_setPowerLogPDF_RK2_1, SK_"test_setPowerLogPDF_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_setPowerLogPDF_RK1_1, SK_"test_setPowerLogPDF_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getPowerLogPDF_RK5_1, SK_"test_getPowerLogPDF_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_getPowerLogPDF_RK4_1, SK_"test_getPowerLogPDF_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_getPowerLogPDF_RK3_1, SK_"test_getPowerLogPDF_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_getPowerLogPDF_RK2_1, SK_"test_getPowerLogPDF_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_getPowerLogPDF_RK1_1, SK_"test_getPowerLogPDF_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setPowerLogCDF_RK5_1, SK_"test_setPowerLogCDF_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_setPowerLogCDF_RK4_1, SK_"test_setPowerLogCDF_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_setPowerLogCDF_RK3_1, SK_"test_setPowerLogCDF_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_setPowerLogCDF_RK2_1, SK_"test_setPowerLogCDF_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_setPowerLogCDF_RK1_1, SK_"test_setPowerLogCDF_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getPowerLogCDF_RK5_1, SK_"test_getPowerLogCDF_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_getPowerLogCDF_RK4_1, SK_"test_getPowerLogCDF_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_getPowerLogCDF_RK3_1, SK_"test_getPowerLogCDF_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_getPowerLogCDF_RK2_1, SK_"test_getPowerLogCDF_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_getPowerLogCDF_RK1_1, SK_"test_getPowerLogCDF_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setPowerLogQuan_RK5_1, SK_"test_setPowerLogQuan_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_setPowerLogQuan_RK4_1, SK_"test_setPowerLogQuan_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_setPowerLogQuan_RK3_1, SK_"test_setPowerLogQuan_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_setPowerLogQuan_RK2_1, SK_"test_setPowerLogQuan_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_setPowerLogQuan_RK1_1, SK_"test_setPowerLogQuan_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getPowerLogQuan_RK5_1, SK_"test_getPowerLogQuan_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_getPowerLogQuan_RK4_1, SK_"test_getPowerLogQuan_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_getPowerLogQuan_RK3_1, SK_"test_getPowerLogQuan_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_getPowerLogQuan_RK2_1, SK_"test_getPowerLogQuan_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_getPowerLogQuan_RK1_1, SK_"test_getPowerLogQuan_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setPowerLogRand_RK5_1, SK_"test_setPowerLogRand_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_setPowerLogRand_RK4_1, SK_"test_setPowerLogRand_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_setPowerLogRand_RK3_1, SK_"test_setPowerLogRand_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_setPowerLogRand_RK2_1, SK_"test_setPowerLogRand_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_setPowerLogRand_RK1_1, SK_"test_setPowerLogRand_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getPowerLogRand_RK5_1, SK_"test_getPowerLogRand_RK5_1")
#endif
#if     RK4_ENABLED
        call test%run(test_getPowerLogRand_RK4_1, SK_"test_getPowerLogRand_RK4_1")
#endif
#if     RK3_ENABLED
        call test%run(test_getPowerLogRand_RK3_1, SK_"test_getPowerLogRand_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_getPowerLogRand_RK2_1, SK_"test_getPowerLogRand_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_getPowerLogRand_RK1_1, SK_"test_getPowerLogRand_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_distPower