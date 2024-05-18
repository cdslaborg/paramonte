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
!>  This module contains tests of the module [pm_distUnif](@ref pm_distUnif).
!>
!>  \final
!>
!>  \author 
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_distUnif

    use pm_distUnif
    use pm_err, only: err_type
    use pm_test, only: test_type, LK

    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     IK5_ENABLED
        module function test_getUnifCDF_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_getUnifCDF_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_getUnifCDF_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_getUnifCDF_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_getUnifCDF_IK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     CK5_ENABLED
        module function test_getUnifCDF_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_getUnifCDF_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_getUnifCDF_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_getUnifCDF_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_getUnifCDF_CK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK5_ENABLED
        module function test_getUnifCDF_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getUnifCDF_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getUnifCDF_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getUnifCDF_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getUnifCDF_RK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK5_ENABLED && IK5_ENABLED
        module function test_setUnifCDF_RK5_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK5_ENABLED && IK4_ENABLED
        module function test_setUnifCDF_RK5_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK5_ENABLED && IK3_ENABLED
        module function test_setUnifCDF_RK5_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK5_ENABLED && IK2_ENABLED
        module function test_setUnifCDF_RK5_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK5_ENABLED && IK1_ENABLED
        module function test_setUnifCDF_RK5_IK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK4_ENABLED && IK5_ENABLED
        module function test_setUnifCDF_RK4_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED && IK4_ENABLED
        module function test_setUnifCDF_RK4_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED && IK3_ENABLED
        module function test_setUnifCDF_RK4_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED && IK2_ENABLED
        module function test_setUnifCDF_RK4_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED && IK1_ENABLED
        module function test_setUnifCDF_RK4_IK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK3_ENABLED && IK5_ENABLED
        module function test_setUnifCDF_RK3_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED && IK4_ENABLED
        module function test_setUnifCDF_RK3_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED && IK3_ENABLED
        module function test_setUnifCDF_RK3_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED && IK2_ENABLED
        module function test_setUnifCDF_RK3_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED && IK1_ENABLED
        module function test_setUnifCDF_RK3_IK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK2_ENABLED && IK5_ENABLED
        module function test_setUnifCDF_RK2_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED && IK4_ENABLED
        module function test_setUnifCDF_RK2_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED && IK3_ENABLED
        module function test_setUnifCDF_RK2_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED && IK2_ENABLED
        module function test_setUnifCDF_RK2_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED && IK1_ENABLED
        module function test_setUnifCDF_RK2_IK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK1_ENABLED && IK5_ENABLED
        module function test_setUnifCDF_RK1_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED && IK4_ENABLED
        module function test_setUnifCDF_RK1_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED && IK3_ENABLED
        module function test_setUnifCDF_RK1_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED && IK2_ENABLED
        module function test_setUnifCDF_RK1_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED && IK1_ENABLED
        module function test_setUnifCDF_RK1_IK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     CK5_ENABLED
        module function test_setUnifCDF_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setUnifCDF_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setUnifCDF_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setUnifCDF_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setUnifCDF_CK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK5_ENABLED
        module function test_setUnifCDF_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setUnifCDF_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setUnifCDF_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setUnifCDF_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setUnifCDF_RK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     SK5_ENABLED
        module function test_getUnifRand_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_getUnifRand_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_getUnifRand_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_getUnifRand_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_getUnifRand_SK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     IK5_ENABLED
        module function test_getUnifRand_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_getUnifRand_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_getUnifRand_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_getUnifRand_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_getUnifRand_IK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     LK5_ENABLED
        module function test_getUnifRand_LK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_getUnifRand_LK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_getUnifRand_LK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_getUnifRand_LK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_getUnifRand_LK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     CK5_ENABLED
        module function test_getUnifRand_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_getUnifRand_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_getUnifRand_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_getUnifRand_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_getUnifRand_CK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK5_ENABLED
        module function test_getUnifRand_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_getUnifRand_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_getUnifRand_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_getUnifRand_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_getUnifRand_RK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     SK5_ENABLED
        module function test_setUnifRand_SK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK4_ENABLED
        module function test_setUnifRand_SK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK3_ENABLED
        module function test_setUnifRand_SK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK2_ENABLED
        module function test_setUnifRand_SK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     SK1_ENABLED
        module function test_setUnifRand_SK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     IK5_ENABLED
        module function test_setUnifRand_IK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED
        module function test_setUnifRand_IK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED
        module function test_setUnifRand_IK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED
        module function test_setUnifRand_IK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED
        module function test_setUnifRand_IK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     LK5_ENABLED
        module function test_setUnifRand_LK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK4_ENABLED
        module function test_setUnifRand_LK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK3_ENABLED
        module function test_setUnifRand_LK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK2_ENABLED
        module function test_setUnifRand_LK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     LK1_ENABLED
        module function test_setUnifRand_LK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     CK5_ENABLED
        module function test_setUnifRand_CK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK4_ENABLED
        module function test_setUnifRand_CK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK3_ENABLED
        module function test_setUnifRand_CK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK2_ENABLED
        module function test_setUnifRand_CK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     CK1_ENABLED
        module function test_setUnifRand_CK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
#if     RK5_ENABLED
        module function test_setUnifRand_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setUnifRand_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setUnifRand_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setUnifRand_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setUnifRand_RK1() result(assertion); logical(LK) :: assertion; end function
#endif
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_getUnifCDF_IK5, SK_"test_getUnifCDF_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_getUnifCDF_IK4, SK_"test_getUnifCDF_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_getUnifCDF_IK3, SK_"test_getUnifCDF_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_getUnifCDF_IK2, SK_"test_getUnifCDF_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_getUnifCDF_IK1, SK_"test_getUnifCDF_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_getUnifCDF_CK5, SK_"test_getUnifCDF_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_getUnifCDF_CK4, SK_"test_getUnifCDF_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_getUnifCDF_CK3, SK_"test_getUnifCDF_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_getUnifCDF_CK2, SK_"test_getUnifCDF_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_getUnifCDF_CK1, SK_"test_getUnifCDF_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getUnifCDF_RK5, SK_"test_getUnifCDF_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getUnifCDF_RK4, SK_"test_getUnifCDF_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getUnifCDF_RK3, SK_"test_getUnifCDF_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getUnifCDF_RK2, SK_"test_getUnifCDF_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getUnifCDF_RK1, SK_"test_getUnifCDF_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED && IK5_ENABLED
        call test%run(test_setUnifCDF_RK5_IK5, SK_"test_setUnifCDF_RK5_IK5")
#endif
#if     RK5_ENABLED && IK4_ENABLED
        call test%run(test_setUnifCDF_RK5_IK4, SK_"test_setUnifCDF_RK5_IK4")
#endif
#if     RK5_ENABLED && IK3_ENABLED
        call test%run(test_setUnifCDF_RK5_IK3, SK_"test_setUnifCDF_RK5_IK3")
#endif
#if     RK5_ENABLED && IK2_ENABLED
        call test%run(test_setUnifCDF_RK5_IK2, SK_"test_setUnifCDF_RK5_IK2")
#endif
#if     RK5_ENABLED && IK1_ENABLED
        call test%run(test_setUnifCDF_RK5_IK1, SK_"test_setUnifCDF_RK5_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK4_ENABLED && IK5_ENABLED
        call test%run(test_setUnifCDF_RK4_IK5, SK_"test_setUnifCDF_RK4_IK5")
#endif
#if     RK4_ENABLED && IK4_ENABLED
        call test%run(test_setUnifCDF_RK4_IK4, SK_"test_setUnifCDF_RK4_IK4")
#endif
#if     RK4_ENABLED && IK3_ENABLED
        call test%run(test_setUnifCDF_RK4_IK3, SK_"test_setUnifCDF_RK4_IK3")
#endif
#if     RK4_ENABLED && IK2_ENABLED
        call test%run(test_setUnifCDF_RK4_IK2, SK_"test_setUnifCDF_RK4_IK2")
#endif
#if     RK4_ENABLED && IK1_ENABLED
        call test%run(test_setUnifCDF_RK4_IK1, SK_"test_setUnifCDF_RK4_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK3_ENABLED && IK5_ENABLED
        call test%run(test_setUnifCDF_RK3_IK5, SK_"test_setUnifCDF_RK3_IK5")
#endif
#if     RK3_ENABLED && IK4_ENABLED
        call test%run(test_setUnifCDF_RK3_IK4, SK_"test_setUnifCDF_RK3_IK4")
#endif
#if     RK3_ENABLED && IK3_ENABLED
        call test%run(test_setUnifCDF_RK3_IK3, SK_"test_setUnifCDF_RK3_IK3")
#endif
#if     RK3_ENABLED && IK2_ENABLED
        call test%run(test_setUnifCDF_RK3_IK2, SK_"test_setUnifCDF_RK3_IK2")
#endif
#if     RK3_ENABLED && IK1_ENABLED
        call test%run(test_setUnifCDF_RK3_IK1, SK_"test_setUnifCDF_RK3_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK2_ENABLED && IK5_ENABLED
        call test%run(test_setUnifCDF_RK2_IK5, SK_"test_setUnifCDF_RK2_IK5")
#endif
#if     RK2_ENABLED && IK4_ENABLED
        call test%run(test_setUnifCDF_RK2_IK4, SK_"test_setUnifCDF_RK2_IK4")
#endif
#if     RK2_ENABLED && IK3_ENABLED
        call test%run(test_setUnifCDF_RK2_IK3, SK_"test_setUnifCDF_RK2_IK3")
#endif
#if     RK2_ENABLED && IK2_ENABLED
        call test%run(test_setUnifCDF_RK2_IK2, SK_"test_setUnifCDF_RK2_IK2")
#endif
#if     RK2_ENABLED && IK1_ENABLED
        call test%run(test_setUnifCDF_RK2_IK1, SK_"test_setUnifCDF_RK2_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK1_ENABLED && IK5_ENABLED
        call test%run(test_setUnifCDF_RK1_IK5, SK_"test_setUnifCDF_RK1_IK5")
#endif
#if     RK1_ENABLED && IK4_ENABLED
        call test%run(test_setUnifCDF_RK1_IK4, SK_"test_setUnifCDF_RK1_IK4")
#endif
#if     RK1_ENABLED && IK3_ENABLED
        call test%run(test_setUnifCDF_RK1_IK3, SK_"test_setUnifCDF_RK1_IK3")
#endif
#if     RK1_ENABLED && IK2_ENABLED
        call test%run(test_setUnifCDF_RK1_IK2, SK_"test_setUnifCDF_RK1_IK2")
#endif
#if     RK1_ENABLED && IK1_ENABLED
        call test%run(test_setUnifCDF_RK1_IK1, SK_"test_setUnifCDF_RK1_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setUnifCDF_CK5, SK_"test_setUnifCDF_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_setUnifCDF_CK4, SK_"test_setUnifCDF_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_setUnifCDF_CK3, SK_"test_setUnifCDF_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_setUnifCDF_CK2, SK_"test_setUnifCDF_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_setUnifCDF_CK1, SK_"test_setUnifCDF_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setUnifCDF_RK5, SK_"test_setUnifCDF_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setUnifCDF_RK4, SK_"test_setUnifCDF_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setUnifCDF_RK3, SK_"test_setUnifCDF_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setUnifCDF_RK2, SK_"test_setUnifCDF_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setUnifCDF_RK1, SK_"test_setUnifCDF_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_getUnifRand_SK5, SK_"test_getUnifRand_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_getUnifRand_SK4, SK_"test_getUnifRand_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_getUnifRand_SK3, SK_"test_getUnifRand_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_getUnifRand_SK2, SK_"test_getUnifRand_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_getUnifRand_SK1, SK_"test_getUnifRand_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_getUnifRand_IK5, SK_"test_getUnifRand_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_getUnifRand_IK4, SK_"test_getUnifRand_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_getUnifRand_IK3, SK_"test_getUnifRand_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_getUnifRand_IK2, SK_"test_getUnifRand_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_getUnifRand_IK1, SK_"test_getUnifRand_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_getUnifRand_LK5, SK_"test_getUnifRand_LK5")
#endif
#if     LK4_ENABLED
        call test%run(test_getUnifRand_LK4, SK_"test_getUnifRand_LK4")
#endif
#if     LK3_ENABLED
        call test%run(test_getUnifRand_LK3, SK_"test_getUnifRand_LK3")
#endif
#if     LK2_ENABLED
        call test%run(test_getUnifRand_LK2, SK_"test_getUnifRand_LK2")
#endif
#if     LK1_ENABLED
        call test%run(test_getUnifRand_LK1, SK_"test_getUnifRand_LK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_getUnifRand_CK5, SK_"test_getUnifRand_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_getUnifRand_CK4, SK_"test_getUnifRand_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_getUnifRand_CK3, SK_"test_getUnifRand_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_getUnifRand_CK2, SK_"test_getUnifRand_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_getUnifRand_CK1, SK_"test_getUnifRand_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_getUnifRand_RK5, SK_"test_getUnifRand_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_getUnifRand_RK4, SK_"test_getUnifRand_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_getUnifRand_RK3, SK_"test_getUnifRand_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_getUnifRand_RK2, SK_"test_getUnifRand_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_getUnifRand_RK1, SK_"test_getUnifRand_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK5_ENABLED
        call test%run(test_setUnifRand_SK5, SK_"test_setUnifRand_SK5")
#endif
#if     SK4_ENABLED
        call test%run(test_setUnifRand_SK4, SK_"test_setUnifRand_SK4")
#endif
#if     SK3_ENABLED
        call test%run(test_setUnifRand_SK3, SK_"test_setUnifRand_SK3")
#endif
#if     SK2_ENABLED
        call test%run(test_setUnifRand_SK2, SK_"test_setUnifRand_SK2")
#endif
#if     SK1_ENABLED
        call test%run(test_setUnifRand_SK1, SK_"test_setUnifRand_SK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED
        call test%run(test_setUnifRand_IK5, SK_"test_setUnifRand_IK5")
#endif
#if     IK4_ENABLED
        call test%run(test_setUnifRand_IK4, SK_"test_setUnifRand_IK4")
#endif
#if     IK3_ENABLED
        call test%run(test_setUnifRand_IK3, SK_"test_setUnifRand_IK3")
#endif
#if     IK2_ENABLED
        call test%run(test_setUnifRand_IK2, SK_"test_setUnifRand_IK2")
#endif
#if     IK1_ENABLED
        call test%run(test_setUnifRand_IK1, SK_"test_setUnifRand_IK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK5_ENABLED
        call test%run(test_setUnifRand_LK5, SK_"test_setUnifRand_LK5")
#endif
#if     LK4_ENABLED
        call test%run(test_setUnifRand_LK4, SK_"test_setUnifRand_LK4")
#endif
#if     LK3_ENABLED
        call test%run(test_setUnifRand_LK3, SK_"test_setUnifRand_LK3")
#endif
#if     LK2_ENABLED
        call test%run(test_setUnifRand_LK2, SK_"test_setUnifRand_LK2")
#endif
#if     LK1_ENABLED
        call test%run(test_setUnifRand_LK1, SK_"test_setUnifRand_LK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK5_ENABLED
        call test%run(test_setUnifRand_CK5, SK_"test_setUnifRand_CK5")
#endif
#if     CK4_ENABLED
        call test%run(test_setUnifRand_CK4, SK_"test_setUnifRand_CK4")
#endif
#if     CK3_ENABLED
        call test%run(test_setUnifRand_CK3, SK_"test_setUnifRand_CK3")
#endif
#if     CK2_ENABLED
        call test%run(test_setUnifRand_CK2, SK_"test_setUnifRand_CK2")
#endif
#if     CK1_ENABLED
        call test%run(test_setUnifRand_CK1, SK_"test_setUnifRand_CK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        call test%run(test_setUnifRand_RK5, SK_"test_setUnifRand_RK5")
#endif
#if     RK4_ENABLED
        call test%run(test_setUnifRand_RK4, SK_"test_setUnifRand_RK4")
#endif
#if     RK3_ENABLED
        call test%run(test_setUnifRand_RK3, SK_"test_setUnifRand_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setUnifRand_RK2, SK_"test_setUnifRand_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setUnifRand_RK1, SK_"test_setUnifRand_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_distUnif