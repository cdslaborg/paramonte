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
!>  This module contains tests of the module [pm_sampleECDF](@ref pm_sampleECDF).
!>
!>  \author
!>  \FatemehBagheri, Thursday 12:45 AM, August 20, 2021, Dallas, TX

module test_pm_sampleECDF

    use pm_sampleECDF
    use pm_kind, only: IK, RK
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK5_ENABLED
        module function test_setECDF_D1_RK5_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK4_ENABLED
        module function test_setECDF_D1_RK4_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK3_ENABLED
        module function test_setECDF_D1_RK3_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK2_ENABLED
        module function test_setECDF_D1_RK2_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     RK1_ENABLED
        module function test_setECDF_D1_RK1_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK5_ENABLED && RK5_ENABLED
        module function test_setECDF_D1_IK5_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK5_ENABLED && RK4_ENABLED
        module function test_setECDF_D1_IK5_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK5_ENABLED && RK3_ENABLED
        module function test_setECDF_D1_IK5_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK5_ENABLED && RK2_ENABLED
        module function test_setECDF_D1_IK5_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK5_ENABLED && RK1_ENABLED
        module function test_setECDF_D1_IK5_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK4_ENABLED && RK5_ENABLED
        module function test_setECDF_D1_IK4_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED && RK4_ENABLED
        module function test_setECDF_D1_IK4_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED && RK3_ENABLED
        module function test_setECDF_D1_IK4_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED && RK2_ENABLED
        module function test_setECDF_D1_IK4_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK4_ENABLED && RK1_ENABLED
        module function test_setECDF_D1_IK4_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK3_ENABLED && RK5_ENABLED
        module function test_setECDF_D1_IK3_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED && RK4_ENABLED
        module function test_setECDF_D1_IK3_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED && RK3_ENABLED
        module function test_setECDF_D1_IK3_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED && RK2_ENABLED
        module function test_setECDF_D1_IK3_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK3_ENABLED && RK1_ENABLED
        module function test_setECDF_D1_IK3_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK2_ENABLED && RK5_ENABLED
        module function test_setECDF_D1_IK2_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED && RK4_ENABLED
        module function test_setECDF_D1_IK2_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED && RK3_ENABLED
        module function test_setECDF_D1_IK2_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED && RK2_ENABLED
        module function test_setECDF_D1_IK2_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK2_ENABLED && RK1_ENABLED
        module function test_setECDF_D1_IK2_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK1_ENABLED && RK5_ENABLED
        module function test_setECDF_D1_IK1_RK5() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED && RK4_ENABLED
        module function test_setECDF_D1_IK1_RK4() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED && RK3_ENABLED
        module function test_setECDF_D1_IK1_RK3() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED && RK2_ENABLED
        module function test_setECDF_D1_IK1_RK2() result(assertion); logical(LK) :: assertion; end function
#endif
#if     IK1_ENABLED && RK1_ENABLED
        module function test_setECDF_D1_IK1_RK1() result(assertion); logical(LK) :: assertion; end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK3_ENABLED
        call test%run(test_setECDF_D1_RK3_RK3 , SK_"test_setECDF_D1_RK3_RK3")
#endif
#if     RK2_ENABLED
        call test%run(test_setECDF_D1_RK2_RK2   , SK_"test_setECDF_D1_RK2_RK2")
#endif
#if     RK1_ENABLED
        call test%run(test_setECDF_D1_RK1_RK1   , SK_"test_setECDF_D1_RK1_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK4_ENABLED && RK3_ENABLED
        call test%run(test_setECDF_D1_IK4_RK3  , SK_"test_setECDF_D1_IK4_RK3")
#endif
#if     IK4_ENABLED && RK2_ENABLED
        call test%run(test_setECDF_D1_IK4_RK2   , SK_"test_setECDF_D1_IK4_RK2")
#endif
#if     IK4_ENABLED && RK1_ENABLED
        call test%run(test_setECDF_D1_IK4_RK1   , SK_"test_setECDF_D1_IK4_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK3_ENABLED && RK3_ENABLED
        call test%run(test_setECDF_D1_IK3_RK3  , SK_"test_setECDF_D1_IK3_RK3")
#endif
#if     IK3_ENABLED && RK2_ENABLED
        call test%run(test_setECDF_D1_IK3_RK2   , SK_"test_setECDF_D1_IK3_RK2")
#endif
#if     IK3_ENABLED && RK1_ENABLED
        call test%run(test_setECDF_D1_IK3_RK1   , SK_"test_setECDF_D1_IK3_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK2_ENABLED && RK3_ENABLED
        call test%run(test_setECDF_D1_IK2_RK3  , SK_"test_setECDF_D1_IK2_RK3")
#endif
#if     IK2_ENABLED && RK2_ENABLED
        call test%run(test_setECDF_D1_IK2_RK2   , SK_"test_setECDF_D1_IK2_RK2")
#endif
#if     IK2_ENABLED && RK1_ENABLED
        call test%run(test_setECDF_D1_IK2_RK1   , SK_"test_setECDF_D1_IK2_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK1_ENABLED && RK3_ENABLED
        call test%run(test_setECDF_D1_IK1_RK3   , SK_"test_setECDF_D1_IK1_RK3")
#endif
#if     IK1_ENABLED && RK2_ENABLED
        call test%run(test_setECDF_D1_IK1_RK2    , SK_"test_setECDF_D1_IK1_RK2")
#endif
#if     IK1_ENABLED && RK1_ENABLED
        call test%run(test_setECDF_D1_IK1_RK1    , SK_"test_setECDF_D1_IK1_RK1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_sampleECDF