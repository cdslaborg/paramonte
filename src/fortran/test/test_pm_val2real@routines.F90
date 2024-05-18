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

!>  \brief This file contains the implementations of the tests of module [pm_val2Real](@ref pm_val2real).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi

submodule (test_pm_val2real) routines

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED
    module procedure test_getReal128_LK_1
        use pm_val2real, only: getReal => getReal128
        use pm_kind, only: IK, RK => RK3
#define test_getReal128_LK_ENABLED 1
#include "test_pm_val2real@routines.inc.F90"
#undef test_getReal128_LK_ENABLED
    end procedure

    module procedure test_getReal128_SK_1
        use pm_val2real, only: getReal => getReal128
        use pm_kind, only: IK, RK => RK3
#define test_getReal128_SK_ENABLED 1
#include "test_pm_val2real@routines.inc.F90"
#undef test_getReal128_SK_ENABLED
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getReal64_LK_1
        use pm_val2real, only: getReal => getReal64
        use pm_kind, only: IK, RK => RK2
#define test_getReal64_LK_ENABLED 1
#include "test_pm_val2real@routines.inc.F90"
#undef test_getReal64_LK_ENABLED
    end procedure

    module procedure test_getReal64_SK_1
        use pm_val2real, only: getReal => getReal64
        use pm_kind, only: IK, RK => RK2
#define test_getReal64_SK_ENABLED 1
#include "test_pm_val2real@routines.inc.F90"
#undef test_getReal64_SK_ENABLED
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getReal32_LK_1
        use pm_val2real, only: getReal => getReal32
        use pm_kind, only: IK, RK => RK1
#define test_getReal32_LK_ENABLED 1
#include "test_pm_val2real@routines.inc.F90"
#undef test_getReal32_LK_ENABLED
    end procedure

    module procedure test_getReal32_SK_1
        use pm_val2real, only: getReal => getReal32
        use pm_kind, only: IK, RK => RK1
#define test_getReal32_SK_ENABLED 1
#include "test_pm_val2real@routines.inc.F90"
#undef test_getReal32_SK_ENABLED
    end procedure
#endif

    module procedure test_getReal_LK_1
        use pm_kind, only: IK, RK
        use pm_val2real, only: getReal
#define test_getReal_LK_ENABLED 1
#include "test_pm_val2real@routines.inc.F90"
#undef test_getReal_LK_ENABLED
    end procedure

    module procedure test_getReal_SK_1
        use pm_kind, only: IK, RK
        use pm_val2real, only: getReal
#define test_getReal_SK_ENABLED 1
#include "test_pm_val2real@routines.inc.F90"
#undef test_getReal_SK_ENABLED
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE
