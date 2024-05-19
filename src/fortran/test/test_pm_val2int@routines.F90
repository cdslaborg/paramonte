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

!>  \brief This file contains the implementations of the tests of module [val2int_pmod](@ref pm_val2int).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi

submodule (test_pm_val2int) routines

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK4_ENABLED
    module procedure test_getInt64_LK_1
        use pm_val2int, only: getInt => getInt64
        use pm_kind, only: IKG => IK4
#define test_getInt64_LK_ENABLED 1
#include "test_pm_val2int@routines.inc.F90"
#undef test_getInt64_LK_ENABLED
    end procedure

    module procedure test_getInt64_SK_1
        use pm_val2int, only: getInt => getInt64
        use pm_kind, only: IKG => IK4
#define test_getInt64_SK_ENABLED 1
#include "test_pm_val2int@routines.inc.F90"
#undef test_getInt64_SK_ENABLED
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK3_ENABLED
    module procedure test_getInt32_LK_1
        use pm_val2int, only: getInt => getInt32
        use pm_kind, only: IKG => IK3
#define test_getInt32_LK_ENABLED 1
#include "test_pm_val2int@routines.inc.F90"
#undef test_getInt32_LK_ENABLED
    end procedure

    module procedure test_getInt32_SK_1
        use pm_val2int, only: getInt => getInt32
        use pm_kind, only: IKG => IK3
#define test_getInt32_SK_ENABLED 1
#include "test_pm_val2int@routines.inc.F90"
#undef test_getInt32_SK_ENABLED
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK2_ENABLED
    module procedure test_getInt16_LK_1
        use pm_val2int, only: getInt => getInt16
        use pm_kind, only: IKG => IK2
#define test_getInt16_LK_ENABLED 1
#include "test_pm_val2int@routines.inc.F90"
#undef test_getInt16_LK_ENABLED
    end procedure

    module procedure test_getInt16_SK_1
        use pm_val2int, only: getInt => getInt16
        use pm_kind, only: IKG => IK2
#define test_getInt16_SK_ENABLED 1
#include "test_pm_val2int@routines.inc.F90"
#undef test_getInt16_SK_ENABLED
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK1_ENABLED
    module procedure test_getInt8_LK_1
        use pm_val2int, only: getInt => getInt8
        use pm_kind, only: IKG => IK1
#define test_getInt8_LK_ENABLED 1
#include "test_pm_val2int@routines.inc.F90"
#undef test_getInt8_LK_ENABLED
    end procedure

    module procedure test_getInt8_SK_1
        use pm_val2int, only: getInt => getInt8
        use pm_kind, only: IKG => IK1
#define test_getInt8_SK_ENABLED 1
#include "test_pm_val2int@routines.inc.F90"
#undef test_getInt8_SK_ENABLED
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure test_getInt_LK_1
        use pm_kind, only: IKG => IK
        use pm_val2int, only: getInt
#define test_getInt_LK_ENABLED 1
#include "test_pm_val2int@routines.inc.F90"
#undef test_getInt_LK_ENABLED
    end procedure

    module procedure test_getInt_SK_1
        use pm_kind, only: IKG => IK
        use pm_val2int, only: getInt
#define test_getInt_SK_ENABLED 1
#include "test_pm_val2int@routines.inc.F90"
#undef test_getInt_SK_ENABLED
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE
