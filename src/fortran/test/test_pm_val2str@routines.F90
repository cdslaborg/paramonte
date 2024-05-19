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
!>  This file contains procedure implementations of [pm_val2Str](@ref pm_val2str).
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

submodule (test_pm_val2str) routines ! LCOV_EXCL_LINE

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getStr_SK_ENABLED 1

    module procedure test_getStr_SK_1
        use pm_kind, only: SK, IK
#include "test_pm_val2str@routines.inc.F90"
    end procedure

#undef getStr_SK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getStr_LK_ENABLED 1

    module procedure test_getStr_LK_1
        use pm_kind, only: IK, LK
#include "test_pm_val2str@routines.inc.F90"
    end procedure

#undef getStr_LK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getStr_RK_ENABLED 1

#if RK3_ENABLED
    module procedure test_getStr_RK3_1
        use pm_kind, only: IK, RK => RK3
#include "test_pm_val2str@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getStr_RK2_1
        use pm_kind, only: IK, RK => RK2
#include "test_pm_val2str@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getStr_RK1_1
        use pm_kind, only: IK, RK => RK1
#include "test_pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef getStr_RK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getStr_CK_ENABLED 1

#if CK3_ENABLED
    module procedure test_getStr_CK3_1
        use pm_kind, only: IK, CK => CK3
#include "test_pm_val2str@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getStr_CK2_1
        use pm_kind, only: IK, CK => CK2
#include "test_pm_val2str@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getStr_CK1_1
        use pm_kind, only: IK, CK => CK1
#include "test_pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef getStr_CK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getStr_IK_ENABLED 1

#if IK4_ENABLED
    module procedure test_getStr_IK4_1
        use pm_kind, only: IK, IKG => IK4
#include "test_pm_val2str@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getStr_IK3_1
        use pm_kind, only: IK, IKG => IK3
#include "test_pm_val2str@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getStr_IK2_1
        use pm_kind, only: IK, IKG => IK2
#include "test_pm_val2str@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getStr_IK1_1
        use pm_kind, only: IK, IKG => IK1
#include "test_pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef getStr_IK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines
