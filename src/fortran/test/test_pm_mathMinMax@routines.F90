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

!>  \brief This file contains the implementations of the tests of module [pm_mathMinMax](@ref pm_mathMinMax).
!>
!>  \author
!>  \AmirShahmoradi

submodule (test_pm_mathMinMax) routines

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMinMax_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMinMax_SK_ENABLED 1

    module procedure test_getMinMax_SK_1
        use pm_kind, only: SK, IK => IK1
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure

#undef getMinMax_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMinMax_CK_ENABLED 1

#if CK3_ENABLED
    module procedure test_getMinMax_CK3_1
        use pm_kind, only: IK, CK => CK3
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getMinMax_CK2_1
        use pm_kind, only: IK, CK => CK2
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getMinMax_CK1_1
        use pm_kind, only: IK, CK => CK1
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef getMinMax_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMinMax_RK_ENABLED 1

#if RK3_ENABLED
    module procedure test_getMinMax_RK3_1
        use pm_kind, only: IK, RK => RK3
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getMinMax_RK2_1
        use pm_kind, only: IK, RK => RK2
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getMinMax_RK1_1
        use pm_kind, only: IK, RK => RK1
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef getMinMax_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMinMax_IK_ENABLED 1

#if IK4_ENABLED
    module procedure test_getMinMax_IK4_1
        use pm_kind, only: IK => IK4
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getMinMax_IK3_1
        use pm_kind, only: IK => IK3
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getMinMax_IK2_1
        use pm_kind, only: IK => IK2
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getMinMax_IK1_1
        use pm_kind, only: IK => IK1
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef getMinMax_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMinMax_SSK_ENABLED 1

    module procedure test_getMinMax_SSK_1
        use pm_kind, only: SK, IK => IK1
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure

#undef getMinMax_SSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMinMax_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMinMax_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMinMax_SK_ENABLED 1

    module procedure test_setMinMax_SK_1
        use pm_kind, only: SK, IK => IK1
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure

#undef setMinMax_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMinMax_CK_ENABLED 1

#if CK3_ENABLED
    module procedure test_setMinMax_CK3_1
        use pm_kind, only: IK, CK => CK3
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setMinMax_CK2_1
        use pm_kind, only: IK, CK => CK2
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setMinMax_CK1_1
        use pm_kind, only: IK, CK => CK1
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef setMinMax_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMinMax_RK_ENABLED 1

#if RK3_ENABLED
    module procedure test_setMinMax_RK3_1
        use pm_kind, only: IK, RK => RK3
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setMinMax_RK2_1
        use pm_kind, only: IK, RK => RK2
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setMinMax_RK1_1
        use pm_kind, only: IK, RK => RK1
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef setMinMax_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMinMax_IK_ENABLED 1

#if IK4_ENABLED
    module procedure test_setMinMax_IK4_1
        use pm_kind, only: IK => IK4
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setMinMax_IK3_1
        use pm_kind, only: IK => IK3
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setMinMax_IK2_1
        use pm_kind, only: IK => IK2
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setMinMax_IK1_1
        use pm_kind, only: IK => IK1
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef setMinMax_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMinMax_SSK_ENABLED 1

    module procedure test_setMinMax_SSK_1
        use pm_kind, only: SK, IK => IK1
#include "test_pm_mathMinMax@routines.inc.F90"
    end procedure

#undef setMinMax_SSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMinMax_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE
