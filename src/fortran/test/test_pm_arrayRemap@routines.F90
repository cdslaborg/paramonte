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
!>  This file contains procedure implementations of [test_pm_arrayRemap](@ref test_pm_arrayRemap).
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_arrayRemap) routines ! LCOV_EXCL_LINE

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRemapped_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getRemapped_D0_SK5_1
        use pm_kind, only: SKC => SK5
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRemapped_D0_SK4_1
        use pm_kind, only: SKC => SK4
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRemapped_D0_SK3_1
        use pm_kind, only: SKC => SK3
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRemapped_D0_SK2_1
        use pm_kind, only: SKC => SK2
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRemapped_D0_SK1_1
        use pm_kind, only: SKC => SK1
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getRemapped_D1_SK5_1
        use pm_kind, only: SKC => SK5
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRemapped_D1_SK4_1
        use pm_kind, only: SKC => SK4
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRemapped_D1_SK3_1
        use pm_kind, only: SKC => SK3
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRemapped_D1_SK2_1
        use pm_kind, only: SKC => SK2
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRemapped_D1_SK1_1
        use pm_kind, only: SKC => SK1
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getRemapped_D1_IK5_1
        use pm_kind, only: IKC => IK5
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getRemapped_D1_IK4_1
        use pm_kind, only: IKC => IK4
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getRemapped_D1_IK3_1
        use pm_kind, only: IKC => IK3
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getRemapped_D1_IK2_1
        use pm_kind, only: IKC => IK2
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getRemapped_D1_IK1_1
        use pm_kind, only: IKC => IK1
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_getRemapped_D1_LK5_1
        use pm_kind, only: LKC => LK5
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_getRemapped_D1_LK4_1
        use pm_kind, only: LKC => LK4
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_getRemapped_D1_LK3_1
        use pm_kind, only: LKC => LK3
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_getRemapped_D1_LK2_1
        use pm_kind, only: LKC => LK2
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_getRemapped_D1_LK1_1
        use pm_kind, only: LKC => LK1
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getRemapped_D1_CK5_1
        use pm_kind, only: CKC => CK5
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getRemapped_D1_CK4_1
        use pm_kind, only: CKC => CK4
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getRemapped_D1_CK3_1
        use pm_kind, only: CKC => CK3
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getRemapped_D1_CK2_1
        use pm_kind, only: CKC => CK2
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getRemapped_D1_CK1_1
        use pm_kind, only: CKC => CK1
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getRemapped_D1_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getRemapped_D1_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getRemapped_D1_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getRemapped_D1_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getRemapped_D1_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRemapped_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRemapped_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setRemapped_D0_SK5_1
        use pm_kind, only: SKC => SK5
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRemapped_D0_SK4_1
        use pm_kind, only: SKC => SK4
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRemapped_D0_SK3_1
        use pm_kind, only: SKC => SK3
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRemapped_D0_SK2_1
        use pm_kind, only: SKC => SK2
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRemapped_D0_SK1_1
        use pm_kind, only: SKC => SK1
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setRemapped_D1_SK5_1
        use pm_kind, only: SKC => SK5
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRemapped_D1_SK4_1
        use pm_kind, only: SKC => SK4
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRemapped_D1_SK3_1
        use pm_kind, only: SKC => SK3
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRemapped_D1_SK2_1
        use pm_kind, only: SKC => SK2
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRemapped_D1_SK1_1
        use pm_kind, only: SKC => SK1
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setRemapped_D1_IK5_1
        use pm_kind, only: IKC => IK5
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setRemapped_D1_IK4_1
        use pm_kind, only: IKC => IK4
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setRemapped_D1_IK3_1
        use pm_kind, only: IKC => IK3
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setRemapped_D1_IK2_1
        use pm_kind, only: IKC => IK2
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setRemapped_D1_IK1_1
        use pm_kind, only: IKC => IK1
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setRemapped_D1_LK5_1
        use pm_kind, only: LKC => LK5
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setRemapped_D1_LK4_1
        use pm_kind, only: LKC => LK4
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setRemapped_D1_LK3_1
        use pm_kind, only: LKC => LK3
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setRemapped_D1_LK2_1
        use pm_kind, only: LKC => LK2
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setRemapped_D1_LK1_1
        use pm_kind, only: LKC => LK1
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setRemapped_D1_CK5_1
        use pm_kind, only: CKC => CK5
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setRemapped_D1_CK4_1
        use pm_kind, only: CKC => CK4
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setRemapped_D1_CK3_1
        use pm_kind, only: CKC => CK3
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setRemapped_D1_CK2_1
        use pm_kind, only: CKC => CK2
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setRemapped_D1_CK1_1
        use pm_kind, only: CKC => CK1
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setRemapped_D1_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setRemapped_D1_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setRemapped_D1_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setRemapped_D1_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setRemapped_D1_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRemapped_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines
