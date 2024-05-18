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
!>  This file contains procedure implementations of [pm_arrayMerge](@ref pm_arrayMerge).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_arrayMerge) routines ! LCOV_EXCL_LINE

    use pm_distUnif, only: getUnifRand
    use pm_distUnif, only: setUnifRand
    use pm_arraySort, only: setSorted
    use pm_val2str, only: getStr
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMerged_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMerged_D0_SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getMerged_D0_SK5_1
        use pm_kind, only: SKC => SK5
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getMerged_D0_SK4_1
        use pm_kind, only: SKC => SK4
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getMerged_D0_SK3_1
        use pm_kind, only: SKC => SK3
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getMerged_D0_SK2_1
        use pm_kind, only: SKC => SK2
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getMerged_D0_SK1_1
        use pm_kind, only: SKC => SK1
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef getMerged_D0_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMerged_D1_SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getMerged_D1_SK5_1
        use pm_kind, only: SKC => SK5
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getMerged_D1_SK4_1
        use pm_kind, only: SKC => SK4
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getMerged_D1_SK3_1
        use pm_kind, only: SKC => SK3
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getMerged_D1_SK2_1
        use pm_kind, only: SKC => SK2
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getMerged_D1_SK1_1
        use pm_kind, only: SKC => SK1
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef getMerged_D1_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMerged_D1_IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getMerged_D1_IK5_1
        use pm_kind, only: IKC => IK5
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getMerged_D1_IK4_1
        use pm_kind, only: IKC => IK4
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getMerged_D1_IK3_1
        use pm_kind, only: IKC => IK3
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getMerged_D1_IK2_1
        use pm_kind, only: IKC => IK2
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getMerged_D1_IK1_1
        use pm_kind, only: IKC => IK1
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef getMerged_D1_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMerged_D1_LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_getMerged_D1_LK5_1
        use pm_kind, only: LKC => LK5
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_getMerged_D1_LK4_1
        use pm_kind, only: LKC => LK4
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_getMerged_D1_LK3_1
        use pm_kind, only: LKC => LK3
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_getMerged_D1_LK2_1
        use pm_kind, only: LKC => LK2
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_getMerged_D1_LK1_1
        use pm_kind, only: LKC => LK1
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef getMerged_D1_LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMerged_D1_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getMerged_D1_CK5_1
        use pm_kind, only: CKC => CK5
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getMerged_D1_CK4_1
        use pm_kind, only: CKC => CK4
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getMerged_D1_CK3_1
        use pm_kind, only: CKC => CK3
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getMerged_D1_CK2_1
        use pm_kind, only: CKC => CK2
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getMerged_D1_CK1_1
        use pm_kind, only: CKC => CK1
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef getMerged_D1_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMerged_D1_RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getMerged_D1_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getMerged_D1_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getMerged_D1_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getMerged_D1_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getMerged_D1_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef getMerged_D1_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if 0
#define getMerged_D1_PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getMerged_D1_PSSK5_1
        use pm_kind, only: SKC => SK5
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getMerged_D1_PSSK4_1
        use pm_kind, only: SKC => SK4
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getMerged_D1_PSSK3_1
        use pm_kind, only: SKC => SK3
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getMerged_D1_PSSK2_1
        use pm_kind, only: SKC => SK2
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getMerged_D1_PSSK1_1
        use pm_kind, only: SKC => SK1
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef getMerged_D1_PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMerged_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMerged_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMerged_D0_SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setMerged_D0_SK5_1
        use pm_kind, only: SKC => SK5
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setMerged_D0_SK4_1
        use pm_kind, only: SKC => SK4
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setMerged_D0_SK3_1
        use pm_kind, only: SKC => SK3
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setMerged_D0_SK2_1
        use pm_kind, only: SKC => SK2
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setMerged_D0_SK1_1
        use pm_kind, only: SKC => SK1
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef setMerged_D0_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMerged_D1_SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setMerged_D1_SK5_1
        use pm_kind, only: SKC => SK5
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setMerged_D1_SK4_1
        use pm_kind, only: SKC => SK4
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setMerged_D1_SK3_1
        use pm_kind, only: SKC => SK3
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setMerged_D1_SK2_1
        use pm_kind, only: SKC => SK2
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setMerged_D1_SK1_1
        use pm_kind, only: SKC => SK1
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef setMerged_D1_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMerged_D1_IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setMerged_D1_IK5_1
        use pm_kind, only: IKC => IK5
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setMerged_D1_IK4_1
        use pm_kind, only: IKC => IK4
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setMerged_D1_IK3_1
        use pm_kind, only: IKC => IK3
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setMerged_D1_IK2_1
        use pm_kind, only: IKC => IK2
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setMerged_D1_IK1_1
        use pm_kind, only: IKC => IK1
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef setMerged_D1_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMerged_D1_LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setMerged_D1_LK5_1
        use pm_kind, only: LKC => LK5
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setMerged_D1_LK4_1
        use pm_kind, only: LKC => LK4
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setMerged_D1_LK3_1
        use pm_kind, only: LKC => LK3
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setMerged_D1_LK2_1
        use pm_kind, only: LKC => LK2
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setMerged_D1_LK1_1
        use pm_kind, only: LKC => LK1
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef setMerged_D1_LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMerged_D1_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setMerged_D1_CK5_1
        use pm_kind, only: CKC => CK5
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setMerged_D1_CK4_1
        use pm_kind, only: CKC => CK4
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setMerged_D1_CK3_1
        use pm_kind, only: CKC => CK3
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setMerged_D1_CK2_1
        use pm_kind, only: CKC => CK2
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setMerged_D1_CK1_1
        use pm_kind, only: CKC => CK1
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef setMerged_D1_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMerged_D1_RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setMerged_D1_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setMerged_D1_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setMerged_D1_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setMerged_D1_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setMerged_D1_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef setMerged_D1_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if 0
#define setMerged_D1_PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setMerged_D1_PSSK5_1
        use pm_kind, only: SKC => SK5
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setMerged_D1_PSSK4_1
        use pm_kind, only: SKC => SK4
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setMerged_D1_PSSK3_1
        use pm_kind, only: SKC => SK3
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setMerged_D1_PSSK2_1
        use pm_kind, only: SKC => SK2
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setMerged_D1_PSSK1_1
        use pm_kind, only: SKC => SK1
#include "test_pm_arrayMerge@routines.inc.F90"
    end procedure
#endif

#undef setMerged_D1_PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMerged_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines
