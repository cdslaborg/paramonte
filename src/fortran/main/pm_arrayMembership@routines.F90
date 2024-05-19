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
!>  This file contains procedure implementations of [pm_arrayMembership](@ref pm_arrayMembership).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_arrayMembership) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define in_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure in_D0_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure in_D0_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure in_D0_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure in_D0_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure in_D0_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure in_D0_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure in_D0_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure in_D0_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure in_D0_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure in_D0_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure in_D0_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure in_D0_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure in_D0_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure in_D0_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure in_D0_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure in_D0_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure in_D0_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure in_D0_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure in_D0_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure in_D0_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure in_D0_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure in_D0_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure in_D0_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure in_D0_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure in_D0_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure in_D0_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure in_D0_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure in_D0_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure in_D0_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure in_D0_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure in_D1_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure in_D1_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure in_D1_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure in_D1_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure in_D1_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure in_D1_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure in_D1_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure in_D1_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure in_D1_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure in_D1_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure in_D1_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure in_D1_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure in_D1_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure in_D1_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure in_D1_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure in_D1_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure in_D1_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure in_D1_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure in_D1_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure in_D1_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure in_D1_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure in_D1_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure in_D1_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure in_D1_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure in_D1_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef in_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define inrange_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure inrange_D0_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure inrange_D0_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure inrange_D0_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure inrange_D0_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure inrange_D0_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure inrange_D0_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure inrange_D0_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure inrange_D0_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure inrange_D0_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure inrange_D0_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure inrange_D0_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure inrange_D0_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure inrange_D0_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure inrange_D0_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure inrange_D0_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure inrange_D0_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure inrange_D0_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure inrange_D0_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure inrange_D0_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure inrange_D0_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure inrange_D0_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure inrange_D0_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure inrange_D0_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure inrange_D0_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure inrange_D0_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure inrange_D0_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure inrange_D0_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure inrange_D0_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure inrange_D0_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure inrange_D0_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure inrange_D1_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure inrange_D1_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure inrange_D1_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure inrange_D1_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure inrange_D1_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure inrange_D1_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure inrange_D1_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure inrange_D1_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure inrange_D1_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure inrange_D1_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure inrange_D1_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure inrange_D1_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure inrange_D1_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure inrange_D1_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure inrange_D1_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure inrange_D1_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure inrange_D1_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure inrange_D1_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure inrange_D1_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure inrange_D1_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure inrange_D1_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure inrange_D1_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure inrange_D1_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure inrange_D1_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure inrange_D1_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef inrange_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define allin_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure allin_D0_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure allin_D0_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure allin_D0_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure allin_D0_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure allin_D0_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure allin_D1_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure allin_D1_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure allin_D1_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure allin_D1_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure allin_D1_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure allin_D1_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure allin_D1_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure allin_D1_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure allin_D1_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure allin_D1_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure allin_D1_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure allin_D1_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure allin_D1_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure allin_D1_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure allin_D1_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure allin_D1_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure allin_D1_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure allin_D1_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure allin_D1_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure allin_D1_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure allin_D1_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure allin_D1_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure allin_D1_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure allin_D1_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure allin_D1_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef allin_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define allinrange_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure allinrange_D0_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure allinrange_D0_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure allinrange_D0_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure allinrange_D0_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure allinrange_D0_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure allinrange_D1_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure allinrange_D1_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure allinrange_D1_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure allinrange_D1_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure allinrange_D1_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure allinrange_D1_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure allinrange_D1_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure allinrange_D1_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure allinrange_D1_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure allinrange_D1_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure allinrange_D1_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure allinrange_D1_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure allinrange_D1_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure allinrange_D1_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure allinrange_D1_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure allinrange_D1_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure allinrange_D1_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure allinrange_D1_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure allinrange_D1_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure allinrange_D1_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure allinrange_D1_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure allinrange_D1_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure allinrange_D1_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure allinrange_D1_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure allinrange_D1_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef allinrange_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define anyin_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure anyin_D0_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure anyin_D0_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure anyin_D0_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure anyin_D0_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure anyin_D0_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure anyin_D1_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure anyin_D1_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure anyin_D1_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure anyin_D1_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure anyin_D1_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure anyin_D1_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure anyin_D1_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure anyin_D1_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure anyin_D1_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure anyin_D1_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure anyin_D1_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure anyin_D1_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure anyin_D1_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure anyin_D1_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure anyin_D1_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure anyin_D1_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure anyin_D1_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure anyin_D1_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure anyin_D1_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure anyin_D1_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure anyin_D1_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure anyin_D1_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure anyin_D1_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure anyin_D1_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure anyin_D1_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef anyin_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define anyinrange_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure anyinrange_D0_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure anyinrange_D0_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure anyinrange_D0_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure anyinrange_D0_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure anyinrange_D0_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure anyinrange_D1_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure anyinrange_D1_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure anyinrange_D1_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure anyinrange_D1_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure anyinrange_D1_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure anyinrange_D1_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure anyinrange_D1_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure anyinrange_D1_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure anyinrange_D1_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure anyinrange_D1_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure anyinrange_D1_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure anyinrange_D1_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure anyinrange_D1_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure anyinrange_D1_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure anyinrange_D1_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure anyinrange_D1_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure anyinrange_D1_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure anyinrange_D1_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure anyinrange_D1_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure anyinrange_D1_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure anyinrange_D1_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure anyinrange_D1_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure anyinrange_D1_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure anyinrange_D1_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure anyinrange_D1_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayMembership@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef anyinrange_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines