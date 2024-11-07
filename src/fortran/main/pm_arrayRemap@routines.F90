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
!>  This file contains procedure implementations of [pm_arrayRemap](@ref pm_arrayRemap).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_arrayRemap) routines ! LCOV_EXCL_LINE

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

#define getRemapped_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define For_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getRemappedFor_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRemappedFor_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRemappedFor_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRemappedFor_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRemappedFor_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRemap@routines.inc.F90"
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
    module procedure getRemappedFor_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRemappedFor_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRemappedFor_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRemappedFor_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRemappedFor_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRemappedFor_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRemappedFor_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRemappedFor_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRemappedFor_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRemappedFor_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getRemappedFor_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getRemappedFor_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getRemappedFor_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getRemappedFor_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getRemappedFor_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getRemappedFor_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getRemappedFor_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getRemappedFor_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getRemappedFor_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getRemappedFor_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRemappedFor_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRemappedFor_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRemappedFor_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRemappedFor_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRemappedFor_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef For_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Rev_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getRemappedRev_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRemappedRev_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRemappedRev_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRemappedRev_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRemappedRev_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRemap@routines.inc.F90"
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
    module procedure getRemappedRev_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRemappedRev_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRemappedRev_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRemappedRev_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRemappedRev_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRemappedRev_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRemappedRev_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRemappedRev_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRemappedRev_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRemappedRev_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getRemappedRev_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getRemappedRev_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getRemappedRev_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getRemappedRev_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getRemappedRev_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getRemappedRev_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getRemappedRev_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getRemappedRev_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getRemappedRev_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getRemappedRev_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRemappedRev_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRemappedRev_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRemappedRev_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRemappedRev_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRemappedRev_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Rev_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef New_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRemapped_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRemapped_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define For_ENABLED 1
#define Old_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRemappedForOld_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRemappedForOld_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRemappedForOld_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRemappedForOld_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRemappedForOld_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRemap@routines.inc.F90"
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
    module procedure setRemappedForOld_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRemappedForOld_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRemappedForOld_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRemappedForOld_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRemappedForOld_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRemappedForOld_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRemappedForOld_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRemappedForOld_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRemappedForOld_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRemappedForOld_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setRemappedForOld_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setRemappedForOld_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setRemappedForOld_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setRemappedForOld_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setRemappedForOld_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setRemappedForOld_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setRemappedForOld_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setRemappedForOld_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setRemappedForOld_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setRemappedForOld_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRemappedForOld_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRemappedForOld_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRemappedForOld_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRemappedForOld_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRemappedForOld_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef For_ENABLED
#undef Old_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Rev_ENABLED 1
#define Old_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRemappedRevOld_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRemappedRevOld_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRemappedRevOld_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRemappedRevOld_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRemappedRevOld_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRemap@routines.inc.F90"
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
    module procedure setRemappedRevOld_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRemappedRevOld_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRemappedRevOld_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRemappedRevOld_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRemappedRevOld_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRemappedRevOld_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRemappedRevOld_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRemappedRevOld_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRemappedRevOld_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRemappedRevOld_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setRemappedRevOld_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setRemappedRevOld_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setRemappedRevOld_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setRemappedRevOld_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setRemappedRevOld_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setRemappedRevOld_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setRemappedRevOld_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setRemappedRevOld_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setRemappedRevOld_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setRemappedRevOld_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRemappedRevOld_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRemappedRevOld_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRemappedRevOld_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRemappedRevOld_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRemappedRevOld_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Rev_ENABLED
#undef Old_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define For_ENABLED 1
#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRemappedForNew_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRemappedForNew_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRemappedForNew_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRemappedForNew_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRemappedForNew_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRemap@routines.inc.F90"
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
    module procedure setRemappedForNew_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRemappedForNew_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRemappedForNew_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRemappedForNew_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRemappedForNew_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRemappedForNew_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRemappedForNew_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRemappedForNew_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRemappedForNew_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRemappedForNew_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setRemappedForNew_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setRemappedForNew_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setRemappedForNew_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setRemappedForNew_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setRemappedForNew_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setRemappedForNew_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setRemappedForNew_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setRemappedForNew_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setRemappedForNew_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setRemappedForNew_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRemappedForNew_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRemappedForNew_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRemappedForNew_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRemappedForNew_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRemappedForNew_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef For_ENABLED
#undef New_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Rev_ENABLED 1
#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRemappedRevNew_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRemappedRevNew_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRemappedRevNew_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRemappedRevNew_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRemappedRevNew_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRemap@routines.inc.F90"
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
    module procedure setRemappedRevNew_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRemappedRevNew_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRemappedRevNew_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRemappedRevNew_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRemappedRevNew_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRemappedRevNew_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRemappedRevNew_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRemappedRevNew_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRemappedRevNew_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRemappedRevNew_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setRemappedRevNew_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setRemappedRevNew_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setRemappedRevNew_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setRemappedRevNew_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setRemappedRevNew_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setRemappedRevNew_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setRemappedRevNew_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setRemappedRevNew_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setRemappedRevNew_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setRemappedRevNew_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRemappedRevNew_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRemappedRevNew_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRemappedRevNew_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRemappedRevNew_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRemappedRevNew_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_arrayRemap@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Rev_ENABLED
#undef New_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRemapped_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines