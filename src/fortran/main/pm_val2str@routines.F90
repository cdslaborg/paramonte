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
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_val2str) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,SK_"@file::"//__FILE__//SK_"@line::"//getStr(LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_arrayResize, only: setResized

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getStr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getStr_D0_SK5_SK
        use pm_kind, only: SKO => SK, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getStr_D0_SK4_SK
        use pm_kind, only: SKO => SK, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getStr_D0_SK3_SK
        use pm_kind, only: SKO => SK, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getStr_D0_SK2_SK
        use pm_kind, only: SKO => SK, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getStr_D0_SK1_SK
        use pm_kind, only: SKO => SK, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getStr_D0_IK5_SK
        use pm_kind, only: SKO => SK, IKG => IK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getStr_D0_IK4_SK
        use pm_kind, only: SKO => SK, IKG => IK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getStr_D0_IK3_SK
        use pm_kind, only: SKO => SK, IKG => IK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getStr_D0_IK2_SK
        use pm_kind, only: SKO => SK, IKG => IK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getStr_D0_IK1_SK
        use pm_kind, only: SKO => SK, IKG => IK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getStr_D0_LK5_SK
        use pm_kind, only: SKO => SK, LKG => LK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getStr_D0_LK4_SK
        use pm_kind, only: SKO => SK, LKG => LK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getStr_D0_LK3_SK
        use pm_kind, only: SKO => SK, LKG => LK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getStr_D0_LK2_SK
        use pm_kind, only: SKO => SK, LKG => LK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getStr_D0_LK1_SK
        use pm_kind, only: SKO => SK, LKG => LK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getStr_D0_CK5_SK
        use pm_kind, only: SKO => SK, CKG => CK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getStr_D0_CK4_SK
        use pm_kind, only: SKO => SK, CKG => CK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getStr_D0_CK3_SK
        use pm_kind, only: SKO => SK, CKG => CK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getStr_D0_CK2_SK
        use pm_kind, only: SKO => SK, CKG => CK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getStr_D0_CK1_SK
        use pm_kind, only: SKO => SK, CKG => CK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getStr_D0_RK5_SK
        use pm_kind, only: SKO => SK, RKG => RK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getStr_D0_RK4_SK
        use pm_kind, only: SKO => SK, RKG => RK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getStr_D0_RK3_SK
        use pm_kind, only: SKO => SK, RKG => RK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getStr_D0_RK2_SK
        use pm_kind, only: SKO => SK, RKG => RK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getStr_D0_RK1_SK
        use pm_kind, only: SKO => SK, RKG => RK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getStr_D0_PSSK5_SK
        use pm_kind, only: SKO => SK, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getStr_D0_PSSK4_SK
        use pm_kind, only: SKO => SK, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getStr_D0_PSSK3_SK
        use pm_kind, only: SKO => SK, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getStr_D0_PSSK2_SK
        use pm_kind, only: SKO => SK, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getStr_D0_PSSK1_SK
        use pm_kind, only: SKO => SK, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1
    module procedure getStr_D0_BSSK_SK
        use pm_kind, only: SKO => SK, SKG => SK
#include "pm_val2str@routines.inc.F90"
    end procedure
#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getStr_D1_SK5_SK
        use pm_kind, only: SKO => SK, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getStr_D1_SK4_SK
        use pm_kind, only: SKO => SK, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getStr_D1_SK3_SK
        use pm_kind, only: SKO => SK, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getStr_D1_SK2_SK
        use pm_kind, only: SKO => SK, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getStr_D1_SK1_SK
        use pm_kind, only: SKO => SK, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getStr_D1_IK5_SK
        use pm_kind, only: SKO => SK, IKG => IK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getStr_D1_IK4_SK
        use pm_kind, only: SKO => SK, IKG => IK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getStr_D1_IK3_SK
        use pm_kind, only: SKO => SK, IKG => IK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getStr_D1_IK2_SK
        use pm_kind, only: SKO => SK, IKG => IK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getStr_D1_IK1_SK
        use pm_kind, only: SKO => SK, IKG => IK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getStr_D1_LK5_SK
        use pm_kind, only: SKO => SK, LKG => LK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getStr_D1_LK4_SK
        use pm_kind, only: SKO => SK, LKG => LK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getStr_D1_LK3_SK
        use pm_kind, only: SKO => SK, LKG => LK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getStr_D1_LK2_SK
        use pm_kind, only: SKO => SK, LKG => LK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getStr_D1_LK1_SK
        use pm_kind, only: SKO => SK, LKG => LK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getStr_D1_CK5_SK
        use pm_kind, only: SKO => SK, CKG => CK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getStr_D1_CK4_SK
        use pm_kind, only: SKO => SK, CKG => CK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getStr_D1_CK3_SK
        use pm_kind, only: SKO => SK, CKG => CK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getStr_D1_CK2_SK
        use pm_kind, only: SKO => SK, CKG => CK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getStr_D1_CK1_SK
        use pm_kind, only: SKO => SK, CKG => CK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getStr_D1_RK5_SK
        use pm_kind, only: SKO => SK, RKG => RK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getStr_D1_RK4_SK
        use pm_kind, only: SKO => SK, RKG => RK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getStr_D1_RK3_SK
        use pm_kind, only: SKO => SK, RKG => RK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getStr_D1_RK2_SK
        use pm_kind, only: SKO => SK, RKG => RK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getStr_D1_RK1_SK
        use pm_kind, only: SKO => SK, RKG => RK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getStr_D1_PSSK5_SK
        use pm_kind, only: SKO => SK, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getStr_D1_PSSK4_SK
        use pm_kind, only: SKO => SK, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getStr_D1_PSSK3_SK
        use pm_kind, only: SKO => SK, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getStr_D1_PSSK2_SK
        use pm_kind, only: SKO => SK, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getStr_D1_PSSK1_SK
        use pm_kind, only: SKO => SK, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1
    module procedure getStr_D1_BSSK_SK
        use pm_kind, only: SKO => SK, SKG => SK
#include "pm_val2str@routines.inc.F90"
    end procedure
#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getStr_D2_SK5_SK
        use pm_kind, only: SKO => SK, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getStr_D2_SK4_SK
        use pm_kind, only: SKO => SK, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getStr_D2_SK3_SK
        use pm_kind, only: SKO => SK, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getStr_D2_SK2_SK
        use pm_kind, only: SKO => SK, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getStr_D2_SK1_SK
        use pm_kind, only: SKO => SK, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getStr_D2_IK5_SK
        use pm_kind, only: SKO => SK, IKG => IK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getStr_D2_IK4_SK
        use pm_kind, only: SKO => SK, IKG => IK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getStr_D2_IK3_SK
        use pm_kind, only: SKO => SK, IKG => IK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getStr_D2_IK2_SK
        use pm_kind, only: SKO => SK, IKG => IK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getStr_D2_IK1_SK
        use pm_kind, only: SKO => SK, IKG => IK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getStr_D2_LK5_SK
        use pm_kind, only: SKO => SK, LKG => LK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getStr_D2_LK4_SK
        use pm_kind, only: SKO => SK, LKG => LK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getStr_D2_LK3_SK
        use pm_kind, only: SKO => SK, LKG => LK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getStr_D2_LK2_SK
        use pm_kind, only: SKO => SK, LKG => LK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getStr_D2_LK1_SK
        use pm_kind, only: SKO => SK, LKG => LK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getStr_D2_CK5_SK
        use pm_kind, only: SKO => SK, CKG => CK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getStr_D2_CK4_SK
        use pm_kind, only: SKO => SK, CKG => CK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getStr_D2_CK3_SK
        use pm_kind, only: SKO => SK, CKG => CK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getStr_D2_CK2_SK
        use pm_kind, only: SKO => SK, CKG => CK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getStr_D2_CK1_SK
        use pm_kind, only: SKO => SK, CKG => CK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getStr_D2_RK5_SK
        use pm_kind, only: SKO => SK, RKG => RK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getStr_D2_RK4_SK
        use pm_kind, only: SKO => SK, RKG => RK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getStr_D2_RK3_SK
        use pm_kind, only: SKO => SK, RKG => RK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getStr_D2_RK2_SK
        use pm_kind, only: SKO => SK, RKG => RK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getStr_D2_RK1_SK
        use pm_kind, only: SKO => SK, RKG => RK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getStr_D2_PSSK5_SK
        use pm_kind, only: SKO => SK, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getStr_D2_PSSK4_SK
        use pm_kind, only: SKO => SK, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getStr_D2_PSSK3_SK
        use pm_kind, only: SKO => SK, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getStr_D2_PSSK2_SK
        use pm_kind, only: SKO => SK, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getStr_D2_PSSK1_SK
        use pm_kind, only: SKO => SK, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1
    module procedure getStr_D2_BSSK_SK
        use pm_kind, only: SKO => SK, SKG => SK
#include "pm_val2str@routines.inc.F90"
    end procedure
#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getStr_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setStr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK1_ENABLED && SK5_ENABLED
    module procedure setStr_D0_SK5_SK1
        use pm_kind, only: SKO => SK1, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK4_ENABLED
    module procedure setStr_D0_SK4_SK1
        use pm_kind, only: SKO => SK1, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK3_ENABLED
    module procedure setStr_D0_SK3_SK1
        use pm_kind, only: SKO => SK1, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK2_ENABLED
    module procedure setStr_D0_SK2_SK1
        use pm_kind, only: SKO => SK1, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK1_ENABLED
    module procedure setStr_D0_SK1_SK1
        use pm_kind, only: SKO => SK1, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK2_ENABLED && SK5_ENABLED
    module procedure setStr_D0_SK5_SK2
        use pm_kind, only: SKO => SK2, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK4_ENABLED
    module procedure setStr_D0_SK4_SK2
        use pm_kind, only: SKO => SK2, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK3_ENABLED
    module procedure setStr_D0_SK3_SK2
        use pm_kind, only: SKO => SK2, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK2_ENABLED
    module procedure setStr_D0_SK2_SK2
        use pm_kind, only: SKO => SK2, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK1_ENABLED
    module procedure setStr_D0_SK1_SK2
        use pm_kind, only: SKO => SK2, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK3_ENABLED && SK5_ENABLED
    module procedure setStr_D0_SK5_SK3
        use pm_kind, only: SKO => SK3, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK4_ENABLED
    module procedure setStr_D0_SK4_SK3
        use pm_kind, only: SKO => SK3, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK3_ENABLED
    module procedure setStr_D0_SK3_SK3
        use pm_kind, only: SKO => SK3, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK2_ENABLED
    module procedure setStr_D0_SK2_SK3
        use pm_kind, only: SKO => SK3, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK1_ENABLED
    module procedure setStr_D0_SK1_SK3
        use pm_kind, only: SKO => SK3, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK4_ENABLED && SK5_ENABLED
    module procedure setStr_D0_SK5_SK4
        use pm_kind, only: SKO => SK4, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK4_ENABLED
    module procedure setStr_D0_SK4_SK4
        use pm_kind, only: SKO => SK4, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK3_ENABLED
    module procedure setStr_D0_SK3_SK4
        use pm_kind, only: SKO => SK4, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK2_ENABLED
    module procedure setStr_D0_SK2_SK4
        use pm_kind, only: SKO => SK4, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK1_ENABLED
    module procedure setStr_D0_SK1_SK4
        use pm_kind, only: SKO => SK4, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED && SK5_ENABLED
    module procedure setStr_D0_SK5_SK5
        use pm_kind, only: SKO => SK5, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK4_ENABLED
    module procedure setStr_D0_SK4_SK5
        use pm_kind, only: SKO => SK5, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK3_ENABLED
    module procedure setStr_D0_SK3_SK5
        use pm_kind, only: SKO => SK5, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK2_ENABLED
    module procedure setStr_D0_SK2_SK5
        use pm_kind, only: SKO => SK5, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK1_ENABLED
    module procedure setStr_D0_SK1_SK5
        use pm_kind, only: SKO => SK5, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if SK1_ENABLED && IK5_ENABLED
    module procedure setStr_D0_IK5_SK1
        use pm_kind, only: SKO => SK1, IKG => IK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && IK4_ENABLED
    module procedure setStr_D0_IK4_SK1
        use pm_kind, only: SKO => SK1, IKG => IK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && IK3_ENABLED
    module procedure setStr_D0_IK3_SK1
        use pm_kind, only: SKO => SK1, IKG => IK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && IK2_ENABLED
    module procedure setStr_D0_IK2_SK1
        use pm_kind, only: SKO => SK1, IKG => IK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && IK1_ENABLED
    module procedure setStr_D0_IK1_SK1
        use pm_kind, only: SKO => SK1, IKG => IK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if SK2_ENABLED && IK5_ENABLED
    module procedure setStr_D0_IK5_SK2
        use pm_kind, only: SKO => SK2, IKG => IK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && IK4_ENABLED
    module procedure setStr_D0_IK4_SK2
        use pm_kind, only: SKO => SK2, IKG => IK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && IK3_ENABLED
    module procedure setStr_D0_IK3_SK2
        use pm_kind, only: SKO => SK2, IKG => IK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && IK2_ENABLED
    module procedure setStr_D0_IK2_SK2
        use pm_kind, only: SKO => SK2, IKG => IK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && IK1_ENABLED
    module procedure setStr_D0_IK1_SK2
        use pm_kind, only: SKO => SK2, IKG => IK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if SK3_ENABLED && IK5_ENABLED
    module procedure setStr_D0_IK5_SK3
        use pm_kind, only: SKO => SK3, IKG => IK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && IK4_ENABLED
    module procedure setStr_D0_IK4_SK3
        use pm_kind, only: SKO => SK3, IKG => IK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && IK3_ENABLED
    module procedure setStr_D0_IK3_SK3
        use pm_kind, only: SKO => SK3, IKG => IK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && IK2_ENABLED
    module procedure setStr_D0_IK2_SK3
        use pm_kind, only: SKO => SK3, IKG => IK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && IK1_ENABLED
    module procedure setStr_D0_IK1_SK3
        use pm_kind, only: SKO => SK3, IKG => IK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if SK4_ENABLED && IK5_ENABLED
    module procedure setStr_D0_IK5_SK4
        use pm_kind, only: SKO => SK4, IKG => IK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && IK4_ENABLED
    module procedure setStr_D0_IK4_SK4
        use pm_kind, only: SKO => SK4, IKG => IK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && IK3_ENABLED
    module procedure setStr_D0_IK3_SK4
        use pm_kind, only: SKO => SK4, IKG => IK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && IK2_ENABLED
    module procedure setStr_D0_IK2_SK4
        use pm_kind, only: SKO => SK4, IKG => IK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && IK1_ENABLED
    module procedure setStr_D0_IK1_SK4
        use pm_kind, only: SKO => SK4, IKG => IK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if SK5_ENABLED && IK5_ENABLED
    module procedure setStr_D0_IK5_SK5
        use pm_kind, only: SKO => SK5, IKG => IK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && IK4_ENABLED
    module procedure setStr_D0_IK4_SK5
        use pm_kind, only: SKO => SK5, IKG => IK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && IK3_ENABLED
    module procedure setStr_D0_IK3_SK5
        use pm_kind, only: SKO => SK5, IKG => IK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && IK2_ENABLED
    module procedure setStr_D0_IK2_SK5
        use pm_kind, only: SKO => SK5, IKG => IK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && IK1_ENABLED
    module procedure setStr_D0_IK1_SK5
        use pm_kind, only: SKO => SK5, IKG => IK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if SK1_ENABLED && LK5_ENABLED
    module procedure setStr_D0_LK5_SK1
        use pm_kind, only: SKO => SK1, LKG => LK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && LK4_ENABLED
    module procedure setStr_D0_LK4_SK1
        use pm_kind, only: SKO => SK1, LKG => LK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && LK3_ENABLED
    module procedure setStr_D0_LK3_SK1
        use pm_kind, only: SKO => SK1, LKG => LK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && LK2_ENABLED
    module procedure setStr_D0_LK2_SK1
        use pm_kind, only: SKO => SK1, LKG => LK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && LK1_ENABLED
    module procedure setStr_D0_LK1_SK1
        use pm_kind, only: SKO => SK1, LKG => LK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if SK2_ENABLED && LK5_ENABLED
    module procedure setStr_D0_LK5_SK2
        use pm_kind, only: SKO => SK2, LKG => LK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && LK4_ENABLED
    module procedure setStr_D0_LK4_SK2
        use pm_kind, only: SKO => SK2, LKG => LK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && LK3_ENABLED
    module procedure setStr_D0_LK3_SK2
        use pm_kind, only: SKO => SK2, LKG => LK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && LK2_ENABLED
    module procedure setStr_D0_LK2_SK2
        use pm_kind, only: SKO => SK2, LKG => LK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && LK1_ENABLED
    module procedure setStr_D0_LK1_SK2
        use pm_kind, only: SKO => SK2, LKG => LK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if SK3_ENABLED && LK5_ENABLED
    module procedure setStr_D0_LK5_SK3
        use pm_kind, only: SKO => SK3, LKG => LK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && LK4_ENABLED
    module procedure setStr_D0_LK4_SK3
        use pm_kind, only: SKO => SK3, LKG => LK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && LK3_ENABLED
    module procedure setStr_D0_LK3_SK3
        use pm_kind, only: SKO => SK3, LKG => LK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && LK2_ENABLED
    module procedure setStr_D0_LK2_SK3
        use pm_kind, only: SKO => SK3, LKG => LK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && LK1_ENABLED
    module procedure setStr_D0_LK1_SK3
        use pm_kind, only: SKO => SK3, LKG => LK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if SK4_ENABLED && LK5_ENABLED
    module procedure setStr_D0_LK5_SK4
        use pm_kind, only: SKO => SK4, LKG => LK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && LK4_ENABLED
    module procedure setStr_D0_LK4_SK4
        use pm_kind, only: SKO => SK4, LKG => LK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && LK3_ENABLED
    module procedure setStr_D0_LK3_SK4
        use pm_kind, only: SKO => SK4, LKG => LK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && LK2_ENABLED
    module procedure setStr_D0_LK2_SK4
        use pm_kind, only: SKO => SK4, LKG => LK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && LK1_ENABLED
    module procedure setStr_D0_LK1_SK4
        use pm_kind, only: SKO => SK4, LKG => LK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if SK5_ENABLED && LK5_ENABLED
    module procedure setStr_D0_LK5_SK5
        use pm_kind, only: SKO => SK5, LKG => LK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && LK4_ENABLED
    module procedure setStr_D0_LK4_SK5
        use pm_kind, only: SKO => SK5, LKG => LK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && LK3_ENABLED
    module procedure setStr_D0_LK3_SK5
        use pm_kind, only: SKO => SK5, LKG => LK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && LK2_ENABLED
    module procedure setStr_D0_LK2_SK5
        use pm_kind, only: SKO => SK5, LKG => LK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && LK1_ENABLED
    module procedure setStr_D0_LK1_SK5
        use pm_kind, only: SKO => SK5, LKG => LK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if SK1_ENABLED && CK5_ENABLED
    module procedure setStr_D0_CK5_SK1
        use pm_kind, only: SKO => SK1, CKG => CK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && CK4_ENABLED
    module procedure setStr_D0_CK4_SK1
        use pm_kind, only: SKO => SK1, CKG => CK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && CK3_ENABLED
    module procedure setStr_D0_CK3_SK1
        use pm_kind, only: SKO => SK1, CKG => CK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && CK2_ENABLED
    module procedure setStr_D0_CK2_SK1
        use pm_kind, only: SKO => SK1, CKG => CK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && CK1_ENABLED
    module procedure setStr_D0_CK1_SK1
        use pm_kind, only: SKO => SK1, CKG => CK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if SK2_ENABLED && CK5_ENABLED
    module procedure setStr_D0_CK5_SK2
        use pm_kind, only: SKO => SK2, CKG => CK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && CK4_ENABLED
    module procedure setStr_D0_CK4_SK2
        use pm_kind, only: SKO => SK2, CKG => CK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && CK3_ENABLED
    module procedure setStr_D0_CK3_SK2
        use pm_kind, only: SKO => SK2, CKG => CK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && CK2_ENABLED
    module procedure setStr_D0_CK2_SK2
        use pm_kind, only: SKO => SK2, CKG => CK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && CK1_ENABLED
    module procedure setStr_D0_CK1_SK2
        use pm_kind, only: SKO => SK2, CKG => CK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if SK3_ENABLED && CK5_ENABLED
    module procedure setStr_D0_CK5_SK3
        use pm_kind, only: SKO => SK3, CKG => CK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && CK4_ENABLED
    module procedure setStr_D0_CK4_SK3
        use pm_kind, only: SKO => SK3, CKG => CK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && CK3_ENABLED
    module procedure setStr_D0_CK3_SK3
        use pm_kind, only: SKO => SK3, CKG => CK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && CK2_ENABLED
    module procedure setStr_D0_CK2_SK3
        use pm_kind, only: SKO => SK3, CKG => CK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && CK1_ENABLED
    module procedure setStr_D0_CK1_SK3
        use pm_kind, only: SKO => SK3, CKG => CK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if SK4_ENABLED && CK5_ENABLED
    module procedure setStr_D0_CK5_SK4
        use pm_kind, only: SKO => SK4, CKG => CK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && CK4_ENABLED
    module procedure setStr_D0_CK4_SK4
        use pm_kind, only: SKO => SK4, CKG => CK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && CK3_ENABLED
    module procedure setStr_D0_CK3_SK4
        use pm_kind, only: SKO => SK4, CKG => CK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && CK2_ENABLED
    module procedure setStr_D0_CK2_SK4
        use pm_kind, only: SKO => SK4, CKG => CK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && CK1_ENABLED
    module procedure setStr_D0_CK1_SK4
        use pm_kind, only: SKO => SK4, CKG => CK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if SK5_ENABLED && CK5_ENABLED
    module procedure setStr_D0_CK5_SK5
        use pm_kind, only: SKO => SK5, CKG => CK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && CK4_ENABLED
    module procedure setStr_D0_CK4_SK5
        use pm_kind, only: SKO => SK5, CKG => CK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && CK3_ENABLED
    module procedure setStr_D0_CK3_SK5
        use pm_kind, only: SKO => SK5, CKG => CK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && CK2_ENABLED
    module procedure setStr_D0_CK2_SK5
        use pm_kind, only: SKO => SK5, CKG => CK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && CK1_ENABLED
    module procedure setStr_D0_CK1_SK5
        use pm_kind, only: SKO => SK5, CKG => CK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if SK1_ENABLED && RK5_ENABLED
    module procedure setStr_D0_RK5_SK1
        use pm_kind, only: SKO => SK1, RKG => RK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && RK4_ENABLED
    module procedure setStr_D0_RK4_SK1
        use pm_kind, only: SKO => SK1, RKG => RK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && RK3_ENABLED
    module procedure setStr_D0_RK3_SK1
        use pm_kind, only: SKO => SK1, RKG => RK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && RK2_ENABLED
    module procedure setStr_D0_RK2_SK1
        use pm_kind, only: SKO => SK1, RKG => RK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && RK1_ENABLED
    module procedure setStr_D0_RK1_SK1
        use pm_kind, only: SKO => SK1, RKG => RK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if SK2_ENABLED && RK5_ENABLED
    module procedure setStr_D0_RK5_SK2
        use pm_kind, only: SKO => SK2, RKG => RK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && RK4_ENABLED
    module procedure setStr_D0_RK4_SK2
        use pm_kind, only: SKO => SK2, RKG => RK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && RK3_ENABLED
    module procedure setStr_D0_RK3_SK2
        use pm_kind, only: SKO => SK2, RKG => RK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && RK2_ENABLED
    module procedure setStr_D0_RK2_SK2
        use pm_kind, only: SKO => SK2, RKG => RK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && RK1_ENABLED
    module procedure setStr_D0_RK1_SK2
        use pm_kind, only: SKO => SK2, RKG => RK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if SK3_ENABLED && RK5_ENABLED
    module procedure setStr_D0_RK5_SK3
        use pm_kind, only: SKO => SK3, RKG => RK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && RK4_ENABLED
    module procedure setStr_D0_RK4_SK3
        use pm_kind, only: SKO => SK3, RKG => RK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && RK3_ENABLED
    module procedure setStr_D0_RK3_SK3
        use pm_kind, only: SKO => SK3, RKG => RK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && RK2_ENABLED
    module procedure setStr_D0_RK2_SK3
        use pm_kind, only: SKO => SK3, RKG => RK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && RK1_ENABLED
    module procedure setStr_D0_RK1_SK3
        use pm_kind, only: SKO => SK3, RKG => RK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if SK4_ENABLED && RK5_ENABLED
    module procedure setStr_D0_RK5_SK4
        use pm_kind, only: SKO => SK4, RKG => RK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && RK4_ENABLED
    module procedure setStr_D0_RK4_SK4
        use pm_kind, only: SKO => SK4, RKG => RK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && RK3_ENABLED
    module procedure setStr_D0_RK3_SK4
        use pm_kind, only: SKO => SK4, RKG => RK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && RK2_ENABLED
    module procedure setStr_D0_RK2_SK4
        use pm_kind, only: SKO => SK4, RKG => RK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && RK1_ENABLED
    module procedure setStr_D0_RK1_SK4
        use pm_kind, only: SKO => SK4, RKG => RK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if SK5_ENABLED && RK5_ENABLED
    module procedure setStr_D0_RK5_SK5
        use pm_kind, only: SKO => SK5, RKG => RK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && RK4_ENABLED
    module procedure setStr_D0_RK4_SK5
        use pm_kind, only: SKO => SK5, RKG => RK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && RK3_ENABLED
    module procedure setStr_D0_RK3_SK5
        use pm_kind, only: SKO => SK5, RKG => RK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && RK2_ENABLED
    module procedure setStr_D0_RK2_SK5
        use pm_kind, only: SKO => SK5, RKG => RK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && RK1_ENABLED
    module procedure setStr_D0_RK1_SK5
        use pm_kind, only: SKO => SK5, RKG => RK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK1_ENABLED && SK5_ENABLED
    module procedure setStr_D0_PSSK5_SK1
        use pm_kind, only: SKO => SK1, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK4_ENABLED
    module procedure setStr_D0_PSSK4_SK1
        use pm_kind, only: SKO => SK1, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK3_ENABLED
    module procedure setStr_D0_PSSK3_SK1
        use pm_kind, only: SKO => SK1, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK2_ENABLED
    module procedure setStr_D0_PSSK2_SK1
        use pm_kind, only: SKO => SK1, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK1_ENABLED
    module procedure setStr_D0_PSSK1_SK1
        use pm_kind, only: SKO => SK1, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK2_ENABLED && SK5_ENABLED
    module procedure setStr_D0_PSSK5_SK2
        use pm_kind, only: SKO => SK2, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK4_ENABLED
    module procedure setStr_D0_PSSK4_SK2
        use pm_kind, only: SKO => SK2, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK3_ENABLED
    module procedure setStr_D0_PSSK3_SK2
        use pm_kind, only: SKO => SK2, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK2_ENABLED
    module procedure setStr_D0_PSSK2_SK2
        use pm_kind, only: SKO => SK2, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK1_ENABLED
    module procedure setStr_D0_PSSK1_SK2
        use pm_kind, only: SKO => SK2, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK3_ENABLED && SK5_ENABLED
    module procedure setStr_D0_PSSK5_SK3
        use pm_kind, only: SKO => SK3, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK4_ENABLED
    module procedure setStr_D0_PSSK4_SK3
        use pm_kind, only: SKO => SK3, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK3_ENABLED
    module procedure setStr_D0_PSSK3_SK3
        use pm_kind, only: SKO => SK3, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK2_ENABLED
    module procedure setStr_D0_PSSK2_SK3
        use pm_kind, only: SKO => SK3, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK1_ENABLED
    module procedure setStr_D0_PSSK1_SK3
        use pm_kind, only: SKO => SK3, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK4_ENABLED && SK5_ENABLED
    module procedure setStr_D0_PSSK5_SK4
        use pm_kind, only: SKO => SK4, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK4_ENABLED
    module procedure setStr_D0_PSSK4_SK4
        use pm_kind, only: SKO => SK4, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK3_ENABLED
    module procedure setStr_D0_PSSK3_SK4
        use pm_kind, only: SKO => SK4, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK2_ENABLED
    module procedure setStr_D0_PSSK2_SK4
        use pm_kind, only: SKO => SK4, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK1_ENABLED
    module procedure setStr_D0_PSSK1_SK4
        use pm_kind, only: SKO => SK4, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED && SK5_ENABLED
    module procedure setStr_D0_PSSK5_SK5
        use pm_kind, only: SKO => SK5, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK4_ENABLED
    module procedure setStr_D0_PSSK4_SK5
        use pm_kind, only: SKO => SK5, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK3_ENABLED
    module procedure setStr_D0_PSSK3_SK5
        use pm_kind, only: SKO => SK5, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK2_ENABLED
    module procedure setStr_D0_PSSK2_SK5
        use pm_kind, only: SKO => SK5, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK1_ENABLED
    module procedure setStr_D0_PSSK1_SK5
        use pm_kind, only: SKO => SK5, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setStr_D0_BSSK5_SK5
        use pm_kind, only: SKO => SK5, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED
    module procedure setStr_D0_BSSK4_SK5
        use pm_kind, only: SKO => SK5, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED
    module procedure setStr_D0_BSSK3_SK5
        use pm_kind, only: SKO => SK5, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED
    module procedure setStr_D0_BSSK2_SK5
        use pm_kind, only: SKO => SK5, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED
    module procedure setStr_D0_BSSK1_SK5
        use pm_kind, only: SKO => SK5, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK1_ENABLED && SK5_ENABLED
    module procedure setStr_D1_SK5_SK1
        use pm_kind, only: SKO => SK1, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK4_ENABLED
    module procedure setStr_D1_SK4_SK1
        use pm_kind, only: SKO => SK1, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK3_ENABLED
    module procedure setStr_D1_SK3_SK1
        use pm_kind, only: SKO => SK1, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK2_ENABLED
    module procedure setStr_D1_SK2_SK1
        use pm_kind, only: SKO => SK1, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK1_ENABLED
    module procedure setStr_D1_SK1_SK1
        use pm_kind, only: SKO => SK1, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK2_ENABLED && SK5_ENABLED
    module procedure setStr_D1_SK5_SK2
        use pm_kind, only: SKO => SK2, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK4_ENABLED
    module procedure setStr_D1_SK4_SK2
        use pm_kind, only: SKO => SK2, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK3_ENABLED
    module procedure setStr_D1_SK3_SK2
        use pm_kind, only: SKO => SK2, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK2_ENABLED
    module procedure setStr_D1_SK2_SK2
        use pm_kind, only: SKO => SK2, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK1_ENABLED
    module procedure setStr_D1_SK1_SK2
        use pm_kind, only: SKO => SK2, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK3_ENABLED && SK5_ENABLED
    module procedure setStr_D1_SK5_SK3
        use pm_kind, only: SKO => SK3, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK4_ENABLED
    module procedure setStr_D1_SK4_SK3
        use pm_kind, only: SKO => SK3, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK3_ENABLED
    module procedure setStr_D1_SK3_SK3
        use pm_kind, only: SKO => SK3, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK2_ENABLED
    module procedure setStr_D1_SK2_SK3
        use pm_kind, only: SKO => SK3, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK1_ENABLED
    module procedure setStr_D1_SK1_SK3
        use pm_kind, only: SKO => SK3, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK4_ENABLED && SK5_ENABLED
    module procedure setStr_D1_SK5_SK4
        use pm_kind, only: SKO => SK4, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK4_ENABLED
    module procedure setStr_D1_SK4_SK4
        use pm_kind, only: SKO => SK4, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK3_ENABLED
    module procedure setStr_D1_SK3_SK4
        use pm_kind, only: SKO => SK4, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK2_ENABLED
    module procedure setStr_D1_SK2_SK4
        use pm_kind, only: SKO => SK4, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK1_ENABLED
    module procedure setStr_D1_SK1_SK4
        use pm_kind, only: SKO => SK4, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED && SK5_ENABLED
    module procedure setStr_D1_SK5_SK5
        use pm_kind, only: SKO => SK5, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK4_ENABLED
    module procedure setStr_D1_SK4_SK5
        use pm_kind, only: SKO => SK5, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK3_ENABLED
    module procedure setStr_D1_SK3_SK5
        use pm_kind, only: SKO => SK5, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK2_ENABLED
    module procedure setStr_D1_SK2_SK5
        use pm_kind, only: SKO => SK5, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK1_ENABLED
    module procedure setStr_D1_SK1_SK5
        use pm_kind, only: SKO => SK5, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if SK1_ENABLED && IK5_ENABLED
    module procedure setStr_D1_IK5_SK1
        use pm_kind, only: SKO => SK1, IKG => IK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && IK4_ENABLED
    module procedure setStr_D1_IK4_SK1
        use pm_kind, only: SKO => SK1, IKG => IK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && IK3_ENABLED
    module procedure setStr_D1_IK3_SK1
        use pm_kind, only: SKO => SK1, IKG => IK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && IK2_ENABLED
    module procedure setStr_D1_IK2_SK1
        use pm_kind, only: SKO => SK1, IKG => IK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && IK1_ENABLED
    module procedure setStr_D1_IK1_SK1
        use pm_kind, only: SKO => SK1, IKG => IK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if SK2_ENABLED && IK5_ENABLED
    module procedure setStr_D1_IK5_SK2
        use pm_kind, only: SKO => SK2, IKG => IK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && IK4_ENABLED
    module procedure setStr_D1_IK4_SK2
        use pm_kind, only: SKO => SK2, IKG => IK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && IK3_ENABLED
    module procedure setStr_D1_IK3_SK2
        use pm_kind, only: SKO => SK2, IKG => IK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && IK2_ENABLED
    module procedure setStr_D1_IK2_SK2
        use pm_kind, only: SKO => SK2, IKG => IK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && IK1_ENABLED
    module procedure setStr_D1_IK1_SK2
        use pm_kind, only: SKO => SK2, IKG => IK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if SK3_ENABLED && IK5_ENABLED
    module procedure setStr_D1_IK5_SK3
        use pm_kind, only: SKO => SK3, IKG => IK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && IK4_ENABLED
    module procedure setStr_D1_IK4_SK3
        use pm_kind, only: SKO => SK3, IKG => IK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && IK3_ENABLED
    module procedure setStr_D1_IK3_SK3
        use pm_kind, only: SKO => SK3, IKG => IK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && IK2_ENABLED
    module procedure setStr_D1_IK2_SK3
        use pm_kind, only: SKO => SK3, IKG => IK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && IK1_ENABLED
    module procedure setStr_D1_IK1_SK3
        use pm_kind, only: SKO => SK3, IKG => IK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if SK4_ENABLED && IK5_ENABLED
    module procedure setStr_D1_IK5_SK4
        use pm_kind, only: SKO => SK4, IKG => IK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && IK4_ENABLED
    module procedure setStr_D1_IK4_SK4
        use pm_kind, only: SKO => SK4, IKG => IK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && IK3_ENABLED
    module procedure setStr_D1_IK3_SK4
        use pm_kind, only: SKO => SK4, IKG => IK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && IK2_ENABLED
    module procedure setStr_D1_IK2_SK4
        use pm_kind, only: SKO => SK4, IKG => IK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && IK1_ENABLED
    module procedure setStr_D1_IK1_SK4
        use pm_kind, only: SKO => SK4, IKG => IK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if SK5_ENABLED && IK5_ENABLED
    module procedure setStr_D1_IK5_SK5
        use pm_kind, only: SKO => SK5, IKG => IK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && IK4_ENABLED
    module procedure setStr_D1_IK4_SK5
        use pm_kind, only: SKO => SK5, IKG => IK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && IK3_ENABLED
    module procedure setStr_D1_IK3_SK5
        use pm_kind, only: SKO => SK5, IKG => IK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && IK2_ENABLED
    module procedure setStr_D1_IK2_SK5
        use pm_kind, only: SKO => SK5, IKG => IK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && IK1_ENABLED
    module procedure setStr_D1_IK1_SK5
        use pm_kind, only: SKO => SK5, IKG => IK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if SK1_ENABLED && LK5_ENABLED
    module procedure setStr_D1_LK5_SK1
        use pm_kind, only: SKO => SK1, LKG => LK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && LK4_ENABLED
    module procedure setStr_D1_LK4_SK1
        use pm_kind, only: SKO => SK1, LKG => LK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && LK3_ENABLED
    module procedure setStr_D1_LK3_SK1
        use pm_kind, only: SKO => SK1, LKG => LK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && LK2_ENABLED
    module procedure setStr_D1_LK2_SK1
        use pm_kind, only: SKO => SK1, LKG => LK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && LK1_ENABLED
    module procedure setStr_D1_LK1_SK1
        use pm_kind, only: SKO => SK1, LKG => LK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if SK2_ENABLED && LK5_ENABLED
    module procedure setStr_D1_LK5_SK2
        use pm_kind, only: SKO => SK2, LKG => LK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && LK4_ENABLED
    module procedure setStr_D1_LK4_SK2
        use pm_kind, only: SKO => SK2, LKG => LK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && LK3_ENABLED
    module procedure setStr_D1_LK3_SK2
        use pm_kind, only: SKO => SK2, LKG => LK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && LK2_ENABLED
    module procedure setStr_D1_LK2_SK2
        use pm_kind, only: SKO => SK2, LKG => LK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && LK1_ENABLED
    module procedure setStr_D1_LK1_SK2
        use pm_kind, only: SKO => SK2, LKG => LK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if SK3_ENABLED && LK5_ENABLED
    module procedure setStr_D1_LK5_SK3
        use pm_kind, only: SKO => SK3, LKG => LK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && LK4_ENABLED
    module procedure setStr_D1_LK4_SK3
        use pm_kind, only: SKO => SK3, LKG => LK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && LK3_ENABLED
    module procedure setStr_D1_LK3_SK3
        use pm_kind, only: SKO => SK3, LKG => LK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && LK2_ENABLED
    module procedure setStr_D1_LK2_SK3
        use pm_kind, only: SKO => SK3, LKG => LK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && LK1_ENABLED
    module procedure setStr_D1_LK1_SK3
        use pm_kind, only: SKO => SK3, LKG => LK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if SK4_ENABLED && LK5_ENABLED
    module procedure setStr_D1_LK5_SK4
        use pm_kind, only: SKO => SK4, LKG => LK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && LK4_ENABLED
    module procedure setStr_D1_LK4_SK4
        use pm_kind, only: SKO => SK4, LKG => LK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && LK3_ENABLED
    module procedure setStr_D1_LK3_SK4
        use pm_kind, only: SKO => SK4, LKG => LK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && LK2_ENABLED
    module procedure setStr_D1_LK2_SK4
        use pm_kind, only: SKO => SK4, LKG => LK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && LK1_ENABLED
    module procedure setStr_D1_LK1_SK4
        use pm_kind, only: SKO => SK4, LKG => LK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if SK5_ENABLED && LK5_ENABLED
    module procedure setStr_D1_LK5_SK5
        use pm_kind, only: SKO => SK5, LKG => LK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && LK4_ENABLED
    module procedure setStr_D1_LK4_SK5
        use pm_kind, only: SKO => SK5, LKG => LK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && LK3_ENABLED
    module procedure setStr_D1_LK3_SK5
        use pm_kind, only: SKO => SK5, LKG => LK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && LK2_ENABLED
    module procedure setStr_D1_LK2_SK5
        use pm_kind, only: SKO => SK5, LKG => LK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && LK1_ENABLED
    module procedure setStr_D1_LK1_SK5
        use pm_kind, only: SKO => SK5, LKG => LK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if SK1_ENABLED && CK5_ENABLED
    module procedure setStr_D1_CK5_SK1
        use pm_kind, only: SKO => SK1, CKG => CK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && CK4_ENABLED
    module procedure setStr_D1_CK4_SK1
        use pm_kind, only: SKO => SK1, CKG => CK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && CK3_ENABLED
    module procedure setStr_D1_CK3_SK1
        use pm_kind, only: SKO => SK1, CKG => CK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && CK2_ENABLED
    module procedure setStr_D1_CK2_SK1
        use pm_kind, only: SKO => SK1, CKG => CK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && CK1_ENABLED
    module procedure setStr_D1_CK1_SK1
        use pm_kind, only: SKO => SK1, CKG => CK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if SK2_ENABLED && CK5_ENABLED
    module procedure setStr_D1_CK5_SK2
        use pm_kind, only: SKO => SK2, CKG => CK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && CK4_ENABLED
    module procedure setStr_D1_CK4_SK2
        use pm_kind, only: SKO => SK2, CKG => CK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && CK3_ENABLED
    module procedure setStr_D1_CK3_SK2
        use pm_kind, only: SKO => SK2, CKG => CK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && CK2_ENABLED
    module procedure setStr_D1_CK2_SK2
        use pm_kind, only: SKO => SK2, CKG => CK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && CK1_ENABLED
    module procedure setStr_D1_CK1_SK2
        use pm_kind, only: SKO => SK2, CKG => CK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if SK3_ENABLED && CK5_ENABLED
    module procedure setStr_D1_CK5_SK3
        use pm_kind, only: SKO => SK3, CKG => CK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && CK4_ENABLED
    module procedure setStr_D1_CK4_SK3
        use pm_kind, only: SKO => SK3, CKG => CK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && CK3_ENABLED
    module procedure setStr_D1_CK3_SK3
        use pm_kind, only: SKO => SK3, CKG => CK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && CK2_ENABLED
    module procedure setStr_D1_CK2_SK3
        use pm_kind, only: SKO => SK3, CKG => CK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && CK1_ENABLED
    module procedure setStr_D1_CK1_SK3
        use pm_kind, only: SKO => SK3, CKG => CK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if SK4_ENABLED && CK5_ENABLED
    module procedure setStr_D1_CK5_SK4
        use pm_kind, only: SKO => SK4, CKG => CK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && CK4_ENABLED
    module procedure setStr_D1_CK4_SK4
        use pm_kind, only: SKO => SK4, CKG => CK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && CK3_ENABLED
    module procedure setStr_D1_CK3_SK4
        use pm_kind, only: SKO => SK4, CKG => CK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && CK2_ENABLED
    module procedure setStr_D1_CK2_SK4
        use pm_kind, only: SKO => SK4, CKG => CK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && CK1_ENABLED
    module procedure setStr_D1_CK1_SK4
        use pm_kind, only: SKO => SK4, CKG => CK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if SK5_ENABLED && CK5_ENABLED
    module procedure setStr_D1_CK5_SK5
        use pm_kind, only: SKO => SK5, CKG => CK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && CK4_ENABLED
    module procedure setStr_D1_CK4_SK5
        use pm_kind, only: SKO => SK5, CKG => CK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && CK3_ENABLED
    module procedure setStr_D1_CK3_SK5
        use pm_kind, only: SKO => SK5, CKG => CK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && CK2_ENABLED
    module procedure setStr_D1_CK2_SK5
        use pm_kind, only: SKO => SK5, CKG => CK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && CK1_ENABLED
    module procedure setStr_D1_CK1_SK5
        use pm_kind, only: SKO => SK5, CKG => CK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if SK1_ENABLED && RK5_ENABLED
    module procedure setStr_D1_RK5_SK1
        use pm_kind, only: SKO => SK1, RKG => RK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && RK4_ENABLED
    module procedure setStr_D1_RK4_SK1
        use pm_kind, only: SKO => SK1, RKG => RK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && RK3_ENABLED
    module procedure setStr_D1_RK3_SK1
        use pm_kind, only: SKO => SK1, RKG => RK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && RK2_ENABLED
    module procedure setStr_D1_RK2_SK1
        use pm_kind, only: SKO => SK1, RKG => RK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && RK1_ENABLED
    module procedure setStr_D1_RK1_SK1
        use pm_kind, only: SKO => SK1, RKG => RK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if SK2_ENABLED && RK5_ENABLED
    module procedure setStr_D1_RK5_SK2
        use pm_kind, only: SKO => SK2, RKG => RK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && RK4_ENABLED
    module procedure setStr_D1_RK4_SK2
        use pm_kind, only: SKO => SK2, RKG => RK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && RK3_ENABLED
    module procedure setStr_D1_RK3_SK2
        use pm_kind, only: SKO => SK2, RKG => RK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && RK2_ENABLED
    module procedure setStr_D1_RK2_SK2
        use pm_kind, only: SKO => SK2, RKG => RK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && RK1_ENABLED
    module procedure setStr_D1_RK1_SK2
        use pm_kind, only: SKO => SK2, RKG => RK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if SK3_ENABLED && RK5_ENABLED
    module procedure setStr_D1_RK5_SK3
        use pm_kind, only: SKO => SK3, RKG => RK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && RK4_ENABLED
    module procedure setStr_D1_RK4_SK3
        use pm_kind, only: SKO => SK3, RKG => RK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && RK3_ENABLED
    module procedure setStr_D1_RK3_SK3
        use pm_kind, only: SKO => SK3, RKG => RK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && RK2_ENABLED
    module procedure setStr_D1_RK2_SK3
        use pm_kind, only: SKO => SK3, RKG => RK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && RK1_ENABLED
    module procedure setStr_D1_RK1_SK3
        use pm_kind, only: SKO => SK3, RKG => RK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if SK4_ENABLED && RK5_ENABLED
    module procedure setStr_D1_RK5_SK4
        use pm_kind, only: SKO => SK4, RKG => RK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && RK4_ENABLED
    module procedure setStr_D1_RK4_SK4
        use pm_kind, only: SKO => SK4, RKG => RK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && RK3_ENABLED
    module procedure setStr_D1_RK3_SK4
        use pm_kind, only: SKO => SK4, RKG => RK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && RK2_ENABLED
    module procedure setStr_D1_RK2_SK4
        use pm_kind, only: SKO => SK4, RKG => RK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && RK1_ENABLED
    module procedure setStr_D1_RK1_SK4
        use pm_kind, only: SKO => SK4, RKG => RK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if SK5_ENABLED && RK5_ENABLED
    module procedure setStr_D1_RK5_SK5
        use pm_kind, only: SKO => SK5, RKG => RK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && RK4_ENABLED
    module procedure setStr_D1_RK4_SK5
        use pm_kind, only: SKO => SK5, RKG => RK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && RK3_ENABLED
    module procedure setStr_D1_RK3_SK5
        use pm_kind, only: SKO => SK5, RKG => RK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && RK2_ENABLED
    module procedure setStr_D1_RK2_SK5
        use pm_kind, only: SKO => SK5, RKG => RK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && RK1_ENABLED
    module procedure setStr_D1_RK1_SK5
        use pm_kind, only: SKO => SK5, RKG => RK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK1_ENABLED && SK5_ENABLED
    module procedure setStr_D1_PSSK5_SK1
        use pm_kind, only: SKO => SK1, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK4_ENABLED
    module procedure setStr_D1_PSSK4_SK1
        use pm_kind, only: SKO => SK1, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK3_ENABLED
    module procedure setStr_D1_PSSK3_SK1
        use pm_kind, only: SKO => SK1, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK2_ENABLED
    module procedure setStr_D1_PSSK2_SK1
        use pm_kind, only: SKO => SK1, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK1_ENABLED
    module procedure setStr_D1_PSSK1_SK1
        use pm_kind, only: SKO => SK1, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK2_ENABLED && SK5_ENABLED
    module procedure setStr_D1_PSSK5_SK2
        use pm_kind, only: SKO => SK2, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK4_ENABLED
    module procedure setStr_D1_PSSK4_SK2
        use pm_kind, only: SKO => SK2, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK3_ENABLED
    module procedure setStr_D1_PSSK3_SK2
        use pm_kind, only: SKO => SK2, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK2_ENABLED
    module procedure setStr_D1_PSSK2_SK2
        use pm_kind, only: SKO => SK2, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK1_ENABLED
    module procedure setStr_D1_PSSK1_SK2
        use pm_kind, only: SKO => SK2, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK3_ENABLED && SK5_ENABLED
    module procedure setStr_D1_PSSK5_SK3
        use pm_kind, only: SKO => SK3, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK4_ENABLED
    module procedure setStr_D1_PSSK4_SK3
        use pm_kind, only: SKO => SK3, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK3_ENABLED
    module procedure setStr_D1_PSSK3_SK3
        use pm_kind, only: SKO => SK3, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK2_ENABLED
    module procedure setStr_D1_PSSK2_SK3
        use pm_kind, only: SKO => SK3, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK1_ENABLED
    module procedure setStr_D1_PSSK1_SK3
        use pm_kind, only: SKO => SK3, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK4_ENABLED && SK5_ENABLED
    module procedure setStr_D1_PSSK5_SK4
        use pm_kind, only: SKO => SK4, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK4_ENABLED
    module procedure setStr_D1_PSSK4_SK4
        use pm_kind, only: SKO => SK4, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK3_ENABLED
    module procedure setStr_D1_PSSK3_SK4
        use pm_kind, only: SKO => SK4, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK2_ENABLED
    module procedure setStr_D1_PSSK2_SK4
        use pm_kind, only: SKO => SK4, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK1_ENABLED
    module procedure setStr_D1_PSSK1_SK4
        use pm_kind, only: SKO => SK4, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED && SK5_ENABLED
    module procedure setStr_D1_PSSK5_SK5
        use pm_kind, only: SKO => SK5, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK4_ENABLED
    module procedure setStr_D1_PSSK4_SK5
        use pm_kind, only: SKO => SK5, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK3_ENABLED
    module procedure setStr_D1_PSSK3_SK5
        use pm_kind, only: SKO => SK5, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK2_ENABLED
    module procedure setStr_D1_PSSK2_SK5
        use pm_kind, only: SKO => SK5, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK1_ENABLED
    module procedure setStr_D1_PSSK1_SK5
        use pm_kind, only: SKO => SK5, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setStr_D1_BSSK5_SK5
        use pm_kind, only: SKO => SK5, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED
    module procedure setStr_D1_BSSK4_SK5
        use pm_kind, only: SKO => SK5, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED
    module procedure setStr_D1_BSSK3_SK5
        use pm_kind, only: SKO => SK5, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED
    module procedure setStr_D1_BSSK2_SK5
        use pm_kind, only: SKO => SK5, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED
    module procedure setStr_D1_BSSK1_SK5
        use pm_kind, only: SKO => SK5, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK1_ENABLED && SK5_ENABLED
    module procedure setStr_D2_SK5_SK1
        use pm_kind, only: SKO => SK1, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK4_ENABLED
    module procedure setStr_D2_SK4_SK1
        use pm_kind, only: SKO => SK1, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK3_ENABLED
    module procedure setStr_D2_SK3_SK1
        use pm_kind, only: SKO => SK1, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK2_ENABLED
    module procedure setStr_D2_SK2_SK1
        use pm_kind, only: SKO => SK1, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK1_ENABLED
    module procedure setStr_D2_SK1_SK1
        use pm_kind, only: SKO => SK1, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK2_ENABLED && SK5_ENABLED
    module procedure setStr_D2_SK5_SK2
        use pm_kind, only: SKO => SK2, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK4_ENABLED
    module procedure setStr_D2_SK4_SK2
        use pm_kind, only: SKO => SK2, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK3_ENABLED
    module procedure setStr_D2_SK3_SK2
        use pm_kind, only: SKO => SK2, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK2_ENABLED
    module procedure setStr_D2_SK2_SK2
        use pm_kind, only: SKO => SK2, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK1_ENABLED
    module procedure setStr_D2_SK1_SK2
        use pm_kind, only: SKO => SK2, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK3_ENABLED && SK5_ENABLED
    module procedure setStr_D2_SK5_SK3
        use pm_kind, only: SKO => SK3, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK4_ENABLED
    module procedure setStr_D2_SK4_SK3
        use pm_kind, only: SKO => SK3, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK3_ENABLED
    module procedure setStr_D2_SK3_SK3
        use pm_kind, only: SKO => SK3, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK2_ENABLED
    module procedure setStr_D2_SK2_SK3
        use pm_kind, only: SKO => SK3, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK1_ENABLED
    module procedure setStr_D2_SK1_SK3
        use pm_kind, only: SKO => SK3, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK4_ENABLED && SK5_ENABLED
    module procedure setStr_D2_SK5_SK4
        use pm_kind, only: SKO => SK4, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK4_ENABLED
    module procedure setStr_D2_SK4_SK4
        use pm_kind, only: SKO => SK4, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK3_ENABLED
    module procedure setStr_D2_SK3_SK4
        use pm_kind, only: SKO => SK4, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK2_ENABLED
    module procedure setStr_D2_SK2_SK4
        use pm_kind, only: SKO => SK4, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK1_ENABLED
    module procedure setStr_D2_SK1_SK4
        use pm_kind, only: SKO => SK4, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED && SK5_ENABLED
    module procedure setStr_D2_SK5_SK5
        use pm_kind, only: SKO => SK5, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK4_ENABLED
    module procedure setStr_D2_SK4_SK5
        use pm_kind, only: SKO => SK5, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK3_ENABLED
    module procedure setStr_D2_SK3_SK5
        use pm_kind, only: SKO => SK5, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK2_ENABLED
    module procedure setStr_D2_SK2_SK5
        use pm_kind, only: SKO => SK5, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK1_ENABLED
    module procedure setStr_D2_SK1_SK5
        use pm_kind, only: SKO => SK5, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if SK1_ENABLED && IK5_ENABLED
    module procedure setStr_D2_IK5_SK1
        use pm_kind, only: SKO => SK1, IKG => IK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && IK4_ENABLED
    module procedure setStr_D2_IK4_SK1
        use pm_kind, only: SKO => SK1, IKG => IK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && IK3_ENABLED
    module procedure setStr_D2_IK3_SK1
        use pm_kind, only: SKO => SK1, IKG => IK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && IK2_ENABLED
    module procedure setStr_D2_IK2_SK1
        use pm_kind, only: SKO => SK1, IKG => IK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && IK1_ENABLED
    module procedure setStr_D2_IK1_SK1
        use pm_kind, only: SKO => SK1, IKG => IK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if SK2_ENABLED && IK5_ENABLED
    module procedure setStr_D2_IK5_SK2
        use pm_kind, only: SKO => SK2, IKG => IK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && IK4_ENABLED
    module procedure setStr_D2_IK4_SK2
        use pm_kind, only: SKO => SK2, IKG => IK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && IK3_ENABLED
    module procedure setStr_D2_IK3_SK2
        use pm_kind, only: SKO => SK2, IKG => IK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && IK2_ENABLED
    module procedure setStr_D2_IK2_SK2
        use pm_kind, only: SKO => SK2, IKG => IK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && IK1_ENABLED
    module procedure setStr_D2_IK1_SK2
        use pm_kind, only: SKO => SK2, IKG => IK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if SK3_ENABLED && IK5_ENABLED
    module procedure setStr_D2_IK5_SK3
        use pm_kind, only: SKO => SK3, IKG => IK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && IK4_ENABLED
    module procedure setStr_D2_IK4_SK3
        use pm_kind, only: SKO => SK3, IKG => IK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && IK3_ENABLED
    module procedure setStr_D2_IK3_SK3
        use pm_kind, only: SKO => SK3, IKG => IK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && IK2_ENABLED
    module procedure setStr_D2_IK2_SK3
        use pm_kind, only: SKO => SK3, IKG => IK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && IK1_ENABLED
    module procedure setStr_D2_IK1_SK3
        use pm_kind, only: SKO => SK3, IKG => IK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if SK4_ENABLED && IK5_ENABLED
    module procedure setStr_D2_IK5_SK4
        use pm_kind, only: SKO => SK4, IKG => IK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && IK4_ENABLED
    module procedure setStr_D2_IK4_SK4
        use pm_kind, only: SKO => SK4, IKG => IK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && IK3_ENABLED
    module procedure setStr_D2_IK3_SK4
        use pm_kind, only: SKO => SK4, IKG => IK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && IK2_ENABLED
    module procedure setStr_D2_IK2_SK4
        use pm_kind, only: SKO => SK4, IKG => IK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && IK1_ENABLED
    module procedure setStr_D2_IK1_SK4
        use pm_kind, only: SKO => SK4, IKG => IK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if SK5_ENABLED && IK5_ENABLED
    module procedure setStr_D2_IK5_SK5
        use pm_kind, only: SKO => SK5, IKG => IK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && IK4_ENABLED
    module procedure setStr_D2_IK4_SK5
        use pm_kind, only: SKO => SK5, IKG => IK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && IK3_ENABLED
    module procedure setStr_D2_IK3_SK5
        use pm_kind, only: SKO => SK5, IKG => IK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && IK2_ENABLED
    module procedure setStr_D2_IK2_SK5
        use pm_kind, only: SKO => SK5, IKG => IK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && IK1_ENABLED
    module procedure setStr_D2_IK1_SK5
        use pm_kind, only: SKO => SK5, IKG => IK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if SK1_ENABLED && LK5_ENABLED
    module procedure setStr_D2_LK5_SK1
        use pm_kind, only: SKO => SK1, LKG => LK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && LK4_ENABLED
    module procedure setStr_D2_LK4_SK1
        use pm_kind, only: SKO => SK1, LKG => LK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && LK3_ENABLED
    module procedure setStr_D2_LK3_SK1
        use pm_kind, only: SKO => SK1, LKG => LK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && LK2_ENABLED
    module procedure setStr_D2_LK2_SK1
        use pm_kind, only: SKO => SK1, LKG => LK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && LK1_ENABLED
    module procedure setStr_D2_LK1_SK1
        use pm_kind, only: SKO => SK1, LKG => LK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if SK2_ENABLED && LK5_ENABLED
    module procedure setStr_D2_LK5_SK2
        use pm_kind, only: SKO => SK2, LKG => LK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && LK4_ENABLED
    module procedure setStr_D2_LK4_SK2
        use pm_kind, only: SKO => SK2, LKG => LK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && LK3_ENABLED
    module procedure setStr_D2_LK3_SK2
        use pm_kind, only: SKO => SK2, LKG => LK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && LK2_ENABLED
    module procedure setStr_D2_LK2_SK2
        use pm_kind, only: SKO => SK2, LKG => LK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && LK1_ENABLED
    module procedure setStr_D2_LK1_SK2
        use pm_kind, only: SKO => SK2, LKG => LK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if SK3_ENABLED && LK5_ENABLED
    module procedure setStr_D2_LK5_SK3
        use pm_kind, only: SKO => SK3, LKG => LK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && LK4_ENABLED
    module procedure setStr_D2_LK4_SK3
        use pm_kind, only: SKO => SK3, LKG => LK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && LK3_ENABLED
    module procedure setStr_D2_LK3_SK3
        use pm_kind, only: SKO => SK3, LKG => LK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && LK2_ENABLED
    module procedure setStr_D2_LK2_SK3
        use pm_kind, only: SKO => SK3, LKG => LK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && LK1_ENABLED
    module procedure setStr_D2_LK1_SK3
        use pm_kind, only: SKO => SK3, LKG => LK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if SK4_ENABLED && LK5_ENABLED
    module procedure setStr_D2_LK5_SK4
        use pm_kind, only: SKO => SK4, LKG => LK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && LK4_ENABLED
    module procedure setStr_D2_LK4_SK4
        use pm_kind, only: SKO => SK4, LKG => LK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && LK3_ENABLED
    module procedure setStr_D2_LK3_SK4
        use pm_kind, only: SKO => SK4, LKG => LK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && LK2_ENABLED
    module procedure setStr_D2_LK2_SK4
        use pm_kind, only: SKO => SK4, LKG => LK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && LK1_ENABLED
    module procedure setStr_D2_LK1_SK4
        use pm_kind, only: SKO => SK4, LKG => LK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if SK5_ENABLED && LK5_ENABLED
    module procedure setStr_D2_LK5_SK5
        use pm_kind, only: SKO => SK5, LKG => LK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && LK4_ENABLED
    module procedure setStr_D2_LK4_SK5
        use pm_kind, only: SKO => SK5, LKG => LK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && LK3_ENABLED
    module procedure setStr_D2_LK3_SK5
        use pm_kind, only: SKO => SK5, LKG => LK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && LK2_ENABLED
    module procedure setStr_D2_LK2_SK5
        use pm_kind, only: SKO => SK5, LKG => LK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && LK1_ENABLED
    module procedure setStr_D2_LK1_SK5
        use pm_kind, only: SKO => SK5, LKG => LK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if SK1_ENABLED && CK5_ENABLED
    module procedure setStr_D2_CK5_SK1
        use pm_kind, only: SKO => SK1, CKG => CK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && CK4_ENABLED
    module procedure setStr_D2_CK4_SK1
        use pm_kind, only: SKO => SK1, CKG => CK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && CK3_ENABLED
    module procedure setStr_D2_CK3_SK1
        use pm_kind, only: SKO => SK1, CKG => CK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && CK2_ENABLED
    module procedure setStr_D2_CK2_SK1
        use pm_kind, only: SKO => SK1, CKG => CK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && CK1_ENABLED
    module procedure setStr_D2_CK1_SK1
        use pm_kind, only: SKO => SK1, CKG => CK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if SK2_ENABLED && CK5_ENABLED
    module procedure setStr_D2_CK5_SK2
        use pm_kind, only: SKO => SK2, CKG => CK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && CK4_ENABLED
    module procedure setStr_D2_CK4_SK2
        use pm_kind, only: SKO => SK2, CKG => CK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && CK3_ENABLED
    module procedure setStr_D2_CK3_SK2
        use pm_kind, only: SKO => SK2, CKG => CK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && CK2_ENABLED
    module procedure setStr_D2_CK2_SK2
        use pm_kind, only: SKO => SK2, CKG => CK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && CK1_ENABLED
    module procedure setStr_D2_CK1_SK2
        use pm_kind, only: SKO => SK2, CKG => CK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if SK3_ENABLED && CK5_ENABLED
    module procedure setStr_D2_CK5_SK3
        use pm_kind, only: SKO => SK3, CKG => CK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && CK4_ENABLED
    module procedure setStr_D2_CK4_SK3
        use pm_kind, only: SKO => SK3, CKG => CK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && CK3_ENABLED
    module procedure setStr_D2_CK3_SK3
        use pm_kind, only: SKO => SK3, CKG => CK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && CK2_ENABLED
    module procedure setStr_D2_CK2_SK3
        use pm_kind, only: SKO => SK3, CKG => CK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && CK1_ENABLED
    module procedure setStr_D2_CK1_SK3
        use pm_kind, only: SKO => SK3, CKG => CK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if SK4_ENABLED && CK5_ENABLED
    module procedure setStr_D2_CK5_SK4
        use pm_kind, only: SKO => SK4, CKG => CK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && CK4_ENABLED
    module procedure setStr_D2_CK4_SK4
        use pm_kind, only: SKO => SK4, CKG => CK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && CK3_ENABLED
    module procedure setStr_D2_CK3_SK4
        use pm_kind, only: SKO => SK4, CKG => CK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && CK2_ENABLED
    module procedure setStr_D2_CK2_SK4
        use pm_kind, only: SKO => SK4, CKG => CK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && CK1_ENABLED
    module procedure setStr_D2_CK1_SK4
        use pm_kind, only: SKO => SK4, CKG => CK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if SK5_ENABLED && CK5_ENABLED
    module procedure setStr_D2_CK5_SK5
        use pm_kind, only: SKO => SK5, CKG => CK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && CK4_ENABLED
    module procedure setStr_D2_CK4_SK5
        use pm_kind, only: SKO => SK5, CKG => CK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && CK3_ENABLED
    module procedure setStr_D2_CK3_SK5
        use pm_kind, only: SKO => SK5, CKG => CK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && CK2_ENABLED
    module procedure setStr_D2_CK2_SK5
        use pm_kind, only: SKO => SK5, CKG => CK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && CK1_ENABLED
    module procedure setStr_D2_CK1_SK5
        use pm_kind, only: SKO => SK5, CKG => CK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if SK1_ENABLED && RK5_ENABLED
    module procedure setStr_D2_RK5_SK1
        use pm_kind, only: SKO => SK1, RKG => RK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && RK4_ENABLED
    module procedure setStr_D2_RK4_SK1
        use pm_kind, only: SKO => SK1, RKG => RK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && RK3_ENABLED
    module procedure setStr_D2_RK3_SK1
        use pm_kind, only: SKO => SK1, RKG => RK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && RK2_ENABLED
    module procedure setStr_D2_RK2_SK1
        use pm_kind, only: SKO => SK1, RKG => RK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && RK1_ENABLED
    module procedure setStr_D2_RK1_SK1
        use pm_kind, only: SKO => SK1, RKG => RK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if SK2_ENABLED && RK5_ENABLED
    module procedure setStr_D2_RK5_SK2
        use pm_kind, only: SKO => SK2, RKG => RK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && RK4_ENABLED
    module procedure setStr_D2_RK4_SK2
        use pm_kind, only: SKO => SK2, RKG => RK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && RK3_ENABLED
    module procedure setStr_D2_RK3_SK2
        use pm_kind, only: SKO => SK2, RKG => RK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && RK2_ENABLED
    module procedure setStr_D2_RK2_SK2
        use pm_kind, only: SKO => SK2, RKG => RK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && RK1_ENABLED
    module procedure setStr_D2_RK1_SK2
        use pm_kind, only: SKO => SK2, RKG => RK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if SK3_ENABLED && RK5_ENABLED
    module procedure setStr_D2_RK5_SK3
        use pm_kind, only: SKO => SK3, RKG => RK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && RK4_ENABLED
    module procedure setStr_D2_RK4_SK3
        use pm_kind, only: SKO => SK3, RKG => RK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && RK3_ENABLED
    module procedure setStr_D2_RK3_SK3
        use pm_kind, only: SKO => SK3, RKG => RK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && RK2_ENABLED
    module procedure setStr_D2_RK2_SK3
        use pm_kind, only: SKO => SK3, RKG => RK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && RK1_ENABLED
    module procedure setStr_D2_RK1_SK3
        use pm_kind, only: SKO => SK3, RKG => RK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if SK4_ENABLED && RK5_ENABLED
    module procedure setStr_D2_RK5_SK4
        use pm_kind, only: SKO => SK4, RKG => RK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && RK4_ENABLED
    module procedure setStr_D2_RK4_SK4
        use pm_kind, only: SKO => SK4, RKG => RK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && RK3_ENABLED
    module procedure setStr_D2_RK3_SK4
        use pm_kind, only: SKO => SK4, RKG => RK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && RK2_ENABLED
    module procedure setStr_D2_RK2_SK4
        use pm_kind, only: SKO => SK4, RKG => RK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && RK1_ENABLED
    module procedure setStr_D2_RK1_SK4
        use pm_kind, only: SKO => SK4, RKG => RK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if SK5_ENABLED && RK5_ENABLED
    module procedure setStr_D2_RK5_SK5
        use pm_kind, only: SKO => SK5, RKG => RK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && RK4_ENABLED
    module procedure setStr_D2_RK4_SK5
        use pm_kind, only: SKO => SK5, RKG => RK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && RK3_ENABLED
    module procedure setStr_D2_RK3_SK5
        use pm_kind, only: SKO => SK5, RKG => RK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && RK2_ENABLED
    module procedure setStr_D2_RK2_SK5
        use pm_kind, only: SKO => SK5, RKG => RK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && RK1_ENABLED
    module procedure setStr_D2_RK1_SK5
        use pm_kind, only: SKO => SK5, RKG => RK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK1_ENABLED && SK5_ENABLED
    module procedure setStr_D2_PSSK5_SK1
        use pm_kind, only: SKO => SK1, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK4_ENABLED
    module procedure setStr_D2_PSSK4_SK1
        use pm_kind, only: SKO => SK1, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK3_ENABLED
    module procedure setStr_D2_PSSK3_SK1
        use pm_kind, only: SKO => SK1, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK2_ENABLED
    module procedure setStr_D2_PSSK2_SK1
        use pm_kind, only: SKO => SK1, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED && SK1_ENABLED
    module procedure setStr_D2_PSSK1_SK1
        use pm_kind, only: SKO => SK1, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK2_ENABLED && SK5_ENABLED
    module procedure setStr_D2_PSSK5_SK2
        use pm_kind, only: SKO => SK2, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK4_ENABLED
    module procedure setStr_D2_PSSK4_SK2
        use pm_kind, only: SKO => SK2, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK3_ENABLED
    module procedure setStr_D2_PSSK3_SK2
        use pm_kind, only: SKO => SK2, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK2_ENABLED
    module procedure setStr_D2_PSSK2_SK2
        use pm_kind, only: SKO => SK2, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED && SK1_ENABLED
    module procedure setStr_D2_PSSK1_SK2
        use pm_kind, only: SKO => SK2, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK3_ENABLED && SK5_ENABLED
    module procedure setStr_D2_PSSK5_SK3
        use pm_kind, only: SKO => SK3, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK4_ENABLED
    module procedure setStr_D2_PSSK4_SK3
        use pm_kind, only: SKO => SK3, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK3_ENABLED
    module procedure setStr_D2_PSSK3_SK3
        use pm_kind, only: SKO => SK3, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK2_ENABLED
    module procedure setStr_D2_PSSK2_SK3
        use pm_kind, only: SKO => SK3, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED && SK1_ENABLED
    module procedure setStr_D2_PSSK1_SK3
        use pm_kind, only: SKO => SK3, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK4_ENABLED && SK5_ENABLED
    module procedure setStr_D2_PSSK5_SK4
        use pm_kind, only: SKO => SK4, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK4_ENABLED
    module procedure setStr_D2_PSSK4_SK4
        use pm_kind, only: SKO => SK4, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK3_ENABLED
    module procedure setStr_D2_PSSK3_SK4
        use pm_kind, only: SKO => SK4, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK2_ENABLED
    module procedure setStr_D2_PSSK2_SK4
        use pm_kind, only: SKO => SK4, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED && SK1_ENABLED
    module procedure setStr_D2_PSSK1_SK4
        use pm_kind, only: SKO => SK4, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED && SK5_ENABLED
    module procedure setStr_D2_PSSK5_SK5
        use pm_kind, only: SKO => SK5, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK4_ENABLED
    module procedure setStr_D2_PSSK4_SK5
        use pm_kind, only: SKO => SK5, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK3_ENABLED
    module procedure setStr_D2_PSSK3_SK5
        use pm_kind, only: SKO => SK5, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK2_ENABLED
    module procedure setStr_D2_PSSK2_SK5
        use pm_kind, only: SKO => SK5, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED && SK1_ENABLED
    module procedure setStr_D2_PSSK1_SK5
        use pm_kind, only: SKO => SK5, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setStr_D2_BSSK5_SK5
        use pm_kind, only: SKO => SK5, SKG => SK5
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED
    module procedure setStr_D2_BSSK4_SK5
        use pm_kind, only: SKO => SK5, SKG => SK4
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED
    module procedure setStr_D2_BSSK3_SK5
        use pm_kind, only: SKO => SK5, SKG => SK3
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED
    module procedure setStr_D2_BSSK2_SK5
        use pm_kind, only: SKO => SK5, SKG => SK2
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#if SK5_ENABLED
    module procedure setStr_D2_BSSK1_SK5
        use pm_kind, only: SKO => SK5, SKG => SK1
#include "pm_val2str@routines.inc.F90"
    end procedure
#endif

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setStr_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines
