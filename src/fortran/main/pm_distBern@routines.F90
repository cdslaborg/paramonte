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
!>  This file contains procedure implementations of [pm_distBern](@ref pm_distBern).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distBern) routines ! LCOV_EXCL_LINE

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

#define isHead_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    module procedure isHeadDD
        rand = isHead(p = .5)
    end procedure

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DS_ENABLED 1

    module procedure isHeadDS
        rand = isHead(p = .5, size = size)
    end procedure

#undef DS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PD_ENABLED 1

#if RK5_ENABLED
    module procedure isHeadPD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isHeadPD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isHeadPD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isHeadPD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isHeadPD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#undef PD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PS_ENABLED 1

#if RK5_ENABLED
    module procedure isHeadPS_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isHeadPS_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isHeadPS_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isHeadPS_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isHeadPS_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#undef PS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isHead_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getBernRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getBernRandPD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getBernRandPD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getBernRandPD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getBernRandPD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getBernRandPD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getBernRandPS_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getBernRandPS_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getBernRandPS_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getBernRandPS_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getBernRandPS_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getBernRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setBernRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RUP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED && RK5_ENABLED
    module procedure setBernRandRUP_IK5_RK5
        use pm_kind, only: IKC => IK5, RKC => RK5
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && RK4_ENABLED
    module procedure setBernRandRUP_IK5_RK4
        use pm_kind, only: IKC => IK5, RKC => RK4
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && RK3_ENABLED
    module procedure setBernRandRUP_IK5_RK3
        use pm_kind, only: IKC => IK5, RKC => RK3
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && RK2_ENABLED
    module procedure setBernRandRUP_IK5_RK2
        use pm_kind, only: IKC => IK5, RKC => RK2
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && RK1_ENABLED
    module procedure setBernRandRUP_IK5_RK1
        use pm_kind, only: IKC => IK5, RKC => RK1
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK4_ENABLED && RK5_ENABLED
    module procedure setBernRandRUP_IK4_RK5
        use pm_kind, only: IKC => IK4, RKC => RK5
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && RK4_ENABLED
    module procedure setBernRandRUP_IK4_RK4
        use pm_kind, only: IKC => IK4, RKC => RK4
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && RK3_ENABLED
    module procedure setBernRandRUP_IK4_RK3
        use pm_kind, only: IKC => IK4, RKC => RK3
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && RK2_ENABLED
    module procedure setBernRandRUP_IK4_RK2
        use pm_kind, only: IKC => IK4, RKC => RK2
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && RK1_ENABLED
    module procedure setBernRandRUP_IK4_RK1
        use pm_kind, only: IKC => IK4, RKC => RK1
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK3_ENABLED && RK5_ENABLED
    module procedure setBernRandRUP_IK3_RK5
        use pm_kind, only: IKC => IK3, RKC => RK5
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && RK4_ENABLED
    module procedure setBernRandRUP_IK3_RK4
        use pm_kind, only: IKC => IK3, RKC => RK4
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && RK3_ENABLED
    module procedure setBernRandRUP_IK3_RK3
        use pm_kind, only: IKC => IK3, RKC => RK3
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && RK2_ENABLED
    module procedure setBernRandRUP_IK3_RK2
        use pm_kind, only: IKC => IK3, RKC => RK2
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && RK1_ENABLED
    module procedure setBernRandRUP_IK3_RK1
        use pm_kind, only: IKC => IK3, RKC => RK1
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK2_ENABLED && RK5_ENABLED
    module procedure setBernRandRUP_IK2_RK5
        use pm_kind, only: IKC => IK2, RKC => RK5
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && RK4_ENABLED
    module procedure setBernRandRUP_IK2_RK4
        use pm_kind, only: IKC => IK2, RKC => RK4
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && RK3_ENABLED
    module procedure setBernRandRUP_IK2_RK3
        use pm_kind, only: IKC => IK2, RKC => RK3
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && RK2_ENABLED
    module procedure setBernRandRUP_IK2_RK2
        use pm_kind, only: IKC => IK2, RKC => RK2
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && RK1_ENABLED
    module procedure setBernRandRUP_IK2_RK1
        use pm_kind, only: IKC => IK2, RKC => RK1
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK1_ENABLED && RK5_ENABLED
    module procedure setBernRandRUP_IK1_RK5
        use pm_kind, only: IKC => IK1, RKC => RK5
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && RK4_ENABLED
    module procedure setBernRandRUP_IK1_RK4
        use pm_kind, only: IKC => IK1, RKC => RK4
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && RK3_ENABLED
    module procedure setBernRandRUP_IK1_RK3
        use pm_kind, only: IKC => IK1, RKC => RK3
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && RK2_ENABLED
    module procedure setBernRandRUP_IK1_RK2
        use pm_kind, only: IKC => IK1, RKC => RK2
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && RK1_ENABLED
    module procedure setBernRandRUP_IK1_RK1
        use pm_kind, only: IKC => IK1, RKC => RK1
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED && RK5_ENABLED
    module procedure setBernRandRUP_LK5_RK5
        use pm_kind, only: LKC => LK5, RKC => RK5
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK5_ENABLED && RK4_ENABLED
    module procedure setBernRandRUP_LK5_RK4
        use pm_kind, only: LKC => LK5, RKC => RK4
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK5_ENABLED && RK3_ENABLED
    module procedure setBernRandRUP_LK5_RK3
        use pm_kind, only: LKC => LK5, RKC => RK3
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK5_ENABLED && RK2_ENABLED
    module procedure setBernRandRUP_LK5_RK2
        use pm_kind, only: LKC => LK5, RKC => RK2
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK5_ENABLED && RK1_ENABLED
    module procedure setBernRandRUP_LK5_RK1
        use pm_kind, only: LKC => LK5, RKC => RK1
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK4_ENABLED && RK5_ENABLED
    module procedure setBernRandRUP_LK4_RK5
        use pm_kind, only: LKC => LK4, RKC => RK5
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED && RK4_ENABLED
    module procedure setBernRandRUP_LK4_RK4
        use pm_kind, only: LKC => LK4, RKC => RK4
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED && RK3_ENABLED
    module procedure setBernRandRUP_LK4_RK3
        use pm_kind, only: LKC => LK4, RKC => RK3
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED && RK2_ENABLED
    module procedure setBernRandRUP_LK4_RK2
        use pm_kind, only: LKC => LK4, RKC => RK2
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED && RK1_ENABLED
    module procedure setBernRandRUP_LK4_RK1
        use pm_kind, only: LKC => LK4, RKC => RK1
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK3_ENABLED && RK5_ENABLED
    module procedure setBernRandRUP_LK3_RK5
        use pm_kind, only: LKC => LK3, RKC => RK5
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED && RK4_ENABLED
    module procedure setBernRandRUP_LK3_RK4
        use pm_kind, only: LKC => LK3, RKC => RK4
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED && RK3_ENABLED
    module procedure setBernRandRUP_LK3_RK3
        use pm_kind, only: LKC => LK3, RKC => RK3
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED && RK2_ENABLED
    module procedure setBernRandRUP_LK3_RK2
        use pm_kind, only: LKC => LK3, RKC => RK2
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED && RK1_ENABLED
    module procedure setBernRandRUP_LK3_RK1
        use pm_kind, only: LKC => LK3, RKC => RK1
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK2_ENABLED && RK5_ENABLED
    module procedure setBernRandRUP_LK2_RK5
        use pm_kind, only: LKC => LK2, RKC => RK5
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED && RK4_ENABLED
    module procedure setBernRandRUP_LK2_RK4
        use pm_kind, only: LKC => LK2, RKC => RK4
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED && RK3_ENABLED
    module procedure setBernRandRUP_LK2_RK3
        use pm_kind, only: LKC => LK2, RKC => RK3
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED && RK2_ENABLED
    module procedure setBernRandRUP_LK2_RK2
        use pm_kind, only: LKC => LK2, RKC => RK2
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED && RK1_ENABLED
    module procedure setBernRandRUP_LK2_RK1
        use pm_kind, only: LKC => LK2, RKC => RK1
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK1_ENABLED && RK5_ENABLED
    module procedure setBernRandRUP_LK1_RK5
        use pm_kind, only: LKC => LK1, RKC => RK5
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED && RK4_ENABLED
    module procedure setBernRandRUP_LK1_RK4
        use pm_kind, only: LKC => LK1, RKC => RK4
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED && RK3_ENABLED
    module procedure setBernRandRUP_LK1_RK3
        use pm_kind, only: LKC => LK1, RKC => RK3
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED && RK2_ENABLED
    module procedure setBernRandRUP_LK1_RK2
        use pm_kind, only: LKC => LK1, RKC => RK2
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED && RK1_ENABLED
    module procedure setBernRandRUP_LK1_RK1
        use pm_kind, only: LKC => LK1, RKC => RK1
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module procedure setBernRandRUP_RK5_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBernRandRUP_RK4_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBernRandRUP_RK3_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBernRandRUP_RK2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBernRandRUP_RK1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RUP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setBernRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines
