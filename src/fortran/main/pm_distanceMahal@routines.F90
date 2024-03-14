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
!>  This file contains procedure implementations of [pm_distanceMahal](@ref pm_distanceMahal).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distanceMahal) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_matrixClass, only: isMatClass, posdefmat
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMahalSq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvDef_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMahalSqEleInvDef_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMahalSqEleInvDef_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMahalSqEleInvDef_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMahalSqEleInvDef_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMahalSqEleInvDef_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMahalSqEleInvDef_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMahalSqEleInvDef_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMahalSqEleInvDef_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMahalSqEleInvDef_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMahalSqEleInvDef_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvDef_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvCen_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMahalSqEleInvCen_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMahalSqEleInvCen_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMahalSqEleInvCen_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMahalSqEleInvCen_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMahalSqEleInvCen_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMahalSqEleInvCen_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMahalSqEleInvCen_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMahalSqEleInvCen_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMahalSqEleInvCen_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMahalSqEleInvCen_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvCen_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define One_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvDef_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMahalSqOneInvDef_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMahalSqOneInvDef_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMahalSqOneInvDef_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMahalSqOneInvDef_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMahalSqOneInvDef_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMahalSqOneInvDef_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMahalSqOneInvDef_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMahalSqOneInvDef_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMahalSqOneInvDef_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMahalSqOneInvDef_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvDef_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvCen_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMahalSqOneInvCen_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMahalSqOneInvCen_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMahalSqOneInvCen_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMahalSqOneInvCen_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMahalSqOneInvCen_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMahalSqOneInvCen_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMahalSqOneInvCen_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMahalSqOneInvCen_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMahalSqOneInvCen_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMahalSqOneInvCen_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvCen_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvDef_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMahalSqOneInvDef_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMahalSqOneInvDef_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMahalSqOneInvDef_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMahalSqOneInvDef_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMahalSqOneInvDef_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMahalSqOneInvDef_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMahalSqOneInvDef_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMahalSqOneInvDef_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMahalSqOneInvDef_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMahalSqOneInvDef_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvDef_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvCen_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMahalSqOneInvCen_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMahalSqOneInvCen_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMahalSqOneInvCen_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMahalSqOneInvCen_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMahalSqOneInvCen_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMahalSqOneInvCen_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMahalSqOneInvCen_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMahalSqOneInvCen_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMahalSqOneInvCen_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMahalSqOneInvCen_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvCen_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef One_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Mix_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvDef_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMahalSqMixInvDef_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMahalSqMixInvDef_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMahalSqMixInvDef_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMahalSqMixInvDef_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMahalSqMixInvDef_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMahalSqMixInvDef_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMahalSqMixInvDef_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMahalSqMixInvDef_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMahalSqMixInvDef_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMahalSqMixInvDef_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvDef_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvCen_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMahalSqMixInvCen_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMahalSqMixInvCen_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMahalSqMixInvCen_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMahalSqMixInvCen_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMahalSqMixInvCen_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMahalSqMixInvCen_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMahalSqMixInvCen_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMahalSqMixInvCen_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMahalSqMixInvCen_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMahalSqMixInvCen_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvCen_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvDef_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMahalSqMixInvDef_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMahalSqMixInvDef_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMahalSqMixInvDef_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMahalSqMixInvDef_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMahalSqMixInvDef_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMahalSqMixInvDef_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMahalSqMixInvDef_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMahalSqMixInvDef_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMahalSqMixInvDef_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMahalSqMixInvDef_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvDef_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvCen_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMahalSqMixInvCen_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMahalSqMixInvCen_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMahalSqMixInvCen_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMahalSqMixInvCen_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMahalSqMixInvCen_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMahalSqMixInvCen_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMahalSqMixInvCen_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMahalSqMixInvCen_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMahalSqMixInvCen_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMahalSqMixInvCen_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvCen_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Mix_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMahalSq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMahalSq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvDef_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMahalSqEleInvDef_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMahalSqEleInvDef_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMahalSqEleInvDef_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMahalSqEleInvDef_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMahalSqEleInvDef_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMahalSqEleInvDef_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMahalSqEleInvDef_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMahalSqEleInvDef_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMahalSqEleInvDef_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMahalSqEleInvDef_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvDef_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvCen_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMahalSqEleInvCen_D0_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMahalSqEleInvCen_D0_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMahalSqEleInvCen_D0_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMahalSqEleInvCen_D0_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMahalSqEleInvCen_D0_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMahalSqEleInvCen_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMahalSqEleInvCen_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMahalSqEleInvCen_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMahalSqEleInvCen_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMahalSqEleInvCen_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvCen_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define One_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvDef_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMahalSqOneInvDef_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMahalSqOneInvDef_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMahalSqOneInvDef_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMahalSqOneInvDef_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMahalSqOneInvDef_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMahalSqOneInvDef_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMahalSqOneInvDef_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMahalSqOneInvDef_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMahalSqOneInvDef_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMahalSqOneInvDef_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvDef_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvCen_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMahalSqOneInvCen_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMahalSqOneInvCen_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMahalSqOneInvCen_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMahalSqOneInvCen_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMahalSqOneInvCen_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMahalSqOneInvCen_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMahalSqOneInvCen_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMahalSqOneInvCen_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMahalSqOneInvCen_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMahalSqOneInvCen_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvCen_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvDef_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMahalSqOneInvDef_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMahalSqOneInvDef_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMahalSqOneInvDef_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMahalSqOneInvDef_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMahalSqOneInvDef_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMahalSqOneInvDef_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMahalSqOneInvDef_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMahalSqOneInvDef_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMahalSqOneInvDef_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMahalSqOneInvDef_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvDef_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvCen_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMahalSqOneInvCen_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMahalSqOneInvCen_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMahalSqOneInvCen_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMahalSqOneInvCen_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMahalSqOneInvCen_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMahalSqOneInvCen_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMahalSqOneInvCen_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMahalSqOneInvCen_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMahalSqOneInvCen_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMahalSqOneInvCen_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvCen_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef One_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Mix_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvDef_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMahalSqMixInvDef_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMahalSqMixInvDef_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMahalSqMixInvDef_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMahalSqMixInvDef_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMahalSqMixInvDef_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMahalSqMixInvDef_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMahalSqMixInvDef_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMahalSqMixInvDef_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMahalSqMixInvDef_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMahalSqMixInvDef_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvDef_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvCen_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMahalSqMixInvCen_D1_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMahalSqMixInvCen_D1_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMahalSqMixInvCen_D1_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMahalSqMixInvCen_D1_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMahalSqMixInvCen_D1_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMahalSqMixInvCen_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMahalSqMixInvCen_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMahalSqMixInvCen_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMahalSqMixInvCen_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMahalSqMixInvCen_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvCen_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvDef_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMahalSqMixInvDef_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMahalSqMixInvDef_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMahalSqMixInvDef_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMahalSqMixInvDef_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMahalSqMixInvDef_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMahalSqMixInvDef_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMahalSqMixInvDef_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMahalSqMixInvDef_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMahalSqMixInvDef_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMahalSqMixInvDef_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvDef_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define InvCen_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMahalSqMixInvCen_D2_CK5
        use pm_kind, only: CKC => CK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMahalSqMixInvCen_D2_CK4
        use pm_kind, only: CKC => CK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMahalSqMixInvCen_D2_CK3
        use pm_kind, only: CKC => CK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMahalSqMixInvCen_D2_CK2
        use pm_kind, only: CKC => CK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMahalSqMixInvCen_D2_CK1
        use pm_kind, only: CKC => CK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMahalSqMixInvCen_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMahalSqMixInvCen_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMahalSqMixInvCen_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMahalSqMixInvCen_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMahalSqMixInvCen_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceMahal@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef InvCen_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Mix_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMahalSq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines