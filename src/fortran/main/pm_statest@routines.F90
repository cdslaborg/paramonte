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
!>  This file contains procedure implementations of [pm_statest](@ref pm_statest).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_statest) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_distKolm, only: setKolmCDF

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getProbKS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WIX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getProbKS_WIX_D0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getProbKS_WIX_D0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getProbKS_WIX_D0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getProbKS_WIX_D0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getProbKS_WIX_D0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WIX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WRX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getProbKS_WRX_D0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getProbKS_WRX_D0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getProbKS_WRX_D0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getProbKS_WRX_D0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getProbKS_WRX_D0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WRX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WII_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getProbKS_WII_D0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getProbKS_WII_D0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getProbKS_WII_D0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getProbKS_WII_D0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getProbKS_WII_D0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WII_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WRI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getProbKS_WRI_D0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getProbKS_WRI_D0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getProbKS_WRI_D0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getProbKS_WRI_D0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getProbKS_WRI_D0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WRI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WRR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getProbKS_WRR_D0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getProbKS_WRR_D0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getProbKS_WRR_D0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getProbKS_WRR_D0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getProbKS_WRR_D0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WRR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getProbKS_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setProbKS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WIX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setProbKS_WIX_D0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setProbKS_WIX_D0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setProbKS_WIX_D0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setProbKS_WIX_D0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setProbKS_WIX_D0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WIX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WRX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setProbKS_WRX_D0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setProbKS_WRX_D0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setProbKS_WRX_D0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setProbKS_WRX_D0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setProbKS_WRX_D0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WRX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WII_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setProbKS_WII_D0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setProbKS_WII_D0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setProbKS_WII_D0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setProbKS_WII_D0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setProbKS_WII_D0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WII_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WRI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setProbKS_WRI_D0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setProbKS_WRI_D0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setProbKS_WRI_D0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setProbKS_WRI_D0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setProbKS_WRI_D0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WRI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WRR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setProbKS_WRR_D0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setProbKS_WRR_D0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setProbKS_WRR_D0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setProbKS_WRR_D0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setProbKS_WRR_D0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_statest@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WRR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setProbKS_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines