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
!>  This file contains procedure implementations of [pm_distanceKolm](@ref pm_distanceKolm).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distanceKolm) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_arraySort, only: isAscending
    use pm_arraySort, only: setSorted
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDisKolm_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSSD_WDD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSSD_WDD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSSD_WDD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSSD_WDD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSSD_WDD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSSA_WDD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSSA_WDD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSSA_WDD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSSA_WDD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSSA_WDD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WID_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSSD_WID_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSSD_WID_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSSD_WID_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSSD_WID_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSSD_WID_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSSA_WID_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSSA_WID_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSSA_WID_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSSA_WID_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSSA_WID_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WID_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WRD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSSD_WRD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSSD_WRD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSSD_WRD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSSD_WRD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSSD_WRD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSSA_WRD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSSA_WRD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSSA_WRD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSSA_WRD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSSA_WRD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WRD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WII_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSSD_WII_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSSD_WII_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSSD_WII_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSSD_WII_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSSD_WII_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSSA_WII_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSSA_WII_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSSA_WII_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSSA_WII_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSSA_WII_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WII_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WRR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSSD_WRR_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSSD_WRR_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSSD_WRR_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSSD_WRR_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSSD_WRR_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSSA_WRR_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSSA_WRR_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSSA_WRR_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSSA_WRR_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSSA_WRR_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WRR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSXD_WDD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSXD_WDD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSXD_WDD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSXD_WDD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSXD_WDD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SXA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSXA_WDD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSXA_WDD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSXA_WDD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSXA_WDD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSXA_WDD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SXA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WID_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSXD_WID_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSXD_WID_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSXD_WID_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSXD_WID_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSXD_WID_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SXA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSXA_WID_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSXA_WID_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSXA_WID_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSXA_WID_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSXA_WID_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SXA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WID_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WRD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSXD_WRD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSXD_WRD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSXD_WRD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSXD_WRD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSXD_WRD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SXA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSXA_WRD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSXA_WRD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSXA_WRD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSXA_WRD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSXA_WRD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SXA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WRD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SCD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSCD_WDD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSCD_WDD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSCD_WDD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSCD_WDD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSCD_WDD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SCD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SCA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSCA_WDD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSCA_WDD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSCA_WDD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSCA_WDD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSCA_WDD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SCA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WID_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SCD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSCD_WID_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSCD_WID_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSCD_WID_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSCD_WID_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSCD_WID_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SCD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SCA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSCA_WID_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSCA_WID_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSCA_WID_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSCA_WID_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSCA_WID_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SCA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WID_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WRD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SCD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSCD_WRD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSCD_WRD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSCD_WRD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSCD_WRD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSCD_WRD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SCD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SCA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisKolmSCA_WRD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisKolmSCA_WRD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisKolmSCA_WRD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisKolmSCA_WRD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisKolmSCA_WRD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SCA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WRD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDisKolm_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setDisKolm_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSSD_WDD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSSD_WDD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSSD_WDD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSSD_WDD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSSD_WDD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSSA_WDD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSSA_WDD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSSA_WDD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSSA_WDD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSSA_WDD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WID_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSSD_WID_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSSD_WID_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSSD_WID_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSSD_WID_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSSD_WID_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSSA_WID_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSSA_WID_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSSA_WID_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSSA_WID_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSSA_WID_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WID_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WRD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSSD_WRD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSSD_WRD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSSD_WRD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSSD_WRD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSSD_WRD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSSA_WRD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSSA_WRD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSSA_WRD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSSA_WRD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSSA_WRD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WRD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WII_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSSD_WII_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSSD_WII_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSSD_WII_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSSD_WII_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSSD_WII_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSSA_WII_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSSA_WII_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSSA_WII_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSSA_WII_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSSA_WII_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WII_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WRR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSSD_WRR_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSSD_WRR_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSSD_WRR_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSSD_WRR_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSSD_WRR_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SSA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSSA_WRR_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSSA_WRR_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSSA_WRR_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSSA_WRR_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSSA_WRR_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SSA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WRR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSXD_WDD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSXD_WDD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSXD_WDD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSXD_WDD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSXD_WDD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SXA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSXA_WDD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSXA_WDD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSXA_WDD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSXA_WDD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSXA_WDD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SXA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WID_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSXD_WID_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSXD_WID_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSXD_WID_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSXD_WID_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSXD_WID_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SXA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSXA_WID_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSXA_WID_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSXA_WID_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSXA_WID_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSXA_WID_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SXA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WID_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WRD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSXD_WRD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSXD_WRD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSXD_WRD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSXD_WRD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSXD_WRD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SXA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSXA_WRD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSXA_WRD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSXA_WRD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSXA_WRD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSXA_WRD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SXA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WRD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SCD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSCD_WDD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSCD_WDD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSCD_WDD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSCD_WDD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSCD_WDD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SCD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SCA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSCA_WDD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSCA_WDD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSCA_WDD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSCA_WDD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSCA_WDD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SCA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WID_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SCD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSCD_WID_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSCD_WID_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSCD_WID_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSCD_WID_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSCD_WID_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SCD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SCA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSCA_WID_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSCA_WID_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSCA_WID_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSCA_WID_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSCA_WID_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SCA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WID_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WRD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SCD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSCD_WRD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSCD_WRD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSCD_WRD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSCD_WRD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSCD_WRD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SCD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SCA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDisKolmSCA_WRD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDisKolmSCA_WRD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDisKolmSCA_WRD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDisKolmSCA_WRD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDisKolmSCA_WRD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distanceKolm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SCA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WRD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setDisKolm_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines